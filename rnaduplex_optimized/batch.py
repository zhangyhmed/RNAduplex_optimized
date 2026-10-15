"""Ordered, bounded, occurrence-preserving JSONL batch processing."""
import concurrent.futures
import gzip
import hashlib
import itertools
import json
import re
import subprocess
import sys
from collections import deque
from pathlib import Path

PAIR_NAMES=('repeat_pair','stem_internal_pair','repeat_stem_internal_pair',
            'stem_external_pair','repeat_stem_external_pair')
PATTERN=re.compile(r'^([().&]+)\s+(-?\d+),(-?\d+)\s*:\s*(-?\d+),(-?\d+)\s+\(\s*(-?\d+\.\d{2})\)\s*$')
REFERENCE=Path(__file__).with_name('_RNAduplex_2_7_2')


def open_text(path,mode):
    return gzip.open(path,mode,compresslevel=1) if str(path).endswith('.gz') else open(path,mode)


def input_hash(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def expanded_records(path):
    digest=input_hash(path)
    with open_text(path,'rt') as f:
        for ordinal,line in enumerate(f,1):
            record=json.loads(line)
            if not isinstance(record,dict):
                raise ValueError(f'line {ordinal}: record must be an object')
            if not isinstance(record.get('pairs'),dict):
                raise ValueError(f'line {ordinal}: missing/malformed pairs object')
            for field in ('junction_id','tract_id'):
                if field not in record:raise ValueError(f'line {ordinal}: missing {field}')
            unexpected=set(record['pairs'])-set(PAIR_NAMES)
            if unexpected:raise ValueError(f'line {ordinal}: unexpected groups {sorted(unexpected)}')
            missing=[name for name in PAIR_NAMES if name not in record['pairs']]
            if missing:
                print(f'RNAduplex_optimized: line {ordinal}: missing pair groups: '+', '.join(missing),file=sys.stderr)
            metadata={k:v for k,v in record.items() if k!='pairs'}
            for name in PAIR_NAMES:
                if name not in record['pairs']:continue
                pair=record['pairs'][name]
                for field in ('left_rna','right_rna'):
                    if not isinstance(pair,dict) or not isinstance(pair.get(field),str):
                        raise ValueError(f'line {ordinal} {name}: missing/malformed {field}')
                left,right=pair['left_rna'],pair['right_rna']
                for field,sequence in [('left_rna',left),('right_rna',right)]:
                    if sequence and (not sequence.isascii() or not sequence.isalpha()):
                        raise ValueError(f'line {ordinal} {name}: {field} must contain ASCII letters only')
                # Do not reconstruct, swap, or transform supplied ordered RNA strings.
                identity=dict(key=f'{digest}:{ordinal}:{name}',record_ordinal=ordinal,
                              pair_name=name,junction_id=record['junction_id'],tract_id=record['tract_id'],
                              metadata=metadata,left_length=len(left),right_length=len(right),
                              ordered_sequence_sha256=hashlib.sha256((left+'\0'+right).encode()).hexdigest(),
                              missing_groups=missing)
                yield identity,left,right


def parse_prediction(line):
    match=PATTERN.fullmatch(line)
    if not match:raise ValueError(f'unchanged backend produced malformed prediction: {line!r}')
    return dict(structure=match[1],left_start=int(match[2]),left_end=int(match[3]),
                right_start=int(match[4]),right_end=int(match[5]),energy=match[6],text=line+'\n')


def predict_chunk(payload):
    chunk,reference_args=payload
    results=[]
    if reference_args is None:
        from . import _native
        for identity,left,right in chunk:
            result=dict(identity,backend='optimized-native',status='ok',predictions=[])
            if not left or not right:result['status']='empty_sequence'
            else:result['predictions']=[_native.predict(left.upper(),right.upper())]
            results.append(result)
    else:
        # Header boundaries retain multiplicity and result order for suboptimal output.
        text=''.join(f'>__rdo_{i}\n{a}\n{b}\n' for i,(_,a,b) in enumerate(chunk) if a and b)
        p=subprocess.run([str(REFERENCE),*reference_args],input=text,text=True,capture_output=True)
        if p.returncode:raise RuntimeError(f'unchanged backend exited {p.returncode}: {p.stderr.strip()}')
        if p.stderr:sys.stderr.write(p.stderr)
        predictions={};index=None
        for line in p.stdout.splitlines():
            if line.startswith('>__rdo_'):
                index=int(line[len('>__rdo_'):]);predictions[index]=[]
            elif line:
                if index is None:raise RuntimeError('unchanged backend output lacks identity header')
                predictions[index].append(parse_prediction(line))
        for i,(identity,left,right) in enumerate(chunk):
            result=dict(identity,backend='unchanged-2.7.2',status='ok',predictions=predictions.get(i,[]))
            if not left or not right:result['status']='empty_sequence'
            elif not result['predictions']:result['status']='no_prediction'
            results.append(result)
    return results


def chunks(records,size):
    records=iter(records)
    while chunk:=list(itertools.islice(records,size)):yield chunk


def run_batch(args,reference_args):
    if args.threads<1 or args.chunk_size<1:raise ValueError('threads and chunk-size must be positive')
    if args.output!='-':
        source,target=Path(args.input),Path(args.output)
        if source.resolve()==target.resolve() or (target.exists() and source.samefile(target)):
            raise ValueError('input and output must be different files')
    if args.backend=='reference' or reference_args:
        backend_args=reference_args
        print('RNAduplex_optimized: batch backend unchanged-2.7.2',file=sys.stderr)
    else:backend_args=None
    encoder=json.JSONEncoder(separators=(',',':'))
    jobs=((chunk,backend_args) for chunk in chunks(expanded_records(args.input),args.chunk_size))
    output=sys.stdout if args.output=='-' else open_text(args.output,'wt')
    try:
        if args.threads==1:
            for job in jobs:
                output.write(''.join(encoder.encode(row)+'\n' for row in predict_chunk(job)))
        else:
            # At most 2*workers chunks are live; reassemble strictly in submission order.
            with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as pool:
                pending=deque()
                for job in jobs:
                    pending.append(pool.submit(predict_chunk,job))
                    if len(pending)>=2*args.threads:
                        output.write(''.join(encoder.encode(row)+'\n' for row in pending.popleft().result()))
                while pending:
                    output.write(''.join(encoder.encode(row)+'\n' for row in pending.popleft().result()))
    finally:
        if output is not sys.stdout:output.close()
