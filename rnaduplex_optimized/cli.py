"""RNAduplex-compatible streaming and identity-preserving batch command."""
import os
import sys
from . import __version__, COMPATIBILITY_VERSION

REFERENCE=os.path.join(os.path.dirname(__file__),'_RNAduplex_2_7_2')


def parser():
    import argparse
    try:
        from . import _native
        revision=_native.build_info()['revision']
    except (ImportError,AttributeError):
        revision='unbuilt'
    p=argparse.ArgumentParser(prog='RNAduplex_optimized',description=(
        'Exact RNAduplex 2.7.2 predictions. Default MFE uses the native core; '
        'model/suboptimal options use the packaged unchanged 2.7.2 backend.'),
        epilog='Examples: printf "ACGU\\nUGCA\\n" | RNAduplex_optimized\n'
               'RNAduplex_optimized batch --input pairs.jsonl.gz --output results.jsonl.gz --threads 2',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--version','-V',action='version',version=f'RNAduplex_optimized {__version__}; implementation {revision}; ViennaRNA compatibility {COMPATIBILITY_VERSION}')
    p.add_argument('--buffered',action='store_true',help='buffer optimized default text results (default flushes each result)')
    p.add_argument('--reference-help',action='store_true',help='show all model/suboptimal flags of the unchanged backend')
    p.add_argument('--temp','-T',metavar='C',help='temperature; unchanged 2.7.2 backend')
    p.add_argument('--paramFile','-P',metavar='FILE',help='parameter file or DNA; unchanged backend')
    p.add_argument('--salt',metavar='M',help='salt option with original 2.7.2 behavior; unchanged backend')
    p.add_argument('--dangles','-d',metavar='N',help='dangling-end model flag; unchanged backend')
    p.add_argument('--deltaEnergy','-e',metavar='KCAL',help='suboptimal range; unchanged backend')
    p.add_argument('--sorted','-s',action='store_true',help='sort suboptimal results; unchanged backend')
    p.add_argument('--noLP',action='store_true',help='original noLP flag; unchanged backend')
    p.add_argument('--noGU',action='store_true',help='original noGU flag; unchanged backend')
    p.add_argument('--noClosingGU',action='store_true',help='original GU closure flag; unchanged backend')
    p.add_argument('--noconv',action='store_true',help='original T conversion flag; unchanged backend')
    p.add_argument('mode',nargs='?',choices=['batch'],help='omit for stdin sequence pairs (one sequence per line); FASTA headers accepted')
    return p


def text_input():
    # RNAlib removes LF only; CR is a sequence character in the original CLI.
    # Reading bytes avoids Python's universal-newline conversion.
    for raw in sys.stdin.buffer:
        line=raw.rstrip(b'\n')
        if not line or line.startswith(b'*'):continue
        if line.startswith(b'@'):return
        # Preserve the distinction between initially blank (skipped above) and
        # whitespace-only (upstream trims it into an empty strand after skipping).
        line=line.rstrip(b' \t') or b' '
        if line.startswith(b'>'):
            words=line[1:].split()
            if not words:return
            sys.stdout.buffer.write(b'>'+words[0]+b'\n')
            sys.stdout.buffer.flush()
        else:yield line


def stream(buffered):
    from . import _native
    lines=iter(text_input())
    while True:
        a=next(lines,None)
        if a is None:return 0
        b=next(lines,None)
        if b is None:return 0
        if a.isalpha() and b.isalpha() and a.isascii() and b.isascii():
            result=_native.predict(a.upper().replace(b'T',b'U'),b.upper().replace(b'T',b'U'))
            sys.stdout.buffer.write(result['text'].encode('ascii'))
        else:
            # Preserve unusual upstream input behavior in an isolated subprocess.
            # A crash is an error, never a successful native prediction.
            import subprocess
            print('RNAduplex_optimized: backend unchanged-2.7.2 for nonletter text input',file=sys.stderr)
            p=subprocess.run([str(REFERENCE)],input=a+b'\n'+b+b'\n',capture_output=True)
            if p.stderr:sys.stderr.buffer.write(p.stderr)
            if p.returncode:raise RuntimeError(f'unchanged backend exited {p.returncode} on nonletter input')
            sys.stdout.buffer.write(p.stdout)
        if not buffered:sys.stdout.flush()


def main(argv=None):
    argv=list(sys.argv[1:] if argv is None else argv)
    try:
        if argv and argv[0]=='batch':
            import argparse
            from .batch import run_batch
            p=argparse.ArgumentParser(prog='RNAduplex_optimized batch',description='Expand supplied RNA pairs without strand transformations; emit every occurrence in physical-line/canonical-group order.')
            p.add_argument('--input',required=True,help='JSONL or JSONL.gz with nested pairs and left_rna/right_rna')
            p.add_argument('--output',default='-',help='JSONL/JSONL.gz output or - for stdout')
            p.add_argument('--threads',type=int,default=1,help='independent worker processes; default 1')
            p.add_argument('--chunk-size',type=int,default=128,help='pairs per bounded worker task; default 128')
            p.add_argument('--backend',choices=['optimized','reference'],default='optimized')
            args,reference_args=p.parse_known_args(argv[1:])
            run_batch(args,reference_args);return 0
        if '--reference-help' in argv:
            import subprocess
            return subprocess.call([str(REFERENCE),'--detailed-help'])
        if any(x in argv for x in ('--help','-h','--version','-V')):
            parser().parse_args(argv);return 0
        buffered='--buffered' in argv
        reference_args=[x for x in argv if x!='--buffered']
        # Upstream colors terminal stdout even when stdin is a pipe. Its
        # unchanged backend independently prompts only when both are terminals.
        if reference_args or sys.stdout.isatty():
            import subprocess
            if buffered:
                print('RNAduplex_optimized: --buffered applies only to optimized default text; unchanged backend retains its own buffering',file=sys.stderr)
            print('RNAduplex_optimized: backend unchanged-2.7.2',file=sys.stderr)
            return subprocess.call([str(REFERENCE),*reference_args])
        return stream(buffered)
    except (ValueError,RuntimeError,OSError) as error:
        print(f'RNAduplex_optimized: {error}',file=sys.stderr);return 2
    except BrokenPipeError:return 0
