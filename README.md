# RNAduplex_optimized

RNAduplex_optimized is a native RNAduplex implementation rewritten and optimized
from ViennaRNA 2.7.2 as a dependency for **CalibRNA**. It preserves the reference
thermodynamic model, integer energies, selected structures, coordinates, and
deterministic tie-breaking.

## What is faster

- **Compatible-pair indexing:** The native search visits compatible base-pair positions directly, avoiding impossible transitions and unnecessary work for sequences without compatible pairs.
- **Specialized energy evaluation:** Inlined energy calculations reuse shared terms and reduce repeated work while keeping the original integer arithmetic.
- **Exact dynamic-programming bounds:** Bounds derived from the thermodynamic parameters skip transitions that cannot improve the minimum. Traceback retains the original search order and tie-breaking.
- **Reusable memory:** Each worker reuses its dynamic-programming, sequence-encoding, and result buffers across predictions, reducing allocation overhead.
- **Batch serialization:** A reusable JSON encoder and bounded, combined writes reduce output overhead. Independent worker processes can evaluate pairs in parallel while preserving output order.
- **Native compilation:** Link-time optimization, loop unrolling, and host CPU tuning improve native execution without unsafe floating-point optimizations.

## Installation

Requires Python 3.10+, a GCC-compatible C compiler, and an existing
**ViennaRNA 2.7.2** installation containing `bin/RNAduplex`,
`include/ViennaRNA/`, and the static library `lib/libRNA.a`.

Activate the intended Python or Conda environment, then run:

```bash
git clone https://github.com/zhangyhmed/RNAduplex_optimized.git
cd RNAduplex_optimized
VIENNA_PREFIX=/path/to/viennarna-prefix CC=gcc python -m pip install .
```

If ViennaRNA is in the active Conda environment, use `VIENNA_PREFIX="$CONDA_PREFIX"`.
Compilation builds both the optimized extension and an unchanged compatibility
backend; it does not replace the original `RNAduplex` installation.

Builds target the build CPU by default. Set `RDO_NATIVE_FLAGS="-O3 -flto"` for a
generic instruction set. Wheels require compatible C/C++ and OpenMP runtime
libraries and retain the configured ViennaRNA prefix's library search path.

## Usage

```bash
RNAduplex_optimized --help
printf 'ACGU\nUGCA\n' | RNAduplex_optimized
RNAduplex_optimized batch --input sequence_pairs.jsonl.gz --output results.jsonl.gz --threads 2
```

Text input contains two sequence lines per pair. Batch input accepts JSONL or
gzip JSONL records with a `pairs` object whose entries contain ordered
`left_rna` and `right_rna` strings. `--threads` selects worker processes;
results retain input order, without deduplication or prediction caching.

The optimized backend handles the default minimum-free-energy model. Original
model options and suboptimal output use the bundled unchanged 2.7.2 backend;
`RNAduplex_optimized --reference-help` lists those options. The native interface
accepts nonempty ASCII-letter strands of at most 32,760 bases each; batch output
records empty strands explicitly.

## Attribution and license

Derived from ViennaRNA 2.7.2 by Ivo L. Hofacker and the ViennaRNA contributors,
Institute for Theoretical Chemistry, University of Vienna. See
[NOTICE](NOTICE) for provenance and [COPYING.ViennaRNA](COPYING.ViennaRNA) for
the upstream license and redistribution terms.
