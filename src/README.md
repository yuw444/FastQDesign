# fastF Source Code

C tool for generating pseudo-design datasets from FastQ references.

## Files

| File | Description |
|------|-------------|
| main.c | Entry point, CLI commands |
| filter.h/.c | FASTQ parsing, whitelist filtering, binary tree |
| bam2db.h/.c | BAM to SQLite DB conversion, DNA encoding |
| utils.h/.c | Sampling, sequencing, BAM extraction functions |
| argparse.h/.c | CLI argument parsing |
| hashtable.h/.c | Hash table implementation |
| mt19937ar.h/.c | Mersenne Twister RNG |

## Build

```bash
mkdir -p build && cd build
ml cmake
cmake ..
make
```

## Test

```bash
cd build
ml cmake
cmake ..
make
ctest --output-on-failure
```

Test results: 4/4 test suites passing (test_filter, test_bam2db, test_utils, test_integration)
