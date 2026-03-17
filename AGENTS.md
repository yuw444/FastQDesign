# PROJECT KNOWLEDGE BASE

**Generated:** 2026-03-13 **Commit:** a912cd3 **Branch:** main

## OVERVIEW

R package for scRNA-seq experiment design. Guides investigators in
designing single-cell RNA sequencing experiments by shifting focus from
UMI matrix back to raw FastQ reads. Includes C tool `fastF` for
pseudo-design dataset generation.

## STRUCTURE

    FastQDesign/
    ├── R/              # R package functions (8 files)
    ├── src/            # C source for fastF tool (11 files)
    ├── tests/          # testthat tests
    ├── data/           # package data
    ├── man/            # R documentation
    ├── CMakeLists.txt  # fastF build config
    └── NAMESPACE       # 16 exported functions

## WHERE TO LOOK

| Task             | Location                   | Notes                   |
|------------------|----------------------------|-------------------------|
| Main R functions | `R/`                       | design.R largest (15kb) |
| fastF C tool     | `src/`                     | main.c entry point      |
| Tests            | `tests/testthat/`          | 3 test files            |
| Package metadata | `NAMESPACE`, `DESCRIPTION` | roxygen2 generated      |

## CODE MAP

| Symbol         | Type     | Location       | Role                     |
|----------------|----------|----------------|--------------------------|
| FastQDesign    | function | R/design.R     | Main design optimization |
| SamplePrep     | function | R/design.R     | Reference preparation    |
| SampleMatch    | function | R/design.R     | Sample comparison        |
| DownSample     | function | R/downsample.R | Downsampling logic       |
| FindPseudotime | function | R/pseudotime.R | Pseudotime analysis      |
| fastF          | binary   | src/main.c     | C entry point            |

## CONVENTIONS

- **Roxygen2**: Documentation generated via roxygen2 (NAMESPACE
  auto-generated)
- **R CMD build**: Standard R package build
- **CMake**: C tool build (Debug by default, requires htslib, sqlite)
- **Pipe operator**: Uses `%>%` from magrittr

## ANTI-PATTERNS (THIS PROJECT)

- No explicit anti-patterns found in code comments

## UNIQUE STYLES

- Heavy dependency on Seurat, monocle3, dplyr for single-cell analysis
- Mixed R+C architecture (fastF as separate C tool)
- Power analysis integration for experiment design

## COMMANDS

``` bash
# R package build/install
R CMD build .
R CMD install FastQDesign_*.tar.gz

# Run tests
Rscript -e "devtools::test()"

# Build fastF C tool
cd /path/to/FastQDesign
mkdir build && cd build
cmake ..
make
```

## NOTES

- fastF subproject has separate CMake build (not integrated with R CMD
  build)
- Depends on: glmGamPoi, Seurat, monocle3, ggplot2, dplyr, patchwork
- Published: Nature Communications Biology (PMID: 40175506)
