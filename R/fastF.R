#' @keywords internal
"_PACKAGE"

get_fastF_bin <- function() {
  # First check in build directory (for development)
  dev_bin <- "/scratch/g/chlin/Yu/FastQDesign/build/fastF"
  if (file.exists(dev_bin)) {
    return(dev_bin)
  }
  
  # Then check system PATH
  bin <- Sys.which("fastF")
  if (bin == "") {
    # Finally check package inst/bin
    bin <- system.file("bin", "fastF", package = "FastQDesign")
    if (bin == "") {
      stop("fastF binary not found. Please reinstall FastQDesign package.")
    }
  }
  bin
}

get_htslib_path <- function() {
  # Try to find htslib via multiple methods
  # 1. Try pkg-config if available (suppress warnings)
  pc_path <- ""
  if (Sys.which("pkg-config") != "") {
    pc_path <- tryCatch({
      suppressWarnings(system("pkg-config --variable=libdir htslib", intern = TRUE))
    }, error = function(e) "")
  }
  
  # 2. Check common htslib locations
  std_paths <- c(
    "/hpc/apps/htslib/1.22.1/lib",
    "/usr/lib", "/usr/local/lib",
    "/usr/lib64", "/usr/local/lib64",
    "/opt/htslib/lib", "/opt/lib"
  )
  
  # Use first valid path found
  paths <- c(pc_path, std_paths)
  for (p in paths) {
    if (nzchar(p) && dir.exists(p) && file.exists(file.path(p, "libhts.so"))) {
      return(p)
    }
  }
  
  # Fallback - let system find it via LD_LIBRARY_PATH
  NULL
}

run_fastF <- function(args, must_succeed = TRUE) {
  bin <- get_fastF_bin()
  
  # Build environment - use system() for better control
  htslib <- get_htslib_path()
  
  if (!is.null(htslib)) {
    # Use system() with explicit LD_LIBRARY_PATH
    cmd <- sprintf("LD_LIBRARY_PATH=%s %s %s", htslib, bin, paste(args, collapse = " "))
    result <- system(cmd, intern = TRUE, ignore.stderr = FALSE)
    status <- attr(result, "status")
  } else {
    # No htslib path, just run directly
    cmd <- sprintf("%s %s", bin, paste(args, collapse = " "))
    result <- system(cmd, intern = TRUE, ignore.stderr = FALSE)
    status <- attr(result, "status")
  }
  
  if (must_succeed && !is.null(status) && status != 0) {
    stop(paste(result, collapse = "\n"))
  }
  
  list(output = result, status = status)
}

#' Find Cell Barcode Whitelist from FastQ
#'
#' @param R1 Path to R1 fastq files (required)
#' @param out Output directory for whitelist (default: ".")
#' @param len Length of cell barcode (default: 16)
#' @param umi Length of UMI (default: 10)
#' @return Invisible NULL. Creates whitelist.txt in output directory.
#' @export
fastF_freq <- function(R1, out = ".", len = 16L, umi = 10L) {
  if (!file.exists(R1)) stop("R1 file does not exist: ", R1)
  
  args <- c("freq",
            "-R", R1,
            "-o", out,
            "-l", as.integer(len),
            "-u", as.integer(umi))
  
  run_fastF(args)
  invisible(NULL)
}

#' Filter FastQ Files Using Whitelist
#'
#' @param R1 Path to R1 fastq files (required)
#' @param R2 Path to R2 fastq files (optional)
#' @param I1 Path to I1 fastq files (optional)
#' @param out Output directory (default: ".")
#' @param whitelist Path to whitelist file (required unless allcells=TRUE)
#' @param len Length of cell barcode (default: 16)
#' @param seed Random seed (default: 926)
#' @param rate Rate of reads to keep (default: 0)
#' @param allcells Keep all reads with cell barcode (default: FALSE)
#' @return Invisible NULL. Creates I1.fastq.gz, R1.fastq.gz, R2.fastq.gz in output.
#' @export
fastF_filter <- function(R1, R2 = NULL, I1 = NULL, out = ".",
                         whitelist = NULL, len = 16L, seed = 926L,
                         rate = 0, allcells = FALSE) {
  if (!file.exists(R1)) stop("R1 file does not exist: ", R1)
  if (is.null(whitelist) && !allcells) {
    stop("Either whitelist or allcells=TRUE is required")
  }
  
  args <- c("filter",
            "-R", R1,
            "-o", out,
            "-l", as.integer(len),
            "-s", as.integer(seed),
            "-t", as.numeric(rate))
  
  if (!is.null(I1)) args <- c(args, "-I", I1)
  if (!is.null(R2)) args <- c(args, "-r", R2)
  if (!is.null(whitelist)) args <- c(args, "-w", whitelist)
  if (allcells) args <- c(args, "-a")
  
  run_fastF(args)
  invisible(NULL)
}

#' Extract CR/CB Tags from BAM
#'
#' @param bam Path to BAM file (required)
#' @param out Output TSV file (default: "crb_output.tsv")
#' @return Invisible NULL. Creates TSV file with CB->CR->count mapping.
#' @export
fastF_crb <- function(bam, out = "crb_output.tsv") {
  args <- c("crb",
            "-b", bam,
            "-o", out)
  
  run_fastF(args)
  invisible(NULL)
}

#' Convert BAM to SQLite UMI Matrix
#'
#' @param bam Path to BAM file (required)
#' @param feature Path to feature list file (required)
#' @param barcode Path to barcode list file (required)
#' @param out Output directory (default: ".")
#' @param dbname Database name (default: ":memory:")
#' @param cell Rate of cell barcode (default: 1.0)
#' @param depth Rate of depth (default: 1.0)
#' @param seed Random seed (default: 926)
#' @param umicopies Whether to store umi.tsv.gz (default: FALSE)
#' @return Invisible NULL. Creates SQLite DB in output directory.
#' @export
fastF_bam2db <- function(bam, feature, barcode, out = ".",
                          dbname = ":memory:", cell = 1.0, depth = 1.0,
                          seed = 926L, umicopies = FALSE) {
  if (!file.exists(feature)) stop("Feature file does not exist: ", feature)
  if (!file.exists(barcode)) stop("Barcode file does not exist: ", barcode)
  
  args <- c("bam2db",
            "-b", bam,
            "-f", feature,
            "-a", barcode,
            "-o", out,
            "-d", dbname,
            "-c", as.numeric(cell),
            "-r", as.numeric(depth),
            "-s", as.integer(seed))
  
  if (umicopies) args <- c(args, "-u")
  
  run_fastF(args)
  invisible(NULL)
}

#' Extract Custom Tags from BAM
#'
#' @param bam Path to BAM file (required)
#' @param tag Tag to extract, e.g., "CB" (required)
#' @param type Tag type: 0=string, 1=integer (default: 0)
#' @return Invisible NULL. Creates tag_summary.csv in current directory.
#' @export
fastF_extract <- function(bam, tag, type = 0L) {
  if (is.null(tag)) stop("--tag is required")
  
  args <- c("extract",
            "-b", bam,
            "-t", tag,
            "-T", as.integer(type))
  
  run_fastF(args)
  invisible(NULL)
}
