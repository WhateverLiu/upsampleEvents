

getRcppRelatedIncludePath <- function(pkgNames)
{
  ps = .libPaths()
  sapply(pkgNames, function(pkgName)
  {
    for (p in ps)
    {
      rst = paste0(p, '/', pkgName, '/include')
      if (dir.exists(rst)) break
    }; rst
  })
}




#'
#' Compile C++ script and load functions into environment
#' 
#' Built upon Rcpp, the function makes the process of compiling and linking
#' C++ scripts simpler and more flexible. Users supply compiler flags,
#' library include paths, shared library paths instead of relying on the
#' \code{R CMD} toolchain. The function also accepts sanitizer flags.
#' Compiling with sanitizers requires a new R session invoked on Linux like this: 
#' \code{LD_PRELOAD="/lib64/libasan.so.5 /lib64/libubsan.so.1"  R}
#' 
#' @param file  Path to the C++ script.
#' 
#' @param env  What environment should the function be loaded into. Default globalenv().
#' 
#' @param compilerPath  Full path to the compiler. If the compiler does not auto
#' search locations that contain standard headers, add the paths in \code{includePaths}.
#' 
#' @param optFlag  Compiler optmization flag. Default '-O2'.
#' 
#' @param flags  Compiler flags to be added.
#' 
#' @param includePaths  A vector of paths to C++ libraries to be included.
#' Libraries referenced using \code{// [[Rcpp::depends(XXX)]]} do not need to
#' included.
#' 
#' @param dllLinkFilePaths  A vector of paths to shared library to be linked.
#' Internally, the function will use \code{-L} to link the directories of 
#' these libraries, and then use \code{-l} to link the shared libraries without
#'  their extension names.
#'  
#' @param sanitize  A boolean. \code{TRUE} will add the sanitizer flags.
#' 
#' @param sanitizerFlags  A string. Sanitizer flags.
#' 
#' @param rebuild  A boolean. \code{TRUE} enforces a rebuild.
#' 
#' @param cacheDir  Where should the shared libraries be stored. Default
#' \code{'../tempFiles/CharlieRcpp'}. Folder will be created if nonexistent.
#' 
#' @param verbose  \code{TRUE} prints the full build command.
#' 
#' @return A list. The build context returned from calling 
#' \code{.Call("sourceCppContext", PACKAGE = "Rcpp", file, NULL, rebuild, cacheDir, .Platform)}.
#'
#'
CharlieSourceCpp <- function (
    file = "", 
    # code = NULL, 
    env = globalenv(), 
    
    compilerPath = Sys.which('g++'),
    
    optFlag = '-O2',
    
    flags = '-std=gnu++17 -shared -DNDEBUG -Wall -fpic -m64 -march=native -mfpmath=sse -msse2 -mstackrealign',
    
    includePaths = c(
      R.home('include'), 
      getRcppRelatedIncludePath('Rcpp')
      # getRcppRelatedIncludePath('RcppEigen')
    ),
    
    dllLinkFilePaths = c(
      # R.home('lib')
    ),
    
    sanitize = FALSE,
    sanitizerFlags = '-g -fno-omit-frame-pointer -fsanitize=address,undefined',
    
    
    # embeddedR = TRUE, 
    rebuild = FALSE, 
    cacheDir = '../tempFiles',
    # cleanupCacheDir = FALSE, 
    verbose = TRUE
    # dryRun = FALSE
    # windowsDebugDLL = FALSE, 
    # echo = TRUE
    )
{
  # rebuild = rebuild || sanitize
  if (sanitize) cacheDir = paste0(cacheDir, '/sanitized')
  dir.create(cacheDir, showWarnings = F, recursive = TRUE)
  cacheDir <- path.expand(cacheDir)
  cacheDir <- Rcpp:::.sourceCppPlatformCacheDir(cacheDir)
  cacheDir <- normalizePath(cacheDir)
  # if (!missing(code)) {
  #   rWorkingDir <- getwd()
  #   file <- tempfile(fileext = ".cpp", tmpdir = cacheDir)
  #   con <- file(file, open = "w")
  #   writeLines(code, con)
  #   close(con)
  # }
  # else {
  rWorkingDir <- normalizePath(dirname(file))
  # }
  file <- normalizePath(file, winslash = "/")
  
  
  # if (!tools::file_ext(file) %in% c("cc", "cpp")) {
  #   stop("The filename '", basename(file), "' does not have an ", 
  #        "extension of .cc or .cpp so cannot be compiled.")
  # }
  
  
  # if (.Platform$OS.type == "windows") {
  #   if (grepl(" ", basename(file), fixed = TRUE)) {
  #     stop("The filename '", basename(file), "' contains spaces. This ", 
  #          "is not permitted.")
  #   }
  # }
  # else {
  #   if (windowsDebugDLL) {
  #     if (verbose) {
  #       message("The 'windowsDebugDLL' toggle is ignored on ", 
  #               "non-Windows platforms.")
  #     }
  #     windowsDebugDLL <- FALSE
  #   }
  # }
  
  
  
  context <- .Call("sourceCppContext", PACKAGE = "Rcpp", file, 
                   NULL, rebuild, cacheDir, .Platform)
  
  
  # ============================================================================
  # Modify the cpp file to be compiled.
  # ============================================================================
  srcFile = paste0(context$buildDirectory, '/', context$cppSourceFilename)
  srcFileContent = readLines(srcFile)
  srcFileContent = srcFileContent[-(1:length(readLines(context$cppSourcePath)))]
  srcFileContent = c(paste0('#include "', context$cppSourcePath, '"'), srcFileContent)
  writeLines(srcFileContent, srcFile)
  
  
  if (context$buildRequired || rebuild) 
  {
    # if (verbose) Rcpp:::.printVerboseOutput(context)
    succeeded <- FALSE
    output <- NULL
    depends <- Rcpp:::.getSourceCppDependencies(context$depends, file)
    Rcpp:::.validatePackages(depends, context$cppSourceFilename)
    envRestore <- Rcpp:::.setupBuildEnvironment(depends, context$plugins, file)
    cwd <- getwd()
    setwd(context$buildDirectory)
    fromCode = FALSE
    if (!Rcpp:::.callBuildHook(context$cppSourcePath, fromCode, verbose)) 
    {
      Rcpp:::.restoreEnvironment(envRestore)
      setwd(cwd)
      return(invisible(NULL))
    }
    on.exit({
      if (!succeeded) Rcpp:::.showBuildFailureDiagnostics()
      Rcpp:::.callBuildCompleteHook(succeeded, output)
      setwd(cwd)
      Rcpp:::.restoreEnvironment(envRestore)
    })
    if (file.exists(context$previousDynlibPath)) 
    {
      try(silent = TRUE, dyn.unload(context$previousDynlibPath))
      file.remove(context$previousDynlibPath)
    }
    
    
    # r <- file.path(R.home("bin"), "R")
    # lib <- context$dynlibFilename
    # deps <- context$cppDependencySourcePaths
    # src <- context$cppSourceFilename
    
    # ==========================================================================
    # Double ensure the build directory only contains .cpp or .R files.
    # ==========================================================================
    filesInBuildDir = list.files(dirname(context$dynlibPath), full.names = TRUE)
    files2delete = filesInBuildDir[!(
      tools::file_ext(filesInBuildDir) %in% c("cpp", "R"))]
    file.remove(files2delete)
    
    
    # ==========================================================================
    # Fill in include paths, link paths and link libraries.
    # ==========================================================================
    includePaths = c(includePaths, getRcppRelatedIncludePath(context$depends))
    includePaths = includePaths[includePaths != '']
    includePaths = paste0('-I"', includePaths, '"')
    includePaths = paste0(includePaths[includePaths != ''], collapse = '  ')
    if (length(dllLinkFilePaths) != 0)
    {
      Lpaths = dirname(dllLinkFilePaths)
      lnames = tools::file_path_sans_ext(basename(dllLinkFilePaths))
      Lpaths = paste0(paste0('-L"', Lpaths, '"'), collapse = '  ')
      lnames = paste0(paste0('-l', lnames), collapse = '  ')
    } else 
    {
      Lpaths = ''
      lnames = ''
    }
    
    
    flags = paste0(flags, ' ', optFlag)
    if (sanitize) flags = paste0(flags, '  ', sanitizerFlags)
    
    
    # return(list(compilerPath, flags, includePaths, Lpaths))
    
    cmd = paste0(compilerPath, '  ', 
                 flags, '  ', 
                 includePaths, '  ',
                 Lpaths, '  ',
                 lnames)
    cmd = paste0(cmd, '  ', srcFile, '  -o  ', context$dynlibPath)
    if (sanitize) cmd = paste0("unset LD_PRELOAD;  ", cmd)
    if (verbose) cat(cmd, "\n")
    result = system(cmd)
    
    
    if (result != 0) stop(paste0(
      "Error occured building shared library. Build directory = ", 
      context$buildDirectory))
    if (!file.exists(context$dynlibFilename)) stop(paste0(
      'Shared library is not built. Build directory = ', 
      context$buildDirectory ))
    
    
    # args <- c("CMD", "SHLIB", 
    #           if (windowsDebugDLL) "-d", 
    #           if (rebuild) "--preclean", 
    #           if (dryRun) "--dry-run", 
    #           "-o", shQuote(lib), 
    #           if (length(deps)) paste(shQuote(deps), collapse = " "), 
    #           shQuote(src))
    
    
    # if (verbose) 
    #   cat(paste(c(r, args), collapse = " "), "\n")
    # so <- if (verbose) "" else TRUE
    # 
    # 
    # result <- suppressWarnings(system2(
    #   compilerPath, args, stdout = so, stderr = so))
    
    
    # if (!verbose) {
    #   output <- result
    #   attributes(output) <- NULL
    #   status <- attr(result, "status")
    #   if (!is.null(status)) {
    #     cat(result, sep = "\n")
    #     succeeded <- FALSE
    #     stop("Error ", status, " occurred building shared library.")
    #   }
    #   else if (!file.exists(context$dynlibFilename)) {
    #     cat(result, sep = "\n")
    #     succeeded <- FALSE
    #     stop("Error occurred building shared library.")
    #   }
    #   else {
    #     succeeded <- TRUE
    #   }
    # }
    # else if (!identical(as.character(result), "0")) {
    #   succeeded <- FALSE
    #   stop("Error ", result, " occurred building shared library.")
    # }
    # else {
    #   succeeded <- TRUE
    # }
  }
  else 
  {
    cwd <- getwd()
    on.exit({setwd(cwd)})
    if (verbose) cat("\nNo rebuild required (use rebuild = TRUE to ", 
                     "force a rebuild)\n\n", sep = "")
  }
  
  
  # if (dryRun) return(invisible(NULL))
  
  
  if (length(context$exportedFunctions) > 0 || length(context$modules) > 0) 
  {
    exports <- c(context$exportedFunctions, context$modules)
    removeObjs <- exports[exports %in% ls(envir = env, all.names = T)]
    remove(list = removeObjs, envir = env)
    scriptPath <- file.path(context$buildDirectory, context$rSourceFilename)
    source(scriptPath, local = env)
  }
  else if (getOption("rcpp.warnNoExports", default = TRUE)) 
  {
    warning("No Rcpp::export attributes or RCPP_MODULE declarations ", 
            "found in source")
  }
  
  
  # if (embeddedR && (length(context$embeddedR) > 0)) {
  #   srcConn <- textConnection(context$embeddedR)
  #   setwd(rWorkingDir)
  #   source(file = srcConn, local = env, echo = echo)
  # }
  # if (cleanupCacheDir) 
  #   Rcpp:::cleanupSourceCppCache(cacheDir, context$cppSourcePath, 
  #                         context$buildDirectory)
  
  
  # invisible(list(functions = context$exportedFunctions, 
  #                modules = context$modules, 
  #                cppSourcePath = context$cppSourcePath, 
  #                buildDirectory = context$buildDirectory))
  invisible(context)
}











