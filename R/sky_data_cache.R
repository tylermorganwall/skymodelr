#' List cached sky data files
#'
#' @description
#' Lists files cached by [download_sky_data()] under
#' `tools::R_user_dir("skymodelr", "data")`.
#'
#' @returns A data frame with columns `file`, `path`, `size`, and `modified`.
#' @export
#'
#' @examples
#' list_sky_data()
list_sky_data = function() {
  cache_dir = prague_coef_cache_dir()
  empty = data.frame(
    file = character(),
    path = character(),
    size = numeric(),
    modified = as.POSIXct(character()),
    stringsAsFactors = FALSE
  )

  if (!dir.exists(cache_dir)) {
    return(empty)
  }

  paths = list.files(cache_dir, full.names = TRUE, no.. = TRUE)
  if (!length(paths)) {
    return(empty)
  }

  info = file.info(paths)
  keep = !is.na(info$isdir) & !info$isdir
  if (!any(keep)) {
    return(empty)
  }

  paths = paths[keep]
  info = info[keep, , drop = FALSE]

  data.frame(
    file = basename(paths),
    path = normalizePath(paths, winslash = "/", mustWork = FALSE),
    size = info$size,
    modified = info$mtime,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' Clear cached sky data files
#'
#' @description
#' Removes files cached by [download_sky_data()] under
#' `tools::R_user_dir("skymodelr", "data")`.
#'
#' @param files Default `NULL`. Character vector of cached file basenames to remove.
#' @param ask Default `interactive()`. Whether to ask for confirmation before deleting files.
#'
#' @returns Invisibly, the paths successfully removed.
#' @export
#'
#' @examplesIf interactive()
#' clear_sky_data()
clear_sky_data = function(files = NULL, ask = interactive()) {
  cache_dir = prague_coef_cache_dir()
  if (!dir.exists(cache_dir)) {
    return(invisible(character()))
  }

  cached = list_sky_data()
  if (is.null(files)) {
    paths = cached$path
  } else {
    requested = unique(basename(files))
    requested = requested[nzchar(requested) & !requested %in% c(".", "..")]
    paths = cached$path[cached$file %in% requested]
  }

  if (!length(paths)) {
    message("No matching cached sky data files found.")
    return(invisible(character()))
  }

  if (isTRUE(ask)) {
    prompt = sprintf(
      "Remove %d cached sky data file%s?",
      length(paths),
      if (length(paths) == 1) "" else "s"
    )
    if (!isTRUE(utils::askYesNo(prompt, default = FALSE))) {
      return(invisible(character()))
    }
  }

  unlink(paths)
  removed = paths[!file.exists(paths)]
  invisible(removed)
}
