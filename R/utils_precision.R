warn_precision_loss = function(path) {
  ext = tools::file_ext(path)
  if (nzchar(ext) && tolower(ext) != "exr") {
    message("Note: writing to non-EXR formats may introduce precision loss and unphysically clamp luminosity values.")
  }
}
