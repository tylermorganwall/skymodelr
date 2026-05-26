clamp_value = function(x, min_value = 0, max_value = 1) {
  pmin(pmax(x, min_value), max_value)
}

get_moon_texture = function() {
  texture_file = system.file(
    "textures",
    "lroc_color_poles_1k.jpg",
    package = "skymodelr"
  )
  if (!nzchar(texture_file)) {
    stop(
      "Moon texture file 'lroc_color_poles_1k.jpg' was not found in the installed package.",
      call. = FALSE
    )
  }
  texture_file
}
