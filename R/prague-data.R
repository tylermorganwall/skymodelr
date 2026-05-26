prague_coef_cache_dir = function() {
  tools::R_user_dir("skymodelr", "data")
}

prague_coef_info = function(
  sea_level = TRUE,
  wide_spectrum = FALSE,
  altitude = NULL
) {
  if (!is.null(altitude)) {
    sea_level = all(altitude == 0)
  }
  if (wide_spectrum && !sea_level) {
    stop(
      "`wide_spectrum = TRUE` is only valid when `altitude == 0`.",
      call. = FALSE
    )
  }

  if (sea_level && wide_spectrum) {
    list(
      name = "PragueSkyModelDatasetGroundInfra.dat",
      url = "https://skydata.tylermw.com/PragueSkyModelDatasetGroundInfra.dat",
      md5 = "33b97729ce8cd7fbdfc60317ff2805e7",
      size = "574MB",
      sea_level = TRUE,
      wide_spectrum = TRUE
    )
  } else if (sea_level) {
    list(
      name = "SkyModelDatasetGround.dat",
      url = "https://skydata.tylermw.com/SkyModelDatasetGround.dat",
      md5 = "953864f62c9434269b7f909cdff6c7cd",
      size = "107MB",
      sea_level = TRUE,
      wide_spectrum = FALSE
    )
  } else {
    list(
      name = "SkyModelDataset.dat",
      url = "https://skydata.tylermw.com/SkyModelDataset.dat",
      md5 = "dcfa2bd165cf803c1791c469e05d134c",
      size = "2.4GB",
      sea_level = FALSE,
      wide_spectrum = FALSE
    )
  }
}

prague_download_call = function(info) {
  sprintf(
    "download_sky_data(sea_level = %s, wide_spectrum = %s)",
    if (info$sea_level) "TRUE" else "FALSE",
    if (info$wide_spectrum) "TRUE" else "FALSE"
  )
}

stop_missing_prague_coef = function(info, model = "Prague sky model") {
  stop(
    sprintf(
      "The %s coefficient file '%s' is not installed. Run %s explicitly to download it.",
      model,
      info$name,
      prague_download_call(info)
    ),
    call. = FALSE
  )
}

resolve_prague_coef_file = function(
  altitude = 0,
  wide_spectrum = FALSE,
  allow_download = interactive(),
  model = "Prague sky model"
) {
  info = prague_coef_info(altitude = altitude, wide_spectrum = wide_spectrum)
  coef_file = file.path(prague_coef_cache_dir(), info$name)
  if (file.exists(coef_file)) {
    return(coef_file)
  }

  if (!allow_download) {
    stop_missing_prague_coef(info, model = model)
  }

  answer = readline(
    prompt = sprintf(
      "Coefficient file '%s' for the %s is not installed. This is a large file (%s). Download now? [y/n] ",
      info$name,
      model,
      info$size
    )
  )
  answer = tolower(trimws(answer))

  if (answer %in% c("y", "yes")) {
    return(download_sky_data(
      sea_level = info$sea_level,
      wide_spectrum = info$wide_spectrum
    ))
  }
  if (answer %in% c("n", "no")) {
    stop_missing_prague_coef(info, model = model)
  }
  stop("Input not recognized.", call. = FALSE)
}
