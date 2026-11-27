testthat::test_that("constant env is resolution invariant", {
  make_env = function(resolution) {
    env = array(1, dim = c(resolution, resolution * 2, 3))
    env_path = tempfile(fileext = ".exr")
    rayimage::ray_write_image(env, env_path, clamp = FALSE)
    env_path
  }

  render_plane = function(env_path) {
    scene = rayrender::xz_rect(
      x = 0,
      xwidth = 2,
      z = 0,
      zwidth = 2,
      y = 0,
      material = rayrender::diffuse(color = "#ffffff")
    )
    set.seed(1)
    render = rayrender::render_scene(
      scene,
      width = 32,
      height = 32,
      fov = 0,
      ortho_dimensions = c(2, 2),
      lookfrom = c(0, 1, 0),
      lookat = c(0, 0, 0),
      camera_up = c(0, 0, 1),
      samples = 64,
      ambient_light = FALSE,
      backgroundhigh = "#000000",
      backgroundlow = "#000000",
      environment_light = env_path,
      intensity_env = 1,
      tonemap = "raw",
      bloom = FALSE,
      preview = FALSE,
      interactive = FALSE,
      plot_scene = FALSE,
      progress = FALSE
    )
    mean(render[,, 1:3])
  }

  env_low = make_env(32)
  env_high = make_env(256)

  mean_low = render_plane(env_low)
  mean_high = render_plane(env_high)

  testthat::expect_lt(abs(mean_low - mean_high), 0.05)
})
