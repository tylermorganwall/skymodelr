// [[Rcpp::depends(RcppThread)]]
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <RcppThread.h>
#include <algorithm>
#include <mutex>
#include <cmath>

inline double mag_to_luminance(double Vmagnitude,
                               double zeroPoint = 1.0) {
  return zeroPoint * std::pow(10.0, (-14.18 - Vmagnitude)/2.5);
}


static inline int wrap_col(int x, int W) {
  int r = x % W; return r < 0 ? r + W : r;
}

static inline double jd_to_gmst(double jd)
{
    double d  = jd - 2451545.0;                  // days since J2000
    double t  = d / 36525.0;                     // Julian centuries
    double gmst =
        280.46061837 + 360.98564736629 * d +
        0.000387933 * t * t - t * t * t / 38710000.0; // IAU 2006
    gmst = fmod(gmst, 360.0);
    if (gmst < 0) gmst += 360.0;
    return gmst * M_PI / 180.0;                  // rad
}

// [[Rcpp::export]]
Rcpp::NumericVector make_starfield_rcpp(Rcpp::DataFrame stars,
                         unsigned int    resolution = 2048,
                         double          zero_point = 1.0,
            						 double          lon_deg   = 0.0,
            						 double          lat_deg   = 0.0,
            						 double          jd        = 2451545.0, // Julian Date of the frame (UTC noon 1 Jan 2000 = 2451545.0)
                         double          turbidity = 3.0,
            						 double          ozone_du  = 300.0,
            						 double          altitude  = 0.0,
            						 double          star_width = 1.0,
            						 bool            use_rgb   = true,
            						 bool            atmosphere_effects = true,
            						 bool            upper_hemisphere_only = true,
            						 unsigned int    numbercores  = 1)
{
  if (resolution == 0) Rcpp::stop("resolution must be greater than or equal to 1");
  if(!upper_hemisphere_only & atmosphere_effects) {
    atmosphere_effects = false;
    Rcpp::message(Rcpp::CharacterVector("When rendering full sphere (`upper_hemisphere_only = false`), atmosphere effects are not rendered (setting `atmosphere_effects = false`)"));
  }
  const int nTheta = resolution;        // rows (latitude)
  const int nPhi   = 2 * resolution;    // cols (longitude)
  const int nPix   = nTheta * nPhi;
  const double dtheta = M_PI / double(nTheta);
  const double dphi   = 2.0 * M_PI / double(nPhi);

  auto omega_row = [&](int t) {
    double theta_c = (double(t) + 0.5) * dtheta;
    double s = std::sin(theta_c);
    if (s < 1e-8) s = 1e-8;
    return dtheta * dphi * s;
  };

  std::vector<double> img(3 * nPix, 0.0);

  // pull columns out of the data frame
  Rcpp::NumericVector ra   = stars["ra_rad"];
  Rcpp::NumericVector dec  = stars["dec_rad"];
  Rcpp::NumericVector mag  = stars["v_mag"];
  Rcpp::NumericVector r_vec  = stars["r"];
  Rcpp::NumericVector g_vec  = stars["g"];
  Rcpp::NumericVector b_vec  = stars["b"];

  const std::size_t   N    = ra.size();

  double psf_fwhm_pix = star_width;
  double psf_trunc_sigma = 3.0;  // truncate at 3 sigma

	double gmst_rad = jd_to_gmst(jd);
	double lst_rad  = gmst_rad + lon_deg * M_PI / 180.0;
	double sin_lat  = std::sin(lat_deg * M_PI / 180.0);
	double cos_lat  = std::cos(lat_deg * M_PI / 180.0);
	std::vector<std::mutex> row_locks(nTheta);

  // parallelise over stars
  RcppThread::parallelFor(0u, (unsigned)N, [&](unsigned i)
    {
  		double ha = lst_rad - ra[i];

  		double sin_alt = std::sin(dec[i]) * sin_lat +
                   std::cos(dec[i]) * cos_lat * std::cos(ha);
  		double alt_rad = std::asin( std::clamp(sin_alt, -1.0, 1.0) );
  		double cos_alt = std::cos(alt_rad);

  		if(cos_alt < 1e-6) cos_alt = 1e-6;

  		double sin_az  =  std::cos(dec[i]) * std::sin(ha) / cos_alt;
  		double cos_az  = (std::sin(dec[i]) * cos_lat -
  						std::cos(dec[i]) * std::cos(ha) * sin_lat) / cos_alt;

  		double x_h =  -cos_alt *  cos_az;
  		double y_h =  sin_alt;
  		double z_h =  -cos_alt *  sin_az;

  		//map to lat/long image
  		double theta = std::acos(std::clamp(y_h, -1.0, 1.0));
  		double phi   = std::atan2(-z_h, x_h);
  		if (phi < 0.0) phi += 2.0 * M_PI;

  		double phi_img = phi + M_PI;
  		if (phi_img >= 2.0 * M_PI) phi_img -= 2.0 * M_PI;

  		double r = 1.0, g = 1.0, b = 1.0;
  		if(use_rgb) {
  			r = r_vec[i];
  			g = g_vec[i];
  			b = b_vec[i];
  		}

  		double T_r = 1.0, T_g = 1.0, T_b = 1.0;
  		constexpr double tauR[9] = {0.00835,0.00673,0.00595,0.01750,0.01350,
  									0.01100,0.02690,0.02050,0.01880};
  		constexpr double tauM0[9]= {0.00220,0.00190,0.00175,0.00300,0.00265,
  									0.00245,0.00360,0.00320,0.00310};
  		constexpr double tauO[9] = {0.0030,0.0048,0.0060,0.0014,0.0010,
  									0.0009,0.0006,0.0007,0.0008};
  		if (theta < M_PI_2) {           // above horizon only
  			// ---------- airmass m(theta) -----------------------------------------
  			double secz = 1.0 / std::max(1e-6, std::cos(theta));
  			double m    = secz < 12.0
  						? secz - 0.001816*pow(secz-1,1)
  						 	   - 0.002875*pow(secz-1,2)
  							   - 0.000808*pow(secz-1,3)
  						: secz;

  			// density scale with observer altitude
  			double rho = std::exp(-altitude / 8000.0);          // rho/rho0  (simple barometric)

  			// coefficient scalers
  			double kR = rho;                                    // Rayleigh air density
  			double kM = rho * (turbidity - 1.0);                // Mie 0 at T=1
  			double kO = ozone_du / 300.0;                       // 1.0 at 300 DU

  			double T_band[9];
  			for (int c=0; c<9; ++c) {
  				T_band[c] = std::exp( -(kR*tauR[c] + kM*tauM0[c] + kO*tauO[c]) * m );
  			}
  			T_r = (T_band[0] + T_band[1] + T_band[2]) / 3.0;
  			T_g = (T_band[3] + T_band[4] + T_band[5]) / 3.0;
  			T_b = (T_band[6] + T_band[7] + T_band[8]) / 3.0;
  		} else {
			  T_r = T_g = T_b = 0.0;      // star is below horizon
  		}

  		if(!upper_hemisphere_only) {
  		  T_r = T_g = T_b = 1.0;
  		}

  		if (T_r==0.0 && T_g==0.0 && T_b==0.0) return;

  		// ----- continuous image coords -----
  		double u = (phi_img / (2.0 * M_PI)) * nPhi;    // column in [0, nPhi)
  		double v = (theta   / M_PI)         * nTheta;  // row    in [0, nTheta)

  		// avoid rare u == nPhi or v == nTheta due to floating point
  		if (u >= nPhi)   u = std::nextafter((double)nPhi, 0.0);
  		if (v >= nTheta) v = std::nextafter((double)nTheta, 0.0);

  		// subpixel center and integer anchors
  		int    p0 = (int)std::floor(u);
  		int    t0 = (int)std::floor(v);

  		// PSF parameters
  		double sigma = psf_fwhm_pix / 2.355;                   // FWHM→σ
  		int    rad   = (int)std::ceil(psf_trunc_sigma * sigma);
  		double inv_2s2 = 1.0 / (2.0 * sigma * sigma);

  		// Pre-compute 1D weights using distance from the pixel center (p + 0.5, t + 0.5)
  		std::vector<double> wx(2*rad + 1), wy(2*rad + 1);

  		for (int dx = -rad; dx <= rad; ++dx) {
 		  double px_center = (p0 + dx) + 0.5;
 		  double ddx = px_center - u;                   // subpixel offset in x
 		  double w   = std::exp(-ddx*ddx * inv_2s2);
		  wx[dx + rad] = w;
  		}
  		for (int dy = -rad; dy <= rad; ++dy) {
 		  int ty = std::clamp(t0 + dy, 0, nTheta - 1);
 		  double py_center = ty + 0.5;
 		  double ddy = py_center - v;                   // subpixel offset in y
 		  double w   = std::exp(-ddy*ddy * inv_2s2);
		  wy[dy + rad] = w;
  		}

  		// Previous norm was image-space only; now include solid angle for env-map consistency
  		double norm = 0.0;
  		for (int dy = -rad; dy <= rad; ++dy) {
  		  int t = std::clamp(t0 + dy, 0, nTheta - 1);
  		  double wyj = wy[dy + rad];
  		  double omega = omega_row(t);
  		  for (int dx = -rad; dx <= rad; ++dx) {
  		    double w = wyj * wx[dx + rad];
  		    norm += w * omega;
  		  }
  		}
  		if (norm <= 0.0) return;

  		// Star radiance in channels
  		double  L  = mag_to_luminance(mag[i], double(zero_point));
  		double R0 = L * r * T_r / norm;
  		double G0 = L * g * T_g / norm;
  		double B0 = L * b * T_b / norm;

  		// ----- accumulate with row locks (one lock per affected row) -----
  		for (int dy = -rad; dy <= rad; ++dy) {
  		  int t = std::clamp(t0 + dy, 0, nTheta - 1);
  		  double wyj = wy[dy + rad];

  		  // lock the row once; update all columns for this row
  		  {
  		    std::scoped_lock lock(row_locks[t]);
  		    std::size_t base = 3ULL * ( (std::size_t)t * (std::size_t)nPhi );

  		    for (int dx = -rad; dx <= rad; ++dx) {
  		      int p = wrap_col(p0 + dx, nPhi);
  		      double w = wyj * wx[dx + rad];   // separable Gaussian

  		      std::size_t idx = base + 3ULL * (std::size_t)p;
  		      img[idx    ] += double(w * R0);
  		      img[idx + 1] += double(w * G0);
  		      img[idx + 2] += double(w * B0);
  		    }
  		  }
  		}
    }, numbercores);

    Rcpp::NumericVector result(nPix * 3);
    result.attr("dim") = Rcpp::IntegerVector::create(nTheta, nPhi, 3);

    for (int t = 0; t < nTheta; ++t) {
      for (int p = 0; p < nPhi; ++p) {
        const int pix_index = t * nPhi + p;
        for (int c = 0; c < 3; ++c) {
          const int arr_index = t + nTheta * (p + nPhi * c);
          result[arr_index] = img[3 * pix_index + c];
        }
      }
    }

    return result;
}
