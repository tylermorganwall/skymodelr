// [[Rcpp::depends(RcppThread)]]
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <RcppThread.h>
#include <OpenEXR/ImfRgbaFile.h>
#include <OpenEXR/ImfRgba.h>
#include <Imath/ImathVec.h>

inline float mag_to_luminance(float Vmagnitude,
                              float zeroPoint = 1.0f)
{
    // Physical scaling constant not needed for relative EXR; zeroPoint
    // lets you tune exposure in post.
    return zeroPoint * std::pow(10.0f, -0.4f * Vmagnitude);
}

using Imf::Rgba;
using Imf::RgbaOutputFile;
using Imath::Vec3;

static inline int clamp_int(int x, int lo, int hi)
{
	return x < lo ? lo : (x > hi ? hi : x);
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

//' Render a star-field equirectangular EXR
//'
//' @param outfile Output `.exr` path.
//' @param stars   Data frame produced by `read_bsc5()` (must contain
//'                `ra_rad`, `dec_rad`, `v_mag`).
//' @param resolution Horizontal *half-resolution*; final size is
//'        2 × resolution × resolution.
//' @param zero_point Default `1`. Exposure scaling (see docs).
//' @param ncores Number of CPU threads.
// [[Rcpp::export]]
void make_starfield_rcpp(std::string     outfile,
                         Rcpp::DataFrame stars,
                         unsigned int    resolution = 2048,
                         double          zero_point = 1.0,
						 double          lon_deg   = 0.0,
						 double          lat_deg   = 0.0,
						 double          jd        = 2451545.0,   // Julian Date of the frame (UTC noon 1 Jan 2000 = 2451545.0)
                         double          turbidity = 3.0,
						 double          ozone_du  = 300.0,
						 double          altitude  = 0.0,
						 bool            use_rgb   = true,
						 unsigned int    ncores     = 1)
{
    if (resolution == 0) Rcpp::stop("resolution must be ≥1");

    const int nTheta = resolution;        // rows (latitude)
    const int nPhi   = 2 * resolution;    // cols (longitude)
    const int nPix   = nTheta * nPhi;

    std::vector<float> img(3 * nPix, 0.0f);

    // pull columns out of the data frame
    Rcpp::NumericVector ra   = stars["ra_rad"];
    Rcpp::NumericVector dec  = stars["dec_rad"];
    Rcpp::NumericVector mag  = stars["v_mag"];
	Rcpp::NumericVector r_vec  = stars["r"];
    Rcpp::NumericVector g_vec  = stars["g"];
    Rcpp::NumericVector b_vec  = stars["b"];

    const std::size_t   N    = ra.size();

	double gmst_rad = jd_to_gmst(jd);
	double lst_rad  = gmst_rad + lon_deg * M_PI / 180.0;
	double sin_lat  = std::sin(lat_deg * M_PI / 180.0);
	double cos_lat  = std::cos(lat_deg * M_PI / 180.0);

    // parallelise over stars – cheap, embarrassingly parallel
    RcppThread::parallelFor(0u, (unsigned)N, [&](unsigned i)
    {
        // equirectangular: phi = RA in [0,2π), theta = π/2 - Dec  in [0,π]
        // double phi   = std::fmod(ra[i] + 2*M_PI, 2*M_PI);
        // double theta = M_PI_2 - dec[i];
		// 1. hour angle (east +)
		double ha = lst_rad - ra[i];

		// 2. equatorial unit vector (right-handed X = south)
		// double cos_dec = std::cos(dec[i]);
		// double x_eq =  cos_dec * std::cos(ha);
		// double y_eq =  cos_dec * std::sin(ha);
		// double z_eq =  std::sin(dec[i]);

		double sin_alt = std::sin(dec[i]) * sin_lat +
                 std::cos(dec[i]) * cos_lat * std::cos(ha);
		double alt_rad     = std::asin( std::clamp(sin_alt, -1.0, 1.0) );
		double cos_alt = std::cos(alt_rad);

		if(cos_alt < 1e-6) cos_alt = 1e-6;

		double sin_az  =  std::cos(dec[i]) * std::sin(ha) / cos_alt;
		double cos_az  = (std::sin(dec[i]) * cos_lat -
						std::cos(dec[i]) * std::cos(ha) * sin_lat) / cos_alt;

		/*  coordinate frame used by makesky:
			+X  = South
			+Y  = Up (zenith)
			+Z  = West
			azimuth counted positive to the East from South            */

		double x_h =  cos_alt *  cos_az;   // South (+)  /  North (−)
		double y_h =  sin_alt;             // Up
		double z_h =  cos_alt *  sin_az;   // West (+)   /  East (−)

		// 3. rotate to horizontal frame
		// double x_h =  x_eq * sin_lat - z_eq * cos_lat;   // points south
		// double y_h =  y_eq;                              // points up
		// double z_h =  x_eq * cos_lat + z_eq * sin_lat;   // points west

		// 4. map to lat/long image
		double theta = std::acos(std::clamp(y_h, -1.0, 1.0));      // 0..π
		double phi   = std::atan2(-z_h, x_h);                      // -π..π
		if (phi < 0.0) phi += 2.0 * M_PI;                          // 0..2π

        // ignore stars that fall outside north-hemisphere canvas if user
        // feeds non-northern declinations and keeps theta∈[0,π]
        // if (theta < 0.0 || theta > M_PI) return;

        int p = int(phi / (2*M_PI) * nPhi);        // column
        int t = int(theta / M_PI * nTheta);        // row

        p = clamp_int(p, 0, nPhi-1);
        t = clamp_int(t, 0, nTheta-1);

        // pixel index
        std::size_t idx = 3 * (t*nPhi + p);

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
			double rho = std::exp(-altitude / 8000.0);          // ρ/ρ0  (simple barometric)

			// coefficient scalers
			double kR = rho;                                    // Rayleigh ∝ air density
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

        float L = mag_to_luminance(mag[i], float(zero_point));
        img[idx    ] += L * r * T_r;
        img[idx + 1] += L * g * T_g;
        img[idx + 2] += L * b * T_b;
    }, ncores);

    // convert to OpenEXR RGBA array
    std::vector<Rgba> px(nPix);
    for (int k = 0; k < nPix; ++k)
        px[k] = Rgba(img[3*k], img[3*k+1], img[3*k+2]);

    RgbaOutputFile out(outfile.c_str(), nPhi, nTheta, Imf::WRITE_RGB);
    out.setFrameBuffer(px.data(), 1, nPhi);
    out.writePixels(nTheta);
}
