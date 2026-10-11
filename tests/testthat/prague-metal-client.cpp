// [[Rcpp::depends(skymodelr)]]
#include <Rcpp.h>
#include <skymodelr/prague.hpp>
#include <cstring>
// [[Rcpp::export]]
Rcpp::List prague_metal_contract(std::string filename) {
  std::shared_ptr<const skymodelr::PragueMetalSnapshot> snapshot;
  bool shared=false;
  {
    skymodelr::PragueSkyModel model(filename,50,false,false,0);
    snapshot=model.metalSnapshot();
    shared=snapshot==model.metalSnapshot();
  } // Export must survive CPU model destruction.
  size_t count; const float *data=snapshot->data(count);
  auto m=snapshot->layout();
  bool bounds=true;
  for(auto axis:{m.visibility,m.albedo,m.altitude,m.elevation,m.gamma,m.alpha,m.emphasis,
                 m.trans_altitude,m.trans_visibility,m.sun}) {
    bounds=bounds && axis.count>0 && size_t(axis.offset)+axis.count<=count;
  }
  return Rcpp::List::create(Rcpp::Named("shared")=shared,Rcpp::Named("bounds")=bounds,
    Rcpp::Named("finite")=std::isfinite(data[m.radiance]),
    Rcpp::Named("source")=std::strstr(snapshot->source(),"prague_spectrum_batch")!=nullptr,
    Rcpp::Named("parameters")=std::strstr(snapshot->source(),"prague_parameters_batch")!=nullptr,
    Rcpp::Named("bytes")=double(count*sizeof(float)));
}
// [[Rcpp::export]]
Rcpp::NumericVector prague_short_transmission(std::string filename) {
  skymodelr::PragueSkyModel model(filename,50,false,false,0);
  auto p=model.computeParameters({0,0,1000},{1,0,.3},.5,0,50,.3);
  // The fitted value need not be one, but sub-millimetre changes must be smooth.
  Rcpp::NumericVector output;
  for(double distance:{1e-6,1.001e-6,1.002e-6,1.003e-6,1.004e-6})
    output.push_back(model.transmittance(p,560,distance));
  return output;
}
