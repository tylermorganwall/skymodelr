R"PRAGUE(
// Reconstruction derived from PragueSkyModel (Copyright 2022 Charles
// University, Apache-2.0); distributed within the GPL-3 skymodelr package.
// Same compressed spectral fit as the CPU evaluator. No RGB/angular LUT.
namespace skymodelr_gpu {
struct Axis { uint offset, count; };
struct Layout {
  Axis visibility, albedo, altitude, elevation, gamma, alpha, emphasis;
  Axis trans_altitude, trans_visibility, sun;
  uint radiance, trans_u, trans_v, rank, sun_offset, sun_stride;
  uint zenith_offset, zenith_stride, emphasis_offset, config_stride;
  uint trans_rank, a_dim, d_dim, channels;
  float channel_start, channel_width;
};
struct Parameters { float theta, gamma, shadow, zero, elevation, altitude, visibility, albedo; };
struct Interpolation { uint index; float factor; };
constant float radius = 6378000.f;
constant float pi = 3.14159265358979323846f;
// Preserve metre/submetre heights: rationalize |p+(0,0,R)|-R instead of
// subtracting two approximately 6,378,000 m floats. atan2 preserves small angles.
float height(float3 p) {
  float n = dot(p,p) + 2*radius*p.z;
  return n / (length(p+float3(0,0,radius)) + radius);
}
float angle(float3 a, float3 b) { return atan2(length(cross(a,b)),dot(a,b)); }
Parameters parameters(float3 p, float3 w, float elevation, float azimuth, float visibility, float albedo) {
  float3 up = normalize(p+float3(0,0,radius));
  w = normalize(w);
  float h = max(0.f,height(p)+50.f);
  float correction = sqrt(h*(2*radius+h))/(radius+h);
  // Algebraic cancellation of the two shifted Earth-space points in the CPU
  // formula, before rounding, preserves the one-metre viewing direction.
  float3 corrected = normalize(w + up*correction);
  float3 sun(cos(azimuth)*cos(elevation),sin(azimuth)*cos(elevation),sin(elevation));
  float3 shadow(-sin(elevation)*cos(azimuth),-sin(elevation)*sin(azimuth),cos(elevation));
  return {angle(w,up), angle(w,sun), angle(corrected,shadow), angle(corrected,up),
          asin(clamp(dot(up,sun),-1.f,1.f)), h, visibility, albedo};
}
Interpolation lookup(float value, Axis axis, device const float *data) {
  float v=clamp(value,data[axis.offset],data[axis.offset+axis.count-1]);
  uint lo=1,hi=axis.count;
  while(lo<hi) { uint mid=(lo+hi)/2; if(data[axis.offset+mid]<=v) lo=mid+1; else hi=mid; }
  uint i=lo-1;
  return {i,lo==axis.count ? 0.f : (v-data[axis.offset+i])/(data[axis.offset+lo]-data[axis.offset+i])};
}
float linear(device const float *data,uint offset,float factor) {
  if(factor==0) return data[offset]; // Last knot need not have a successor.
  return (data[offset+1]-data[offset])*factor+data[offset];
}
struct SkyQuery { Interpolation gamma, alpha, emphasis, config[4]; };
// Prepare wavelength-independent coordinates once for a spectral integration.
SkyQuery prepare_sky(Parameters p,device const Layout &m,device const float *data) {
  return {lookup(p.gamma,m.gamma,data),lookup(p.elevation<0 ? p.shadow : p.zero,m.alpha,data),
    lookup(p.zero,m.emphasis,data),{lookup(p.visibility,m.visibility,data),lookup(p.albedo,m.albedo,data),
    lookup(p.altitude,m.altitude,data),lookup(p.elevation*180/pi,m.elevation,data)}};
}
float sky_prepared(SkyQuery query,float wavelength,device const Layout &m,device const float *data) {
  if(wavelength<m.channel_start || wavelength>=m.channel_start+m.channels*m.channel_width) return 0;
  uint channel=uint((wavelength-m.channel_start)/m.channel_width);
  auto g=query.gamma,a=query.alpha,e=query.emphasis;
  auto config=query.config;
  uint counts[4]={m.visibility.count,m.albedo.count,m.altitude.count,m.elevation.count};
  // Weighted corner accumulation avoids a dynamically indexed 16-value local
  // array. It retains the CPU's zero-weight threshold and per-corner clamping;
  // summation order differs, so numerical agreement is tolerance-based.
  float result=0;
  for(uint i=0;i<16;++i) {
    uint index=0; float weight=1;
    for(uint j=0;j<4;++j) {
      uint high=(i>>(3-j))&1;
      float f=config[j].factor<1e-6f ? 0.f : config[j].factor;
      weight*=high ? f : 1-f;
      index=index*counts[j]+min(config[j].index+high,counts[j]-1);
    }
    float value=0;
    if(weight>0) {
      uint base=m.radiance+(index*m.channels+channel)*m.config_stride;
      for(uint r=0;r<m.rank;++r)
        value+=linear(data,base+m.sun_offset+r*m.sun_stride+g.index,g.factor)*
               linear(data,base+m.zenith_offset+r*m.zenith_stride+a.index,a.factor);
      value=max(0.f,value*linear(data,base+m.emphasis_offset+e.index,e.factor));
    }
    result+=weight*value;
  }
  return result;
}
float sky(Parameters p,float wavelength,device const Layout &m,device const float *data) {
  return sky_prepared(prepare_sky(p,m,data),wavelength,m,data);
}
// Nearest forward intersection. c is factored to preserve altitude close to
// Earth; the near root uses c/(-b+sqrt(discriminant)) to avoid cancellation.
float circle(float cosine,float altitude,float shell) {
  float r=radius+altitude, b=r*cosine;
  float c=(altitude-shell)*(2*radius+altitude+shell);
  float discriminant=b*b-c;
  if(discriminant<=0) return -1;
  float root=sqrt(discriminant), far=-b+root;
  float near=c/far;
  return near>0 && far>0 ? min(near,far) : max(near,far);
}
Interpolation trans_lookup(float value,uint n,uint power) {
  uint i=min(uint(value*n),n-1); float f=0;
  if(i<n-1) {
    float a=float(i)/n,b=float(i+1)/n;
    float va=pow(a,float(power)),vb=pow(b,float(power));
    f=clamp((pow(value,float(power))-va)/(vb-va),0.f,1.f);
  }
  return {i,f};
}
struct TransParameters { Interpolation altitude, distance; };
TransParameters trans_parameters(Parameters p,float distance,device const Layout &m) {
  float cosine=cos(p.theta),sine=sin(p.theta),t,shell=0;
  if(p.altitude<.3f) {
    if(p.theta<=pi/2) {shell=100000;t=circle(cosine,p.altitude,shell);} else t=0;
  } else {
    t=circle(cosine,p.altitude,0);
    if(t<0) {shell=100000;t=circle(cosine,p.altitude,shell);}
  }
  bool boundary=t<=distance;
  t=min(max(0.f,t),distance);
  float x=sine*t,z=p.altitude+cosine*t;
  float h=boundary ? shell : clamp(height(float3(x,0,z)),0.f,100000.f);
  float a=pow(h/100000.f,1.f/3.f);
  // acos(y/length) rounds to zero for short paths. atan2 is equivalent here.
  float d=min(1.f,sqrt(sqrt(atan2(abs(x),radius+z)*radius/1571524.413613f)));
  return {trans_lookup(a,m.a_dim,3),trans_lookup(d,m.d_dim,4)};
}
float reconstruct_trans(uint visibility,uint altitude,uint channel,TransParameters p,
                        device const Layout &m,device const float *data) {
  uint coefficients=m.trans_v+((visibility*m.trans_altitude.count+altitude)*m.channels+channel)*m.trans_rank;
  float values[4]={0,0,0,0}; uint index=0;
  for(uint a=p.altitude.index;a<=p.altitude.index+1;++a) if(a<m.a_dim)
    for(uint d=p.distance.index;d<=p.distance.index+1;++d) if(d<m.d_dim) {
      uint base=m.trans_u+((altitude*m.d_dim+d)*m.a_dim+a)*m.trans_rank;
      float value=0;
      for(uint r=0;r<m.trans_rank;++r) value+=data[base+r]*data[coefficients+r];
      values[index++]=value;
    }
  if(p.distance.factor>0) {
    values[0]=mix(values[0],values[1],p.distance.factor);
    values[1]=mix(values[2],values[3],p.distance.factor);
  }
  return p.altitude.factor>0 ? mix(max(0.f,values[0]),max(0.f,values[1]),p.altitude.factor) : max(0.f,values[0]);
}
struct TransQuery { Interpolation visibility, altitude; TransParameters distance; };
TransQuery prepare_transmission(Parameters p,float distance,device const Layout &m,device const float *data) {
  return {lookup(p.visibility,m.trans_visibility,data),lookup(p.altitude,m.trans_altitude,data),trans_parameters(p,distance,m)};
}
float transmission_prepared(TransQuery query,float wavelength,device const Layout &m,device const float *data) {
  if(wavelength<m.channel_start || wavelength>=m.channel_start+m.channels*m.channel_width) return 0;
  uint channel=uint((wavelength-m.channel_start)/m.channel_width);
  auto v=query.visibility,h=query.altitude;auto t=query.distance;
  float result[2];
  for(uint j=0;j<2;++j) {
    uint vi=min(v.index+j,m.trans_visibility.count-1);
    float value=reconstruct_trans(vi,h.index,channel,t,m,data);
    if(h.factor>0) value=mix(value,reconstruct_trans(vi,h.index+1,channel,t,m,data),h.factor);
    result[j]=value;
    if(v.factor==0) break;
  }
  float value=v.factor>0 ? mix(result[0],result[1],v.factor) : result[0];
  return value*value;
}
float transmission(Parameters p,float wavelength,float distance,device const Layout &m,device const float *data) {
  return transmission_prepared(prepare_transmission(p,distance,m,data),wavelength,m,data);
}
float sun(Parameters p,float wavelength,bool attenuate,device const Layout &m,device const float *data) {
  if(wavelength<280 || wavelength>=280+m.sun.count-1 || p.gamma>.004654793f) return 0;
  float index=wavelength-280;
  float value=linear(data,m.sun.offset+uint(index),fract(index));
  return attenuate ? value*transmission(p,wavelength,INFINITY,m,data) : value;
}
// Batch input uses already computed parameters; consumers may also call
// parameters() directly inside their own shaders. One thread per wavelength.
struct Query { Parameters parameters; float wavelength, distance; uint quantity, attenuate_sun; };
struct Geometry { float4 position, direction, settings; };
}
kernel void prague_parameters_batch(device const skymodelr_gpu::Geometry *input [[buffer(0)]],
    device skymodelr_gpu::Parameters *output [[buffer(1)]],constant uint &count [[buffer(2)]],
    uint i [[thread_position_in_grid]]) {
  if(i>=count) return;
  auto q=input[i];
  output[i]=skymodelr_gpu::parameters(q.position.xyz,q.direction.xyz,q.settings.x,q.settings.y,q.settings.z,q.settings.w);
}
kernel void prague_spectrum_batch(device const skymodelr_gpu::Layout &model [[buffer(0)]],
    device const float *data [[buffer(1)]], device const skymodelr_gpu::Query *queries [[buffer(2)]],
    device float *output [[buffer(3)]],constant uint &count [[buffer(4)]],uint i [[thread_position_in_grid]]) {
  if(i>=count) return;
  auto q=queries[i];
  output[i]=q.quantity==0 ? skymodelr_gpu::sky(q.parameters,q.wavelength,model,data) :
    (q.quantity==1 ? skymodelr_gpu::sun(q.parameters,q.wavelength,q.attenuate_sun,model,data) :
     skymodelr_gpu::transmission(q.parameters,q.wavelength,q.distance,model,data));
}
)PRAGUE"
