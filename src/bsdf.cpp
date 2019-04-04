#include "bsdf.h"

#include <iostream>
#include <algorithm>
#include <utility>

using std::min;
using std::max;
using std::swap;

namespace CGL {

void make_coord_space(Matrix3x3& o2w, const Vector3D& n) {
  Vector3D z = Vector3D(n.x, n.y, n.z);
  Vector3D h = z;
  if (fabs(h.x) <= fabs(h.y) && fabs(h.x) <= fabs(h.z)) h.x = 1.0;
  else if (fabs(h.y) <= fabs(h.x) && fabs(h.y) <= fabs(h.z)) h.y = 1.0;
  else h.z = 1.0;

  z.normalize();
  Vector3D y = cross(h, z);
  y.normalize();
  Vector3D x = cross(z, y);
  x.normalize();

  o2w[0] = x;
  o2w[1] = y;
  o2w[2] = z;
}

// Mirror BSDF //

Spectrum MirrorBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum MirrorBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO: 1.2
  // Using BSDF::reflect(), implement sample_f for a mirror surface
    * pdf = 1.0;
    reflect(wo, wi);
    return reflectance / abs_cos_theta(*wi);
}

// Microfacet BSDF //

double MicrofacetBSDF::G(const Vector3D& wo, const Vector3D& wi) {
    return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
}

double MicrofacetBSDF::D(const Vector3D& h) {
  // TODO: 2.2
  // Compute Beckmann normal distribution function (NDF) here.
  // You will need the roughness alpha.
    double tan_theta = sin_theta(h) / cos_theta(h);
    double num = exp(-pow(tan_theta, 2) / pow(alpha, 2));
    double den = PI * pow(alpha, 2) * pow(cos_theta(h), 4);
    return num / den;
//  return std::pow(cos_theta(h), 100.0);;
}

Spectrum MicrofacetBSDF::F(const Vector3D& wi) {
  // TODO: 2.3
  // Compute Fresnel term for reflection on dielectric-conductor interface.
  // You will need both eta and etaK, both of which are Spectrum.
    double cost = cos(getTheta(wi));
    Spectrum Rs = (((eta*eta) + (k*k)) - (2 * eta * cost) + (cost*cost)) / (((eta*eta) + (k*k)) + (2 * eta * cost) + (cost * cost));
    Spectrum Rp = (((eta*eta) + (k*k))*(cost*cost) - (2 * eta * cost) + 1 ) / ( ((eta*eta) + (k*k))*(cost * cost) + (2 * eta * cost) + 1);
    return ( Rs + Rp ) / 2.0 ;
//  return Spectrum();
}

Spectrum MicrofacetBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  // TODO: 2.1
  // Implement microfacet model here
    Vector3D n = Vector3D(0.0, 0.0, 1.0);
    Vector3D h = (wo + wi);
    h.normalize();
    if (dot(n, wi) > 0 && dot(n, wo) > 0){
        return (F(wi) * G(wo, wi) * D(h)) / (4.0 * (dot(n, wo)) * (dot(n, wi)));
    }
    return Spectrum(0,0,0);
}

Spectrum MicrofacetBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO: 2.4
  // *Importance* sample Beckmann normal distribution function (NDF) here.
  // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
  //       and return the sampled BRDF value.
    Vector2D r = sampler.get_sample();
    double r1 = r.x;
    double r2 = r.y;

    double theta = atan(sqrt(-(alpha*alpha)*log(1.0-r1)));
    double phi = 2 * PI * r2;

    double e = exp( - pow(tan(theta), 2) / pow(alpha,2));
    double thetaPDF = (2 * sin(theta) * e) / (pow(alpha,2) * pow(cos(theta), 3));
    double phiPDF = 1 / (2 * PI);

    Vector3D h = Vector3D(cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta));
    * wi = (2 * dot(wo, h) * h) - wo;
    if (dot(Vector3D(0,0,1), *wi) <= 0){
        *pdf = 0;
        return Spectrum(0,0,0);
    }
    double hPDF = (thetaPDF*phiPDF) / sin(theta);
    double wiPDF = hPDF / (4* dot(* wi, h));
    * pdf = wiPDF;

    return MicrofacetBSDF::f(wo, *wi);
//    
//  *wi = cosineHemisphereSampler.get_sample(pdf); //placeholder
//  return MicrofacetBSDF::f(wo, *wi);
}

// Refraction BSDF //

Spectrum RefractionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum RefractionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  return Spectrum();
}

// Glass BSDF //

Spectrum GlassBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum GlassBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // TODO: 1.4
  // Compute Fresnel coefficient and either reflect or refract based on it.
    if (refract(wo, wi, ior) == false){
        reflect(wo, wi);
        * pdf = 1.0;
        return reflectance / abs_cos_theta(* wi);
    }
    else{
        double R0 = pow((1.0 - ior)/(1.0 + ior), 2);
        double R = R0 + (1.0 - R0)*pow(1.0-abs(wo.z),5); //wo.z = +/- cos theta
        if (coin_flip(R)){
            reflect(wo, wi);
            * pdf = R;
            return R * reflectance / abs_cos_theta(*wi);
        }
        else{
            refract(wo, wi, ior);
            * pdf = 1.0 - R;
            double n;
            if (wo.z > 0){ // entering non-air material
                n = 1.0 / ior;
            }
            else{ // exiting
                n = ior;
            }
            return (1.0 - R) * transmittance / abs_cos_theta(*wi) / pow(n, 2);
        }
    }
  return Spectrum();
}

void BSDF::reflect(const Vector3D& wo, Vector3D* wi) {

  // TODO: 1.1
  // Implement reflection of wo about normal (0,0,1) and store result in wi.
    * wi = Vector3D(-wo.x, -wo.y, wo.z);
}

bool BSDF::refract(const Vector3D& wo, Vector3D* wi, float ior) {

  // TODO: 1.3
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
  // ray entering the surface through vacuum.
    double n;
    if (wo.z > 0){ // entering non-air material
        n = 1.0 / ior;
    }
    else{ // exiting
        n = ior;
    }
    
    double totIntRefl = 1.0 - pow(n, 2)*(1.0 - pow(wo.z, 2));
    if (totIntRefl < 0.0){
        return false;
    }
    else{
        double sign = -1.0;
        if (wo.z < 0.0){
            sign = 1.0;
        }
        * wi =  Vector3D( -n*wo.x, -n*wo.y, sign* sqrt(totIntRefl));
        return true;
    }
  return true;
}

// Emission BSDF //

Spectrum EmissionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum EmissionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *pdf = 1.0 / PI;
  *wi  = sampler.get_sample(pdf);
  return Spectrum();
}

} // namespace CGL
