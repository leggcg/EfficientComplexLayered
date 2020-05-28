#ifndef __PATH_HPP__
#define __PATH_HPP__

#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/hw/basicshader.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/scene.h>
#include <Eigen/Core>

#include "microfacet.h"
#include "ior.h"

/* Note that base material (diffuse, dielectric and conductor) and HG_phase functions are modified versions
                of Mitsuba's implementation, to avoid loading external plugins.
Some parts of the code are based on Belcour2018 work.
*/

MTS_NAMESPACE_BEGIN
inline Float DiffusePDF(const Vector3 &wi,const Vector3 &wo)
{
	return warp::squareToCosineHemispherePdf(wo);
}
Spectrum interactWithDiffuse(const Vector3 &wi, const Spectrum &reflectance, const Point2 &uv, Vector3 &wo) 
{
	wo = warp::squareToCosineHemisphere(uv);
	return reflectance;
}

Spectrum interactWithDiffusePDF(const Vector3 &wi, const Spectrum &reflectance, const Point2 &uv, 
				Vector3 &wo, Float &pdf) 
{
	wo = warp::squareToCosineHemisphere(uv);
	pdf = DiffusePDF(wi,wo);
	return reflectance;
}
inline Spectrum evalDiffuse(const Vector3 &wi,const Vector3 &wo, const Spectrum &reflectance)
{
	return reflectance * (INV_PI * Frame::cosTheta(wo));
}
Spectrum evalDiffusePDF(const Vector3 &wi,const Vector3 &wo, const Spectrum &reflectance, Float &pdf)
{
	pdf = DiffusePDF(wi,wo);
	return reflectance * (INV_PI * Frame::cosTheta(wo));
}

/* Interact with a conducting layer
 *
 * Inputs:
 *    + `wi` input direction
 *    + `eta   =   eta2 / eta1` upper index of refaction
 *    + `kappa = kappa2 / eta1` complex index of refraction
 *    + `alpha` roughness of the interface
 *    + `uv` two random numbers in [0, 1]
 *
 * Outputs:
 *    + `wo` outgoing direction
 *    + `return value` reflectance of the layer
 */
Spectrum interactWithConductor(const Vector3& wi,
                               const Spectrum& eta,
                               const Spectrum& kappa,
                               const Float &alphaU,
                               const Float &alphaV,
                               const Vector2& uv,
                               Vector3& wo,
                               const MicrofacetDistribution::EType type) {

   /* Construct the microfacet distribution matching the
      roughness values at the current surface position. */
   MicrofacetDistribution distr(type, alphaU, alphaV, true);

   /* Sample M, the microfacet normal */
   Float pdf; Point2 sample(uv.x, uv.y);
   Normal m = distr.sample(wi, sample, pdf);

   if (pdf == 0)
      return Spectrum(0.0f);

   /* Perfect specular reflection based on the microfacet normal */
   wo = reflect(wi, m);

   /* Side check */
   if (Frame::cosTheta(wo) <= 0)
      return Spectrum(0.0f);

   Spectrum F   = fresnelConductorExact(dot(wi, m), eta, kappa);
   //F = Spectrum(1.0f);
   Float weight = distr.smithG1(wo, m);

   return F * weight;
}
Spectrum interactWithConductorPDF(const Vector3& wi,
                               const Spectrum& eta,
                               const Spectrum& kappa,
                               const Float &alphaU,
                               const Float &alphaV,
                               const Vector2& uv,
                               Vector3& wo,
			       Float &pdf,
                               const MicrofacetDistribution::EType type) {

	/* Construct the microfacet distribution matching the
	   roughness values at the current surface position. */
	MicrofacetDistribution distr(type, alphaU, alphaV, true);

	/* Sample M, the microfacet normal */
	Point2 sample(uv.x, uv.y);
	Normal m = distr.sample(wi, sample, pdf);

	if (pdf == 0)
		return Spectrum(0.0f);

	/* Perfect specular reflection based on the microfacet normal */
	wo = reflect(wi, m);

	/* Side check */
	if (Frame::cosTheta(wo) <= 0)
		return Spectrum(0.0f);

	Spectrum F   = fresnelConductorExact(dot(wi, m), eta, kappa);
	//F = Spectrum(1.0f);
	Float weight = distr.smithG1(wo, m);

	pdf /= 4.0f * dot(wo, m);

	return F * weight;
}
Spectrum evalConductor(const Vector3& wi,const Vector3& wo,
		       const Spectrum& eta,
		       const Spectrum& kappa,
		       const Float &alphaU,
		       const Float &alphaV,
		       const MicrofacetDistribution::EType type) 
{
	/* Construct the microfacet distribution matching the
	   roughness values at the current surface position. */
	MicrofacetDistribution distr(type, alphaU, alphaV, true);
	
	/* Calculate the reflection half-vector */
	Vector H = normalize(wo+wi);

	/* Evaluate the microfacet normal distribution */
	const Float D = distr.eval(H);
	if (D == 0)
		return Spectrum(0.0f);
	
	/* Fresnel factor */
	const Spectrum F = fresnelConductorExact(dot(wi, H), eta, kappa);
	
	/* Smith's shadow-masking function */
	const Float G = distr.G(wi, wo, H);

	/* Calculate the total amount of reflection */
	Float model = D * G / (4.0f * Frame::cosTheta(wi));

	return F * model;
}

Spectrum evalConductorPDF(const Vector3& wi,const Vector3& wo,
		       const Spectrum& eta,
		       const Spectrum& kappa,
		       const Float &alphaU,
		       const Float &alphaV,
		       Float &pdf,
		       const MicrofacetDistribution::EType type) 
{
	/* Construct the microfacet distribution matching the
	   roughness values at the current surface position. */
	MicrofacetDistribution distr(type, alphaU, alphaV, true);
	
	/* Calculate the reflection half-vector */
	Vector H = normalize(wo+wi);

	/* Evaluate the microfacet normal distribution */
	const Float D = distr.eval(H);
	if (D == 0)
		return Spectrum(0.0f);
	
	/* Fresnel factor */
	const Spectrum F = fresnelConductorExact(dot(wi, H), eta, kappa);
	
	/* Smith's shadow-masking function */
	const Float G = distr.G(wi, wo, H);

	/* Calculate the total amount of reflection */
	Float denom = 4.0f * Frame::cosTheta(wi);
	Float model = D * G / denom;
	pdf = D * distr.smithG1(wi, H) / denom;

	return F * model;
}

Float ConductorPDF(const Vector3& wi,const Vector3& wo,
		       const Float &alphaU,
		       const Float &alphaV,
		       const MicrofacetDistribution::EType type) 
{
	//if(Frame::cosTheta(wi) <= 0 || Frame::cosTheta(wo) <= 0)
	//	return 0.0;

	/* Construct the microfacet distribution matching the
	   roughness values at the current surface position. */
	MicrofacetDistribution distr(type, alphaU, alphaV, true);
	
	/* Calculate the reflection half-vector */
	Vector H = normalize(wo+wi);

	return distr.eval(H) * distr.smithG1(wi, H) / (4.0f * Frame::cosTheta(wi));
}

/* Interact with a dielectric layer
 *
 * Inputs:
 *    + `wi` input direction
 *    + `eta   =   eta2 / eta1` upper index of refaction
 *    + `kappa = kappa2 / eta1` complex index of refraction
 *    + `alpha` roughness of the interface
 *    + `uv` two random numbers in [0, 1]
 *
 * Outputs:
 *    + `wo` outgoing direction
 *    + `is_reflected` light flux direction
 *    + `return value` reflectance of the layer
 */
Spectrum interactWithDielectric(const Vector3& wi,
                                const Float &eta,
                                const Float &alphaU,
                                const Float &alphaV,
                                const Vector3& uvw,
                                Vector3& wo,
                                bool& is_reflected,
                                const MicrofacetDistribution::EType type,
				const ETransportMode mode = ERadiance)
{

    /* Construct the microfacet distribution matching the
        roughness values at the current surface position. */
    MicrofacetDistribution distr(type, alphaU, alphaV, true);

    /* Sample M, the microfacet normal */
    Float microfacetPDF;
    Point2 sample(uvw.x, uvw.y);
    const Normal m = distr.sample(wi, sample, microfacetPDF);
    if (microfacetPDF == 0)
        return Spectrum(0.0f);

    Float cosThetaT;
    Float F = fresnelDielectricExt(dot(wi, m), cosThetaT, eta);
    Spectrum weight(1.0f);

    is_reflected = uvw.z <= F;

    if(is_reflected) {
      /* Perfect specular reflection based on the microfacet normal */
      wo = reflect(wi, m);

      /* Side check */
      if (Frame::cosTheta(wo) <= 0)
         return Spectrum(0.0f);

    } else {
        if (cosThetaT == 0)
            return Spectrum(0.0f);

        /* Perfect specular transmission based on the microfacet normal */
        wo = refract(wi, m, eta, cosThetaT);

        /* Side check */
        if (Frame::cosTheta(wo) >= 0)
            return Spectrum(0.0f);

	Float factor = (mode == ERadiance) ? (cosThetaT < 0 ? 1.0/eta : eta) : 1.0f;

	weight *= (factor * factor);

        //weight /= wo.z/wi.z;
   }
		
	if(mode == ETransportModes)
	{
		//Float auxeta = cosThetaT < 0 ? eta : 1.0/eta;
		//Float J = std::abs(wo.z / wi.z) * (auxeta * auxeta);
		//weight /= J;
		Float auxeta = cosThetaT < 0 ? 1.0/eta : eta;
		// we multiply by wi.z/wo.z so that the costheta is correct
		//	i.e. we traced in the other direction, so we have costheta for wo
		//	but we want it the costheta for wi as outgoing
		Float J = std::abs(wi.z / wo.z) * (auxeta * auxeta);
		weight *= J;
	}

	weight *= distr.smithG1(wo, m);
	if(!is_reflected)
		wo=-wo;
	return weight;
}
Spectrum interactWithDielectricPDF(const Vector3& wi,
                                const Float &eta,
                                const Float &alphaU,
                                const Float &alphaV,
                                const Vector3& uvw,
                                Vector3& wo,
				Float &pdf,
                                bool& is_reflected,
                                const MicrofacetDistribution::EType type,
				const ETransportMode mode = ERadiance) 
{

	/* Construct the microfacet distribution matching the
	   roughness values at the current surface position. */
	MicrofacetDistribution distr(type, alphaU, alphaV, true);

	/* Sample M, the microfacet normal */
	Float microfacetPDF;
	Point2 sample(uvw.x, uvw.y);
	const Normal m = distr.sample(wi, sample, microfacetPDF);
	if (microfacetPDF == 0)
		return Spectrum(0.0f);
	pdf = microfacetPDF;

	Float cosThetaT;
	Float F = fresnelDielectricExt(dot(wi, m), cosThetaT, eta);
	Spectrum weight(1.0f);

	is_reflected = uvw.z <= F;

	Float dwh_dwo;
	if(is_reflected) {
		/* Perfect specular reflection based on the microfacet normal */
		wo = reflect(wi, m);

		/* Side check */
		if (Frame::cosTheta(wo) <= 0)
			return Spectrum(0.0f);
      
		/* Jacobian of the half-direction mapping */
		dwh_dwo = 1.0f / (4.0f * dot(wo, m));
		pdf *= F;
	} else {
		if (cosThetaT == 0)
			return Spectrum(0.0f);

		/* Perfect specular transmission based on the microfacet normal */
		wo = refract(wi, m, eta, cosThetaT);

		/* Side check */
		if (Frame::cosTheta(wo) >= 0)
			return Spectrum(0.0f);

		Float factor = (mode == ERadiance) ? (cosThetaT < 0 ? 1.0/eta : eta) : 1.0f;

		weight *= (factor * factor);
	
		//weight /= wo.z/wi.z;

		Float auxeta = cosThetaT < 0 ? eta : 1.0/eta;
		Float sqrtDenom = dot(wi, m) + auxeta * dot(wo, m);
		dwh_dwo = (auxeta*auxeta * dot(wo, m)) / (sqrtDenom*sqrtDenom);
		pdf *= 1-F;
		
		if(mode == ETransportModes)
		{
			//Float J = std::abs(wo.z / wi.z) * (auxeta * auxeta);
			//weight /= J;
			const Float invauxeta = 1.0/auxeta;
			// we multiply by wi.z/wo.z so that the costheta is correct
			//	i.e. we traced in the other direction, so we have costheta for wo
			//	but we want it the costheta for wi as outgoing
			Float J = std::abs(wi.z / wo.z) * (invauxeta * invauxeta);
			weight *= J;
		}
	}
	pdf *= std::abs(dwh_dwo);

	weight *= distr.smithG1(wo, m);
	if(!is_reflected)
		wo=-wo;
	
	return weight;
}
Spectrum evalDielectric(const Vector3& wi, const Vector3& wo,
			const Float &eta,
			const Float &alphaU,
			const Float &alphaV,
			const MicrofacetDistribution::EType type,
			const ETransportMode mode = ERadiance) 
{
	/* Determine the type of interaction */
	bool reflect = Frame::cosTheta(wi) * Frame::cosTheta(wo) > 0;

	Vector H;
	if (reflect) {
		/* Calculate the reflection half-vector */
		H = normalize(wo+wi);
	} else {
		/* Calculate the transmission half-vector */
		Float eta_ = Frame::cosTheta(wi) > 0 ? eta : 1.0/eta;
		H = normalize(wi + wo*eta_);
	}

	/* Ensure that the half-vector points into the
	   same hemisphere as the macrosurface normal */
	H *= math::signum(Frame::cosTheta(H));

	/* Construct the microfacet distribution matching the
	   roughness values at the current surface position. */
	MicrofacetDistribution distr(type, alphaU, alphaV, true);

	/* Evaluate the microfacet normal distribution */
	const Float D = distr.eval(H);
	if (D == 0)
		return Spectrum(0.0f);

	/* Fresnel factor */
	const Float F = fresnelDielectricExt(dot(wi, H), eta);

	/* Smith's shadow-masking function */
	const Float G = distr.G(wi, wo, H);

	if (reflect) {
		/* Calculate the total amount of reflection */
		Float value = F * D * G / (4.0f * std::abs(Frame::cosTheta(wi)));

		return Spectrum(value);
	} else {
		Float eta_ = Frame::cosTheta(wi) > 0.0f ? eta : 1.0/eta;

		/* Calculate the total amount of transmission */
		Float sqrtDenom = dot(wi, H) + eta_ * dot(wo, H);
		Float value = ((1 - F) * D * G * eta_ * eta_ * dot(wi, H) * dot(wo, H)) /
			(Frame::cosTheta(wi) * sqrtDenom * sqrtDenom);

		/* Missing term in the original paper: account for the solid angle
		   compression when tracing radiance -- this is necessary for
		   bidirectional methods */
		Float factor = (mode == ERadiance) ? (Frame::cosTheta(wi) > 0 ? 1.0/eta : eta) : 1.0f;

		return Spectrum(std::abs(value * factor * factor));
	}
}
Spectrum evalDielectricPDF(const Vector3& wi, const Vector3& wo,
			const Float &eta,
			const Float &alphaU,
			const Float &alphaV,
			Float &pdf,
			const MicrofacetDistribution::EType type,
			const ETransportMode mode = ERadiance) 
{
	/* Determine the type of interaction */
	bool reflect = Frame::cosTheta(wi) * Frame::cosTheta(wo) > 0;

	Vector H;
	Float dwh_dwo;
	if (reflect) {
		/* Calculate the reflection half-vector */
		H = normalize(wo+wi);
		/* Jacobian of the half-direction mapping */
		dwh_dwo = 1.0f / (4.0f * dot(wo, H));
	} else {
		/* Calculate the transmission half-vector */
		Float eta_ = Frame::cosTheta(wi) > 0 ? eta : 1.0/eta;
		H = normalize(wi + wo*eta_);
		/* Jacobian of the half-direction mapping */
		Float sqrtDenom = dot(wi, H) + eta_ * dot(wo, H);
		dwh_dwo = (eta_*eta_ * dot(wo, H)) / (sqrtDenom*sqrtDenom);
	}

	/* Ensure that the half-vector points into the
	   same hemisphere as the macrosurface normal */
	H *= math::signum(Frame::cosTheta(H));

	/* Construct the microfacet distribution matching the
	   roughness values at the current surface position. */
	MicrofacetDistribution distr(type, alphaU, alphaV, true);

	/* Evaluate the microfacet normal distribution */
	const Float D = distr.eval(H);
	if (D == 0)
	{
		pdf = 0.0;
		return Spectrum(0.0f);
	}

	/* Fresnel factor */
	const Float F = fresnelDielectricExt(dot(wi, H), eta);

	/* Smith's shadow-masking function */
	const Float G = distr.G(wi, wo, H);

	/* Evaluate the microfacet model sampling density function */
	Float prob = distr.pdf(math::signum(Frame::cosTheta(wi)) * wi, H);
	prob *= reflect ? F : (1-F);
	pdf = std::abs(prob * dwh_dwo);

	if (reflect) {
		/* Calculate the total amount of reflection */
		Float value = F * D * G / (4.0f * std::abs(Frame::cosTheta(wi)));

		return Spectrum(value);
	} else {
		Float eta_ = Frame::cosTheta(wi) > 0.0f ? eta : 1.0/eta;

		/* Calculate the total amount of transmission */
		Float sqrtDenom = dot(wi, H) + eta_ * dot(wo, H);
		Float value = ((1 - F) * D * G * eta_ * eta_ * dot(wi, H) * dot(wo, H)) /
			(Frame::cosTheta(wi) * sqrtDenom * sqrtDenom);

		/* Missing term in the original paper: account for the solid angle
		   compression when tracing radiance -- this is necessary for
		   bidirectional methods */
		Float factor = (mode == ERadiance) ? (Frame::cosTheta(wi) > 0 ? 1.0/eta : eta) : 1.0f;
		Float J;
		if(mode == ETransportModes)
		{
			Float auxeta = wo.z > 0 ? eta : Float(1.0) / eta;
			J=std::abs(wi.z / wo.z) * (auxeta*auxeta);
		}
		else
		{
			J=1.0f;
		}

		return Spectrum(std::abs(value * factor * factor * J));
	}
}
Float DielectricPDF(const Vector3& wi, const Vector3& wo,
			const Float &eta,
			const Float &alphaU,
			const Float &alphaV,
			const MicrofacetDistribution::EType type)
{
	/* Determine the type of interaction */
	bool reflect = Frame::cosTheta(wi) * Frame::cosTheta(wo) > 0;

	Vector H;
	Float dwh_dwo;
	if (reflect) {
		/* Calculate the reflection half-vector */
		H = normalize(wo+wi);
		/* Jacobian of the half-direction mapping */
		dwh_dwo = 1.0f / (4.0f * dot(wo, H));
	} else {
		/* Calculate the transmission half-vector */
		Float eta_ = Frame::cosTheta(wi) > 0 ? eta : 1.0/eta;
		H = normalize(wi + wo*eta_);
		/* Jacobian of the half-direction mapping */
		Float sqrtDenom = dot(wi, H) + eta_ * dot(wo, H);
		dwh_dwo = (eta_*eta_ * dot(wo, H)) / (sqrtDenom*sqrtDenom);
	}

	/* Ensure that the half-vector points into the
	   same hemisphere as the macrosurface normal */
	H *= math::signum(Frame::cosTheta(H));

	/* Construct the microfacet distribution matching the
	   roughness values at the current surface position. */
	MicrofacetDistribution distr(type, alphaU, alphaV, true);

	/* Fresnel factor */
	const Float F = fresnelDielectricExt(dot(wi, H), eta);

	/* Evaluate the microfacet model sampling density function */
	Float prob = distr.pdf(math::signum(Frame::cosTheta(wi)) * wi, H);
	prob *= reflect ? F : (1-F);
	return std::abs(prob * dwh_dwo);
}

/* Sample the HG phase function
 */
// same code as src/phase/hg.cpp:74
Vector3 HG_Sample(const Vector3& wi, Float g, const Vector2& uv) {

   Float cosTheta;
   if (std::abs(g) < Epsilon) {
      cosTheta = 1 - 2*uv.x;
   } else {
      Float sqrTerm = (1 - g * g) / (1 - g + 2 * g * uv.x);
      cosTheta = (1 + g * g - sqrTerm * sqrTerm) / (2 * g);
   }

   Float sinTheta = math::safe_sqrt(1.0f-cosTheta*cosTheta),
         sinPhi, cosPhi;

   math::sincos(2*M_PI*uv.y, &sinPhi, &cosPhi);

   const Vector ws(sinTheta * cosPhi, sinTheta * sinPhi, cosTheta);
   return Frame(wi).toWorld(ws);
}
Vector3 HG_SamplePDF(const Vector3& wi, Float g, const Vector2& uv, Float &pdf) {

   Float cosTheta;
   if (std::abs(g) < Epsilon) {
      cosTheta = 1 - 2*uv.x;
   } else {
      Float sqrTerm = (1 - g * g) / (1 - g + 2 * g * uv.x);
      cosTheta = (1 + g * g - sqrTerm * sqrTerm) / (2 * g);
   }

   Float sinTheta = math::safe_sqrt(1.0f-cosTheta*cosTheta),
         sinPhi, cosPhi;

   math::sincos(2*M_PI*uv.y, &sinPhi, &cosPhi);

   const Vector ws(sinTheta * cosPhi, sinTheta * sinPhi, cosTheta);
   Vector wo( Frame(wi).toWorld(ws));

   Float temp = 1.0f + g*g + 2.0f * g * dot(wi, wo);
   pdf = INV_FOURPI * (1 - g*g) / (temp * std::sqrt(temp));

   return wo;
}
Float HG_Eval(const Vector3& wi, const Vector3& wo, const Float &g)
{
	Float temp = 1.0f + g*g + 2.0f * g * dot(wi, wo);
	return INV_FOURPI * (1 - g*g) / (temp * std::sqrt(temp));
}

Spectrum evalTr(const Float &depth, const Spectrum &sigma_t, const Vector &wi)
{
	const Float opt_depth = depth/Frame::cosTheta(wi);
	if(std::isnan(opt_depth) || std::isinf(opt_depth)) { return Spectrum(0.0); }

	return (sigma_t * (-opt_depth)).exp();
}
Spectrum evalTr(const Float &distance, const Spectrum &sigma_t)
{
	return (sigma_t * (-distance)).exp();
}
Float MediumPdfFailure(const Spectrum &sigma_t, const Float &depth, const Vector &wi)
{
	const Float opt_depth = depth/std::abs(Frame::cosTheta(wi));
	if(std::isnan(opt_depth) || std::isinf(opt_depth)) { return 0.0; }

	// Sampling pdf failure
	//const float pdf = exp(-sigma_f * opt_depth);
	return (sigma_t * (-opt_depth)).exp().average();
}
/// Helper function: reflect \c wi with respect to a given surface normal
/*inline Vector reflect(const Vector &wi, const Normal &m) 
{
	return 2 * dot(wi, m) * Vector(m) - wi;
}*/

Frame getFrame(const Intersection &its, const Vector &n) 
{
	Frame result;

	//Frame frame;
	//BSDF::getFrame(its,frame);
	//Frame frame=BSDF::getFrame(its);
	Frame frame;
	computeShadingFrame(its.shFrame.n, its.dpdu, frame);

	result.n = normalize(frame.toWorld(n));

	result.s = normalize(its.dpdu - result.n
	         * dot(result.n, its.dpdu));

	result.t = cross(result.n, result.s);

	return result;
}


struct PathInfo
{
	// position
	enum EType {
		EBSDF  = 0,
		EPhase = 1//,
		//EVolume = 2
	};
	EType type;
	int layer;
	Float depth;

	// path info
	int bounce; 
	Spectrum thruBefore;
	Spectrum thruAfter;

	// direction
	Vector wi;
	bool is_up;
	bool is_reflected;

	// sample event
	Float pdf;
	Vector wo;
	//Float d;

	std::string toString() const
	{
		std::ostringstream oss;
		oss << "PathInfo[" << endl
		    << "  Type = \"" 
		       << (type==EType::EBSDF ? "BSDF" : "Phase") << endl
		    << "  Layer = " << layer << endl
		    << "  Depth = " << depth << endl
		    << "  Bounce = " << bounce << endl
		    << "  thruBefore = " << thruBefore.toString() << endl
		    << "  thruAfter  = " << thruAfter.toString() << endl
		    << "  wi = " << wi.toString() << endl
		    << "  wo = " << wo.toString() << endl
		    << "  is_up = " <<  (is_up ? "true" : "false") << endl
		    << "  is_reflected = " <<  (is_reflected ? "true" : "false") << endl
		    << "  pdf = " << pdf << endl
		    //<< "  d = " << d << endl
		    << "]";
		return oss.str();
	}

	PathInfo(const EType &type_,const int &bounce_, const int &layer_, const bool &is_up_, const bool &is_reflctd_,	
			const Vector &wi_, const Vector &wo_, const Float &pdf_,
			const Spectrum &thruBefore_, const Spectrum &thruAfter_,
			const Float &depth_ = 0.0)
	{
		type=type_;
		bounce=bounce_;
		layer=layer_;
		is_up=is_up_;
		is_reflected=is_reflctd_;
		wi=wi_;
		wo=wo_;
		pdf=pdf_;
		thruBefore=thruBefore_;
		thruAfter=thruAfter_;
		if(depth_==0.0) // BSDF
			depth=layer_;//layer/2;
		else
			depth=depth_;
/*		else // Phase function
		{
			if(is_up)
				depth=(layer+1)/2 - depth_;
			else
				depth=layer/2 + depth_;
		}*/
	}

/*	PathInfo(const int &bounce_, const int &layer_, const bool &is_up_, 
			const Vector &wi_, const Float &pdf_, const Float &d_,
			const Spectrum &thruBefore_, const Spectrum &thruAfter_)
	{
		type=EVolume;
		bounce=bounce_;
		layer=layer_;
		is_up=is_up_;
		wi=wi_;
		wo=wi_;
		pdf=pdf_;
		d=d_;
		thruBefore=thruBefore_;
		thruAfter=thruAfter_;
	}*/

	void setVertex(const EType &type_,const int &bounce_, const int &layer_, const bool &is_up_, 
			const Vector &wi_, const Vector &wo_, const Float &pdf_,
			const Spectrum &thruBefore_, const Spectrum &thruAfter_)
	{
		type=type_;
		bounce=bounce_;
		layer=layer_;
		is_up=is_up_;
		wi=wi_;
		wo=wo_;
		pdf=pdf_;
		thruBefore=thruBefore_;
		thruAfter=thruAfter_;
	}

/*	void setSegment(const int &bounce_, const int &layer_, const bool &is_up_, 
			const Vector &wi_, const Float &pdf_, const Float &d_,
			const Spectrum &thruBefore_, const Spectrum &thruAfter_)
	{
		type=EVolume;
		bounce=bounce_;
		layer=layer_;
		is_up=is_up_;
		wi=wi_;
		wo=wi_;
		pdf=pdf_;
		d=d_;
		thruBefore=thruBefore_;
		thruAfter=thruAfter_;
	}*/
};


class Path 
{
public:
Path(BSDFSamplingRecord &bRec_, const int &max_bounces_ = 100) : bRec(bRec_), max_bounces(max_bounces_)
{
}

~Path()
{
	for(unsigned int i=0;i<m_path.size();i++)
		delete m_path[i];
}

// WE MUST KEEP THIS CODE IN SYNC WITH sampleDirection()
// -- all changes should go in this version first then copy paste
// -- USE ifdef/ifndef SAMPLEDIRECTION to code the portions relevant to one or the other
bool generatePath(const int &nb_layers, const std::vector<Float> &m_depths, const std::vector<Float> &m_gs,
		const MicrofacetDistribution::EType &m_type,
		const Spectrum *m_etas, const Spectrum *m_kappas, const Float *m_alphasU, const Float *m_alphasV,
		const Spectrum *m_reflectances, const Vector *m_normals,
		const Spectrum *m_sigmas, const Spectrum *m_sigmaa,
		const ETransportMode mode = ERadiance)
{
	if (Frame::cosTheta(bRec.wi) < 0 ||
		((bRec.component != -1 && bRec.component != 0) ||
		!(bRec.typeMask & BSDF::EGlossyReflection)))
			return false;

	Vector3 wi   = bRec.wi;
	int index    = 0;
	int bounces  = 0;
	Spectrum e   = Spectrum(1.0);
	#ifndef SAMPLEDIRECTION
	Spectrum epre= Spectrum(1.0);
	#endif
	bool is_up   = false;

	// Track the current uv position
	const Intersection &its=bRec.its;
	//Point2 uvposition=its.uv;
	//const Vector &dpdu=its.dpdu;
	//const Vector &dpdv=its.dpdv;

	/* Start the iteration loop */
	while(bounces < max_bounces) 
	{
		/* Retreive the layer depth. If set, we assume that the layer is a scattering medium
		 * and sample the transmittance. We can as well do scattering in the medium
		 */
		Float depth = m_depths[index];
		if(depth > 0.0f) 
		{
			/* Obtain the meidum information */
			const Spectrum &sigma_a = m_sigmaa[index];
			const Spectrum &sigma_s = m_sigmas[index];
			const Spectrum &sigma_t = (sigma_a+sigma_s);
			const Float    &g       = m_gs[index];
			// use balance strategy for sampling
			int channel = std::min((int) (nextRandom()*SPECTRUM_SAMPLES), SPECTRUM_SAMPLES-1);
			const Float    &sigma_f = sigma_t[channel];

			/* Remaining depth along the layer dimension */
			Float rmn_depth = depth;

			//#define DEBUG_LAYERINFO
			#ifdef DEBUG_LAYERINFO		
			printf("Bounce: %d Layer: %d(depth=%f isup=%d) uv=%s e=%s ", bounces, index, depth, is_up, bRec.its.uv.toString().c_str(), e.toString().c_str());
			printf("sigma_a=%s, sigma_s=%s, g=%f\n", sigma_a.toString().c_str(), sigma_s.toString().c_str(), g);
			exit(-1);
			#endif


			/* Perform multiple scattering in the layer */
			while(true) 
			{
				/* Sample the transmittance */
				Float d  = - math::fastlog(1.0 - nextRandom()) / sigma_f;
				Float rd = d * Frame::cosTheta(wi); // Depth along the layer dimension


				/* Should we stop the ray and sample the phase function? */
				bool is_stopped   = rd < rmn_depth;
				if(is_stopped && sigma_f>0.0) 
				{
					Spectrum transmittance=(sigma_t * (-d)).exp();
					if(transmittance.max() < 1e-10)
						return false;

					// Sampling distance pdf (success) balance
					Float pdf = (sigma_t*transmittance).average();//sigma_f * exp(-sigma_f *d);
					if(pdf < 1e-10)
						return false;

					// sigma_s * transmittance/pdfsuccess (volpath.cpp:113)
					#ifndef SAMPLEDIRECTION
					epre=e;
					#endif
					e *= sigma_s * transmittance * (1.0 / pdf);
					//e *= sigma_s*(- sigma_t * d).exp() / pdf;
					//Bench Assert(e.isValid());

					// Add segment to the path
					//m_path.push_back(new PathInfo(bounces, index, is_up, wi, pdf, d, epre, e));

					// Sample the HG function
					Vector2 uv = Vector2(nextRandom(), nextRandom());
					#ifndef SAMPLEDIRECTION
					Vector3 wo = HG_SamplePDF(wi, g, uv, pdf);
					#else
					Vector3 wo = HG_Sample(wi, g, uv);
					#endif
					// no need to update throughput, as hg phase sample returns 1.0

					// Check if the new ray is up
					bool is_wo_up = wo.z >= 0.0;
					bool is_wi_up = wi.z >= 0.0;
					bool chg_dir  = is_wo_up != is_wi_up;
					if(chg_dir) {
						wo = -wo;
					}
					if(wo.z == 0.0f) return false;

					#ifndef SAMPLEDIRECTION
					// Add vertex to the Path
					//m_path.push_back(new PathInfo(PathInfo::EPhase, bounces, index, is_up, chg_dir, wi, wo, pdf, e, e, (chg_dir) ? rd : rmn_depth - rd));
					m_path.push_back(new PathInfo(PathInfo::EPhase, bounces, index, is_up, chg_dir, wi, wo, pdf, e, e, rmn_depth - rd));
					#endif

					// Update transport information
					bounces  += 1;
					wi        = wo;
					is_up     = (chg_dir) ? !is_up : is_up;
					//rmn_depth = (chg_dir) ? rd : rmn_depth - rd;
					rmn_depth = (chg_dir) ? depth - (rmn_depth - rd) : rmn_depth - rd;

					// update theoretical uv coordinate
					//uvposition=getUV(uvposition, wo, d, dpdu, dpdv);

					/* Leave to the next interface */
				} 
				else 
				{
					const Float opt_depth = rmn_depth/Frame::cosTheta(wi);
					if(std::isnan(opt_depth) || std::isinf(opt_depth)) { return false; }

					// Sampling pdf failure
					//const float pdf = exp(-sigma_f * opt_depth);
					const Spectrum &transmittance = (sigma_t * (-opt_depth)).exp();
					const float &pdf = transmittance.average();
					if(transmittance.max() < 1e-10)
						return false;
					if(pdf < 1e-10)
						return false;

					// transmittance / pdffailure (volpath:189)
					#ifndef SAMPLEDIRECTION
					epre=e;
					#endif
					e *= transmittance * (1.0 / pdf);
					//e *= (sigma_t * (-opt_depth)).exp() / pdf;
					//Bench Assert(e.isValid());
					
					// Add segment to the path
					//m_path.push_back(new PathInfo(bounces, index, is_up, wi, pdf, opt_depth, epre, e));
			
					//uvposition=getUV(uvposition, wi, opt_depth, dpdu, dpdv);
					break;
				}
				/* Short stop if the ray has no energy */
				if(e.isZero() || !e.isValid()) { return false; }
				//Bench Assert(e.isValid());
			}
		} 
		else 
		{
			/* Obtain the material information */
			const float &alphaU    = m_alphasU[index];
			const float &alphaV    = m_alphasV[index];
			const bool is_diffuse  = alphaU < 1e-10;
			const Spectrum &eta1   = (!is_up) ?   m_etas[index]   :   m_etas[index+1] ;
			const Spectrum &eta2   = (!is_up) ?   m_etas[index+1] :   m_etas[index]   ;
			const Spectrum &kappa2 = (!is_up) ? m_kappas[index+1] : m_kappas[index]   ;
			const Spectrum &reflectance = m_reflectances[index];

			// get the normal from the normal map, keep in mind that we need to track  the "actual" 
			//	intersection within the texture not just the intersection of the surface (bRec.its)
			//	----> we have the actual position tracked in uvposition
			const Vector &n=m_normals[index];//getNormal(index,uvposition,is_up);
			const Frame &perturbed=getFrame(its,n);
			const Vector &wiperturbed=perturbed.toLocal(its.toWorld(wi));


			const Spectrum &eta    = eta2   / eta1;
			const Spectrum &kappa  = kappa2 / eta1;
		
			#ifdef DEBUG_LAYERINFO		
			printf("Bounce: %d Layer: %d(isup=%d) uv=%s e=%s", bounces, index, is_up,
					bRec.its.uv.toString().c_str(),e.toString().c_str());
			printf("alpha=%f, eta1=%s, eta2=%s, kappa2=%s\n", alpha, eta1.toString().c_str(), 
					eta2.toString().c_str(), kappa2.toString().c_str());
			#endif

			/* Get the next interaction with the layer */
			Vector3 wo;
			if(is_diffuse) // lambertian
			{
				Point2 uv(nextRandom(), nextRandom());

				#ifndef SAMPLEDIRECTION
				Float pdf;
				epre=e;
				e *= interactWithDiffusePDF(wiperturbed, reflectance, uv, wo, pdf);
				#else
				e *= interactWithDiffuse(wiperturbed, reflectance, uv, wo);
				#endif
				//Bench Assert(e.isValid());
				wo = its.toLocal(perturbed.toWorld(wo));

				#ifndef SAMPLEDIRECTION
				// Add vertex to the Path
				m_path.push_back(new PathInfo(PathInfo::EBSDF, bounces, index, is_up, true, wi, wo, pdf, epre, e));
				#endif

				is_up = !is_up;
			}
			else if(kappa.isZero()) // dielectric
			{	
				Vector3 uvw(nextRandom(), nextRandom(), nextRandom());
				bool is_reflected = true;
				#ifndef SAMPLEDIRECTION
				Float pdf;
				epre=e;
				e *= interactWithDielectricPDF(wiperturbed, eta.average(), alphaU, alphaV, uvw, wo, pdf, is_reflected, m_type, mode);
				#else
				e *= interactWithDielectric(wiperturbed, eta.average(), alphaU, alphaV, uvw, wo, is_reflected, m_type, mode);
				#endif
				//Bench Assert(e.isValid());
				wo = its.toLocal(perturbed.toWorld(wo));

				#ifndef SAMPLEDIRECTION
				// Add vertex to the Path
				m_path.push_back(new PathInfo(PathInfo::EBSDF, bounces, index, is_up, is_reflected, wi, wo, pdf, epre, e));
				#endif
				
				if(is_reflected) {
					is_up = !is_up;
				}
			} 
			else // conductor
			{
				Vector2 uv(nextRandom(), nextRandom());
				#ifndef SAMPLEDIRECTION
				Float pdf;
				epre=e;
				e *= interactWithConductorPDF(wiperturbed, eta, kappa, alphaU, alphaV, uv, wo, pdf, m_type) * reflectance;
				#else
				e *= interactWithConductor(wiperturbed, eta, kappa, alphaU, alphaV, uv, wo, m_type) * reflectance;
				#endif
				//Bench Assert(e.isValid());
				wo = its.toLocal(perturbed.toWorld(wo));

				#ifndef SAMPLEDIRECTION
				// Add vertex to the Path
				m_path.push_back(new PathInfo(PathInfo::EBSDF, bounces, index, is_up, true, wi, wo, pdf, epre, e));
				#endif

				is_up = !is_up;
			}

			if(e.isZero() || !e.isValid()) { return false; }

			if(wo.z<0.0) // Avoid numerical errors or change of side due to normal map
				return false;

			/* If the ray continues its travel, updates its parameters */
			wi       = wo;
			bounces += 1;
		}

		/* Update transport information */
		index   += (is_up) ? -1 : 1;

		if(index < 0) {
			// BEGIN COMMENTED TO COMPILE CAREFUL
			/*bRec.wo = wi;
			bRec.eta = 1.0f;
			bRec.sampledComponent = 0;
			bRec.sampledType = EGlossyReflection;
			pdf = this->pdf(bRec, ESolidAngle);*/
			// END COMMENTED TO COMPILE CAREFUL

			// last node of m_path contains:
			// Total path throughput: thruAfter
			// sampled direction: wo
			//Bench Assert(e.isValid());
			#ifdef SAMPLEDIRECTION
			wo_=wi;
			e_=e;
			#endif
			return true;
		} 
		else if(index >= nb_layers) {
			return false;
		}
	}

	return false;
}
// returns Path throughput and sampled direction wo
// WARNING: only works if invoked after a buildPath() call 
Spectrum getPathSample(Vector &wo) const
{
	const PathInfo *pinfo=m_path.back();
	wo=pinfo->wo;
	return pinfo->thruAfter;
}

// WE MUST KEEP THIS CODE IN SYNC WITH buildPath()
// EASIEST WAY: we just need to remove the push_backs and epre, final if do e_=e, wo_=wo
//		- change calls to non PDF version
// this is used for BSDF::sample(), we must do all path tracing
//	we ONLY care for path throughput e (sampling weight) and outgoing direction wo
//	no need to store path information in m_path
#define SAMPLEDIRECTION
bool sampleDirection(Spectrum &e_, Vector &wo_, const int &nb_layers, 
		const std::vector<Float> &m_depths, const std::vector<Float> &m_gs,
		const MicrofacetDistribution::EType &m_type,
		const Spectrum *m_etas, const Spectrum *m_kappas, const Float *m_alphasU, const Float *m_alphasV,
		const Spectrum *m_reflectances, const Vector *m_normals,
		const Spectrum *m_sigmas, const Spectrum *m_sigmaa,
		const ETransportMode mode = ERadiance)
{
	if (Frame::cosTheta(bRec.wi) < 0 ||
		((bRec.component != -1 && bRec.component != 0) ||
		!(bRec.typeMask & BSDF::EGlossyReflection)))
			return false;

	Vector3 wi   = bRec.wi;
	int index    = 0;
	int bounces  = 0;
	Spectrum e   = Spectrum(1.0);
	#ifndef SAMPLEDIRECTION
	Spectrum epre= Spectrum(1.0);
	#endif
	bool is_up   = false;

	// Track the current uv position
	const Intersection &its=bRec.its;
	//Point2 uvposition=its.uv;
	//const Vector &dpdu=its.dpdu;
	//const Vector &dpdv=its.dpdv;

	/* Start the iteration loop */
	while(bounces < max_bounces) 
	{
		/* Retreive the layer depth. If set, we assume that the layer is a scattering medium
		 * and sample the transmittance. We can as well do scattering in the medium
		 */
		Float depth = m_depths[index];
		if(depth > 0.0f) 
		{
			/* Obtain the meidum information */
			const Spectrum &sigma_a = m_sigmaa[index];
			const Spectrum &sigma_s = m_sigmas[index];
			const Spectrum &sigma_t = (sigma_a+sigma_s);
			const Float    &g       = m_gs[index];
			// use balance strategy for sampling
			int channel = std::min((int) (nextRandom()*SPECTRUM_SAMPLES), SPECTRUM_SAMPLES-1);
			const Float    &sigma_f = sigma_t[channel];

			/* Remaining depth along the layer dimension */
			Float rmn_depth = depth;

			//#define DEBUG_LAYERINFO
			#ifdef DEBUG_LAYERINFO		
			printf("Bounce: %d Layer: %d(depth=%f isup=%d) uv=%s e=%s ", bounces, index, depth, is_up, bRec.its.uv.toString().c_str(), e.toString().c_str());
			printf("sigma_a=%s, sigma_s=%s, g=%f\n", sigma_a.toString().c_str(), sigma_s.toString().c_str(), g);
			exit(-1);
			#endif


			/* Perform multiple scattering in the layer */
			while(true) 
			{
				/* Sample the transmittance */
				Float d  = - math::fastlog(1.0 - nextRandom()) / sigma_f;
				Float rd = d * Frame::cosTheta(wi); // Depth along the layer dimension


				/* Should we stop the ray and sample the phase function? */
				bool is_stopped   = rd < rmn_depth;
				if(is_stopped && sigma_f>0.0) 
				{
					Spectrum transmittance=(sigma_t * (-d)).exp();
					if(transmittance.max() < 1e-10)
						return false;

					// Sampling distance pdf (success) balance
					Float pdf = (sigma_t*transmittance).average();//sigma_f * exp(-sigma_f *d);
					if(pdf < 1e-10)
						return false;

					// sigma_s * transmittance/pdfsuccess (volpath.cpp:113)
					#ifndef SAMPLEDIRECTION
					epre=e;
					#endif
					e *= sigma_s * transmittance * (1.0 / pdf);
					//e *= sigma_s*(- sigma_t * d).exp() / pdf;
					//Bench Assert(e.isValid());

					// Add segment to the path
					//m_path.push_back(new PathInfo(bounces, index, is_up, wi, pdf, d, epre, e));

					// Sample the HG function
					Vector2 uv = Vector2(nextRandom(), nextRandom());
					#ifndef SAMPLEDIRECTION
					Vector3 wo = HG_SamplePDF(wi, g, uv, pdf);
					#else
					Vector3 wo = HG_Sample(wi, g, uv);
					#endif
					// no need to update throughput, as hg phase sample returns 1.0

					// Check if the new ray is up
					bool is_wo_up = wo.z >= 0.0;
					bool is_wi_up = wi.z >= 0.0;
					bool chg_dir  = is_wo_up != is_wi_up;
					if(chg_dir) {
						wo = -wo;
					}
					if(wo.z == 0.0f) return false;

					#ifndef SAMPLEDIRECTION
					// Add vertex to the Path
					//m_path.push_back(new PathInfo(PathInfo::EPhase, bounces, index, is_up, chg_dir, wi, wo, pdf, e, e, (chg_dir) ? rd : rmn_depth - rd));
					m_path.push_back(new PathInfo(PathInfo::EPhase, bounces, index, is_up, chg_dir, wi, wo, pdf, e, e, rmn_depth - rd));
					#endif

					// Update transport information
					bounces  += 1;
					wi        = wo;
					is_up     = (chg_dir) ? !is_up : is_up;
					//rmn_depth = (chg_dir) ? rd : rmn_depth - rd;
					rmn_depth = (chg_dir) ? depth - (rmn_depth - rd) : rmn_depth - rd;

					// update theoretical uv coordinate
					//uvposition=getUV(uvposition, wo, d, dpdu, dpdv);

					/* Leave to the next interface */
				} 
				else 
				{
					const Float opt_depth = rmn_depth/Frame::cosTheta(wi);
					if(std::isnan(opt_depth) || std::isinf(opt_depth)) { return false; }

					// Sampling pdf failure
					//const float pdf = exp(-sigma_f * opt_depth);
					const Spectrum &transmittance = (sigma_t * (-opt_depth)).exp();
					const float &pdf = transmittance.average();
					if(transmittance.max() < 1e-10)
						return false;
					if(pdf < 1e-10)
						return false;

					// transmittance / pdffailure (volpath:189)
					#ifndef SAMPLEDIRECTION
					epre=e;
					#endif
					e *= transmittance * (1.0 / pdf);
					//e *= (sigma_t * (-opt_depth)).exp() / pdf;
					//Bench Assert(e.isValid());
					
					// Add segment to the path
					//m_path.push_back(new PathInfo(bounces, index, is_up, wi, pdf, opt_depth, epre, e));
			
					//uvposition=getUV(uvposition, wi, opt_depth, dpdu, dpdv);
					break;
				}
				/* Short stop if the ray has no energy */
				if(e.isZero() || !e.isValid()) { return false; }
				//Bench Assert(e.isValid());
			}
		} 
		else 
		{
			/* Obtain the material information */
			const float &alphaU    = m_alphasU[index];
			const float &alphaV    = m_alphasV[index];
			const bool is_diffuse  = alphaU < 1e-10;
			const Spectrum &eta1   = (!is_up) ?   m_etas[index]   :   m_etas[index+1] ;
			const Spectrum &eta2   = (!is_up) ?   m_etas[index+1] :   m_etas[index]   ;
			const Spectrum &kappa2 = (!is_up) ? m_kappas[index+1] : m_kappas[index]   ;
			const Spectrum &reflectance = m_reflectances[index];

			// get the normal from the normal map, keep in mind that we need to track  the "actual" 
			//	intersection within the texture not just the intersection of the surface (bRec.its)
			//	----> we have the actual position tracked in uvposition
			const Vector &n=m_normals[index];//getNormal(index,uvposition,is_up);
			const Frame &perturbed=getFrame(its,n);
			const Vector &wiperturbed=perturbed.toLocal(its.toWorld(wi));


			const Spectrum &eta    = eta2   / eta1;
			const Spectrum &kappa  = kappa2 / eta1;
		
			#ifdef DEBUG_LAYERINFO		
			printf("Bounce: %d Layer: %d(isup=%d) uv=%s e=%s", bounces, index, is_up,
					bRec.its.uv.toString().c_str(),e.toString().c_str());
			printf("alpha=%f, eta1=%s, eta2=%s, kappa2=%s\n", alpha, eta1.toString().c_str(), 
					eta2.toString().c_str(), kappa2.toString().c_str());
			#endif

			/* Get the next interaction with the layer */
			Vector3 wo;
			if(is_diffuse) // lambertian
			{
				Point2 uv(nextRandom(), nextRandom());

				#ifndef SAMPLEDIRECTION
				Float pdf;
				epre=e;
				e *= interactWithDiffusePDF(wiperturbed, reflectance, uv, wo, pdf);
				#else
				e *= interactWithDiffuse(wiperturbed, reflectance, uv, wo);
				#endif
				//Bench Assert(e.isValid());
				wo = its.toLocal(perturbed.toWorld(wo));

				#ifndef SAMPLEDIRECTION
				// Add vertex to the Path
				m_path.push_back(new PathInfo(PathInfo::EBSDF, bounces, index, is_up, true, wi, wo, pdf, epre, e));
				#endif

				is_up = !is_up;
			}
			else if(kappa.isZero()) // dielectric
			{	
				Vector3 uvw(nextRandom(), nextRandom(), nextRandom());
				bool is_reflected = true;
				#ifndef SAMPLEDIRECTION
				Float pdf;
				epre=e;
				e *= interactWithDielectricPDF(wiperturbed, eta.average(), alphaU, alphaV, uvw, wo, pdf, is_reflected, m_type, mode);
				#else
				e *= interactWithDielectric(wiperturbed, eta.average(), alphaU, alphaV, uvw, wo, is_reflected, m_type, mode);
				#endif
				//Bench Assert(e.isValid());
				wo = its.toLocal(perturbed.toWorld(wo));

				#ifndef SAMPLEDIRECTION
				// Add vertex to the Path
				m_path.push_back(new PathInfo(PathInfo::EBSDF, bounces, index, is_up, is_reflected, wi, wo, pdf, epre, e));
				#endif
				
				if(is_reflected) {
					is_up = !is_up;
				}
			} 
			else // conductor
			{
				Vector2 uv(nextRandom(), nextRandom());
				#ifndef SAMPLEDIRECTION
				Float pdf;
				epre=e;
				e *= interactWithConductorPDF(wiperturbed, eta, kappa, alphaU, alphaV, uv, wo, pdf, m_type) * reflectance;
				#else
				e *= interactWithConductor(wiperturbed, eta, kappa, alphaU, alphaV, uv, wo, m_type) * reflectance;
				#endif
				//Bench Assert(e.isValid());
				wo = its.toLocal(perturbed.toWorld(wo));

				#ifndef SAMPLEDIRECTION
				// Add vertex to the Path
				m_path.push_back(new PathInfo(PathInfo::EBSDF, bounces, index, is_up, true, wi, wo, pdf, epre, e));
				#endif

				is_up = !is_up;
			}

			if(e.isZero() || !e.isValid()) { return false; }

			if(wo.z<0.0) // Avoid numerical errors or change of side due to normal map
				return false;

			/* If the ray continues its travel, updates its parameters */
			wi       = wo;
			bounces += 1;
		}

		/* Update transport information */
		index   += (is_up) ? -1 : 1;

		if(index < 0) {
			// BEGIN COMMENTED TO COMPILE CAREFUL
			/*bRec.wo = wi;
			bRec.eta = 1.0f;
			bRec.sampledComponent = 0;
			bRec.sampledType = EGlossyReflection;
			pdf = this->pdf(bRec, ESolidAngle);*/
			// END COMMENTED TO COMPILE CAREFUL

			// last node of m_path contains:
			// Total path throughput: thruAfter
			// sampled direction: wo
			//Bench Assert(e.isValid());
			#ifdef SAMPLEDIRECTION
			wo_=wi;
			e_=e;
			#endif
			return true;
		} 
		else if(index >= nb_layers) {
			return false;
		}
	}

	return false;
}

// This method connects a full path (this) with a very short path in the reverse direction 
//	the connection is performed in the next layer
//#define DEBUGEVAL
Spectrum EvalPath(const Path &shortpath, const int &nb_layers,
		const std::vector<Float> &m_depths, const std::vector<Float> &m_gs,
		const MicrofacetDistribution::EType &m_type,
		const Spectrum *m_etas, const Spectrum *m_kappas, const Float *m_alphasU, const Float *m_alphasV,
		const Spectrum *m_reflectances, const Vector *m_normals,
		const Spectrum *m_sigmas, const Spectrum *m_sigmaa) 
{
	/*enum EType {
		EInterface = 0,
		EVolume = 1,
		EScattering = 2
	};
	EType type;
	int layer;
	Float depth;*/

	Spectrum Li(0.0);

	const unsigned int &fullsize=m_path.size();
	const unsigned int &shortsize=shortpath.m_path.size();
	if(fullsize==0 || shortsize==0)
	{
		printf("Empty paths? This should never happen...\n");
		return Spectrum(0.0);
	}

	const Intersection &its=bRec.its;

	const PathInfo *qs=shortpath.m_path[0];
	// connection at the top layer
	if(true)
	{
		if(m_path[0]->layer!=0 || qs->layer!=0)
		{
			printf("Bottom layer not supported yet\n");
			return Spectrum(0.0);
		}
		const float &alphaU    = m_alphasU[0];
		const float &alphaV    = m_alphasV[0];
		const Spectrum &eta1   = m_etas[0];
		const Spectrum &eta2   = m_etas[1];
		//const Spectrum &reflectance = m_reflectances[0];

		// get the normal from the normal map, keep in mind that we need to track  the "actual" 
		//	intersection within the texture not just the intersection of the surface (bRec.its)
		//	----> we have the actual position tracked in uvposition
		const Vector &n=m_normals[0];//getNormal(index,uvposition,is_up);
		const Frame &perturbed=getFrame(its,n);
		const Vector &wi=perturbed.toLocal(its.toWorld(m_path[0]->wi));
		const Vector &wo=perturbed.toLocal(its.toWorld(qs->wi));

		const Spectrum &eta    = eta2   / eta1;

		//return evalDielectric(wi, wo, eta.average(), alphaU, alphaV, 
		//			m_type, ERadiance) * reflectance;
		Li += evalDielectric(wi, wo, eta.average(), alphaU, alphaV,
					m_type, ERadiance);// * reflectance;

	}
	if(fullsize<2 || qs->is_reflected)
		return Li;

	// TODO: code only works for wi,wo on the +z side

	// start with fullpath (this) to find a suitable event 
	unsigned int fullindex;
	for(fullindex=1; fullindex < fullsize; fullindex++)
	{
		const PathInfo *qf=m_path[fullindex];
		if(qf->thruBefore.isZero())
			break;

		const int &layer  = qf->layer;
		switch(layer)
		{
			case 0: // reached top interface
			{
				break;

				if(qf->is_up==false)
					break;
				// wi is up but remember that wi is always in the z+, we need to flip it
				Vector qfwi= -qf->wi;

				const float &alphaU  = m_alphasU[0];
				const float &alphaV  = m_alphasV[0];
				const Spectrum &eta1 = m_etas[1];
				const Spectrum &eta2 = m_etas[0];
				//const Spectrum &reflectance = m_reflectances[0];

				// get the normal from the normal map, keep in mind that we need to track  the "actual" 
				//	intersection within the texture not just the intersection of the surface (bRec.its)
				//	----> we have the actual position tracked in uvposition
				const Vector &n=m_normals[0];//getNormal(index,uvposition,is_up);
				const Frame &perturbed=getFrame(its,n);
				const Vector &wi=perturbed.toLocal(its.toWorld(qfwi));
				const Vector &wo=perturbed.toLocal(its.toWorld(qs->wi));

				const Spectrum &eta    = eta2   / eta1;

				Spectrum e = qf->thruBefore;
				//return evalDielectric(wi, wo, eta.average(), alpha, alpha,
				//			m_type, ERadiance) * reflectance;
				e *= evalDielectric(wi, wo, eta.average(), alphaU, alphaV,
							m_type, ERadiance);// * reflectance;

				Li += e;
			}break;
			case 1: // medium
			{
				const Float &depth	= m_depths[layer];
				if(depth<=0.0)
					break;

				const Vector &wi = qf->wi;
				Vector wo = qs->wo;
				Spectrum e = qf->thruBefore * qs->thruAfter;
				const bool &is_up = qf->is_up;

				const Spectrum &sigma_a = m_sigmaa[layer];
				const Spectrum &sigma_s = m_sigmas[layer];
				const Spectrum &sigma_t = (sigma_a+sigma_s);
				const Float    &g       = m_gs[layer];
				
				if(is_up==true)
					wo = -wo;

				e *= HG_Eval(wi,wo,g);

				// remaining depth in wi direction to next interface
				// but we want remaining depth to interface 0
				//Float rmn_depth = qf->depth;
				Float rmn_depth = is_up ? qf->depth : depth - qf->depth;

				/*bool is_wo_up = wo.z >= 0.0;
				bool is_wi_up = wi.z >= 0.0;
				bool chg_dir  = is_wo_up != is_wi_up;*/

				// remaining depth in wo direction
				//rmn_depth = (chg_dir) ? depth - rmn_depth : rmn_depth;
#ifdef DEBUGEVAL
				if(rmn_depth<0.0)
				{
					//printf("chd_dir=%d qf->depth=%f depth=%f rmn_depth=%f wi.z=%f wo.z=%f\n",chg_dir,qf->depth,depth,rmn_depth,wi.z,wo.z);
					printf("qf->depth=%f depth=%f rmn_depth=%f wi.z=%f wo.z=%f\n",qf->depth,depth,rmn_depth,wi.z,wo.z);
					for(unsigned int bla=0;bla<fullindex;bla++)
						cout << m_path[bla]->toString() << endl;
					cout << "Current node*****=" << endl;
					cout << qf->toString() << endl;
				}
#endif
				// exp(-sigma_t * distance)
				e *= (-sigma_t * rmn_depth/fabs(wo.z)).exp();
				
				Li += e;
			}break;
			case 2: // bsdf
			{
				Vector wo = qs->wo;
				Spectrum e;//(0.0);

				const bool &is_up = qf->is_up;

				const float &alphaU    = m_alphasU[layer];
				const float &alphaV    = m_alphasV[layer];
				const bool is_diffuse  = alphaU < 1e-10;
				const Spectrum &eta1   = (!is_up) ?   m_etas[layer]   :   m_etas[layer+1] ;
				const Spectrum &eta2   = (!is_up) ?   m_etas[layer+1] :   m_etas[layer]   ;
				const Spectrum &kappa2 = (!is_up) ? m_kappas[layer+1] : m_kappas[layer]   ;
				const Spectrum &reflectance = m_reflectances[layer];

				const Vector &n=m_normals[layer];//getNormal(index,uvposition,is_up);
				const Frame &perturbed=getFrame(its,n);
				const Vector &wi=perturbed.toLocal(its.toWorld(qf->wi));

				const Spectrum &eta    = eta2   / eta1;
				const Spectrum &kappa  = kappa2 / eta1;
				
				/* Obtain the meidum information */
				const Float &depth = m_depths[layer-1];
				const Spectrum &sigma_a = m_sigmaa[layer-1];
				const Spectrum &sigma_s = m_sigmas[layer-1];
				const Spectrum &sigma_t = (sigma_a+sigma_s);

				// perform connection using evaluation of BSDF(2)
				if(is_diffuse) // lambertian
				{
					// NO MIS BEGIN
					e = qf->thruBefore * qs->thruAfter;
					wo=perturbed.toLocal(its.toWorld(wo));
					e *= evalDiffuse(wi, wo, reflectance);
					if(depth > 0.0f)
					{
						const Float opt_depth = depth/std::abs(Frame::cosTheta(qs->wo));
						if(std::isnan(opt_depth) || std::isinf(opt_depth)) 
						{ 
							return Spectrum(0.0); 
						}
						e *= evalTr(opt_depth, sigma_t);
					}
					// NO MIS END
					/*wo=perturbed.toLocal(its.toWorld(wo));
					Float pdf1 = qs->pdf;
					Float pdf2;
					Spectrum e1 = evalDiffusePDF(wi, wo, reflectance, pdf2);

					Spectrum Tr;
					Float opt_depth = depth/std::abs(Frame::cosTheta(qs->wo));
					if(std::isnan(opt_depth) || std::isinf(opt_depth)) 
						Tr = Spectrum(0.0); 
					else if(depth > 0.0)
						Tr = evalTr(opt_depth, sigma_t);
					else
						Tr = Spectrum(1.0);

					e = qf->thruBefore * qs->thruAfter * e1 * Tr * miWeight(pdf1,pdf2);
					
					if(!qf->thruAfter.isZero())
					{ // using the sample from bottom interface (i.e. stored in qf)
						// qf->wo is always stored in the z+ side
						const Frame perttop=getFrame(its,m_normals[0]);
						const Vector qfwo = perttop.toLocal(its.toWorld(-qf->wo));
						const Vector qswi = perttop.toLocal(its.toWorld(qs->wi));

						pdf2 = qf->pdf;
						Spectrum e2 = evalDielectricPDF(qswi, qfwo, 
							  (m_etas[1]/m_etas[0]).average(),
							  m_alphasU[0], m_alphasV[0], pdf1,
							  m_type, ETransportModes);
							// * m_reflectances[0];
					
						opt_depth = depth/std::abs(Frame::cosTheta(qf->wo));
						if(std::isnan(opt_depth) || std::isinf(opt_depth)) 
							Tr = Spectrum(0.0); 
						else if(depth > 0.0)
							Tr = evalTr(opt_depth, sigma_t);
						else
							Tr = Spectrum(1.0);
						
						// no need to multiply for qs->thruBefore == 1
					 	e += qf->thruAfter * e2 * Tr * miWeight(pdf2,pdf1);
					} */
				}
				else if(kappa.isZero()) // dielectric
				{ 
					// we know wo is_up==false (because qs->is_reflected is false)
					// if wi is_up==false, we do nothing
					// wi is_up==true we keep wi in the z+ but change wo side 
					//	(eta values already consider this above at initialization)
					if(is_up==true)
						wo = -wo;
					
					// BEGIN WARNING CODE FOR TREATING EXPLICITLY LAYERS 3 and 4
					// if the ray is coming from the -z side of the interface do nothing
					//    this is handled in cases 3 and 4
					if(is_up==true)
						break;
					// END WARNING CODE FOR TREATING EXPLICITLY LAYERS 3 and 4


					wo=perturbed.toLocal(its.toWorld(wo));

					Spectrum Tr;
					Float opt_depth;

					// NO MIS Version
					/*e = qf->thruBefore * qs->thruAfter;
					e *= evalDielectric(wi, wo, eta.average(), alphaU, alphaV, 
						m_type, EImportance) * reflectance;
					opt_depth = depth/std::abs(Frame::cosTheta(qs->wo));
					if(std::isnan(opt_depth) || std::isinf(opt_depth)) 
						Tr = Spectrum(0.0); 
					else if(depth > 0.0)
						Tr = evalTr(opt_depth, sigma_t);
					else
						Tr = Spectrum(1.0);
					e *= Tr; */
					// END NO MIS Version

					// MIS
					Float pdf1 = qs->pdf;
					Float pdf2;
					Spectrum e1 = evalDielectricPDF(wi, wo, eta.average(), alphaU, alphaV, pdf2,
						  m_type, EImportance); 
						//* reflectance;
					
					opt_depth = depth/std::abs(Frame::cosTheta(qs->wo));
					if(std::isnan(opt_depth) || std::isinf(opt_depth)) 
						Tr = Spectrum(0.0); 
					else if(depth > 0.0)
						Tr = evalTr(opt_depth, sigma_t);
					else
						Tr = Spectrum(1.0);

					e = qf->thruBefore * qs->thruAfter * e1 * Tr * miWeight(pdf1,pdf2);
				
					if(((is_up && !qf->is_reflected) || (!is_up && qf->is_reflected)) && !qf->thruAfter.isZero())
					{
						// qf->wo is always stored in the z+ side
						const Frame perttop=getFrame(its,m_normals[0]);
						const Vector qfwo = perttop.toLocal(its.toWorld(-qf->wo));
						const Vector qswi = perttop.toLocal(its.toWorld(qs->wi));

						pdf2 = qf->pdf;
						Spectrum e2 = evalDielectricPDF(qswi, qfwo, 
							  (m_etas[1]/m_etas[0]).average(),
							  m_alphasU[0], m_alphasV[0], pdf1,
							  m_type, ETransportModes);
							//* m_reflectances[0];
					
						opt_depth = depth/std::abs(Frame::cosTheta(qf->wo));
						if(std::isnan(opt_depth) || std::isinf(opt_depth)) 
							Tr = Spectrum(0.0); 
						else if(depth > 0.0)
							Tr = evalTr(opt_depth, sigma_t);
						else
							Tr = Spectrum(1.0);
						
						// no need to multiply for qs->thruBefore == 1
					 	e += qf->thruAfter * e2 * Tr * miWeight(pdf2,pdf1);
					} 
				}
				else // conductor
				{ 
					// wo is_up==false (because qs->is_reflected is false)
					// wi is_up SHOULD BE false
					/*if(is_up==true)
					{
						printf("wi is_up\n");
						return Spectrum(0.0);
					}*/
					// NON MIS VERSION
					/*e = qf->thruBefore * qs->thruAfter;
					wo=perturbed.toLocal(its.toWorld(wo));
					e *= evalConductor(wi, wo, eta, kappa, alphaU, alphaV, 
						m_type) * reflectance;
					if(depth > 0.0f)
					{
						const Float opt_depth = depth/std::abs(Frame::cosTheta(qs->wo));
						if(std::isnan(opt_depth) || std::isinf(opt_depth)) 
						{ 
							return Spectrum(0.0); 
						}
						e *= evalTr(opt_depth, sigma_t);
					}*/
					// MIS VERSION
					// Using the sample from top interface (i.e. stored in qs)
					wo=perturbed.toLocal(its.toWorld(wo));
					Float pdf1 = qs->pdf;
					Float pdf2;
					Spectrum e1 = evalConductorPDF(wi, wo, eta, kappa, alphaU, alphaV, pdf2, 
						m_type) * reflectance;

					Spectrum Tr;
					Float opt_depth = depth/std::abs(Frame::cosTheta(qs->wo));
					if(std::isnan(opt_depth) || std::isinf(opt_depth)) 
						Tr = Spectrum(0.0); 
					else if(depth > 0.0)
						Tr = evalTr(opt_depth, sigma_t);
					else
						Tr = Spectrum(1.0);

					e = qf->thruBefore * qs->thruAfter * e1 * Tr * miWeight(pdf1,pdf2);
					
					if(!qf->thruAfter.isZero())
					{ // using the sample from bottom interface (i.e. stored in qf)
						// qf->wo is always stored in the z+ side
						const Frame perttop=getFrame(its,m_normals[0]);
						const Vector qfwo = perttop.toLocal(its.toWorld(-qf->wo));
						const Vector qswi = perttop.toLocal(its.toWorld(qs->wi));

						pdf2 = qf->pdf;
						Spectrum e2 = evalDielectricPDF(qswi, qfwo, 
							  (m_etas[1]/m_etas[0]).average(),
							  m_alphasU[0], m_alphasV[0], pdf1,
							  m_type, ETransportModes);
							//* m_reflectances[0];
					
						opt_depth = depth/std::abs(Frame::cosTheta(qf->wo));
						if(std::isnan(opt_depth) || std::isinf(opt_depth)) 
							Tr = Spectrum(0.0); 
						else if(depth > 0.0)
							Tr = evalTr(opt_depth, sigma_t);
						else
							Tr = Spectrum(1.0);
						
						// no need to multiply for qs->thruBefore == 1
					 	e += qf->thruAfter * e2 * Tr * miWeight(pdf2,pdf1);
					} 
				}

				Li += e;
			}break;
			case 3: // second medium
			{	
				// get above interface parameters
				const float &alphaU    = m_alphasU[2];
				const float &alphaV    = m_alphasV[2];
				const Spectrum &eta1   = m_etas[2];
				const Spectrum &eta2   = m_etas[3];
				const Spectrum &eta    = eta2   / eta1;

				const Vector &n=m_normals[2];//getNormal(index,uvposition,is_up);
				const Frame &perturbed=getFrame(its,n);

				// direction on the other side of the interface
				const Vector &reversewi=perturbed.toLocal(its.toWorld(qs->wo));
				Vector wo;
				Vector uvw(nextRandom(), nextRandom(), nextRandom());
				bool is_reflected;
				Spectrum thruInterface = interactWithDielectric(reversewi, eta.average(), 
						alphaU, alphaV, uvw, wo, is_reflected, m_type, ETransportModes);
				if(is_reflected)
					break;
#ifdef DEBUGEVAL
				Vector woaux=wo;
#endif
				wo = its.toLocal(perturbed.toWorld(wo));
				
				// Tr of medium above the interface above (i.e. between layers 0 and 2)
				Spectrum Tr1;
				Float opt_depth=m_depths[1]/Frame::cosTheta(qs->wo);
				if(std::isnan(opt_depth) || std::isinf(opt_depth)) 
					break; 
				else 
					Tr1 = evalTr(opt_depth, m_sigmaa[1]+m_sigmas[1]);

				// thru interface layer 0 * Tr medium layer 1 * thru interface layer 2
				Spectrum e = qs->thruAfter * Tr1 * thruInterface;
				
				// current medium properties
				const Float &depth	= m_depths[layer];
				if(depth<=0.0)
					break;

				const Vector &wi = qf->wi;
				const bool &is_up = qf->is_up;

				const Spectrum &sigma_a = m_sigmaa[layer];
				const Spectrum &sigma_s = m_sigmas[layer];
				const Spectrum &sigma_t = (sigma_a+sigma_s);
				const Float    &g       = m_gs[layer];
				
				if(is_up==true)
					wo = -wo;
				
				// thru up to current vertex * phase(wi,wo,g)
				e *= qf->thruBefore * HG_Eval(wi,wo,g);
				if(!e.isValid())
				{
#ifdef DEBUGEVAL
					cout << "\nEval case 3 fail:" << endl;
					cout << "\tHG_Eval(" << wi.toString() << "," << wo.toString() << "," << g << ")=" << HG_Eval(wi,wo,g) << endl;
					cout << "\tqs->wo=" << qs->wo.toString() << endl;
					cout << "\treversewi=" << reversewi.toString() << endl;
					cout << "\twoaux=" << woaux.toString() << endl;
					cout << "\tthruBefore=" << qf->thruBefore.toString() << endl;
#endif
					break;
				}

				// remaining depth to next interface in wi direction
				// change it to depth to upper interface
				Float rmn_depth = is_up ? qf->depth : depth - qf->depth;
				//Float rmn_depth = qf->depth;
#ifdef DEBUGEVAL
				if(rmn_depth<0.0)
				{
					//printf("chd_dir=%d qf->depth=%f depth=%f rmn_depth=%f wi.z=%f wo.z=%f\n",chg_dir,qf->depth,depth,rmn_depth,wi.z,wo.z);
					printf("qf->depth=%f depth=%f rmn_depth=%f wi.z=%f wo.z=%f\n",qf->depth,depth,rmn_depth,wi.z,wo.z);
					for(unsigned int bla=0;bla<fullindex;bla++)
						cout << m_path[bla]->toString() << endl;
					cout << "Current node*****=" << endl;
					cout << qf->toString() << endl;
				}
#endif
				/*bool is_wo_up = wo.z >= 0.0;
				bool is_wi_up = wi.z >= 0.0;
				bool chg_dir  = is_wo_up != is_wi_up;

				// remaining depth in wo direction
				rmn_depth = (chg_dir) ? depth - rmn_depth : rmn_depth;*/
				// exp(-sigma_t * distance)
				
				// Tr medium layer 2 from current depth to interface layer 1
				e *= (-sigma_t * rmn_depth/fabs(wo.z)).exp();
			
				// Li+= thru up to current vertex before sampling event 
				// * phase evaluation
				// * Tr medium layer 2 from vertex to interface layer 1	
				// * thru interface layer 1
				// * Tr medium layer 1 from interface layer 2 to layer 0
				// * thru interface layer 0
				Li += e;
			}break;
			case 4: // third interface
			{
				const float &alphaUint2    = m_alphasU[2];
				const float &alphaVint2    = m_alphasV[2];
				const Spectrum &eta1int2   = m_etas[2];
				const Spectrum &eta2int2   = m_etas[3];
				const Spectrum &etaint2    = eta2int2   / eta1int2;

				const Vector &nint2=m_normals[2];//getNormal(index,uvposition,is_up);
				const Frame &perturbedint2=getFrame(its,nint2);

				// direction on the other side of the interface
				const Vector &reversewi=perturbedint2.toLocal(its.toWorld(qs->wo));
				Vector wo;
				Vector uvw(nextRandom(), nextRandom(), nextRandom());
				bool is_reflected;
				Float pdf1;
				Spectrum thruInterface = interactWithDielectricPDF(reversewi, etaint2.average(), 
						alphaUint2, alphaVint2, uvw, wo, pdf1, is_reflected, m_type, ETransportModes);
				//thruInterface *= std::abs(reversewi.z / wo.z);
				if(is_reflected)
					break;
				wo = its.toLocal(perturbedint2.toWorld(wo));
				
				// Tr of medium above the interface above (i.e. between layers 0 and 2)
				Spectrum Tr1;
				Float opt_depth=m_depths[1]/Frame::cosTheta(qs->wo);
				if(std::isnan(opt_depth) || std::isinf(opt_depth)) 
					break; 
				else 
					Tr1 = evalTr(opt_depth, m_sigmaa[1]+m_sigmas[1]);

				
				// Medium above the interface properties
				opt_depth=m_depths[3]/Frame::cosTheta(wo);
				Spectrum Tr3;
				if(std::isnan(opt_depth) || std::isinf(opt_depth)) 
					break; 
				else 
					Tr3 = evalTr(opt_depth, m_sigmaa[3]+m_sigmas[3]);
				// thru interface layer 0 * Tr medium layer 1 * thru interface layer 2 * Tr mlayer 3
				Spectrum e = qs->thruAfter * Tr1 * thruInterface * Tr3;
				
				// now do the actual interface where we're at
				const bool &is_up = qf->is_up;

				const float &alphaU    = m_alphasU[layer];
				const float &alphaV    = m_alphasV[layer];
				const bool is_diffuse  = alphaU < 1e-10;
				const Spectrum &eta1   = (!is_up) ?   m_etas[layer]   :   m_etas[layer+1] ;
				const Spectrum &eta2   = (!is_up) ?   m_etas[layer+1] :   m_etas[layer]   ;
				const Spectrum &kappa2 = (!is_up) ? m_kappas[layer+1] : m_kappas[layer]   ;
				const Spectrum &reflectance = m_reflectances[layer];

				const Vector &n=m_normals[layer];//getNormal(index,uvposition,is_up);
				const Frame &perturbed=getFrame(its,n);
				const Vector &wi=perturbed.toLocal(its.toWorld(qf->wi));
				const Float &depth = m_depths[layer-1];

				const Spectrum &eta    = eta2   / eta1;
				const Spectrum &kappa  = kappa2 / eta1;
				if(is_diffuse) // lambertian
				{
					// NO MIS BEGIN
					wo=perturbed.toLocal(its.toWorld(wo));
					e *= qf->thruBefore * evalDiffuse(wi, wo, reflectance);
				}
				else if(kappa.isZero()) // dielectric
				{ 
					// we know wo is_up==false (because qs->is_reflected is false)
					// if wi is_up==false, we do nothing
					// wi is_up==true we keep wi in the z+ but change wo side 
					//	(eta values already consider this above at initialization)
					if(is_up==true)
						wo = -wo;

					wo=perturbed.toLocal(its.toWorld(wo));

					// NO MIS Version
					//e *= qf->thruBefore * evalDielectric(wi, wo, eta.average(), alphaU, alphaV, 
					//	m_type, EImportance) * reflectance;

					// MIS
					Float pdf2;
					Spectrum e1 = evalDielectricPDF(wi, wo, eta.average(), alphaU, alphaV, pdf2,
						  m_type, EImportance); 
					e *= qf->thruBefore * e1 * miWeight(pdf1,pdf2);
					
					if(((is_up && !qf->is_reflected) || (!is_up && qf->is_reflected)) && !qf->thruAfter.isZero())
					{
						// qf->wo is always stored in the z+ side
						const Vector qfwo = perturbedint2.toLocal(its.toWorld(-qf->wo));

						pdf2 = qf->pdf;
						Spectrum e2 = evalDielectricPDF(reversewi, qfwo, 
							  etaint2.average(),
							  alphaUint2, alphaVint2, pdf1,
							  m_type, ETransportModes) * qs->thruAfter * Tr1;
					
						opt_depth = depth/std::abs(Frame::cosTheta(qf->wo));
						if(std::isnan(opt_depth) || std::isinf(opt_depth)) 
							Tr3 = Spectrum(0.0); 
						else if(depth > 0.0)
							Tr3 = evalTr(opt_depth, m_sigmaa[3]+m_sigmas[3]);
						else
							Tr3 = Spectrum(1.0);
						
						// no need to multiply for qs->thruBefore == 1
					 	e += qf->thruAfter * e2 * Tr3 * miWeight(pdf2,pdf1);
					} 
				}
				else
				{ // conductor
					wo=perturbed.toLocal(its.toWorld(wo));
					// NON MIS VERSION
					//e *= qf->thruBefore * evalConductor(wi, wo, eta, kappa, alphaU, alphaV, 
					//	m_type) * reflectance;
					
					// MIS
					Float pdf2;
					Spectrum e1 = evalConductorPDF(wi, wo, eta, kappa, alphaU, alphaV, pdf2,
						  m_type) * reflectance; 
					e *= qf->thruBefore * e1 * miWeight(pdf1,pdf2);
					
					if(!qf->thruAfter.isZero())
					{
						// qf->wo is always stored in the z+ side
						const Vector qfwo = perturbedint2.toLocal(its.toWorld(-qf->wo));

						pdf2 = qf->pdf;
						Spectrum e2 = evalDielectricPDF(reversewi, qfwo, 
							  etaint2.average(),
							  alphaUint2, alphaVint2, pdf1,
							  m_type, ETransportModes) * qs->thruAfter * Tr1;
					
						opt_depth = depth/std::abs(Frame::cosTheta(qf->wo));
						if(std::isnan(opt_depth) || std::isinf(opt_depth)) 
							Tr3 = Spectrum(0.0); 
						else if(depth > 0.0)
							Tr3 = evalTr(opt_depth, m_sigmaa[3]+m_sigmas[3]);
						else
							Tr3 = Spectrum(1.0);
						
						// no need to multiply for qs->thruBefore == 1
					 	e += qf->thruAfter * e2 * Tr3 * miWeight(pdf2,pdf1);
					} 
				}
				Li+=e;
			}break;
		}
	}

	if(Li.isValid())
	{
		return Li;
		//return clamp(Li);
	}
	return Spectrum(0.0);
}

Float EvalPDFFull(const Path &shortpath, const int &nb_layers,
		const std::vector<Float> &m_depths, const std::vector<Float> &m_gs,
		const MicrofacetDistribution::EType &m_type,
		const Spectrum *m_etas, const Spectrum *m_kappas, const Float *m_alphasU, const Float *m_alphasV,
		const Spectrum *m_reflectances, const Vector *m_normals,
		const Spectrum *m_sigmas, const Spectrum *m_sigmaa) const
{
	/*enum EType {
		EInterface = 0,
		EVolume = 1,
		EScattering = 2
	};
	EType type;
	int layer;
	Float depth;*/

	Float pdf = 0.0;

	const unsigned int &fullsize=m_path.size();
	const unsigned int &shortsize=shortpath.m_path.size();
	if(fullsize==0 || shortsize==0)
	{
		printf("Empty paths? This should never happen...\n");
		return 0.0;
	}

	const Intersection &its=bRec.its;

	const PathInfo *qs=shortpath.m_path[0];
	// connection at the top layer
	Float pdfqs_re;
	if(true)
	{
		if(m_path[0]->layer!=0 || qs->layer!=0)
		{
			printf("Bottom layer not supported yet\n");
			return 0.0;
		}
		const float &alphaU    = m_alphasU[0];
		const float &alphaV    = m_alphasV[0];
		const Spectrum &eta1   = m_etas[0];
		const Spectrum &eta2   = m_etas[1];

		const Vector &n=m_normals[0];//getNormal(index,uvposition,is_up);
		const Frame &perturbed=getFrame(its,n);
		const Vector &wi=perturbed.toLocal(its.toWorld(m_path[0]->wi));
		const Vector &wo=perturbed.toLocal(its.toWorld(qs->wi));

		const Spectrum &eta    = eta2   / eta1;

		//return evalDielectric(wi, wo, eta.average(), alphaU, alphaV, 
		//			m_type, ERadiance) * reflectance;
		pdf += DielectricPDF(wi, wo, eta.average(), alphaU, alphaV, m_type);

		const Vector &wirevqs = perturbed.toLocal(its.toWorld(-qs->wo));
		//pdfqs_re = DielectricPDF(-qs->wo, qs->wi, eta.average(), alphaU, alphaV, m_type);
		pdfqs_re = DielectricPDF(wirevqs, wo, eta.average(), alphaU, alphaV, m_type);
		/*printf("Begin ------\n");
		//cout << toString() << endl;
		printf("toplayerPDF=%f\n",pdf);
		//cout << "Reverse -qs->wo=" << wirevqs.toString() << " qs->wi=" << wo.toString() << endl;*/

	}
	if(fullsize<2 || qs->is_reflected)
	{
		//printf("END---------------pdf=%f\n\n\n",pdf);
		return pdf;
	}

	// TODO: code only works for wi,wo on the +z side

	// start with fullpath (this) to find a suitable event 
	unsigned int fullindex;
	//Float pdfaccum=m_path[0]->pdf;
	Float pdfqs=pdfqs_re / qs->pdf;
	for(fullindex=1; fullindex < fullsize; fullindex++)
	{
		const PathInfo *qf=m_path[fullindex];
		if(qf->thruBefore.isZero())
			break;

		const int &layer  = qf->layer;
		switch(layer)
		{
			case 0: // reached top interface
			{
				break;
				if(qf->is_up==false)
					break;
				// wi is up but remember that wi is always in the z+, we need to flip it
				Vector qfwi= -qf->wi;

				const float &alphaU     = m_alphasU[0];
				const float &alphaV     = m_alphasV[0];
				const Spectrum &eta1   = m_etas[1];
				const Spectrum &eta2   = m_etas[0];

				// get the normal from the normal map, keep in mind that we need to track  the "actual" 
				//	intersection within the texture not just the intersection of the surface (bRec.its)
				//	----> we have the actual position tracked in uvposition
				const Vector &n=m_normals[0];//getNormal(index,uvposition,is_up);
				const Frame &perturbed=getFrame(its,n);
				const Vector &wi=perturbed.toLocal(its.toWorld(qfwi));
				const Vector &wo=perturbed.toLocal(its.toWorld(qs->wi));

				const Spectrum &eta    = eta2   / eta1;

				//return evalDielectric(wi, wo, eta.average(), alphaU, alphaV, 
				//			m_type, ERadiance) * reflectance;
				//pdf += pdfaccum*DielectricPDF(wi, wo, eta.average(), alphaU, alphaV, m_type)*pdfqs;
				pdf += DielectricPDF(wi, wo, eta.average(), alphaU, alphaV, m_type);
				//pdf += DielectricPDF(wi, wo, eta.average(), alphaU, alphaV, m_type);
			}break;
			case 2: // bsdf: change to 1 for no medium, 2 for medium
			{
				Vector wo = qs->wo;

				const bool &is_up = qf->is_up;

				const Float &alphaU    = m_alphasU[layer];
				const Float &alphaV    = m_alphasV[layer];
				const bool is_diffuse  = alphaU < 1e-10;
				const Spectrum &eta1   = (!is_up) ?   m_etas[layer]   :   m_etas[layer+1] ;
				const Spectrum &eta2   = (!is_up) ?   m_etas[layer+1] :   m_etas[layer]   ;
				const Spectrum &kappa2 = (!is_up) ? m_kappas[layer+1] : m_kappas[layer]   ;
				//const Spectrum &sigma_t= m_sigmaa[layer] + m_sigmas[layer];
				const Spectrum &sigma_t= m_sigmaa[layer-1] + m_sigmas[layer-1];

				const Vector &n=m_normals[layer];//getNormal(index,uvposition,is_up);
				const Frame &perturbed=getFrame(its,n);
				const Vector &wi=perturbed.toLocal(its.toWorld(qf->wi));

				const Spectrum &eta    = eta2   / eta1;
				const Spectrum &kappa  = kappa2 / eta1;
				
				// perform connection using evaluation of BSDF(2)
				if(is_diffuse)
				{
					if(is_up==true)
					{
						printf("wi is_up\n");
						return 0.0;
					}
					wo=perturbed.toLocal(its.toWorld(wo));
					const Float pdfFailure = MediumPdfFailure(sigma_t, m_depths[layer-1], wo);
					pdf += DiffusePDF(wi, wo) * pdfFailure * pdfqs;
				}
				else if(kappa.isZero()) // dielectric
				{ 
					// we know wo is_up==false (because qs->is_reflected is false)
					// if wi is_up==false, we do nothing
					// wi is_up==true we keep wi in the z+ but change wo side 
					//	(eta values already consider this above at initialization)
					if(is_up==true)
						wo = -wo;

					wo=perturbed.toLocal(its.toWorld(wo));
					//const Float pdfFailure = MediumPdfFailure(sigma_t, m_depths[layer], wo);
					const Float pdfFailure = MediumPdfFailure(sigma_t, m_depths[layer-1], wo);

					//pdf += pdfaccum * DielectricPDF(wi, wo, eta.average(), alphaU, alphaV, m_type)
					//	* pdfqs;
					Float bsdfpdf = DielectricPDF(wi, wo, eta.average(), alphaU, alphaV, m_type);
					// NON MIS
					//pdf +=  bsdfpdf * pdfFailure * pdfqs;
					// NON MIS END

					// MIS BEGIN
					pdf +=  bsdfpdf	* pdfFailure * pdfqs * miWeight(qs->pdf,bsdfpdf);

					if((!is_up && qf->is_reflected) && !qf->thruAfter.isZero())
					{
						// qf->wo is always stored in the z+ side
						const Frame perttop=getFrame(its,m_normals[0]);
						const Vector qfwo = perttop.toLocal(its.toWorld(-qf->wo));
						const Vector qswi = perttop.toLocal(its.toWorld(qs->wi));

						Float pdfrefract = DielectricPDF(qswi, qfwo, 
							  (m_etas[1]/m_etas[0]).average(),
							  m_alphasU[0], m_alphasV[0], m_type);
							//* m_reflectances[0];
					
						Float pdfrefract_rev = DielectricPDF(qfwo, qswi, 
							  (m_etas[1]/m_etas[0]).average(),
							  m_alphasU[0], m_alphasV[0], m_type);

						const Float pdfFailure = MediumPdfFailure(sigma_t, m_depths[layer-1], qfwo);

						pdf += pdfrefract_rev * pdfFailure * miWeight(qf->pdf, pdfrefract);

					} 
					// MIS END
					//pdf += DielectricPDF(wi, wo, eta.average(), alphaU, alphaV, m_type)*pdfqs;
				}
				else // conductor
				{ 
					// wo is_up==false (because qs->is_reflected is false)
					// wi is_up SHOULD BE false
					if(is_up==true)
					{
						printf("wi is_up\n");
						return 0.0;
					}
					wo=perturbed.toLocal(its.toWorld(wo));
					//const Float pdfFailure = MediumPdfFailure(sigma_t, m_depths[layer], wo);
					const Float pdfFailure = MediumPdfFailure(sigma_t, m_depths[layer-1], wo);
					//pdf += pdfaccum*ConductorPDF(wi, wo, alphaU, alphaV, m_type)*pdfqs;
					pdf += ConductorPDF(wi, wo, alphaU, alphaV, m_type)
						* pdfFailure * pdfqs;
					/*cout << "bsdfPdf=" << ConductorPDF(wi, wo, alphaU, alphaV, m_type) << endl;
					cout << "pdfFailure=" << pdfFailure << endl;
					cout << "pdfqs=" << pdfqs << endl;
					cout << "incrementtotal=" << ConductorPDF(wi, wo, alphaU, alphaV, m_type) * pdfFailure * pdfqs << endl;*/
					//pdf += ConductorPDF(wi, wo, alphaU, alphaV)*pdfqs;
				}
			}break;
		}
		//cout << "pdfvalue=" << pdf << endl << endl;
		//pdfaccum *= qf->pdf;
	}

	//printf("END---------------pdf=%f\n\n\n",pdf);
	return pdf;
}

std::string toString() const
{
	std::ostringstream oss;
	oss << "Path[" << endl;
	for(unsigned int i=0;i<m_path.size();i++)
	{
		oss << "[" << i << "] " << m_path[i]->toString() << endl;
	}
	oss << "]";
	return oss.str();
}

private:
inline Float miWeight(Float pdfA, Float pdfB) const {
	pdfA *= pdfA; pdfB *= pdfB;
	return pdfA / (pdfA + pdfB);
}
			    
inline Spectrum clamp(Spectrum &s, const Float &max_value = 3.0f) const
{
	Float r, g, b;

	s.toLinearRGB(r, g, b);

	r = std::min(r, max_value);
	g = std::min(g, max_value);
	b = std::min(b, max_value);
	
	s.fromLinearRGB(r, g, b);
	return s;
}

inline Float nextRandom() 
{
	//return const_cast<Random*>(rng.get())->nextFloat();
	//return rng.get()->nextFloat();
	return bRec.sampler->next1D();
}
//ref<Random> rng;

std::vector<PathInfo *> m_path;

BSDFSamplingRecord &bRec;
int max_bounces;


/*MicrofacetDistribution::EType m_type;
Float *m_alphas;
Spectrum *m_kappas, *m_etas;
Spectrum *m_reflectances;
Spectrum *m_sigmas;
Spectrum *m_sigmaa;
Vector *m_normals;

int nb_layers;
std::vector<Float>    m_depths;
std::vector<Float>    m_gs;*/
};


MTS_NAMESPACE_END

#endif

