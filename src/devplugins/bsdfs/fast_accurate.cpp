// Mitsuba includes
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/hw/basicshader.h>

// STL includes
//#include <random>

// Local includes
#include "microfacet.h"
#include "ior.h"

// Local includes
// #include "brdf.h"
#include "layered_fa_parse.hpp"
#include "layered_fa_path.hpp"

MTS_NAMESPACE_BEGIN

class FastAccurate : public BSDF {
public:

FastAccurate(const Properties &props) : BSDF(props) {
	ref<FileResolver> fResolver = Thread::getThread()->getFileResolver();

	MicrofacetDistribution distr(props);
	m_type = distr.getType();

	max_bounces = props.getInteger("max_bounces", 20);

	/* Parse a layered structure */
	parseLayers(props, nb_layers, m_tex_etas, m_tex_kappas, m_etafactors, 
			m_tex_alphasU, m_tex_alphasV, m_tex_reflectances,
			m_tex_normalmaps, m_tiles,
			m_depths, m_tex_sigmas, m_tex_sigmaa, m_sigmasfactors, m_sigmaafactors, m_gs);
	for(int i=0; i<nb_layers; i++)
	{
		m_depthsEmpty.push_back(0.0);
		if(!(m_depths[i]>0.0))
		{
			if(i>0)
			{
				if(!(m_depths[i-1]>0.0))
					Log(EError, "Missing medium between interfaces?");
				// store medium information (previous layer) in the current layer
				m_depthsPDF.push_back(m_depths[i-1]);
			}
			else
				m_depthsPDF.push_back(0.0);
		}
	}
/*	for(int i=0; i<m_depths.size(); i++)
	{
		cout << "depth[" << i << "]=" << m_depths[i] << endl;
	}
	for(int i=0; i<m_depthsEmpty.size(); i++)
	{
		cout << "empty[" << i << "]=" << m_depthsEmpty[i] << endl;
	}
	for(int i=0; i<m_depthsPDF.size(); i++)
	{
		cout << "PDF[" << i << "]=" << m_depthsPDF[i] << endl;
	}*/
}

FastAccurate(Stream *stream, InstanceManager *manager) : BSDF(stream, manager)  {
	Log(EError, "Not implemented");
	configure();
}

void serialize(Stream *stream, InstanceManager *manager) const {
	BSDF::serialize(stream, manager);
	Log(EError, "Not implemented");
}

void configure() {
	unsigned int extraFlags = 0;

	extraFlags |= EAnisotropic;
	extraFlags |= ESpatiallyVarying;

	m_components.clear();
	m_components.push_back(EGlossyReflection | EFrontSide | extraFlags);


        m_usesRayDifferentials = true;

        // Print the resulting microstructure
	Log(EInfo, "---Material composition---");
	Log(EInfo, "Exterior IOR");
	Log(EInfo, " + n = %s", m_tex_etas[0]->getAverage().toString().c_str());
	Log(EInfo, " + k = %s", m_tex_kappas[0]->getAverage().toString().c_str());
	Log(EInfo, "");

	for(int k=0; k<nb_layers; ++k)
	{
		if(m_depths[k] > 0.0f)
		{
			Log(EInfo, "Layer %d [medium]", k);
			Log(EInfo, " + d  = %f", m_depths[k]);
			Log(EInfo, " + ss = %s", m_tex_sigmas[k]->getAverage().toString().c_str());
			Log(EInfo, " + sa = %s", m_tex_sigmaa[k]->getAverage().toString().c_str());
			Log(EInfo, " + g  = %f", m_gs[k]);
			Log(EInfo, "");
		}
		else if(!m_tex_alphasU[k]->getAverage().isZero())
		{
			if(m_tex_kappas[k+1]->getAverage().isZero())
				Log(EInfo, "Layer %d [interface dielectric]", k);
			else
				Log(EInfo, "Layer %d [interface conductor]", k);
			Log(EInfo, " + n = %s", m_tex_etas[k+1]->getAverage().toString().c_str());
			Log(EInfo, " + k = %s", m_tex_kappas[k+1]->getAverage().toString().c_str());
			Log(EInfo, " + a = [u=%f, v=%f]", m_tex_alphasU[k]->getAverage().average(),
							m_tex_alphasV[k]->getAverage().average());
			Log(EInfo, "");
		}
		else
		{
			Log(EInfo, "Layer %d [interface diffuse]", k);
			Log(EInfo, " + albedo = %s", m_tex_reflectances[k]->getAverage().toString().c_str());
			Log(EInfo, "");
		}
	}



        BSDF::configure();
}


/* Fetch the material properties at the intersection point */
inline void evalMaps(const BSDFSamplingRecord &bRec,
			Spectrum* m_etas, Spectrum* m_kappas,
			Float* m_alphasU, Float* m_alphasV, Spectrum* m_reflectances,
			Vector* m_normals, Spectrum* m_sigmas, Spectrum* m_sigmaa) const 
{

	// Fetch air values
	m_etas[0]   = m_tex_etas[0]->eval(bRec.its, false);
	m_kappas[0] = m_tex_kappas[0]->eval(bRec.its, false);

	// Fetch textured values
	for(int i=0; i<nb_layers; ++i) 
	{
		m_etas[i+1]   = m_tex_etas[i+1]->eval(bRec.its, false)*m_etafactors[i];
		m_kappas[i+1] = m_tex_kappas[i+1]->eval(bRec.its, false)*m_etafactors[i];
		m_alphasU[i]   = m_tex_alphasU[i]->eval(bRec.its, false).average();
		m_alphasV[i]   = m_tex_alphasV[i]->eval(bRec.its, false).average();
		m_reflectances[i] = m_tex_reflectances[i]->eval(bRec.its);
		m_sigmas[i]   = m_tex_sigmas[i]->eval(bRec.its, false)*m_sigmasfactors[i];
		m_sigmaa[i]   = m_tex_sigmaa[i]->eval(bRec.its, false)*m_sigmaafactors[i];
		Vector n(0.0,0.0,1.0);
		if(m_tex_normalmaps[i]!=NULL)
		{
			// Tile the normal map
			Float u=bRec.its.uv.x*m_tiles[i];
			u = u - floor(u);
			Float v=bRec.its.uv.y*m_tiles[i];
			v = v - floor(v);
			Point2 uv(u,v);
			m_tex_normalmaps[i]->eval(uv).toLinearRGB(n.x, n.y, n.z);
			for (int j=0; j<3; ++j)
				n[j] = 2 * n[j] - 1;
			// WARNING: forbid normals on the other side
			if(n.z<0.0)
				n.z=0.0;
			n=normalize(n);
		}
		m_normals[i]=n;
	}
}


Spectrum eval(const BSDFSamplingRecord &_bRec, EMeasure measure) const {
	BSDFSamplingRecord bRec(_bRec);
	if (measure != ESolidAngle ||
		Frame::cosTheta(bRec.wi) <= 0 ||
		Frame::cosTheta(bRec.wo) <= 0 ||
		((bRec.component != -1 && bRec.component != 0) ||
		!(bRec.typeMask & EGlossyReflection)))
			return Spectrum(0.0f);
	
	// use integrator pathwsampler
	Assert(bRec.sampler);

	Float m_alphasU[nb_layers];
	Float m_alphasV[nb_layers];
	Spectrum m_reflectances[nb_layers];
	Vector m_normals[nb_layers];
	Spectrum m_kappas[nb_layers+1], m_etas[nb_layers+1];
	Spectrum m_sigmas[nb_layers], m_sigmaa[nb_layers];
	evalMaps(bRec, m_etas, m_kappas, m_alphasU, m_alphasV, m_reflectances, m_normals, m_sigmas, m_sigmaa);


	Path forwardpath(bRec, max_bounces);
	forwardpath.generatePath(nb_layers, m_depths, m_gs, m_type, 
		m_etas, m_kappas, m_alphasU, m_alphasV,
		m_reflectances, m_normals, m_sigmas, m_sigmaa, EImportance);

	BSDFSamplingRecord bRecback(bRec);
	bRecback.wi=bRec.wo;
	Path backwardpath(bRecback, 1);
	backwardpath.generatePath(nb_layers, m_depths, m_gs, m_type, 
		m_etas, m_kappas, m_alphasU, m_alphasV,
		m_reflectances, m_normals, m_sigmas, m_sigmaa, ETransportModes);

	return forwardpath.EvalPath(backwardpath, nb_layers, m_depths, m_gs, m_type, 
				m_etas, m_kappas, m_alphasU, m_alphasV,
	                        m_reflectances, m_normals, m_sigmas, m_sigmaa);
}

Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
	if (Frame::cosTheta(bRec.wi) < 0 ||
		((bRec.component != -1 && bRec.component != 0) ||
		!(bRec.typeMask & EGlossyReflection)))
			return Spectrum(0.0f);

	Spectrum e;

	Float m_alphasU[nb_layers];
	Float m_alphasV[nb_layers];
	Spectrum m_reflectances[nb_layers];
	Vector m_normals[nb_layers];
	Spectrum m_kappas[nb_layers+1], m_etas[nb_layers+1];
	Spectrum m_sigmas[nb_layers], m_sigmaa[nb_layers];
	evalMaps(bRec, m_etas, m_kappas, m_alphasU, m_alphasV, m_reflectances, m_normals, m_sigmas, m_sigmaa);
	
	Path path(bRec, max_bounces);
	//bool success = path.generatePath(nb_layers, m_depths, m_gs, m_type, 
	bool success = path.sampleDirection(e, bRec.wo , nb_layers, m_depths, m_gs, m_type, 
				m_etas, m_kappas, m_alphasU, m_alphasV,
	                        m_reflectances, m_normals, m_sigmas, m_sigmaa);
	if(success)
	{
		bRec.eta = 1.0f;
		bRec.sampledComponent = 0;
		bRec.sampledType = EGlossyReflection;

		//e = path.getPathSample(bRec.wo);
		if(e.isValid())
			return e;
	}
	return Spectrum(0.0f);
}

Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
	if (Frame::cosTheta(bRec.wi) < 0 ||
		((bRec.component != -1 && bRec.component != 0) ||
		!(bRec.typeMask & EGlossyReflection)))
			return Spectrum(0.0f);

	Spectrum value = this->sample(bRec, sample);
	pdf = this->pdf(bRec, ESolidAngle);

	return value;
}

Float pdf(const BSDFSamplingRecord &_bRec, EMeasure measure) const 
{
	BSDFSamplingRecord bRec(_bRec);
	if (measure != ESolidAngle ||
		Frame::cosTheta(bRec.wi) <= 0 ||
		Frame::cosTheta(bRec.wo) <= 0 ||
		((bRec.component != -1 && bRec.component != 0) ||
		!(bRec.typeMask & EGlossyReflection)))
		return 0.0f;
	
	// use integrator pathwsampler
	Assert(bRec.sampler);

	return pdfEvalNoScatMedium(bRec) + 0.1;
}

Float pdfEvalNoScatMedium(BSDFSamplingRecord &bRec) const
{
	Float m_alphasU[nb_layers];
	Float m_alphasV[nb_layers];
	Spectrum m_reflectances[nb_layers];
	Vector m_normals[nb_layers];
	Spectrum m_kappas[nb_layers+1], m_etas[nb_layers+1];
	Spectrum m_sigmas[nb_layers], m_sigmaa[nb_layers];

	evalMaps(bRec, m_etas, m_kappas, m_alphasU, m_alphasV, m_reflectances, m_normals, m_sigmas, m_sigmaa);
	// all absorption, no scattering
	for(int i=0;i<nb_layers;i++)
	{
		m_sigmaa[i]+=m_sigmas[i];
		m_sigmas[i]=Spectrum(0.0);
	}

	BSDFSamplingRecord bRecback(bRec);
	bRecback.wi=bRec.wo;

	int numPdfSamples=1;
	Float weight = 1.0/((Float)numPdfSamples);
	Float value = 0.0;
	for(int i=0; i<numPdfSamples; i++)
	{
		Path forwardpath(bRec, nb_layers+1);
		forwardpath.generatePath(nb_layers, m_depths, m_gs, m_type, 
			m_etas, m_kappas, m_alphasU, m_alphasV,
			m_reflectances, m_normals, m_sigmas, m_sigmaa, EImportance);
	
		Path backwardpath(bRecback, 1);
		backwardpath.generatePath(nb_layers, m_depths, m_gs, m_type, 
			m_etas, m_kappas, m_alphasU, m_alphasV,
			m_reflectances, m_normals, m_sigmas, m_sigmaa, ETransportModes);

		value += forwardpath.EvalPDFFull(backwardpath, nb_layers, m_depths, m_gs, m_type, 
				m_etas, m_kappas, m_alphasU, m_alphasV, 
	                        m_reflectances, m_normals, m_sigmas, m_sigmaa) * weight;
	}
	if(std::isnormal(value))
		return value;
	return 0.0;
}



void addChild(const std::string &name, ConfigurableObject *child) {
	if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) 
	{
		int n_layers = nb_layers-1;
		std::string substr = name.substr(0, 3);
		if(substr == "alp")
		{
			std::string substr = name.substr(0,6);
			int index = atoi(name.substr(7).c_str());
			if (substr == "alpha_")
			{
				index = atoi(name.substr(6).c_str());
				m_tex_alphasU[index] = m_tex_alphasV[index] = static_cast<Texture *>(child);
			}
			else if (substr == "alphaU")
			{
				m_tex_alphasU[index] = static_cast<Texture *>(child);
			}
			else if (substr == "alphaV")
			{
				m_tex_alphasV[index] = static_cast<Texture *>(child);
			}
			Log(EInfo, "Adding texture for %s in layer %d/%d", substr.c_str(), index, n_layers);
		}
		else if (substr == "sig") 
		{
			std::string substr = name.substr(0,6);
			int index = atoi(name.substr(7).c_str());
			if(substr == "sigmas")
			{
				m_tex_sigmas[index] = static_cast<Texture *>(child);
			}
			else if(substr == "sigmaa")
			{
				m_tex_sigmaa[index] = static_cast<Texture *>(child);
			}
			Log(EInfo, "Adding texture for %s in layer %d/%d", substr.c_str(), index, n_layers);
	    	}
		else if(substr == "spe")
		{
			std::string substr = name.substr(0,19);
			int index = atoi(name.substr(20).c_str());
			if(substr == "specularReflectance")
			{
				m_tex_reflectances[index] = static_cast<Texture *>(child);
				Log(EInfo, "Adding texture for %s in layer %d/%d", substr.c_str(), index, n_layers);
			}
		}
		else if(substr == "dif")
		{
			std::string substr = name.substr(0,7);
			int index = atoi(name.substr(8).c_str());
			if(substr == "diffuse")
			{
				m_tex_reflectances[index] = static_cast<Texture *>(child);
				Log(EInfo, "Adding texture for %s in layer %d/%d", substr.c_str(), index, n_layers);
				m_tex_alphasU[index]=new ConstantFloatTexture(0.0);
				m_tex_alphasV[index]=new ConstantFloatTexture(0.0);
			}
		}
		else
			BSDF::addChild(name, child);
        } else {
            BSDF::addChild(name, child);
        }
    }

    Float getRoughness(const Intersection &its, int component) const {
	    // TODO: correct this, either set to min/max alpha or avg alpha
	    return 0.1f;
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "FastAccurate[" << endl
            << "  id = \"" << getID() << "\"," << endl
            << "  distribution = " << MicrofacetDistribution::distributionName(m_type) << "," << endl
/*            << "  sampleVisible = " << m_sampleVisible << "," << endl
            << "  alphaU = " << indent(m_alphaU->toString()) << "," << endl
            << "  alphaV = " << indent(m_alphaV->toString()) << "," << endl
            << "  specularReflectance = " << indent(m_specularReflectance->toString()) << "," << endl
            << "  eta = " << m_eta.toString() << "," << endl
            << "  k = " << m_k.toString() << endl*/
            << "]";
        return oss.str();
    }

    Shader *createShader(Renderer *renderer) const;

    MTS_DECLARE_CLASS()
private:
    MicrofacetDistribution::EType m_type;

    /* Layered structure */
    int nb_layers;
    std::vector<ref<const Texture>> m_tex_etas;
    std::vector<ref<const Texture>> m_tex_kappas;
    std::vector<Spectrum> m_etafactors;
    std::vector<ref<const Texture>> m_tex_alphasU;
    std::vector<ref<const Texture>> m_tex_alphasV;
    std::vector<ref<const Texture> > m_tex_reflectances;
    std::vector<ref<const Texture2D> > m_tex_normalmaps;
    std::vector<Float> m_tiles;

    std::vector<Float>    m_depths;
    std::vector<Float>    m_depthsPDF;
    std::vector<Float>    m_depthsEmpty; // full of zeros for pdf computation
    std::vector<ref<const Texture>> m_tex_sigmaa;
    std::vector<Float> m_sigmaafactors;
    std::vector<ref<const Texture>> m_tex_sigmas;
    std::vector<Float> m_sigmasfactors;
    std::vector<Float>    m_gs;

    /* Constants */
    int max_bounces = 10;
};

/**
 * GLSL port of the rough conductor shader. This version is much more
 * approximate -- it only supports the Ashikhmin-Shirley distribution,
 * does everything in RGB, and it uses the Schlick approximation to the
 * Fresnel reflectance of conductors. When the roughness is lower than
 * \alpha < 0.2, the shader clamps it to 0.2 so that it will still perform
 * reasonably well in a VPL-based preview.
 */
class FastAccurateShader : public Shader {
public:
    FastAccurateShader(Renderer *renderer, const Texture *specularReflectance,
            const Texture *alphaU, const Texture *alphaV, const Spectrum &eta,
            const Spectrum &k) : Shader(renderer, EBSDFShader),
            m_specularReflectance(specularReflectance), m_alphaU(alphaU), m_alphaV(alphaV){
        m_specularReflectanceShader = renderer->registerShaderForResource(m_specularReflectance.get());
        m_alphaUShader = renderer->registerShaderForResource(m_alphaU.get());
        m_alphaVShader = renderer->registerShaderForResource(m_alphaV.get());

        /* Compute the reflectance at perpendicular incidence */
        m_R0 = fresnelConductorExact(1.0f, eta, k);
    }

    bool isComplete() const {
        return m_specularReflectanceShader.get() != NULL &&
               m_alphaUShader.get() != NULL &&
               m_alphaVShader.get() != NULL;
    }

    void putDependencies(std::vector<Shader *> &deps) {
        deps.push_back(m_specularReflectanceShader.get());
        deps.push_back(m_alphaUShader.get());
        deps.push_back(m_alphaVShader.get());
    }

    void cleanup(Renderer *renderer) {
        renderer->unregisterShaderForResource(m_specularReflectance.get());
        renderer->unregisterShaderForResource(m_alphaU.get());
        renderer->unregisterShaderForResource(m_alphaV.get());
    }

    void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
        parameterIDs.push_back(program->getParameterID(evalName + "_R0", false));
    }

    void bind(GPUProgram *program, const std::vector<int> &parameterIDs, int &textureUnitOffset) const {
        program->setParameter(parameterIDs[0], m_R0);
    }

    void generateCode(std::ostringstream &oss,
            const std::string &evalName,
            const std::vector<std::string> &depNames) const {
        oss << "uniform vec3 " << evalName << "_R0;" << endl
            << endl
            << "float " << evalName << "_D(vec3 m, float alphaU, float alphaV) {" << endl
            << "    float ct = cosTheta(m), ds = 1-ct*ct;" << endl
            << "    if (ds <= 0.0)" << endl
            << "        return 0.0f;" << endl
            << "    alphaU = 2 / (alphaU * alphaU) - 2;" << endl
            << "    alphaV = 2 / (alphaV * alphaV) - 2;" << endl
            << "    float exponent = (alphaU*m.x*m.x + alphaV*m.y*m.y)/ds;" << endl
            << "    return sqrt((alphaU+2) * (alphaV+2)) * 0.15915 * pow(ct, exponent);" << endl
            << "}" << endl
            << endl
            << "float " << evalName << "_G(vec3 m, vec3 wi, vec3 wo) {" << endl
            << "    if ((dot(wi, m) * cosTheta(wi)) <= 0 || " << endl
            << "        (dot(wo, m) * cosTheta(wo)) <= 0)" << endl
            << "        return 0.0;" << endl
            << "    float nDotM = cosTheta(m);" << endl
            << "    return min(1.0, min(" << endl
            << "        abs(2 * nDotM * cosTheta(wo) / dot(wo, m))," << endl
            << "        abs(2 * nDotM * cosTheta(wi) / dot(wi, m))));" << endl
            << "}" << endl
            << endl
            << "vec3 " << evalName << "_schlick(float ct) {" << endl
            << "    float ctSqr = ct*ct, ct5 = ctSqr*ctSqr*ct;" << endl
            << "    return " << evalName << "_R0 + (vec3(1.0) - " << evalName << "_R0) * ct5;" << endl
            << "}" << endl
            << endl
            << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "   if (cosTheta(wi) <= 0 || cosTheta(wo) <= 0)" << endl
            << "    	return vec3(0.0);" << endl
            << "   vec3 H = normalize(wi + wo);" << endl
            << "   vec3 reflectance = " << depNames[0] << "(uv);" << endl
            << "   float alphaU = max(0.2, " << depNames[1] << "(uv).r);" << endl
            << "   float alphaV = max(0.2, " << depNames[2] << "(uv).r);" << endl
            << "   float D = " << evalName << "_D(H, alphaU, alphaV)" << ";" << endl
            << "   float G = " << evalName << "_G(H, wi, wo);" << endl
            << "   vec3 F = " << evalName << "_schlick(1-dot(wi, H));" << endl
            << "   return reflectance * F * (D * G / (4*cosTheta(wi)));" << endl
            << "}" << endl
            << endl
            << "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "    if (cosTheta(wi) < 0.0 || cosTheta(wo) < 0.0)" << endl
            << "    	return vec3(0.0);" << endl
            << "    return " << evalName << "_R0 * inv_pi * inv_pi * cosTheta(wo);"<< endl
            << "}" << endl;
    }
    MTS_DECLARE_CLASS()
private:
    ref<const Texture> m_specularReflectance;
    ref<const Texture> m_alphaU;
    ref<const Texture> m_alphaV;
    ref<Shader> m_specularReflectanceShader;
    ref<Shader> m_alphaUShader;
    ref<Shader> m_alphaVShader;
    Spectrum m_R0;
};

Shader *FastAccurate::createShader(Renderer *renderer) const {
    return new FastAccurateShader(renderer,
        new ConstantSpectrumTexture(Spectrum(1.0)), new ConstantFloatTexture(0.05), new ConstantFloatTexture(0.05), 
	Spectrum(1.0), Spectrum(0.0));
        //m_specularReflectance.get(), m_alphaU.get(), m_alphaV.get(), m_eta, m_k);
}

MTS_IMPLEMENT_CLASS(FastAccurateShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(FastAccurate, false, BSDF)
MTS_EXPORT_PLUGIN(FastAccurate, "Fast Accurate Layered Material");
MTS_NAMESPACE_END
