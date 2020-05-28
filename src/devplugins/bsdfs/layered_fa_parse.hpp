// Mitsuba includes
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/logger.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/hw/basicshader.h>

// STL includes
#include <fstream>
#include <sstream>

/* Parsing code is heavily modified and extended from Belcour2018 implementation */

MTS_NAMESPACE_BEGIN

void parseLayers(const Properties &props,
                 int& nb_layers,
                 std::vector<ref<const Texture> >& m_tex_etas,
                 std::vector<ref<const Texture> >& m_tex_kappas,
                 std::vector<Spectrum>& m_eta_factors,
                 std::vector<ref<const Texture> >& m_tex_alphasU,
                 std::vector<ref<const Texture> >& m_tex_alphasV,
                 std::vector<ref<const Texture> >& m_tex_reflectances,
                 std::vector<ref<const Texture2D> >& m_tex_normalmaps,
		 std::vector<Float>& m_tiles,
                 std::vector<Float>& m_depths,
		 std::vector<ref<const Texture> >& m_tex_sigmas,
		 std::vector<ref<const Texture> >& m_tex_sigmaa,
		 std::vector<Float>& m_sigmas_factors,
		 std::vector<Float>& m_sigmaa_factors,
//                 std::vector<Spectrum>& m_sigma_s,
//                 std::vector<Spectrum>& m_sigma_a,
                 std::vector<Float>&    m_gs) {
    /* Parse a layered structure */
    nb_layers = props.getInteger("nb_layers", 0);
    //AssertEx(nb_layers > 0, "layered must have at least one layer");

    // Add the external IOR
    Float extEta = lookupIOR(props, "extEta", "air");
    m_tex_etas.push_back(new ConstantSpectrumTexture(Spectrum(extEta)));
    m_tex_kappas.push_back(new ConstantSpectrumTexture(Spectrum(0.0)));

    // Add the layers IOR and interfaces
    for(int k=0; k<nb_layers; ++k) {
        std::string index = std::to_string(k);
        std::string name;

        // Adding a rough surface parameters
	Spectrum eta_k;
	Spectrum kappa_k;
	name = std::string("material_") + index;
	std::string materialName = props.getString(name, "");
	if(materialName=="")
	{
        	name = std::string("eta_") + index;
	        eta_k = props.getSpectrum(name, Spectrum(1.0f));

	        name = std::string("kappa_") + index;
	        kappa_k = props.getSpectrum(name, Spectrum(0.0f));
	}
	else
	{
		ref<FileResolver> fResolver = Thread::getThread()->getFileResolver();
		eta_k.fromContinuousSpectrum(InterpolatedSpectrum(
				fResolver->resolve("data/ior/" + materialName + ".eta.spd")));
		kappa_k.fromContinuousSpectrum(InterpolatedSpectrum(
		                fResolver->resolve("data/ior/" + materialName + ".k.spd")));
	}

        name = std::string("etafactor_") + index;
        Spectrum etafactor_k = props.getSpectrum(name, Spectrum(1.0f));

        name = std::string("alpha_") + index;
        Float alpha_k = props.getFloat(name, 0.0f);
	Float alphaU_k, alphaV_k;
	if(alpha_k!=0.0)
	{
		alphaU_k = alphaV_k = alpha_k;
	}
	else
	{
		std::string name1 = std::string("alphaU_") + index;
		std::string name2 = std::string("alphaV_") + index;
		alphaU_k = props.getFloat(name1, 0.1f);
		alphaV_k = props.getFloat(name2, 0.1f);
	}
        
	name = std::string("specularReflectance_") + index;
        Spectrum Reflectance = props.getSpectrum(name, Spectrum(1.0f));

	name = std::string("diffuse_") + index;
	if(props.hasProperty(name))
	{
		if(props.hasProperty(std::string("alpha_") + index) || 
				props.hasProperty(std::string("alphaU_") + index))
		{
			SLog(EError, "Layer %d: diffuse surfaces have no roughness values", k);
		}
		else
			alphaU_k = alphaV_k = 0.0; // remove default values
		Reflectance = props.getSpectrum(name);
	}


        name = std::string("normalmap_") + index;
        std::string normalmap = props.getString(name, "");

        name = std::string("tiles_") + index;
        Float tiles = props.getFloat(name, 1.0f);

        // Adding the participating media parameters. This should overloads the
        // already defined parameters in computation. This has to be tested
        // using the depth > 0.0 test.
        name = std::string("depth_") + index;
        Float depth = props.getFloat(name, 0.0f);

        std::string name1 = std::string("sigmas_") + index;
        std::string name2 = std::string("sigmaa_") + index;
	Spectrum sigma_s, sigma_a;
	if(props.hasProperty(name1) || props.hasProperty(name2))
	{
		sigma_s = props.getSpectrum(name1, Spectrum(0.0f));
		sigma_a = props.getSpectrum(name2, Spectrum(0.0f));
	}
	else
	{
		// Match default Guo2018 paremeters
		name = std::string("mediumalbedo_") + index;
	        Spectrum mediumalbedo = props.getSpectrum(name, Spectrum(1.0));

	        name = std::string("sigmat_") + index;
	        Spectrum sigma_t = props.getSpectrum(name, Spectrum(1.0));

		sigma_s = mediumalbedo * sigma_t;
		sigma_a = sigma_t - sigma_s; // materials.h:157
	}
        
	name = std::string("sigmasfactor_") + index;
	Float sigmas_factor = props.getFloat(name, 1.0f);

	name = std::string("sigmaafactor_") + index;
	Float sigmaa_factor = props.getFloat(name, 1.0f);

        name = std::string("g_") + index;
        Float g = props.getFloat(name, 0.9f);

        // IORs when a volume is set
        if(depth > 0.0f) {
            eta_k   = m_tex_etas.back()->getAverage();
            kappa_k = m_tex_kappas.back()->getAverage();
        }

        // Update roughness
        m_tex_etas.push_back(new ConstantSpectrumTexture(eta_k));
        m_tex_kappas.push_back(new ConstantSpectrumTexture(kappa_k));
	m_eta_factors.push_back(etafactor_k);
        m_tex_alphasU.push_back(new ConstantFloatTexture(alphaU_k));
        m_tex_alphasV.push_back(new ConstantFloatTexture(alphaV_k));
        m_tex_reflectances.push_back(new ConstantSpectrumTexture(Reflectance));

	// Update normalmap
	if(normalmap=="")
	{
		m_tex_normalmaps.push_back(NULL);
	}
	else
	{
		Properties normalmapprops("bitmap");
                normalmapprops.setFloat("gamma",1.0);
                normalmapprops.setString("filename",Thread::getThread()->
							getFileResolver()->resolve(normalmap).string());
                normalmapprops.setBoolean("cache",false);
                normalmapprops.setString("filterType","nearest");
                normalmapprops.setString("wrapMode","clamp");
		// load the bitmap texture
		m_tex_normalmaps.push_back(static_cast<Texture2D *>
						(PluginManager::getInstance()->createObject(normalmapprops)));
	}
	m_tiles.push_back(tiles);

        // Update the media
        m_depths.push_back(depth);
        //m_sigma_s.push_back(sigma_s);
        //m_sigma_a.push_back(sigma_a);
	m_tex_sigmas.push_back(new ConstantSpectrumTexture(sigma_s));
	m_tex_sigmaa.push_back(new ConstantSpectrumTexture(sigma_a));
	m_sigmas_factors.push_back(sigmas_factor);
	m_sigmaa_factors.push_back(sigmaa_factor);
        m_gs.push_back(g);

    }

}

MTS_NAMESPACE_END
