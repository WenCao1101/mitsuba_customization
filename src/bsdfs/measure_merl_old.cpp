
#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/hw/basicshader.h>
#include "microfacet.h"
#include "ior.h"

#define POWITACQ_IMPLEMENTATION
#include "powitacq_rgb.h"

#define DJ_BRDF_IMPLEMENTATION 1
#include "dj_brdf.h"

MTS_NAMESPACE_BEGIN


class measure_merl : public BSDF {
public:
	measure_merl(const Properties &props)
		: BSDF(props) {

        fs::path filename = Thread::getThread()->getFileResolver()->resolve(
            props.getString("filename"));
        brdf_merl = new djb::merl(filename.string().c_str());
	//	m_brdf = new powitacq_rgb::measured(filename.string().c_str());
	//	powitacq_rgb::Vector3f tt=m_brdf->get_diffuse();
	//	m_albedo =Color3(tt[0], tt[1], tt[2]);
	//	Float m_diffuse=(m_albedo[0]+m_albedo[1]+m_albedo[2])/3.0f;
	//	m_alpha = m_brdf->get_roughness();
	//	m_specularSamplingWeight=m_specular;
	//	m_albedoReflectance = new ConstantSpectrumTexture(
     //       props.getSpectrum("diffuseReflectance", Spectrum(m_diffuse)));
     //   m_specularReflectance = new ConstantSpectrumTexture(
      //      props.getSpectrum("specularReflectance", Spectrum(m_specular)));
      //  m_params =djb::microfacet::params::isotropic(m_alpha);
	   m_params = djb::tabular::fit_ggx_parameters(djb::tabular(*brdf_merl, 90, false));
	    m_ggx = new djb::ggx();		
		
	}

	measure_merl(Stream *stream, InstanceManager *manager)
		: BSDF(stream, manager) {

		configure();
	}

	~measure_merl()
	{
	//	delete m_brdf;
		delete m_ggx;
		delete brdf_merl;
          
	}

     void configure() {
		/* Verify the input parameter and fix them if necessary */
		m_components.clear();	
		m_components.push_back(EDiffuseReflection | EFrontSide | 0);
		m_usesRayDifferentials = false;
		BSDF::configure();
	}


	Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		if (!(bRec.typeMask & EDiffuseReflection) || measure != ESolidAngle
			|| Frame::cosTheta(bRec.wi) <= 0
			|| Frame::cosTheta(bRec.wo) <= 0)
			return Spectrum(0.0f);
		Color3 result(0.0f);	
	//	powitacq_rgb::Vector2f uu = sample_reverse(bRec);
		
    //    Point2 uum=warp::uniformDiskToSquareConcentric(Point2(bRec.wo.x, bRec.wo.y));
     //   powitacq_rgb::Vector2f uu=powitacq_rgb::Vector2f(uum.x, uum.y); 
	 /*
        powitacq_rgb::Vector3f wi = powitacq_rgb::Vector3f(bRec.wi.x, bRec.wi.y, bRec.wi.z);
        powitacq_rgb::Vector3f wo = powitacq_rgb::Vector3f(bRec.wo.x, bRec.wo.y, bRec.wo.z);
        powitacq_rgb::Vector3f fr = m_brdf->eval(wi, wo,uu);
        result = Color3(fr[0], fr[1], fr[2]);*/
	  /*
      Vector H =normalize(bRec.wi+bRec.wo);
     Float factor1 = 1.0f / (4.0f * M_PI * m_alpha * m_alpha *
                        std::sqrt(Frame::cosTheta(bRec.wi)*Frame::cosTheta(bRec.wo)));
	  Float tan_thetaH=Frame::tanTheta(H);

       Float exponent =-(tan_thetaH* tan_thetaH)/(m_alpha*m_alpha);
	  Float  specRef = factor1 * math::fastexp(exponent);
	     result=m_specularReflectance->eval(bRec.its) * specRef;*/

	//	result = result +m_albedo * INV_PI;
       
	   djb::vec3 o(bRec.wi.x, bRec.wi.y, bRec.wi.z);
		djb::vec3 i(bRec.wo.x, bRec.wo.y, bRec.wo.z);
		djb::vec3 fr_p = brdf_merl->evalp(i, o);
		result = Color3(fr_p.x, fr_p.y, fr_p.z);

		//return result* Frame::cosTheta(bRec.wo);
		return result;
	}
    
    Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        if (!(bRec.typeMask & EDiffuseReflection) || measure != ESolidAngle
			|| Frame::cosTheta(bRec.wi) <= 0
			|| Frame::cosTheta(bRec.wo) <= 0)
			return 0.0f;

     

        Float diffuseProb = 0.0f, specProb = 0.0f;

            
           djb::vec3 o(bRec.wi.x, bRec.wi.y, bRec.wi.z);
		   djb::vec3 i(bRec.wo.x, bRec.wo.y, bRec.wo.z);
       
        //    diffuseProb = warp::squareToCosineHemispherePdf(bRec.wo);
            specProb=m_ggx->pdf(i, o, &m_params);
       
            return  diffuseProb;
        
    }
   
    Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
        if (!(bRec.typeMask & EDiffuseReflection) || Frame::cosTheta(bRec.wi) <= 0)
			return Spectrum(0.0f);

     
       
          /* Sample the tabulated microfacet BRDF */
	     	djb::vec3 o = djb::vec3(bRec.wi.x, bRec.wi.y, bRec.wi.z);
                djb::vec3 i = m_ggx->sample(sample.x, sample.y, o, &m_params);
		

            bRec.wo = Vector(i.x, i.y, i.z);
			
            bRec.sampledComponent = 1;
            bRec.sampledType = EGlossyReflection;

            if (Frame::cosTheta(bRec.wo) <= 0.0f)
                return Spectrum(0.0f);
        
        bRec.eta = 1.0f;
 
       djb::vec3 fr_p = brdf_merl->evalp(i, o) / pdf(bRec, ESolidAngle);
		return Color3(fr_p.x, fr_p.y, fr_p.z);
			
    }
    inline powitacq_rgb::Vector2f sample_reverse(const BSDFSamplingRecord &bRec) const {
        Vector H =normalize(bRec.wi+bRec.wo);
		Point2 sample;
		sample.x = 0.5f;
		sample.y = 0.5f;
		

		powitacq_rgb::Vector2f uu = powitacq_rgb::Vector2f(sample.x,sample.y);

        return uu;
		}
    Spectrum sample(BSDFSamplingRecord &bRec,Float &pdf_,const Point2 &sample_) const {
        Spectrum res = sample(bRec, sample_);
		pdf_ = pdf(bRec, ESolidAngle);
        return res;
    }

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture))
				&& (name == "reflectance" || name == "albedoReflectance")) {

		} else {
			BSDF::addChild(name, child);
		}
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "measure_merl[" << endl
			<< "  id = \"" << getID() << "\"," << endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
private:
	ref<Texture> m_albedoReflectance;
    ref<Texture> m_specularReflectance;
   // std::vector<float> m_brdf;
//	Color3  m_albedo;
//	Float  m_alpha;
//	Float m_specular;
//	 powitacq_rgb::measured *m_brdf;
//	Float m_specularSamplingWeight;
	djb::ggx *m_ggx;
	djb::microfacet::params m_params;
	djb::brdf* brdf_merl;
};

// ================ Hardware shader implementation ================

class measure_merl_shader : public Shader {
public:
	measure_merl_shader(Renderer *renderer, const Texture *m_albedo)
		: Shader(renderer, EBSDFShader), m_albedoReflectance(m_albedo) {
		m_reflectanceShader = renderer->registerShaderForResource(m_albedoReflectance.get());
	}

	bool isComplete() const {
		return m_reflectanceShader.get() != NULL;
	}

	void cleanup(Renderer *renderer) {
		renderer->unregisterShaderForResource(m_albedoReflectance.get());
	}

	void putDependencies(std::vector<Shader *> &deps) {
		deps.push_back(m_reflectanceShader.get());
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    if (cosTheta(wi) < 0.0 || cosTheta(wo) < 0.0)" << endl
			<< "    	return vec3(0.0);" << endl
			<< "    return " << depNames[0] << "(uv) * inv_pi * cosTheta(wo);" << endl
			<< "}" << endl
			<< endl
			<< "vec3 " << evalName << "_albedo(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    return " << evalName << "(uv, wi, wo);" << endl
			<< "}" << endl;
	}

	MTS_DECLARE_CLASS()
private:
	ref<const Texture> m_albedoReflectance;
	ref<Shader> m_reflectanceShader;
	
};

Shader *measure_merl::createShader(Renderer *renderer) const {
	return new measure_merl_shader(renderer,m_albedoReflectance.get());
}

MTS_IMPLEMENT_CLASS(measure_merl_shader, false, Shader)
MTS_IMPLEMENT_CLASS_S(measure_merl, false, BSDF)
MTS_EXPORT_PLUGIN(measure_merl, "measure_merl BRDF")
MTS_NAMESPACE_END
