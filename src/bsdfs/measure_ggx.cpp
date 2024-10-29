
#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/hw/basicshader.h>
#include "microfacet.h"
#include "ior.h"

#define POWITACQ_IMPLEMENTATION
#include "powitacq_rgb.h"


MTS_NAMESPACE_BEGIN


class measure_ggx : public BSDF {
public:
	measure_ggx(const Properties &props)
		: BSDF(props) {

        fs::path filename = Thread::getThread()->getFileResolver()->resolve(
            props.getString("filename"));

       
		m_brdf = new powitacq_rgb::measured(filename.string().c_str());
		powitacq_rgb::Vector3f tt=m_brdf->get_diffuse();
		m_diffuse=(tt[0]+tt[1]+tt[2])/3.0f;
		m_specular=1.0-m_diffuse;
		m_alpha = m_brdf->get_roughness();
		m_specularSamplingWeight=m_specular;
		m_diffuseReflectance = new ConstantSpectrumTexture(
            props.getSpectrum("diffuseReflectance", Spectrum(m_diffuse)));
        m_specularReflectance = new ConstantSpectrumTexture(
            props.getSpectrum("specularReflectance", Spectrum(m_specular)));
		
	}

	measure_ggx(Stream *stream, InstanceManager *manager)
		: BSDF(stream, manager) {

		configure();
	}

	~measure_ggx()
	{
		delete m_brdf;
          
	}

     void configure() {
		/* Verify the input parameter and fix them if necessary */
		m_components.clear();	
		m_components.push_back(EDiffuseReflection | EFrontSide | 1);
		m_components.push_back(EGlossyReflection | EFrontSide | 0);
		m_usesRayDifferentials = false;
		BSDF::configure();
	}


	Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		if (!(bRec.typeMask & EDiffuseReflection) || measure != ESolidAngle
			|| Frame::cosTheta(bRec.wi) <= 0
			|| Frame::cosTheta(bRec.wo) <= 0)
			return Spectrum(0.0f);
		Color3 result(0.0f);	
    	//powitacq_rgb::Vector2f uu = sample_reverse(bRec);
		 
        Point2 uum=warp::uniformDiskToSquareConcentric(Point2(bRec.wo.x, bRec.wo.y));
        powitacq_rgb::Vector2f uu=powitacq_rgb::Vector2f(uum.x, uum.y); 
        powitacq_rgb::Vector3f wi = powitacq_rgb::Vector3f(bRec.wi.x, bRec.wi.y, bRec.wi.z);
        powitacq_rgb::Vector3f wo = powitacq_rgb::Vector3f(bRec.wo.x, bRec.wo.y, bRec.wo.z);
        powitacq_rgb::Vector3f fr = m_brdf->eval(wi, wo,uu);
        result = Color3(fr[0], fr[1], fr[2]);
	  /*
      Vector H =normalize(bRec.wi+bRec.wo);
     Float factor1 = 1.0f / (4.0f * M_PI * m_alpha * m_alpha *
                        std::sqrt(Frame::cosTheta(bRec.wi)*Frame::cosTheta(bRec.wo)));
	  Float tan_thetaH=Frame::tanTheta(H);

       Float exponent =-(tan_thetaH* tan_thetaH)/(m_alpha*m_alpha);
	  Float  specRef = factor1 * math::fastexp(exponent);
	     result=m_specularReflectance->eval(bRec.its) * specRef;*/

		result = result +m_diffuseReflectance->eval(bRec.its) * INV_PI;
		return result* Frame::cosTheta(bRec.wo);
	}
    
    Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        if (Frame::cosTheta(bRec.wi) <= 0 ||
            Frame::cosTheta(bRec.wo) <= 0 || measure != ESolidAngle)
            return 0.0f;

     

        Float diffuseProb = 0.0f, specProb = 0.0f;

            
            Vector H = normalize(bRec.wi+bRec.wo);
			Float factor1=1.0f/(4.0f*M_PI*m_alpha*m_alpha*dot(H,bRec.wi)*std::pow(Frame::cosTheta(H),3));
           // Float factor1 = 1.0f / (4.0f * M_PI * alphaU * alphaV *
            //    dot(H, bRec.wi) * std::pow(Frame::cosTheta(H), 3));
          //  Float factor2 = H.x / alphaU, factor3 = H.y / alphaV;
              Float tan_thetaH=Frame::tanTheta(H);
         //   Float exponent = -(factor2*factor2+factor3*factor3)/(H.z*H.z);
		    Float exponent =-(tan_thetaH* tan_thetaH)/(m_alpha*m_alpha);
            specProb = factor1 * math::fastexp(exponent);
       
            diffuseProb = warp::squareToCosineHemispherePdf(bRec.wo);

       
            return m_specularSamplingWeight * specProb +
                   (1-m_specularSamplingWeight) * diffuseProb;
	 //  return specProb;
	    
        
    }
   
    inline Spectrum sample(BSDFSamplingRecord &bRec, Float &_pdf, const Point2 &_sample) const {
        Point2 sample(_sample);

        bool hasSpecular = (bRec.typeMask & EGlossyReflection)
                && (bRec.component == -1 || bRec.component == 0);
        bool hasDiffuse  = (bRec.typeMask & EDiffuseReflection)
                && (bRec.component == -1 || bRec.component == 1);

        if (!hasSpecular && !hasDiffuse)
            return Spectrum(0.0f);
        Spectrum result(0.0f);
        bool choseSpecular = true;

      //  if (hasDiffuse && hasSpecular) {
		
       //     if (sample.x <= m_specularSamplingWeight) {
               // sample.x /= m_specularSamplingWeight;
	//		   choseSpecular = true;	

     //       } else {
              //  sample.x = (sample.x - m_specularSamplingWeight)
               //     / (1-m_specularSamplingWeight);
        //        choseSpecular = false;
      //      }
    //    }

        if (choseSpecular) {
          //  Float alphaU = m_roughness;
           // Float alphaV = m_roughness;
            Float phiH=sample.y*2.0f*M_PI;
          //  Float phiH = std::atan(alphaV/alphaU
           //     * std::tan(2.0f * M_PI * sample.y));
          //  if (sample.y > 0.5f)
          //      phiH += M_PI;
            Float cosPhiH = std::cos(phiH);
            Float sinPhiH = math::safe_sqrt(1.0f-cosPhiH*cosPhiH);
            Float thetaH = std::atan(m_alpha*math::safe_sqrt(-math::fastlog(sample.x)));
			
         //   Float thetaH = std::atan(math::safe_sqrt(
           //     -math::fastlog(sample.x) / (
           //         (cosPhiH*cosPhiH) / (alphaU*alphaU) +
           //         (sinPhiH*sinPhiH) / (alphaV*alphaV)
           // )));

            //Vector H = sphericalDirection(thetaH, phiH);
			Vector H (std::sin(thetaH)*std::cos(phiH),std::sin(thetaH)*std::sin(phiH),std::cos(thetaH));
            bRec.wo = H * (2.0f * dot(bRec.wi, H)) - bRec.wi;
			
            bRec.sampledComponent = 1;
            bRec.sampledType = EGlossyReflection;

            if (Frame::cosTheta(bRec.wo) <= 0.0f)
                return Spectrum(0.0f);
        } else {
            bRec.wo = warp::squareToCosineHemisphere(sample);
			
            bRec.sampledComponent = 0;
            bRec.sampledType = EDiffuseReflection;
        }
        bRec.eta = 1.0f;

        _pdf = pdf(bRec, ESolidAngle);

        if (_pdf == 0)
            return Spectrum(0.0f);
        else
            return eval(bRec, ESolidAngle) / _pdf;
			
    }
    inline powitacq_rgb::Vector2f sample_reverse(const BSDFSamplingRecord &bRec) const {
        Vector H =normalize(bRec.wi+bRec.wo);
		Point2 sample;
        Float thetaH = std::acos(H.z);
		//Float thetaH=std::asin(.5f * math::safe_sqrt(math::safe_sqrt(H.x) + math::safe_sqrt(H.y) + math::safe_sqrt(H.z - 1.f)))*2.0f;
		Float phiH = std::atan2(H.y, H.x);
		if(phiH>=0.0)
		    phiH=phiH;
		else
			phiH=phiH+2.0f*M_PI;
        sample.y = phiH/(2.0f*M_PI);
		sample.x = math::fastexp(-std::pow((std::tan(thetaH))/m_alpha,2));
		if (sample.x<0.0f)
		{
			sample.x=0.000001f;
		}
		

		powitacq_rgb::Vector2f uu = powitacq_rgb::Vector2f(sample.x,sample.y);

        return uu;
		}
    Spectrum sample(BSDFSamplingRecord &bRec,const Point2 &sample) const {
        Float pdf;
        return measure_ggx::sample(bRec, pdf, sample);
    }

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture))
				&& (name == "reflectance" || name == "diffuseReflectance")) {

		} else {
			BSDF::addChild(name, child);
		}
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "measure_ggx[" << endl
			<< "  id = \"" << getID() << "\"," << endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
private:
	ref<Texture> m_diffuseReflectance;
    ref<Texture> m_specularReflectance;
   // std::vector<float> m_brdf;
	Float  m_diffuse;
	Float  m_alpha;
	Float m_specular;
	 powitacq_rgb::measured *m_brdf;
	Float m_specularSamplingWeight;
};

// ================ Hardware shader implementation ================

class measure_ggx_shader : public Shader {
public:
	measure_ggx_shader(Renderer *renderer, const Texture *m_diffuse)
		: Shader(renderer, EBSDFShader), m_diffuseReflectance(m_diffuse) {
		m_reflectanceShader = renderer->registerShaderForResource(m_diffuseReflectance.get());
	}

	bool isComplete() const {
		return m_reflectanceShader.get() != NULL;
	}

	void cleanup(Renderer *renderer) {
		renderer->unregisterShaderForResource(m_diffuseReflectance.get());
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
			<< "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    return " << evalName << "(uv, wi, wo);" << endl
			<< "}" << endl;
	}

	MTS_DECLARE_CLASS()
private:
	ref<const Texture> m_diffuseReflectance;
	ref<Shader> m_reflectanceShader;
	
};

Shader *measure_ggx::createShader(Renderer *renderer) const {
	return new measure_ggx_shader(renderer,m_diffuseReflectance.get());
}

MTS_IMPLEMENT_CLASS(measure_ggx_shader, false, Shader)
MTS_IMPLEMENT_CLASS_S(measure_ggx, false, BSDF)
MTS_EXPORT_PLUGIN(measure_ggx, "measure_ggx BRDF")
MTS_NAMESPACE_END
