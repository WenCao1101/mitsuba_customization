#include <mitsuba/core/properties.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/tensor.h>

#include <mitsuba/core/warp.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/ior.h>
#include <mitsuba/render/microfacet.h>
#include <mitsuba/render/texture.h>

#define DJ_BRDF_IMPLEMENTATION 1
#include "merl_brdf.h"


NAMESPACE_BEGIN(mitsuba)

/**!

.. _bsdf-measured:

Measured material (:monosp:`measured`)
--------------------------------------

.. pluginparameters::

 * - filename
   - |string|
   - Filename of the material data file to be loaded


*/
template <typename Float, typename Spectrum>
class MERLDEC final : public BSDF<Float, Spectrum> {
public:
    MI_IMPORT_BASE(BSDF, m_flags, m_components)
    MI_IMPORT_TYPES(Texture, MicrofacetDistribution)

    MERLDEC(const Properties &props) : Base(props) {
        
        auto fs            = Thread::thread()->file_resolver();
        fs::path file_path = fs->resolve(props.string("filename"));
        m_name             = file_path.filename().string();
        m_diffuse_reflectance  = props.texture<Texture>("diffuse_reflectance",  .5f);
        if (props.has_property("specular_reflectance"))
            m_specular_reflectance = props.texture<Texture>("specular_reflectance", 1.f);

//////////////////////////////////////////
     mitsuba::MicrofacetDistribution<ScalarFloat, Spectrum> distr(props);
       // m_params= distr.params(); 
     m_type =MicrofacetType::GGX;
     m_sample_visible = true;

    if (distr.is_anisotropic())
            Throw("The 'roughplastic' plugin currently does not support "
                  "anisotropic microfacet distributions!");

      //  m_alpha = distr.alpha();
        
////////////////////////////////////////////////
        m_brdf= new merl<Float>(file_path);
        // load tabulated
       tabular<Float>  *temp=new tabular<Float>(*m_brdf, 90, false);
       m_alpha = tabular<Float>::fit_ggx_parameter(*temp);

     

     //   m_ggx = new ggx<Float>();

//////////////////////////////////////////////////////////////////////////////////////
     m_components.push_back(BSDFFlags::GlossyReflection | BSDFFlags::FrontSide);
     m_components.push_back(BSDFFlags::DiffuseReflection | BSDFFlags::FrontSide);
        m_flags =  m_components[0] | m_components[1];
        dr::set_attr(this, "flags", m_flags);
        parameters_changed();
    }

    void traverse(TraversalCallback *callback) override {
        callback->put_object("diffuse_reflectance", m_diffuse_reflectance.get(), +ParamFlags::Differentiable);
        callback->put_parameter("alpha",            m_alpha,                     ParamFlags::Differentiable | ParamFlags::Discontinuous);
        
        if (m_specular_reflectance)
            callback->put_object("specular_reflectance", m_specular_reflectance.get(), +ParamFlags::Differentiable);
    }
  void parameters_changed(const std::vector<std::string> &keys = {}) override {
        

        /* Compute weights that further steer samples towards
           the specular or diffuse components */
        Float d_mean = m_diffuse_reflectance->mean(),
              s_mean = 1.f;

        if (m_specular_reflectance)
            s_mean = m_specular_reflectance->mean();

         m_specular_sampling_weight = s_mean / (d_mean + s_mean);

        // Precompute rough reflectance (vectorized)
       
    }
    

    std::pair<BSDFSample3f, Spectrum> sample(const BSDFContext &ctx,
                                             const SurfaceInteraction3f &si,
                                             Float sample1,
                                             const Point2f &sample2,
                                             Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);

        bool has_specular = ctx.is_enabled(BSDFFlags::GlossyReflection, 0),
             has_diffuse  = ctx.is_enabled(BSDFFlags::DiffuseReflection, 1);

        Float cos_theta_i = Frame3f::cos_theta(si.wi);
        active &= cos_theta_i > 0.f;

        BSDFSample3f bs = dr::zeros<BSDFSample3f>();
        Spectrum result(0.f);
        if (unlikely((!has_specular && !has_diffuse) || dr::none_or<false>(active)))
            return { bs, result };

      
         // Determine which component should be sampled
        Float prob_specular =  m_specular_sampling_weight,
              prob_diffuse  =  (1.f - m_specular_sampling_weight);

        if (unlikely(has_specular != has_diffuse))
            prob_specular = has_specular ? 1.f : 0.f;
        else
            prob_specular = prob_specular / (prob_specular + prob_diffuse);
        prob_diffuse = 1.f - prob_specular;

        Mask sample_specular = active && (sample1 < prob_specular),
             sample_diffuse = active && !sample_specular;

        bs.eta = 1.f;

        if (dr::any_or<true>(sample_specular)) {
            MicrofacetDistribution distr(m_type, m_alpha, m_sample_visible);
            Normal3f m = std::get<0>(distr.sample(si.wi, sample2));

            dr::masked(bs.wo, sample_specular) = reflect(si.wi, m);
            dr::masked(bs.sampled_component, sample_specular) = 0;
            dr::masked(bs.sampled_type, sample_specular) = +BSDFFlags::GlossyReflection;
        }

        if (dr::any_or<true>(sample_diffuse)) {
            dr::masked(bs.wo, sample_diffuse) = warp::square_to_cosine_hemisphere(sample2);
            dr::masked(bs.sampled_component, sample_diffuse) = 1;
            dr::masked(bs.sampled_type, sample_diffuse) = +BSDFFlags::DiffuseReflection;
        }

        bs.pdf = pdf(ctx, si, bs.wo, active);
        active &= bs.pdf > 0.f;
        result = eval(ctx, si, bs.wo, active);

        return { bs, (UnpolarizedSpectrum(result) / bs.pdf) & active };
    }

     Spectrum eval(const BSDFContext &ctx, const SurfaceInteraction3f &si,
                  const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        Vector3f wi = si.wi;

        active &= Frame3f::cos_theta(wi) > 0.f &&
                  Frame3f::cos_theta(wo) > 0.f;

        if (!ctx.is_enabled(BSDFFlags::GlossyReflection) || dr::none_or<false>(active))
            return Spectrum(0.f);
      //  UnpolarizedSpectrum spec(0.0);
     //   merl_brdf<Float>* merl_brdf;

        Vector3f m = dr::normalize(wo + wi);
    //    djb::vec3 o(float(wi.x()), float(wi.y()), float(wi.z()));
	//	djb::vec3 i(float(wo.x()), float(wo.y()), float(wo.z()));
		Vector3f fr_p = m_brdf->evalp(wi, wo);
      // Color3f fr = Color3f(fr_p.x(), fr_p.y(), fr_p.z());
        return UnpolarizedSpectrum (fr_p) & active;
    }

    Float pdf(const BSDFContext &ctx, const SurfaceInteraction3f &si,
              const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        bool has_specular = ctx.is_enabled(BSDFFlags::GlossyReflection, 0),
             has_diffuse = ctx.is_enabled(BSDFFlags::DiffuseReflection, 1);

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);

        active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

        if (unlikely((!has_specular && !has_diffuse) || dr::none_or<false>(active)))
            return 0.f;

    

        // Determine which component should be sampled
        Float prob_specular = m_specular_sampling_weight,
              prob_diffuse  =  (1.f - m_specular_sampling_weight);

        if (unlikely(has_specular != has_diffuse))
            prob_specular = has_specular ? 1.f : 0.f;
        else
            prob_specular = prob_specular / (prob_specular + prob_diffuse);
        prob_diffuse = 1.f - prob_specular;

        Vector3f H = dr::normalize(wo + si.wi);

        MicrofacetDistribution distr(m_type, m_alpha, m_sample_visible);
        Float result = 0.f;
        if (m_sample_visible)
            result = distr.eval(H) * distr.smith_g1(si.wi, H) /
                     (4.f * cos_theta_i);
        else
            result = distr.pdf(si.wi, H) / (4.f * dr::dot(wo, H));
        result *= prob_specular;

        result += prob_diffuse * warp::square_to_cosine_hemisphere_pdf(wo);
        return result;
 // MI_SAMPLE_DIFFUSE
    }
 std::pair<Spectrum, Float> eval_pdf(const BSDFContext &ctx,
                                        const SurfaceInteraction3f &si,
                                        const Vector3f &wo,
                                        Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        bool has_specular = ctx.is_enabled(BSDFFlags::GlossyReflection, 0),
             has_diffuse  = ctx.is_enabled(BSDFFlags::DiffuseReflection, 1);

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);

        active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

        if (unlikely((!has_specular && !has_diffuse) || dr::none_or<false>(active)))
            return { 0.f, 0.f };

       

        // Determine which component should be sampled
        Float prob_specular =  m_specular_sampling_weight,
              prob_diffuse  =  (1.f - m_specular_sampling_weight);

        if (unlikely(has_specular != has_diffuse))
            prob_specular = has_specular ? 1.f : 0.f;
        else
            prob_specular = prob_specular / (prob_specular + prob_diffuse);
        prob_diffuse = 1.f - prob_specular;

        // Calculate the reflection half-vector
        Vector3f H = dr::normalize(wo + si.wi);

        MicrofacetDistribution distr(m_type, m_alpha, m_sample_visible);

        // Evaluate the microfacet normal distribution
        Float D = distr.eval(H);

        // Evaluate shadow/masking term for incoming direction
        Float smith_g1_wi = distr.smith_g1(si.wi, H);

        Float pdf = 0.f;
        if (m_sample_visible)
            pdf = D * smith_g1_wi / (4.f * cos_theta_i);
        else
            pdf = distr.pdf(si.wi, H) / (4.f * dr::dot(wo, H));
        pdf *= prob_specular;

        pdf += prob_diffuse * warp::square_to_cosine_hemisphere_pdf(wo);

        UnpolarizedSpectrum value(0.f);
        if (has_specular) {
            // Fresnel term
            

            // Smith's shadow-masking function
            Float G = distr.smith_g1(wo, H) * smith_g1_wi;

            // Calculate the specular reflection component
            value =  D * G / (4.f * cos_theta_i);

            if (m_specular_reflectance)
                value *= m_specular_reflectance->eval(si, active);
        }

        if (has_diffuse) {
             UnpolarizedSpectrum diff = m_diffuse_reflectance->eval(si, active);
            value += diff * (dr::InvPi<Float>  * cos_theta_o);
        }

        return { depolarizer<Spectrum>(value) & active, pdf };
                                        }
    std::string to_string() const override {
        std::ostringstream oss;
        oss << "merl_dec[" << std::endl
        << "  distribution = " << m_type << "," << std::endl
        << "  sample_visible = "           << m_sample_visible                    << "," << std::endl
        << "  alpha = "                    << m_alpha                             << "," << std::endl
        << "  diffuse_reflectance = "      << m_diffuse_reflectance               << "," << std::endl
         << "]";
        return oss.str();
    }

MI_DECLARE_CLASS()

private:
    std::string m_name;
     ref<Texture> m_diffuse_reflectance;
    ref<Texture> m_specular_reflectance;
     MicrofacetType m_type;
    Float m_alpha;
    Float m_specular_sampling_weight;
	merl<Float>* m_brdf;
    tabular<Float> * m_tabular;
     bool m_sample_visible;
  //  ggx<Float> *m_ggx;
 // typename microfacet<Float>::params m_params;
};

MI_IMPLEMENT_CLASS_VARIANT(MERLDEC, BSDF)
MI_EXPORT_PLUGIN(MERLDEC, "merl_dec")
NAMESPACE_END(mitsuba)
