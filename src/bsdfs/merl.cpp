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
        if (!fs::exists(file_path))
            Throw("File \"%s\" not found!", file_path);
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
     //  m_alpha = tabular<Float>::fit_ggx_parameter(*temp);
       m_alpha=12.5;
     //  Throw("Alpha\"%f\" calculated!", m_alpha);

//////////////////////////////////////////////////////////////////////////////////////
     m_components.push_back(BSDFFlags::GlossyReflection | BSDFFlags::FrontSide);
   //  m_components.push_back(BSDFFlags::DiffuseReflection | BSDFFlags::FrontSide);
        m_flags =  m_components[0] | m_components[1];
        dr::set_attr(this, "flags", m_flags);
        parameters_changed();
    }

    void traverse(TraversalCallback *callback) override {
     //   callback->put_object("diffuse_reflectance", m_diffuse_reflectance.get(), +ParamFlags::Differentiable);
        callback->put_parameter("alpha",            m_alpha,                     ParamFlags::Differentiable | ParamFlags::Discontinuous);
        
    }
  void parameters_changed(const std::vector<std::string> &keys = {}) override {
        
       
    }
    

    std::pair<BSDFSample3f, Spectrum> sample(const BSDFContext &ctx,
                                             const SurfaceInteraction3f &si,
                                             Float sample1,
                                             const Point2f &sample2,
                                             Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);

        bool has_specular = ctx.is_enabled(BSDFFlags::GlossyReflection, 0);
         //    has_diffuse  = ctx.is_enabled(BSDFFlags::DiffuseReflection, 1);

        Float cos_theta_i = Frame3f::cos_theta(si.wi);
        active &= cos_theta_i > 0.f;

        BSDFSample3f bs = dr::zeros<BSDFSample3f>();
        Spectrum result(0.f);
        if (!ctx.is_enabled(BSDFFlags::GlossyReflection) ||dr::none_or<false>(active))
            return {bs, result};

        bs.eta = 1.f;

        MicrofacetDistribution distr(m_type, m_alpha, m_sample_visible);
        Normal3f m = std::get<0>(distr.sample(si.wi, sample2));

        dr::masked(bs.wo, has_specular) = reflect(si.wi, m);
           
        
        bs.pdf = pdf(ctx, si, bs.wo, active);
        active &= bs.pdf > 0.f;
        result = eval(ctx, si, bs.wo, active);

        return {bs, (UnpolarizedSpectrum(result) / bs.pdf) & active };
    }

     Spectrum eval(const BSDFContext &ctx, const SurfaceInteraction3f &si,
                  const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        Vector3f wi = si.wi;

        active &= Frame3f::cos_theta(wi) > 0.f &&
                  Frame3f::cos_theta(wo) > 0.f;

        if (!ctx.is_enabled(BSDFFlags::GlossyReflection) || dr::none_or<false>(active))
            return Spectrum(0.f);
 
        Vector3f m = dr::normalize(wo + wi);
  
		Vector3f fr_p = m_brdf->evalp(wi, wo);
     
        return UnpolarizedSpectrum (fr_p) & active;
    }

    Float pdf(const BSDFContext &ctx, const SurfaceInteraction3f &si,
              const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        bool has_specular = ctx.is_enabled(BSDFFlags::GlossyReflection, 0);
       

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);

        active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

        if (unlikely((!has_specular) || dr::none_or<false>(active)))
            return 0.f;

        Vector3f H = dr::normalize(wo + si.wi);

        MicrofacetDistribution distr(m_type, m_alpha, m_sample_visible);
        Float result = 0.f;
        if (m_sample_visible)
            result = distr.eval(H) * distr.smith_g1(si.wi, H) /
                     (4.f * cos_theta_i);
        else
            result = distr.pdf(si.wi, H) / (4.f * dr::dot(wo, H));

        return result;
 // MI_SAMPLE_DIFFUSE
    }
 
    std::string to_string() const override {
        std::ostringstream oss;
        oss << "merl_dec[" << std::endl
        << "  distribution = " << m_type << "," << std::endl
        << "  sample_visible = "           << m_sample_visible                    << "," << std::endl
        << "  alpha = "                    << m_alpha                             << "," << std::endl
         << "]";
        return oss.str();
    }

MI_DECLARE_CLASS()

private:
    std::string m_name;

     MicrofacetType m_type;
    Float m_alpha;

	merl<Float>* m_brdf;
    tabular<Float> * m_tabular;
     bool m_sample_visible;
};

MI_IMPLEMENT_CLASS_VARIANT(MERLDEC, BSDF)
MI_EXPORT_PLUGIN(MERLDEC, "merl_dec")
NAMESPACE_END(mitsuba)
