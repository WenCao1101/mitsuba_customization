#include <mitsuba/core/properties.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/microfacet.h>
#include <mitsuba/render/texture.h>

NAMESPACE_BEGIN(mitsuba)

/**!

.. _bsdf-diffuse:

Smooth diffuse material (:monosp:`diffuse`)
-------------------------------------------

.. pluginparameters::

 * - reflectance
   - |spectrum| or |texture|
   - Specifies the diffuse albedo of the material (Default: 0.5)
   - |exposed|, |differentiable|

The smooth diffuse material (also referred to as *Lambertian*)
represents an ideally diffuse material with a user-specified amount of
reflectance. Any received illumination is scattered so that the surface
looks the same independently of the direction of observation.

.. subfigstart::
.. subfigure:: ../../resources/data/docs/images/render/bsdf_diffuse_plain.jpg
   :caption: Homogeneous reflectance
.. subfigure:: ../../resources/data/docs/images/render/bsdf_diffuse_textured.jpg
   :caption: Textured reflectance
.. subfigend::
   :label: fig-diffuse

Apart from a homogeneous reflectance value, the plugin can also accept
a nested or referenced texture map to be used as the source of reflectance
information, which is then mapped onto the shape based on its UV
parameterization. When no parameters are specified, the model uses the default
of 50% reflectance.

Note that this material is one-sided---that is, observed from the
back side, it will be completely black. If this is undesirable,
consider using the :ref:`twosided <bsdf-twosided>` BRDF adapter plugin.
The following XML snippet describes a diffuse material,
whose reflectance is specified as an sRGB color:

.. tabs::
    .. code-tab:: xml
        :name: diffuse-srgb

        <bsdf type="diffuse">
            <rgb name="reflectance" value="0.2, 0.25, 0.7"/>
        </bsdf>

    .. code-tab:: python

        'type': 'diffuse',
        'reflectance': {
            'type': 'rgb',
            'value': [0.2, 0.25, 0.7]
        }

Alternatively, the reflectance can be textured:

.. tabs::
    .. code-tab:: xml
        :name: diffuse-texture

        <bsdf type="diffuse">
            <texture type="bitmap" name="reflectance">
                <string name="filename" value="wood.jpg"/>
            </texture>
        </bsdf>

    .. code-tab:: python

        'type': 'diffuse',
        'reflectance': {
            'type': 'bitmap',
            'filename': 'wood.jpg'
        }
*/
template <typename Float, typename Spectrum>
class Spec_Test final : public BSDF<Float, Spectrum> {
public:
    MI_IMPORT_BASE(BSDF, m_flags, m_components)
    MI_IMPORT_TYPES(Texture, MicrofacetDistribution)

    Spec_Test(const Properties &props) : Base(props) {
        std::string material = props.string("material", "none");
        m_type = MicrofacetType::GGX;

        m_sample_visible = props.get<bool>("sample_visible", true);
        
        m_alpha_u = props.texture<Texture>("alpha", 0.1f);
 
        if (props.has_property("albedo")) 
            albedo = props.texture<Texture>("albedo", 0.7f);
        
        d_mean = albedo->mean();
        s_mean=1.0f-d_mean;
        m_specular_sampling_weight = s_mean / (d_mean + s_mean);

        m_components.push_back(BSDFFlags::GlossyReflection | BSDFFlags::FrontSide);
        m_components.push_back(BSDFFlags::DiffuseReflection | BSDFFlags::FrontSide);
        m_flags =  m_components[0] | m_components[1];
        dr::set_attr(this, "flags", m_flags);
        m_components.push_back(m_flags);
    }

    void traverse(TraversalCallback *callback) override {
        callback->put_object("albedo",albedo.get(), +ParamFlags::Differentiable);
        callback->put_object("alpha",      m_alpha_u.get(),  ParamFlags::Differentiable | ParamFlags::Discontinuous);
    }

    std::pair<BSDFSample3f, Spectrum> sample(const BSDFContext &ctx,
                                             const SurfaceInteraction3f &si,
                                             Float  sample1 ,
                                             const Point2f &sample2,
                                             Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);
        bool has_specular = ctx.is_enabled(BSDFFlags::GlossyReflection, 0),
             has_diffuse  = ctx.is_enabled(BSDFFlags::DiffuseReflection, 1);
        Float cos_theta_i = Frame3f::cos_theta(si.wi);
        BSDFSample3f bs = dr::zeros<BSDFSample3f>();

        active &= cos_theta_i > 0.f;
       if (unlikely((!has_specular && !has_diffuse) || dr::none_or<false>(active)))
            return { bs, 0.f };
        bs.eta = 1.f;
        Spectrum result(0.f);
        Mask sample_specular = active && (sample1 <m_specular_sampling_weight),
             sample_diffuse = active && !sample_specular; 

        Normal3f h;
        if (dr::any_or<true>(sample_specular)) {  
           
        MicrofacetDistribution distr(m_type,
                                    m_alpha_u->eval_1(si, active),
                                    m_sample_visible);   
       
        std::tie(h, bs.pdf) = distr.sample(si.wi, sample2); 
        dr::masked(bs.wo, sample_specular) = reflect(si.wi, h);

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

        return { bs, (depolarizer<Spectrum>(result)/ bs.pdf) & active };
    }

    Spectrum eval(const BSDFContext &ctx, const SurfaceInteraction3f &si,
                  const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        bool has_specular = ctx.is_enabled(BSDFFlags::GlossyReflection, 0),
             has_diffuse  = ctx.is_enabled(BSDFFlags::DiffuseReflection, 1);
        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);

        active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

         if (unlikely((!has_specular && !has_diffuse) || dr::none_or<false>(active)))
            return 0.f;

        // Calculate the half-direction vector
        Vector3f H = dr::normalize(wo + si.wi);

        /* Construct a microfacet distribution matching the
           roughness values at the current surface position. */
        MicrofacetDistribution distr(m_type,
                                     m_alpha_u->eval_1(si, active),
                                     m_sample_visible);

        // Evaluate the microfacet normal distribution
        Float D = distr.eval(H);

        active &= dr::neq(D, 0.f);

        // Evaluate Smith's shadow-masking function
        Float G = distr.G(si.wi, wo, H);

        // Evaluate the full microfacet model (except Fresnel)
        UnpolarizedSpectrum result = D * G / (4.f * Frame3f::cos_theta(si.wi));

        Spectrum F;
  
        F = UnpolarizedSpectrum(dr::dot(si.wi, H));
        

        /* If requested, include the specular reflectance component */
       // if (m_specular_reflectance)
       //     result *= m_specular_reflectance->eval(si, active);
       result=F*result* s_mean;
       result+=albedo->eval(si,active)* dr::InvPi<Float> * cos_theta_o;

        return depolarizer<Spectrum>(result) & active;
    }

    Float pdf(const BSDFContext &ctx, const SurfaceInteraction3f &si,
              const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        bool has_specular = ctx.is_enabled(BSDFFlags::GlossyReflection, 0),
             has_diffuse = ctx.is_enabled(BSDFFlags::DiffuseReflection, 1);
        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);

        // Calculate the half-direction vector
        Vector3f h = dr::normalize(wo + si.wi);
        Float kh=dr::dot(wo,h);
 
        active &= cos_theta_i > 0.f && cos_theta_o > 0.f &&
                  dr::dot(si.wi, h) > 0.f && dr::dot(wo, h) > 0.f;

        if (unlikely((!has_specular && !has_diffuse) || dr::none_or<false>(active)))
            return 0.f;
       Float result=0.0; 
        /* Construct a microfacet distribution matching the
           roughness values at the current surface position. */
          
        MicrofacetDistribution distr(m_type,
                                     m_alpha_u->eval_1(si, active),
                                     m_sample_visible);

        if (likely(m_sample_visible))
            result = distr.eval(h) * distr.smith_g1(si.wi, h) /
                     (4.f * cos_theta_i);
        else
            result = distr.pdf(si.wi, h) / (4.f * dr::dot(wo, h));

        result*=s_mean;
        result+=d_mean* warp::square_to_cosine_hemisphere_pdf(wo);

        return dr::select(active, result, 0.f);
    }

    std::pair<Spectrum, Float> eval_pdf(const BSDFContext &ctx,
                                        const SurfaceInteraction3f &si,
                                        const Vector3f &wo,
                                        Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);

        // Calculate the half-direction vector
        Vector3f H = dr::normalize(wo + si.wi);

        /* Filter cases where the micro/macro-surface don't agree on the side.
           This logic is evaluated in smith_g1() called as part of the eval()
           and sample() methods and needs to be replicated in the probability
           density computation as well. */
        active &= cos_theta_i > 0.f && cos_theta_o > 0.f &&
                  dr::dot(si.wi, H) > 0.f && dr::dot(wo, H) > 0.f;

        if (unlikely(!ctx.is_enabled(BSDFFlags::GlossyReflection) || dr::none_or<false>(active)))
            return { 0.f, 0.f };

        /* Construct a microfacet distribution matching the
           roughness values at the current surface position. */
        MicrofacetDistribution distr(m_type,
                                     m_alpha_u->eval_1(si, active),
                                     m_sample_visible);

        // Evaluate the microfacet normal distribution
        Float D = distr.eval(H);

        active &= dr::neq(D, 0.f);

        // Evaluate Smith's shadow-masking function
        Float smith_g1_wi = distr.smith_g1(si.wi, H);
        Float G = smith_g1_wi * distr.smith_g1(wo, H);

        // Evaluate the full microfacet model (except Fresnel)
        UnpolarizedSpectrum value = D * G / (4.f * Frame3f::cos_theta(si.wi));

        Spectrum F;
     
        F = UnpolarizedSpectrum(dr::dot(si.wi, H));
        value=F*value* s_mean;
        value+=albedo->eval(si,active)* dr::InvPi<Float> * Frame3f::cos_theta(wo);  

        // If requested, include the specular reflectance component
     //   if (m_specular_reflectance)
     //       value *= m_specular_reflectance->eval(si, active);

        Float pdf;
        if (likely(m_sample_visible))
            pdf = D * smith_g1_wi / (4.f * cos_theta_i);
        else
            pdf = distr.pdf(si.wi, H) / (4.f * dr::dot(wo, H));
        pdf*=s_mean;
        pdf+=d_mean* warp::square_to_cosine_hemisphere_pdf(wo);

        return { depolarizer<Spectrum>(value) & active, dr::select(active, pdf, 0.f) };
    }
    std::string to_string() const override {
        std::ostringstream oss;
        oss << "Spec_Test[" << std::endl
            << "albedo = " << string::indent(albedo) << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS()
private:
     MicrofacetType m_type;
    /// Anisotropic roughness values
    ref<Texture> m_alpha_u;
    /// Importance sample the distribution of visible normals?
    bool m_sample_visible;
    /// Relative refractive index (real component)
    Float m_specular_sampling_weight;
    Float s_mean;
    Float d_mean;

    ref<Texture> albedo;
};

MI_IMPLEMENT_CLASS_VARIANT(Spec_Test, BSDF)
MI_EXPORT_PLUGIN(Spec_Test, "Spec_Test")
NAMESPACE_END(mitsuba)
