#include <mitsuba/core/string.h>
#include <mitsuba/core/fwd.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/ior.h>
#include <mitsuba/render/microfacet.h>
#include <mitsuba/render/texture.h>

NAMESPACE_BEGIN(mitsuba)

/**!
.. _bsdf-roughconductor:

Rough conductor material (:monosp:`roughconductor`)
---------------------------------------------------

.. pluginparameters::

 * - material
   - |string|
   - Name of the material preset, see :num:`conductor-ior-list`. (Default: none)

 * - eta, k
   - |spectrum| or |texture|
   - Real and imaginary components of the material's index of refraction. (Default: based on the value of :monosp:`material`)
   - |exposed|, |differentiable|, |discontinuous|

 * - specular_reflectance
   - |spectrum| or |texture|
   - Optional factor that can be used to modulate the specular reflection component.
     Note that for physical realism, this parameter should never be touched. (Default: 1.0)
   - |exposed|, |differentiable|

 * - distribution
   - |string|
   - Specifies the type of microfacet normal distribution used to model the surface roughness.

     - :monosp:`beckmann`: Physically-based distribution derived from Gaussian random surfaces.
       This is the default.
     - :monosp:`ggx`: The GGX :cite:`Walter07Microfacet` distribution (also known as Trowbridge-Reitz
       :cite:`Trowbridge19975Average` distribution) was designed to better approximate the long
       tails observed in measurements of ground surfaces, which are not modeled by the Beckmann
       distribution.

 * - alpha, alpha_u, alpha_v
   - |texture| or |float|
   - Specifies the roughness of the unresolved surface micro-geometry along the tangent and
     bitangent directions. When the Beckmann distribution is used, this parameter is equal to the
     **root mean square** (RMS) slope of the microfacets. :monosp:`alpha` is a convenience
     parameter to initialize both :monosp:`alpha_u` and :monosp:`alpha_v` to the same value. (Default: 0.1)
   - |exposed|, |differentiable|, |discontinuous|

 * - sample_visible
   - |bool|
   - Enables a sampling technique proposed by Heitz and D'Eon :cite:`Heitz1014Importance`, which
     focuses computation on the visible parts of the microfacet normal distribution, considerably
     reducing variance in some cases. (Default: |true|, i.e. use visible normal sampling)

This plugin implements a realistic microfacet scattering model for rendering
rough conducting materials, such as metals.

.. subfigstart::
.. subfigure:: ../../resources/data/docs/images/render/bsdf_roughconductor_copper.jpg
   :caption: Rough copper (Beckmann, :math:`\alpha=0.1`)
.. subfigure:: ../../resources/data/docs/images/render/bsdf_roughconductor_anisotropic_aluminium.jpg
   :caption: Vertically brushed aluminium (Anisotropic Beckmann, :math:`\alpha_u=0.05,\ \alpha_v=0.3`)
.. subfigure:: ../../resources/data/docs/images/render/bsdf_roughconductor_textured_carbon.jpg
   :caption: Carbon fiber using two inverted checkerboard textures for ``alpha_u`` and ``alpha_v``
.. subfigend::
    :label: fig-bsdf-roughconductor


Microfacet theory describes rough surfaces as an arrangement of unresolved
and ideally specular facets, whose normal directions are given by a
specially chosen *microfacet distribution*. By accounting for shadowing
and masking effects between these facets, it is possible to reproduce the
important off-specular reflections peaks observed in real-world measurements
of such materials.

This plugin is essentially the *roughened* equivalent of the (smooth) plugin
:ref:`conductor <bsdf-conductor>`. For very low values of :math:`\alpha`, the two will
be identical, though scenes using this plugin will take longer to render
due to the additional computational burden of tracking surface roughness.

The implementation is based on the paper *Microfacet Models
for Refraction through Rough Surfaces* by Walter et al.
:cite:`Walter07Microfacet` and it supports two different types of microfacet
distributions.

To facilitate the tedious task of specifying spectrally-varying index of
refraction information, this plugin can access a set of measured materials
for which visible-spectrum information was publicly available
(see the corresponding table in the :ref:`conductor <bsdf-conductor>` reference).

When no parameters are given, the plugin activates the default settings,
which describe a 100% reflective mirror with a medium amount of roughness modeled
using a Beckmann distribution.

To get an intuition about the effect of the surface roughness parameter
:math:`\alpha`, consider the following approximate classification: a value of
:math:`\alpha=0.001-0.01` corresponds to a material with slight imperfections
on an otherwise smooth surface finish, :math:`\alpha=0.1` is relatively rough,
and :math:`\alpha=0.3-0.7` is **extremely** rough (e.g. an etched or ground
finish). Values significantly above that are probably not too realistic.


The following XML snippet describes a material definition for brushed aluminium:

.. tabs::
    .. code-tab:: xml
        :name: lst-roughconductor-aluminium

        <bsdf type="roughconductor">
            <string name="material" value="Al"/>
            <string name="distribution" value="ggx"/>
            <float name="alpha_u" value="0.05"/>
            <float name="alpha_v" value="0.3"/>
        </bsdf>

    .. code-tab:: python

        'type': 'roughconductor',
        'material': 'Al',
        'distribution': 'ggx',
        'alpha_u': 0.05,
        'alpha_v': 0.3

Technical details
*****************

All microfacet distributions allow the specification of two distinct
roughness values along the tangent and bitangent directions. This can be
used to provide a material with a *brushed* appearance. The alignment
of the anisotropy will follow the UV parameterization of the underlying
mesh. This means that such an anisotropic material cannot be applied to
triangle meshes that are missing texture coordinates.

Since Mitsuba 0.5.1, this plugin uses a new importance sampling technique
contributed by Eric Heitz and Eugene D'Eon, which restricts the sampling
domain to the set of visible (unmasked) microfacet normals. The previous
approach of sampling all normals is still available and can be enabled
by setting :monosp:`sample_visible` to :monosp:`false`. However this will lead
to significantly slower convergence.

When using this plugin, you should ideally compile Mitsuba with support for
spectral rendering to get the most accurate results. While it also works
in RGB mode, the computations will be more approximate in nature.
Also note that this material is one-sided---that is, observed from the
back side, it will be completely black. If this is undesirable,
consider using the :ref:`twosided <bsdf-twosided>` BRDF adapter.

In *polarized* rendering modes, the material automatically switches to a polarized
implementation of the underlying Fresnel equations.

 */

template <typename Float, typename Spectrum>
class Spec_Merl final : public BSDF<Float, Spectrum> {
public:
    MI_IMPORT_BASE(BSDF, m_flags, m_components)
    MI_IMPORT_TYPES(Texture, MicrofacetDistribution)

    Spec_Merl(const Properties &props) : Base(props) {
        std::string material = props.string("material", "none");
        m_type = MicrofacetType::GGX;

        m_sample_visible = props.get<bool>("sample_visible", true);
        
        m_alpha_u = props.texture<Texture>("alpha", 0.1f);
        specular = props.texture<Texture>("specular", 0.1f);
        s_mean=specular->mean();
        
 

        if (props.has_property("albedo")) 
            albedo = props.texture<Texture>("albedo", 0.7f);
        
       // d_mean = albedo->mean();
       d_mean=1.0f-s_mean;
       // s_mean=1.0f-d_mean;
       // m_specular_sampling_weight = s_mean / (d_mean + s_mean);

        m_components.push_back(BSDFFlags::GlossyReflection | BSDFFlags::FrontSide);
        m_components.push_back(BSDFFlags::DiffuseReflection | BSDFFlags::FrontSide);
        m_flags =  m_components[0] | m_components[1];
        dr::set_attr(this, "flags", m_flags);
        parameters_changed();
        m_sqrt_one_minus_rho_sqr= dr::safe_sqrt(1.0-m_rho*m_rho);
    }

    void traverse(TraversalCallback *callback) override {
        if (albedo)
            callback->put_object("albedo", albedo.get(), +ParamFlags::Differentiable);
        callback->put_object("alpha",                m_alpha_u.get(),              ParamFlags::Differentiable | ParamFlags::Discontinuous);
        callback->put_object("specular",                specular.get(),              ParamFlags::Differentiable | ParamFlags::Discontinuous);
    }
 void parameters_changed(const std::vector<std::string> &keys = {}) override {
        // Compute inverse of eta squared
    Float    phi_a=0.0;
    Float alpha=m_alpha_u->mean();
   Float cos_phi_a=dr::cos(phi_a);
  Float  sin_phi_a=dr::sin(phi_a);
  Float  cos_phi_a_sq=2.0*cos_phi_a*cos_phi_a-1.0;
  Float  a1_sqr=alpha*alpha;
  Float  a2_sqr=alpha*alpha;
  Float  tmp1=a1_sqr+a2_sqr;
   Float tmp2=a1_sqr-a2_sqr;
    m_ax=dr::sqrt(0.5*(tmp1+tmp2*cos_phi_a_sq));
   m_ay=dr::sqrt(0.5*(tmp1-tmp2*cos_phi_a_sq));
   m_rho=(a2_sqr-a1_sqr)*sin_phi_a*cos_phi_a/(m_ax*m_ay);

 }
    std::pair<BSDFSample3f, Spectrum> sample(const BSDFContext &ctx,
                                             const SurfaceInteraction3f &si,
                                             Float  sample1 ,
                                             const Point2f &sample2,
                                             Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);
        bool has_specular = ctx.is_enabled(BSDFFlags::GlossyReflection, 0),
             has_diffuse  = ctx.is_enabled(BSDFFlags::DiffuseReflection, 1);
        BSDFSample3f bs = dr::zeros<BSDFSample3f>();
        Float cos_theta_i = Frame3f::cos_theta(si.wi);
        active &= cos_theta_i > 0.f;

        /* Construct a microfacet distribution matching the
           roughness values at the current surface position. */
        if (unlikely((!has_specular && !has_diffuse) || dr::none_or<false>(active)))
            return { bs, 0.f };

        // Sample M, the microfacet normal
        

        // Perfect specular reflection based on the microfacet normal
       
        bs.eta = 1.f;
        Spectrum result(0.f);

        // Ensure that this is a valid sample
    
        Mask sample_specular = active && (sample1 <s_mean),
             sample_diffuse = active && !sample_specular; 

    Normal3f h;
     
        if (dr::any_or<true>(sample_specular)) {  
        /*    
        MicrofacetDistribution distr(m_type,
                                    m_alpha_u->eval_1(si, active),
                                    m_sample_visible);   
       
        std::tie(h, bs.pdf) = distr.sample(si.wi, sample2); 
         dr::masked(bs.wo, sample_specular) = reflect(si.wi, h);
*/

        Float u1=sample2.x()*0.99998+0.00001;
       Float u2=sample2.y()*0.99998+0.00001;
      Float a=si.wi.x()*m_ax+si.wi.y()*m_ay*m_rho;
       Float b=si.wi.y()*m_ay*m_sqrt_one_minus_rho_sqr;
       Float c=si.wi.z()-si.wi.x()*m_tx_n-si.wi.y()*m_ty_n;
      Vector3f o_std(a,b,c);
        o_std=dr::normalize(o_std);
        Float tx_m;
        Float ty_m;
     std::tie(tx_m,ty_m)=sample_vp22_std(u1,u2,o_std);
     
      Float  tx_h=m_ax*tx_m+m_tx_n;
      Float  choleski=m_rho*tx_m+m_sqrt_one_minus_rho_sqr*ty_m;
       Float ty_h=m_ay*choleski+m_ty_n;
     Normal3f h=Vector3f(-tx_h,-ty_h,1.0);
        h=dr::normalize(h);
        
        dr::masked(bs.wo, sample_specular) = dr::select(o_std.z()>0.0,reflect(si.wi, h),Vector3f(0.0f,0.0f,1.0f));
        
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
        
        return { bs, (depolarizer<Spectrum>(result)/ bs.pdf)  & active };
    }
std::pair<Float,Float> sample_vp22_std(Float u1,Float u2,const Vector3f o_std) const
{
 Float cos_theta_k=o_std.z();
 
 Float sin_theta_k=dr::select(cos_theta_k<1.0,dr::safe_sqrt(1.0-cos_theta_k*cos_theta_k),Float(0.0));

 Float tx=qf2_radial(u1,cos_theta_k,sin_theta_k);
 Float       ty=q3_radial(u2,tx);

  Float      xslope1=tx;
   Float     yslope1=ty;

   Float     nrm=1.0/dr::safe_sqrt(o_std.x()*o_std.x()+o_std.y()*o_std.y());
  Float      cos_phi_k=o_std.x()*nrm;
   Float     sin_phi_k=o_std.y()*nrm;
   Float     xslope2=cos_phi_k*tx-sin_phi_k*ty;
   Float     yslope2=sin_phi_k*tx+cos_phi_k*ty;
 Float  xslope=dr::select(sin_theta_k==0.0,xslope1,xslope2);
 Float   yslope=dr::select(sin_theta_k==0.0,yslope1,yslope2);
    return {xslope,yslope};
  
}
Float qf2_radial(Float u,Float cos_theta_k,Float sin_theta_k)const
{
Float sin_theta=u*(1.0+cos_theta_k)-1.0;
 Float   cos_theta=dr::safe_sqrt(1.0-sin_theta*sin_theta);
        
  Float      tan_theta=sin_theta/cos_theta;
   Float     tan_theta_k=sin_theta_k/cos_theta_k;
   Float     res1=-(tan_theta+tan_theta_k)/(1.0-tan_theta*tan_theta_k);
        
  Float      cot_theta_k2=cos_theta_k/sin_theta_k;
   Float     res2=(1.0+tan_theta*cot_theta_k2)/(tan_theta-cot_theta_k2);
   Float     value1=dr::select(sin_theta_k<0.707107,res1,res2);

   Float     cot_theta=cos_theta/sin_theta;
   Float     tan_theta_k2=sin_theta_k/cos_theta_k;
  Float      res3=(1.0+tan_theta_k2*cot_theta)/(tan_theta_k2-cot_theta);

   Float     cot_theta_k=cos_theta_k/sin_theta_k;
   Float     res4=(cot_theta+cot_theta_k)/(1.0-cot_theta*cot_theta_k);
  Float  value2=dr::select(sin_theta_k<0.707107,res3,res4);

    Float    value=dr::select(cos_theta>0.707107,value1,value2);
    return value;

}
Float q3_radial(Float u,Float qf2) const
{
  Float alpha=dr::safe_sqrt(1.0+qf2*qf2);
     Float   S=0;
     
    Float    u1=2.0*(0.5-u);
    Float    S1=-1.0;

    Float    u2=2.0*(u-0.5);
   Float     S2=1.0;
        u=dr::select(u<0.5,u1,u2);
       S=dr::select(u<0.5,S1,S2);
    Float    p=u * (u * (u * (-0.365728915865723) 
	          + 0.790235037209296) - 0.424965825137544) + 0.000152998850436920;
    Float    q=u * (u * (u * (u * 0.169507819808272 - 0.397203533833404) 
	          - 0.232500544458471) + 1) - 0.539825872510702;
    Float    value=S*alpha*(p/q);


return value;

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
       return result & active;
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
        /*   
        MicrofacetDistribution distr(m_type,
                                     m_alpha_u->eval_1(si, active),
                                     m_sample_visible);

        if (likely(m_sample_visible))
            result = distr.eval(h) * distr.smith_g1(si.wi, h) /
                     (4.f * cos_theta_i);
        else
            result = distr.pdf(si.wi, h) / (4.f * dr::dot(wo, h));
         */

       Float  g1_o=g1(h,wo);
        
        Float g1_i=g1(h,si.wi);
      Float  tmp=g1_o*g1_i;

        active &=tmp>0.0;
      Float  G=tmp/(g1_i+g1_o-tmp);
        active &=G>0.0;
        active &=kh>0.0;
        active &=h.z()>0.0;

      Float cos_therta_h_sqr=h.z()*h.z();
      Float  cos_theta_h_sqr_sqr=cos_therta_h_sqr*cos_therta_h_sqr;
      Float  xslope=-h.x()/h.z();
       Float yslope=-h.y()/h.z();
      Float  nrm_p22=m_ax*m_ay*m_sqrt_one_minus_rho_sqr;
       Float x=xslope-m_tx_n;
       Float y=yslope-m_ty_n;
       Float x_=x/m_ax;
       Float tmp1=m_ax*y-m_rho*m_ay*x;
      Float  tmp2=m_ax*m_ay*m_sqrt_one_minus_rho_sqr;
      Float  y_=tmp1/tmp2;
      Float  r_sqr=x_*x_+y_*y_;
        tmp=1.0+r_sqr;
      Float  p22_std=(1.0/(dr::Pi<Float>*tmp*tmp))/nrm_p22/cos_theta_h_sqr_sqr;
      Float  vndf=p22_std*kh/sigma(wo);
        result=vndf/(4.0*dr::dot(si.wi,h));



        result*=s_mean;
        result+=d_mean* warp::square_to_cosine_hemisphere_pdf(wo);
        return dr::select(active, result, 0.f);
    }
    Float g1(Vector3f wi,Vector3f wh)const
    {
        Float g1_loacal=dr::dot(wh,m_n);
        Float ss=sigma(wh);
         Float value=wh.z()/ss;
        return dr::select(g1_loacal>0.0,value,0.0f);
    }
    Float sigma(Vector3f wh) const
    { Float  a=wh.x()*m_ax+wh.y()*m_ay*m_rho;
      Float  b=wh.y()*m_ay*m_sqrt_one_minus_rho_sqr;
     Float   c=wh.z()-wh.x()*m_tx_n-wh.y()*m_ty_n;
     Float   nrm=dr::safe_sqrt(a*a+b*b+c*c);
      Vector3f ttt=Vector3f(a,b,c)/nrm;

      Float  result=nrm*((1.0+ttt.z())/2.0);

        return result;
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

        return {depolarizer<Spectrum>(value)& active, pdf };
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "RoughConductor[" << std::endl
            << "  distribution = " << m_type << "," << std::endl
            << "  sample_visible = " << m_sample_visible << "," << std::endl
            << "  alpha_u = " << string::indent(m_alpha_u) << "," << std::endl;
           oss << "  albedo = " << string::indent(albedo) << "," << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS()
private:
    /// Specifies the type of microfacet distribution
    MicrofacetType m_type;
    /// Anisotropic roughness values
    ref<Texture> m_alpha_u;
    ref<Texture> specular;
    /// Importance sample the distribution of visible normals?
    bool m_sample_visible;
    /// Relative refractive index (real component)
 //Float m_specular_sampling_weight;
 Float s_mean;
    Float d_mean;
    Float m_ax;
   Float m_ay;
    Float m_rho;
    ref<Texture> albedo;
   Vector3f m_n=Vector3f(0.0f,0.0f,0.0f);
    Float m_sqrt_one_minus_rho_sqr;
    Float m_tx_n=0.0f;
    Float m_ty_n=0.0f;
    
};

MI_IMPLEMENT_CLASS_VARIANT(Spec_Merl, BSDF)
MI_EXPORT_PLUGIN(Spec_Merl, "Spec_Merl")
NAMESPACE_END(mitsuba)
