#pragma once


#ifndef DJB_INCLUDE_DJ_BRDF_H
#define DJB_INCLUDE_DJ_BRDF_H


#include <mitsuba/core/properties.h>
#include <mitsuba/core/warp.h>



#ifndef DJB_EPSILON
#define DJB_EPSILON 1e-4f
#endif

#ifndef DJB_ASSERT
#	include <assert.h>
#	define DJB_ASSERT(x) assert(x)
#endif
/*
#ifndef DJB_LOG
#	include <stdio.h>
#	define DJB_LOG(format, ...) fprintf(stdout, format, ##__VA_ARGS__)
#endif
*/
NAMESPACE_BEGIN(mitsuba)
/* Exception API */
struct exc : public std::exception {
	exc(const char *fmt, ...);
	virtual ~exc() throw() {}
	const char *what() const throw() {return m_str.c_str();}
	std::string m_str;
};
template<typename Float> static Float sat(const Float& x) {return dr::minimum(dr::maximum(x,0.0f),1.0f);}

template <typename Float>

Float intensity( Vector<Float, 3> inty) { return 0.2126 * inty.x() + 0.7152 * inty.y() + 0.0722 * inty.z();}
template <typename Float>
static void
uniform_to_concentric(Float u1, Float u2, Float *r1, Float *r2)
{   using Mask   = dr::mask_t<Float>;
	/* Concentric map code with less branching (by Dave Cline), see
	   http://psgraphics.blogspot.ch/2011/01/improved-code-for-concentric-map.html */
//	Float r1 = 2.0 * u1 - 1.0;
//	Float r2 = 2.0 * u2 - 1.0;
	/*
	Float phi, r;

	if (r1 == 0 && r2 == 0) {
		r = phi = 0;
	} else if (dr::abs(r1 * r1) > dr::abs(r2 * r2)) {
		r = r1;
		phi = (M_PI / 4.0) * (r2 / r1);
	} else {
		r = r2;
		phi = (M_PI / 2.0) - (r1 / r2) * (M_PI / 4.0);
	}

	(*x) = r * dr::cos(phi);
	(*y) = r * dr::sin(phi);*/
	Float x = dr::fmsub(2.f, u1, 1.f),
          y = dr::fmsub(2.f, u2, 1.f);
    Mask is_zero         = dr::eq(x, 0.f) &&
                           dr::eq(y, 0.f),
         quadrant_1_or_3 = dr::abs(x) < dr::abs(y);

    Float r  = dr::select(quadrant_1_or_3, y, x),
          rp = dr::select(quadrant_1_or_3, x, y);

    Float phi = 0.25f * dr::Pi<Float> * rp / r;
    dr::masked(phi, quadrant_1_or_3) = 0.5f * dr::Pi<Float> - phi;
    dr::masked(phi, is_zero) = 0.f;

    auto [s, c] = dr::sincos(phi);
	(*r1) = r * c;
	(*r2) = r * s;


}

/* BRDF interface */
template <typename Float>
class brdf {
public:
     MI_IMPORT_CORE_TYPES()
   
	// evaluate f_r
	virtual  Vector3f eval(const Vector3f& i, const Vector3f& o,
	                  const void *user_param = NULL) const = 0;
	virtual Vector3f eval_hd(const Vector3f& h, const Vector3f& d,
	                     const void *user_param = NULL) const
						 {   
	Vector3f i, o;
	hd_to_io(h, d, &i, &o);

	return eval(i, o, user_param);
}
	// evaluate f_r * cos
	virtual Vector3f evalp(const Vector3f& i, const Vector3f& o,
	                   const void *user_param = NULL) const{
	return eval(i, o, user_param) * i.z();
}
	virtual Vector3f evalp_hd(const Vector3f& h, const Vector3f& d,
	                      const void *user_param = NULL) const{
	Vector3f i, o;
	hd_to_io(h, d, &i, &o);

	return eval(i, o, user_param) * i.z();
}
	// evaluate f_r * cos / pdf
	virtual Vector3f evalp_is(Float u1, Float u2,
	                      const Vector3f& o,
	                      Vector3f *i, Float *pdf,
	                      const void *user_param = NULL) const{
	const Vector3f i_ = sample(u1, u2, o, user_param);
	Float pdf_ = this->pdf(i_, o);
	if (i) (*i) = i_;
	if (pdf) (*pdf) = pdf_;

	return (evalp(i_, o, user_param) / pdf_);
}
	// importance sample f_r * cos using two uniform numbers
	virtual Vector3f sample(Float u1, Float u2,
	                    const Vector3f& o,
	                    const void *user_param = NULL) const{
	Float x, y, z;

	uniform_to_concentric(u1, u2, &x, &y);
	z = dr::sqrt(1.0 - x * x - y * y);
	DJB_ASSERT(z >= 0.0);

	return Vector3f(x, y, z);
}
	// evaluate the PDF of a sample
	virtual Float pdf(const Vector3f& i, const Vector3f& o,
	                    const void *user_param = NULL) const{
	return /* cos(theta_i) */i.z() / dr::Pi<Float>;
}
	// utilities
	static void io_to_hd(const Vector3f& i, const Vector3f& o, Vector3f *h, Vector3f *d);
	static void hd_to_io(const Vector3f& h, const Vector3f& d, Vector3f *i, Vector3f *o);
	// ctor / dtor
	brdf() {}
	virtual ~brdf() {}
#if 1 // noncopyable
private:
	brdf(const brdf& fr);
	brdf& operator=(const brdf& fr);
#endif
};

/* Lambertian BRDF */
template <typename Float>
class lambert : public brdf<Float> {
public:
     MI_IMPORT_CORE_TYPES()
  
	/* Lambertian Parameters */
	class params {
	public:
		params(const Vector3f& reflectance = Vector3f(1)):m_reflectance(reflectance)
{}
		Vector3f m_reflectance;
	};
	/* Implementation */
	Vector3f eval(const Vector3f& i, const Vector3f& o,
	          const void *user_param = NULL) const{
	const lambert::params params = 
		user_param ? *reinterpret_cast<const lambert::params *>(user_param)
		           : lambert::params();

	return (params.m_reflectance / dr::Pi<Float>);
};

};

//NAMESPACE_BEGIN(mitsuba)
/* MERL BRDF */

//NAMESPACE_END(mitsuba)
#define DJB__UTIA_STEP_T  15.0
#define DJB__UTIA_STEP_P   7.5
#define DJB__UTIA_NTI      6
#define DJB__UTIA_NPI     48
#define DJB__UTIA_NTV      6
#define DJB__UTIA_NPV     48
#define DJB__UTIA_PLANES   3
/* UTIA BRDF */
template <typename Float>
class utia : public brdf<Float> {
	std::vector<double> m_samples;
	double m_norm;
public:
    MI_IMPORT_CORE_TYPES()
	
	utia(const char *filename);
	Vector3f eval(const Vector3f& i, const Vector3f& o,
	          const void *user_param = NULL) const{
	Float r2d = 180.0 / dr::Pi<Float>;
	Float theta_i = r2d * dr::acos(i.z());
	Float theta_o = r2d * dr::acos(o.z());
	Float phi_i = r2d * dr::atan2(i.y(), i.x());
	Float phi_o = r2d * dr::atan2(o.y(), o.x());

	// make sure we're above horizon
	if (theta_i >= 90.0 || theta_o >= 90.0)
		return Vector3f(0);

	// make sure phi is in [0, 360)
	while (phi_i < 0.0) {phi_i+= 360.0;}
	while (phi_o < 0.0) {phi_o+= 360.0;}
	while (phi_i >= 360) {phi_i-= 360.0;}
	while (phi_o >= 360) {phi_o-= 360.0;}

	uint32_t iti[2], itv[2], ipi[2], ipv[2];
	//iti[0] = dr::floor(theta_i / DJB__UTIA_STEP_T);
	iti[0]=1;
	iti[1] = iti[0] + 1;
	if(iti[0] > DJB__UTIA_NTI - 2) {
		iti[0] = DJB__UTIA_NTI - 2;
		iti[1] = DJB__UTIA_NTI - 1;
	}
	itv[0] = dr::floor2int<uint32_t>(theta_o / DJB__UTIA_STEP_T);
	//itv[0]=1;
	itv[1] = itv[0] + 1;
	if(itv[0] > DJB__UTIA_NTV - 2) {
		itv[0] = DJB__UTIA_NTV - 2;
		itv[1] = DJB__UTIA_NTV - 1;
	}

	ipi[0] = dr::floor2int<uint32_t>(phi_i / DJB__UTIA_STEP_P);
	//ipi[0] = 1;
	ipi[1] = ipi[0] + 1;
	ipv[0] = dr::floor2int<uint32_t>(phi_o / DJB__UTIA_STEP_P);
	//ipv[0] = 1;
	ipv[1] = ipv[0] + 1;

	Float sum;
	Float wti[2], wtv[2], wpi[2], wpv[2];
	wti[1] = theta_i - (Float)(DJB__UTIA_STEP_T * iti[0]);
	wti[0] = (Float)(DJB__UTIA_STEP_T * iti[1]) - theta_i;
	sum = wti[0] + wti[1];
	wti[0]/= sum;
	wti[1]/= sum;
	wtv[1] = theta_o - (Float)(DJB__UTIA_STEP_T * itv[0]);
	wtv[0] = (Float)(DJB__UTIA_STEP_T * itv[1]) - theta_o;
	sum = wtv[0] + wtv[1];
	wtv[0]/= sum;
	wtv[1]/= sum;

	wpi[1] = phi_i - (Float)(DJB__UTIA_STEP_P * ipi[0]);
	wpi[0] = (Float)(DJB__UTIA_STEP_P * ipi[1]) - phi_i;
	sum = wpi[0] + wpi[1];
	wpi[0]/= sum;
	wpi[1]/= sum;
	wpv[1] = phi_o - (Float)(DJB__UTIA_STEP_P * ipv[0]);
	wpv[0] = (Float)(DJB__UTIA_STEP_P * ipv[1]) - phi_o;
	sum = wpv[0] + wpv[1];
	wpv[0]/= sum;
	wpv[1]/= sum;

	if(ipi[1] == DJB__UTIA_NPI)
		ipi[1] = 0;
	if(ipv[1] == DJB__UTIA_NPV)
		ipv[1] = 0;

	uint32_t nc = DJB__UTIA_NPV * DJB__UTIA_NTV;
	uint32_t nr = DJB__UTIA_NPI * DJB__UTIA_NTI;
	Float RGB[DJB__UTIA_PLANES];
	for(uint32_t isp = 0; isp < DJB__UTIA_PLANES; ++isp) {
		uint32_t i, j, k, l;

		RGB[isp] = 0.0;
		for(i = 0; i < 2; ++i)
		for(j = 0; j < 2; ++j)
		for(k = 0; k < 2; ++k)
		for(l = 0; l < 2; ++l) {
			Float w = wti[i] * wtv[j] * wpi[k] * wpv[l];
			uint32_t idx = isp * nr * nc + nc * (DJB__UTIA_NPI * iti[i] + ipi[k]) 
			        + DJB__UTIA_NPV * itv[j] + ipv[l];

			RGB[isp]+= w * (Float)m_samples[idx];
		}
		if (RGB[isp] > 0.0375)
			RGB[isp] = pow((Float)(RGB[isp] + 0.055)/1.055, (Float)2.4); 
		else
			RGB[isp]/= (Float)12.92;
		RGB[isp]*= (Float)100.0; 
	}
	return Vector3f(
		max((Float)0, RGB[0]),
		max((Float)0, RGB[1]),
		max((Float)0, RGB[2])
	);
}

	const std::vector<double>& get_samples() const {return m_samples;}
private:
	void normalize();
};

/* Fresnel API */


namespace fresnel {
	/* Utilities */
	template <typename Float>
	
	class utils {
	public:
	MI_IMPORT_CORE_TYPES()
    
	void ior_to_f0(Float ior, Float *f0);
	void ior_to_f0(const Vector3f& ior, Vector3f *f0);
	void f0_to_ior(Float f0, Float *ior);
	void f0_to_ior(const Vector3f& f0, Vector3f *ior);
	};
	/* Interface */
	template <typename Float>
	class impl {
	public:
      MI_IMPORT_CORE_TYPES()
	
		virtual ~impl() {}
		virtual Vector3f eval(Float cos_theta_d) const = 0;
		virtual impl *copy() const = 0;
	};

	/* Ideal Specular Reflection */
	template <typename Float>
	class ideal : public impl<Float> {
	public:
      MI_IMPORT_CORE_TYPES()

		Vector3f eval(Float cos_theta_d) const {return Vector3f(1);}
		impl <Float> *copy() const {return new ideal();}
	};

	/* Fresnel for Unpolarized Light */
	template <typename Float>
	class unpolarized : public impl<Float> {
		MI_IMPORT_CORE_TYPES()
		Vector3f ior; // index of refraction
	public:
        
        
		unpolarized(const Vector3f& ior): ior(ior) {}
		Vector3f eval(Float cos_theta_d) const{
		DJB_ASSERT(cos_theta_d >= 0.f && cos_theta_d <= 1.0 && "Invalid Angle");
	Vector3f  F;

	//for (int i = 0; i < 3; ++i)
		F.x() = unpolarized__eval(cos_theta_d, Vector3f::to_raw(ior).x());
		F.y() = unpolarized__eval(cos_theta_d, Vector3f::to_raw(ior).y());
		F.z() = unpolarized__eval(cos_theta_d, Vector3f::to_raw(ior).z());

	return {F.x(), F.y(), F.z()};}
		impl<Float> *copy() const {return new unpolarized(*this);}
	};

	/* Schlick's Fresnel */
	template <typename Float>
	class schlick : public impl<Float> {
		 MI_IMPORT_CORE_TYPES()
		Vector3f f0; // backscattering fresnel
	public:
       
		schlick(const Vector3f& f0): f0(f0) {}
		Vector3f eval(Float cos_theta_d) const{
        DJB_ASSERT(cos_theta_d >= 0.f && cos_theta_d <= 1.0 && "Invalid Angle");
	    Float c1 = 1.0 - cos_theta_d;
	    Float c2 = c1 * c1;
	    Float c5 = c2 * c2 * c1;

	    return f0 + c5 * (Vector3f(1) - f0);

		}


		impl<Float> *copy() const {return new schlick(*this);}
	};

	/* SGD's Fresnel */
	template <typename Float>
	class sgd : public impl<Float> {
		MI_IMPORT_CORE_TYPES()
		Vector3f f0, f1;
	public:
       
       
		sgd(const Vector3f& f0, const Vector3f& f1) : f0(f0), f1(f1) {}
		Vector3f eval(Float cos_theta_d) const{
	DJB_ASSERT(cos_theta_d >= 0.f && cos_theta_d <= 1.0 && "Invalid Angle");
	Float c = cos_theta_d;

	return f0 - c * f1 + pow(1.0 - c, 5.0) * (Vector3f(1) - f0);
}
		impl<Float> *copy() const {return new sgd(*this);}
	};

	/* Arbitrary Fresnel Function */
	template <typename Float>
	class spline : public impl<Float> {
		MI_IMPORT_CORE_TYPES()
		std::vector<Vector3f >m_points;
	public:
        
		explicit spline(const std::vector<Vector3f>& points): m_points(points) {}
		const std::vector<Vector3f>& get_points() const {return m_points;}
		Vector3f eval(Float cos_theta_d) const;
		impl<Float> *copy() const {return new spline(*this);}
	};
} // namespace fresnel

/* Microfacet API */
template <typename Float>
class microfacet : public brdf<Float> {
public:
    MI_IMPORT_CORE_TYPES()
	/* microfacet parameters */
	
	class params {
		friend class microfacet;
	public:
		// Factories
		static params standard();
		static params isotropic(Float a);
		static params elliptic(Float a1, Float a2, Float phi_a = 0.0);
		static params pdfparams(Float ax, Float ay, Float rho = 0.0,
		                        Float tx_n = 0.0, Float ty_n = 0.0);
		// mutators
		void set_ellipse(Float a1, Float a2, Float phi_a = 0.0);
		void set_pdfparams(Float ax, Float ay, Float rho = 0.0,
		                   Float tx_n = 0.0, Float ty_n = 0.0);
		void set_location(Float tx_n, Float ty_n);
		void set_location(const Vector3f& n);
		// accessors
		void get_ellipse(Float *a1, Float *a2, Float *phi_a = NULL) const;
		void get_pdfparams(Float *ax, Float *ay, Float *rho = NULL,
		                   Float *tx_n = NULL, Float *ty_n = NULL) const;
		void get_location(Float *tx_n, Float *ty_n) const;
		void get_location(Vector3f *n) const;
		// Ctors (prefer factories for construction)
		params(Float a1 = 1.0, Float a2 = 1.0, Float phi_a = 0.0);
		params(Float ax, Float ay, Float rho, Float tx_n, Float ty_n);
	private:
		Vector3f m_n; // mean normal
		Float m_a1, m_a2, m_phi_a; // ellipse parameters
		Float m_ax, m_ay; // scale parameters
		Float m_rho, m_sqrt_one_minus_rho_sqr; // correlation
		Float m_tx_n, m_ty_n; // location parameters
	};
	// Dtor
	virtual ~microfacet() {delete m_fresnel;}
	// BRDF interface
	Vector3f eval(const Vector3f& i, const Vector3f& o,
	          const void *user_param = NULL) const{
	return (evalp(i, o, user_param) / i.z());
}

	Vector3f evalp(const Vector3f& i, const Vector3f& o,
	           const void *user_param = NULL) const{
	const microfacet::params params = 
		user_param ? *reinterpret_cast<const microfacet::params *>(user_param)
		           : microfacet::params::standard();

	Vector3f h = normalize(i + o);
	Float G = gaf(h, i, o, params);
    Mask active=true;
	active &= G > 0.0;
	if (dr::none_or<true>(active)) {
		Float cos_theta_d = sat(dr::dot(o, h));
		Vector3f F = fresnel(cos_theta_d);
		Float D = ndf(h, params);

		return (F * ((D * G) / (4.0 * o.z())));
	}
	return Vector3f(0);
}
	Vector3f sample(Float u1, Float u2, const Vector3f& o,
	            const void *user_param = NULL) const{
	// get params
	const microfacet<Float>::params params = 
		user_param ? *reinterpret_cast<const microfacet<Float>::params *>(user_param)
		           : microfacet<Float>::params::standard();

#if 1 // clamp samples
	u1 = sat(u1) * (Float)0.99998 + (Float)0.00001;
	u2 = sat(u2) * (Float)0.99998 + (Float)0.00001;
#endif

	// warp the receiver direction
	Float a = o.x() * params.m_ax + o.y() * params.m_ay * params.m_rho;
	Float b = o.y() * params.m_ay * params.m_sqrt_one_minus_rho_sqr;
	Float c = o.z() - o.x() * params.m_tx_n - o.y() * params.m_ty_n;
	Vector3f o_std = dr::normalize(Vector3f(a, b, c));
    Mask active=true;
	active &= o_std.z() > 0.0;
	if (dr::none_or<true>(active)) {
		// create a standard variate
		Float tx_m, ty_m;
		sample_vp22_std_smith(u1, u2, o_std, &tx_m, &ty_m);

		// warp the variate with the microfacet parameters
		Float tx_h = params.m_ax * tx_m + params.m_tx_n;
		Float choleski = params.m_rho * tx_m
			             + params.m_sqrt_one_minus_rho_sqr * ty_m;
		Float ty_h = params.m_ay * choleski + params.m_ty_n;

		// get the associated normal
		Vector3f h = dr::normalize(Vector3f(-tx_h, -ty_h, 1));

		// return the reflected direction around the micronormal
		return (2.0 * dr::dot(o, h) * h - o);
	}
	return Vector3f(0, 0, 1);
}
	Float pdf(const Vector3f& i, const Vector3f& o,
	            const void *user_param = NULL) const;
	Vector3f evalp_is(Float u1, Float u2, const Vector3f& o,
	              Vector3f *i, Float *pdf, const void *user_param = NULL) const{
	const microfacet::params params = 
		user_param ? *reinterpret_cast<const microfacet::params *>(user_param)
		           : microfacet::params::standard();
	Vector3f i_ = sample(u1, u2, o, user_param);
	Vector3f h = dr::normalize(i_ + o);
	Float G = gaf(h, i_, o, params);
	if (pdf) (*pdf) = 0.f;

    Mask active=true;
	active &= G > 0.0;
	if (dr::none_or<true>(active)) {
		Float cos_theta_d = sat(dot(o, h));

		if (i) (*i) = i_;
		if (!supports_smith_vndf_sampling()) {
			Float pdf_ = (h.z() * ndf(h, params) / (4.0 * cos_theta_d));
			if (pdf) (*pdf) = pdf_;
			return evalp(i_, o, user_param) / pdf_;
		} else {
			Vector3f F = fresnel(cos_theta_d);
			Float G1 = g1(h, o, params);
			if (pdf) (*pdf) = (vndf(h, o, params) / (4.0 * cos_theta_d));
			return (F * (G / G1));
		}
	}
	return Vector3f(0);
}
	// eval queries
	Vector3f fresnel(Float cos_theta_d) const {return m_fresnel->eval(cos_theta_d);}
	Float ndf(const Vector3f& h,
	            const params& params = params::standard()) const;
	Float gaf(const Vector3f& h, const Vector3f& i, const Vector3f& o,
	            const params& params = params::standard()) const;
	Float g1(const Vector3f& h, const Vector3f& k,
	           const params& params = params::standard()) const;
	Float sigma(const Vector3f& k,
	              const params& params = params::standard()) const;
	Float p22(Float x, Float y,
	            const params& params = params::standard()) const;
	Float vp22(Float x, Float y, const Vector3f& k,
	             const params& params = params::standard()) const;
	Float vndf(const Vector3f& h, const Vector3f& k,
	             const params& params = params::standard()) const;
	// sampling queries
	virtual bool supports_smith_vndf_sampling() const = 0;
	virtual Float qf2(Float u, const Vector3f& k) const;
	virtual Float qf3(Float u, const Vector3f& k, Float qf2) const;
	// mutators
	void set_shadow(bool shadow) {m_shadow = shadow;}
	void set_fresnel(const fresnel::impl<Float>& f){
	delete m_fresnel;
	m_fresnel = f.copy();
}
	// accessors
	uint32_t get_shadow() const {return m_shadow;}
	const fresnel::impl<Float>& get_fresnel() const {return *m_fresnel;}
protected:
	// ctor
	microfacet(const fresnel::impl<Float>& f = fresnel::ideal<Float>(),
	           bool shadow = true):m_fresnel(f.copy()), m_shadow(shadow){}
	// microfacet interfaces
	virtual Float sigma_std(const Vector3f& k) const = 0; // standard masking
	virtual Float p22_std(Float x, Float y) const = 0; // standard slope pdf
	// sampling functions
	virtual void sample_vp22_std_smith(Float u1, Float u2, const Vector3f& k,
	                                   Float *xslope, Float *yslope) const;
	virtual void sample_vp22_std_nmap(Float u1, Float u2, const Vector3f& k,
	                                  Float *xslope, Float *yslope) const = 0;
	// members
	const fresnel::impl<Float> *m_fresnel;
	bool m_shadow;
};

/* Radial microfacets */
template <typename Float>
class radial : public microfacet<Float> {
public:
	MI_IMPORT_CORE_TYPES()
	radial(const fresnel::impl<Float>& fresnel = fresnel::ideal<Float>(),
	       bool shadow = true): microfacet<Float>(fresnel, shadow)
	{}
	// queries
	virtual Float p22_radial(Float r_sqr) const = 0;
	virtual Float sigma_std_radial(Float cos_theta_k) const = 0;
	virtual Float cdf_radial(Float r) const = 0;
	virtual Float qf_radial(Float u) const = 0;
	virtual Float qf2_radial(Float u, 
	                           Float cos_theta_k,
	                           Float sin_theta_k) const;
	virtual Float qf3_radial(Float u, Float qf2) const;
	bool supports_smith_vndf_sampling() const {return true;}
private:
	// eval
	Float p22_std(Float x, Float y) const;
	Float sigma_std(const Vector3f& k) const;
	// sample
	void sample_vp22_std_nmap(Float u1, Float u2, const Vector3f& o,
	                          Float *xslope, Float *yslope) const;
	void sample_vp22_std_smith(Float u1, Float u2, const Vector3f& o,
	                           Float *xslope, Float *yslope) const;
};

/* Beckmann Microfacet NDF */
template <typename Float>		
class beckmann : public radial<Float> {
public:
   
    MI_IMPORT_CORE_TYPES()
	// Linear Representation
	class lrep {
		friend class beckmann;
		public:
		// Ctor
		lrep(Float E1 = 0, Float E2 = 0,
		     Float E3 = 1, Float E4 = 1,
		     Float E5 = 0);
		// linear operators
		lrep operator+(const lrep& r) const // combine
		{
	     return beckmann::lrep(m_E1 + r.m_E1,
	                      m_E2 + r.m_E2,
	                      m_E3 + r.m_E3 + (Float)2.0 * m_E1 * r.m_E1,
	                      m_E4 + r.m_E4 + (Float)2.0 * m_E2 * r.m_E2,
	                      m_E5 + r.m_E5 + m_E1 * r.m_E2 + m_E2 * r.m_E1);
          }
		lrep operator*(Float sc) const{
	DJB_ASSERT(sc >= (Float)0.0 && "Invalid scale");
	Float sc_sqr = sc * sc;

	return beckmann::lrep(m_E1 * sc, m_E2 * sc,
	                      m_E3 * sc_sqr, m_E4 * sc_sqr,
	                      m_E5 * sc_sqr);
} // scale
		lrep& operator+=(const lrep& rep){
	m_E1+= rep.m_E1;
	m_E2+= rep.m_E2;
	m_E3+= rep.m_E3 + (Float)2.0 * m_E1 * rep.m_E1;
	m_E4+= rep.m_E4 + (Float)2.0 * m_E2 * rep.m_E2;
	m_E5+= rep.m_E5 + m_E1 * rep.m_E2 + m_E2 * rep.m_E1;

	return (*this);
} // compound combine
		lrep& operator*=(Float sc){
	DJB_ASSERT(sc >= (Float)0.0 && "Invalid scale");
	Float sc_sqr = sc * sc;
	m_E1*= sc;
	m_E2*= sc;
	m_E3*= sc_sqr;
	m_E4*= sc_sqr;
	m_E5*= sc_sqr;

	return (*this);
} // compound scale
		// utilities
		void scale(Float sc);
		void scale(Float x, Float y);
		void shear(Float x, Float y);
		void reparameterize(Float dudx, Float dvdx,
		                    Float dudy, Float dvdy);
	private:
		// members
		Float m_E1, m_E2; // first order slope moments
		Float m_E3, m_E4; // second order slope moments
		Float m_E5; // first order joint slope moment
	};
	// utilities
	static void params_to_lrep(const typename microfacet<Float>::params& params, lrep *lrep){
	DJB_ASSERT(lrep && "Null output ptr");
	Float ax, ay, rho, tx_n, ty_n;
	params.get_pdfparams(&ax, &ay, &rho, &tx_n, &ty_n);
	(*lrep) = beckmann::lrep(tx_n, ty_n,
	                         (Float)0.5 * ax * ax + tx_n * tx_n,
	                         (Float)0.5 * ay * ay + ty_n * ty_n,
	                         (Float)0.5 * rho * ax * ay + tx_n * ty_n);
}
	static void lrep_to_params(const lrep& lrep, typename microfacet<Float>::params *params){
	DJB_ASSERT(params && "Null output ptr");
	Float tx_n = lrep.m_E1;
	Float ty_n = lrep.m_E2;
	Float tmp1 = max((Float)0.0, lrep.m_E3 - lrep.m_E1 * lrep.m_E1);
	Float tmp2 = max((Float)0.0, lrep.m_E4 - lrep.m_E2 * lrep.m_E2);
	/* clamp the parameters to valid values */
	Float ax = (Float)max(1e-5, sqrt(2.0 * tmp1));
	Float ay = (Float)max(1e-5, sqrt(2.0 * tmp2));
	Float rho = (Float)2.0 * (lrep.m_E5 - lrep.m_E1 * lrep.m_E2) / (ax * ay);
	rho = min((Float)0.99, max((Float)-0.99, rho));

	(*params) = microfacet<Float>::params::pdfparams(ax, ay, rho, tx_n, ty_n);
}
	// ctor
	beckmann(const fresnel::impl<Float>& fresnel = fresnel::ideal<Float>(),
	         bool shadow = true): radial<Float>(fresnel, shadow)
	{}
	// queries
	Float p22_radial(Float r_sqr) const;
	Float sigma_std_radial(Float cos_theta_k) const;
	Float cdf_radial(Float r) const;
	Float qf_radial(Float u) const;
	Float qf1(Float u) const;
	Float qf2_radial(Float u,
	                   Float cos_theta_k, Float sin_theta_k) const;
	Float qf3_radial(Float u, Float qf2) const;
	bool supports_smith_vndf_sampling() const {return true;}
};

/* GGX Microfacet NDF */
template <typename Float>	
class ggx : public radial<Float> {
public:
    MI_IMPORT_CORE_TYPES()
	
	ggx(const fresnel::impl<Float>& fresnel = fresnel::ideal<Float>(),
	    bool shadow = true): radial<Float>(fresnel, shadow)
	{}
	// queries
	Float p22_radial(Float r_sqr) const;
	Float sigma_std_radial(Float cos_theta_k) const;
	Float cdf_radial(Float r) const;
	Float qf_radial(Float u) const;
	Float qf1(Float u) const;
	Float qf2_radial(Float u,
	                   Float cos_theta_k, Float sin_theta_k) const;
	Float qf3_radial(Float u, Float qf2) const;
	bool supports_smith_vndf_sampling() const {return true;}
private:
	Float qf3_rational_approx(Float u) const;
};

/* Tabulated Microfacet NDF */
template <typename Float>
class tabular : public radial<Float> {
	std::vector<Float> m_p22;
	std::vector<Float> m_sigma;
	// tables for Nmap VNDF sampling
	std::vector<Float> m_cdf;
	std::vector<Float> m_qf;
public:
    MI_IMPORT_CORE_TYPES()
	tabular(const brdf<Float>& brdf, uint32_t resolution, bool shadow = true):radial<Float>(){
	DJB_ASSERT(resolution > 2 && "Invalid Resolution");
//	void set_shadow(bool m_shadow){shadow = m_shadow;}

	// allocate memory
	m_p22.reserve(resolution);
	m_sigma.reserve(resolution);
	m_cdf.reserve(resolution);
	m_qf.reserve(resolution);

	// eval
	
	compute_p22_smith(brdf, resolution);
	normalize_p22();
	compute_sigma();
	//compute_fresnel(brdf, resolution);

	// sample
	//compute_cdf();
	//compute_qf();
	}

	static typename microfacet<Float>::params fit_beckmann_parameters(const tabular<Float>& tabular)
	{const uint32_t ntheta = 128;
	Float dtheta = dr::Pi<Float> / (Float)ntheta;
	Float nint = 0.0;
	Float alpha;

	for (uint32_t i = 0; i < ntheta; ++i) {
		Float u = (Float)i / (Float)ntheta; // in [0,1)
		Float theta_h = u * u * dr::Pi<Float> * 0.5;
		Float cos_theta_h = cos(theta_h);
		Float r_h = tan(theta_h);
		Float r_h_sqr = r_h * r_h;
		Float p22_r = tabular.p22_radial(r_h_sqr);

		nint+= (u * r_h_sqr * r_h * p22_r) / (cos_theta_h * cos_theta_h);
	}
	nint*= dtheta * /* int_0^2pi cos^2 phi dphi */dr::Pi<Float>;
	alpha = sqrt(2.0 * nint);

#ifndef NVERBOSE
    Throw("Invalid file structure: %f", alpha);
	//DJB_LOG("djb_verbose: Beckmann_alpha = %.9f\n", alpha);
#endif

	return microfacet<Float>::params::isotropic(alpha);}

	static typename microfacet<Float>::params fit_ggx_parameters(tabular<Float>& tabular)
	{
	const int32_t ntheta = 128;
	Float dtheta = dr::Pi<Float> / (Float)ntheta;
	Float nint = 0.0;
	Float alpha;

	for (uint32_t i = 0; i < ntheta; ++i) {
		Float u = i / ntheta; // in [0,1)
		Float theta_h = u * u * dr::Pi<Float> * 0.5;
		Float cos_theta_h = dr::cos(theta_h);
		Float r_h = dr::tan(theta_h);
		Float r_h_sqr = r_h * r_h;
		Float p22_r = tabular.p22_radial(r_h_sqr);

		nint+= (u * r_h_sqr * p22_r) / (cos_theta_h * cos_theta_h);
	}
	nint*= dtheta * /* int_0^2pi fabs(cos(phi)) dphi */4.0;
	alpha = nint;
#ifndef NVERBOSE
    Throw("Invalid file structure: %f", alpha);
	//DJB_LOG("djb_verbose: GGX_alpha = %.9f\n", alpha);
#endif

	return microfacet<Float>::params::isotropic(alpha);
}
	const std::vector<Float>& get_p22v() const {return m_p22;}
	const std::vector<Float>& get_sigmav() const {return m_sigma;}
	const std::vector<Float>& get_cdfv() const {return m_cdf;}
	const std::vector<Float>& get_qfv() const {return m_qf;}
	// queries
	Float p22_radial(Float r_sqr) const;
	Float sigma_std_radial(Float cos_theta_k) const;
	Float cdf_radial(Float r) const;
	Float qf_radial(Float u) const;
	bool supports_smith_vndf_sampling() const {return false;}
private:
	// eval precomputations
	void compute_p22_smith(const brdf<Float>& brdf, uint32_t res);
	void compute_fresnel(const brdf<Float>& brdf, uint32_t res);
	void normalize_p22();
	void compute_sigma();
	// sample precomputations
	void compute_cdf();
	void compute_qf();
	void compute_pdf1();
	void normalize_pdf1();
};

/* Tabulated Anisotropic Microfacet NDF */
template <typename Float>
class tabular_anisotropic : public microfacet<Float> {
	std::vector<Float> m_p22;
	std::vector<Float> m_sigma;
	std::vector<Float> m_pdf1;
	std::vector<Float> m_pdf2;
	std::vector<Float> m_cdf1;
	std::vector<Float> m_cdf2;
	std::vector<Float> m_qf1;
	std::vector<Float> m_qf2;
	uint32_t m_elevation_res;
	uint32_t m_azimuthal_res;
	bool supports_smith_vndf_sampling() const {return false;}
public:
   MI_IMPORT_CORE_TYPES()
	tabular_anisotropic(const brdf<Float>& brdf,
	                    uint32_t elevation_res,
	                    uint32_t azimuthal_res,
	                    bool shadow = true);
	static typename microfacet<Float>::params fit_beckmann_parameters(const tabular_anisotropic<Float>& tabular)
	{
	const uint32_t ntheta = 128;
	const uint32_t nphi   = 512;
	Float dtheta = dr::sqrt(dr::Pi<Float> * 0.5) / (Float)ntheta;
	Float dphi = 2.0 * dr::Pi<Float> / (Float)nphi;
	Float nint[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
	Float mux, muy, alphax, alphay, rho;

	for (uint32_t j = 0; j < nphi; ++j) {
		Float tmp = (Float)j / (Float)nphi; // in [0,1)
		Float phi = tmp * 2.0 * dr::Pi<Float>; // in [0,2pi)
		Float cos_phi = cos(phi);
		Float sin_phi = sin(phi);
		Float cos_phi_sqr = cos_phi * cos_phi;
		Float sin_phi_sqr = sin_phi * sin_phi;

		for (uint32_t i = 0; i < ntheta; ++i) {
			Float tmp1 = (Float)i / (Float)ntheta; // in [0,1)
			Float theta = tmp1 * dr::sqrt(dr::Pi<Float> * 0.5); // in [0, sqrt(pi/2))
			Float theta_sqr = theta * theta; // in [0, pi/2)
			Float p22 = tabular.p22_std_theta_phi(theta_sqr, phi);
			Float tan_theta = tan(theta_sqr);
			Float cos_theta = cos(theta_sqr);
			Float tan_theta_sqr = tan_theta * tan_theta;
			Float cos_theta_sqr = cos_theta * cos_theta;
			Float tmp2 = theta * p22 * tan_theta / cos_theta_sqr;
			Float e1 = -tan_theta * cos_phi; // mux
			Float e2 = -tan_theta * sin_phi; // muy
			Float e3 = tan_theta_sqr * cos_phi_sqr; // var1
			Float e4 = tan_theta_sqr * sin_phi_sqr; // var2
			Float e5 = tan_theta_sqr * cos_phi * sin_phi; // cov

			nint[0]+= tmp2 * e1;
			nint[1]+= tmp2 * e2;
			nint[2]+= tmp2 * e3;
			nint[3]+= tmp2 * e4;
			nint[4]+= tmp2 * e5;
		}
	}
	for (uint32_t i = 0; i < 5; i++)
		nint[i]*= 2.0 * dtheta * dphi;
	mux = nint[0];
	muy = nint[1];
	alphax = sqrt((Float)2.0 * (nint[2] - mux * mux));
	alphay = sqrt((Float)2.0 * (nint[3] - muy * muy));
	rho = 2.0 * (nint[4] - mux * muy) / (alphax * alphay);

#ifndef NVERBOSE
    Throw("Invalid file structure: %.9f %.9f %.9f %.9f %.9f", alphax,
		alphay,
		rho,
		mux,
		muy);
	//DJB_LOG(
	//	"djb_verbose: Beckmann params = %.9f %.9f %.9f %.9f %.9f\n", 
	//	alphax,
	//	alphay,
	//	rho,
	//	mux,
	//	muy
	//);
#endif

	return microfacet<Float>::params::pdfparams(alphax, alphay, rho, mux, muy);
}
	static typename microfacet<Float>::params fit_ggx_parameters(const tabular_anisotropic<Float>& tabular)
	{
	const uint32_t ntheta = 128;
	const uint32_t nphi   = 512;
	Float dtheta = dr::sqrt(dr::Pi<Float> * 0.5) / (Float)ntheta;
	Float dphi = 2.0 * dr::Pi<Float> / (Float)nphi;
	Float nint[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
	Float mux, muy, alphax, alphay, rho;

	for (uint32_t j = 0; j < nphi; ++j) {
		Float tmp = (Float)j / (Float)nphi; // in [0,1)
		Float phi = tmp * 2.0 * dr::Pi<Float>; // in [0,2pi)
		Float cos_phi = cos(phi);
		Float sin_phi = sin(phi);

		for (uint32_t i = 0; i < ntheta; ++i) {
			Float tmp1 = (Float)i / (Float)ntheta; // in [0,1)
			Float theta = tmp1 * sqrt(dr::Pi<Float> * 0.5); // in [0, sqrt(pi/2))
			Float theta_sqr = theta * theta; // in [0, pi/2)
			Float p22 = tabular.p22_std_theta_phi(theta_sqr, phi);
			Float tan_theta = tan(theta_sqr);
			Float cos_theta = cos(theta_sqr);
			Float cos_theta_sqr = cos_theta * cos_theta;
			Float tmp2 = theta * p22 * tan_theta / cos_theta_sqr;
			Float e1 = -tan_theta * cos_phi; // mux
			Float e2 = -tan_theta * sin_phi; // muy
			Float e3 = fabs(e1); // var1
			Float e4 = fabs(e2); // var2
			Float e5 = 0; // TODO

			nint[0]+= tmp2 * e1;
			nint[1]+= tmp2 * e2;
			nint[2]+= tmp2 * e3;
			nint[3]+= tmp2 * e4;
			nint[4]+= tmp2 * e5;
		}
	}
	for (uint32_t i = 0; i < 5; i++)
		nint[i]*= 2.0 * dtheta * dphi;
	mux = nint[0];
	muy = nint[1];
	alphax = sqrt(nint[2] * nint[2] - mux * mux);
	alphay = sqrt(nint[3] * nint[3] - muy * muy);
	rho = 0.0; // TODO

#ifndef NVERBOSE
    Throw("Invalid file structure: %.9f %.9f %.9f %.9f %.9f", alphax,
		alphay,
		rho,
		mux,
		muy);
	//DJB_LOG(
//		"djb_verbose: GGX params = %.9f %.9f %.9f %.9f %.9f\n", 
//		alphax,
//		alphay,
//		rho,
//		mux,
//		muy
//	);
#endif

	return microfacet<Float>::params::pdfparams(alphax, alphay, rho, mux, muy);
}

	const std::vector<Float>& get_p22v(uint32_t *elev_cnt, uint32_t *azim_cnt) const;
	const std::vector<Float>& get_sigmav(uint32_t *elev_cnt, uint32_t *azim_cnt) const;
	// queries
	Float pdf1(Float phi) const;
	Float pdf2(Float theta, Float phi) const;
	Float cdf1(Float phi) const;
	Float cdf2(Float theta, Float phi) const;
	Float qf1(Float u1) const;
	Float qf2(Float u, Float phi) const;
private:
	// eval interface
	Float sigma_std(const Vector3f& k) const; // standard masking
	Float p22_std(Float x, Float y) const; // standard slope pdf
	Float p22_std_theta_phi(Float theta, Float phi) const; // standard slope pdf
	// sample interface
	void sample_vp22_std_nmap(Float u1, Float u2, const Vector3f& k,
	                          Float *xslope, Float *yslope) const;
	// eval precomputations
	void compute_p22_smith(const brdf<Float>& brdf);
	void compute_fresnel(const brdf<Float>& brdf, uint32_t res);
	void normalize_p22();
	void compute_sigma();
	// sample precomputations
	void compute_pdf1();
	void compute_pdf2();
	void normalize_pdf1();
	void normalize_pdf2();
	void compute_cdf1();
	void compute_cdf2();
	void compute_qf1();
	void compute_qf2();
};

/* Shifted Gamma Distribution BRDF */
template <typename Float>
class sgd : public brdf<Float> {
	struct data {
		const char *name;
		const char *otherName;
		double rhoD[3];
		double rhoS[3];
		double alpha[3];
		double p[3];
		double f0[3];
		double f1[3];
		double kap[3];
		double lambda[3];
		double c[3];
		double k[3];
		double theta0[3];
		double error[3];
	};
	static const data s_data[100]= {
	{ "alum-bronze", "alum-bronze", { 0.0478786, 0.0313514, 0.0200638 }, { 0.0364976, 0.664975, 0.268836 }, { 0.014832, 0.0300126, 0.0490339 }, { 0.459076, 0.450056, 0.529272 }, { 6.05524, 0.235756, 0.580647 }, { 5.05524, 0.182842, 0.476088 }, { 46.3841, 24.5961, 14.8261 }, { 2.60672, 2.97371, 2.7827 }, { 1.12717e-07, 1.06401e-07, 5.27952e-08 }, { 47.783, 36.2767, 31.6066 }, { 0.205635, 0.066289, -0.0661091 }, { 0.100735, 0.0878706, 0.0861907 } },
	{ "alumina-oxide", "alumina-oxide", { 0.316358, 0.292248, 0.25416 }, { 0.00863128, 0.00676832, 0.0103309 }, { 0.000159222, 0.000139421, 0.000117714 }, { 0.377727, 0.318496, 0.402598 }, { 0.0300766, 1.70375, 1.96622 }, { -0.713784, 0.70375, 1.16019 }, { 4381.96, 5413.74, 5710.42 }, { 3.31076, 4.93831, 2.84538 }, { 6.72897e-08, 1.15769e-07, 6.32199e-08 }, { 354.275, 367.448, 414.581 }, { 0.52701, 0.531166, 0.53301 }, { 0.213276, 0.147418, 0.27746 } },
	{ "aluminium", "aluminium", { 0.0305166, 0.0358788, 0.0363463 }, { 0.0999739, 0.131797, 0.0830361 }, { 0.0012241, 0.000926487, 0.000991844 }, { 0.537669, 0.474562, 0.435936 }, { 0.977854, 0.503108, 1.77905 }, { -0.0221457, -0.0995445, 0.77905 }, { 449.321, 658.044, 653.86 }, { 8.2832e-07, 9.94692e-08, 6.11887e-08 }, { 3.54592e-07, 16.0175, 15.88 }, { 23.8656, 10.6911, 9.69801 }, { -0.510356, 0.570179, 0.566156 }, { 0.303567, 0.232628, 0.441578 } },
	{ "aventurnine", "aventurnine", { 0.0548217, 0.0621179, 0.0537826 }, { 0.0348169, 0.0872381, 0.111961 }, { 0.000328039, 0.000856166, 0.00145342 }, { 0.387167, 0.504525, 0.652122 }, { 0.252033, 0.133897, 0.087172 }, { 0.130593, 0.0930416, 0.0567429 }, { 2104.51, 676.157, 303.59 }, { 3.12126, 2.50965e-07, 2.45778e-05 }, { 1.03849e-07, 8.53824e-07, 3.20722e-07 }, { 251.265, 24.2886, 29.0236 }, { 0.510125, -0.41764, -0.245097 }, { 0.0359759, 0.0297523, 0.0285881 } },
	{ "beige-fabric", "fabric-beige", { 0.20926, 0.160666, 0.145337 }, { 0.121663, 0.0501577, 0.00177279 }, { 0.39455, 0.15975, 0.110706 }, { 0.474725, 0.0144728, 1.70871e-12 }, { 0.0559459, 0.222268, 8.4764 }, { -0.318718, -0.023826, 7.4764 }, { 3.8249, 7.32453, 10.0904 }, { 2.26283, 2.97144, 3.55311 }, { 0.0375346, 0.073481, 0.0740222 }, { 7.52635, 9.05672, 10.6185 }, { 0.217453, 0.407084, 0.450203 }, { 0.00217528, 0.00195262, 0.00171008 } },
	{ "black-fabric", "fabric-black", { 0.0189017, 0.0112353, 0.0110067 }, { 2.20654e-16, 6.76197e-15, 1.57011e-13 }, { 0.132262, 0.128044, 0.127838 }, { 0.189024, 0.18842, 0.188426 }, { 1, 1, 1 }, { 0, 0, 0 }, { 8.1593, 8.38075, 8.39184 }, { 3.83017, 3.89536, 3.89874 }, { 0.00415117, 0.00368324, 0.00365826 }, { 12.9974, 13.2597, 13.2737 }, { 0.207997, 0.205597, 0.205424 }, { 0.000363154, 0.000272253, 0.000274773 } },
	{ "black-obsidian", "black-obsidian", { 0.00130399, 0.0011376, 0.00107233 }, { 0.133029, 0.125362, 0.126188 }, { 0.000153649, 0.000148939, 0.000179285 }, { 0.186234, 0.227495, 0.25745 }, { 2.42486e-12, 0.0174133, 0.091766 }, { -0.0800755, -0.048671, 0.0406445 }, { 5668.57, 5617.79, 4522.84 }, { 13.7614, 8.59526, 6.44667 }, { 1e+38, 1e+38, 1e+38 }, { 117.224, 120.912, 113.366 }, { 1.19829, 1.19885, 1.19248 }, { 0.0841374, 0.0943085, 0.123514 } },
	{ "black-oxidized-steel", "black-oxidized-steel", { 0.0149963, 0.0120489, 0.0102471 }, { 0.373438, 0.344382, 0.329202 }, { 0.187621, 0.195704, 0.200503 }, { 0.661367, 0.706913, 0.772267 }, { 0.0794166, 0.086518, 0.080815 }, { 0.0470402, 0.0517633, 0.0455037 }, { 5.1496, 4.91636, 4.69009 }, { 4.0681, 3.95489, 3.71052 }, { 1.07364e-07, 1.05341e-07, 1.16556e-07 }, { 20.2383, 20.1786, 20.2553 }, { -0.479617, -0.4885, -0.478388 }, { 0.00129576, 0.0011378, 0.000986163 } },
	{ "black-phenolic", "black-phenolic", { 0.00204717, 0.00196935, 0.00182908 }, { 0.177761, 0.293146, 0.230592 }, { 0.00670804, 0.00652009, 0.00656043 }, { 0.706648, 0.677776, 0.673986 }, { 0.16777, 0.12335, 0.166663 }, { 0.111447, 0.0927321, 0.125663 }, { 65.4189, 70.8936, 70.9951 }, { 1.06318, 1.15283, 1.16529 }, { 1.24286e-07, 3.00039e-08, 9.77334e-08 }, { 74.0711, 75.1165, 73.792 }, { 0.338204, 0.319306, 0.33434 }, { 0.0307129, 0.0531183, 0.0454238 } },
	{ "black-soft-plastic", "black-plastic-soft", { 0.00820133, 0.00777718, 0.00764537 }, { 0.110657, 0.0980322, 0.100579 }, { 0.0926904, 0.0935964, 0.0949975 }, { 0.14163, 0.148703, 0.143694 }, { 0.150251, 0.169418, 0.170457 }, { 0.100065, 0.113089, 0.114468 }, { 11.2419, 11.113, 10.993 }, { 4.3545, 4.3655, 4.31586 }, { 0.00464641, 0.00384785, 0.0046145 }, { 14.6751, 14.8089, 14.5436 }, { 0.275651, 0.262317, 0.271284 }, { 0.000527402, 0.000550255, 0.000527718 } },
	{ "blue-acrylic", "acrylic-blue", { 0.0134885, 0.0373766, 0.10539 }, { 0.0864901, 0.0228191, 0.204042 }, { 0.000174482, 0.000269795, 0.0015211 }, { 0.373948, 0.362425, 0.563636 }, { 0.0185562, 0.399982, 0.0525861 }, { -0.0209713, 0.241543, 0.0169474 }, { 4021.24, 2646.36, 346.898 }, { 3.38722, 3.62885, 1.83684e-06 }, { 9.64334e-08, 9.96105e-08, 3.61787e-07 }, { 338.073, 272.828, 23.5039 }, { 0.526039, 0.515404, -0.526935 }, { 0.0612235, 0.0789826, 0.0461093 } },
	{ "blue-fabric", "fabric-blue", { 0.0267828, 0.0281546, 0.066668 }, { 0.0825614, 0.0853369, 0.0495164 }, { 0.248706, 0.249248, 0.18736 }, { 9.23066e-13, 1.66486e-12, 2.27218e-12 }, { 0.201626, 0.213723, 0.56548 }, { 0.225891, 0.226267, 0.638493 }, { 5.15615, 5.14773, 6.43713 }, { 2.25846, 2.25536, 2.68382 }, { 0.128037, 0.128363, 0.0944915 }, { 6.95531, 6.94633, 8.17665 }, { 0.407528, 0.407534, 0.411378 }, { 0.000323043, 0.00032064, 0.000653631 } },
	{ "blue-metallic-paint2", "ch-ball-blue-metallic", { 0.010143, 0.0157349, 0.0262717 }, { 0.0795798, 0.0234493, 0.0492337 }, { 0.00149045, 0.00110477, 0.00141008 }, { 0.624615, 0.598721, 0.67116 }, { 9.36434e-14, 3.61858e-15, 1.15633e-14 }, { -0.210234, -1, -1 }, { 314.024, 441.812, 299.726 }, { 1.20935e-05, 7.51792e-06, 3.86474e-05 }, { 3.38901e-07, 2.94502e-07, 3.15718e-07 }, { 27.0491, 28.576, 30.6214 }, { -0.326593, -0.274443, -0.187842 }, { 0.0908879, 0.163236, 0.286541 } },
	{ "blue-metallic-paint", "metallic-blue", { 0.00390446, 0.00337319, 0.00848198 }, { 0.0706771, 0.0415082, 0.104423 }, { 0.155564, 0.139, 0.15088 }, { 1.01719, 1.02602, 1.16153 }, { 0.149347, 0.153181, 1.87241e-14 }, { -0.487331, -0.76557, -1 }, { 4.4222, 4.59265, 3.93929 }, { 2.54345, 2.33884, 2.24405 }, { 6.04906e-08, 5.81858e-08, 1.2419e-07 }, { 23.9533, 25.0641, 24.6856 }, { -0.34053, -0.294595, -0.258117 }, { 0.00202544, 0.00246268, 0.0059725 } },
	{ "blue-rubber", "blue-rubber", { 0.0371302, 0.0732915, 0.146637 }, { 0.384232, 0.412357, 0.612608 }, { 0.218197, 0.2668, 0.478375 }, { 0.815054, 1.00146, 1.24995 }, { 0.0631713, 0.0622636, 0.0399196 }, { 0.0478254, 0.0422186, 0.007015 }, { 4.41586, 3.76795, 3.46276 }, { 3.77807, 3.82679, 3.33186 }, { 1.2941e-07, 1.07194e-07, 0.00045665 }, { 19.8046, 19.3115, 11.4364 }, { -0.499472, -0.557706, -0.172177 }, { 0.000630461, 0.000733287, 0.00171092 } },
	{ "brass", "brass", { 0.0301974, 0.0223812, 0.0139381 }, { 0.0557826, 0.0376687, 0.0775998 }, { 0.0002028, 0.000258468, 0.00096108 }, { 0.362322, 0.401593, 0.776606 }, { 0.639886, 0.12354, 0.0197853 }, { -0.360114, -0.87646, -0.0919344 }, { 3517.61, 2612.49, 331.815 }, { 3.64061, 2.87206, 0.529487 }, { 1.01146e-07, 9.83073e-08, 5.48819e-08 }, { 312.802, 283.431, 183.091 }, { 0.522711, 0.516719, 0.474834 }, { 0.440765, 0.200948, 0.15549 } },
	{ "cherry-235", "cherry-235", { 0.0497502, 0.0211902, 0.0120688 }, { 0.166001, 0.202786, 0.165189 }, { 0.0182605, 0.0277997, 0.0255721 }, { 0.0358348, 0.163231, 0.129135 }, { 0.0713408, 0.0571719, 0.0791809 }, { 0.0200814, 0.00887306, 0.0251675 }, { 54.7448, 33.7294, 37.36 }, { 6.11314, 6.23697, 6.16351 }, { 30.3886, 0.00191869, 0.0495069 }, { 18.8114, 25.8454, 23.3753 }, { 0.816378, 0.387479, 0.522125 }, { 0.00235908, 0.00199931, 0.00155902 } },
	{ "chrome", "chrome", { 0.00697189, 0.00655268, 0.0101854 }, { 0.0930656, 0.041946, 0.104558 }, { 0.00013, 0.0002, 6e-5 }, { 0.3, 0.36, 0.26 }, { 0.256314, 0.819565, 3.22085e-13 }, { -0.743686, -0.180435, -1 }, { 2200, 2150, 4100 }, { 3.8545, 4.44817, 5.40959 }, { 5.30781e-08, 1.04045e-07, 1.8e+10 }, { 354.965, 349.356, 457 }, { 0.526796, 0.528469, 1.00293 }, { 0.802138, 1.29121, 1.06148 } },
	{ "chrome-steel", "chrome-steel", { 0.0206718, 0.0240818, 0.024351 }, { 0.129782, 0.109032, 0.0524555 }, { 5.51292e-05, 3.13288e-05, 4.51944e-05 }, { 0.207979, 0.152758, 0.325431 }, { 1.18818e-12, 2.06813e-11, 0.580895 }, { -0.316807, -0.265326, -0.419105 }, { 15466.8, 28628.9, 16531.2 }, { 12.8988, 68.7898, 4.68237 }, { 1e+10, 1e+10, 44.7025 }, { 510, 665, 618.155 }, { 1.20035, 1.2003, 0.579562 }, { 0.552277, 0.448267, 1.12869 } },
	{ "colonial-maple-223", "colonial-maple-223", { 0.100723, 0.0356306, 0.0162408 }, { 0.059097, 0.0661341, 0.11024 }, { 0.0197628, 0.0279336, 0.0621265 }, { 0.0311867, 0.112022, 0.344348 }, { 0.0576683, 0.0617498, 0.0479061 }, { -0.0503364, -0.0382196, -0.0382636 }, { 50.7952, 34.6835, 14.201 }, { 6.03342, 6.01053, 4.56588 }, { 18.9034, 0.087246, 1.03757e-07 }, { 18.3749, 21.6159, 26.9347 }, { 0.801289, 0.544191, -0.139944 }, { 0.00229125, 0.00193164, 0.0013723 } },
	{ "color-changing-paint1", "", { 0.00513496, 0.00500415, 0.00296872 }, { 1.53167, 0.430731, 1.48308 }, { 0.00320129, 0.0023053, 0.0329464 }, { 0.167301, 0.100003, 1.04868 }, { 0.0208525, 0.109366, 0.0553211 }, { 0.0101418, 0.0801275, 0.0357004 }, { 279.058, 407.716, 9.71061 }, { 6.73768, 7.19333, 1.04729 }, { 2.66398e+08, 8.98323e+10, 4.16363e-08 }, { 36.7257, 27.8164, 45.5388 }, { 1.01208, 1.19549, 0.131148 }, { 0.164144, 0.0515527, 0.116602 } },
	{ "color-changing-paint2", "", { 0.00463528, 0.00544054, 0.0070818 }, { 1.35172, 1.47838, 1.29831 }, { 0.0279961, 0.0267135, 0.0257468 }, { 0.720154, 0.717648, 0.694662 }, { 0.019073, 0.00825302, 0.0301024 }, { -0.0181694, -0.0198592, -0.000522292 }, { 18.6722, 19.4709, 20.7514 }, { 1.56832, 1.554, 1.60402 }, { 1.19918e-07, 9.95358e-08, 9.69838e-08 }, { 41.295, 42.1716, 42.3986 }, { 0.124794, 0.129237, 0.132024 }, { 0.0489296, 0.0499662, 0.054555 } },
	{ "color-changing-paint3", "", { 0.00305737, 0.00257341, 0.00263616 }, { 0.880793, 0.691268, 0.707821 }, { 0.014742, 0.0135513, 0.0108894 }, { 0.537248, 0.572188, 0.457665 }, { 0.055479, 0.0666783, 0.0585094 }, { 0.0315282, 0.0457554, 0.0365558 }, { 42.1798, 43.3684, 62.0624 }, { 2.06002, 1.82703, 2.54086 }, { 5.0531e-08, 9.99088e-08, 5.726e-08 }, { 50.0948, 52.0344, 54.643 }, { 0.197971, 0.229219, 0.241049 }, { 0.0484699, 0.0342594, 0.0347198 } },
	{ "dark-blue-paint", "dark-blue-paint", { 0.00665057, 0.0139696, 0.0472605 }, { 0.231099, 0.18931, 0.12528 }, { 0.130681, 0.112103, 0.0629285 }, { 0.238562, 0.17, 0.0157371 }, { 0.112486, 0.106018, 0.139316 }, { 0.0662142, 0.0682314, 0.10607 }, { 8.0828, 9.4312, 16.8231 }, { 4.05526, 4.09607, 4.54368 }, { 0.00104743, 0.0035546, 0.10813 }, { 14.347, 13.928, 13.5626 }, { 0.11894, 0.226752, 0.524826 }, { 0.000406243, 0.000291986, 0.000555778 } },
	{ "dark-red-paint", "dark-red-paint", { 0.237125, 0.0365577, 0.0106149 }, { 0.227405, 0.111055, 0.150433 }, { 0.67372, 0.0706716, 0.20328 }, { 1.38127, 1.4345e-13, 0.27723 }, { 1.56101e-11, 0.0779523, 0.136566 }, { -0.148315, 0.0560349, 0.0961999 }, { 3.95965, 15.1861, 5.71355 }, { 2.27942, 4.30359, 3.20838 }, { 0.0225467, 0.113039, 0.00667899 }, { 7.45998, 12.679, 10.9071 }, { 0.119045, 0.521154, 0.166378 }, { 0.00155397, 0.000508377, 0.000218865 } },
	{ "dark-specular-fabric", "", { 0.0197229, 0.00949167, 0.00798414 }, { 0.556218, 0.401495, 0.378651 }, { 0.140344, 0.106541, 0.166715 }, { 0.249059, 0.177611, 0.434167 }, { 0.0351133, 0.0387177, 0.0370533 }, { 0.0243153, 0.0293178, 0.0264913 }, { 7.60492, 9.81673, 6.19307 }, { 3.93869, 4.23097, 4.3775 }, { 0.00122421, 0.00238545, 8.47126e-06 }, { 13.889, 14.5743, 17.2049 }, { 0.114655, 0.210179, -0.227628 }, { 0.00158681, 0.000974676, 0.000638865 } },
	{ "delrin", "delrin", { 0.272703, 0.249805, 0.220642 }, { 0.536593, 0.727886, 0.64011 }, { 0.176535, 0.344018, 0.208011 }, { 0.762088, 0.823603, 0.85976 }, { 0.0398465, 0.0430719, 0.0290162 }, { -0.0121332, -0.0255915, -0.0305128 }, { 5.02017, 3.78669, 4.39454 }, { 3.44325, 3.49289, 3.48886 }, { 9.20654e-08, 0.000215686, 1.06338e-07 }, { 21.2534, 12.5556, 20.455 }, { -0.438862, -0.184248, -0.478699 }, { 0.00493658, 0.00655617, 0.00422698 } },
	{ "fruitwood-241", "fruitwood-241", { 0.0580445, 0.0428667, 0.0259801 }, { 0.203894, 0.233494, 0.263882 }, { 0.00824986, 0.0534794, 0.0472951 }, { 0.160382, 1.07206, 0.768335 }, { 0.00129482, 0.00689891, 0.01274 }, { -0.0211778, -0.0140791, -0.00665974 }, { 110.054, 6.98485, 11.6203 }, { 6.74678, 1.31986, 1.76021 }, { 162.121, 1.18843e-07, 1.0331e-07 }, { 32.6966, 36.752, 34.371 }, { 0.765181, 0.0514459, 0.0106852 }, { 0.00493613, 0.00460352, 0.00373027 } },
	{ "gold-metallic-paint2", "ch-ball-gold-metallic2", { 0.0796008, 0.0538361, 0.0649523 }, { 0.633627, 1.77116, 0.0564028 }, { 0.00376608, 0.00871206, 0.000572055 }, { 0.415684, 0.368424, 0.623038 }, { 0.0343265, 0.00330259, 4.15759e-12 }, { -0.00929705, -0.0219437, -0.0596944 }, { 181.769, 85.527, 795.042 }, { 2.77279, 3.50114, 0.985212 }, { 4.78294e-08, 9.8399e-08, 4.45358e-08 }, { 83.871, 57.2291, 216.537 }, { 0.365256, 0.276734, 0.491264 }, { 0.117145, 0.17457, 0.0754722 } },
	{ "gold-metallic-paint3", "ch-ball-gold-metallic", { 0.0579212, 0.0416649, 0.0271208 }, { 0.0729896, 0.0597695, 0.037684 }, { 0.00146432, 0.00156513, 0.000977438 }, { 0.529437, 0.551234, 0.504486 }, { 1.84643e-14, 5.7212e-15, 8.68546e-13 }, { -1, -1, -0.648897 }, { 382.887, 345.197, 593.725 }, { 4.97078e-07, 1.12888e-06, 4.15366e-07 }, { 3.98586e-07, 3.69533e-07, 16.0596 }, { 22.0196, 22.6462, 11.7126 }, { -0.634103, -0.588056, 0.578355 }, { 0.116233, 0.0897222, 0.0711514 } },
	{ "gold-metallic-paint", "metallic-gold", { 0.0178625, 0.00995704, 0.00335044 }, { 0.17127, 0.120714, 0.115473 }, { 0.127954, 0.127825, 0.109623 }, { 0.781093, 0.795517, 0.661313 }, { 1.00776e-12, 1.22243e-15, 6.18432e-13 }, { -1, -1, -0.334497 }, { 5.90039, 5.83556, 7.14985 }, { 2.7548, 2.7057, 2.94348 }, { 9.46481e-08, 1.06951e-07, 1.10733e-07 }, { 23.8811, 23.9059, 24.2972 }, { -0.303345, -0.293778, -0.267019 }, { 0.0095612, 0.00483452, 0.00145599 } },
	{ "gold-paint", "gold-paint", { 0.147708, 0.0806975, 0.033172 }, { 0.160592, 0.217282, 0.236425 }, { 0.122506, 0.108069, 0.12187 }, { 0.795078, 0.637578, 0.936117 }, { 9.16095e-12, 1.81225e-12, 0.0024589 }, { -0.596835, -0.331147, -0.140729 }, { 5.98176, 7.35539, 5.29722 }, { 2.64832, 3.04253, 2.3013 }, { 9.3111e-08, 8.80143e-08, 9.65288e-08 }, { 24.3593, 24.4037, 25.3623 }, { -0.284195, -0.277297, -0.245352 }, { 0.00313716, 0.00203922, 0.00165683 } },
	{ "gray-plastic", "gray-plastic", { 0.103233, 0.104428, 0.0983734 }, { 0.494656, 0.517207, 0.52772 }, { 0.00758705, 0.00848095, 0.00887135 }, { 0.557908, 0.556548, 0.545887 }, { 0.0428175, 0.0438899, 0.0569098 }, { 0.0208304, 0.0221893, 0.0375763 }, { 75.6274, 68.3098, 66.5689 }, { 1.72111, 1.76597, 1.8406 }, { 1.06783e-07, 5.31845e-08, 6.53296e-08 }, { 65.6538, 63.2018, 61.5805 }, { 0.308907, 0.283768, 0.280211 }, { 0.0512156, 0.050844, 0.0505763 } },
	{ "grease-covered-steel", "", { 0.0196306, 0.0200926, 0.0187026 }, { 0.0433721, 0.0311621, 0.0326401 }, { 0.00019081, 0.000173919, 0.000217638 }, { 0.164569, 0.141125, 0.219421 }, { 7.12672e-13, 1.06789e-14, 1.61131e-13 }, { -1, -1, -1 }, { 4655.4, 5210.63, 3877.45 }, { 17.2827, 29.2523, 8.54352 }, { 1e+38, 1e+38, 1e+38 }, { 103.46, 105.572, 99.1097 }, { 1.20022, 1.20678, 1.19961 }, { 0.261184, 0.20789, 0.210192 } },
	{ "green-acrylic", "acrylic-green", { 0.0176527, 0.0761863, 0.0432331 }, { 0.0517555, 0.15899, 0.0754193 }, { 0.000185443, 5.19959e-05, 7.95188e-05 }, { 0.288191, 0.170979, 0.145492 }, { 0.418137, 0.0486445, 0.170328 }, { 0.330999, 0.0294954, 0.110218 }, { 4223.98, 16977.4, 11351.3 }, { 5.6545, 26.1353, 1.54193e+07 }, { 1e+38, 1e+38, 4.23355e+06 }, { 172.888, 204.112, 403.267 }, { 0.995381, 1.19417, 0.646645 }, { 0.138255, 0.168884, 0.116563 } },
	{ "green-fabric", "fabric-green", { 0.0511324, 0.0490447, 0.0577457 }, { 0.043898, 0.108081, 0.118528 }, { 0.0906425, 0.14646, 0.125546 }, { 0.199121, 0.21946, 0.130311 }, { 0.117671, 0.0797822, 0.0840896 }, { 0.107501, 0.0628391, 0.0668466 }, { 11.1681, 7.43507, 8.69777 }, { 4.65909, 3.72793, 3.72472 }, { 0.00055264, 0.00331292, 0.0119365 }, { 16.8276, 12.7802, 12.1538 }, { 0.153958, 0.17378, 0.291679 }, { 0.000577582, 0.000635885, 0.000729507 } },
	{ "green-latex", "fabric-green-latex", { 0.0885476, 0.13061, 0.0637004 }, { 0.177041, 0.16009, 0.101365 }, { 0.241826, 0.21913, 0.2567 }, { 0.175925, 0.162514, 0.326958 }, { 0.0213854, 0.0498004, 0.0677643 }, { -0.0864353, -0.0518848, -0.00668045 }, { 5.17117, 5.55867, 4.84094 }, { 2.61304, 2.75985, 2.84236 }, { 0.0417628, 0.0352726, 0.0128547 }, { 8.4496, 8.91691, 9.57998 }, { 0.296393, 0.295144, 0.18077 }, { 0.000875063, 0.00103454, 0.000568252 } },
	{ "green-metallic-paint2", "ch-ball-green-metallic", { 0.00536389, 0.0147585, 0.0072232 }, { 0.0553207, 0.0656441, 0.0608999 }, { 0.00131834, 0.00140737, 0.00171711 }, { 0.436009, 0.57969, 0.54703 }, { 0.291808, 0.194127, 0.199549 }, { 0.126008, -0.0681937, 0.0148264 }, { 493.87, 362.85, 317.92 }, { 4.96631e-08, 3.41674e-06, 8.72658e-07 }, { 15.5833, 3.60732e-07, 3.78431e-07 }, { 8.31751, 25.0069, 21.7501 }, { 0.561587, -0.431953, -0.657077 }, { 0.019819, 0.0506927, 0.0279691 } },
	{ "green-metallic-paint", "green-metallic-paint", { 0.00368935, 0.0155555, 0.022272 }, { 0.185621, 0.436002, 0.322925 }, { 0.131402, 0.146271, 0.154061 }, { 1.28366, 0.865104, 0.944013 }, { 0.12483, 0.0443223, 0.0955612 }, { 0.028874, -0.118581, -0.139373 }, { 3.62502, 5.13725, 4.71161 }, { 1.97389, 2.72437, 2.63239 }, { 5.5607e-08, 1.1239e-07, 9.76666e-08 }, { 27.5467, 23.1688, 23.2911 }, { -0.204638, -0.326555, -0.334024 }, { 0.00120994, 0.00217841, 0.00218444 } },
	{ "green-plastic", "green-plastic", { 0.015387, 0.0851675, 0.0947402 }, { 0.0607427, 0.156977, 0.125155 }, { 0.000302146, 0.00197038, 0.000690284 }, { 0.373134, 0.751741, 0.505735 }, { 0.116395, 0.0388464, 0.0476683 }, { 0.0123694, -0.0038785, 0.00864048 }, { 2329.55, 181.586, 833.886 }, { 3.39581, 0.000128411, 3.06753e-07 }, { 6.38801e-08, 2.49751e-07, 5.82332e-07 }, { 260.135, 31.8416, 26.4101 }, { 0.510576, -0.156086, -0.337325 }, { 0.0411204, 0.0490786, 0.0385724 } },
	{ "hematite", "hematite", { 0.00948374, 0.0117628, 0.00985037 }, { 0.0705694, 0.118965, 0.115059 }, { 0.000908552, 0.000601576, 0.00184248 }, { 0.515183, 0.498157, 0.73351 }, { 0.235045, 0.0609324, 0.00720315 }, { -0.120178, -0.096611, -0.235146 }, { 626.075, 967.655, 201.871 }, { 3.94735e-07, 3.54792e-07, 0.000101676 }, { 5.04376e-07, 16.4874, 2.78427e-07 }, { 24.7601, 14.7537, 31.4869 }, { -0.431506, 0.577853, -0.162174 }, { 0.0598105, 0.0605776, 0.0805566 } },
	{ "ipswich-pine-221", "ipswich-pine-221", { 0.0560746, 0.0222518, 0.0105117 }, { 0.0991995, 0.106719, 0.110343 }, { 0.014258, 0.0178759, 0.0188163 }, { 0.0625943, 0.113994, 0.118296 }, { 1.55288e-13, 6.595e-12, 3.97788e-13 }, { -0.0675784, -0.0696373, -0.0703103 }, { 68.7248, 53.3521, 50.6148 }, { 6.36482, 6.37738, 6.35839 }, { 111.962, 1.68378, 0.850892 }, { 20.956, 23.9475, 24.0837 }, { 0.842456, 0.66795, 0.641354 }, { 0.00364018, 0.00320748, 0.00300736 } },
	{ "light-brown-fabric", "", { 0.0612259, 0.0263619, 0.0187761 }, { 3.65487e-12, 9.7449e-12, 4.13685e-12 }, { 0.147778, 0.137639, 0.13071 }, { 0.188292, 0.188374, 0.189026 }, { 1, 1, 1 }, { 0, 0, 0 }, { 7.46192, 7.90085, 8.23842 }, { 3.59765, 3.74493, 3.85476 }, { 0.00660661, 0.0049581, 0.00395337 }, { 12.0657, 12.6487, 13.0977 }, { 0.221546, 0.213332, 0.206743 }, { 0.00126964, 0.000811029, 0.000670968 } },
	{ "light-red-paint", "paint-light-red", { 0.391162, 0.0458387, 0.0059411 }, { 0.522785, 0.144057, 0.214076 }, { 0.854048, 0.0627781, 0.205965 }, { 1.40232, 0.101554, 0.73691 }, { 0.0346142, 0.0655657, 0.0725632 }, { -0.0688995, 0.0400933, 0.0252314 }, { 4.68691, 16.2538, 4.71913 }, { 1.72476, 4.88106, 3.96894 }, { 0.0969212, 0.0124952, 1.02202e-07 }, { 5.57044, 16.1561, 19.9884 }, { 0.249637, 0.388724, -0.505514 }, { 0.00271456, 0.000889762, 0.000505975 } },
	{ "maroon-plastic", "maroon-plastic", { 0.189951, 0.0353828, 0.0321504 }, { 0.127693, 0.100703, 0.115731 }, { 0.00160715, 0.00110827, 0.00100127 }, { 0.684406, 0.645917, 0.569111 }, { 0.0479368, 0.0624437, 0.0921161 }, { -0.0134224, 0.00888653, 0.0423901 }, { 257.032, 398.91, 515.363 }, { 4.57088e-05, 2.62643e-05, 3.15876e-06 }, { 3.05979e-07, 2.76112e-07, 3.67984e-07 }, { 29.9497, 31.8569, 27.4899 }, { -0.211122, -0.159036, -0.309239 }, { 0.0285283, 0.0277773, 0.0286259 } },
	{ "natural-209", "natural-209", { 0.0961753, 0.0349012, 0.0121752 }, { 0.0781649, 0.0898869, 0.111321 }, { 0.0137282, 0.0154247, 0.0233645 }, { 0.0491415, 0.0673163, 0.149331 }, { 0.0522205, 0.0456428, 0.0279324 }, { -0.0277529, -0.0308045, -0.0492147 }, { 71.8965, 63.4418, 40.2168 }, { 6.3511, 6.32932, 6.32009 }, { 224.515, 48.774, 0.0171617 }, { 20.2337, 21.0421, 25.7335 }, { 0.874926, 0.811948, 0.483785 }, { 0.00465507, 0.00387214, 0.00295941 } },
	{ "neoprene-rubber", "neoprene-rubber", { 0.259523, 0.220477, 0.184871 }, { 0.275058, 0.391429, 0.0753145 }, { 0.143818, 0.207586, 0.0764912 }, { 0.770284, 0.774203, 0.700644 }, { 0.113041, 0.110436, 0.16895 }, { 0.060346, 0.0565499, 0.0788468 }, { 5.56845, 4.61088, 8.84784 }, { 3.02411, 3.82334, 2.38942 }, { 5.3042e-08, 1.01087e-07, 6.70643e-08 }, { 23.2109, 20.1103, 28.3099 }, { -0.379072, -0.500903, -0.156114 }, { 0.00218754, 0.00314666, 0.00180987 } },
	{ "nickel", "nickel", { 0.0144009, 0.0115339, 0.00989042 }, { 0.157696, 0.293022, 0.450103 }, { 0.00556292, 0.00627392, 0.00660563 }, { 0.171288, 0.168324, 0.161023 }, { 2.21884, 1.61986, 0.931645 }, { 1.21884, 1.22103, 0.698939 }, { 160.907, 143.234, 136.933 }, { 6.76252, 6.76447, 6.76278 }, { 16179.3, 3451.04, 3453.27 }, { 36.1378, 35.1378, 33.614 }, { 0.846825, 0.821209, 0.830948 }, { 0.132612, 0.177506, 0.155358 } },
	{ "nylon", "nylon", { 0.204199, 0.211192, 0.19234 }, { 0.156797, 0.303324, 0.236394 }, { 0.0250344, 0.0436802, 0.0421753 }, { 0.528875, 0.617086, 0.620808 }, { 0.240279, 0.117115, 0.127421 }, { 0.219191, 0.0940848, 0.101339 }, { 26.4669, 14.8484, 15.213 }, { 2.32486, 2.22618, 2.18663 }, { 1.17309e-07, 5.34378e-08, 6.68165e-08 }, { 39.9707, 33.9774, 34.3152 }, { 0.11793, -0.0190508, -0.00222626 }, { 0.00830775, 0.00885046, 0.00645556 } },
	{ "orange-paint", "orange-paint", { 0.368088, 0.147113, 0.00692426 }, { 0.524979, 0.116386, 0.199437 }, { 0.818115, 0.064743, 0.229391 }, { 1.44385, 0.0709512, 0.483597 }, { 6.92565e-13, 0.106161, 0.102279 }, { -0.174318, 0.0934385, 0.0625648 }, { 4.57466, 16.0185, 4.96427 }, { 1.84547, 4.70387, 3.6232 }, { 0.072629, 0.0299825, 0.000333551 }, { 5.96872, 14.9466, 13.2194 }, { 0.222125, 0.438216, -0.0759733 }, { 0.00178442, 0.000789226, 0.000301814 } },
	{ "pearl-paint", "pearl-paint", { 0.181967, 0.159068, 0.143348 }, { 0.105133, 0.0928717, 0.0802367 }, { 0.0724063, 0.0808503, 0.0596139 }, { 0.194454, 0.203296, 0.091958 }, { 0.168966, 0.297431, 0.401185 }, { -0.831034, -0.702569, -0.226546 }, { 13.6335, 12.3113, 17.1256 }, { 5.07525, 4.9097, 4.92468 }, { 0.000252814, 0.000258641, 0.0180128 }, { 18.9547, 18.2047, 16.1366 }, { 0.156731, 0.135808, 0.415562 }, { 0.00244402, 0.00167924, 0.0012928 } },
	{ "pickled-oak-260", "pickled-oak-260", { 0.181735, 0.14142, 0.125486 }, { 0.0283411, 0.0296418, 0.025815 }, { 0.0105853, 0.0102771, 0.0101188 }, { 2.31337e-14, 2.35272e-14, 1.99762e-14 }, { 5.38184e-13, 2.15933e-13, 3.55496e-12 }, { -0.309259, -0.291046, -0.329625 }, { 95.4759, 98.3089, 99.831 }, { 6.37433, 6.39032, 6.39857 }, { 4641.15, 5970.15, 6818.76 }, { 18.1707, 18.2121, 18.2334 }, { 1.00563, 1.01274, 1.01645 }, { 0.00336548, 0.0032891, 0.00330829 } },
	{ "pink-fabric2", "", { 0.24261, 0.0829238, 0.0751196 }, { 0.161823, 0.0591236, 0.00907967 }, { 0.220011, 0.148623, 0.111966 }, { 6.44517e-12, 3.24286e-13, 3.83556e-11 }, { 0.242032, 0.456181, 3.11925 }, { -0.0582931, 0.295844, 2.84268 }, { 5.66376, 7.80657, 9.98941 }, { 2.43742, 3.05804, 3.53392 }, { 0.111404, 0.079067, 0.0738938 }, { 7.47166, 9.23588, 10.5653 }, { 0.408014, 0.422779, 0.448866 }, { 0.00145188, 0.000903357, 0.00103814 } },
	{ "pink-fabric", "fabric-pink", { 0.270553, 0.223977, 0.240993 }, { 0.299998, 0.418074, 4.07112e-13 }, { 0.787023, 0.234345, 0.17346 }, { 1.77629, 4.64947e-15, 0.203846 }, { 0.121124, 0.0271041, 1 }, { 0.0151107, 0.00273969, 0 }, { 4.71724, 5.3941, 6.55218 }, { 2.23471, 2.34423, 3.31887 }, { 0.0266324, 0.11955, 0.00966123 }, { 7.20596, 7.20326, 11.0696 }, { 0.130136, 0.407562, 0.22282 }, { 0.00129262, 0.00124216, 0.00262093 } },
	{ "pink-felt", "fabric-pink-felt", { 0.259533, 0.192978, 0.185581 }, { 0.359813, 0.533498, 0.0390541 }, { 0.46679, 0.203504, 0.314933 }, { 0.663613, 1.70267e-13, 0.778919 }, { 0.0530603, 0.0108612, 1.12976 }, { -0.124124, -0.0455941, 1.05206 }, { 3.6266, 6.02293, 3.92458 }, { 2.22042, 2.55548, 3.68735 }, { 0.0333497, 0.102524, 8.10122e-05 }, { 7.40195, 7.81026, 13.5108 }, { 0.185696, 0.409228, -0.236166 }, { 0.000823299, 0.00122612, 0.00141316 } },
	{ "pink-jasper", "pink-jasper", { 0.226234, 0.138929, 0.110785 }, { 0.0846118, 0.0984038, 0.078693 }, { 0.00223592, 0.00203213, 0.0013737 }, { 0.53138, 0.4995, 0.436111 }, { 0.172825, 0.129714, 0.154397 }, { 0.0911362, 0.0599408, 0.0851167 }, { 253.003, 292.741, 474.2 }, { 1.64528, 1.86649, 5.29965e-08 }, { 5.71153e-08, 6.24765e-08, 15.5452 }, { 110.886, 113.953, 8.17872 }, { 0.416268, 0.422653, 0.561605 }, { 0.0166667, 0.0194625, 0.0142451 } },
	{ "pink-plastic", "pink-plastic", { 0.354572, 0.0905002, 0.0696372 }, { 0.0316585, 0.0444153, 6.36158e-15 }, { 0.0566727, 0.0369011, 0.142628 }, { 0.215149, 0.0621896, 0.189286 }, { 0.317781, 0.0925394, 1 }, { 0.266863, 0.0652821, 0 }, { 16.788, 27.3043, 7.67451 }, { 5.62766, 5.4919, 3.67439 }, { 1.71214e-05, 0.187553, 0.00561921 }, { 23.0822, 17.5369, 12.3737 }, { 0.0833775, 0.577298, 0.215985 }, { 0.00219041, 0.000905537, 0.00219415 } },
	{ "polyethylene", "polyethylene", { 0.228049, 0.239339, 0.240326 }, { 0.0420869, 0.134269, 0.0867928 }, { 0.0472725, 0.260465, 0.0719615 }, { 2.09548e-13, 0.743064, 2.27838e-13 }, { 0.489907, 0.316434, 0.181688 }, { 0.424297, 0.244434, 0.0772152 }, { 22.178, 4.24558, 14.9332 }, { 4.92064, 4.37149, 4.27422 }, { 0.330549, 9.78063e-07, 0.109488 }, { 14.3439, 17.178, 12.5991 }, { 0.61045, -0.470479, 0.517652 }, { 0.00285999, 0.00147486, 0.00191841 } },
	{ "polyurethane-foam", "polyurethane-foam", { 0.0898318, 0.0428583, 0.0340984 }, { 4.0852e-12, 8.0217e-14, 4.05682e-14 }, { 0.154984, 0.142104, 0.139418 }, { 0.188586, 0.188095, 0.188124 }, { 1, 1, 1 }, { 0, 0, 0 }, { 7.18433, 7.70043, 7.81983 }, { 3.50087, 3.67779, 3.7174 }, { 0.0079259, 0.00567317, 0.00525063 }, { 11.6922, 12.3795, 12.537 }, { 0.226848, 0.217315, 0.215132 }, { 0.00119032, 0.000823112, 0.000753708 } },
	{ "pure-rubber", "pure-rubber", { 0.284259, 0.251873, 0.223824 }, { 0.542899, 0.598765, 0.142162 }, { 0.62899, 0.41413, 0.0693873 }, { 1.13687, 0.452063, 0.262952 }, { 0.0185379, 0.00420901, 0.0353132 }, { -0.0265343, -0.035723, 0.0236736 }, { 3.75577, 3.75969, 13.6137 }, { 2.17187, 2.13584, 5.56595 }, { 0.0311422, 0.0526357, 7.8942e-07 }, { 7.14885, 7.0586, 23.9264 }, { 0.152489, 0.245615, -0.0979145 }, { 0.00177602, 0.00171837, 0.0007762 } },
	{ "purple-paint", "purple-paint", { 0.290743, 0.0347118, 0.0339802 }, { 0.301308, 0.205258, 0.280717 }, { 0.0339173, 0.00839425, 0.0378011 }, { 0.703279, 0.20083, 0.797376 }, { 0.0654958, 0.0786512, 0.0696424 }, { 0.0299494, 0.0467651, 0.0373649 }, { 16.3705, 104.698, 13.2101 }, { 1.72814, 6.77684, 1.54679 }, { 1.10548e-07, 0.177314, 6.66798e-08 }, { 38.0847, 41.1712, 38.2781 }, { 0.0789738, 0.576929, 0.0587144 }, { 0.00296394, 0.00277256, 0.00225685 } },
	{ "pvc", "pvc", { 0.0322978, 0.0357449, 0.0403426 }, { 0.28767, 0.317369, 0.310067 }, { 0.0171547, 0.0176681, 0.0213663 }, { 0.769726, 0.730637, 0.797555 }, { 0.0289552, 0.026258, 0.0281305 }, { 0.00815373, 0.00651989, 0.00714485 }, { 25.8231, 26.9026, 20.5927 }, { 1.20095, 1.31883, 1.23655 }, { 1.02768e-07, 1.09767e-07, 5.97648e-08 }, { 51.5975, 50.0184, 48.1953 }, { 0.218049, 0.208814, 0.173894 }, { 0.00471012, 0.00620207, 0.00575443 } },
	{ "red-fabric2", "", { 0.155216, 0.0226757, 0.0116884 }, { 1.80657e-15, 5.51946e-13, 1.35221e-15 }, { 0.16689, 0.135884, 0.128307 }, { 0.184631, 0.18856, 0.1883 }, { 1, 1, 1 }, { 0, 0, 0 }, { 6.78759, 7.98303, 8.36701 }, { 3.33819, 3.77225, 3.89063 }, { 0.0112096, 0.00468671, 0.00372538 }, { 11.0557, 12.7596, 13.2393 }, { 0.240935, 0.21164, 0.206003 }, { 0.00263226, 0.00060129, 0.000538747 } },
	{ "red-fabric", "fabric-red", { 0.201899, 0.0279008, 0.0103965 }, { 0.168669, 0.0486346, 0.040485 }, { 0.324447, 0.228455, 0.109436 }, { 0.787411, 0.821197, 0.279212 }, { 0.0718348, 0.0644687, 0.0206123 }, { -0.0585917, -0.0062547, -0.050402 }, { 3.88129, 4.32067, 9.16355 }, { 3.59825, 3.93046, 4.66379 }, { 0.000130047, 1.04152e-07, 4.59182e-05 }, { 13.0776, 19.649, 17.8852 }, { -0.209387, -0.530789, -0.0242035 }, { 0.00103676, 0.000248481, 0.000220698 } },
	{ "red-metallic-paint", "ch-ball-red-metallic", { 0.0380897, 0.00540095, 0.00281156 }, { 0.0416724, 0.07642, 0.108438 }, { 0.00133258, 0.00106883, 0.00128863 }, { 0.693854, 0.52857, 0.539477 }, { 2.45718e-16, 0.0598671, 0.0633332 }, { -1, -0.08904, -0.0114056 }, { 300.371, 521.418, 425.982 }, { 6.45857e-05, 6.3446e-07, 8.51754e-07 }, { 2.75773e-07, 4.05125e-07, 3.41703e-07 }, { 33.0213, 24.3646, 23.5793 }, { -0.121415, -0.469753, -0.532083 }, { 0.287672, 0.0467824, 0.0286558 } },
	{ "red-phenolic", "red-phenolic", { 0.165227, 0.0256259, 0.00935644 }, { 0.240561, 0.360634, 0.475777 }, { 0.0052844, 0.00467439, 0.00613717 }, { 0.568938, 0.509763, 0.575762 }, { 0.156419, 0.0972193, 0.069671 }, { 0.0752589, 0.0444558, 0.0266428 }, { 104.336, 128.839, 89.6357 }, { 1.57629, 1.92067, 1.57462 }, { 6.38793e-08, 6.71457e-08, 5.25172e-08 }, { 77.4088, 79.3795, 73.0158 }, { 0.343482, 0.352913, 0.324912 }, { 0.0680783, 0.0859239, 0.0936197 } },
	{ "red-plastic", "red-plastic", { 0.247569, 0.049382, 0.0175621 }, { 0.406976, 0.151478, 0.176348 }, { 0.28723, 0.0572489, 0.0624682 }, { 0.939617, 0.0851973, 0.0701483 }, { 0.10036, 0.178468, 0.149441 }, { 0.0512697, 0.13191, 0.102958 }, { 3.80564, 17.8374, 16.5631 }, { 4.37031, 4.96092, 4.75924 }, { 1.03822e-07, 0.0236845, 0.0321546 }, { 18.6088, 16.1371, 15.1243 }, { -0.609047, 0.43547, 0.445882 }, { 0.00134715, 0.000536577, 0.00047104 } },
	{ "red-specular-plastic", "red-plastic-specular", { 0.252589, 0.0397665, 0.0185317 }, { 0.0139957, 0.0343278, 0.0527973 }, { 6.01746e-05, 8.07327e-05, 0.000205705 }, { 0.174569, 0.202455, 0.390522 }, { 0.441328, 0.179378, 0.150221 }, { 0.191312, 0.0691835, 0.0760833 }, { 14623.1, 10620.2, 3332.98 }, { 22.6631, 12.9299, 3.06408 }, { 1e+38, 1e+38, 6.36761e-08 }, { 187.799, 162.245, 315.067 }, { 1.19812, 1.20136, 0.52102 }, { 0.0620224, 0.0541448, 0.0836154 } },
	{ "silicon-nitrade", "silicon-nitrade", { 0.0141611, 0.0115865, 0.00842477 }, { 0.0710113, 0.0670906, 0.015769 }, { 6.40406e-05, 0.000138867, 0.00224354 }, { 0.159422, 0.283527, 0.734323 }, { 0.0516164, 0.10318, 2.36643 }, { -0.0277792, -0.0531505, 1.36643 }, { 13926.7, 5668.82, 167.824 }, { 31.0579, 5.63563, 8.52534e-05 }, { 1e+38, 1e+38, 2.96855e-07 }, { 176.231, 148.745, 29.0282 }, { 1.20757, 1.14062, -0.241968 }, { 0.0712429, 0.0720965, 0.11404 } },
	{ "silver-metallic-paint2", "", { 0.0554792, 0.0573803, 0.0563376 }, { 0.121338, 0.115673, 0.10966 }, { 0.029859, 0.0303706, 0.0358666 }, { 0.144097, 0.104489, 0.158163 }, { 1.03749e-14, 3.52034e-15, 4.41778e-12 }, { -1, -1, -1 }, { 31.9005, 32.1514, 26.5685 }, { 6.08248, 5.89319, 5.95156 }, { 0.00761403, 0.0839948, 0.00138703 }, { 23.5891, 20.6786, 23.3297 }, { 0.435405, 0.54017, 0.34665 }, { 0.021384, 0.0187248, 0.0175762 } },
	{ "silver-metallic-paint", "metallic-silver", { 0.0189497, 0.0205686, 0.0228822 }, { 0.173533, 0.168901, 0.165266 }, { 0.037822, 0.038145, 0.0381908 }, { 0.165579, 0.162955, 0.160835 }, { 5.66903e-12, 1.65276e-14, 4.28399e-14 }, { -1, -1, -1 }, { 25.1551, 24.9957, 25.0003 }, { 5.92591, 5.90225, 5.89007 }, { 0.000684235, 0.000840795, 0.000995926 }, { 23.4725, 23.1898, 23.0144 }, { 0.310575, 0.317996, 0.324938 }, { 0.0122324, 0.0117672, 0.0113004 } },
	{ "silver-paint", "paint-silver", { 0.152796, 0.124616, 0.113375 }, { 0.30418, 0.30146, 0.283174 }, { 0.110819, 0.105318, 0.0785677 }, { 0.640378, 0.641115, 0.445228 }, { 2.37347e-13, 7.68194e-13, 2.9434e-12 }, { -0.350607, -0.355433, -0.359297 }, { 7.21531, 7.46519, 10.8289 }, { 3.06016, 2.97652, 3.78287 }, { 9.71666e-08, 1.09342e-07, 9.82336e-08 }, { 24.1475, 24.5083, 25.6757 }, { -0.281384, -0.257633, -0.200968 }, { 0.00361589, 0.00384995, 0.00444223 } },
	{ "special-walnut-224", "special-walnut-224", { 0.0121712, 0.00732998, 0.00463072 }, { 0.209603, 0.216118, 0.211885 }, { 0.117091, 0.119932, 0.131119 }, { 0.548899, 0.524858, 0.569425 }, { 0.0808859, 0.0802614, 0.0789982 }, { 0.0327605, 0.0324012, 0.0274637 }, { 7.42314, 7.41578, 6.77215 }, { 3.61532, 3.82778, 3.71096 }, { 1.19182e-07, 1.03098e-07, 1.08004e-07 }, { 22.9756, 22.7318, 22.2976 }, { -0.31158, -0.332257, -0.35462 }, { 0.000934753, 0.000878135, 0.000758769 } },
	{ "specular-black-phenolic", "black-bball", { 0.00212164, 0.00308282, 0.00410253 }, { 0.0881574, 0.0923246, 0.0398117 }, { 0.00119167, 0.000641898, 0.000186605 }, { 0.616914, 0.578026, 0.40121 }, { 0.122486, 0.0907984, 0.13855 }, { 0.0482851, 0.0348218, 0.0481813 }, { 395.656, 781.361, 3615.49 }, { 1.18051e-05, 1.22149, 2.87626 }, { 2.88932e-07, 6.09017e-08, 9.79363e-08 }, { 28.9748, 200.529, 331.336 }, { -0.257695, 0.487197, 0.524707 }, { 0.0311874, 0.0437836, 0.0853211 } },
	{ "specular-blue-phenolic", "blue-bball", { 0.00497564, 0.0138836, 0.032815 }, { 0.1077, 0.0898232, 0.175296 }, { 0.000918571, 0.0010348, 0.00176322 }, { 0.570978, 0.639916, 0.709385 }, { 0.0354139, 0.0488958, 0.023159 }, { -0.0434314, -0.0294427, -0.039938 }, { 558.482, 431.708, 222.547 }, { 3.57872e-06, 2.38985e-05, 6.87316e-05 }, { 3.33356e-07, 3.3336e-07, 3.2851e-07 }, { 28.5289, 32.1362, 30.3802 }, { -0.272244, -0.141224, -0.19005 }, { 0.0324614, 0.0338654, 0.057732 } },
	{ "specular-green-phenolic", "green-bball", { 0.00781782, 0.0259654, 0.0233739 }, { 0.0688449, 0.144658, 0.143654 }, { 0.000307494, 0.0010353, 0.00155331 }, { 0.365481, 0.585805, 0.787512 }, { 0.129429, 0.0443676, 5.423e-12 }, { 0.047932, -0.0136055, -0.0592743 }, { 2313.47, 482.847, 206.646 }, { 3.55591, 5.3306e-06, 0.000279403 }, { 1.03804e-07, 3.47154e-07, 3.16888e-07 }, { 256.625, 28.2901, 38.0691 }, { 0.511901, -0.277015, 0.00513058 }, { 0.0373526, 0.0450242, 0.0470867 } },
	{ "specular-maroon-phenolic", "maroon-bball", { 0.152486, 0.0263216, 0.00802748 }, { 0.0761775, 0.098375, 0.165913 }, { 0.000342958, 0.000605578, 0.00144136 }, { 0.4052, 0.553617, 0.65133 }, { 0.0646024, 0.0116325, 0.037551 }, { -0.0555983, -0.0527264, -0.0283844 }, { 1961.35, 868.119, 306.537 }, { 2.81751, 1.361, 2.43015e-05 }, { 5.70302e-08, 5.89096e-08, 3.27578e-07 }, { 248.591, 203.707, 29.053 }, { 0.506511, 0.488814, -0.242746 }, { 0.0465516, 0.0434438, 0.0545231 } },
	{ "specular-orange-phenolic", "orange-bball", { 0.32771, 0.0540131, 0.00883213 }, { 0.051915, 0.0686764, 0.0489478 }, { 7.91913e-05, 0.000139576, 1.62017e-05 }, { 0.253564, 0.354675, 3.55583e-12 }, { 0.0768695, 0.0496641, 0.0223538 }, { -0.015514, -0.023864, -0.0433847 }, { 10274.5, 5159.16, 61722.9 }, { 7.4207, 3.83488, 8.72207e+06 }, { 1e+38, 1.02522e-07, 2.30124e+07 }, { 170.364, 373.255, 855.811 }, { 1.19475, 0.53082, 0.608547 }, { 0.0405677, 0.0498947, 0.0953035 } },
	{ "specular-red-phenolic", "red-bball", { 0.303563, 0.0354891, 0.00899721 }, { 0.151819, 0.0938022, 0.196935 }, { 0.00117843, 0.00056476, 0.00185124 }, { 0.570146, 0.524406, 0.732785 }, { 0.0287445, 0.0672098, 0.0178433 }, { -0.0386172, 0.0144813, -0.0322308 }, { 438.968, 982.441, 201.328 }, { 2.88657e-06, 8.12418e-07, 9.95633e-05 }, { 3.8194e-07, 3.69661e-07, 3.03179e-07 }, { 26.0046, 30.0091, 31.3125 }, { -0.376008, -0.21799, -0.162894 }, { 0.0307055, 0.0199408, 0.0652829 } },
	{ "specular-violet-phenolic", "violet-bball", { 0.0686035, 0.0181856, 0.0210368 }, { 0.108459, 0.0471612, 0.171691 }, { 0.00123271, 0.000443974, 0.00149517 }, { 0.657484, 0.546753, 0.653065 }, { 0.0403569, 0.121081, 0.035323 }, { -0.0295013, 0.0563904, -0.0275623 }, { 351.208, 1193.45, 294.897 }, { 3.17585e-05, 1.3817, 2.44051e-05 }, { 3.02028e-07, 6.19706e-08, 3.40809e-07 }, { 31.3319, 234.879, 28.7237 }, { -0.168991, 0.500354, -0.252626 }, { 0.033584, 0.0495535, 0.0510203 } },
	{ "specular-white-phenolic", "white-bball", { 0.282896, 0.231703, 0.127818 }, { 0.0678467, 0.0683808, 0.032756 }, { 2.00918e-05, 3.22307e-05, 9.59333e-05 }, { 4.35028e-11, 0.131903, 0.390338 }, { 0.0132613, 0.0416359, 0.285115 }, { -0.1025, -0.0828341, 0.0759745 }, { 49772.5, 28318.4, 7131.79 }, { 4.58788e+06, 1.12234e+07, 3.05696 }, { 4.5813e+06, 1.00372e+07, 6.93295e-08 }, { 767.24, 646.93, 455.184 }, { 0.60999, 0.61958, 0.536789 }, { 0.134616, 0.114414, 0.179676 } },
	{ "specular-yellow-phenolic", "yellow-bball", { 0.309395, 0.135278, 0.0159106 }, { 0.0607659, 0.141526, 0.110839 }, { 0.000200174, 0.00179326, 0.00048094 }, { 0.355381, 0.767285, 0.505029 }, { 0.077039, 0.0175574, 0.038586 }, { -0.0182146, -0.0417745, -0.0114057 }, { 3597.35, 190.656, 1192.15 }, { 3.80702, 0.000181166, 1.69208 }, { 1.01181e-07, 3.0765e-07, 1.27637e-07 }, { 313.8, 34.1523, 221.08 }, { 0.523001, -0.0813947, 0.500075 }, { 0.0451868, 0.0491835, 0.031808 } },
	{ "ss440", "ss440", { 0.0229923, 0.0187037, 0.0153204 }, { 0.127809, 0.14899, 0.0473376 }, { 4.06782e-05, 4.86278e-05, 8.15367e-05 }, { 0.0888931, 0.259135, 0.332727 }, { 2.26479e-12, 5.7223e-13, 0.820577 }, { -0.474442, -0.223179, -0.179423 }, { 23197.8, 16625.7, 9083.17 }, { 1.43395e+07, 7.33446, 4.46921 }, { 1.06287e+07, 1e+38, 9.97795e-08 }, { 556.069, 215.021, 479.875 }, { 0.627834, 1.20059, 0.540075 }, { 0.366952, 0.463741, 0.874071 } },
	{ "steel", "steel", { 0.019973, 0.0127074, 0.0246402 }, { 0.0615275, 0.0469644, 0.0402151 }, { 8.63865e-05, 0.000249576, 5.77865e-05 }, { 0.3729, 0.583679, 0.471406 }, { 0.665193, 0.798139, 0.0115189 }, { -0.334807, -0.201861, -0.988481 }, { 8119.74, 1951.31, 10374.2 }, { 3.40501, 1.12086, 1.88549 }, { 1.25096e-07, 9.51373e-08, 1.20412e-07 }, { 474.668, 313.81, 603.334 }, { 0.539696, 0.519433, 0.545491 }, { 1.11996, 1.28262, 1.42004 } },
	{ "teflon", "teflon", { 0.276442, 0.263098, 0.260294 }, { 1.56924, 1.52804, 1.43859 }, { 0.678586, 0.662167, 0.577852 }, { 1.2402, 1.14126, 1.44077 }, { 9.40662e-13, 3.57656e-11, 1.50208e-11 }, { -0.0492032, -0.0587316, -0.0548028 }, { 3.91532, 3.82834, 3.64944 }, { 2.09871, 2.05044, 2.90315 }, { 0.0378525, 0.0435869, 0.0028469 }, { 6.8689, 6.72229, 9.65585 }, { 0.166634, 0.181605, -0.0468321 }, { 0.00462131, 0.00463132, 0.00408386 } },
	{ "tungsten-carbide", "tungsten-carbide", { 0.0151872, 0.0103016, 0.0123192 }, { 0.0504358, 0.075701, 0.0556673 }, { 6.6122e-05, 7.65809e-05, 4.80196e-05 }, { 0.255291, 0.270824, 0.26732 }, { 5.30357e-14, 5.44537e-12, 1.09586e-11 }, { -1, -1, -1 }, { 12280.7, 10423.3, 16683.5 }, { 7.42976, 6.36068, 6.81602 }, { 1e+38, 1e+38, 1e+38 }, { 185.623, 173.962, 218.975 }, { 1.197, 1.19586, 1.19745 }, { 0.472021, 0.583364, 0.850304 } },
	{ "two-layer-gold", "", { 0.0415046, 0.0312801, 0.0253658 }, { 1.58161, 1.18736, 1.63847 }, { 0.0263104, 0.0293804, 0.0241265 }, { 0.355682, 0.354281, 0.36415 }, { 0.117355, 0.0614942, 0.0447004 }, { 0.0411678, -0.0237579, -0.0100488 }, { 30.4478, 27.5367, 32.7278 }, { 3.91606, 3.99968, 3.79114 }, { 9.38513e-08, 6.62806e-08, 6.79906e-08 }, { 36.9091, 35.6578, 38.4915 }, { 0.0815181, 0.0465613, 0.0922721 }, { 0.165383, 0.160769, 0.168742 } },
	{ "two-layer-silver", "", { 0.0657916, 0.0595705, 0.0581288 }, { 1.55275, 2.00145, 1.93045 }, { 0.0149977, 0.0201665, 0.0225062 }, { 0.382631, 0.35975, 0.361657 }, { 4.93242e-13, 1.00098e-14, 0.0103259 }, { -0.0401315, -0.0395054, -0.0312454 }, { 50.1263, 38.8508, 34.9978 }, { 3.41873, 3.77545, 3.78138 }, { 6.09709e-08, 1.02036e-07, 1.01016e-07 }, { 46.6236, 40.8229, 39.1812 }, { 0.183797, 0.139103, 0.117092 }, { 0.170639, 0.189329, 0.21468 } },
	{ "violet-acrylic", "acrylic-violet", { 0.0599875, 0.023817, 0.0379025 }, { 0.134984, 0.13337, 0.295509 }, { 0.0011295, 0.00126481, 0.00186818 }, { 0.523244, 0.551606, 0.478149 }, { 0.100325, 0.0939057, 0.0663939 }, { 0.0603016, 0.0579001, 0.0353439 }, { 498.715, 424.307, 328.834 }, { 4.79917e-07, 1.39397e-06, 2.03376 }, { 4.27008e-07, 3.4655e-07, 1.05578e-07 }, { 23.5924, 24.3831, 116.7 }, { -0.515008, -0.476979, 0.432199 }, { 0.041711, 0.0411698, 0.0691096 } },
	{ "violet-rubber", "violet-rubber", { 0.223179, 0.0553634, 0.113238 }, { 0.547456, 0.0966027, 0.185463 }, { 0.445092, 0.089421, 0.1543 }, { 0.923481, 0.782897, 0.800907 }, { 0.0518927, 0.113508, 0.094549 }, { -0.0208955, 0.0585063, 0.036751 }, { 3.58673, 7.33539, 5.24404 }, { 2.83007, 2.27788, 3.00009 }, { 0.00374266, 9.65628e-08, 1.13746e-07 }, { 9.67534, 27.1898, 22.3519 }, { -0.0037936, -0.174481, -0.363226 }, { 0.00274248, 0.000994259, 0.00101746 } },
	{ "white-acrylic", "acrylic-white", { 0.314106, 0.300008, 0.263648 }, { 0.015339, 0.0736169, 0.0976209 }, { 0.00233914, 0.00168956, 0.00117717 }, { 0.570501, 0.51337, 0.440543 }, { 2.651, 0.38473, 0.249988 }, { 2.27139, 0.341493, 0.218093 }, { 226.134, 342.423, 548.396 }, { 1.41218, 2.17179e-07, 4.52325e-08 }, { 9.83895e-08, 5.84338e-06, 15.7173 }, { 110.209, 18.8595, 8.78916 }, { 0.419684, -0.61884, 0.56243 }, { 0.0537471, 0.0401864, 0.0400114 } },
	{ "white-diffuse-bball", "white-diffuse-bball", { 0.284726, 0.239199, 0.177227 }, { 0.680658, 0.828508, 0.502128 }, { 0.689731, 0.601478, 0.209084 }, { 1.23085, 0.956833, 0.800488 }, { 1.06465e-11, 7.45307e-19, 0.089058 }, { -0.109723, -0.090129, 0.0545089 }, { 3.94153, 3.65714, 4.53186 }, { 2.04845, 2.06125, 3.72838 }, { 0.0434345, 0.0439287, 1.06879e-07 }, { 6.69695, 6.79535, 20.1466 }, { 0.178983, 0.18898, -0.494774 }, { 0.007671, 0.0072205, 0.00441788 } },
	{ "white-fabric2", "", { 0.10784, 0.102669, 0.113943 }, { 0.0375359, 0.0296317, 0.0364218 }, { 0.0583526, 0.0520763, 0.063905 }, { 9.1836e-12, 5.41897e-13, 3.47708e-12 }, { 0.0528459, 0.0854122, 0.0860966 }, { 0.0578589, 0.0971296, 0.0985717 }, { 18.1669, 20.2291, 16.6809 }, { 4.60692, 4.77932, 4.46485 }, { 0.172556, 0.23927, 0.138435 }, { 13.5006, 13.9648, 13.1165 }, { 0.561147, 0.587176, 0.54153 }, { 0.00113864, 0.00114867, 0.00119185 } },
	{ "white-fabric", "", { 0.290107, 0.219835, 0.160654 }, { 0.230066, 0.156787, 2.19005e-12 }, { 0.479844, 0.196767, 0.162438 }, { 0.669662, 2.22894e-13, 0.194582 }, { 7.17773e-14, 0.0364632, 1 }, { -0.209237, -0.0528699, 0 }, { 3.60967, 6.18732, 6.91147 }, { 2.17095, 2.60737, 3.42418 }, { 0.0379118, 0.0990846, 0.00860197 }, { 7.221, 7.95861, 11.4289 }, { 0.195346, 0.409992, 0.225019 }, { 0.000728813, 0.000806025, 0.00177083 } },
	{ "white-marble", "", { 0.236183, 0.221746, 0.192889 }, { 0.24075, 0.221456, 0.23389 }, { 0.00430204, 0.00429388, 0.00405907 }, { 0.681306, 0.676061, 0.620701 }, { 0.118569, 0.167963, 0.140415 }, { 0.0620704, 0.102032, 0.071968 }, { 103.053, 104.312, 122.169 }, { 1.0396, 1.05866, 1.26372 }, { 1.08288e-07, 6.82995e-08, 5.76023e-08 }, { 88.7766, 89.0261, 88.95 }, { 0.378101, 0.372681, 0.371975 }, { 0.0493991, 0.0689784, 0.0584648 } },
	{ "white-paint", "white-paint", { 0.356024, 0.3536, 0.324889 }, { 4.14785, 0.255488, 0.108438 }, { 0.0654126, 0.0841905, 0.538225 }, { 0.796927, 0.7778, 0.888627 }, { 0.0770738, 0.25087, 1.75144e-13 }, { 0.0782297, 0.226285, -0.367441 }, { 8.89985, 7.66809, 3.59697 }, { 1.96428, 2.23069, 2.23051 }, { 5.9674e-08, 1.07673e-07, 0.0280623 }, { 31.0512, 27.6872, 7.41492 }, { -0.0896936, -0.151411, 0.154604 }, { 0.160631, 0.0657026, 0.0710504 } },
	{ "yellow-matte-plastic", "yellow-matte-plastic", { 0.276745, 0.108557, 0.0203686 }, { 0.806628, 1.99624, 0.977002 }, { 0.0237573, 0.0265756, 0.0305873 }, { 0.558101, 0.537647, 0.710679 }, { 0.0624565, 0.0605404, 0.0541077 }, { 0.049915, 0.0567466, 0.0473569 }, { 26.7398, 24.8437, 17.6137 }, { 2.1481, 2.32064, 1.65758 }, { 5.99553e-08, 5.96103e-08, 6.38097e-08 }, { 41.7723, 39.7008, 40.1979 }, { 0.11569, 0.0891281, 0.0886802 }, { 0.0469972, 0.10918, 0.0466045 } },
	{ "yellow-paint", "paint-yellow", { 0.288876, 0.195348, 0.0314583 }, { 0.449392, 0.412812, 0.168707 }, { 0.650734, 0.190849, 0.16131 }, { 1.21986, 0.0333524, 0.137577 }, { 0.00201415, 7.31639e-13, 0.13481 }, { -0.0897316, -0.0639292, 0.103378 }, { 3.83308, 6.31195, 7.06429 }, { 2.18358, 2.72018, 3.26561 }, { 0.0298708, 0.0764778, 0.0198928 }, { 7.1694, 8.39507, 10.5682 }, { 0.146537, 0.387991, 0.292185 }, { 0.00152911, 0.00123385, 0.000424783 } },
	{ "yellow-phenolic", "phenolic-yellow", { 0.26924, 0.190177, 0.0858303 }, { 0.0861694, 0.0960246, 0.122709 }, { 0.00126171, 0.00197611, 0.00187166 }, { 0.444222, 0.512534, 0.504223 }, { 0.0988294, 0.120249, 0.0929719 }, { 0.0295258, 0.0539044, 0.0471952 }, { 509.34, 294.418, 314.704 }, { 4.59604e-08, 1.76251, 1.81658 }, { 15.6472, 6.5671e-08, 9.96466e-08 }, { 8.52565, 116.04, 118.064 }, { 0.561753, 0.425536, 0.432459 }, { 0.0238647, 0.0190162, 0.0173945 } },
	{ "yellow-plastic", "yellow-plastic", { 0.221083, 0.193042, 0.0403393 }, { 0.265199, 0.340361, 0.0670333 }, { 0.280789, 0.146396, 0.0248514 }, { 0.920018, 0.0356237, 7.80502e-12 }, { 0.144471, 0.10845, 0.125678 }, { 0.113592, 0.0663735, 0.103644 }, { 3.85808, 7.84753, 41.2517 }, { 4.34432, 3.16869, 5.71439 }, { 1.04213e-07, 0.055929, 7.84273 }, { 18.6813, 9.73614, 16.4495 }, { -0.60265, 0.393386, 0.781122 }, { 0.00125483, 0.00140731, 0.00101943 } }
};
	const fresnel::impl<Float> *m_fresnel;
	const data *m_data;
public:
   MI_IMPORT_CORE_TYPES()
	sgd(const char *name);
	~sgd() {delete m_fresnel;}
	Vector3f eval(const Vector3f& i, const Vector3f& o,
	          const void *user_param = NULL) const{
	if (i.z > 0.0 && o.z > 0.0) {
		Vector3f h = normalize(i + o);
		const Vector3f Ks = Vector3f::from_raw(m_data->rhoS);
		const Vector3f Kd = Vector3f::from_raw(m_data->rhoD);
		Vector3f F = fresnel(sat(dot(i, h)));
		Vector3f G = gaf(h, i, o);
		Vector3f D = ndf(h);

		return ((Kd + Ks * (F * D * G) / (i.z * o.z)) / dr::Pi<Float>);
	}
	return Vector3f(0);
}
			  
	Vector3f ndf(const Vector3f& h) const{
	double cos_theta_h = h.z;
	double ndf[3];

	for (uint32_t i = 0; i < 3; ++i)
		ndf[i] = sgd__ndf(cos_theta_h, m_data->alpha[i],
	                      m_data->p[i], m_data->kap[i]);

	return Vector3f::from_raw(ndf);
}

	Vector3f gaf(const Vector3f& h, const Vector3f& i, const Vector3f& o) const{
	return g1(i) * g1(o);
}
	Vector3f g1(const Vector3f& k) const{
	double g1[3];

	for (uint32_t i = 0; i < 3; ++i)
		g1[i] = sgd__g1(k, m_data->theta0[i], m_data->c[i],
	                    m_data->k[i], m_data->lambda[i]);

	return Vector3f::from_raw(g1);
}
	Vector3f fresnel(Float cos_theta_d) const {return m_fresnel->eval(cos_theta_d);}
	const fresnel::impl<Float> &get_fresnel() const {return *m_fresnel;}
};

/* ABC Distribution BRDF */
template <typename Float>
class abc : public brdf<Float> {
	struct data {
		const char *name;
		double kD[3];
		double A[3];
		double B;
		double C;
		double ior;
	};
	static const data s_data[100] = {
	{ "alum-bronze", {0.024514, 0.022346, 0.016220}, {51.714554, 37.932693, 27.373095}, 10482.133789, 0.816737, 2.236525 },
	{ "alumina-oxide", {0.310854, 0.288744, 0.253423}, {2191.683838, 1851.525513, 2064.303467}, 97381.554688, 1.393850, 1.437706 },
	{ "aluminium", {0.000000, 0.002017, 0.005639}, {1137.956787, 1090.751831, 1113.879150}, 218716.593750, 1.092006, 6.282688 },
	{ "aventurnine", {0.052933, 0.060724, 0.052465}, {1835.591309, 1594.715820, 1585.083984}, 324913.375000, 1.110341, 1.403899 },
	{ "beige-fabric", {0.202460, 0.132107, 0.110584}, {0.146439, 0.133698, 0.120796}, 3442009.500000, 0.083210, 10.406549 },
	{ "black-fabric", {0.014945, 0.008790, 0.008178}, {3.109159, 2.772291, 2.723374}, 871444.875000, 0.292300, 1.057058 },
	{ "black-obsidian", {0.000000, 0.000000, 0.000000}, {1447.612183, 1281.373169, 1313.115601}, 223549.906250, 1.196914, 1.485313 },
	{ "black-oxidized-steel", {0.012978, 0.011743, 0.011219}, {1.100971, 1.057305, 1.015528}, 30.661055, 1.869546, 1.660413 },
	{ "black-phenolic", {0.001445, 0.001674, 0.001885}, {104.792976, 99.654602, 100.537247}, 8696.112305, 1.524949, 1.839669 },
	{ "black-soft-plastic", {0.000666, 0.000594, 0.000723}, {0.449065, 0.442911, 0.443035}, 124.908409, 0.618811, 1.902569 },
	{ "blue-acrylic", {0.011106, 0.035057, 0.104254}, {1191.397339, 1099.257202, 1082.125366}, 255868.953125, 1.084736, 1.475061 },
	{ "blue-fabric", {0.008485, 0.009937, 0.034319}, {0.262721, 0.290408, 0.465468}, 3115416.500000, 0.138407, 2.581633 },
	{ "blue-metallic-paint", {0.003770, 0.001030, 0.009177}, {0.495104, 0.469843, 1.097868}, 48.572060, 1.702520, 9.024582 },
	{ "blue-metallic-paint2", {0.000000, 0.000000, 0.000000}, {399.504974, 589.572571, 1115.380737}, 68107.390625, 1.124658, 2.756175 },
	{ "blue-rubber", {0.035403, 0.074082, 0.154563}, {1.295714, 1.279542, 1.248921}, 37.265717, 1.565906, 1.448611 },
	{ "brass", {0.007722, 0.012504, 0.010908}, {2975.156982, 1931.137207, 987.877380}, 244596.531250, 1.136163, 2.698691 },
	{ "cherry-235", {0.037077, 0.014564, 0.007235}, {2.240743, 2.211937, 2.200735}, 123.542030, 1.498179, 1.741225 },
	{ "chrome-steel", {0.006340, 0.008416, 0.009893}, {795.266418, 629.395203, 620.750549}, 122360.437500, 1.486198, 99.999619 },
	{ "chrome", {0.002692, 0.002520, 0.003088}, {986.846741, 771.035522, 696.134583}, 64661.257812, 1.905420, 99.999413 },
	{ "colonial-maple-223", {0.081471, 0.023697, 0.008944}, {1.003711, 0.932037, 0.910946}, 24.846296, 2.121997, 1.805962 },
	{ "color-changing-paint1", {0.002224, 0.002598, 0.005303}, {32.486813, 34.621838, 36.775242}, 2456.161865, 1.388325, 2.066753 },
	{ "color-changing-paint2", {0.000368, 0.000022, 0.001103}, {16.119473, 13.414398, 13.766625}, 1099.414307, 1.331730, 3.186646 },
	{ "color-changing-paint3", {0.007087, 0.004775, 0.003475}, {35.495728, 34.895355, 33.480820}, 2193.942139, 1.442522, 1.918318 },
	{ "dark-blue-paint", {0.002488, 0.007395, 0.034947}, {0.541177, 0.465818, 0.482633}, 24.620041, 1.411947, 1.882264 },
	{ "dark-red-paint", {0.251814, 0.030605, 0.009098}, {0.457430, 0.406366, 0.394383}, 34.332043, 1.023539, 1.670872 },
	{ "dark-specular-fabric", {0.020414, 0.009048, 0.007769}, {1.659937, 1.561381, 1.496425}, 130.493088, 0.975626, 1.356871 },
	{ "delrin", {0.296617, 0.284276, 0.247152}, {4.629758, 4.546722, 4.481231}, 526.069458, 0.931238, 1.545662 },
	{ "fruitwood-241", {0.048981, 0.034689, 0.019745}, {5.976665, 5.760932, 5.568172}, 192.199890, 1.573300, 1.483582 },
	{ "gold-metallic-paint", {0.009164, 0.005649, 0.000424}, {1.685630, 1.147257, 0.418757}, 53.738525, 1.511730, 9.031934 },
	{ "gold-metallic-paint2", {0.000000, 0.000000, 0.000000}, {204.653580, 176.726898, 145.123520}, 1705396.875000, 0.621691, 8.599911 },
	{ "gold-metallic-paint3", {0.000000, 0.000000, 0.000000}, {712.306885, 534.449829, 299.714020}, 94549.195312, 0.970436, 3.965659 },
	{ "gold-paint", {0.150463, 0.080224, 0.022239}, {1.010258, 0.811582, 0.504538}, 46.942558, 1.756670, 7.936110 },
	{ "gray-plastic", {0.099830, 0.100258, 0.092187}, {208.924332, 204.230804, 203.086151}, 33710.691406, 1.023289, 1.400441 },
	{ "grease-covered-steel", {0.000000, 0.000000, 0.000000}, {292.647400, 277.053131, 268.833618}, 15585.240234, 1.170786, 2.556668 },
	{ "green-acrylic", {0.016581, 0.074299, 0.040329}, {1279.877197, 1196.967407, 1221.283569}, 59359.058594, 1.407153, 1.354611 },
	{ "green-fabric", {0.040464, 0.035315, 0.041862}, {0.655209, 0.935078, 1.071461}, 1602934.125000, 0.158285, 1.639683 },
	{ "green-latex", {0.071622, 0.106434, 0.047407}, {0.618660, 0.699535, 0.427560}, 701732.250000, 0.230061, 3.772501 },
	{ "green-metallic-paint", {0.000000, 0.027206, 0.035274}, {0.983370, 1.731232, 1.827718}, 43.999603, 2.172307, 3.235001 },
	{ "green-metallic-paint2", {0.000000, 0.000000, 0.000000}, {289.186859, 578.364441, 361.170837}, 117416.273438, 0.980195, 2.199395 },
	{ "green-plastic", {0.012262, 0.079555, 0.087255}, {2172.868896, 1950.600586, 1812.148193}, 365361.062500, 1.158465, 1.428372 },
	{ "hematite", {0.000000, 0.000273, 0.000528}, {1261.585449, 1214.305420, 1255.108643}, 453901.000000, 1.011091, 2.480358 },
	{ "ipswich-pine-221", {0.049757, 0.016449, 0.004662}, {1.680928, 1.649491, 1.647754}, 49.356201, 1.853272, 1.729069 },
	{ "light-brown-fabric", {0.048293, 0.019506, 0.013650}, {1.725058, 1.219095, 1.080752}, 3929432.500000, 0.155265, 1.153189 },
	{ "light-red-paint", {0.410883, 0.038975, 0.006158}, {0.591418, 0.516404, 0.539583}, 12.604274, 2.137619, 1.743125 },
	{ "maroon-plastic", {0.192636, 0.033163, 0.030345}, {2021.205688, 1751.832031, 1723.350952}, 255986.625000, 1.178735, 1.438042 },
	{ "natural-209", {0.080513, 0.023515, 0.003459}, {2.211772, 2.104441, 1.964800}, 91.766762, 1.154241, 1.641148 },
	{ "neoprene-rubber", {0.248773, 0.211228, 0.170494}, {1.327917, 1.299042, 1.233595}, 63.131874, 1.701130, 1.942141 },
	{ "nickel", {0.006227, 0.006663, 0.007587}, {36.614742, 32.745403, 28.906111}, 705.733887, 1.945258, 5.441883 },
	{ "nylon", {0.232688, 0.232487, 0.207358}, {19.215639, 18.194090, 16.938326}, 3334.451904, 0.886752, 1.382529 },
	{ "orange-paint", {0.407889, 0.131526, 0.003770}, {0.451945, 0.373510, 0.347662}, 11.196602, 1.972428, 1.887030 },
	{ "pearl-paint", {0.181452, 0.162150, 0.136895}, {0.843725, 0.714731, 0.493903}, 30.882542, 1.941556, 9.303603 },
	{ "pickled-oak-260", {0.127652, 0.093767, 0.077318}, {0.813234, 0.803229, 0.781871}, 67.082710, 1.625856, 2.933847 },
	{ "pink-fabric", {0.252317, 0.191945, 0.194732}, {0.434443, 0.450356, 0.464137}, 239167.375000, 0.146117, 2.710645 },
	{ "pink-fabric2", {0.210380, 0.046071, 0.034554}, {0.970204, 0.643877, 0.617551}, 1448228.625000, 0.193577, 4.907917 },
	{ "pink-felt", {0.267287, 0.178590, 0.156370}, {0.225846, 0.226874, 0.221785}, 54.853024, 0.344321, 2.938402 },
	{ "pink-jasper", {0.196444, 0.119217, 0.093004}, {372.695374, 357.903107, 363.584320}, 149508.703125, 0.934382, 1.652987 },
	{ "pink-plastic", {0.360193, 0.079077, 0.053058}, {0.826036, 0.561110, 0.579990}, 16122.769531, 0.217444, 1.810468 },
	{ "polyethylene", {0.183218, 0.202033, 0.207279}, {0.839841, 0.740593, 0.637590}, 3082111.000000, 0.118759, 3.010911 },
	{ "polyurethane-foam", {0.057868, 0.017982, 0.010879}, {0.146946, 0.101932, 0.092976}, 5723581.500000, 0.073071, 3.375618 },
	{ "pure-rubber", {0.287230, 0.252144, 0.211784}, {0.741122, 0.716288, 0.661539}, 85.876984, 0.982928, 1.670395 },
	{ "purple-paint", {0.280217, 0.031190, 0.032274}, {14.006363, 13.769886, 13.753972}, 408.684143, 1.776584, 1.483770 },
	{ "pvc", {0.030285, 0.033992, 0.038152}, {40.608692, 39.046635, 37.916531}, 2901.928223, 1.203405, 1.407985 },
	{ "red-fabric", {0.185341, 0.020015, 0.003522}, {0.934393, 0.314632, 0.273278}, 130451.187500, 0.287254, 3.335066 },
	{ "red-fabric2", {0.116616, 0.012705, 0.002205}, {0.586502, 0.141217, 0.134042}, 8092575.000000, 0.084146, 1.976978 },
	{ "red-metallic-paint", {0.000000, 0.000000, 0.000000}, {2077.145508, 423.511139, 278.621643}, 161943.781250, 1.015078, 2.346517 },
	{ "red-phenolic", {0.162631, 0.021944, 0.006764}, {198.859695, 194.122971, 198.083405}, 57814.757812, 1.049066, 1.752476 },
	{ "red-plastic", {0.268000, 0.044432, 0.010675}, {1.081538, 1.045352, 1.037089}, 34.117554, 1.529797, 1.666450 },
	{ "red-specular-plastic", {0.246119, 0.035688, 0.015152}, {2890.193359, 2449.908936, 2543.227783}, 215481.906250, 1.268530, 1.322196 },
	{ "silicon-nitrade", {0.002805, 0.003669, 0.003936}, {1185.676392, 1031.834106, 1048.810913}, 185127.203125, 1.146923, 1.828749 },
	{ "silver-metallic-paint", {0.000342, 0.000025, 0.002221}, {2.419921, 2.373467, 2.319122}, 58.200935, 1.718832, 9.073233 },
	{ "silver-metallic-paint2", {0.036282, 0.048095, 0.055908}, {1.964917, 1.771562, 1.583858}, 84.018394, 1.141540, 8.961876 },
	{ "silver-paint", {0.159599, 0.127910, 0.108671}, {1.151248, 1.207139, 1.248789}, 38.884895, 1.971779, 8.336617 },
	{ "special-walnut-224", {0.008121, 0.004822, 0.003478}, {0.704096, 0.693490, 0.676845}, 15.779606, 2.668091, 1.883156 },
	{ "specular-black-phenolic", {0.001379, 0.001737, 0.001891}, {2866.697021, 2568.570068, 3042.953613}, 291094.218750, 1.252478, 1.415415 },
	{ "specular-blue-phenolic", {0.002601, 0.011165, 0.029456}, {3080.328613, 2655.796143, 2536.662354}, 298749.312500, 1.201234, 1.437517 },
	{ "specular-green-phenolic", {0.004986, 0.023480, 0.020383}, {3187.704102, 2609.277832, 2708.220947}, 267365.562500, 1.224307, 1.411174 },
	{ "specular-maroon-phenolic", {0.150765, 0.024140, 0.006117}, {2987.999268, 2544.122314, 2451.749023}, 261203.859375, 1.225481, 1.424407 },
	{ "specular-orange-phenolic", {0.325956, 0.051359, 0.005049}, {3129.503174, 2678.228516, 2803.262207}, 462669.000000, 1.185505, 1.439050 },
	{ "specular-red-phenolic", {0.301706, 0.032497, 0.006908}, {2943.967773, 2402.380127, 2540.010010}, 261627.484375, 1.209261, 1.410007 },
	{ "specular-violet-phenolic", {0.066982, 0.015686, 0.018350}, {3117.315918, 2600.714111, 2717.788330}, 257952.890625, 1.251725, 1.407573 },
	{ "specular-white-phenolic", {0.279240, 0.225120, 0.122981}, {4513.249023, 3816.075195, 3782.943848}, 265093.843750, 1.318991, 1.426619 },
	{ "specular-yellow-phenolic", {0.305828, 0.133078, 0.012394}, {3565.679932, 2820.262207, 2784.792236}, 275647.718750, 1.221645, 1.393108 },
	{ "ss440", {0.003789, 0.004251, 0.004992}, {758.115112, 613.785278, 553.360596}, 96436.515625, 1.550240, 99.999336 },
	{ "steel", {0.006151, 0.008683, 0.011281}, {5709.654297, 4641.991211, 4419.819824}, 234651.921875, 1.401774, 2.658189 },
	{ "teflon", {0.321152, 0.313597, 0.305526}, {1.787717, 1.816174, 1.846979}, 41.644257, 2.066234, 1.611215 },
	{ "tungsten-carbide", {0.003085, 0.003349, 0.003638}, {808.095520, 620.232666, 606.583435}, 90806.242188, 1.700902, 99.999184 },
	{ "two-layer-gold", {0.000000, 0.000000, 0.000000}, {103.922806, 91.724258, 76.325027}, 481560.750000, 0.631262, 8.633638 },
	{ "two-layer-silver", {0.000000, 0.000000, 0.000000}, {125.852753, 121.561691, 119.139694}, 1048834.875000, 0.598443, 8.718435 },
	{ "violet-acrylic", {0.044077, 0.006396, 0.023810}, {428.049377, 362.420013, 431.304535}, 1148630.875000, 0.752451, 2.462801 },
	{ "violet-rubber", {0.232500, 0.045230, 0.108831}, {0.840866, 0.796264, 0.822774}, 371.841339, 0.646749, 2.143676 },
	{ "white-acrylic", {0.316720, 0.300911, 0.263206}, {1379.722168, 1291.507080, 1357.878052}, 153332.265625, 1.153521, 1.291205 },
	{ "white-diffuse-bball", {0.319956, 0.272448, 0.178031}, {2.104850, 2.082304, 2.011517}, 75.010391, 1.940200, 1.693726 },
	{ "white-fabric", {0.331064, 0.233945, 0.160661}, {3.545784, 3.587269, 3.573341}, 18457.023438, 0.103981, 1.039433 },
	{ "white-fabric2", {0.074490, 0.068798, 0.080224}, {0.078094, 0.081885, 0.083064}, 6166446.000000, 0.023514, 2.997413 },
	{ "white-marble", {0.222139, 0.208435, 0.183381}, {349.574860, 331.314941, 323.574066}, 49167.527344, 1.153413, 1.650249 },
	{ "white-paint", {0.320307, 0.312274, 0.305093}, {159.141998, 146.118881, 139.397903}, 2777.948730, 1.813009, 1.197418 },
	{ "yellow-matte-plastic", {0.290588, 0.113993, 0.019811}, {85.230881, 80.641701, 78.267593}, 6072.900879, 1.145334, 1.315196 },
	{ "yellow-paint", {0.301059, 0.196365, 0.025507}, {0.446085, 0.446681, 0.399212}, 12.881507, 2.137449, 1.822657 },
	{ "yellow-phenolic", {0.291521, 0.210540, 0.097144}, {498.945862, 477.551605, 480.963745}, 114769.390625, 1.006011, 1.543437 },
	{ "yellow-plastic", {0.237134, 0.202702, 0.030379}, {0.997742, 0.967993, 0.828363}, 289.454254, 0.522704, 1.596116 }
};
	const fresnel::impl<Float> *m_fresnel;
	const data *m_data;
public:
   MI_IMPORT_CORE_TYPES()
	abc(const char *name);
	~abc() {delete m_fresnel;}
	Vector3f eval(const Vector3f& i, const Vector3f& o,
	          const void *user_param = NULL) const{
	if (i.z > 0.0 && o.z > 0.0) {
		Vector3f h = normalize(i + o);
		const Vector3f Kd = Vector3f::from_raw(m_data->kD);
		Vector3f F = fresnel(sat(dot(i, h)));
		double G = gaf(h, i, o);
		Vector3f D = ndf(h);

		return (Kd / dr::Pi<Float>  + (F * D * G) / (dr::Pi<Float> * i.z * o.z));
	}
	return Vector3f(0);
}

	Vector3f ndf(const Vector3f& h) const{
	double cos_theta_h = h.z;
	double ndf[3];

	for (uint32_t i = 0; i < 3; ++i)
		ndf[i] = abc__ndf(cos_theta_h, m_data->A[i], m_data->B, m_data->C);

	return Vector3f::from_raw(ndf);
}

	Float gaf(const Vector3f& h, const Vector3f& i, const Vector3f& o) const;
	Vector3f fresnel(Float cos_theta_d) const {return m_fresnel->eval(cos_theta_d);}
	const fresnel::impl<Float> &get_fresnel() const {return *m_fresnel;}
};

 // namespace djb
NAMESPACE_END(mitsuba)

//
//
//// end header file /////////////////////////////////////////////////////
#endif // DJB_INCLUDE_DJ_BRDF_H

#if DJ_BRDF_IMPLEMENTATION
#include <cmath>
#include <cstdarg>
#include <iostream>     // std::ios, std::istream, std::cout
#include <fstream>      // std::filebuf
#include <cstring>      // memcpy
#include <stdint.h>     // uint32_t
#include <mitsuba/core/mmap.h>
#include <mitsuba/core/mstream.h>
#include <drjit/dynamic.h>
#include <array>
#ifndef DJB_ASSERT
#	include <assert.h>
#	define DJB_ASSERT(x) assert(x)
#endif
/*
#ifndef DJB_LOG
#	include <stdio.h>
#	define DJB_LOG(format, ...) fprintf(stdout, format, ##__VA_ARGS__)
#endif
*/
#ifndef M_PI
#	define M_PI 3.1415926535897932384626433832795
#endif

#ifdef _MSC_VER
#	pragma warning(disable: 4244) // possible loss of data
#endif


NAMESPACE_BEGIN(mitsuba)
// *************************************************************************************************
// utility API
/*
exc::exc(const char *fmt, ...)
{
	char buf[256];
	va_list args;
	va_start(args, fmt);
	vsnprintf(buf, 256, fmt, args);
	va_end(args);

	m_str = std::string(buf);
}*/
/*
template <typename Float>
Vector<Float, 3>::Vector<Float, 3>(Float theta, Float phi)
{
	Float s = sin(theta);
	x = s * cos(phi);
	y = s * sin(phi);
	z = cos(theta);
}


template <typename Float>
static Float inversesqrt(Float x)
{
	DJB_ASSERT(x > 0.0);
	return (1.0 / sqrt(x));
}


template<typename T>
static T max3(const T& x, const T& y, const T& z)
{
	T m = x;

	(m < y) && (m = y);
	(m < z) && (m = z);

	return m;
}
*/
template <typename Float>
static Float inversesqrt(Float x)
{
	DJB_ASSERT(x > 0.0);
	return (1.0 / dr::sqrt(x));
}

// *************************************************************************************************
// Special functions

// See: http://www.johndcook.com/blog/cpp_erf/, by John D. Cook
/*
static Float erf(Float x)
{
	// constants
	Float a1 =  (Float)0.254829592;
	Float a2 = -(Float)0.284496736;
	Float a3 =  (Float)1.421413741;
	Float a4 = -(Float)1.453152027;
	Float a5 =  (Float)1.061405429;
	Float p  =  (Float)0.3275911;

	// Save the sign of x
	int sign = 1;
	if (x < 0)
		sign = -1;
	x = fabs(x);

	// A&S formula 7.1.26
	Float t = 1.0/(1.0 + p*x);
	Float y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

	return sign * y;
}
*/
// See: Approximating the erfinv function, by Mike Giles
template <typename Float>
static Float erfinv(Float u)
{
	Float w, p;

	w = -logf(((Float)1.0 - u) * ((Float)1.0 + u));
	if (w < (Float)5.0) {
		w = w - (Float)2.500000;
		p = (Float)2.81022636e-08;
		p = (Float)3.43273939e-07 + p * w;
		p = (Float)-3.5233877e-06 + p * w;
		p = (Float)-4.39150654e-06 + p * w;
		p = (Float)0.00021858087 + p * w;
		p = (Float)-0.00125372503 + p * w;
		p = (Float)-0.00417768164 + p * w;
		p = (Float)0.246640727 + p * w;
		p = (Float)1.50140941 + p * w;
	} else {
		w = sqrt(w) - (Float)3.0;
		p = (Float)-0.000200214257;
		p = (Float)0.000100950558 + p * w;
		p = (Float)0.00134934322 + p * w;
		p = (Float)-0.00367342844 + p * w;
		p = (Float)0.00573950773 + p * w;
		p = (Float)-0.0076224613 + p * w;
		p = (Float)0.00943887047 + p * w;
		p = (Float)1.00167406 + p * w;
		p = (Float)2.83297682 + p * w;
	}

	return p * u;
}

// *************************************************************************************************
// Sample warps
/*
template <typename Float>
static void
uniform_to_concentric(Float u1, Float u2, Float *r1, Float *r2)
{   using Mask   = dr::mask_t<Float>;
	/* Concentric map code with less branching (by Dave Cline), see
	   http://psgraphics.blogspot.ch/2011/01/improved-code-for-concentric-map.html */
//	Float r1 = 2.0 * u1 - 1.0;
//	Float r2 = 2.0 * u2 - 1.0;
	/*
	Float phi, r;

	if (r1 == 0 && r2 == 0) {
		r = phi = 0;
	} else if (dr::abs(r1 * r1) > dr::abs(r2 * r2)) {
		r = r1;
		phi = (M_PI / 4.0) * (r2 / r1);
	} else {
		r = r2;
		phi = (M_PI / 2.0) - (r1 / r2) * (M_PI / 4.0);
	}

	(*x) = r * dr::cos(phi);
	(*y) = r * dr::sin(phi);*//*
	Float x = dr::fmsub(2.f, u1, 1.f),
          y = dr::fmsub(2.f, u2, 1.f);
    Mask is_zero         = dr::eq(x, 0.f) &&
                           dr::eq(y, 0.f),
         quadrant_1_or_3 = dr::abs(x) < dr::abs(y);

    Float r  = dr::select(quadrant_1_or_3, y, x),
          rp = dr::select(quadrant_1_or_3, x, y);

    Float phi = 0.25f * dr::Pi<Float> * rp / r;
    dr::masked(phi, quadrant_1_or_3) = 0.5f * dr::Pi<Float> - phi;
    dr::masked(phi, is_zero) = 0.f;

    auto [s, c] = dr::sincos(phi);
	(*r1) = r * c;
	(*r2) = r * s;


}
*/
// *************************************************************************************************
// Geometric operations

//---------------------------------------------------------------------------
// rotate vector along one axis
template <typename Float>
static Vector<Float, 3> rotate_vector(const Vector<Float, 3>& x, const Vector<Float, 3>& axis, Float angle)
{
	Float cos_angle = dr::cos(angle);
	Float sin_angle = dr::sin(angle);
	Vector<Float, 3> out = cos_angle * x;
	Float tmp1 = dr::dot(axis, x);
	Float tmp2 = tmp1 * (1.0 - cos_angle);
	out+= axis * tmp2;
	out+= sin_angle * dr::cross(axis, x);

	return out;
}


// *************************************************************************************************
// BRDF API
template <typename Float>
void brdf<Float>::io_to_hd(const Vector3f& i, const Vector3f& o, Vector3f *h, Vector3f *d)
{
	const Vector3f y_axis = Vector3f(0, 1, 0);
	const Vector3f z_axis = Vector3f(0, 0, 1);
	Float theta_h, phi_h;

	(*h) = dr::normalize(i + o);
	auto h_p=*h;
	//xyz_to_theta_phi(*h, &theta_h, &phi_h);
	theta_h=dr::acos(h_p.z());	
	phi_h=dr::atan2(h_p.y(), h_p.x());
	Vector3f tmp = rotate_vector(i, z_axis, -phi_h);
	(*d) = dr::normalize(rotate_vector(tmp, y_axis, -theta_h));
}
template <typename Float>
void brdf<Float>::hd_to_io(const Vector3f& h, const Vector3f& d, Vector3f *i, Vector3f *o)
{
	const Vector3f y_axis = Vector3f(0, 1, 0);
	const Vector3f z_axis = Vector3f(0, 0, 1);
	Float theta_h, phi_h;

//	xyz_to_theta_phi(h, &theta_h, &phi_h);
	theta_h=dr::acos(h.z());	
	phi_h=dr::atan2(h.y(), h.x());
	Vector3f tmp = rotate_vector(d, y_axis, theta_h);
	(*i) = dr::normalize(rotate_vector(tmp, z_axis, phi_h));
	(*o) = dr::normalize(2.0 * dr::dot((*i), h) * h - (*i));
}
//template <typename Float>
/*
Vector3f
brdf<Float>::eval_hd(const Vector3f& h, const Vector3f& d, const void *user_param) const
{   
	Vector3f i, o;
	hd_to_io(h, d, &i, &o);

	return eval(i, o, user_param);
}
template <typename Float>
Vector3f brdf<Float>::evalp(const Vector3f& i, const Vector3f& o, const void *user_param) const
{
	return eval(i, o, user_param) * i.z();
}
template <typename Float>
Vector3f brdf<Float>::evalp_hd(const Vector3f& h, const Vector3f& d, const void *user_param) const
{
	Vector3f i, o;
	hd_to_io(h, d, &i, &o);

	return eval(i, o, user_param) * i.z();
}
template <typename Float>
Vector3f
brdf<Float>::evalp_is(
	Float u1, Float u2,
	const Vector3f& o, Vector3f *i, Float *pdf,
	const void *user_param
) const {
	const Vector3f i_ = sample(u1, u2, o, user_param);
	Float pdf_ = this->pdf(i_, o);
	if (i) (*i) = i_;
	if (pdf) (*pdf) = pdf_;

	return (evalp(i_, o, user_param) / pdf_);
}
template <typename Float>
Vector<Float, 3>
brdf<Float>::sample(Float u1, Float u2, const Vector<Float, 3>& o, const void *user_param) const
{
	Float x, y, z;

	uniform_to_concentric(u1, u2, &x, &y);
	z = dr::sqrt(1.0 - x * x - y * y);
	DJB_ASSERT(z >= 0.0);

	return Vector<Float, 3>(x, y, z);
}*/

//template <typename Float>
//Float brdf<Float>::pdf(const Vector<Float, 3>& i, const Vector<Float, 3>& o, const void *user_param) const
//{
//	return /* cos(theta_i) */i.z() / dr::Pi<Float>;
//}

// *************************************************************************************************
// Lambertian API implementation

//---------------------------------------------------------------------------
// Lambertian Params

/* Constructors */
/*
template <typename Float>
lambert<Float>::params::params(const Vector3f & reflectance):
	m_reflectance(reflectance)
{}

//---------------------------------------------------------------------------
// Lambertian Evaluation
template <typename Float>
Vector<Float, 3>
lambert<Float>::eval(const Vector<Float, 3>& i, const Vector<Float, 3>& o, const void *user_param) const
{
	const lambert::params params = 
		user_param ? *reinterpret_cast<const lambert::params *>(user_param)
		           : lambert::params();

	return (params.m_reflectance / dr::Pi<Float>);
}*/

// *************************************************************************************************
// MERL API implementation

// Copyright 2005 Mitsubishi Electric Research Laboratories All Rights Reserved.

// Permission to use, copy and modify this software and its documentation without
// fee for educational, research and non-profit purposes, is hereby granted, provided
// that the above copyright notice and the following three paragraphs appear in all copies.

// To request permission to incorporate this software into commercial products contact:
// Vice President of Marketing and Business Development;
// Mitsubishi Electric Research Laboratories (MERL), 201 Broadway, Cambridge, MA 02139 or 
// <license@merl.com>.

// IN NO EVENT SHALL MERL BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL,
// OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND
// ITS DOCUMENTATION, EVEN IF MERL HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.

// MERL SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE PROVIDED
// HEREUNDER IS ON AN "AS IS" BASIS, AND MERL HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT,
// UPDATES, ENHANCEMENTS OR MODIFICATIONS.

#define MERL_SAMPLING_RES_THETA_H   90
#define MERL_SAMPLING_RES_THETA_D   90
#define MERL_SAMPLING_RES_PHI_D    360

#define MERL_RED_SCALE   (1.00 / 1500.0)
#define MERL_GREEN_SCALE (1.15 / 1500.0)
#define MERL_BLUE_SCALE  (1.66 / 1500.0)

//---------------------------------------------------------------------------
// Lookup theta_half index
// This is a non-linear mapping!
// In:  [0 .. pi/2]
// Out: [0 .. 89]
template <typename Float>
static Float
theta_half_index(Float theta_half)
{  MI_IMPORT_CORE_TYPES();
	Mask active = true;
    active &= theta_half  > 0;
	if (dr::none_or<false>(active))
		return 0;
	Float theta_half_deg = ((theta_half / (M_PI/2.0)) * MERL_SAMPLING_RES_THETA_H);
	Float temp = theta_half_deg * MERL_SAMPLING_RES_THETA_H;
	temp = dr::sqrt(temp);
	//int ret_val =dr::floor2int<int> temp;
	temp=dr::clamp(temp,0.0,Float(MERL_SAMPLING_RES_THETA_H - 1));

	/*if (ret_val < 0)
		ret_val = 0;
	else if (ret_val >= MERL_SAMPLING_RES_THETA_H)
		ret_val = MERL_SAMPLING_RES_THETA_H - 1;*/
	return temp;
}

//---------------------------------------------------------------------------
// Lookup theta_diff index
// In:  [0 .. pi/2]
// Out: [0 .. 89]
template <typename Float>
static Float 
theta_diff_index(Float theta_diff)
{
	Float tmp = theta_diff / (M_PI* 0.5) * MERL_SAMPLING_RES_THETA_D;
	//int tmp1=dr::floor2int<int>(tmp);
	tmp=dr::clamp(tmp,0.0,Float(MERL_SAMPLING_RES_THETA_D - 1));
	
	/*if (tmp1 < 0)
		return 0;
	else if (tmp1 < MERL_SAMPLING_RES_THETA_D - 1)
		return tmp1;
	else
		return MERL_SAMPLING_RES_THETA_D - 1;*/
	return tmp;
}

//---------------------------------------------------------------------------
// Lookup phi_diff index
template <typename Float>
static Float
phi_diff_index(Float phi_diff)
{
	// Because of reciprocity, the BRDF is unchanged under
	// phi_diff -> phi_diff + dr::Pi<Float>
	//if (phi_diff < 0.0)
	//	phi_diff += M_PI;
    phi_diff=dr::select(phi_diff < 0.0, phi_diff + dr::Pi<Float>, phi_diff);
	// In: phi_diff in [0 .. pi]
	// Out: tmp in [0 .. 179]
	Float tmp = phi_diff / M_PI * MERL_SAMPLING_RES_PHI_D / 2;
	//uint32_t tmp1=dr::floor2int<uint32_t>(tmp);

    tmp=dr::clamp(tmp,0.0,Float(MERL_SAMPLING_RES_PHI_D / 2 - 1));
	/*if (tmp1 < 0)
		return 0;
	else if (tmp1 < MERL_SAMPLING_RES_PHI_D / 2 - 1)
		return tmp1;
	else
		return MERL_SAMPLING_RES_PHI_D / 2 - 1;*/
	return tmp;
}
//NAMESPACE_BEGIN(mitsuba)
// XXX End of 
// Copyright 2005 Mitsubishi Electric Research Laboratories All Rights Reserved.
template <typename Float>
class merl_brdf: public brdf<Float>{
	MI_IMPORT_CORE_TYPES()
	 using FloatStorage                = DynamicBuffer<Float>;
	std::vector<double> m_samples;
	FloatStorage m_data;
public:
   
	uint32_t n, dims[3];
    merl_brdf() = default;
	merl_brdf(const fs::path  &filename){
	
	ref<FileStream> file_merl = new FileStream(filename, FileStream::ERead);
	uint32_t file_size = file_merl->size();
//	if(file_size != sizeof(Float)*180*90*90*3+4*3)
//	    Throw("The file \"%f\" is to be a valid MERL\"%f\" file!", file_size,sizeof(Float)*180*90*90*3);
	if (!file_merl->can_read())
		Throw("Could not open the file \"%s\"!", filename.filename().string().c_str());	
	//uint32_t n, dims[3];	
	file_merl->read((char *)dims, 4 * 3);	
	n = dims[0] * dims[1] * dims[2];
	if (n != 180*90*90)
	   Throw ("djb_error:\"%f\" Failed to read MERL header\n",n);
	// allocate brdf and read data
	m_samples.resize(3 * n);
//	std::unique_ptr<double[]> tmp(new double[sizeof(double) * 3 * n]);
//	file_merl->read(tmp.get(), sizeof(double) * 3 * n);  
	file_merl->read((char *)&m_samples[0], sizeof(double) * 3 * n);

	if (file_merl->is_closed())
		Throw ("djb_error: Reading failed \"%s\"\n", filename.filename().string().c_str());
	//Throw ("djb_error:\"%f\" Failed to read MERL header\n",m_samples[100]);
	/*
	std::fstream f(filename.filename().string().c_str(), std::fstream::in | std::fstream::binary);
	uint32_t n, dims[3];
     
	// check file
	if (!f.is_open())
		Throw ("djb_error: Failed to open  \"%s\"\n", filename.filename().string().c_str());

	// read header
	f.read((char *)dims, 4 * 3);/*bytes*/
	/*n = dims[0] * dims[1] * dims[2];

	if (n != 180*90*90*3)
	   Throw ("djb_error:\"%f\" Failed to read MERL header\n",n);

	// allocate brdf and read data
	m_samples.resize(3 * n);
	f.read((char *)&m_samples[0], sizeof(double) * 3 * n);
	if (f.fail())
		Throw ("djb_error: Reading failed \"%s\"\n", filename.filename().string().c_str());*/
	//////////mitsuba::MemoryMappedFile file = mitsuba::////////////////////////////////////////////////	
   /* ref<MemoryMappedFile>  mm_merl=new MemoryMappedFile(filename, false);
	if (!mm_merl)
		Throw("Could not memory-map the file \"%s\"!", filename.filename().string().c_str());
	if (mm_merl->size()< sizeof(Float)*180*90*90*3)
		{//Throw("The file \"%s\" is too small to be a valid MERL file!", filename.filename().string().c_str());	
		Throw("The file \"%f\" is too small to be a valid MERL\"%f\" file!", mm_merl->size(),sizeof(Float)*180*90*90*3);	}
	 ref<MemoryStream> stream = new MemoryStream(mm_merl->data(), mm_merl->size());
	 stream->read((char *)dims, 12);
	 n = dims[0] * dims[1] * dims[2];
	if (n != 180*90*90)
	   Throw ("djb_error:\"%f\" Failed to read MERL header\n",n);
	const char * ptr_merl            = ( const char *) mm_merl->data()+12;   */
	std::unique_ptr<ScalarFloat[]> data_out(new ScalarFloat[n * 3]);
	ScalarFloat *data_out_ptr = data_out.get();
	double *m_samples_ptr = m_samples.data();
	 for (uint32_t k = 0; k < n * 3; ++k)
        *data_out_ptr++ = *m_samples_ptr++ ;
	m_data = dr::load<FloatStorage>(data_out.get(), (n * 3));
    
//	m_data = dr::load<FloatStorage>(m_samples.data(),sal);// have a problem , need to load one by one
	/*UInt32 idx_r=1; UInt32 idx_g=2; UInt32 idx_b=3;
	Vector3f rgb;
	rgb.x()= dr::gather<Float>(m_data, idx_r);
	rgb.y()= dr::gather<Float>(m_data, idx_g);
	rgb.z()= dr::gather<Float>(m_data, idx_b);
	*/
	
}
Vector3f eval(const Vector3f& i, const Vector3f& o,
	          const void *user_param = NULL) const
			  {
	// convert to half / diff angle coordinates
	Vector3f h, d;
	Float theta_h, phi_h, theta_d, phi_d;
	brdf<Float>::io_to_hd(i, o, &h, &d);
	theta_h=dr::acos(h.z());	
	phi_h=dr::atan2(h.y(), h.x());
	theta_d=dr::acos(d.z());	
	phi_d=dr::atan2(d.y(), d.x());
//	xyz_to_theta_phi<Float>(h, &theta_h, &phi_h);
//	xyz_to_theta_phi<Float>(d, &theta_d, &phi_d);
    
	// compute indexes
	UInt32 idx_r = phi_diff_index(phi_d)
	          + theta_diff_index(theta_d)
	          * MERL_SAMPLING_RES_PHI_D / 2
	          + theta_half_index(theta_h)
	          * MERL_SAMPLING_RES_PHI_D / 2
	          * MERL_SAMPLING_RES_THETA_D;
	UInt32 idx_g = idx_r + MERL_SAMPLING_RES_THETA_H
	          * MERL_SAMPLING_RES_THETA_D 
	          * MERL_SAMPLING_RES_PHI_D / 2;
	UInt32 idx_b = idx_r + MERL_SAMPLING_RES_THETA_H
	          * MERL_SAMPLING_RES_THETA_D 
	          * MERL_SAMPLING_RES_PHI_D;
  //  Throw ("djb_error:\"%f\" Failed to read MERL header\n",idx_g);
	// get color
//	idx_r=10;  idx_g=1458000;  idx_b=4370000;
//	uint32_t idx_r1=0; uint32_t idx_g1=1458000; uint32_t idx_b1=4273000;
	Vector3f rgb;
    
    rgb.x()= dr::gather<Float>(m_data, idx_r)* MERL_RED_SCALE;
	rgb.y()= dr::gather<Float>(m_data, idx_g)* MERL_GREEN_SCALE;
	rgb.z()= dr::gather<Float>(m_data, idx_b)* MERL_BLUE_SCALE;
//	rgb.x() = m_samples[4371001] * MERL_RED_SCALE;
//	rgb.y() = m_samples[4371002] * MERL_GREEN_SCALE;
//	rgb.z() = m_samples[4371003] * MERL_BLUE_SCALE;
    Mask active = true;
	active = active && (rgb.x() >= 0.0 && rgb.y() >= 0.0 && rgb.z() >= 0.0);
	if (dr::none_or<false>(active)) {
#ifndef NVERBOSE
    Throw("djb_verbose: below horizon");
	//DJB_LOG("djb_verbose: below horizon\n");
#endif
		return Vector3f(0);
	}

	return rgb;
}
Vector3f evalp(const Vector3f& i, const Vector3f& o,
	          const void *user_param = NULL) const
			  {
				
				return eval(i, o, user_param) * i.z();
				}
const std::vector<double>& get_samples() const {return m_samples;}
};
//NAMESPACE_END(mitsuba)
//---------------------------------------------------------------------------

// MERL Contructor

//merl::merl(const char *filename)
//{
//	std::fstream f(filename, std::fstream::in | std::fstream::binary);
//	int32_t n, dims[3];

	// check file
//	if (!f.is_open())//
//		throw exc("djb_error: Failed to open %s\n", filename);

	// read header
//	f.read((char *)dims, /*bytes*/4 * 3);
//	n = dims[0] * dims[1] * dims[2];
//	if (n <= 0)
//		throw exc("djb_error: Failed to read MERL header\n");

	// allocate brdf and read data
//	m_samples.resize(3 * n);
//	f.read((char *)&m_samples[0], sizeof(double) * 3 * n);
//	if (f.fail())
//		throw exc("djb_error: Reading %s failed\n", filename);
//}

//---------------------------------------------------------------------------
// look up the BRDF.
/*
Vector<Float, 3> merl::eval(const Vector<Float, 3>& i, const Vector<Float, 3>& o, const void *user_param) const
{
	// convert to half / diff angle coordinates
	Vector<Float, 3> h, d;
	Float theta_h, phi_h, theta_d, phi_d;
	io_to_hd(i, o, &h, &d);
	xyz_to_theta_phi(h, &theta_h, &phi_h);
	xyz_to_theta_phi(d, &theta_d, &phi_d);

	// compute indexes
	int idx_r = phi_diff_index(phi_d)
	          + theta_diff_index(theta_d)
	          * MERL_SAMPLING_RES_PHI_D / 2
	          + theta_half_index(theta_h)
	          * MERL_SAMPLING_RES_PHI_D / 2
	          * MERL_SAMPLING_RES_THETA_D;
	int idx_g = idx_r + MERL_SAMPLING_RES_THETA_H
	          * MERL_SAMPLING_RES_THETA_D 
	          * MERL_SAMPLING_RES_PHI_D / 2;
	int idx_b = idx_r + MERL_SAMPLING_RES_THETA_H
	          * MERL_SAMPLING_RES_THETA_D 
	          * MERL_SAMPLING_RES_PHI_D;

	// get color
	Vector<Float, 3> rgb;
	rgb.x = m_samples[idx_r] * MERL_RED_SCALE;
	rgb.y = m_samples[idx_g] * MERL_GREEN_SCALE;
	rgb.z = m_samples[idx_b] * MERL_BLUE_SCALE;

	if (rgb.x < 0.0 || rgb.y < 0.0 || rgb.z < 0.0) {
#ifndef NVERBOSE
	DJB_LOG("djb_verbose: below horizon\n");
#endif
		return Vector<Float, 3>(0);
	}

	return rgb;
}
*/
// *************************************************************************************************
// UTIA API implementation (based on Jiri Filip's implementation)
/*
#define DJB__UTIA_STEP_T  15.0
#define DJB__UTIA_STEP_P   7.5
#define DJB__UTIA_NTI      6
#define DJB__UTIA_NPI     48
#define DJB__UTIA_NTV      6
#define DJB__UTIA_NPV     48
#define DJB__UTIA_PLANES   3
*/
//---------------------------------------------------------------------------
// Read UTIA BRDF data
template <typename Float>
utia<Float>::utia(const char *filename)
{
	// open file
	std::fstream f(filename, std::fstream::in | std::fstream::binary);
	if (!f.is_open())
		throw exc("djb_error: Failed to open %s\n", filename);

	// allocate memory
	uint32_t cnt = DJB__UTIA_PLANES * DJB__UTIA_NTI * DJB__UTIA_NPI
	        * DJB__UTIA_NTV * DJB__UTIA_NPV;
	m_samples.resize(cnt);

	// read data
	f.read((char *)&m_samples[0], sizeof(double) * cnt);

	// normalize
	m_norm = 1.0;
	normalize();
	if (f.fail())
		throw exc("djb_error: Reading %s failed\n", filename);
}

//---------------------------------------------------------------------------
// Look up the UTIA BRDF
/*
template <typename Float>
Vector<Float, 3> utia<Float>::eval(const Vector<Float, 3>& i, const Vector<Float, 3>& o, const void *user_param) const
{
	Float r2d = 180.0 / dr::Pi<Float>;
	Float theta_i = r2d * acos(i.z);
	Float theta_o = r2d * acos(o.z);
	Float phi_i = r2d * atan2(i.y, i.x);
	Float phi_o = r2d * atan2(o.y, o.x);

	// make sure we're above horizon
	if (theta_i >= 90.0 || theta_o >= 90.0)
		return Vector3f(0);

	// make sure phi is in [0, 360)
	while (phi_i < 0.0) {phi_i+= 360.0;}
	while (phi_o < 0.0) {phi_o+= 360.0;}
	while (phi_i >= 360) {phi_i-= 360.0;}
	while (phi_o >= 360) {phi_o-= 360.0;}

	uint32_t iti[2], itv[2], ipi[2], ipv[2];
	iti[0] = (uint32_t)floor(theta_i / DJB__UTIA_STEP_T);
	iti[1] = iti[0] + 1;
	if(iti[0] > DJB__UTIA_NTI - 2) {
		iti[0] = DJB__UTIA_NTI - 2;
		iti[1] = DJB__UTIA_NTI - 1;
	}
	itv[0] = (uint32_t)floor(theta_o / DJB__UTIA_STEP_T);
	itv[1] = itv[0] + 1;
	if(itv[0] > DJB__UTIA_NTV - 2) {
		itv[0] = DJB__UTIA_NTV - 2;
		itv[1] = DJB__UTIA_NTV - 1;
	}

	ipi[0] = (uint32_t)floor(phi_i / DJB__UTIA_STEP_P);
	ipi[1] = ipi[0] + 1;
	ipv[0] = (uint32_t)floor(phi_o / DJB__UTIA_STEP_P);
	ipv[1] = ipv[0] + 1;

	Float sum;
	Float wti[2], wtv[2], wpi[2], wpv[2];
	wti[1] = theta_i - (Float)(DJB__UTIA_STEP_T * iti[0]);
	wti[0] = (Float)(DJB__UTIA_STEP_T * iti[1]) - theta_i;
	sum = wti[0] + wti[1];
	wti[0]/= sum;
	wti[1]/= sum;
	wtv[1] = theta_o - (Float)(DJB__UTIA_STEP_T * itv[0]);
	wtv[0] = (Float)(DJB__UTIA_STEP_T * itv[1]) - theta_o;
	sum = wtv[0] + wtv[1];
	wtv[0]/= sum;
	wtv[1]/= sum;

	wpi[1] = phi_i - (Float)(DJB__UTIA_STEP_P * ipi[0]);
	wpi[0] = (Float)(DJB__UTIA_STEP_P * ipi[1]) - phi_i;
	sum = wpi[0] + wpi[1];
	wpi[0]/= sum;
	wpi[1]/= sum;
	wpv[1] = phi_o - (Float)(DJB__UTIA_STEP_P * ipv[0]);
	wpv[0] = (Float)(DJB__UTIA_STEP_P * ipv[1]) - phi_o;
	sum = wpv[0] + wpv[1];
	wpv[0]/= sum;
	wpv[1]/= sum;

	if(ipi[1] == DJB__UTIA_NPI)
		ipi[1] = 0;
	if(ipv[1] == DJB__UTIA_NPV)
		ipv[1] = 0;

	uint32_t nc = DJB__UTIA_NPV * DJB__UTIA_NTV;
	uint32_t nr = DJB__UTIA_NPI * DJB__UTIA_NTI;
	Float RGB[DJB__UTIA_PLANES];
	for(uint32_t isp = 0; isp < DJB__UTIA_PLANES; ++isp) {
		uint32_t i, j, k, l;

		RGB[isp] = 0.0;
		for(i = 0; i < 2; ++i)
		for(j = 0; j < 2; ++j)
		for(k = 0; k < 2; ++k)
		for(l = 0; l < 2; ++l) {
			Float w = wti[i] * wtv[j] * wpi[k] * wpv[l];
			uint32_t idx = isp * nr * nc + nc * (DJB__UTIA_NPI * iti[i] + ipi[k]) 
			        + DJB__UTIA_NPV * itv[j] + ipv[l];

			RGB[isp]+= w * (Float)m_samples[idx];
		}
		if (RGB[isp] > 0.0375)
			RGB[isp] = pow((Float)(RGB[isp] + 0.055)/1.055, (Float)2.4); 
		else
			RGB[isp]/= (Float)12.92;
		RGB[isp]*= (Float)100.0; 
	}
	return Vector3f(
		max((Float)0, RGB[0]),
		max((Float)0, RGB[1]),
		max((Float)0, RGB[2])
	);
}
*/
//---------------------------------------------------------------------------
// Some UTIA BRDFs contain negative values for some materials. Therefore, we 
// clamp them. The data also needs to be scaled.
template <typename Float>
void utia<Float>::normalize()
{
	// clamp to zero
	for (uint32_t i = 0; i < (uint32_t)m_samples.size(); ++i) {
#ifndef NVERBOSE
		if (m_samples[i] < 0.0)
		    Throw("djb_verbose: negative UTIA BRDF value found, set to 0");
			//DJB_LOG("djb_verbose: negative UTIA BRDF value found, set to 0\n");
#endif
		m_samples[i] = dr::maximum(m_samples[i], 0.0);
	}

	// scale
	Float k = /* Magic constant provided by Jiri Filip */ 1.f / 140.f;
	for (uint32_t i = 0; i < (uint32_t)m_samples.size(); ++i)
		m_samples[i]*= k;
}

// *************************************************************************************************
// Private Spline API
namespace spline {

int32_t uwrap_repeat(int32_t i, int32_t edge)
{
	while (i >= edge) i-= edge;
	while (i < 0    ) i+= edge;

	return i;
}

int32_t uwrap_edge(int32_t i, int32_t edge)
{
	if      (i >= edge) i = edge - 1;
	else if (i < 0    ) i = 0;

	return i;
}

//typedef int32_t (*uwrap_callback)(int32_t, int32_t);
/*
template <typename Float>
Float lerp(const Float& x1, const Float& x2, Float u)
{
	return x1 + u * (x2 - x1);
}
template <typename Float>
Float lerp(const Float& x1, const Float& x2, Float u)
{
	return x1 + u * (x2 - x1);
}*/
/*
template <typename Float>
Float modf( Float x,  uint32_t intpart)
{
	intpart=uint32_t(x);
	return x=x-intpart;
}*/
template <typename Float>
Float eval(std::vector<Float> points, std::string wrap,Float u)
{    MI_IMPORT_CORE_TYPES();
	
	uint32_t edge   = points.size();
	using FloatStorage                = DynamicBuffer<Float>;

	FloatStorage s_data;
    std::unique_ptr<Float[]> data_out(new Float[edge]);
	Float *data_out_ptr = data_out.get();
	Float *points_ptr = points.data();
	 for (uint32_t k = 0; k < edge; ++k)
        *data_out_ptr++ = *points_ptr++ ;
	s_data = dr::load<FloatStorage>(data_out.get(), (edge));
    
	//Float intpart=dr::floor((u * edge)-u);
	UInt32 intpart=UInt32 ((u * edge)-u);
	//int32_t intpart=1;
	Float frac=u * edge -u - dr::floor((u * edge)-u); 
    //int32_t i11=dr::floor2int<int32_t>(intpart);
	
	//Float frac = modf(u * edge - u, &intpart);
	UInt32 i1 =dr::clamp(intpart, 0,edge-1);
	UInt32 i2 = dr::clamp(intpart + 1,0, edge-1);
   if (wrap == "edge") {
	UInt32 i1 = dr::clamp(intpart, 0,edge-1);
	UInt32 i2 = dr::clamp(intpart + 1,0, edge-1);
   }
   else if (wrap == "repeat") {
   UInt32 i1 = dr::clamp(intpart, intpart+edge,intpart-edge);
	UInt32 i2 = dr::clamp(intpart + 1,intpart+edge,intpart-edge);
   }
	const Float &p1 = dr::gather<Float>(s_data, i1);
	const Float &p2 = dr::gather<Float>(s_data, i2);

	return dr::lerp(p1, p2, frac);
}

template <typename Float>
Float
eval2d(
	const std::vector<Float>& points, 
	int32_t w,
	int32_t h,
	Float u1,
	Float u2,
	std::string wrap1,
	std::string wrap2
) {
	// compute weights and indices
	int32_t intpart1= dr::floor2int ((u1 * w)-u1);
	//int32_t intpart1=1;
	Float frac1 = u1 * w -u1- dr::floor((u1 * w)-u1);
//	Float frac1 = modf(u1 * w - u1, &intpart1);
  /*if (wrap1 == "edge") {
	UInt32 i1 = dr::clamp(intpart1, 0,w-1);
	UInt32 i2 = dr::clamp(intpart1 + 1,0, w-1);
   }
   else if (wrap1 == "repeat") {
   UInt32 i1 = dr::clamp(intpart, i+w,i-w);
	UInt32 i2 = dr::clamp(intpart + 1,i+w,i-w);
   }*/
	int32_t i1 = std::min(intpart1    , w);
	int32_t i2 = std::min(intpart1 + 1, w);
	int32_t intpart2=dr::floor2int<int32_t>((u2 * h)-u2);
//	int32_t intpart2=1;
	Float frac2 = u2 * h-u2 - dr::floor((u2 * h)-u2);
//	Float frac2 = modf(u2 * h - u2, &intpart2);
	int32_t j1 = std::min(intpart2    , h);
	int32_t j2 = std::min(intpart2 + 1, h);

	// fetches
	const Float& p1 = points[i1 + w * j1];
	const Float& p2 = points[i2 + w * j1];
	const Float& p3 = points[i1 + w * j2];
	const Float& p4 = points[i2 + w * j2];

	// return bilinear interpolation
	Float tmp1 = dr::lerp(p1, p2, frac1);
	Float tmp2 =dr::lerp(p3, p4, frac1);;
	return dr::lerp(tmp1, tmp2, frac2);
}

} // namespace spline

// *************************************************************************************************
// Fresnel API implementation
namespace fresnell {
template <typename Float>
void ior_to_f0(Float ior, Float *f0)
{
	DJB_ASSERT(ior > 0.0 && "Invalid ior");
	DJB_ASSERT(f0 && "Null output ptr");
	Float tmp = (ior - 1.0) / (ior + 1.0);

	(*f0) = tmp * tmp;
}
template <typename Float>
void ior_to_f0(const Vector<Float, 3>& ior, Vector<Float, 3> *f0)
{
	DJB_ASSERT(f0 && "Null output ptr");
	ior_to_f0(ior.x, &f0->x);
	ior_to_f0(ior.y, &f0->y);
	ior_to_f0(ior.z, &f0->z);
}
template <typename Float>
void f0_to_ior(Float f0, Float *ior)
{
	DJB_ASSERT(ior && "Null output ptr");
	if (f0 == 1.0) {
		(*ior) = 1.0;
	} else {
		Float sqrt_f0 = sqrt(f0);

		(*ior) = (1.0 + sqrt_f0) / (1.0 - sqrt_f0);
	}
}
template <typename Float>
void f0_to_ior(const Vector<Float, 3>& f0, Vector<Float, 3> *ior)
{
	DJB_ASSERT(ior && "Null output ptr");
	f0_to_ior(f0.x, &ior->x);
	f0_to_ior(f0.y, &ior->y);
	f0_to_ior(f0.z, &ior->z);
}
template <typename Float>
static Float unpolarized__eval(Float cos_theta_d, Float ior)
{
	Float c = cos_theta_d;
	Float n = ior;
	Float g = sqrt(n * n + c * c - 1.0);
	Float tmp1 = c * (g + c) - 1.0;
	Float tmp2 = c * (g - c) + 1.0;
	Float tmp3 = (tmp1 * tmp1) / (tmp2 * tmp2);
	Float tmp4 = ((g - c) * (g - c)) / ((g + c) * (g + c));

	return ((0.5 * tmp4) * (1.0 + tmp3));
}
/*
template <typename Float>
Vector<Float, 3> fresnell::unpolarized<Float>::eval(Float cos_theta_d) const
{
	DJB_ASSERT(cos_theta_d >= 0.f && cos_theta_d <= 1.0 && "Invalid Angle");
	Vector<Float, 3>  F;

	//for (int i = 0; i < 3; ++i)
		F.x() = unpolarized__eval(cos_theta_d, Vector<Float, 3>::to_raw(ior).x());
		F.y() = unpolarized__eval(cos_theta_d, Vector<Float, 3>::to_raw(ior).y());
		F.z() = unpolarized__eval(cos_theta_d, Vector<Float, 3>::to_raw(ior).z());

	return {F.x(), F.y(), F.z()};
}*/
/*template <typename Float>
fresnel::schlick<Float>(const Vector<Float, 3>& f0): f0(f0)
{
}
template <typename Float>
Vector<Float, 3> schlick<Float>::eval(Float cos_theta_d) const
{
	DJB_ASSERT(cos_theta_d >= 0.f && cos_theta_d <= 1.0 && "Invalid Angle");
	Float c1 = 1.0 - cos_theta_d;
	Float c2 = c1 * c1;
	Float c5 = c2 * c2 * c1;

	return f0 + c5 * (Vector<Float, 3>(1) - f0);
}
template <typename Float>
Vector<Float, 3> sgd<Float>::eval(Float cos_theta_d) const
{
	DJB_ASSERT(cos_theta_d >= 0.f && cos_theta_d <= 1.0 && "Invalid Angle");
	Float c = cos_theta_d;

	return f0 - c * f1 + pow(1.0 - c, 5.0) * (Vector<Float, 3>(1) - f0);
}*/
/**
Vector<Float, 3> spline::eval(Float cos_theta_d) const
{
	DJB_ASSERT(cos_theta_d >= 0.f && cos_theta_d <= 1.0 && "Invalid Angle");
	Float u = 2.0 * acos(cos_theta_d) / M_PI;

	return spline::eval(m_points, spline::uwrap_edge, u);
}
*/
} // namespace fresnel

// *************************************************************************************************
// Microfacet API implementation

// -------------------------------------------------------------------------------------------------
/**
 * Compute the PDF parameters associated to the parameters of an ellipse
 */
template <typename Float>
static void
ellipse_to_pdfparams(
	Float a1, Float a2, Float phi_a,
	Float *ax, Float *ay, Float *rho
) {
	Float cos_phi_a = cos(phi_a);
	Float sin_phi_a = sin(phi_a);
	Float cos_2phi_a = 2.0 * cos_phi_a * cos_phi_a - 1.0f; // cos(2x) = 2cos(x)^2 - 1
	Float a1_sqr = a1 * a1;
	Float a2_sqr = a2 * a2;
	Float tmp1 = a1_sqr + a2_sqr;
	Float tmp2 = a1_sqr - a2_sqr;

	(*ax) = sqrt(0.5 * (tmp1 + tmp2 * cos_2phi_a));
	(*ay) = sqrt(0.5 * (tmp1 - tmp2 * cos_2phi_a));
	(*rho) = (a2_sqr - a1_sqr) * cos_phi_a * sin_phi_a / ((*ax) * (*ay));
}


// -------------------------------------------------------------------------------------------------
/**
 * Compute the ellipse parameters associated to the parameters of the PDF
 */
template <typename Float>
static void
pdfparams_to_ellipse(
	Float ax, Float ay, Float rho,
	Float *a1, Float *a2, Float *phi_a
) {
	Float ax_sqr = ax * ax;
	Float ay_sqr = ay * ay;
	Float cov = rho * ax * ay * /*saves 1 mul*/2.0;
	Float tmp1 = ax_sqr + ay_sqr;
	Float tmp2 = ax_sqr - ay_sqr;
	Float tmp3 = sqrt(tmp2 * tmp2 + cov * cov);

	(*a1) = sqrt(0.5 * (tmp1 + tmp3));
	(*a2) = sqrt(0.5 * (tmp1 - tmp3));
	(*phi_a) = (cov != 0.0) ? atan((ax_sqr - ay_sqr - tmp3) / cov) : 0.0;
}

//---------------------------------------------------------------------------
// Microfacet Params

/* Constructors */
template <typename Float>
microfacet<Float>::params::params(Float a1, Float a2, Float phi_a)
{
	set_ellipse(a1, a2, phi_a);
	set_location(0, 0);
}
template <typename Float>
microfacet<Float>::params::params(
	Float ax, Float ay, Float rho, Float tx_n, Float ty_n
) {
	set_pdfparams(ax, ay, rho, tx_n, ty_n);
}

/* Factories */
template <typename Float>
typename microfacet<Float>::params microfacet<Float>::params::standard()
{
	return microfacet::params::isotropic(1);
}
template <typename Float>
typename microfacet<Float>::params microfacet<Float>::params::isotropic(Float a)
{
	return microfacet::params::elliptic(a, a, 0);
}
template <typename Float>
typename microfacet<Float>::params
microfacet<Float>::params::elliptic(Float a1, Float a2, Float phi_a)
{
	return microfacet::params(a1, a2, phi_a);
}
template <typename Float>
typename microfacet<Float>::params
microfacet<Float>::params::pdfparams(
	Float ax, Float ay, Float rho,
	Float tx_n, Float ty_n
) {
	return microfacet::params(ax, ay, rho, tx_n, ty_n);
}

/* Mutators */
template <typename Float>
void microfacet<Float>::params::set_location(Float tx_n, Float ty_n)
{
	m_tx_n = tx_n;
	m_ty_n = ty_n;
	m_n = normalize(Vector3f(-tx_n, -ty_n, 1));
}
template <typename Float>
void microfacet<Float>::params::set_location(const Vector3f& n)
{
	m_n = n;
	m_tx_n = -n.x / n.z;
	m_ty_n = -n.y / n.z;
}
template <typename Float>
void microfacet<Float>::params::set_ellipse(Float a1, Float a2, Float phi_a)
{
	DJB_ASSERT(a1 > 0.0 && a2 > 0.0 && "Invalid ellipse radii");
	m_a1 = a1;
	m_a2 = a2;
	m_phi_a = phi_a;
	ellipse_to_pdfparams(a1, a2, phi_a, &m_ax, &m_ay, &m_rho);
	m_sqrt_one_minus_rho_sqr = sqrt(1.0 - m_rho * m_rho);
}
template <typename Float>
void
microfacet<Float>::params::set_pdfparams(
	Float ax, Float ay, Float rho,
	Float tx_n, Float ty_n
) {
	DJB_ASSERT(ax > 0.0 && ay > 0.0 && "Invalid scale parameters");
	DJB_ASSERT(fabs(rho) < 1.0 && "Invalid correlation parameter");
	m_ax = ax;
	m_ay = ay;
	m_rho = rho;
	m_sqrt_one_minus_rho_sqr = sqrt(1.0 - rho * rho);
	pdfparams_to_ellipse(ax, ay, rho, &m_a1, &m_a2, &m_phi_a);
	set_location(tx_n, ty_n);
}

/* Accessors */
template <typename Float>
void microfacet<Float>::params::get_location(Float *tx_n, Float *ty_n) const
{
	if (tx_n) (*tx_n) = m_tx_n;
	if (ty_n) (*ty_n) = m_ty_n;
}
template <typename Float>
void microfacet<Float>::params::get_location(Vector3f *n) const
{
	if (n) (*n) = m_n;
}
template <typename Float>
void
microfacet<Float>::params::get_ellipse(Float *a1, Float *a2, Float *phi_a) const
{
	if (a1) (*a1) = m_a1;
	if (a2) (*a2) = m_a2;
	if (phi_a) (*phi_a) = m_phi_a;
}
template <typename Float>
void
microfacet<Float>::params::get_pdfparams(
	Float *ax, Float *ay, Float *rho,
	Float *tx_n, Float *ty_n
) const {
	if (ax) (*ax) = m_ax;
	if (ay) (*ay) = m_ay;
	if (rho) (*rho) = m_rho;
	if (tx_n) (*tx_n) = m_tx_n;
	if (ty_n) (*ty_n) = m_ty_n;
}

//---------------------------------------------------------------------------
// Facet Ctor
/*
template <typename Float>
microfacet<Float>::microfacet(
	const typename mitsuba::fresnell::impl<Float>& fresnelll,
	bool shadow
):
	m_fresnel(fresnelll.copy()), m_shadow(shadow)
{
	
}
*/
//---------------------------------------------------------------------------
// Mutators
/*
template <typename Float>
void microfacet<Float>::set_fresnel(const fresnell::impl<Float>& f)
{
	delete m_fresnel;
	m_fresnel = f.copy();
}
*/
//---------------------------------------------------------------------------
// Evaluate f_r * cos
/*
template <typename Float>
Vector<Float, 3>
microfacet<Float>::evalp(const Vector<Float, 3>& i, const Vector<Float, 3>& o, const void *user_param) const
{
	const microfacet::params params = 
		user_param ? *reinterpret_cast<const microfacet::params *>(user_param)
		           : microfacet::params::standard();

	Vector3f h = normalize(i + o);
	Float G = gaf(h, i, o, params);

	if (G > 0.0) {
		Float cos_theta_d = sat(dot(o, h));
		Vector3f F = fresnel(cos_theta_d);
		Float D = ndf(h, params);

		return (F * ((D * G) / (4.0 * o.z())));
	}
	return Vector3f(0);
}

//---------------------------------------------------------------------------
// Evaluate f_r
template <typename Float>
Vector<Float, 3>
microfacet<Float>::eval(const Vector<Float, 3>& i, const Vector<Float, 3>& o, const void *user_param) const
{
	return (evalp(i, o, user_param) / i.z());
}
*/
//---------------------------------------------------------------------------
// Microfacet NDF
template <typename Float>
Float microfacet<Float>::ndf(const Vector3f& h, const microfacet<Float>::params& params) const
{ Mask  active= true;

  active &= h.z() > dr::Epsilon<Float>;	
	//if (h.z() > DJB_EPSILON) {
	Float cos_theta_h_sqr = h.z() * h.z();
	Float cos_theta_h_sqr_sqr = cos_theta_h_sqr * cos_theta_h_sqr;
	Float xslope = -h.x() / h.z();
	Float yslope = -h.y() / h.z();

	return dr::select(active,(p22(xslope, yslope, params) / cos_theta_h_sqr_sqr),0.0f);
	//}
	//return 0.0;
}

//---------------------------------------------------------------------------
// Microfacet Slope PDF
template <typename Float>
Float microfacet<Float>::p22(Float x, Float y, const params& params) const
{
	x-= params.m_tx_n;
	y-= params.m_ty_n;

	Float nrm = params.m_ax * params.m_ay * params.m_sqrt_one_minus_rho_sqr;
	Float x_ = x / params.m_ax;
	Float tmp1 = params.m_ax * y - params.m_rho * params.m_ay * x;
	Float tmp2 = params.m_ax * params.m_ay
	             * params.m_sqrt_one_minus_rho_sqr;
	Float y_ = tmp1 / tmp2;

	return (p22_std(x_, y_) / nrm);
}

//---------------------------------------------------------------------------
// Microfacet Visible Slope PDF
template <typename Float>
Float
microfacet<Float>::vp22(Float x, Float y, const Vector3f& k, const params& params) const
{
	Vector3f h = normalize(Vector<Float, 3>(-x, -y, 1));
	Float jacobian = h.z() * h.z() * h.z(); // solid angles to slope space jacobian

	return jacobian * vndf(h, k, params);
}

//---------------------------------------------------------------------------
// Microfacet VNDF
template <typename Float>
Float
microfacet<Float>::vndf(
	const Vector3f& h, const Vector3f& k,
	const microfacet<Float>::params& params
) const {MI_IMPORT_CORE_TYPES();
	Float kh = dot(k, h);
    Mask active= true;
	active &= kh > 0.0;
	if (dr::none_or<true>(active)) {
		Float D = ndf(h, params);

		return (kh * D / sigma(k, params));
	}
	return 0.0;
}

//---------------------------------------------------------------------------
// Smith Height Correlated Microfacet GAF
template <typename Float>
Float
microfacet<Float>::sigma(
	const Vector3f& k,
	const microfacet<Float>::params& params
) const {
	// warp the incident direction
	Float a = k.x() * params.m_ax + k.y() * params.m_ay * params.m_rho;
	Float b = k.y() * params.m_ay * params.m_sqrt_one_minus_rho_sqr;
	Float c = k.z() - k.x() * params.m_tx_n - k.y() * params.m_ty_n;
	Float nrm = sqrt(a * a + b * b + c * c);

	return nrm * sigma_std(Vector3f(a, b, c) / nrm);
}
template <typename Float>
Float
microfacet<Float>::g1(
	const Vector3f& h, const Vector3f& k,
	const microfacet<Float>::params& params
) const {MI_IMPORT_CORE_TYPES();
	Float g1_local_test = dot(k, params.m_n);
	Mask active= true;	
	active &= g1_local_test > 0.0;
	if (dr::none_or<true>(active))
		return k.z() / sigma(k, params);
	return 0.0;
}
template <typename Float>
Float
microfacet<Float>::gaf(
	const Vector3f& h, const Vector3f& i, const Vector3f& o,
	const microfacet<Float>::params& params
) const {MI_IMPORT_CORE_TYPES();
	Float G1_o = g1(h, o, params);
    Mask active= true;
	
	if (m_shadow) {
		Float G1_i = g1(h, i, params);
		Float tmp = G1_i * G1_o;
		active &= tmp > 0.0;
#if 1 // height correlated smith
		if (dr::none_or<true>(active))
			return (tmp / (G1_i + G1_o - tmp));

		return 0.0; // fully attenuated
#else // separable form
		return tmp;
#endif
	}

	return G1_o;
}

//---------------------------------------------------------------------------
// Importance Sample f_r * cos Using Two Uniform Numbers
/*
template <typename Float>
Vector<Float, 3>
microfacet<Float>::sample(
	Float u1, Float u2,
	const Vector3f& o,
	const void *user_param
) const {
	// get params
	const microfacet<Float>::params params = 
		user_param ? *reinterpret_cast<const microfacet<Float>::params *>(user_param)
		           : microfacet<Float>::params::standard();

#if 1 // clamp samples
	u1 = sat(u1) * (Float)0.99998 + (Float)0.00001;
	u2 = sat(u2) * (Float)0.99998 + (Float)0.00001;
#endif

	// warp the receiver direction
	Float a = o.x * params.m_ax + o.y * params.m_ay * params.m_rho;
	Float b = o.y * params.m_ay * params.m_sqrt_one_minus_rho_sqr;
	Float c = o.z - o.x * params.m_tx_n - o.y * params.m_ty_n;
	Vector3f o_std = normalize(Vector3f(a, b, c));

	if (o_std.z > 0.0) {
		// create a standard variate
		Float tx_m, ty_m;
		sample_vp22_std_smith(u1, u2, o_std, &tx_m, &ty_m);

		// warp the variate with the microfacet parameters
		Float tx_h = params.m_ax * tx_m + params.m_tx_n;
		Float choleski = params.m_rho * tx_m
			             + params.m_sqrt_one_minus_rho_sqr * ty_m;
		Float ty_h = params.m_ay * choleski + params.m_ty_n;

		// get the associated normal
		Vector3f h = normalize(Vector3f(-tx_h, -ty_h, 1));

		// return the reflected direction around the micronormal
		return (2.0 * dot(o, h) * h - o);
	}
	return Vector3f(0, 0, 1);
}
*/
//---------------------------------------------------------------------------
// Evaluate the PDF of a Sample
template <typename Float>
Float
microfacet<Float>::pdf(const Vector3f& i, const Vector3f& o, const void *user_param) const
{
	const microfacet<Float>::params params = 
		user_param ? *reinterpret_cast<const microfacet::params *>(user_param)
		           : microfacet::params::standard();
	Vector3f h = dr::normalize(i + o);
	Float G = gaf(h, i, o, params);
    Mask active= true;
	active &= G > 0.0;
	if (dr::none_or<true>(active)) {
		if (!supports_smith_vndf_sampling()) {
			return (h.z() * ndf(h, params) / (4.0 * dr::dot(i, h)));
		} else {
			return (vndf(h, o, params) / (4.0 * dr::dot(i, h)));
		}
	}
	return 0.0;
}

//---------------------------------------------------------------------------
// Evaluate f_r * cos / PDF
/*
template <typename Float>
Vector<Float, 3>
microfacet<Float>::evalp_is(
	Float u1, Float u2,
	const Vector3f& o,
	Vector3f *i, Float *pdf,
	const void *user_param
) const {
	const microfacet::params params = 
		user_param ? *reinterpret_cast<const microfacet::params *>(user_param)
		           : microfacet::params::standard();
	Vector3f i_ = sample(u1, u2, o, user_param);
	Vector3f h = normalize(i_ + o);
	Float G = gaf(h, i_, o, params);
	if (pdf) (*pdf) = 0.f;

	if (G > 0.0) {
		Float cos_theta_d = sat(dot(o, h));

		if (i) (*i) = i_;
		if (!supports_smith_vndf_sampling()) {
			Float pdf_ = (h.z * ndf(h, params) / (4.0 * cos_theta_d));
			if (pdf) (*pdf) = pdf_;
			return evalp(i_, o, user_param) / pdf_;
		} else {
			Vector3f F = fresnel(cos_theta_d);
			Float G1 = g1(h, o, params);
			if (pdf) (*pdf) = (vndf(h, o, params) / (4.0 * cos_theta_d));
			return (F * (G / G1));
		}
	}
	return Vector3f(0);
}
*/
//---------------------------------------------------------------------------
// Importance Sample a Smith VNDF
template <typename Float>
void
microfacet<Float>::sample_vp22_std_smith(
	Float u1, Float u2,
	const Vector3f& o,
	Float *xslope, Float *yslope
) const {
	if (supports_smith_vndf_sampling()) {
		(*xslope) = qf2(u1, o);
		(*yslope) = qf3(u2, o, *xslope);
	} else {
		sample_vp22_std_nmap(u1, u2, o, xslope, yslope);
	}
}
template <typename Float>
Float microfacet<Float>::qf2(Float u, const Vector3f& k) const 
{
	throw exc("djb_error: Not Implemented");
}
template <typename Float>
Float microfacet<Float>::qf3(Float u, const Vector3f& k, Float qf2) const
{
	throw exc("djb_error: Not Implemented");
}

// *************************************************************************************************
// Radial Microfacet API implementation
template <typename Float>
Float radial<Float>::p22_std(Float x, Float y) const
{
	return p22_radial(x * x + y * y);
}
template <typename Float>
Float radial<Float>::sigma_std(const Vector3f& k) const
{
	return sigma_std_radial(k.z());
}
template <typename Float>
void
radial<Float>::sample_vp22_std_nmap(
	Float u1, Float u2, const Vector3f& k,
	Float *xslope, Float *yslope
) const {
	Float phi_h = u1 * dr::Pi<Float> * 2.0;
	Float r_h = qf_radial(u2);

	(*xslope) = r_h * cos(phi_h);
	(*yslope) = r_h * sin(phi_h);
}
template <typename Float>
void
radial<Float>::sample_vp22_std_smith(
	Float u1, Float u2,
	const Vector3f& k,
	Float *xslope, Float *yslope
) const {
	if (supports_smith_vndf_sampling()) {
		Float cos_theta_k = k.z();
		Mask active= true;
		active &= k.z() < 1.0;
	//	Float sin_theta_k = k.z() < 1.0 ? dr::sqrt(1.0 - k.z() * k.z()) : 0.0;
		Float sin_theta_k =dr::select(active, dr::sqrt(1.0 - k.z() * k.z()),0.0);
		Float tx = qf2_radial(u1, cos_theta_k, sin_theta_k);
		Float ty = qf3_radial(u2, tx);

		if (sin_theta_k == 0.0) {
			/* perfect normal incidence */
			(*xslope) = tx;
			(*yslope) = ty;
		} else {
			/* nonnormal incidence */
			Float nrm = inversesqrt(k.x() * k.x() + k.y() * k.y());
			Float cos_phi_k = k.x() * nrm;
			Float sin_phi_k = k.y() * nrm;

			(*xslope) = cos_phi_k * tx - sin_phi_k * ty;
			(*yslope) = sin_phi_k * tx + cos_phi_k * ty;
		}
	} else {
		sample_vp22_std_nmap(u1, u2, k, xslope, yslope);
	}
}
template <typename Float>
Float
radial<Float>::qf2_radial(
	Float u,
	Float cos_theta_k,
	Float sin_theta_k
) const {
	throw exc("djb_error: Not Implemented");
}
template <typename Float>
Float radial<Float>::qf3_radial(Float u, Float qf2) const
{
	throw exc("djb_error: Not Implemented");
}


// *************************************************************************************************
// Beckmann Microfacet API implementation
template <typename Float>
Float beckmann<Float>::p22_radial(Float r_sqr) const
{
	return (dr::exp(-r_sqr) / dr::Pi<Float>);
}
template <typename Float>
Float beckmann<Float>::sigma_std_radial(Float cos_theta_k) const
{
	if (cos_theta_k == 1.0) return 1.0;
	Float sin_theta_k = dr::sqrt(1.0 - cos_theta_k * cos_theta_k);
	Float cot_theta_k = cos_theta_k / sin_theta_k;
	Float nu = cot_theta_k;
	Float tmp = dr::exp(-nu * nu) * inversesqrt(dr::Pi<Float>);
	return (cos_theta_k * (1.0 + erf(nu)) + sin_theta_k * tmp) / 2.0;
}
template <typename Float>
Float beckmann<Float>::cdf_radial(Float r) const
{
	return (1.0 - dr::exp(-r * r));
}
template <typename Float>
Float beckmann<Float>::qf_radial(Float u) const
{
	return dr::sqrt(-dr::log(1.0 - u));
}
template <typename Float>
Float beckmann<Float>::qf1(Float u) const
{
	return erfinv(2.0 * u - 1.0);
}

/* Based on Mitsuba's source code, from Wenzel Jakob */
template <typename Float>
Float beckmann<Float>::qf2_radial(Float u, Float cos_theta_k, Float sin_theta_k) const
{
	const Float sqrt_pi_inv = 1. / dr::sqrt(dr::Pi<Float>);

	/* The original inversion routine from the paper contained
	   discontinuities, which causes issues for QMC integration
	   and techniques like Kelemen-style MLT. The following code
	   performs a numerical inversion with better behavior */
	Float cot_theta_k = cos_theta_k / sin_theta_k;
	Float tan_theta_k = sin_theta_k / cos_theta_k;

	/* Search interval -- everything is parameterized
	   in the erf() domain */
	Float a = -1, c = erf(cot_theta_k);
	u = max(u, (Float)1e-6f);

	/* Start with a good initial guess */
	/* We can do better (inverse of an approximation computed in Mathematica) */
	Float fit = 1 + cos_theta_k 
	            * (-0.876f + cos_theta_k * (0.4265f - 0.0594f * cos_theta_k));
	Float b = c - (1 + c) * std::pow(1 - u, fit);

	/* Normalization factor for the CDF */
	Float normalization = 1 / (1 + c + sqrt_pi_inv *
		tan_theta_k * exp(-cot_theta_k * cot_theta_k));

	uint32_t it = 0;
	while (++it < 10) {
		/* Bisection criterion -- the oddly-looking
		   boolean expression are intentional to check
		   for NaNs at little additional cost */
		if (!(b >= a && b <= c))
			b = 0.5f * (a + c);

		/* Evaluate the CDF and its derivative
		   (i.e. the density function) */
		Float invErf = erfinv(b);
		Float value = normalization * (1 + b + sqrt_pi_inv *
			tan_theta_k * std::exp(-invErf * invErf)) - u;
		Float derivative = normalization * (1 - invErf * tan_theta_k);

		if (fabs(value) < 1e-5f)
			break;

		/* Update bisection intervals */
		if (value > 0)
			c = b;
		else
			a = b;

		b -= value / derivative;
	}

	/* Now convert back into a slope value */
	return erfinv(max(-(Float)0.9999, b));
}
template <typename Float>
Float beckmann<Float>::qf3_radial(Float u, Float qf2) const
{
	return qf1(u);
}
template <typename Float>
beckmann<Float>::lrep::lrep(
	Float E1, Float E2, Float E3, Float E4, Float E5
): m_E1(E1), m_E2(E2), m_E3(E3), m_E4(E4), m_E5(E5) {
	/* empty */
}
/*
template <typename Float>
void beckmann<Float>::params_to_lrep(const microfacet<Float>::params& params, beckmann<Float>::lrep *lrep)
{
	DJB_ASSERT(lrep && "Null output ptr");
	Float ax, ay, rho, tx_n, ty_n;
	params.get_pdfparams(&ax, &ay, &rho, &tx_n, &ty_n);
	(*lrep) = beckmann::lrep(tx_n, ty_n,
	                         (Float)0.5 * ax * ax + tx_n * tx_n,
	                         (Float)0.5 * ay * ay + ty_n * ty_n,
	                         (Float)0.5 * rho * ax * ay + tx_n * ty_n);
}
template <typename Float>*//*
void beckmann<Float>::lrep_to_params(const lrep& lrep, microfacet::params *params)
{
	DJB_ASSERT(params && "Null output ptr");
	Float tx_n = lrep.m_E1;
	Float ty_n = lrep.m_E2;
	Float tmp1 = max((Float)0.0, lrep.m_E3 - lrep.m_E1 * lrep.m_E1);
	Float tmp2 = max((Float)0.0, lrep.m_E4 - lrep.m_E2 * lrep.m_E2);
	/* clamp the parameters to valid values */
	/*
	Float ax = (Float)max(1e-5, sqrt(2.0 * tmp1));
	Float ay = (Float)max(1e-5, sqrt(2.0 * tmp2));
	Float rho = (Float)2.0 * (lrep.m_E5 - lrep.m_E1 * lrep.m_E2) / (ax * ay);
	rho = min((Float)0.99, max((Float)-0.99, rho));

	(*params) = microfacet::params::pdfparams(ax, ay, rho, tx_n, ty_n);*/
//}
/*
template <typename Float>
beckmann::lrep beckmann::lrep::operator+(const lrep& r) const
{
	return beckmann::lrep(m_E1 + r.m_E1,
	                      m_E2 + r.m_E2,
	                      m_E3 + r.m_E3 + (Float)2.0 * m_E1 * r.m_E1,
	                      m_E4 + r.m_E4 + (Float)2.0 * m_E2 * r.m_E2,
	                      m_E5 + r.m_E5 + m_E1 * r.m_E2 + m_E2 * r.m_E1);
}
template <typename Float>
beckmann::lrep beckmann::lrep::operator*(Float sc) const
{
	DJB_ASSERT(sc >= (Float)0.0 && "Invalid scale");
	Float sc_sqr = sc * sc;

	return beckmann::lrep(m_E1 * sc, m_E2 * sc,
	                      m_E3 * sc_sqr, m_E4 * sc_sqr,
	                      m_E5 * sc_sqr);
}
template <typename Float>
beckmann::lrep& beckmann::lrep::operator+=(const lrep& rep)
{
	m_E1+= rep.m_E1;
	m_E2+= rep.m_E2;
	m_E3+= rep.m_E3 + (Float)2.0 * m_E1 * rep.m_E1;
	m_E4+= rep.m_E4 + (Float)2.0 * m_E2 * rep.m_E2;
	m_E5+= rep.m_E5 + m_E1 * rep.m_E2 + m_E2 * rep.m_E1;

	return (*this);
}
template <typename Float>
beckmann::lrep& beckmann::lrep::operator*=(Float sc)
{
	DJB_ASSERT(sc >= (Float)0.0 && "Invalid scale");
	Float sc_sqr = sc * sc;
	m_E1*= sc;
	m_E2*= sc;
	m_E3*= sc_sqr;
	m_E4*= sc_sqr;
	m_E5*= sc_sqr;

	return (*this);
}*/
template <typename Float>
void beckmann<Float> ::lrep::shear(Float tx, Float ty)
{
	m_E1+= tx;
	m_E2+= ty;
	m_E3+= tx * tx;
	m_E4+= ty * ty;
	m_E5+= tx * ty;
}
template <typename Float>
void beckmann<Float> ::lrep::scale(Float x, Float y)
{
	m_E1*= x;
	m_E2*= y;
	m_E3*= x * x;
	m_E4*= y * y;
	m_E5*= x * y;
}

// *************************************************************************************************
// GGX Microfacet API implementation
template <typename Float>
Float ggx<Float> ::p22_radial(Float r_sqr) const
{
	Float tmp = 1.0 + r_sqr;
	return (1.0 / (dr::Pi<Float> * tmp * tmp));
}
template <typename Float>
Float ggx<Float> ::sigma_std_radial(Float cos_theta_k) const
{
	return ((1.0 + cos_theta_k) / 2.0);
}
template <typename Float>
Float ggx<Float> ::cdf_radial(Float r) const
{
	Float tmp = r * r;
	return (tmp / (1.0 + tmp));
}

template <typename Float>

Float ggx<Float> ::qf_radial(Float u) const
{
	return sqrt(u / (1.0 - u));
}

template <typename Float>
Float ggx<Float> ::qf1(Float u) const
{
	if (u < 0.5) {
		u = (0.5 - u) * 2.0;
		return -u * inversesqrt(1.0 - u * u);
	} else {
		u = (u - 0.5) * 2.0;
		return u * inversesqrt(1.0 - u * u);
	}
}
template <typename Float>
Float ggx<Float> ::qf2_radial(Float u, Float cos_theta_k, Float sin_theta_k) const
{   Mask active = true;
	DRJIT_MARK_USED(active);
	Float sin_theta = u * (1.0 + cos_theta_k) - 1.0;
	Float cos_theta = dr::sqrt(1.0 - sin_theta * sin_theta);
    Mask active_cos = cos_theta > 0.707107;   
	if (dr::none_or<true>(active_cos)) {
		Float tan_theta = sin_theta / cos_theta;
        Mask active_sin = sin_theta_k < 0.707107;
		if (dr::none_or<true>(active_sin)) {
			Float tan_theta_k = sin_theta_k / cos_theta_k;

			return -(tan_theta + tan_theta_k) / (1.0 - tan_theta * tan_theta_k);
		} else {
			Float cot_theta_k = cos_theta_k / sin_theta_k;

			return (1.0 + tan_theta * cot_theta_k) / (tan_theta - cot_theta_k);
		}
	} else {
		Float cot_theta = cos_theta / sin_theta;
        Mask active_sin_k = sin_theta_k < 0.707107;
		if (dr::none_or<true>(active_sin_k)) {
			Float tan_theta_k = sin_theta_k / cos_theta_k;

			return (1.0 + tan_theta_k * cot_theta) / (tan_theta_k - cot_theta);
		} else {
			Float cot_theta_k = cos_theta_k / sin_theta_k;

			return (cot_theta + cot_theta_k) / (1.0 - cot_theta * cot_theta_k);
		}
	}
}
template <typename Float>
Float ggx<Float> ::qf3_radial(Float u, Float qf2) const
{
	Float alpha = sqrt(1.0 + qf2 * qf2);
	Float S = 0;
    Mask active=true;
	active &= u < 0.5;
	if (dr::none_or<true>(active)) {
		u = 2.0 * (0.5 - u); // remap to (0,1)
		S = -1.0;
	} else {
		u = 2.0 * (u - 0.5); // remap to (0,1)
		S = 1.0;
	}

	return S * alpha * qf3_rational_approx(u);
}

/* Based on Mitsuba's source code, from Wenzel Jakob */
template <typename Float>
Float ggx<Float> ::qf3_rational_approx(Float u) const
{
	Float p = u * (u * (u * (-0.365728915865723) 
	          + 0.790235037209296) - 0.424965825137544) + 0.000152998850436920;
	Float q = u * (u * (u * (u * 0.169507819808272 - 0.397203533833404) 
	          - 0.232500544458471) + 1) - 0.539825872510702;

	return (p / q);
}

// *************************************************************************************************
// Tabulated Microfacet API implementation


template <typename Float>
Float tabular<Float> ::p22_radial(Float r_sqr) const
{
	Float r = dr::sqrt(r_sqr);
	Float u = dr::sqrt(2.0 * dr::atan(r) / dr::Pi<Float>);
	return spline::eval<Float>(m_p22, "edge", u);
}
template <typename Float>
Float tabular<Float> ::sigma_std_radial(Float cos_theta_k) const
{
	Float u = 2.0 * dr::acos(cos_theta_k) / dr::Pi<Float>;
	return spline::eval<Float>(m_sigma, "edge", u);
}
template <typename Float>
Float tabular<Float> ::cdf_radial(Float r) const
{
	Float u = dr::atan(r) * 2.0 / dr::Pi<Float>;
	
	u=dr::maximum(0.0f, u);
	//if (u < 0.0) u = 0.0;
	return spline::eval<Float>(m_cdf, "edge", dr::sqrt(u));
}
template <typename Float>
Float tabular<Float> ::qf_radial(Float u) const
{
	DJB_ASSERT(u > 0.0 && u < 1.0);
	Float qf = spline::eval<Float>(m_qf,"edge",  u);
	return dr::tan(qf * dr::Pi<Float> / 2.0);
}
template <typename Float>
Float tabular_anisotropic<Float> ::p22_std(Float x, Float y) const
{
	Float theta = atan(sqrt(x * x + y * y));
	Float phi   = atan2(-y, -x);
	return p22_std_theta_phi(theta, phi);
}
template <typename Float>
Float tabular_anisotropic<Float> ::p22_std_theta_phi(Float theta, Float phi) const
{
	if (phi < 0.0) phi+= 2.0 * dr::Pi<Float>;
	Float u1 = theta * 2.0 / dr::Pi<Float>;
	Float u2 = phi * 0.5 / dr::Pi<Float>;
	uint32_t w = m_elevation_res;
	uint32_t h = m_azimuthal_res;

	return spline::eval2d<Float>(m_p22, w, h, 
	                      "edge", u1,
	                      "repeat",u2);
}
template <typename Float>
Float tabular_anisotropic<Float> ::sigma_std(const Vector3f& k) const
{
	Float theta = acos(k.z);
	Float phi = atan2(k.y, k.x);
	if (phi < 0.0) phi+= 2.0 * dr::Pi<Float>;
	Float u1 = theta * 2.0 / dr::Pi<Float>;
	Float u2 = phi * 0.5 / dr::Pi<Float>;
	uint32_t w = m_elevation_res;
	uint32_t h = m_azimuthal_res;

	return spline::eval2d<Float>(m_sigma, w, h, 
	                      "edge", u1,
	                      "repeat", u2);
}

// -------------------------------------------------------------------------------------------------
// Constructor
/*
tabular::tabular(const brdf& brdf, int res, bool shadow):
	radial()
{
	DJB_ASSERT(res > 2 && "Invalid Resolution");
	set_shadow(shadow);

	// allocate memory
	m_p22.reserve(res);
	m_sigma.reserve(res);
	m_cdf.reserve(res);
	m_qf.reserve(res);

	// eval
	compute_p22_smith(brdf, res);
	normalize_p22();
	compute_sigma();
	compute_fresnel(brdf, res);

	// sample
	compute_cdf();
	compute_qf();
}*/
template <typename Float>
tabular_anisotropic<Float>::tabular_anisotropic(
	const brdf<Float>& brdf,
	uint32_t elevation_res,
	uint32_t azimuthal_res,
	bool shadow
): microfacet<Float>() {
	DJB_ASSERT(elevation_res > 1 && azimuthal_res > 1 && "Invalid Resolution");
//	m_shadow = shadow;

	m_elevation_res = elevation_res;
	m_azimuthal_res = azimuthal_res;

	// allocate memory
	m_p22.reserve(elevation_res * azimuthal_res);
	m_sigma.reserve(elevation_res * azimuthal_res);
	m_pdf1.reserve(azimuthal_res);
	m_cdf1.reserve(azimuthal_res);
	m_qf1.reserve(azimuthal_res);
	m_pdf2.reserve(elevation_res * azimuthal_res);
	m_cdf2.reserve(elevation_res * azimuthal_res);
	m_qf2.reserve(elevation_res * azimuthal_res);
	/*
	m_p22=dr::zeros<FloatStorage>(elevation_res * azimuthal_res);
	m_sigma=dr::zeros<FloatStorage>(elevation_res * azimuthal_res);
	m_pdf1=dr::zeros<FloatStorage>(azimuthal_res);
	m_cdf1=dr::zeros<FloatStorage>(azimuthal_res);
	m_qf1=dr::zeros<FloatStorage>(azimuthal_res);
	m_pdf2=dr::zeros<FloatStorage>(elevation_res * azimuthal_res);
	m_cdf2=dr::zeros<FloatStorage>(elevation_res * azimuthal_res);
	m_qf2=dr::zeros<FloatStorage>(elevation_res * azimuthal_res);
*/
	


	// eval
	compute_p22_smith(brdf);
	normalize_p22();
	compute_sigma();
	compute_fresnel(brdf, elevation_res);

	// sample
	compute_pdf1();
	compute_cdf1();
	compute_qf1();
	compute_pdf2();
	compute_cdf2();
	compute_qf2();
}

// -------------------------------------------------------------------------------------------------
// Normalize the Slope PDF
template <typename Float>
void tabular<Float>::normalize_p22()
{
	const int ntheta = 128;
	const float dphi = 2.0 * M_PI;
	const float dtheta = M_PI / (float_t)ntheta;
	Float nint = 0.0;

	for (int i = 0; i < ntheta; ++i) {
		float u = (float)i / (float)ntheta; // in [0,1)
		float theta_h = u * u * M_PI * 0.5;
		float r_h = tan(theta_h);
		float cos_theta_h = cos(theta_h);
		Float p22_r = p22_radial(r_h * r_h);

		nint+= (u * p22_r * r_h) / (cos_theta_h * cos_theta_h);
	}
	nint*= dtheta * dphi;

	// normalize the slope pdf
	DJB_ASSERT(nint > 0.0); // should never fail
	nint = 1.0 / nint;
	for (int i = 0; i < (int)m_p22.size(); ++i)
		m_p22[i]*= nint;

#ifndef NVERBOSE
	Throw("djb_verbose: Slope PDF norm. constant = %.f", nint);
#endif
}
/*template <typename Float>
void tabular<Float>::normalize_p22()
{
	const int ntheta = 128;
	double dphi = 2.0 * M_PI;
	const double dtheta = M_PI / ntheta;
	Float nint = 0.0;

	for (int i = 0; i < ntheta; ++i) {
		double u = i / ntheta; // in [0,1)
		double theta_h = u * u * M_PI * 0.5;
		double r_h = dr::tan(theta_h);
		double cos_theta_h = dr::cos(theta_h);
		Float p22_r = p22_radial(r_h * r_h);
	   
       
		nint= nint+(u * p22_r * r_h) / (cos_theta_h * cos_theta_h);
	} 
	
	nint*= dtheta * dphi;
   
	// normalize the slope pdf
	DJB_ASSERT(nint > 0.0); // should never fail
	nint = 1.0 / nint;
	for (uint32_t i = 0; i < m_p22.size(); ++i)
		m_p22[i]=m_p22[i]* nint;

#ifndef NVERBOSE
  Throw("djb_verbose: Slope PDF norm. constant = %.f", nint);
	//DJB_LOG("djb_verbose: Slope PDF norm. constant = %.9f\n", nint);
#endif
}*/
template <typename Float>
void tabular_anisotropic<Float>::normalize_p22()
{
	const uint32_t ntheta = 128;
	const uint32_t nphi   = 256;
	Float dtheta = sqrt(0.5 * dr::Pi<Float>) / (Float)ntheta;
	Float dphi   = 2.0 * dr::Pi<Float> / (Float)nphi;
	Float k = 0.0;

	for (uint32_t j = 0; j < nphi; ++j) {
		Float u = (Float)j / (Float)nphi;
		Float phi = u * 2.0 * dr::Pi<Float>;
		for (uint32_t i = 0; i < ntheta; ++i) {
			Float u = (Float)i / (Float)ntheta; // in [0,1)
			Float theta = u * sqrt(dr::Pi<Float> * 0.5);
			Float theta_sqr = theta * theta; // in [0, pi/2)
			Float c = cos(theta_sqr);
			Float pdf = p22_std_theta_phi(theta_sqr, phi);
			Float weight = (theta * tan(theta_sqr)) / (c * c);

			k+= (weight * pdf);
		}
	}
	k*= (2.0 * dtheta * dphi);

	// normalize the slope pdf
	k = 1.0 / k;
	for (uint32_t i = 0; i < (uint32_t)m_p22.size(); ++i)
		m_p22[i]*= k;

#ifndef NVERBOSE
	Throw("djb_verbose: Anisotropic slope PDF norm. constant = %.f", k);
	//DJB_LOG("djb_verbose: Anisotropic slope PDF norm. constant = %.9f\n", k);
#endif
}

// -------------------------------------------------------------------------------------------------
/**
 * Compute the Smith Masking Term
 *
 * This requires the computation of an integral.
 * The resolution at which the integral is computed is hard coded to 
 * yield sufficient precision for most distributions.
 */
template <typename Float>
void tabular<Float>::compute_sigma()
{
	const uint32_t ntheta = 90;
	const uint32_t nphi   = 180;
	Float dtheta = dr::Pi<Float> / ntheta;
	Float dphi   = 2.0 * dr::Pi<Float> / nphi;
	uint32_t cnt = m_p22.size() - 1;
   
	for (uint32_t i = 0; i < cnt; ++i) {
		Float tmp = i / cnt; // in [0, 1)
		Float theta_k = tmp * 0.5 * dr::Pi<Float>; // in [0, pi/2)
		Float cos_theta_k = dr::cos(theta_k);
		Float sin_theta_k = dr::sin(theta_k);
		Float nint = 0.0;

		for (uint32_t j2 = 0; j2 < nphi; ++j2) {
			Float u_j = (Float)j2 / (Float)nphi; // in [0, 1)
			Float phi_h = u_j * 2.0 * dr::Pi<Float>;          // in [0, 2pi)
			for (uint32_t j1 = 0; j1 < ntheta; ++j1) {
				Float u_i = (Float)j1 / (Float)ntheta; // in [0, 1)
				Float theta_h = u_i * u_i * dr::Pi<Float> * 0.5; // in [0, sqrt(pi/2))
				Float sin_theta_h = dr::sin(theta_h);
				Float kh = sin_theta_k * sin_theta_h * dr::cos(phi_h)
				           + cos_theta_k * dr::cos(theta_h);
                Float ss=dr::sin(theta_h);
                Vector3f xyz(ss*dr::cos(phi_h),ss*dr::sin(phi_h),dr::cos(theta_h)); 
				nint+= dr::maximum(kh, 0.0f)
				     *microfacet<Float>::ndf(xyz)
				     * u_i * sin_theta_h;
			}
		}
		nint*= dtheta * dphi;
		//m_sigma[i]=dr::maximum((cos_theta_k, nint));
		m_sigma.push_back(dr::maximum(nint,(Float)cos_theta_k));
	}
	m_sigma.push_back(m_sigma.back());
	//m_sigma

#ifndef NVERBOSE
    Throw("djb_verbose: Projected area term ready");
	//DJB_LOG("djb_verbose: Projected area term ready\n");
#endif
}
template <typename Float>
void tabular_anisotropic<Float>::compute_sigma()
{
	const uint32_t ntheta = 45;
	const uint32_t nphi   = 90;
	Float dtheta = dr::sqrt(dr::Pi<Float> * 0.5) / (Float)ntheta;
	Float dphi   = 2.0 * dr::Pi<Float> / (Float)nphi;
	uint32_t w = m_elevation_res - 1;
	uint32_t h = m_azimuthal_res;

	for (uint32_t i2 = 0; i2 < h; ++i2) {
		Float tmp = (Float)i2 / (Float)h; // in [0, 1)
		Float phi_k = tmp * 2.0 * dr::Pi<Float>; // in [0, 2pi)
		for (uint32_t i1 = 0; i1 < w; ++i1) {
			Float tmp = (Float)i1 / (Float)w; // in [0, 1)
			Float theta_k = tmp * 0.5 * dr::Pi<Float>; // in [0, pi/2)
			Float cos_theta_k = cos(theta_k);
			Float nint = 0.0;

			for (uint32_t j2 = 0; j2 < nphi; ++j2) {
				Float tmp = (Float)j2 / (Float)nphi; // in [0, 1)
				Float phi = tmp * 2.0 * dr::Pi<Float>;          // in [0, 2pi)
				for (uint32_t j1 = 0; j1 < ntheta; ++j1) {
					Float tmp = (Float)j1 / (Float)ntheta; // in [0, 1)
					Float theta = tmp * sqrt(dr::Pi<Float> * 0.5); // in [0, sqrt(pi/2))
					Float theta_sqr = theta * theta; // in [0, pi/2)
					Float sin_theta = sin(theta_sqr);
					Float m_dot_k = sin(theta_k) * sin_theta * cos(phi - phi_k)
					               + cos_theta_k * cos(theta_sqr);
					Float weight = theta * sin_theta;
					Float masking = max((Float)0.0, m_dot_k)
					               * ndf(Vector3f(theta_sqr, phi));

					nint+= weight * masking;
				}
			}
			nint*= 2.0 * dtheta * dphi;
			m_sigma.push_back(max((Float)cos_theta_k, nint));
		}
		m_sigma.push_back(m_sigma.back());
	}

#ifndef NVERBOSE
	Throw("djb_verbose: Anisotropic projected area term ready");
//	DJB_LOG("djb_verbose: Anisotropic projected area term ready\n");
#endif
}

// -------------------------------------------------------------------------------------------------
/**
 * Evaluate the Components from a BRDF (Smith)
 *
 * The NDF and the Fresnel term are extracted from the data.
 * This inversion requires exponentiating a matrix, so a small, self-contained
 * row major matrix API is implemented first.
 */
template <typename Float>
class matrix {

	std::vector<Float> mij;
	
	uint32_t size;
public:

  
 
	matrix(uint32_t size);

	Float& operator()(uint32_t i, uint32_t j) {return mij[j*size+i];}
    //Float & operator()(uint32_t i, uint32_t j) {return mij[j*size+i];}
	void transform( std::vector<Float>& v, std::vector<Float>& out);
	std::vector<Float> eigenvector(uint32_t iterations) {
	uint32_t j = 0;
    std::vector<Float> vec[2];
	//vec[0].reserve(size);
	//vec[1].reserve(size);
	vec[0].resize(size, 0.0f);
    vec[1].resize(size, 0.0f);
	for (uint32_t i = 0; i < size; ++i)
		vec[j][i] = 1.0f;
		
	for (uint32_t i = 0; i < iterations; ++i) {
		transform(vec[j], vec[1-j]);
		j = 1 - j;
	}
	return vec[j];
}
};
template <typename Float>
matrix<Float>::matrix(uint32_t size) : mij (size * size, 0.0),size(size)
{}
template <typename Float>
void matrix<Float>::transform( std::vector<Float>& v, std::vector<Float>& out) 
{   /// Const pointer to the start of the tensor 
	out.resize(0);
	for (uint32_t j = 0; j < size; ++j) {
		out.push_back(0);
		for (uint32_t i = 0; i < size; ++i) {
		    out[j]+= (*this)(i, j) * v[i];
		}
	}
}
/*
template <typename Float>
FloatStorage matrix<Float>::eigenvector(uint32_t iterations) const
{
	uint32_t j = 0;
	FloatStorage vec[2];
	vec[0].reserve(size);
	vec[1].reserve(size);
	for (uint32_t i = 0; i < size; ++i)
		vec[j][i] = 1.0;
	for (uint32_t i = 0; i < iterations; ++i) {
		transform(vec[j], vec[1-j]);
		j = 1 - j;
	}
	return vec[j];
}*/

template <typename Float>
void tabular<Float>::compute_p22_smith(const brdf<Float>& brdf, uint32_t res)
{
	uint32_t cnt = res - 1;
	Float dtheta = dr::sqrt(dr::Pi<Float> * 0.5) / cnt;
	matrix<Float> km(cnt);

	for (uint32_t i = 0; i < cnt; ++i) {
		Float tmp = i / cnt;
		Float theta = tmp * dr::sqrt(dr::Pi<Float> * 0.5);
		Float theta_o = theta * theta;
		Float cos_theta_o = dr::cos(theta_o);
		Float tan_theta_o = dr::tan(theta_o);
        Float ss=dr::sin(theta_o);
        Vector3f xyz(ss*dr::cos(0.0),ss*dr::sin(0.0),dr::cos(theta_o));
		auto fr = brdf.eval(xyz, xyz);
		Float fr_i = intensity(fr);
		Float kji_tmp = (dtheta * dr::pow(cos_theta_o, 6.0)) * (8.0 * fr_i);

		for (uint32_t j = 0; j < cnt; ++j) {
			const Float dphi_h = dr::Pi<Float> / 180.0;
			Float tmp = j / cnt;
			Float theta = tmp * dr::sqrt(dr::Pi<Float> * 0.5);
			Float theta_h = theta * theta;
			Float cos_theta_h = dr::cos(theta_h);
			Float tan_theta_h = dr::tan(theta_h);
			Float tan_product = tan_theta_h * tan_theta_o;
			Float nint = 0.0f;
            uint32_t phi_dh=360;
			Float phi_hh = 0.0f;
			for(uint32_t jj = 0; jj < phi_dh; ++jj)
			   { phi_hh += dphi_h;
				nint+= dr::maximum(tan_product * dr::cos(phi_hh),1.0f);}
		//	for ( Float phi_h = 0.0f; phi_h < 2.0f * dr::Pi<Float>; phi_h+= dphi_h)
		//		nint+= dr::maximum(tan_product * dr::cos(phi_h),1.0f);
			nint*= dphi_h;

			km(j, i) = theta * kji_tmp * nint * tan_theta_h
			         / (cos_theta_h * cos_theta_h);
		}
	}

	// compute slope pdf
	const std::vector<Float> v = km.eigenvector(4);
//    const std::vector<Float> v={0.0f, 1.0f, 2.0f, 3.0f};
	for (uint32_t i = 0; i < (uint32_t)v.size(); ++i)
		m_p22.push_back(1e-2 * v[i]);
	m_p22.push_back(0);
}

template <typename Float>
void tabular_anisotropic<Float>::compute_p22_smith(const brdf<Float>& brdf)
{
	uint32_t w = m_elevation_res - 1;
	uint32_t h = m_azimuthal_res;
	Float dtheta = dr::sqrt(dr::Pi<Float> * 0.5) / (Float)w;
	Float dphi = 2.0 * dr::Pi<Float> / (Float)h;
	std::vector<double> kv(w * h);
	matrix<Float> km(w * h);

	// compute kernel matrix terms
	// compute Smith matrix
	for (uint32_t i2 = 0; i2 < h; ++i2)
	for (uint32_t i1 = 0; i1 < w; ++i1) {
		Float tmp1 = i1 / w; // in [0, 1)
		Float tmp2 = i2 / h; // in [0, 1)
		Float theta = tmp1 * 0.5 * dr::Pi<Float>; // in [0, pi/2)
		Float phi = tmp2 * 2.0 * dr::Pi<Float>; // in [0, 2 pi)
		Float sin_theta = dr::sin(theta);
		Float zo = dr::cos(theta);
		Float xo = sin_theta * dr::cos(phi);
		Float yo = sin_theta * dr::sin(phi);
		Float ss=dr::sin(theta);
        Vector3f xyz(ss*dr::cos(phi),ss*dr::sin(phi),dr::cos(theta));
		Vector3f fr = brdf.eval(xyz, xyz);
		Float fr_i = fr.intensity();
		Float kji_tmp1 = (dtheta * dphi) * (4.0 * fr_i * dr::pow(zo, 5.0));

		for (uint32_t j2 = 0; j2 < h; ++j2)
		for (uint32_t j1 = 0; j1 < w; ++j1) {
			Float tmp1 = j1 / w; // in [0, 1)
			Float tmp2 = j2 / h; // in [0, 1)
			Float theta = tmp1 * 0.5 * dr::Pi<Float>; // in [0, pi/2)
			Float phi = tmp2 * 2.0 * dr::Pi<Float>; // in [0, 2 pi)
			Float cos_theta = dr::cos(theta);
			Float tan_theta = dr::tan(theta);
			Float slope1 = -tan_theta * dr::cos(phi); // x slope
			Float slope2 = -tan_theta * dr::sin(phi); // y slope
			Float m_dot_o = zo - xo * slope1 - yo * slope2;
			Float den = cos_theta * cos_theta;
			Float kji_tmp2 = tan_theta * dr::maximum(m_dot_o,0.0) / den;

			km(j2 * w + j1, i2 * w + i1) = kji_tmp1 * kji_tmp2;
		}
	}

	// compute slope pdf
	kv = km.eigenvector(4);

	// compute pdf
	for (uint32_t j = 0; j < h; ++j) {
		for (uint32_t i = 0; i < w; ++i) {
			Float pdf = kv[j * w + i];
			m_p22.push_back(pdf);
		}
		m_p22.push_back(0.0);
	}
}

// -------------------------------------------------------------------------------------------------
// Fresnel computation
template <typename Float>
void tabular<Float>::compute_fresnel(const brdf<Float>& brdf, uint32_t res)
{
	std::vector<Vector3f> fresnel(res);
	uint32_t cnt = res - 1;

	// compute average ratio between input and Torrance Sparrow equation
	for (uint32_t i = 0; i < cnt; ++i) {
		const Float phi_d = dr::Pi<Float> * 0.5; // this can NOT be tweaked
		const Float phi_h = 0.0;        // this can be tweaked (no impact on MERL)
		float tmp = i /cnt;
		float theta_d = tmp * M_PI * 0.5; // linear parameterization in theta_d
		Vector3f f = Vector3f(0); // Fresnel value at R, G, B wavelengths
		uint32_t count[3] = {0, 0, 0};
		float theta_h = 0.0;
		
        Mask active=true;
		DRJIT_MARK_USED(active);
		for (uint32_t j = 0; theta_h < M_PI * 0.5 - theta_d; ++j) {
			float tmp1 = j / cnt;
			theta_h = tmp1 * tmp1 * M_PI * 0.5; // in [0, pi/2)
         
			if (theta_h > M_PI * 0.5) continue;
            Float ss=dr::sin(theta_h);
            Vector3f xyz(ss*dr::cos(phi_h),ss*dr::sin(phi_h),dr::cos(theta_h));
			Vector3f dir_h = xyz;
			Float ss1=dr::sin(theta_d);
            Vector3f xyzz(ss1*dr::cos(phi_d),ss1*dr::sin(phi_d),dr::cos(theta_d));
			Vector3f dir_d = xyzz;
			Vector3f dir_i, dir_o;
	        mitsuba::brdf<Float>::hd_to_io(dir_h, dir_d, &dir_i, &dir_o);

			dir_i = Vector3f(0, 0, 1); //XXX: hack to reproduce my EGSR fits
			Vector3f fr1 = brdf.eval(dir_i, dir_o);
			Vector3f fr2 =microfacet<Float>::eval(dir_i, dir_o);
             Mask active0 =active && fr2.x() > /* epsilon */dr::Epsilon<Float>;
			if (dr::none_or<true>(active0)) {
				Float ratio0 = fr1.x() / fr2.x();
				f.x()+= ratio0;
				++count[0];
			}
            Mask active1 =active && fr2.y() > /* epsilon */dr::Epsilon<Float>;
			if (dr::none_or<true>(active1)) {
				Float ratio1 = fr1.y() / fr2.y();
				f.y()+= ratio1;
				++count[1];
			}
			Mask active2 =active && fr2.z() > /* epsilon */dr::Epsilon<Float>;
			if (dr::none_or<true>(active2)) {
				Float ratio2 = fr1.z() / fr2.z();
				f.z()+= ratio2;
				++count[2];
			}
		}

		// compute average
		fresnel[i].x() = count[0] == 0 ? 1.0 : dr::minimum(f.x() / (Float)count[0],1.0f);
		fresnel[i].y() = count[1] == 0 ? 1.0 :dr::minimum( f.y() / (Float)count[1],1.0f);
		fresnel[i].z() = count[2] == 0 ? 1.0 : dr::minimum( f.z() / (Float)count[2],1.0f);
	}
	// copy last value
	fresnel[res - 1] = fresnel[res - 2];
	microfacet<Float>::set_fresnel(fresnel::spline<Float>(fresnel));
#ifndef NVERBOSE
	Throw("djb_verbose: Fresnel function ready");
//	DJB_LOG("djb_verbose: Fresnel function ready\n");
#endif
}
template <typename Float>
void tabular_anisotropic<Float>::compute_fresnel(const brdf<Float>& brdf, uint32_t res)
{
	std::vector<Vector3f> fresnel(res);
	uint32_t cnt = res - 1;

	// compute average ratio between input and Torrance Sparrow equation
	for (uint32_t i = 0; i < cnt; ++i) {
		const Float phi_d = dr::Pi<Float> * 0.5; // this can NOT be tweaked
		const Float phi_h = 0.0;        // this can be tweaked (no impact on MERL)
		Float tmp = (Float)i / (Float)cnt;
		Float theta_d = tmp * dr::Pi<Float> * 0.5; // linear parameterization in theta_d
		Vector3f f = Vector3f(0); // Fresnel value at R, G, B wavelengths
		uint32_t count[3] = {0, 0, 0};
		Float theta_h = 0.0;

		for (uint32_t j = 0; theta_h < dr::Pi<Float> * 0.5 - theta_d; ++j) {
			Float tmp1 = (Float)j / (Float)cnt;
			theta_h = tmp1 * tmp1 * dr::Pi<Float> * 0.5; // in [0, pi/2)

			if (theta_h > dr::Pi<Float> * 0.5) continue;

			Vector3f dir_h = Vector3f(theta_h, phi_h);
			Vector3f dir_d = Vector3f(theta_d, phi_d);
			Vector3f dir_i, dir_o;
			hd_to_io(dir_h, dir_d, &dir_i, &dir_o);

			dir_i = Vector3f(0, 0, 1); //XXX: hack to reproduce my EGSR fits
			Vector3f fr1 = brdf.eval(dir_i, dir_o);
			Vector3f fr2 = eval(dir_i, dir_o);

			if (fr2.x > /* epsilon */dr::Epsilon<Float>) {
				Float ratio = fr1.x() / fr2.x();
				f.x()+= ratio;
				++count[0];
			}
			if (fr2.y > /* epsilon */dr::Epsilon<Float>) {
				Float ratio = fr1.y() / fr2.y();
				f.y()+= ratio;
				++count[1];
			}
			if (fr2.z > /* epsilon */dr::Epsilon<Float>) {
				Float ratio = fr1.z() / fr2.z();
				f.z()+= ratio;
				++count[2];
			}
		}

		// compute average
		fresnel[i].x() = count[0] == 0 ? 1.0 : min((Float)1.0, f.x() / (Float)count[0]);
		fresnel[i].y() = count[1] == 0 ? 1.0 : min((Float)1.0, f.y() / (Float)count[1]);
		fresnel[i].z() = count[2] == 0 ? 1.0 : min((Float)1.0, f.z() / (Float)count[2]);
	}
	// copy last value
	fresnel[res - 1] = fresnel[res - 2];
	set_fresnel(fresnel::spline(fresnel));
#ifndef NVERBOSE
	Throw("djb_verbose: Fresnel function ready");
	//DJB_LOG("djb_verbose: Fresnel function ready\n");
#endif
}

// -------------------------------------------------------------------------------------------------
// Compute the CDF (stored for debugging)
template <typename Float>
void tabular<Float>::compute_cdf()
{
	uint32_t cnt = (uint32_t)m_p22.size() - 1;
	Float dtheta = dr::Pi<Float> / (Float)cnt;
	Float nint = 0.0;

	m_cdf.resize(0);
	for (uint32_t i = 0; i < cnt; ++i) {
		Float u = i / cnt;
		Float theta_h = u * u * dr::Pi<Float> * 0.5;
		Float cos_theta_h = dr::cos(theta_h);
		Float r_h = dr::tan(theta_h);
		Float p22_r = p22_radial(r_h * r_h);

		nint+= (u * r_h * p22_r) / (cos_theta_h * cos_theta_h);
		m_cdf.push_back(nint * dtheta * /* normalize */(2.0 * dr::Pi<Float>));
	}
	m_cdf.push_back(1);

#ifndef NVERBOSE
	Throw("djb_verbose: Slope CDF ready");
//	DJB_LOG("djb_verbose: Slope CDF ready\n");
#endif
}

// -------------------------------------------------------------------------------------------------
// Compute the Quantile Function
template <typename Float>
void tabular<Float>::compute_qf()
{    MI_IMPORT_CORE_TYPES();
	uint32_t cnt = (uint32_t)m_p22.size() - 1;
	uint32_t res = cnt * 8; // resolution of inversion
	uint32_t j = 0;
	Mask active=true;
    DRJIT_MARK_USED(active);   
	m_qf.resize(0);
	m_qf.push_back(0);
	for (uint32_t i = 1; i < cnt; ++i) {
		Float cdf = (Float)i / (Float)cnt;

		for (uint32_t j=0; j < res; ++j) {
			Float u = j /res;
			Float theta_h = u * dr::Pi<Float> * 0.5;
			Float qf = cdf_radial(dr::tan(theta_h));

			// lerp lookup
			active &= qf >= cdf;
			if (dr::none_or<true>(active)) {
				m_qf.push_back(u);
				break;
			} else if (j == res) {
				m_qf.push_back(1.0);
				break;
			}
		}
	}

	m_qf.push_back(1.0);
#ifndef NVERBOSE
	Throw("djb_verbose: Slope QF ready");
//	DJB_LOG("djb_verbose: Slope QF ready\n");
#endif
}

// -------------------------------------------------------------------------------------------------
// Tabular Anisotropic Fetches
template <typename Float>
Float tabular_anisotropic<Float>::pdf1(Float phi) const
{
	Float u = phi * 0.5 / dr::Pi<Float>;
	return spline::eval<Float>(m_pdf1, spline::uwrap_repeat, u);
}
template <typename Float>
Float tabular_anisotropic<Float>::cdf1(Float phi) const
{
	Float u = phi * 0.5 / dr::Pi<Float>;
	return spline::eval<Float>(m_cdf1, spline::uwrap_repeat, u);
}
template <typename Float>
Float tabular_anisotropic<Float>::qf1(Float u1) const
{
	DJB_ASSERT(u1 >= 0.0 && u1 <= 1.0 && "Invalid Variate");
	return spline::eval<Float>(m_qf1, spline::uwrap_edge, u1) * 2.0 * dr::Pi<Float>;
}
template <typename Float>
Float tabular_anisotropic<Float>::pdf2(Float theta, Float phi) const
{
	DJB_ASSERT(theta >= 0.0 && "Invalid Angle");
	if (theta >= 0.5 * dr::Pi<Float>) return 0.0;
	uint32_t w = m_elevation_res;
	uint32_t h = m_azimuthal_res;
	Float u1 = theta * 2.0 / dr::Pi<Float>; // in [0, 1]
	Float u2 = phi * 0.5 / dr::Pi<Float>; // in [0, 1]

	return spline::eval2d<Float>(m_pdf2, w, h,
	                      spline::uwrap_edge, u1,
	                      spline::uwrap_repeat, u2);
}
template <typename Float>
Float tabular_anisotropic<Float>::cdf2(Float theta, Float phi) const
{
	DJB_ASSERT(theta >= 0.0 && "Invalid Angle");
	if (theta >= 0.5 * dr::Pi<Float>) return 1.0;
	uint32_t w = m_elevation_res;
	uint32_t h = m_azimuthal_res;
	Float u1 = theta * 2.0 / dr::Pi<Float>; // in [0, 1]
	Float u2 = phi * 0.5 / dr::Pi<Float>; // in [0, 1]

	return spline::eval2d<Float>(m_cdf2, w, h,
	                      spline::uwrap_edge, u1,
	                      spline::uwrap_repeat, u2);
}
template <typename Float>
Float tabular_anisotropic<Float>::qf2(Float u, Float phi) const
{
	DJB_ASSERT(u >= 0.0 && u <= 1.0 && "Invalid Variate");
	uint32_t w = m_elevation_res;
	uint32_t h = m_azimuthal_res;
	Float u1 = phi / (2.0 * dr::Pi<Float>);

	return spline::eval2d<Float>(m_qf2, w, h,
	                      spline::uwrap_edge, u,
	                      spline::uwrap_repeat, u1) * 0.5 * dr::Pi<Float>;
}

// -------------------------------------------------------------------------------------------------
// Importance sample an anisotropic microfacet NDF
template <typename Float>
void
tabular_anisotropic<Float>::sample_vp22_std_nmap(
	Float u1, Float u2, const Vector3f& k,
	Float *xslope, Float *yslope
) const {
	Float phi = qf1(u1);
	Float theta = qf2(u2, phi);
	Float tan_theta = tan(theta);

	(*xslope) = -tan_theta * cos(phi);
	(*yslope) = -tan_theta * sin(phi);
}

// -------------------------------------------------------------------------------------------------
/**
 * Compute the Marginal PDF (stored for debug)
 *
 * We compute the marginal slope PDF in (hemi-)spherical space.
 * The elevation parameter is marginalized.
 * Spherical space is more conveniant for numerical work than slope space:
 * it is free from singularities if we marginalize the elevation parameter.
 */
template <typename Float>
void tabular_anisotropic<Float>::compute_pdf1()
{
	uint32_t ntheta = 256;
	uint32_t nphi   = m_azimuthal_res;
	Float dtheta = 0.5 * dr::Pi<Float> / (Float)ntheta;

	m_pdf1.resize(0);
	for (uint32_t i = 0; i < nphi; ++i) {
		Float u = (Float)i / (Float)nphi; // in (0,1)
		Float phi = u * 2.0 * dr::Pi<Float>; // in [0,2pi)
		Float nint = 0.0;

		for (uint32_t j = 0; j < ntheta; ++j) {
			Float u = (Float)j / (Float)ntheta; // in (0,1)
			Float theta = u * 0.5 * dr::Pi<Float>; // in [0,pi/2)
			Float p22 = this->p22_std_theta_phi(theta, phi);
			Float cos_theta = cos(theta);

			nint+= (p22 * tan(theta)) / (cos_theta * cos_theta);
		}
		nint*= dtheta;
		m_pdf1.push_back(nint);
	}
	normalize_pdf1();
#ifndef NVERBOSE
	Throw("djb_verbose: PDF_1 ready");
//	DJB_LOG("djb_verbose: PDF_1 ready\n");
#endif
}

// -------------------------------------------------------------------------------------------------
// Compute CDF 1 (stored for debug)
template <typename Float>
void tabular_anisotropic<Float>::compute_cdf1()
{
	uint32_t cnt = (uint32_t)m_azimuthal_res - 1;
	Float dphi = 2.0 * dr::Pi<Float> / (Float)cnt;

	m_cdf1.resize(0);
	Float nint = 0.0;

	m_cdf1.push_back(0);
	for (uint32_t i = 1; i < cnt; ++i) {
		Float u = (Float)i / (Float)cnt;
		Float phi = u * 2.0 * dr::Pi<Float>;
		Float pdf = pdf1(phi);

		nint+= pdf;
		m_cdf1.push_back(nint * dphi);
	}

	m_cdf1.push_back(1.0);
#ifndef NVERBOSE
	Throw("djb_verbose: CDF_1 ready");
//	DJB_LOG("djb_verbose: CDF_1 ready\n");
#endif
}

// -------------------------------------------------------------------------------------------------
// Compute the Quantile Function
template <typename Float>
void tabular_anisotropic<Float>::compute_qf1()
{
	uint32_t cnt = (uint32_t)m_cdf1.size() - 1;
	uint32_t res = cnt * 8; // resolution for inversion
	uint32_t j = 0;

	m_qf1.resize(0);
	m_qf1.push_back(0);
	for (uint32_t i = 1; i < cnt; ++i) {
		Float cdf = (Float)i / (Float)cnt;

		for (; j < res; ++j) {
			Float u = (Float)j / (Float)res;
			Float phi = u * 2.0 * dr::Pi<Float>;
			Float qf = cdf1(phi);

			// lerp lookup
			if (qf >= cdf) {
				m_qf1.push_back(u);
				break;
			} else if (j == res) {
				m_qf1.push_back(1.0);
				break;
			}
		}
	}

	m_qf1.push_back(1.0);
#ifndef NVERBOSE
	Throw("djb_verbose: QF_1 ready");
//	DJB_LOG("djb_verbose: QF_1 ready\n");
#endif
}

// -------------------------------------------------------------------------------------------------
/**
 * Compute the Conditional PDF (stored for debug)
 *
 * We compute the conditional PDF in (hemi-)spherical space.
 * This is the PDF of elevation angle conditioned on azimuthal angle.
 */
template <typename Float>
void tabular_anisotropic<Float>::compute_pdf2()
{
	uint32_t ntheta = m_elevation_res - 1;
	uint32_t nphi = m_azimuthal_res;

	m_pdf2.resize(0);
	for (uint32_t i = 0; i < nphi; ++i) {
		Float u = (Float)i / (Float)nphi; // in [0,1)
		Float phi = u * 2.0 * dr::Pi<Float>; // in [0,2pi)

		for (uint32_t j = 0; j < ntheta; ++j) {
			Float u = (Float)j / (Float)ntheta; // in [0,1)
			Float theta = u * 0.5 * dr::Pi<Float>; // in [0,pi/2)
			Float p22 = this->p22_std_theta_phi(theta, phi);
			Float p1 = this->pdf1(phi);

			m_pdf2.push_back(p22 / p1);
		}
		m_pdf2.push_back(0);
	}

	normalize_pdf2();
#ifndef NVERBOSE
	Throw("djb_verbose: PDF_2 ready");
//	DJB_LOG("djb_verbose: PDF_2 ready\n");
#endif
}

// -------------------------------------------------------------------------------------------------
// Compute CDF 2 (stored for debug)
template <typename Float>
void tabular_anisotropic<Float>::compute_cdf2()
{
	uint32_t ntheta = m_elevation_res - 1;
	uint32_t nphi = m_azimuthal_res;
	Float dtheta = 0.5 * dr::Pi<Float> / (Float)ntheta;

	m_cdf2.resize(0);
	for (uint32_t i = 0; i < nphi; ++i) {
		Float u = (Float)i / (Float)nphi; // in [0,1)
		Float phi = u * 2.0 * dr::Pi<Float>; // in [0,2pi)
		Float nint = 0.0;

		for (uint32_t j = 0; j < ntheta; ++j) {
			Float u = (Float)j / (Float)ntheta; // in [0,1)
			Float theta = u * 0.5 * dr::Pi<Float>; // in [0,pi/2)
			Float pdf2 = this->pdf2(theta, phi);
			Float cos_theta = cos(theta);

			nint+= (pdf2 * tan(theta)) / (cos_theta * cos_theta);
			m_cdf2.push_back(nint * dtheta);
		}
		m_cdf2.push_back(1);
	}

#ifndef NVERBOSE
	Throw("djb_verbose: CDF_2 ready");
//	DJB_LOG("djb_verbose: CDF_2 ready\n");
#endif
}

// -------------------------------------------------------------------------------------------------
// Compute the Quantile Function
template <typename Float>
void tabular_anisotropic<Float>::compute_qf2()
{
	uint32_t ntheta = m_elevation_res - 1;
	uint32_t nphi = m_azimuthal_res;
	uint32_t res = ntheta * 8; // resolution of inversion

	m_qf2.resize(0);
	for (uint32_t k = 0; k < nphi; ++k) {
		Float u = (Float)k / (Float)nphi;
		Float phi = u * 2.0 * dr::Pi<Float>;
		uint32_t j = 0;

		m_qf2.push_back(0);
		for (uint32_t i = 1; i < ntheta; ++i) {
			Float cdf = (Float)i / (Float)ntheta;

			for (; j < res; ++j) {
				Float u = (Float)j / (Float)res;
				Float theta = u * 0.5 * dr::Pi<Float>; // in [0, pi/2)
				Float qf = cdf2(theta, phi);

				// lerp lookup
				if (qf >= cdf) {
					m_qf2.push_back(u);
					break;
				} else if (j == res) {
					m_qf2.push_back(1.0);
					break;
				}
			}
		}
		m_qf2.push_back(1.0);
	}

#ifndef NVERBOSE
	Throw("djb_verbose: QF_2 ready");
//	DJB_LOG("djb_verbose: QF_2 ready\n");
#endif
}

// -------------------------------------------------------------------------------------------------
// Normalize pdf1
template <typename Float>
void tabular_anisotropic<Float>::normalize_pdf1()
{
	uint32_t cnt = 512;
	Float dphi = 2.0 * dr::Pi<Float> / (Float)cnt;
	Float nint = 0.0;
	Float k;

	for (uint32_t i = 0; i < cnt; ++i) {
		Float tmp = (Float)i / (Float)cnt; // in (0,1)
		Float phi = tmp * 2.0 * dr::Pi<Float>; // in [0,2pi)
		Float p1 = this->pdf1(phi);

		nint+= p1;
	}
	nint*= dphi;
	k = 1.0 / nint;
	for (uint32_t i = 0; i < (uint32_t)m_pdf1.size(); ++i)
		m_pdf1[i]*= k;
#ifndef NVERBOSE 
    Throw("djb_verbose: PDF P1 norm. cst: %f", k);
//	DJB_LOG("djb_verbose: PDF P1 norm. cst: %.9f\n", k);
#endif
}

// -------------------------------------------------------------------------------------------------
// Normalize pdf2
template <typename Float>
void tabular_anisotropic<Float>::normalize_pdf2()
{
	uint32_t ntheta = 256;
	uint32_t nphi = m_azimuthal_res;
	Float dtheta = 0.5 * dr::Pi<Float> / (Float)ntheta;
	std::vector<Float> k(nphi, 0);

	k.resize(0);
	for (uint32_t j = 0; j < nphi; ++j) {
		Float tmp = (Float)j / (Float)nphi;
		Float phi = tmp * 2.0 * dr::Pi<Float>; // in [0,2pi)
		Float nint = 0.0;

		for (uint32_t i = 0; i < ntheta; ++i) {
			Float tmp = (Float)i / (Float)ntheta; // in [0,1)
			Float theta = tmp * 0.5 * dr::Pi<Float>; // in [0,pi/2)
			Float pdf2 = this->pdf2(theta, phi);
			Float cos_theta = cos(theta);

			nint+= (pdf2 * tan(theta)) / (cos_theta * cos_theta);
		}
		nint*= dtheta;
		k.push_back(1.0 / nint);
	}

	for (uint32_t j = 0; j < nphi; ++j) {
		for (uint32_t i = 0; i < m_elevation_res; ++i)
			m_pdf2[i + m_elevation_res * j]*= k[j];
#ifndef NVERBOSE 
		// DJB_LOG("djb_verbose: PDF P2 norm. cst: %.9f\n", k[j]);
#endif
	}
}

// -------------------------------------------------------------------------------------------------
// Accessors
template <typename Float>
const std::vector<Float>&
tabular_anisotropic<Float>::get_p22v(uint32_t *elev_cnt, uint32_t *azim_cnt) const
{
	if (elev_cnt) *elev_cnt = m_elevation_res;
	if (azim_cnt) *azim_cnt = m_azimuthal_res;

	return m_p22;
}
template <typename Float>
const std::vector<Float>&
tabular_anisotropic<Float>::get_sigmav(uint32_t *elev_cnt, uint32_t *azim_cnt) const
{
	if (elev_cnt) *elev_cnt = m_elevation_res;
	if (azim_cnt) *azim_cnt = m_azimuthal_res;

	return m_sigma;
}

// -------------------------------------------------------------------------------------------------
/**
 * Direct Conversion to Parametric BRDFs
 *
 * The tabulated BRDFs may be converted directly to either a 
 * Gaussian or a GGX BRDF. For this, moments are computed from the 
 * tabulated microfacet slope distribution in order to retrieve
 * the scale and correlation parameters.
 */
/*
template <typename Float>
typename microfacet<Float>::params tabular<Float>::fit_beckmann_parameters(const tabular<Float>& tab)
{
	const int ntheta = 128;
	Float dtheta = dr::Pi<Float> / (Float)ntheta;
	Float nint = 0.0;
	Float alpha;

	for (int i = 0; i < ntheta; ++i) {
		Float u = (Float)i / (Float)ntheta; // in [0,1)
		Float theta_h = u * u * M_PI * 0.5;
		Float cos_theta_h = cos(theta_h);
		Float r_h = tan(theta_h);
		Float r_h_sqr = r_h * r_h;
		Float p22_r = tab.p22_radial(r_h_sqr);

		nint+= (u * r_h_sqr * r_h * p22_r) / (cos_theta_h * cos_theta_h);
	}/*
	////nint*= dtheta * /* int_0^2pi cos^2 phi dphi *///M_PI;
	//alpha = sqrt(2.0 * nint);

//#ifndef NVERBOSE
//	DJB_LOG("djb_verbose: Beckmann_alpha = %.9f\n", alpha);
//#endif

//	return microfacet::params::isotropic(alpha);
//}
/*
template <typename Float>
typename microfacet::params tabular<Float>::fit_ggx_parameters(const tabular<Float>& tab)
{
	const int ntheta = 128;
	Float dtheta = M_PI / (Float)ntheta;
	Float nint = 0.0;
	Float alpha;

	for (int i = 0; i < ntheta; ++i) {
		Float u = (Float)i / (Float)ntheta; // in [0,1)
		Float theta_h = u * u * M_PI * 0.5;
		Float cos_theta_h = cos(theta_h);
		Float r_h = tan(theta_h);
		Float r_h_sqr = r_h * r_h;
		Float p22_r = tab.p22_radial(r_h_sqr);

		nint+= (u * r_h_sqr * p22_r) / (cos_theta_h * cos_theta_h);
	}
	nint*= dtheta * /* int_0^2pi fabs(cos(phi)) dphi *///4.0;/*
	/*
	alpha = nint;
#ifndef NVERBOSE
	DJB_LOG("djb_verbose: GGX_alpha = %.9f\n", alpha);
#endif

	return microfacet::params::isotropic(alpha);
}*/
/*
template <typename Float>
microfacet::params
tabular_anisotropic::fit_beckmann_parameters(const tabular_anisotropic& tab)
{
	const int ntheta = 128;
	const int nphi   = 512;
	Float dtheta = sqrt(M_PI * 0.5) / (Float)ntheta;
	Float dphi = 2.0 * M_PI / (Float)nphi;
	Float nint[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
	Float mux, muy, alphax, alphay, rho;

	for (int j = 0; j < nphi; ++j) {
		Float tmp = (Float)j / (Float)nphi; // in [0,1)
		Float phi = tmp * 2.0 * M_PI; // in [0,2pi)
		Float cos_phi = cos(phi);
		Float sin_phi = sin(phi);
		Float cos_phi_sqr = cos_phi * cos_phi;
		Float sin_phi_sqr = sin_phi * sin_phi;

		for (int i = 0; i < ntheta; ++i) {
			Float tmp1 = (Float)i / (Float)ntheta; // in [0,1)
			Float theta = tmp1 * sqrt(M_PI * 0.5); // in [0, sqrt(pi/2))
			Float theta_sqr = theta * theta; // in [0, pi/2)
			Float p22 = tab.p22_std_theta_phi(theta_sqr, phi);
			Float tan_theta = tan(theta_sqr);
			Float cos_theta = cos(theta_sqr);
			Float tan_theta_sqr = tan_theta * tan_theta;
			Float cos_theta_sqr = cos_theta * cos_theta;
			Float tmp2 = theta * p22 * tan_theta / cos_theta_sqr;
			Float e1 = -tan_theta * cos_phi; // mux
			Float e2 = -tan_theta * sin_phi; // muy
			Float e3 = tan_theta_sqr * cos_phi_sqr; // var1
			Float e4 = tan_theta_sqr * sin_phi_sqr; // var2
			Float e5 = tan_theta_sqr * cos_phi * sin_phi; // cov

			nint[0]+= tmp2 * e1;
			nint[1]+= tmp2 * e2;
			nint[2]+= tmp2 * e3;
			nint[3]+= tmp2 * e4;
			nint[4]+= tmp2 * e5;
		}
	}
	for (int i = 0; i < 5; i++)
		nint[i]*= 2.0 * dtheta * dphi;
	mux = nint[0];
	muy = nint[1];
	alphax = sqrt((Float)2.0 * (nint[2] - mux * mux));
	alphay = sqrt((Float)2.0 * (nint[3] - muy * muy));
	rho = 2.0 * (nint[4] - mux * muy) / (alphax * alphay);

#ifndef NVERBOSE
	DJB_LOG(
		"djb_verbose: Beckmann params = %.9f %.9f %.9f %.9f %.9f\n", 
		alphax,
		alphay,
		rho,
		mux,
		muy
	);
#endif

	return microfacet::params::pdfparams(alphax, alphay, rho, mux, muy);
}
template <typename Float>
microfacet::params
tabular_anisotropic::fit_ggx_parameters(const tabular_anisotropic& tab)
{
	const int ntheta = 128;
	const int nphi   = 512;
	Float dtheta = sqrt(M_PI * 0.5) / (Float)ntheta;
	Float dphi = 2.0 * M_PI / (Float)nphi;
	Float nint[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
	Float mux, muy, alphax, alphay, rho;

	for (int j = 0; j < nphi; ++j) {
		Float tmp = (Float)j / (Float)nphi; // in [0,1)
		Float phi = tmp * 2.0 * M_PI; // in [0,2pi)
		Float cos_phi = cos(phi);
		Float sin_phi = sin(phi);

		for (int i = 0; i < ntheta; ++i) {
			Float tmp1 = (Float)i / (Float)ntheta; // in [0,1)
			Float theta = tmp1 * sqrt(M_PI * 0.5); // in [0, sqrt(pi/2))
			Float theta_sqr = theta * theta; // in [0, pi/2)
			Float p22 = tab.p22_std_theta_phi(theta_sqr, phi);
			Float tan_theta = tan(theta_sqr);
			Float cos_theta = cos(theta_sqr);
			Float cos_theta_sqr = cos_theta * cos_theta;
			Float tmp2 = theta * p22 * tan_theta / cos_theta_sqr;
			Float e1 = -tan_theta * cos_phi; // mux
			Float e2 = -tan_theta * sin_phi; // muy
			Float e3 = fabs(e1); // var1
			Float e4 = fabs(e2); // var2
			Float e5 = 0; // TODO

			nint[0]+= tmp2 * e1;
			nint[1]+= tmp2 * e2;
			nint[2]+= tmp2 * e3;
			nint[3]+= tmp2 * e4;
			nint[4]+= tmp2 * e5;
		}
	}
	for (int i = 0; i < 5; i++)
		nint[i]*= 2.0 * dtheta * dphi;
	mux = nint[0];
	muy = nint[1];
	alphax = sqrt(nint[2] * nint[2] - mux * mux);
	alphay = sqrt(nint[3] * nint[3] - muy * muy);
	rho = 0.0; // TODO

#ifndef NVERBOSE
	DJB_LOG(
		"djb_verbose: GGX params = %.9f %.9f %.9f %.9f %.9f\n", 
		alphax,
		alphay,
		rho,
		mux,
		muy
	);
#endif

	return microfacet::params::pdfparams(alphax, alphay, rho, mux, muy);
}*/

// *************************************************************************************************
// Shifted Gamma Distribution API implementation
/*
template <typename Float>
const sgd<Float> ::data sgd<Float> ::s_data[] = {
	{ "alum-bronze", "alum-bronze", { 0.0478786, 0.0313514, 0.0200638 }, { 0.0364976, 0.664975, 0.268836 }, { 0.014832, 0.0300126, 0.0490339 }, { 0.459076, 0.450056, 0.529272 }, { 6.05524, 0.235756, 0.580647 }, { 5.05524, 0.182842, 0.476088 }, { 46.3841, 24.5961, 14.8261 }, { 2.60672, 2.97371, 2.7827 }, { 1.12717e-07, 1.06401e-07, 5.27952e-08 }, { 47.783, 36.2767, 31.6066 }, { 0.205635, 0.066289, -0.0661091 }, { 0.100735, 0.0878706, 0.0861907 } },
	{ "alumina-oxide", "alumina-oxide", { 0.316358, 0.292248, 0.25416 }, { 0.00863128, 0.00676832, 0.0103309 }, { 0.000159222, 0.000139421, 0.000117714 }, { 0.377727, 0.318496, 0.402598 }, { 0.0300766, 1.70375, 1.96622 }, { -0.713784, 0.70375, 1.16019 }, { 4381.96, 5413.74, 5710.42 }, { 3.31076, 4.93831, 2.84538 }, { 6.72897e-08, 1.15769e-07, 6.32199e-08 }, { 354.275, 367.448, 414.581 }, { 0.52701, 0.531166, 0.53301 }, { 0.213276, 0.147418, 0.27746 } },
	{ "aluminium", "aluminium", { 0.0305166, 0.0358788, 0.0363463 }, { 0.0999739, 0.131797, 0.0830361 }, { 0.0012241, 0.000926487, 0.000991844 }, { 0.537669, 0.474562, 0.435936 }, { 0.977854, 0.503108, 1.77905 }, { -0.0221457, -0.0995445, 0.77905 }, { 449.321, 658.044, 653.86 }, { 8.2832e-07, 9.94692e-08, 6.11887e-08 }, { 3.54592e-07, 16.0175, 15.88 }, { 23.8656, 10.6911, 9.69801 }, { -0.510356, 0.570179, 0.566156 }, { 0.303567, 0.232628, 0.441578 } },
	{ "aventurnine", "aventurnine", { 0.0548217, 0.0621179, 0.0537826 }, { 0.0348169, 0.0872381, 0.111961 }, { 0.000328039, 0.000856166, 0.00145342 }, { 0.387167, 0.504525, 0.652122 }, { 0.252033, 0.133897, 0.087172 }, { 0.130593, 0.0930416, 0.0567429 }, { 2104.51, 676.157, 303.59 }, { 3.12126, 2.50965e-07, 2.45778e-05 }, { 1.03849e-07, 8.53824e-07, 3.20722e-07 }, { 251.265, 24.2886, 29.0236 }, { 0.510125, -0.41764, -0.245097 }, { 0.0359759, 0.0297523, 0.0285881 } },
	{ "beige-fabric", "fabric-beige", { 0.20926, 0.160666, 0.145337 }, { 0.121663, 0.0501577, 0.00177279 }, { 0.39455, 0.15975, 0.110706 }, { 0.474725, 0.0144728, 1.70871e-12 }, { 0.0559459, 0.222268, 8.4764 }, { -0.318718, -0.023826, 7.4764 }, { 3.8249, 7.32453, 10.0904 }, { 2.26283, 2.97144, 3.55311 }, { 0.0375346, 0.073481, 0.0740222 }, { 7.52635, 9.05672, 10.6185 }, { 0.217453, 0.407084, 0.450203 }, { 0.00217528, 0.00195262, 0.00171008 } },
	{ "black-fabric", "fabric-black", { 0.0189017, 0.0112353, 0.0110067 }, { 2.20654e-16, 6.76197e-15, 1.57011e-13 }, { 0.132262, 0.128044, 0.127838 }, { 0.189024, 0.18842, 0.188426 }, { 1, 1, 1 }, { 0, 0, 0 }, { 8.1593, 8.38075, 8.39184 }, { 3.83017, 3.89536, 3.89874 }, { 0.00415117, 0.00368324, 0.00365826 }, { 12.9974, 13.2597, 13.2737 }, { 0.207997, 0.205597, 0.205424 }, { 0.000363154, 0.000272253, 0.000274773 } },
	{ "black-obsidian", "black-obsidian", { 0.00130399, 0.0011376, 0.00107233 }, { 0.133029, 0.125362, 0.126188 }, { 0.000153649, 0.000148939, 0.000179285 }, { 0.186234, 0.227495, 0.25745 }, { 2.42486e-12, 0.0174133, 0.091766 }, { -0.0800755, -0.048671, 0.0406445 }, { 5668.57, 5617.79, 4522.84 }, { 13.7614, 8.59526, 6.44667 }, { 1e+38, 1e+38, 1e+38 }, { 117.224, 120.912, 113.366 }, { 1.19829, 1.19885, 1.19248 }, { 0.0841374, 0.0943085, 0.123514 } },
	{ "black-oxidized-steel", "black-oxidized-steel", { 0.0149963, 0.0120489, 0.0102471 }, { 0.373438, 0.344382, 0.329202 }, { 0.187621, 0.195704, 0.200503 }, { 0.661367, 0.706913, 0.772267 }, { 0.0794166, 0.086518, 0.080815 }, { 0.0470402, 0.0517633, 0.0455037 }, { 5.1496, 4.91636, 4.69009 }, { 4.0681, 3.95489, 3.71052 }, { 1.07364e-07, 1.05341e-07, 1.16556e-07 }, { 20.2383, 20.1786, 20.2553 }, { -0.479617, -0.4885, -0.478388 }, { 0.00129576, 0.0011378, 0.000986163 } },
	{ "black-phenolic", "black-phenolic", { 0.00204717, 0.00196935, 0.00182908 }, { 0.177761, 0.293146, 0.230592 }, { 0.00670804, 0.00652009, 0.00656043 }, { 0.706648, 0.677776, 0.673986 }, { 0.16777, 0.12335, 0.166663 }, { 0.111447, 0.0927321, 0.125663 }, { 65.4189, 70.8936, 70.9951 }, { 1.06318, 1.15283, 1.16529 }, { 1.24286e-07, 3.00039e-08, 9.77334e-08 }, { 74.0711, 75.1165, 73.792 }, { 0.338204, 0.319306, 0.33434 }, { 0.0307129, 0.0531183, 0.0454238 } },
	{ "black-soft-plastic", "black-plastic-soft", { 0.00820133, 0.00777718, 0.00764537 }, { 0.110657, 0.0980322, 0.100579 }, { 0.0926904, 0.0935964, 0.0949975 }, { 0.14163, 0.148703, 0.143694 }, { 0.150251, 0.169418, 0.170457 }, { 0.100065, 0.113089, 0.114468 }, { 11.2419, 11.113, 10.993 }, { 4.3545, 4.3655, 4.31586 }, { 0.00464641, 0.00384785, 0.0046145 }, { 14.6751, 14.8089, 14.5436 }, { 0.275651, 0.262317, 0.271284 }, { 0.000527402, 0.000550255, 0.000527718 } },
	{ "blue-acrylic", "acrylic-blue", { 0.0134885, 0.0373766, 0.10539 }, { 0.0864901, 0.0228191, 0.204042 }, { 0.000174482, 0.000269795, 0.0015211 }, { 0.373948, 0.362425, 0.563636 }, { 0.0185562, 0.399982, 0.0525861 }, { -0.0209713, 0.241543, 0.0169474 }, { 4021.24, 2646.36, 346.898 }, { 3.38722, 3.62885, 1.83684e-06 }, { 9.64334e-08, 9.96105e-08, 3.61787e-07 }, { 338.073, 272.828, 23.5039 }, { 0.526039, 0.515404, -0.526935 }, { 0.0612235, 0.0789826, 0.0461093 } },
	{ "blue-fabric", "fabric-blue", { 0.0267828, 0.0281546, 0.066668 }, { 0.0825614, 0.0853369, 0.0495164 }, { 0.248706, 0.249248, 0.18736 }, { 9.23066e-13, 1.66486e-12, 2.27218e-12 }, { 0.201626, 0.213723, 0.56548 }, { 0.225891, 0.226267, 0.638493 }, { 5.15615, 5.14773, 6.43713 }, { 2.25846, 2.25536, 2.68382 }, { 0.128037, 0.128363, 0.0944915 }, { 6.95531, 6.94633, 8.17665 }, { 0.407528, 0.407534, 0.411378 }, { 0.000323043, 0.00032064, 0.000653631 } },
	{ "blue-metallic-paint2", "ch-ball-blue-metallic", { 0.010143, 0.0157349, 0.0262717 }, { 0.0795798, 0.0234493, 0.0492337 }, { 0.00149045, 0.00110477, 0.00141008 }, { 0.624615, 0.598721, 0.67116 }, { 9.36434e-14, 3.61858e-15, 1.15633e-14 }, { -0.210234, -1, -1 }, { 314.024, 441.812, 299.726 }, { 1.20935e-05, 7.51792e-06, 3.86474e-05 }, { 3.38901e-07, 2.94502e-07, 3.15718e-07 }, { 27.0491, 28.576, 30.6214 }, { -0.326593, -0.274443, -0.187842 }, { 0.0908879, 0.163236, 0.286541 } },
	{ "blue-metallic-paint", "metallic-blue", { 0.00390446, 0.00337319, 0.00848198 }, { 0.0706771, 0.0415082, 0.104423 }, { 0.155564, 0.139, 0.15088 }, { 1.01719, 1.02602, 1.16153 }, { 0.149347, 0.153181, 1.87241e-14 }, { -0.487331, -0.76557, -1 }, { 4.4222, 4.59265, 3.93929 }, { 2.54345, 2.33884, 2.24405 }, { 6.04906e-08, 5.81858e-08, 1.2419e-07 }, { 23.9533, 25.0641, 24.6856 }, { -0.34053, -0.294595, -0.258117 }, { 0.00202544, 0.00246268, 0.0059725 } },
	{ "blue-rubber", "blue-rubber", { 0.0371302, 0.0732915, 0.146637 }, { 0.384232, 0.412357, 0.612608 }, { 0.218197, 0.2668, 0.478375 }, { 0.815054, 1.00146, 1.24995 }, { 0.0631713, 0.0622636, 0.0399196 }, { 0.0478254, 0.0422186, 0.007015 }, { 4.41586, 3.76795, 3.46276 }, { 3.77807, 3.82679, 3.33186 }, { 1.2941e-07, 1.07194e-07, 0.00045665 }, { 19.8046, 19.3115, 11.4364 }, { -0.499472, -0.557706, -0.172177 }, { 0.000630461, 0.000733287, 0.00171092 } },
	{ "brass", "brass", { 0.0301974, 0.0223812, 0.0139381 }, { 0.0557826, 0.0376687, 0.0775998 }, { 0.0002028, 0.000258468, 0.00096108 }, { 0.362322, 0.401593, 0.776606 }, { 0.639886, 0.12354, 0.0197853 }, { -0.360114, -0.87646, -0.0919344 }, { 3517.61, 2612.49, 331.815 }, { 3.64061, 2.87206, 0.529487 }, { 1.01146e-07, 9.83073e-08, 5.48819e-08 }, { 312.802, 283.431, 183.091 }, { 0.522711, 0.516719, 0.474834 }, { 0.440765, 0.200948, 0.15549 } },
	{ "cherry-235", "cherry-235", { 0.0497502, 0.0211902, 0.0120688 }, { 0.166001, 0.202786, 0.165189 }, { 0.0182605, 0.0277997, 0.0255721 }, { 0.0358348, 0.163231, 0.129135 }, { 0.0713408, 0.0571719, 0.0791809 }, { 0.0200814, 0.00887306, 0.0251675 }, { 54.7448, 33.7294, 37.36 }, { 6.11314, 6.23697, 6.16351 }, { 30.3886, 0.00191869, 0.0495069 }, { 18.8114, 25.8454, 23.3753 }, { 0.816378, 0.387479, 0.522125 }, { 0.00235908, 0.00199931, 0.00155902 } },
	{ "chrome", "chrome", { 0.00697189, 0.00655268, 0.0101854 }, { 0.0930656, 0.041946, 0.104558 }, { 0.00013, 0.0002, 6e-5 }, { 0.3, 0.36, 0.26 }, { 0.256314, 0.819565, 3.22085e-13 }, { -0.743686, -0.180435, -1 }, { 2200, 2150, 4100 }, { 3.8545, 4.44817, 5.40959 }, { 5.30781e-08, 1.04045e-07, 1.8e+10 }, { 354.965, 349.356, 457 }, { 0.526796, 0.528469, 1.00293 }, { 0.802138, 1.29121, 1.06148 } },
	{ "chrome-steel", "chrome-steel", { 0.0206718, 0.0240818, 0.024351 }, { 0.129782, 0.109032, 0.0524555 }, { 5.51292e-05, 3.13288e-05, 4.51944e-05 }, { 0.207979, 0.152758, 0.325431 }, { 1.18818e-12, 2.06813e-11, 0.580895 }, { -0.316807, -0.265326, -0.419105 }, { 15466.8, 28628.9, 16531.2 }, { 12.8988, 68.7898, 4.68237 }, { 1e+10, 1e+10, 44.7025 }, { 510, 665, 618.155 }, { 1.20035, 1.2003, 0.579562 }, { 0.552277, 0.448267, 1.12869 } },
	{ "colonial-maple-223", "colonial-maple-223", { 0.100723, 0.0356306, 0.0162408 }, { 0.059097, 0.0661341, 0.11024 }, { 0.0197628, 0.0279336, 0.0621265 }, { 0.0311867, 0.112022, 0.344348 }, { 0.0576683, 0.0617498, 0.0479061 }, { -0.0503364, -0.0382196, -0.0382636 }, { 50.7952, 34.6835, 14.201 }, { 6.03342, 6.01053, 4.56588 }, { 18.9034, 0.087246, 1.03757e-07 }, { 18.3749, 21.6159, 26.9347 }, { 0.801289, 0.544191, -0.139944 }, { 0.00229125, 0.00193164, 0.0013723 } },
	{ "color-changing-paint1", "", { 0.00513496, 0.00500415, 0.00296872 }, { 1.53167, 0.430731, 1.48308 }, { 0.00320129, 0.0023053, 0.0329464 }, { 0.167301, 0.100003, 1.04868 }, { 0.0208525, 0.109366, 0.0553211 }, { 0.0101418, 0.0801275, 0.0357004 }, { 279.058, 407.716, 9.71061 }, { 6.73768, 7.19333, 1.04729 }, { 2.66398e+08, 8.98323e+10, 4.16363e-08 }, { 36.7257, 27.8164, 45.5388 }, { 1.01208, 1.19549, 0.131148 }, { 0.164144, 0.0515527, 0.116602 } },
	{ "color-changing-paint2", "", { 0.00463528, 0.00544054, 0.0070818 }, { 1.35172, 1.47838, 1.29831 }, { 0.0279961, 0.0267135, 0.0257468 }, { 0.720154, 0.717648, 0.694662 }, { 0.019073, 0.00825302, 0.0301024 }, { -0.0181694, -0.0198592, -0.000522292 }, { 18.6722, 19.4709, 20.7514 }, { 1.56832, 1.554, 1.60402 }, { 1.19918e-07, 9.95358e-08, 9.69838e-08 }, { 41.295, 42.1716, 42.3986 }, { 0.124794, 0.129237, 0.132024 }, { 0.0489296, 0.0499662, 0.054555 } },
	{ "color-changing-paint3", "", { 0.00305737, 0.00257341, 0.00263616 }, { 0.880793, 0.691268, 0.707821 }, { 0.014742, 0.0135513, 0.0108894 }, { 0.537248, 0.572188, 0.457665 }, { 0.055479, 0.0666783, 0.0585094 }, { 0.0315282, 0.0457554, 0.0365558 }, { 42.1798, 43.3684, 62.0624 }, { 2.06002, 1.82703, 2.54086 }, { 5.0531e-08, 9.99088e-08, 5.726e-08 }, { 50.0948, 52.0344, 54.643 }, { 0.197971, 0.229219, 0.241049 }, { 0.0484699, 0.0342594, 0.0347198 } },
	{ "dark-blue-paint", "dark-blue-paint", { 0.00665057, 0.0139696, 0.0472605 }, { 0.231099, 0.18931, 0.12528 }, { 0.130681, 0.112103, 0.0629285 }, { 0.238562, 0.17, 0.0157371 }, { 0.112486, 0.106018, 0.139316 }, { 0.0662142, 0.0682314, 0.10607 }, { 8.0828, 9.4312, 16.8231 }, { 4.05526, 4.09607, 4.54368 }, { 0.00104743, 0.0035546, 0.10813 }, { 14.347, 13.928, 13.5626 }, { 0.11894, 0.226752, 0.524826 }, { 0.000406243, 0.000291986, 0.000555778 } },
	{ "dark-red-paint", "dark-red-paint", { 0.237125, 0.0365577, 0.0106149 }, { 0.227405, 0.111055, 0.150433 }, { 0.67372, 0.0706716, 0.20328 }, { 1.38127, 1.4345e-13, 0.27723 }, { 1.56101e-11, 0.0779523, 0.136566 }, { -0.148315, 0.0560349, 0.0961999 }, { 3.95965, 15.1861, 5.71355 }, { 2.27942, 4.30359, 3.20838 }, { 0.0225467, 0.113039, 0.00667899 }, { 7.45998, 12.679, 10.9071 }, { 0.119045, 0.521154, 0.166378 }, { 0.00155397, 0.000508377, 0.000218865 } },
	{ "dark-specular-fabric", "", { 0.0197229, 0.00949167, 0.00798414 }, { 0.556218, 0.401495, 0.378651 }, { 0.140344, 0.106541, 0.166715 }, { 0.249059, 0.177611, 0.434167 }, { 0.0351133, 0.0387177, 0.0370533 }, { 0.0243153, 0.0293178, 0.0264913 }, { 7.60492, 9.81673, 6.19307 }, { 3.93869, 4.23097, 4.3775 }, { 0.00122421, 0.00238545, 8.47126e-06 }, { 13.889, 14.5743, 17.2049 }, { 0.114655, 0.210179, -0.227628 }, { 0.00158681, 0.000974676, 0.000638865 } },
	{ "delrin", "delrin", { 0.272703, 0.249805, 0.220642 }, { 0.536593, 0.727886, 0.64011 }, { 0.176535, 0.344018, 0.208011 }, { 0.762088, 0.823603, 0.85976 }, { 0.0398465, 0.0430719, 0.0290162 }, { -0.0121332, -0.0255915, -0.0305128 }, { 5.02017, 3.78669, 4.39454 }, { 3.44325, 3.49289, 3.48886 }, { 9.20654e-08, 0.000215686, 1.06338e-07 }, { 21.2534, 12.5556, 20.455 }, { -0.438862, -0.184248, -0.478699 }, { 0.00493658, 0.00655617, 0.00422698 } },
	{ "fruitwood-241", "fruitwood-241", { 0.0580445, 0.0428667, 0.0259801 }, { 0.203894, 0.233494, 0.263882 }, { 0.00824986, 0.0534794, 0.0472951 }, { 0.160382, 1.07206, 0.768335 }, { 0.00129482, 0.00689891, 0.01274 }, { -0.0211778, -0.0140791, -0.00665974 }, { 110.054, 6.98485, 11.6203 }, { 6.74678, 1.31986, 1.76021 }, { 162.121, 1.18843e-07, 1.0331e-07 }, { 32.6966, 36.752, 34.371 }, { 0.765181, 0.0514459, 0.0106852 }, { 0.00493613, 0.00460352, 0.00373027 } },
	{ "gold-metallic-paint2", "ch-ball-gold-metallic2", { 0.0796008, 0.0538361, 0.0649523 }, { 0.633627, 1.77116, 0.0564028 }, { 0.00376608, 0.00871206, 0.000572055 }, { 0.415684, 0.368424, 0.623038 }, { 0.0343265, 0.00330259, 4.15759e-12 }, { -0.00929705, -0.0219437, -0.0596944 }, { 181.769, 85.527, 795.042 }, { 2.77279, 3.50114, 0.985212 }, { 4.78294e-08, 9.8399e-08, 4.45358e-08 }, { 83.871, 57.2291, 216.537 }, { 0.365256, 0.276734, 0.491264 }, { 0.117145, 0.17457, 0.0754722 } },
	{ "gold-metallic-paint3", "ch-ball-gold-metallic", { 0.0579212, 0.0416649, 0.0271208 }, { 0.0729896, 0.0597695, 0.037684 }, { 0.00146432, 0.00156513, 0.000977438 }, { 0.529437, 0.551234, 0.504486 }, { 1.84643e-14, 5.7212e-15, 8.68546e-13 }, { -1, -1, -0.648897 }, { 382.887, 345.197, 593.725 }, { 4.97078e-07, 1.12888e-06, 4.15366e-07 }, { 3.98586e-07, 3.69533e-07, 16.0596 }, { 22.0196, 22.6462, 11.7126 }, { -0.634103, -0.588056, 0.578355 }, { 0.116233, 0.0897222, 0.0711514 } },
	{ "gold-metallic-paint", "metallic-gold", { 0.0178625, 0.00995704, 0.00335044 }, { 0.17127, 0.120714, 0.115473 }, { 0.127954, 0.127825, 0.109623 }, { 0.781093, 0.795517, 0.661313 }, { 1.00776e-12, 1.22243e-15, 6.18432e-13 }, { -1, -1, -0.334497 }, { 5.90039, 5.83556, 7.14985 }, { 2.7548, 2.7057, 2.94348 }, { 9.46481e-08, 1.06951e-07, 1.10733e-07 }, { 23.8811, 23.9059, 24.2972 }, { -0.303345, -0.293778, -0.267019 }, { 0.0095612, 0.00483452, 0.00145599 } },
	{ "gold-paint", "gold-paint", { 0.147708, 0.0806975, 0.033172 }, { 0.160592, 0.217282, 0.236425 }, { 0.122506, 0.108069, 0.12187 }, { 0.795078, 0.637578, 0.936117 }, { 9.16095e-12, 1.81225e-12, 0.0024589 }, { -0.596835, -0.331147, -0.140729 }, { 5.98176, 7.35539, 5.29722 }, { 2.64832, 3.04253, 2.3013 }, { 9.3111e-08, 8.80143e-08, 9.65288e-08 }, { 24.3593, 24.4037, 25.3623 }, { -0.284195, -0.277297, -0.245352 }, { 0.00313716, 0.00203922, 0.00165683 } },
	{ "gray-plastic", "gray-plastic", { 0.103233, 0.104428, 0.0983734 }, { 0.494656, 0.517207, 0.52772 }, { 0.00758705, 0.00848095, 0.00887135 }, { 0.557908, 0.556548, 0.545887 }, { 0.0428175, 0.0438899, 0.0569098 }, { 0.0208304, 0.0221893, 0.0375763 }, { 75.6274, 68.3098, 66.5689 }, { 1.72111, 1.76597, 1.8406 }, { 1.06783e-07, 5.31845e-08, 6.53296e-08 }, { 65.6538, 63.2018, 61.5805 }, { 0.308907, 0.283768, 0.280211 }, { 0.0512156, 0.050844, 0.0505763 } },
	{ "grease-covered-steel", "", { 0.0196306, 0.0200926, 0.0187026 }, { 0.0433721, 0.0311621, 0.0326401 }, { 0.00019081, 0.000173919, 0.000217638 }, { 0.164569, 0.141125, 0.219421 }, { 7.12672e-13, 1.06789e-14, 1.61131e-13 }, { -1, -1, -1 }, { 4655.4, 5210.63, 3877.45 }, { 17.2827, 29.2523, 8.54352 }, { 1e+38, 1e+38, 1e+38 }, { 103.46, 105.572, 99.1097 }, { 1.20022, 1.20678, 1.19961 }, { 0.261184, 0.20789, 0.210192 } },
	{ "green-acrylic", "acrylic-green", { 0.0176527, 0.0761863, 0.0432331 }, { 0.0517555, 0.15899, 0.0754193 }, { 0.000185443, 5.19959e-05, 7.95188e-05 }, { 0.288191, 0.170979, 0.145492 }, { 0.418137, 0.0486445, 0.170328 }, { 0.330999, 0.0294954, 0.110218 }, { 4223.98, 16977.4, 11351.3 }, { 5.6545, 26.1353, 1.54193e+07 }, { 1e+38, 1e+38, 4.23355e+06 }, { 172.888, 204.112, 403.267 }, { 0.995381, 1.19417, 0.646645 }, { 0.138255, 0.168884, 0.116563 } },
	{ "green-fabric", "fabric-green", { 0.0511324, 0.0490447, 0.0577457 }, { 0.043898, 0.108081, 0.118528 }, { 0.0906425, 0.14646, 0.125546 }, { 0.199121, 0.21946, 0.130311 }, { 0.117671, 0.0797822, 0.0840896 }, { 0.107501, 0.0628391, 0.0668466 }, { 11.1681, 7.43507, 8.69777 }, { 4.65909, 3.72793, 3.72472 }, { 0.00055264, 0.00331292, 0.0119365 }, { 16.8276, 12.7802, 12.1538 }, { 0.153958, 0.17378, 0.291679 }, { 0.000577582, 0.000635885, 0.000729507 } },
	{ "green-latex", "fabric-green-latex", { 0.0885476, 0.13061, 0.0637004 }, { 0.177041, 0.16009, 0.101365 }, { 0.241826, 0.21913, 0.2567 }, { 0.175925, 0.162514, 0.326958 }, { 0.0213854, 0.0498004, 0.0677643 }, { -0.0864353, -0.0518848, -0.00668045 }, { 5.17117, 5.55867, 4.84094 }, { 2.61304, 2.75985, 2.84236 }, { 0.0417628, 0.0352726, 0.0128547 }, { 8.4496, 8.91691, 9.57998 }, { 0.296393, 0.295144, 0.18077 }, { 0.000875063, 0.00103454, 0.000568252 } },
	{ "green-metallic-paint2", "ch-ball-green-metallic", { 0.00536389, 0.0147585, 0.0072232 }, { 0.0553207, 0.0656441, 0.0608999 }, { 0.00131834, 0.00140737, 0.00171711 }, { 0.436009, 0.57969, 0.54703 }, { 0.291808, 0.194127, 0.199549 }, { 0.126008, -0.0681937, 0.0148264 }, { 493.87, 362.85, 317.92 }, { 4.96631e-08, 3.41674e-06, 8.72658e-07 }, { 15.5833, 3.60732e-07, 3.78431e-07 }, { 8.31751, 25.0069, 21.7501 }, { 0.561587, -0.431953, -0.657077 }, { 0.019819, 0.0506927, 0.0279691 } },
	{ "green-metallic-paint", "green-metallic-paint", { 0.00368935, 0.0155555, 0.022272 }, { 0.185621, 0.436002, 0.322925 }, { 0.131402, 0.146271, 0.154061 }, { 1.28366, 0.865104, 0.944013 }, { 0.12483, 0.0443223, 0.0955612 }, { 0.028874, -0.118581, -0.139373 }, { 3.62502, 5.13725, 4.71161 }, { 1.97389, 2.72437, 2.63239 }, { 5.5607e-08, 1.1239e-07, 9.76666e-08 }, { 27.5467, 23.1688, 23.2911 }, { -0.204638, -0.326555, -0.334024 }, { 0.00120994, 0.00217841, 0.00218444 } },
	{ "green-plastic", "green-plastic", { 0.015387, 0.0851675, 0.0947402 }, { 0.0607427, 0.156977, 0.125155 }, { 0.000302146, 0.00197038, 0.000690284 }, { 0.373134, 0.751741, 0.505735 }, { 0.116395, 0.0388464, 0.0476683 }, { 0.0123694, -0.0038785, 0.00864048 }, { 2329.55, 181.586, 833.886 }, { 3.39581, 0.000128411, 3.06753e-07 }, { 6.38801e-08, 2.49751e-07, 5.82332e-07 }, { 260.135, 31.8416, 26.4101 }, { 0.510576, -0.156086, -0.337325 }, { 0.0411204, 0.0490786, 0.0385724 } },
	{ "hematite", "hematite", { 0.00948374, 0.0117628, 0.00985037 }, { 0.0705694, 0.118965, 0.115059 }, { 0.000908552, 0.000601576, 0.00184248 }, { 0.515183, 0.498157, 0.73351 }, { 0.235045, 0.0609324, 0.00720315 }, { -0.120178, -0.096611, -0.235146 }, { 626.075, 967.655, 201.871 }, { 3.94735e-07, 3.54792e-07, 0.000101676 }, { 5.04376e-07, 16.4874, 2.78427e-07 }, { 24.7601, 14.7537, 31.4869 }, { -0.431506, 0.577853, -0.162174 }, { 0.0598105, 0.0605776, 0.0805566 } },
	{ "ipswich-pine-221", "ipswich-pine-221", { 0.0560746, 0.0222518, 0.0105117 }, { 0.0991995, 0.106719, 0.110343 }, { 0.014258, 0.0178759, 0.0188163 }, { 0.0625943, 0.113994, 0.118296 }, { 1.55288e-13, 6.595e-12, 3.97788e-13 }, { -0.0675784, -0.0696373, -0.0703103 }, { 68.7248, 53.3521, 50.6148 }, { 6.36482, 6.37738, 6.35839 }, { 111.962, 1.68378, 0.850892 }, { 20.956, 23.9475, 24.0837 }, { 0.842456, 0.66795, 0.641354 }, { 0.00364018, 0.00320748, 0.00300736 } },
	{ "light-brown-fabric", "", { 0.0612259, 0.0263619, 0.0187761 }, { 3.65487e-12, 9.7449e-12, 4.13685e-12 }, { 0.147778, 0.137639, 0.13071 }, { 0.188292, 0.188374, 0.189026 }, { 1, 1, 1 }, { 0, 0, 0 }, { 7.46192, 7.90085, 8.23842 }, { 3.59765, 3.74493, 3.85476 }, { 0.00660661, 0.0049581, 0.00395337 }, { 12.0657, 12.6487, 13.0977 }, { 0.221546, 0.213332, 0.206743 }, { 0.00126964, 0.000811029, 0.000670968 } },
	{ "light-red-paint", "paint-light-red", { 0.391162, 0.0458387, 0.0059411 }, { 0.522785, 0.144057, 0.214076 }, { 0.854048, 0.0627781, 0.205965 }, { 1.40232, 0.101554, 0.73691 }, { 0.0346142, 0.0655657, 0.0725632 }, { -0.0688995, 0.0400933, 0.0252314 }, { 4.68691, 16.2538, 4.71913 }, { 1.72476, 4.88106, 3.96894 }, { 0.0969212, 0.0124952, 1.02202e-07 }, { 5.57044, 16.1561, 19.9884 }, { 0.249637, 0.388724, -0.505514 }, { 0.00271456, 0.000889762, 0.000505975 } },
	{ "maroon-plastic", "maroon-plastic", { 0.189951, 0.0353828, 0.0321504 }, { 0.127693, 0.100703, 0.115731 }, { 0.00160715, 0.00110827, 0.00100127 }, { 0.684406, 0.645917, 0.569111 }, { 0.0479368, 0.0624437, 0.0921161 }, { -0.0134224, 0.00888653, 0.0423901 }, { 257.032, 398.91, 515.363 }, { 4.57088e-05, 2.62643e-05, 3.15876e-06 }, { 3.05979e-07, 2.76112e-07, 3.67984e-07 }, { 29.9497, 31.8569, 27.4899 }, { -0.211122, -0.159036, -0.309239 }, { 0.0285283, 0.0277773, 0.0286259 } },
	{ "natural-209", "natural-209", { 0.0961753, 0.0349012, 0.0121752 }, { 0.0781649, 0.0898869, 0.111321 }, { 0.0137282, 0.0154247, 0.0233645 }, { 0.0491415, 0.0673163, 0.149331 }, { 0.0522205, 0.0456428, 0.0279324 }, { -0.0277529, -0.0308045, -0.0492147 }, { 71.8965, 63.4418, 40.2168 }, { 6.3511, 6.32932, 6.32009 }, { 224.515, 48.774, 0.0171617 }, { 20.2337, 21.0421, 25.7335 }, { 0.874926, 0.811948, 0.483785 }, { 0.00465507, 0.00387214, 0.00295941 } },
	{ "neoprene-rubber", "neoprene-rubber", { 0.259523, 0.220477, 0.184871 }, { 0.275058, 0.391429, 0.0753145 }, { 0.143818, 0.207586, 0.0764912 }, { 0.770284, 0.774203, 0.700644 }, { 0.113041, 0.110436, 0.16895 }, { 0.060346, 0.0565499, 0.0788468 }, { 5.56845, 4.61088, 8.84784 }, { 3.02411, 3.82334, 2.38942 }, { 5.3042e-08, 1.01087e-07, 6.70643e-08 }, { 23.2109, 20.1103, 28.3099 }, { -0.379072, -0.500903, -0.156114 }, { 0.00218754, 0.00314666, 0.00180987 } },
	{ "nickel", "nickel", { 0.0144009, 0.0115339, 0.00989042 }, { 0.157696, 0.293022, 0.450103 }, { 0.00556292, 0.00627392, 0.00660563 }, { 0.171288, 0.168324, 0.161023 }, { 2.21884, 1.61986, 0.931645 }, { 1.21884, 1.22103, 0.698939 }, { 160.907, 143.234, 136.933 }, { 6.76252, 6.76447, 6.76278 }, { 16179.3, 3451.04, 3453.27 }, { 36.1378, 35.1378, 33.614 }, { 0.846825, 0.821209, 0.830948 }, { 0.132612, 0.177506, 0.155358 } },
	{ "nylon", "nylon", { 0.204199, 0.211192, 0.19234 }, { 0.156797, 0.303324, 0.236394 }, { 0.0250344, 0.0436802, 0.0421753 }, { 0.528875, 0.617086, 0.620808 }, { 0.240279, 0.117115, 0.127421 }, { 0.219191, 0.0940848, 0.101339 }, { 26.4669, 14.8484, 15.213 }, { 2.32486, 2.22618, 2.18663 }, { 1.17309e-07, 5.34378e-08, 6.68165e-08 }, { 39.9707, 33.9774, 34.3152 }, { 0.11793, -0.0190508, -0.00222626 }, { 0.00830775, 0.00885046, 0.00645556 } },
	{ "orange-paint", "orange-paint", { 0.368088, 0.147113, 0.00692426 }, { 0.524979, 0.116386, 0.199437 }, { 0.818115, 0.064743, 0.229391 }, { 1.44385, 0.0709512, 0.483597 }, { 6.92565e-13, 0.106161, 0.102279 }, { -0.174318, 0.0934385, 0.0625648 }, { 4.57466, 16.0185, 4.96427 }, { 1.84547, 4.70387, 3.6232 }, { 0.072629, 0.0299825, 0.000333551 }, { 5.96872, 14.9466, 13.2194 }, { 0.222125, 0.438216, -0.0759733 }, { 0.00178442, 0.000789226, 0.000301814 } },
	{ "pearl-paint", "pearl-paint", { 0.181967, 0.159068, 0.143348 }, { 0.105133, 0.0928717, 0.0802367 }, { 0.0724063, 0.0808503, 0.0596139 }, { 0.194454, 0.203296, 0.091958 }, { 0.168966, 0.297431, 0.401185 }, { -0.831034, -0.702569, -0.226546 }, { 13.6335, 12.3113, 17.1256 }, { 5.07525, 4.9097, 4.92468 }, { 0.000252814, 0.000258641, 0.0180128 }, { 18.9547, 18.2047, 16.1366 }, { 0.156731, 0.135808, 0.415562 }, { 0.00244402, 0.00167924, 0.0012928 } },
	{ "pickled-oak-260", "pickled-oak-260", { 0.181735, 0.14142, 0.125486 }, { 0.0283411, 0.0296418, 0.025815 }, { 0.0105853, 0.0102771, 0.0101188 }, { 2.31337e-14, 2.35272e-14, 1.99762e-14 }, { 5.38184e-13, 2.15933e-13, 3.55496e-12 }, { -0.309259, -0.291046, -0.329625 }, { 95.4759, 98.3089, 99.831 }, { 6.37433, 6.39032, 6.39857 }, { 4641.15, 5970.15, 6818.76 }, { 18.1707, 18.2121, 18.2334 }, { 1.00563, 1.01274, 1.01645 }, { 0.00336548, 0.0032891, 0.00330829 } },
	{ "pink-fabric2", "", { 0.24261, 0.0829238, 0.0751196 }, { 0.161823, 0.0591236, 0.00907967 }, { 0.220011, 0.148623, 0.111966 }, { 6.44517e-12, 3.24286e-13, 3.83556e-11 }, { 0.242032, 0.456181, 3.11925 }, { -0.0582931, 0.295844, 2.84268 }, { 5.66376, 7.80657, 9.98941 }, { 2.43742, 3.05804, 3.53392 }, { 0.111404, 0.079067, 0.0738938 }, { 7.47166, 9.23588, 10.5653 }, { 0.408014, 0.422779, 0.448866 }, { 0.00145188, 0.000903357, 0.00103814 } },
	{ "pink-fabric", "fabric-pink", { 0.270553, 0.223977, 0.240993 }, { 0.299998, 0.418074, 4.07112e-13 }, { 0.787023, 0.234345, 0.17346 }, { 1.77629, 4.64947e-15, 0.203846 }, { 0.121124, 0.0271041, 1 }, { 0.0151107, 0.00273969, 0 }, { 4.71724, 5.3941, 6.55218 }, { 2.23471, 2.34423, 3.31887 }, { 0.0266324, 0.11955, 0.00966123 }, { 7.20596, 7.20326, 11.0696 }, { 0.130136, 0.407562, 0.22282 }, { 0.00129262, 0.00124216, 0.00262093 } },
	{ "pink-felt", "fabric-pink-felt", { 0.259533, 0.192978, 0.185581 }, { 0.359813, 0.533498, 0.0390541 }, { 0.46679, 0.203504, 0.314933 }, { 0.663613, 1.70267e-13, 0.778919 }, { 0.0530603, 0.0108612, 1.12976 }, { -0.124124, -0.0455941, 1.05206 }, { 3.6266, 6.02293, 3.92458 }, { 2.22042, 2.55548, 3.68735 }, { 0.0333497, 0.102524, 8.10122e-05 }, { 7.40195, 7.81026, 13.5108 }, { 0.185696, 0.409228, -0.236166 }, { 0.000823299, 0.00122612, 0.00141316 } },
	{ "pink-jasper", "pink-jasper", { 0.226234, 0.138929, 0.110785 }, { 0.0846118, 0.0984038, 0.078693 }, { 0.00223592, 0.00203213, 0.0013737 }, { 0.53138, 0.4995, 0.436111 }, { 0.172825, 0.129714, 0.154397 }, { 0.0911362, 0.0599408, 0.0851167 }, { 253.003, 292.741, 474.2 }, { 1.64528, 1.86649, 5.29965e-08 }, { 5.71153e-08, 6.24765e-08, 15.5452 }, { 110.886, 113.953, 8.17872 }, { 0.416268, 0.422653, 0.561605 }, { 0.0166667, 0.0194625, 0.0142451 } },
	{ "pink-plastic", "pink-plastic", { 0.354572, 0.0905002, 0.0696372 }, { 0.0316585, 0.0444153, 6.36158e-15 }, { 0.0566727, 0.0369011, 0.142628 }, { 0.215149, 0.0621896, 0.189286 }, { 0.317781, 0.0925394, 1 }, { 0.266863, 0.0652821, 0 }, { 16.788, 27.3043, 7.67451 }, { 5.62766, 5.4919, 3.67439 }, { 1.71214e-05, 0.187553, 0.00561921 }, { 23.0822, 17.5369, 12.3737 }, { 0.0833775, 0.577298, 0.215985 }, { 0.00219041, 0.000905537, 0.00219415 } },
	{ "polyethylene", "polyethylene", { 0.228049, 0.239339, 0.240326 }, { 0.0420869, 0.134269, 0.0867928 }, { 0.0472725, 0.260465, 0.0719615 }, { 2.09548e-13, 0.743064, 2.27838e-13 }, { 0.489907, 0.316434, 0.181688 }, { 0.424297, 0.244434, 0.0772152 }, { 22.178, 4.24558, 14.9332 }, { 4.92064, 4.37149, 4.27422 }, { 0.330549, 9.78063e-07, 0.109488 }, { 14.3439, 17.178, 12.5991 }, { 0.61045, -0.470479, 0.517652 }, { 0.00285999, 0.00147486, 0.00191841 } },
	{ "polyurethane-foam", "polyurethane-foam", { 0.0898318, 0.0428583, 0.0340984 }, { 4.0852e-12, 8.0217e-14, 4.05682e-14 }, { 0.154984, 0.142104, 0.139418 }, { 0.188586, 0.188095, 0.188124 }, { 1, 1, 1 }, { 0, 0, 0 }, { 7.18433, 7.70043, 7.81983 }, { 3.50087, 3.67779, 3.7174 }, { 0.0079259, 0.00567317, 0.00525063 }, { 11.6922, 12.3795, 12.537 }, { 0.226848, 0.217315, 0.215132 }, { 0.00119032, 0.000823112, 0.000753708 } },
	{ "pure-rubber", "pure-rubber", { 0.284259, 0.251873, 0.223824 }, { 0.542899, 0.598765, 0.142162 }, { 0.62899, 0.41413, 0.0693873 }, { 1.13687, 0.452063, 0.262952 }, { 0.0185379, 0.00420901, 0.0353132 }, { -0.0265343, -0.035723, 0.0236736 }, { 3.75577, 3.75969, 13.6137 }, { 2.17187, 2.13584, 5.56595 }, { 0.0311422, 0.0526357, 7.8942e-07 }, { 7.14885, 7.0586, 23.9264 }, { 0.152489, 0.245615, -0.0979145 }, { 0.00177602, 0.00171837, 0.0007762 } },
	{ "purple-paint", "purple-paint", { 0.290743, 0.0347118, 0.0339802 }, { 0.301308, 0.205258, 0.280717 }, { 0.0339173, 0.00839425, 0.0378011 }, { 0.703279, 0.20083, 0.797376 }, { 0.0654958, 0.0786512, 0.0696424 }, { 0.0299494, 0.0467651, 0.0373649 }, { 16.3705, 104.698, 13.2101 }, { 1.72814, 6.77684, 1.54679 }, { 1.10548e-07, 0.177314, 6.66798e-08 }, { 38.0847, 41.1712, 38.2781 }, { 0.0789738, 0.576929, 0.0587144 }, { 0.00296394, 0.00277256, 0.00225685 } },
	{ "pvc", "pvc", { 0.0322978, 0.0357449, 0.0403426 }, { 0.28767, 0.317369, 0.310067 }, { 0.0171547, 0.0176681, 0.0213663 }, { 0.769726, 0.730637, 0.797555 }, { 0.0289552, 0.026258, 0.0281305 }, { 0.00815373, 0.00651989, 0.00714485 }, { 25.8231, 26.9026, 20.5927 }, { 1.20095, 1.31883, 1.23655 }, { 1.02768e-07, 1.09767e-07, 5.97648e-08 }, { 51.5975, 50.0184, 48.1953 }, { 0.218049, 0.208814, 0.173894 }, { 0.00471012, 0.00620207, 0.00575443 } },
	{ "red-fabric2", "", { 0.155216, 0.0226757, 0.0116884 }, { 1.80657e-15, 5.51946e-13, 1.35221e-15 }, { 0.16689, 0.135884, 0.128307 }, { 0.184631, 0.18856, 0.1883 }, { 1, 1, 1 }, { 0, 0, 0 }, { 6.78759, 7.98303, 8.36701 }, { 3.33819, 3.77225, 3.89063 }, { 0.0112096, 0.00468671, 0.00372538 }, { 11.0557, 12.7596, 13.2393 }, { 0.240935, 0.21164, 0.206003 }, { 0.00263226, 0.00060129, 0.000538747 } },
	{ "red-fabric", "fabric-red", { 0.201899, 0.0279008, 0.0103965 }, { 0.168669, 0.0486346, 0.040485 }, { 0.324447, 0.228455, 0.109436 }, { 0.787411, 0.821197, 0.279212 }, { 0.0718348, 0.0644687, 0.0206123 }, { -0.0585917, -0.0062547, -0.050402 }, { 3.88129, 4.32067, 9.16355 }, { 3.59825, 3.93046, 4.66379 }, { 0.000130047, 1.04152e-07, 4.59182e-05 }, { 13.0776, 19.649, 17.8852 }, { -0.209387, -0.530789, -0.0242035 }, { 0.00103676, 0.000248481, 0.000220698 } },
	{ "red-metallic-paint", "ch-ball-red-metallic", { 0.0380897, 0.00540095, 0.00281156 }, { 0.0416724, 0.07642, 0.108438 }, { 0.00133258, 0.00106883, 0.00128863 }, { 0.693854, 0.52857, 0.539477 }, { 2.45718e-16, 0.0598671, 0.0633332 }, { -1, -0.08904, -0.0114056 }, { 300.371, 521.418, 425.982 }, { 6.45857e-05, 6.3446e-07, 8.51754e-07 }, { 2.75773e-07, 4.05125e-07, 3.41703e-07 }, { 33.0213, 24.3646, 23.5793 }, { -0.121415, -0.469753, -0.532083 }, { 0.287672, 0.0467824, 0.0286558 } },
	{ "red-phenolic", "red-phenolic", { 0.165227, 0.0256259, 0.00935644 }, { 0.240561, 0.360634, 0.475777 }, { 0.0052844, 0.00467439, 0.00613717 }, { 0.568938, 0.509763, 0.575762 }, { 0.156419, 0.0972193, 0.069671 }, { 0.0752589, 0.0444558, 0.0266428 }, { 104.336, 128.839, 89.6357 }, { 1.57629, 1.92067, 1.57462 }, { 6.38793e-08, 6.71457e-08, 5.25172e-08 }, { 77.4088, 79.3795, 73.0158 }, { 0.343482, 0.352913, 0.324912 }, { 0.0680783, 0.0859239, 0.0936197 } },
	{ "red-plastic", "red-plastic", { 0.247569, 0.049382, 0.0175621 }, { 0.406976, 0.151478, 0.176348 }, { 0.28723, 0.0572489, 0.0624682 }, { 0.939617, 0.0851973, 0.0701483 }, { 0.10036, 0.178468, 0.149441 }, { 0.0512697, 0.13191, 0.102958 }, { 3.80564, 17.8374, 16.5631 }, { 4.37031, 4.96092, 4.75924 }, { 1.03822e-07, 0.0236845, 0.0321546 }, { 18.6088, 16.1371, 15.1243 }, { -0.609047, 0.43547, 0.445882 }, { 0.00134715, 0.000536577, 0.00047104 } },
	{ "red-specular-plastic", "red-plastic-specular", { 0.252589, 0.0397665, 0.0185317 }, { 0.0139957, 0.0343278, 0.0527973 }, { 6.01746e-05, 8.07327e-05, 0.000205705 }, { 0.174569, 0.202455, 0.390522 }, { 0.441328, 0.179378, 0.150221 }, { 0.191312, 0.0691835, 0.0760833 }, { 14623.1, 10620.2, 3332.98 }, { 22.6631, 12.9299, 3.06408 }, { 1e+38, 1e+38, 6.36761e-08 }, { 187.799, 162.245, 315.067 }, { 1.19812, 1.20136, 0.52102 }, { 0.0620224, 0.0541448, 0.0836154 } },
	{ "silicon-nitrade", "silicon-nitrade", { 0.0141611, 0.0115865, 0.00842477 }, { 0.0710113, 0.0670906, 0.015769 }, { 6.40406e-05, 0.000138867, 0.00224354 }, { 0.159422, 0.283527, 0.734323 }, { 0.0516164, 0.10318, 2.36643 }, { -0.0277792, -0.0531505, 1.36643 }, { 13926.7, 5668.82, 167.824 }, { 31.0579, 5.63563, 8.52534e-05 }, { 1e+38, 1e+38, 2.96855e-07 }, { 176.231, 148.745, 29.0282 }, { 1.20757, 1.14062, -0.241968 }, { 0.0712429, 0.0720965, 0.11404 } },
	{ "silver-metallic-paint2", "", { 0.0554792, 0.0573803, 0.0563376 }, { 0.121338, 0.115673, 0.10966 }, { 0.029859, 0.0303706, 0.0358666 }, { 0.144097, 0.104489, 0.158163 }, { 1.03749e-14, 3.52034e-15, 4.41778e-12 }, { -1, -1, -1 }, { 31.9005, 32.1514, 26.5685 }, { 6.08248, 5.89319, 5.95156 }, { 0.00761403, 0.0839948, 0.00138703 }, { 23.5891, 20.6786, 23.3297 }, { 0.435405, 0.54017, 0.34665 }, { 0.021384, 0.0187248, 0.0175762 } },
	{ "silver-metallic-paint", "metallic-silver", { 0.0189497, 0.0205686, 0.0228822 }, { 0.173533, 0.168901, 0.165266 }, { 0.037822, 0.038145, 0.0381908 }, { 0.165579, 0.162955, 0.160835 }, { 5.66903e-12, 1.65276e-14, 4.28399e-14 }, { -1, -1, -1 }, { 25.1551, 24.9957, 25.0003 }, { 5.92591, 5.90225, 5.89007 }, { 0.000684235, 0.000840795, 0.000995926 }, { 23.4725, 23.1898, 23.0144 }, { 0.310575, 0.317996, 0.324938 }, { 0.0122324, 0.0117672, 0.0113004 } },
	{ "silver-paint", "paint-silver", { 0.152796, 0.124616, 0.113375 }, { 0.30418, 0.30146, 0.283174 }, { 0.110819, 0.105318, 0.0785677 }, { 0.640378, 0.641115, 0.445228 }, { 2.37347e-13, 7.68194e-13, 2.9434e-12 }, { -0.350607, -0.355433, -0.359297 }, { 7.21531, 7.46519, 10.8289 }, { 3.06016, 2.97652, 3.78287 }, { 9.71666e-08, 1.09342e-07, 9.82336e-08 }, { 24.1475, 24.5083, 25.6757 }, { -0.281384, -0.257633, -0.200968 }, { 0.00361589, 0.00384995, 0.00444223 } },
	{ "special-walnut-224", "special-walnut-224", { 0.0121712, 0.00732998, 0.00463072 }, { 0.209603, 0.216118, 0.211885 }, { 0.117091, 0.119932, 0.131119 }, { 0.548899, 0.524858, 0.569425 }, { 0.0808859, 0.0802614, 0.0789982 }, { 0.0327605, 0.0324012, 0.0274637 }, { 7.42314, 7.41578, 6.77215 }, { 3.61532, 3.82778, 3.71096 }, { 1.19182e-07, 1.03098e-07, 1.08004e-07 }, { 22.9756, 22.7318, 22.2976 }, { -0.31158, -0.332257, -0.35462 }, { 0.000934753, 0.000878135, 0.000758769 } },
	{ "specular-black-phenolic", "black-bball", { 0.00212164, 0.00308282, 0.00410253 }, { 0.0881574, 0.0923246, 0.0398117 }, { 0.00119167, 0.000641898, 0.000186605 }, { 0.616914, 0.578026, 0.40121 }, { 0.122486, 0.0907984, 0.13855 }, { 0.0482851, 0.0348218, 0.0481813 }, { 395.656, 781.361, 3615.49 }, { 1.18051e-05, 1.22149, 2.87626 }, { 2.88932e-07, 6.09017e-08, 9.79363e-08 }, { 28.9748, 200.529, 331.336 }, { -0.257695, 0.487197, 0.524707 }, { 0.0311874, 0.0437836, 0.0853211 } },
	{ "specular-blue-phenolic", "blue-bball", { 0.00497564, 0.0138836, 0.032815 }, { 0.1077, 0.0898232, 0.175296 }, { 0.000918571, 0.0010348, 0.00176322 }, { 0.570978, 0.639916, 0.709385 }, { 0.0354139, 0.0488958, 0.023159 }, { -0.0434314, -0.0294427, -0.039938 }, { 558.482, 431.708, 222.547 }, { 3.57872e-06, 2.38985e-05, 6.87316e-05 }, { 3.33356e-07, 3.3336e-07, 3.2851e-07 }, { 28.5289, 32.1362, 30.3802 }, { -0.272244, -0.141224, -0.19005 }, { 0.0324614, 0.0338654, 0.057732 } },
	{ "specular-green-phenolic", "green-bball", { 0.00781782, 0.0259654, 0.0233739 }, { 0.0688449, 0.144658, 0.143654 }, { 0.000307494, 0.0010353, 0.00155331 }, { 0.365481, 0.585805, 0.787512 }, { 0.129429, 0.0443676, 5.423e-12 }, { 0.047932, -0.0136055, -0.0592743 }, { 2313.47, 482.847, 206.646 }, { 3.55591, 5.3306e-06, 0.000279403 }, { 1.03804e-07, 3.47154e-07, 3.16888e-07 }, { 256.625, 28.2901, 38.0691 }, { 0.511901, -0.277015, 0.00513058 }, { 0.0373526, 0.0450242, 0.0470867 } },
	{ "specular-maroon-phenolic", "maroon-bball", { 0.152486, 0.0263216, 0.00802748 }, { 0.0761775, 0.098375, 0.165913 }, { 0.000342958, 0.000605578, 0.00144136 }, { 0.4052, 0.553617, 0.65133 }, { 0.0646024, 0.0116325, 0.037551 }, { -0.0555983, -0.0527264, -0.0283844 }, { 1961.35, 868.119, 306.537 }, { 2.81751, 1.361, 2.43015e-05 }, { 5.70302e-08, 5.89096e-08, 3.27578e-07 }, { 248.591, 203.707, 29.053 }, { 0.506511, 0.488814, -0.242746 }, { 0.0465516, 0.0434438, 0.0545231 } },
	{ "specular-orange-phenolic", "orange-bball", { 0.32771, 0.0540131, 0.00883213 }, { 0.051915, 0.0686764, 0.0489478 }, { 7.91913e-05, 0.000139576, 1.62017e-05 }, { 0.253564, 0.354675, 3.55583e-12 }, { 0.0768695, 0.0496641, 0.0223538 }, { -0.015514, -0.023864, -0.0433847 }, { 10274.5, 5159.16, 61722.9 }, { 7.4207, 3.83488, 8.72207e+06 }, { 1e+38, 1.02522e-07, 2.30124e+07 }, { 170.364, 373.255, 855.811 }, { 1.19475, 0.53082, 0.608547 }, { 0.0405677, 0.0498947, 0.0953035 } },
	{ "specular-red-phenolic", "red-bball", { 0.303563, 0.0354891, 0.00899721 }, { 0.151819, 0.0938022, 0.196935 }, { 0.00117843, 0.00056476, 0.00185124 }, { 0.570146, 0.524406, 0.732785 }, { 0.0287445, 0.0672098, 0.0178433 }, { -0.0386172, 0.0144813, -0.0322308 }, { 438.968, 982.441, 201.328 }, { 2.88657e-06, 8.12418e-07, 9.95633e-05 }, { 3.8194e-07, 3.69661e-07, 3.03179e-07 }, { 26.0046, 30.0091, 31.3125 }, { -0.376008, -0.21799, -0.162894 }, { 0.0307055, 0.0199408, 0.0652829 } },
	{ "specular-violet-phenolic", "violet-bball", { 0.0686035, 0.0181856, 0.0210368 }, { 0.108459, 0.0471612, 0.171691 }, { 0.00123271, 0.000443974, 0.00149517 }, { 0.657484, 0.546753, 0.653065 }, { 0.0403569, 0.121081, 0.035323 }, { -0.0295013, 0.0563904, -0.0275623 }, { 351.208, 1193.45, 294.897 }, { 3.17585e-05, 1.3817, 2.44051e-05 }, { 3.02028e-07, 6.19706e-08, 3.40809e-07 }, { 31.3319, 234.879, 28.7237 }, { -0.168991, 0.500354, -0.252626 }, { 0.033584, 0.0495535, 0.0510203 } },
	{ "specular-white-phenolic", "white-bball", { 0.282896, 0.231703, 0.127818 }, { 0.0678467, 0.0683808, 0.032756 }, { 2.00918e-05, 3.22307e-05, 9.59333e-05 }, { 4.35028e-11, 0.131903, 0.390338 }, { 0.0132613, 0.0416359, 0.285115 }, { -0.1025, -0.0828341, 0.0759745 }, { 49772.5, 28318.4, 7131.79 }, { 4.58788e+06, 1.12234e+07, 3.05696 }, { 4.5813e+06, 1.00372e+07, 6.93295e-08 }, { 767.24, 646.93, 455.184 }, { 0.60999, 0.61958, 0.536789 }, { 0.134616, 0.114414, 0.179676 } },
	{ "specular-yellow-phenolic", "yellow-bball", { 0.309395, 0.135278, 0.0159106 }, { 0.0607659, 0.141526, 0.110839 }, { 0.000200174, 0.00179326, 0.00048094 }, { 0.355381, 0.767285, 0.505029 }, { 0.077039, 0.0175574, 0.038586 }, { -0.0182146, -0.0417745, -0.0114057 }, { 3597.35, 190.656, 1192.15 }, { 3.80702, 0.000181166, 1.69208 }, { 1.01181e-07, 3.0765e-07, 1.27637e-07 }, { 313.8, 34.1523, 221.08 }, { 0.523001, -0.0813947, 0.500075 }, { 0.0451868, 0.0491835, 0.031808 } },
	{ "ss440", "ss440", { 0.0229923, 0.0187037, 0.0153204 }, { 0.127809, 0.14899, 0.0473376 }, { 4.06782e-05, 4.86278e-05, 8.15367e-05 }, { 0.0888931, 0.259135, 0.332727 }, { 2.26479e-12, 5.7223e-13, 0.820577 }, { -0.474442, -0.223179, -0.179423 }, { 23197.8, 16625.7, 9083.17 }, { 1.43395e+07, 7.33446, 4.46921 }, { 1.06287e+07, 1e+38, 9.97795e-08 }, { 556.069, 215.021, 479.875 }, { 0.627834, 1.20059, 0.540075 }, { 0.366952, 0.463741, 0.874071 } },
	{ "steel", "steel", { 0.019973, 0.0127074, 0.0246402 }, { 0.0615275, 0.0469644, 0.0402151 }, { 8.63865e-05, 0.000249576, 5.77865e-05 }, { 0.3729, 0.583679, 0.471406 }, { 0.665193, 0.798139, 0.0115189 }, { -0.334807, -0.201861, -0.988481 }, { 8119.74, 1951.31, 10374.2 }, { 3.40501, 1.12086, 1.88549 }, { 1.25096e-07, 9.51373e-08, 1.20412e-07 }, { 474.668, 313.81, 603.334 }, { 0.539696, 0.519433, 0.545491 }, { 1.11996, 1.28262, 1.42004 } },
	{ "teflon", "teflon", { 0.276442, 0.263098, 0.260294 }, { 1.56924, 1.52804, 1.43859 }, { 0.678586, 0.662167, 0.577852 }, { 1.2402, 1.14126, 1.44077 }, { 9.40662e-13, 3.57656e-11, 1.50208e-11 }, { -0.0492032, -0.0587316, -0.0548028 }, { 3.91532, 3.82834, 3.64944 }, { 2.09871, 2.05044, 2.90315 }, { 0.0378525, 0.0435869, 0.0028469 }, { 6.8689, 6.72229, 9.65585 }, { 0.166634, 0.181605, -0.0468321 }, { 0.00462131, 0.00463132, 0.00408386 } },
	{ "tungsten-carbide", "tungsten-carbide", { 0.0151872, 0.0103016, 0.0123192 }, { 0.0504358, 0.075701, 0.0556673 }, { 6.6122e-05, 7.65809e-05, 4.80196e-05 }, { 0.255291, 0.270824, 0.26732 }, { 5.30357e-14, 5.44537e-12, 1.09586e-11 }, { -1, -1, -1 }, { 12280.7, 10423.3, 16683.5 }, { 7.42976, 6.36068, 6.81602 }, { 1e+38, 1e+38, 1e+38 }, { 185.623, 173.962, 218.975 }, { 1.197, 1.19586, 1.19745 }, { 0.472021, 0.583364, 0.850304 } },
	{ "two-layer-gold", "", { 0.0415046, 0.0312801, 0.0253658 }, { 1.58161, 1.18736, 1.63847 }, { 0.0263104, 0.0293804, 0.0241265 }, { 0.355682, 0.354281, 0.36415 }, { 0.117355, 0.0614942, 0.0447004 }, { 0.0411678, -0.0237579, -0.0100488 }, { 30.4478, 27.5367, 32.7278 }, { 3.91606, 3.99968, 3.79114 }, { 9.38513e-08, 6.62806e-08, 6.79906e-08 }, { 36.9091, 35.6578, 38.4915 }, { 0.0815181, 0.0465613, 0.0922721 }, { 0.165383, 0.160769, 0.168742 } },
	{ "two-layer-silver", "", { 0.0657916, 0.0595705, 0.0581288 }, { 1.55275, 2.00145, 1.93045 }, { 0.0149977, 0.0201665, 0.0225062 }, { 0.382631, 0.35975, 0.361657 }, { 4.93242e-13, 1.00098e-14, 0.0103259 }, { -0.0401315, -0.0395054, -0.0312454 }, { 50.1263, 38.8508, 34.9978 }, { 3.41873, 3.77545, 3.78138 }, { 6.09709e-08, 1.02036e-07, 1.01016e-07 }, { 46.6236, 40.8229, 39.1812 }, { 0.183797, 0.139103, 0.117092 }, { 0.170639, 0.189329, 0.21468 } },
	{ "violet-acrylic", "acrylic-violet", { 0.0599875, 0.023817, 0.0379025 }, { 0.134984, 0.13337, 0.295509 }, { 0.0011295, 0.00126481, 0.00186818 }, { 0.523244, 0.551606, 0.478149 }, { 0.100325, 0.0939057, 0.0663939 }, { 0.0603016, 0.0579001, 0.0353439 }, { 498.715, 424.307, 328.834 }, { 4.79917e-07, 1.39397e-06, 2.03376 }, { 4.27008e-07, 3.4655e-07, 1.05578e-07 }, { 23.5924, 24.3831, 116.7 }, { -0.515008, -0.476979, 0.432199 }, { 0.041711, 0.0411698, 0.0691096 } },
	{ "violet-rubber", "violet-rubber", { 0.223179, 0.0553634, 0.113238 }, { 0.547456, 0.0966027, 0.185463 }, { 0.445092, 0.089421, 0.1543 }, { 0.923481, 0.782897, 0.800907 }, { 0.0518927, 0.113508, 0.094549 }, { -0.0208955, 0.0585063, 0.036751 }, { 3.58673, 7.33539, 5.24404 }, { 2.83007, 2.27788, 3.00009 }, { 0.00374266, 9.65628e-08, 1.13746e-07 }, { 9.67534, 27.1898, 22.3519 }, { -0.0037936, -0.174481, -0.363226 }, { 0.00274248, 0.000994259, 0.00101746 } },
	{ "white-acrylic", "acrylic-white", { 0.314106, 0.300008, 0.263648 }, { 0.015339, 0.0736169, 0.0976209 }, { 0.00233914, 0.00168956, 0.00117717 }, { 0.570501, 0.51337, 0.440543 }, { 2.651, 0.38473, 0.249988 }, { 2.27139, 0.341493, 0.218093 }, { 226.134, 342.423, 548.396 }, { 1.41218, 2.17179e-07, 4.52325e-08 }, { 9.83895e-08, 5.84338e-06, 15.7173 }, { 110.209, 18.8595, 8.78916 }, { 0.419684, -0.61884, 0.56243 }, { 0.0537471, 0.0401864, 0.0400114 } },
	{ "white-diffuse-bball", "white-diffuse-bball", { 0.284726, 0.239199, 0.177227 }, { 0.680658, 0.828508, 0.502128 }, { 0.689731, 0.601478, 0.209084 }, { 1.23085, 0.956833, 0.800488 }, { 1.06465e-11, 7.45307e-19, 0.089058 }, { -0.109723, -0.090129, 0.0545089 }, { 3.94153, 3.65714, 4.53186 }, { 2.04845, 2.06125, 3.72838 }, { 0.0434345, 0.0439287, 1.06879e-07 }, { 6.69695, 6.79535, 20.1466 }, { 0.178983, 0.18898, -0.494774 }, { 0.007671, 0.0072205, 0.00441788 } },
	{ "white-fabric2", "", { 0.10784, 0.102669, 0.113943 }, { 0.0375359, 0.0296317, 0.0364218 }, { 0.0583526, 0.0520763, 0.063905 }, { 9.1836e-12, 5.41897e-13, 3.47708e-12 }, { 0.0528459, 0.0854122, 0.0860966 }, { 0.0578589, 0.0971296, 0.0985717 }, { 18.1669, 20.2291, 16.6809 }, { 4.60692, 4.77932, 4.46485 }, { 0.172556, 0.23927, 0.138435 }, { 13.5006, 13.9648, 13.1165 }, { 0.561147, 0.587176, 0.54153 }, { 0.00113864, 0.00114867, 0.00119185 } },
	{ "white-fabric", "", { 0.290107, 0.219835, 0.160654 }, { 0.230066, 0.156787, 2.19005e-12 }, { 0.479844, 0.196767, 0.162438 }, { 0.669662, 2.22894e-13, 0.194582 }, { 7.17773e-14, 0.0364632, 1 }, { -0.209237, -0.0528699, 0 }, { 3.60967, 6.18732, 6.91147 }, { 2.17095, 2.60737, 3.42418 }, { 0.0379118, 0.0990846, 0.00860197 }, { 7.221, 7.95861, 11.4289 }, { 0.195346, 0.409992, 0.225019 }, { 0.000728813, 0.000806025, 0.00177083 } },
	{ "white-marble", "", { 0.236183, 0.221746, 0.192889 }, { 0.24075, 0.221456, 0.23389 }, { 0.00430204, 0.00429388, 0.00405907 }, { 0.681306, 0.676061, 0.620701 }, { 0.118569, 0.167963, 0.140415 }, { 0.0620704, 0.102032, 0.071968 }, { 103.053, 104.312, 122.169 }, { 1.0396, 1.05866, 1.26372 }, { 1.08288e-07, 6.82995e-08, 5.76023e-08 }, { 88.7766, 89.0261, 88.95 }, { 0.378101, 0.372681, 0.371975 }, { 0.0493991, 0.0689784, 0.0584648 } },
	{ "white-paint", "white-paint", { 0.356024, 0.3536, 0.324889 }, { 4.14785, 0.255488, 0.108438 }, { 0.0654126, 0.0841905, 0.538225 }, { 0.796927, 0.7778, 0.888627 }, { 0.0770738, 0.25087, 1.75144e-13 }, { 0.0782297, 0.226285, -0.367441 }, { 8.89985, 7.66809, 3.59697 }, { 1.96428, 2.23069, 2.23051 }, { 5.9674e-08, 1.07673e-07, 0.0280623 }, { 31.0512, 27.6872, 7.41492 }, { -0.0896936, -0.151411, 0.154604 }, { 0.160631, 0.0657026, 0.0710504 } },
	{ "yellow-matte-plastic", "yellow-matte-plastic", { 0.276745, 0.108557, 0.0203686 }, { 0.806628, 1.99624, 0.977002 }, { 0.0237573, 0.0265756, 0.0305873 }, { 0.558101, 0.537647, 0.710679 }, { 0.0624565, 0.0605404, 0.0541077 }, { 0.049915, 0.0567466, 0.0473569 }, { 26.7398, 24.8437, 17.6137 }, { 2.1481, 2.32064, 1.65758 }, { 5.99553e-08, 5.96103e-08, 6.38097e-08 }, { 41.7723, 39.7008, 40.1979 }, { 0.11569, 0.0891281, 0.0886802 }, { 0.0469972, 0.10918, 0.0466045 } },
	{ "yellow-paint", "paint-yellow", { 0.288876, 0.195348, 0.0314583 }, { 0.449392, 0.412812, 0.168707 }, { 0.650734, 0.190849, 0.16131 }, { 1.21986, 0.0333524, 0.137577 }, { 0.00201415, 7.31639e-13, 0.13481 }, { -0.0897316, -0.0639292, 0.103378 }, { 3.83308, 6.31195, 7.06429 }, { 2.18358, 2.72018, 3.26561 }, { 0.0298708, 0.0764778, 0.0198928 }, { 7.1694, 8.39507, 10.5682 }, { 0.146537, 0.387991, 0.292185 }, { 0.00152911, 0.00123385, 0.000424783 } },
	{ "yellow-phenolic", "phenolic-yellow", { 0.26924, 0.190177, 0.0858303 }, { 0.0861694, 0.0960246, 0.122709 }, { 0.00126171, 0.00197611, 0.00187166 }, { 0.444222, 0.512534, 0.504223 }, { 0.0988294, 0.120249, 0.0929719 }, { 0.0295258, 0.0539044, 0.0471952 }, { 509.34, 294.418, 314.704 }, { 4.59604e-08, 1.76251, 1.81658 }, { 15.6472, 6.5671e-08, 9.96466e-08 }, { 8.52565, 116.04, 118.064 }, { 0.561753, 0.425536, 0.432459 }, { 0.0238647, 0.0190162, 0.0173945 } },
	{ "yellow-plastic", "yellow-plastic", { 0.221083, 0.193042, 0.0403393 }, { 0.265199, 0.340361, 0.0670333 }, { 0.280789, 0.146396, 0.0248514 }, { 0.920018, 0.0356237, 7.80502e-12 }, { 0.144471, 0.10845, 0.125678 }, { 0.113592, 0.0663735, 0.103644 }, { 3.85808, 7.84753, 41.2517 }, { 4.34432, 3.16869, 5.71439 }, { 1.04213e-07, 0.055929, 7.84273 }, { 18.6813, 9.73614, 16.4495 }, { -0.60265, 0.393386, 0.781122 }, { 0.00125483, 0.00140731, 0.00101943 } }
};*/
template <typename Float>
static double 
sgd__g1(const Vector<Float, 3>& k, double theta0, double c, double k_, double lambda)
{
	double tmp1 = max(0.0, acos(k.z) - theta0);
	double tmp2 = 1.0 - exp(c * pow(tmp1, k_));
	double tmp3 = 1.0 + lambda * tmp2;
	return dr::minimum(dr::maximum(tmp3, 0.0f),1.0f);
}

static double sgd__ndf(double cos_theta_h, double alpha, double p, double kap)
{
	const double inv_pi = 1.0 / M_PI;
	double c2 = cos_theta_h * cos_theta_h;
	double t2 = (1.0 - c2) / c2;
	double ax = alpha + t2 / alpha;

	return ((kap * exp(-ax) * inv_pi) / (pow(ax, p) * c2 * c2));
}

// -------------------------------------------------------------------------------------------------
// SGD Ctor
template <typename Float>
sgd<Float>::sgd(const char *name): m_fresnel(NULL), m_data(NULL)
{
	bool found = false;
	for (uint32_t i = 0; i < (uint32_t)(sizeof(sgd::s_data) / sizeof(sgd::data)); ++i) {
		if (!strcmp(s_data[i].name, name)
		    || !strcmp(s_data[i].otherName, name)) {
			m_data = &s_data[i];
			m_fresnel = new fresnel::sgd(Vector3f::from_raw(m_data->f0),
			                             Vector3f::from_raw(m_data->f1));
			found = true;
			break;
		}
	}
	if (!found) throw exc("djb_error: No SGD parameters for %s\n", name);
}

// -------------------------------------------------------------------------------------------------
// SGD eval
/*
template <typename Float>
Vector<Float, 3> sgd<Float>::eval(const Vector<Float, 3>& i, const Vector<Float, 3>& o, const void *user_param) const
{
	if (i.z > 0.0 && o.z > 0.0) {
		Vector3f h = normalize(i + o);
		const Vector3f Ks = Vector3f::from_raw(m_data->rhoS);
		const Vector3f Kd = Vector3f::from_raw(m_data->rhoD);
		Vector3f F = fresnel(sat(dot(i, h)));
		Vector3f G = gaf(h, i, o);
		Vector3f D = ndf(h);

		return ((Kd + Ks * (F * D * G) / (i.z * o.z)) / dr::Pi<Float>);
	}
	return Vector3f(0);
}

// -------------------------------------------------------------------------------------------------
// SGD Shadowing
template <typename Float>
Vector<Float, 3> sgd<Float>::gaf(const Vector<Float, 3>& h, const Vector<Float, 3>& i, const Vector<Float, 3>& o) const
{
	return g1(i) * g1(o);
}
template <typename Float>
Vector<Float, 3> sgd<Float>::g1(const Vector<Float, 3>& k) const
{
	double g1[3];

	for (uint32_t i = 0; i < 3; ++i)
		g1[i] = sgd__g1(k, m_data->theta0[i], m_data->c[i],
	                    m_data->k[i], m_data->lambda[i]);

	return Vector3f::from_raw(g1);
}

// -------------------------------------------------------------------------------------------------
// SGD NDF
template <typename Float>
Vector<Float, 3> sgd<Float>::ndf(const Vector<Float, 3>& h) const
{
	double cos_theta_h = h.z;
	double ndf[3];

	for (uint32_t i = 0; i < 3; ++i)
		ndf[i] = sgd__ndf(cos_theta_h, m_data->alpha[i],
	                      m_data->p[i], m_data->kap[i]);

	return Vector3f::from_raw(ndf);
}
*/

// *************************************************************************************************
// ABC Distribution API implementation
/*
template <typename Float>
const abc<Float>::data abc<Float>::s_data[] = {
	{ "alum-bronze", {0.024514, 0.022346, 0.016220}, {51.714554, 37.932693, 27.373095}, 10482.133789, 0.816737, 2.236525 },
	{ "alumina-oxide", {0.310854, 0.288744, 0.253423}, {2191.683838, 1851.525513, 2064.303467}, 97381.554688, 1.393850, 1.437706 },
	{ "aluminium", {0.000000, 0.002017, 0.005639}, {1137.956787, 1090.751831, 1113.879150}, 218716.593750, 1.092006, 6.282688 },
	{ "aventurnine", {0.052933, 0.060724, 0.052465}, {1835.591309, 1594.715820, 1585.083984}, 324913.375000, 1.110341, 1.403899 },
	{ "beige-fabric", {0.202460, 0.132107, 0.110584}, {0.146439, 0.133698, 0.120796}, 3442009.500000, 0.083210, 10.406549 },
	{ "black-fabric", {0.014945, 0.008790, 0.008178}, {3.109159, 2.772291, 2.723374}, 871444.875000, 0.292300, 1.057058 },
	{ "black-obsidian", {0.000000, 0.000000, 0.000000}, {1447.612183, 1281.373169, 1313.115601}, 223549.906250, 1.196914, 1.485313 },
	{ "black-oxidized-steel", {0.012978, 0.011743, 0.011219}, {1.100971, 1.057305, 1.015528}, 30.661055, 1.869546, 1.660413 },
	{ "black-phenolic", {0.001445, 0.001674, 0.001885}, {104.792976, 99.654602, 100.537247}, 8696.112305, 1.524949, 1.839669 },
	{ "black-soft-plastic", {0.000666, 0.000594, 0.000723}, {0.449065, 0.442911, 0.443035}, 124.908409, 0.618811, 1.902569 },
	{ "blue-acrylic", {0.011106, 0.035057, 0.104254}, {1191.397339, 1099.257202, 1082.125366}, 255868.953125, 1.084736, 1.475061 },
	{ "blue-fabric", {0.008485, 0.009937, 0.034319}, {0.262721, 0.290408, 0.465468}, 3115416.500000, 0.138407, 2.581633 },
	{ "blue-metallic-paint", {0.003770, 0.001030, 0.009177}, {0.495104, 0.469843, 1.097868}, 48.572060, 1.702520, 9.024582 },
	{ "blue-metallic-paint2", {0.000000, 0.000000, 0.000000}, {399.504974, 589.572571, 1115.380737}, 68107.390625, 1.124658, 2.756175 },
	{ "blue-rubber", {0.035403, 0.074082, 0.154563}, {1.295714, 1.279542, 1.248921}, 37.265717, 1.565906, 1.448611 },
	{ "brass", {0.007722, 0.012504, 0.010908}, {2975.156982, 1931.137207, 987.877380}, 244596.531250, 1.136163, 2.698691 },
	{ "cherry-235", {0.037077, 0.014564, 0.007235}, {2.240743, 2.211937, 2.200735}, 123.542030, 1.498179, 1.741225 },
	{ "chrome-steel", {0.006340, 0.008416, 0.009893}, {795.266418, 629.395203, 620.750549}, 122360.437500, 1.486198, 99.999619 },
	{ "chrome", {0.002692, 0.002520, 0.003088}, {986.846741, 771.035522, 696.134583}, 64661.257812, 1.905420, 99.999413 },
	{ "colonial-maple-223", {0.081471, 0.023697, 0.008944}, {1.003711, 0.932037, 0.910946}, 24.846296, 2.121997, 1.805962 },
	{ "color-changing-paint1", {0.002224, 0.002598, 0.005303}, {32.486813, 34.621838, 36.775242}, 2456.161865, 1.388325, 2.066753 },
	{ "color-changing-paint2", {0.000368, 0.000022, 0.001103}, {16.119473, 13.414398, 13.766625}, 1099.414307, 1.331730, 3.186646 },
	{ "color-changing-paint3", {0.007087, 0.004775, 0.003475}, {35.495728, 34.895355, 33.480820}, 2193.942139, 1.442522, 1.918318 },
	{ "dark-blue-paint", {0.002488, 0.007395, 0.034947}, {0.541177, 0.465818, 0.482633}, 24.620041, 1.411947, 1.882264 },
	{ "dark-red-paint", {0.251814, 0.030605, 0.009098}, {0.457430, 0.406366, 0.394383}, 34.332043, 1.023539, 1.670872 },
	{ "dark-specular-fabric", {0.020414, 0.009048, 0.007769}, {1.659937, 1.561381, 1.496425}, 130.493088, 0.975626, 1.356871 },
	{ "delrin", {0.296617, 0.284276, 0.247152}, {4.629758, 4.546722, 4.481231}, 526.069458, 0.931238, 1.545662 },
	{ "fruitwood-241", {0.048981, 0.034689, 0.019745}, {5.976665, 5.760932, 5.568172}, 192.199890, 1.573300, 1.483582 },
	{ "gold-metallic-paint", {0.009164, 0.005649, 0.000424}, {1.685630, 1.147257, 0.418757}, 53.738525, 1.511730, 9.031934 },
	{ "gold-metallic-paint2", {0.000000, 0.000000, 0.000000}, {204.653580, 176.726898, 145.123520}, 1705396.875000, 0.621691, 8.599911 },
	{ "gold-metallic-paint3", {0.000000, 0.000000, 0.000000}, {712.306885, 534.449829, 299.714020}, 94549.195312, 0.970436, 3.965659 },
	{ "gold-paint", {0.150463, 0.080224, 0.022239}, {1.010258, 0.811582, 0.504538}, 46.942558, 1.756670, 7.936110 },
	{ "gray-plastic", {0.099830, 0.100258, 0.092187}, {208.924332, 204.230804, 203.086151}, 33710.691406, 1.023289, 1.400441 },
	{ "grease-covered-steel", {0.000000, 0.000000, 0.000000}, {292.647400, 277.053131, 268.833618}, 15585.240234, 1.170786, 2.556668 },
	{ "green-acrylic", {0.016581, 0.074299, 0.040329}, {1279.877197, 1196.967407, 1221.283569}, 59359.058594, 1.407153, 1.354611 },
	{ "green-fabric", {0.040464, 0.035315, 0.041862}, {0.655209, 0.935078, 1.071461}, 1602934.125000, 0.158285, 1.639683 },
	{ "green-latex", {0.071622, 0.106434, 0.047407}, {0.618660, 0.699535, 0.427560}, 701732.250000, 0.230061, 3.772501 },
	{ "green-metallic-paint", {0.000000, 0.027206, 0.035274}, {0.983370, 1.731232, 1.827718}, 43.999603, 2.172307, 3.235001 },
	{ "green-metallic-paint2", {0.000000, 0.000000, 0.000000}, {289.186859, 578.364441, 361.170837}, 117416.273438, 0.980195, 2.199395 },
	{ "green-plastic", {0.012262, 0.079555, 0.087255}, {2172.868896, 1950.600586, 1812.148193}, 365361.062500, 1.158465, 1.428372 },
	{ "hematite", {0.000000, 0.000273, 0.000528}, {1261.585449, 1214.305420, 1255.108643}, 453901.000000, 1.011091, 2.480358 },
	{ "ipswich-pine-221", {0.049757, 0.016449, 0.004662}, {1.680928, 1.649491, 1.647754}, 49.356201, 1.853272, 1.729069 },
	{ "light-brown-fabric", {0.048293, 0.019506, 0.013650}, {1.725058, 1.219095, 1.080752}, 3929432.500000, 0.155265, 1.153189 },
	{ "light-red-paint", {0.410883, 0.038975, 0.006158}, {0.591418, 0.516404, 0.539583}, 12.604274, 2.137619, 1.743125 },
	{ "maroon-plastic", {0.192636, 0.033163, 0.030345}, {2021.205688, 1751.832031, 1723.350952}, 255986.625000, 1.178735, 1.438042 },
	{ "natural-209", {0.080513, 0.023515, 0.003459}, {2.211772, 2.104441, 1.964800}, 91.766762, 1.154241, 1.641148 },
	{ "neoprene-rubber", {0.248773, 0.211228, 0.170494}, {1.327917, 1.299042, 1.233595}, 63.131874, 1.701130, 1.942141 },
	{ "nickel", {0.006227, 0.006663, 0.007587}, {36.614742, 32.745403, 28.906111}, 705.733887, 1.945258, 5.441883 },
	{ "nylon", {0.232688, 0.232487, 0.207358}, {19.215639, 18.194090, 16.938326}, 3334.451904, 0.886752, 1.382529 },
	{ "orange-paint", {0.407889, 0.131526, 0.003770}, {0.451945, 0.373510, 0.347662}, 11.196602, 1.972428, 1.887030 },
	{ "pearl-paint", {0.181452, 0.162150, 0.136895}, {0.843725, 0.714731, 0.493903}, 30.882542, 1.941556, 9.303603 },
	{ "pickled-oak-260", {0.127652, 0.093767, 0.077318}, {0.813234, 0.803229, 0.781871}, 67.082710, 1.625856, 2.933847 },
	{ "pink-fabric", {0.252317, 0.191945, 0.194732}, {0.434443, 0.450356, 0.464137}, 239167.375000, 0.146117, 2.710645 },
	{ "pink-fabric2", {0.210380, 0.046071, 0.034554}, {0.970204, 0.643877, 0.617551}, 1448228.625000, 0.193577, 4.907917 },
	{ "pink-felt", {0.267287, 0.178590, 0.156370}, {0.225846, 0.226874, 0.221785}, 54.853024, 0.344321, 2.938402 },
	{ "pink-jasper", {0.196444, 0.119217, 0.093004}, {372.695374, 357.903107, 363.584320}, 149508.703125, 0.934382, 1.652987 },
	{ "pink-plastic", {0.360193, 0.079077, 0.053058}, {0.826036, 0.561110, 0.579990}, 16122.769531, 0.217444, 1.810468 },
	{ "polyethylene", {0.183218, 0.202033, 0.207279}, {0.839841, 0.740593, 0.637590}, 3082111.000000, 0.118759, 3.010911 },
	{ "polyurethane-foam", {0.057868, 0.017982, 0.010879}, {0.146946, 0.101932, 0.092976}, 5723581.500000, 0.073071, 3.375618 },
	{ "pure-rubber", {0.287230, 0.252144, 0.211784}, {0.741122, 0.716288, 0.661539}, 85.876984, 0.982928, 1.670395 },
	{ "purple-paint", {0.280217, 0.031190, 0.032274}, {14.006363, 13.769886, 13.753972}, 408.684143, 1.776584, 1.483770 },
	{ "pvc", {0.030285, 0.033992, 0.038152}, {40.608692, 39.046635, 37.916531}, 2901.928223, 1.203405, 1.407985 },
	{ "red-fabric", {0.185341, 0.020015, 0.003522}, {0.934393, 0.314632, 0.273278}, 130451.187500, 0.287254, 3.335066 },
	{ "red-fabric2", {0.116616, 0.012705, 0.002205}, {0.586502, 0.141217, 0.134042}, 8092575.000000, 0.084146, 1.976978 },
	{ "red-metallic-paint", {0.000000, 0.000000, 0.000000}, {2077.145508, 423.511139, 278.621643}, 161943.781250, 1.015078, 2.346517 },
	{ "red-phenolic", {0.162631, 0.021944, 0.006764}, {198.859695, 194.122971, 198.083405}, 57814.757812, 1.049066, 1.752476 },
	{ "red-plastic", {0.268000, 0.044432, 0.010675}, {1.081538, 1.045352, 1.037089}, 34.117554, 1.529797, 1.666450 },
	{ "red-specular-plastic", {0.246119, 0.035688, 0.015152}, {2890.193359, 2449.908936, 2543.227783}, 215481.906250, 1.268530, 1.322196 },
	{ "silicon-nitrade", {0.002805, 0.003669, 0.003936}, {1185.676392, 1031.834106, 1048.810913}, 185127.203125, 1.146923, 1.828749 },
	{ "silver-metallic-paint", {0.000342, 0.000025, 0.002221}, {2.419921, 2.373467, 2.319122}, 58.200935, 1.718832, 9.073233 },
	{ "silver-metallic-paint2", {0.036282, 0.048095, 0.055908}, {1.964917, 1.771562, 1.583858}, 84.018394, 1.141540, 8.961876 },
	{ "silver-paint", {0.159599, 0.127910, 0.108671}, {1.151248, 1.207139, 1.248789}, 38.884895, 1.971779, 8.336617 },
	{ "special-walnut-224", {0.008121, 0.004822, 0.003478}, {0.704096, 0.693490, 0.676845}, 15.779606, 2.668091, 1.883156 },
	{ "specular-black-phenolic", {0.001379, 0.001737, 0.001891}, {2866.697021, 2568.570068, 3042.953613}, 291094.218750, 1.252478, 1.415415 },
	{ "specular-blue-phenolic", {0.002601, 0.011165, 0.029456}, {3080.328613, 2655.796143, 2536.662354}, 298749.312500, 1.201234, 1.437517 },
	{ "specular-green-phenolic", {0.004986, 0.023480, 0.020383}, {3187.704102, 2609.277832, 2708.220947}, 267365.562500, 1.224307, 1.411174 },
	{ "specular-maroon-phenolic", {0.150765, 0.024140, 0.006117}, {2987.999268, 2544.122314, 2451.749023}, 261203.859375, 1.225481, 1.424407 },
	{ "specular-orange-phenolic", {0.325956, 0.051359, 0.005049}, {3129.503174, 2678.228516, 2803.262207}, 462669.000000, 1.185505, 1.439050 },
	{ "specular-red-phenolic", {0.301706, 0.032497, 0.006908}, {2943.967773, 2402.380127, 2540.010010}, 261627.484375, 1.209261, 1.410007 },
	{ "specular-violet-phenolic", {0.066982, 0.015686, 0.018350}, {3117.315918, 2600.714111, 2717.788330}, 257952.890625, 1.251725, 1.407573 },
	{ "specular-white-phenolic", {0.279240, 0.225120, 0.122981}, {4513.249023, 3816.075195, 3782.943848}, 265093.843750, 1.318991, 1.426619 },
	{ "specular-yellow-phenolic", {0.305828, 0.133078, 0.012394}, {3565.679932, 2820.262207, 2784.792236}, 275647.718750, 1.221645, 1.393108 },
	{ "ss440", {0.003789, 0.004251, 0.004992}, {758.115112, 613.785278, 553.360596}, 96436.515625, 1.550240, 99.999336 },
	{ "steel", {0.006151, 0.008683, 0.011281}, {5709.654297, 4641.991211, 4419.819824}, 234651.921875, 1.401774, 2.658189 },
	{ "teflon", {0.321152, 0.313597, 0.305526}, {1.787717, 1.816174, 1.846979}, 41.644257, 2.066234, 1.611215 },
	{ "tungsten-carbide", {0.003085, 0.003349, 0.003638}, {808.095520, 620.232666, 606.583435}, 90806.242188, 1.700902, 99.999184 },
	{ "two-layer-gold", {0.000000, 0.000000, 0.000000}, {103.922806, 91.724258, 76.325027}, 481560.750000, 0.631262, 8.633638 },
	{ "two-layer-silver", {0.000000, 0.000000, 0.000000}, {125.852753, 121.561691, 119.139694}, 1048834.875000, 0.598443, 8.718435 },
	{ "violet-acrylic", {0.044077, 0.006396, 0.023810}, {428.049377, 362.420013, 431.304535}, 1148630.875000, 0.752451, 2.462801 },
	{ "violet-rubber", {0.232500, 0.045230, 0.108831}, {0.840866, 0.796264, 0.822774}, 371.841339, 0.646749, 2.143676 },
	{ "white-acrylic", {0.316720, 0.300911, 0.263206}, {1379.722168, 1291.507080, 1357.878052}, 153332.265625, 1.153521, 1.291205 },
	{ "white-diffuse-bball", {0.319956, 0.272448, 0.178031}, {2.104850, 2.082304, 2.011517}, 75.010391, 1.940200, 1.693726 },
	{ "white-fabric", {0.331064, 0.233945, 0.160661}, {3.545784, 3.587269, 3.573341}, 18457.023438, 0.103981, 1.039433 },
	{ "white-fabric2", {0.074490, 0.068798, 0.080224}, {0.078094, 0.081885, 0.083064}, 6166446.000000, 0.023514, 2.997413 },
	{ "white-marble", {0.222139, 0.208435, 0.183381}, {349.574860, 331.314941, 323.574066}, 49167.527344, 1.153413, 1.650249 },
	{ "white-paint", {0.320307, 0.312274, 0.305093}, {159.141998, 146.118881, 139.397903}, 2777.948730, 1.813009, 1.197418 },
	{ "yellow-matte-plastic", {0.290588, 0.113993, 0.019811}, {85.230881, 80.641701, 78.267593}, 6072.900879, 1.145334, 1.315196 },
	{ "yellow-paint", {0.301059, 0.196365, 0.025507}, {0.446085, 0.446681, 0.399212}, 12.881507, 2.137449, 1.822657 },
	{ "yellow-phenolic", {0.291521, 0.210540, 0.097144}, {498.945862, 477.551605, 480.963745}, 114769.390625, 1.006011, 1.543437 },
	{ "yellow-plastic", {0.237134, 0.202702, 0.030379}, {0.997742, 0.967993, 0.828363}, 289.454254, 0.522704, 1.596116 }
};*/

static double abc__ndf(double cos_theta_h, double A, double B, double C)
{
	double tmp = 1.0 - cos_theta_h;

	return (A / pow(1.0 + B * tmp, C));
}

// -------------------------------------------------------------------------------------------------
// ABC Ctor
template <typename Float>
abc<Float>::abc(const char *name): m_fresnel(NULL), m_data(NULL)
{
	bool found = false;
	for (uint32_t i = 0; i < (uint32_t)(sizeof(abc::s_data) / sizeof(abc::data)); ++i) {
		if (!strcmp(s_data[i].name, name)) {
			m_data = &s_data[i];
			m_fresnel = new fresnel::unpolarized(Vector3f(m_data->ior)),
			found = true;
			break;
		}
	}
	if (!found) throw exc("djb_error: No ABC parameters for %s\n", name);
}

// -------------------------------------------------------------------------------------------------
// ABC eval
/*
template <typename Float>
Vector<Float, 3> abc<Float>::eval(const Vector3f& i, const Vector3f& o, const void *user_param) const
{
	if (i.z > 0.0 && o.z > 0.0) {
		Vector3f h = normalize(i + o);
		const Vector3f Kd = Vector3f::from_raw(m_data->kD);
		Vector3f F = fresnel(sat(dot(i, h)));
		double G = gaf(h, i, o);
		Vector3f D = ndf(h);

		return (Kd / dr::Pi<Float>  + (F * D * G) / (dr::Pi<Float> * i.z * o.z));
	}
	return Vector3f(0);
}*/

// -------------------------------------------------------------------------------------------------
// ABC Shadowing
template <typename Float>
Float abc<Float>::gaf(const Vector3f& h, const Vector3f& i, const Vector3f& o) const
{
	Float g1_i = min((Float)1.0, (Float)2.0 * (h.z * i.z / dot(h, i)));
	Float g1_o = min((Float)1.0, (Float)2.0 * (h.z * o.z / dot(h, o)));

	return min(g1_i, g1_o);
}

// -------------------------------------------------------------------------------------------------
// ABC NDF
/*
template <typename Float>
Vector<Float, 3> abc<Float>::ndf(const Vector3f& h) const
{
	double cos_theta_h = h.z;
	double ndf[3];

	for (uint32_t i = 0; i < 3; ++i)
		ndf[i] = abc__ndf(cos_theta_h, m_data->A[i], m_data->B, m_data->C);

	return Vector3f::from_raw(ndf);
}*/

 // namespace djb
NAMESPACE_END(mitsuba)

#endif // DJ_BRDF_IMPLEMENTATION

