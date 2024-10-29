#pragma once


#ifndef DJB_INCLUDE_DJ_BRDF_H
#define DJB_INCLUDE_DJ_BRDF_H

#include <mitsuba/core/properties.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/core/mmap.h>
#include <mitsuba/core/mstream.h>
#include <drjit/dynamic.h>
#ifndef DJB_ASSERT
#	include <assert.h>
#	define DJB_ASSERT(x) assert(x)
#endif

NAMESPACE_BEGIN(mitsuba)

template <typename Float>
class brdf {
public:
	// evaluate f_r
    MI_IMPORT_CORE_TYPES();
	virtual Vector3f eval(const Vector3f& i, const Vector3f& o,
	                  const void *user_param = NULL) const = 0;
	
	// evaluate f_r * cos
	virtual Vector3f evalp(const Vector3f& i, const Vector3f& o,
	                   const void *user_param = NULL) const;
	
	// importance sample f_r * cos using two uniform numbers
    // evaluate the PDF of a sample
	virtual Float pdf(const Vector3f& i, const Vector3f& o,
	                    const void *user_param = NULL) const;
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
class merl: public brdf<Float>{
	MI_IMPORT_CORE_TYPES()
	 using FloatStorage                = DynamicBuffer<Float>;
	std::vector<double> m_samples;
	FloatStorage m_data;
public:
   
	uint32_t n, dims[3];
    merl() = default;
	merl(const fs::path  &filename){
	
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
Float pdf(const Vector3f& i, const Vector3f& o,
	                    const void *user_param = NULL) const{}
};

template <typename Float>
class microfacet : public brdf<Float> {
public:
MI_IMPORT_CORE_TYPES();
	/* microfacet parameters */
	class params {
		friend class microfacet;
	public:
		// Factories
        params(Float a) {}; 
        static params isotropic(Float a);
        params(Float a1 = 1.0, Float a2 = 1.0, Float phi_a = 0.0);
        static params elliptic(Float a1, Float a2, Float phi_a = 0.0);
        void set_location(Float tx_n, Float ty_n);
        void set_ellipse(Float a1, Float a2, Float phi_a = 0.0);
        
		static params standard();
		/*
		static params elliptic(float_t a1, float_t a2, float_t phi_a = 0.0);
		static params pdfparams(float_t ax, float_t ay, float_t rho = 0.0,
		                        float_t tx_n = 0.0, float_t ty_n = 0.0);
		// mutators
		void set_ellipse(float_t a1, float_t a2, float_t phi_a = 0.0);
		void set_pdfparams(float_t ax, float_t ay, float_t rho = 0.0,
		                   float_t tx_n = 0.0, float_t ty_n = 0.0);
		void set_location(float_t tx_n, float_t ty_n);
		void set_location(const vec3& n);
		// accessors
		void get_ellipse(float_t *a1, float_t *a2, float_t *phi_a = NULL) const;
		void get_pdfparams(float_t *ax, float_t *ay, float_t *rho = NULL,
		                   float_t *tx_n = NULL, float_t *ty_n = NULL) const;
		void get_location(float_t *tx_n, float_t *ty_n) const;
		void get_location(vec3 *n) const;
		// Ctors (prefer factories for construction)
		params(float_t a1 = 1.0, float_t a2 = 1.0, float_t phi_a = 0.0);
		params(float_t ax, float_t ay, float_t rho, float_t tx_n, float_t ty_n);*/
	private:
		Vector3f m_n; // mean normal
		Float m_a1, m_a2, m_phi_a; // ellipse parameters
		float m_ax, m_ay; // scale parameters
		Float m_rho, m_sqrt_one_minus_rho_sqr; // correlation
		Float m_tx_n, m_ty_n; // location parameters
	};
	// Dtor
    microfacet() = default;
    microfacet(bool shadow=true): m_shadow(shadow) {}
   virtual ~microfacet() {}
//	virtual ~microfacet() {delete m_fresnel;}
	// BRDF interface
    
	Vector3f eval(const Vector3f& i, const Vector3f& o,
	          const void *user_param = NULL) const{}
	Vector3f evalp(const Vector3f& i, const Vector3f& o,
	           const void *user_param = NULL) const{}
    Float pdf(const Vector3f& i, const Vector3f& o,
	            const void *user_param = NULL);

    Float gaf(const Vector3f& h, const Vector3f& i, const Vector3f& o,
	            const params& params = params::standard());
	Float g1(const Vector3f& h, const Vector3f& k,
	           const params& params = params::standard()) ;
	Float sigma(const Vector3f& k,
	              const params& params = params::standard()) ;         

    Float ndf(const Vector3f& h,
	            const params& params = params::standard()) ;     
    Float vndf(const Vector3f& h, const Vector3f& k,
	             const params& params = params::standard()) ;      
    Float p22(Float x, Float y,
	            const params& params = params::standard());     
    virtual bool supports_smith_vndf_sampling() const = 0;            
    protected:
    virtual Float p22_std(Float x, Float y); // standard slope pdf    
    virtual Float sigma_std(const Vector3f k); // standard masking                     
	/*vec3 sample(float_t u1, float_t u2, const vec3& o,
	            const void *user_param = NULL) const;
	float_t pdf(const vec3& i, const vec3& o,
	            const void *user_param = NULL) const;
	vec3 evalp_is(float_t u1, float_t u2, const vec3& o,
	              vec3 *i, float_t *pdf, const void *user_param = NULL) const;
	// eval queries
	vec3 fresnel(float cos_theta_d) const {return m_fresnel->eval(cos_theta_d);}
	float_t ndf(const vec3& h,
	            const params& params = params::standard()) const;
	float_t gaf(const vec3& h, const vec3& i, const vec3& o,
	            const params& params = params::standard()) const;
	float_t g1(const vec3& h, const vec3& k,
	           const params& params = params::standard()) const;
	float_t sigma(const vec3& k,
	              const params& params = params::standard()) const;
	float_t p22(float_t x, float_t y,
	            const params& params = params::standard()) const;
	float_t vp22(float_t x, float_t y, const vec3& k,
	             const params& params = params::standard()) const;
	float_t vndf(const vec3& h, const vec3& k,
	             const params& params = params::standard()) const;
	// sampling queries
	virtual bool supports_smith_vndf_sampling() const = 0;
	virtual float_t qf2(float u, const vec3& k) const;
	virtual float_t qf3(float u, const vec3& k, float qf2) const;
	// mutators
	void set_shadow(bool shadow) {m_shadow = shadow;}
	void set_fresnel(const fresnel::impl& f);
	// accessors
	int get_shadow() const {return m_shadow;}
	const fresnel::impl& get_fresnel() const {return *m_fresnel;}
protected:
	// ctor
	microfacet(const fresnel::impl& f = fresnel::ideal(),
	           bool shadow = true);
	// microfacet interfaces
	virtual float_t sigma_std(const vec3& k) const = 0; // standard masking
	virtual float_t p22_std(float_t x, float_t y) const = 0; // standard slope pdf
	// sampling functions
	virtual void sample_vp22_std_smith(float_t u1, float_t u2, const vec3& k,
	                                   float_t *xslope, float_t *yslope) const;
	virtual void sample_vp22_std_nmap(float_t u1, float_t u2, const vec3& k,
	                                  float_t *xslope, float_t *yslope) const = 0;*/
	// members
//	const fresnel::impl *m_fresnel;
	bool m_shadow;
};

template <typename Float>
class radial : public microfacet<Float> {
public:
    MI_IMPORT_CORE_TYPES();
	radial(bool shadow = true): microfacet<Float>(shadow)
	{}
	// queries
	virtual Float p22_radial(Float r_sqr);
	virtual Float sigma_std_radial(Float cos_theta_k);
private:
	// eval
	Float p22_std(Float x, Float y){
	return p22_radial(x * x + y * y);
      }
    Float sigma_std(Vector3f k){
	return sigma_std_radial(k.z());
}
    };

template <typename Float>
class ggx : public radial<Float> {
public:
ggx(bool shadow = true): radial<Float>(shadow)
{}
bool supports_smith_vndf_sampling() const {return true;}
Float p22_radial(Float r_sqr) {
	Float tmp = 1.0 + r_sqr;
	return (1.0 / (M_PI * tmp * tmp));
}
Float sigma_std_radial(Float cos_theta_k);
};



template <typename Float>
class tabular : public radial<Float> 
{
	std::vector<Float> m_p22;
	std::vector<float> m_sigma;
	// tables for Nmap VNDF sampling
	std::vector<float> m_cdf;
	std::vector<float> m_qf;
public:
   MI_IMPORT_CORE_TYPES();
	tabular(const brdf<Float>& brdf, int32_t resolution, bool shadow = true):radial<Float>()
{
	DJB_ASSERT(res > 2 && "Invalid Resolution");
	//set_shadow(shadow);

	// allocate memory
//	m_p22.resize(resolution,0.0f);
//	m_sigma.reserve(res);
//	m_cdf.reserve(res);
//	m_qf.reserve(res);

	// eval
	//compute_p22_smith(brdf, res);
	//normalize_p22();
	//compute_sigma();
	//compute_fresnel(brdf, res);

	// sample
//	compute_cdf();
//	compute_qf();
}


	static typename microfacet<Float>::params fit_ggx_parameters(const tabular<Float>& tabular);
   static Float fit_ggx_parameter(const tabular<Float>& tabular);
	// queries
	Float p22_radial(Float r_sqr);
	Float sigma_std_radial(Float cos_theta_k)const;
	float cdf_radial(float r) const;
	float qf_radial(float u) const;
	bool supports_smith_vndf_sampling() const {return false;}
private:
	// eval precomputations
	void compute_p22_smith(const brdf<Float>& brdf, int res);
	void compute_fresnel(const brdf<Float>& brdf, int res);
	void normalize_p22();
	void compute_sigma();
	// sample precomputations
	void compute_cdf();
	void compute_qf();
	void compute_pdf1();
	void normalize_pdf1();
};



NAMESPACE_END(mitsuba)

//
//
//// end header file /////////////////////////////////////////////////////
#endif // DJB_INCLUDE_DJ_BRDF_H

#if DJ_BRDF_IMPLEMENTATION
#include <cmath>
#include <cstdarg>
//#include <iostream>     // std::ios, std::istream, std::cout
//#include <fstream>      // std::filebuf
//#include <cstring>      // memcpy
#include <stdint.h>     // uint32_t
#include <mitsuba/core/warp.h>
#include <mitsuba/core/util.h>
#include <drjit/dynamic.h>
#include <array>
#ifndef DJB_ASSERTm_tabular
#	include <assert.h>
#	define DJB_ASSERT(x) assert(x)
#endif

#ifndef DJB_EPSILON
#define DJB_EPSILON (float_t)1e-4
#endif

#ifndef M_PI
#	define M_PI 3.1415926535897932384626433832795
#endif

#ifdef _MSC_VER
#	pragma warning(disable: 4244) // possible loss of data
#endif


NAMESPACE_BEGIN(mitsuba)
// Constructor

template <typename Float>

Float intensity( Vector<Float, 3> inty) { return 0.2126 * inty.x() + 0.7152 * inty.y() + 0.0722 * inty.z();}

template <typename Float>
class matrix {
	std::vector<Float> mij;
	int size;
public:
	matrix(int size);
	Float& operator()(int i, int j) {return mij[j*size+i];}
	const Float& operator()(int i, int j) const {return mij[j*size+i];}
	void transform(const std::vector<Float>& v, std::vector<Float>& out) const;
	std::vector<Float> eigenvector(int iterations) const;
};
template <typename Float>
matrix<Float>::matrix(int size) : mij(size * size, 0), size(size)
{}
template <typename Float>
void matrix<Float>::transform(const std::vector<Float>& v, std::vector<Float>& out) const
{
	out.resize(0);
	for (int j = 0; j < size; ++j) {
		out.push_back(0);
		for (int i = 0; i < size; ++i) {
			out[j]+= (*this)(i, j) * v[i];
		}
	}
}
template <typename Float>
std::vector<Float> matrix<Float>::eigenvector(int iterations) const
{
	int j = 0;
	std::vector<Float> vec[2];
	//vec[0].reserve(size);
	//vec[1].reserve(size);
    vec[0].resize(size,0.0f);
    vec[1].resize(size,0.0f);
	for (int i = 0; i < size; ++i)
		vec[j][i] = 1.0;
	for (int i = 0; i < iterations; ++i) {
		transform(vec[j], vec[1-j]);
		j = 1 - j;
	}
	return vec[j];
}

template <typename Float>
Float
microfacet<Float>::pdf(const Vector3f& i, const Vector3f& o, const void *user_param) 
{
	const microfacet<Float>::params params = 
		user_param ? *reinterpret_cast<const microfacet<Float>::params *>(user_param)
		           : microfacet<Float>::params::standard();
	Vector3f h = dr::normalize(i + o);
	Float G = gaf(h, i, o, params);
    Mask active = true;
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

template <typename Float>
Float
microfacet<Float>::sigma(
	const Vector3f& k,
	const microfacet<Float>::params& params
) {
	// warp the incident direction
	Float a = k.x() * params.m_ax + k.y() * params.m_ay * params.m_rho;
	Float b = k.y() * params.m_ay * params.m_sqrt_one_minus_rho_sqr;
	Float c = k.z() - k.x() * params.m_tx_n - k.y() * params.m_ty_n;
	Float nrm = dr::sqrt(a * a + b * b + c * c);
    Vector3f nrml=Vector3f(a, b, c) / nrm;    
	return nrm * sigma_std(nrml);
    
}

template <typename Float>
Float
microfacet<Float>::g1(
	const Vector3f& h, const Vector3f& k,
	const microfacet<Float>::params& params
)  {
	Float g1_local_test = dr::dot(k, params.m_n);
   Mask active = true;
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
) {
	Float G1_o = g1(h, o, params);

	if (m_shadow) {
		Float G1_i = g1(h, i, params);
		Float tmp = G1_i * G1_o;
#if 1 // height correlated smith
        return dr::maximum(0.0, tmp / (G1_i + G1_o - tmp));
		/*if (tmp > 0.0)
			return (tmp / (G1_i + G1_o - tmp));

		return 0.0; // fully attenuated*/
#else // separable form
		return tmp;
#endif
	}

	return G1_o;
}
// Microfacet NDF
template <typename Float>
Float microfacet<Float>::ndf(const Vector3f& h, const microfacet::params& params)
{ Mask active = true;
    active &= h.z() > DJB_EPSILON;
	if (dr::none_or<true>(active)) {
		Float cos_theta_h_sqr = h.z() * h.z();
		Float cos_theta_h_sqr_sqr = cos_theta_h_sqr * cos_theta_h_sqr;
		Float xslope = -h.x() / h.z();
		Float yslope = -h.y() / h.z();

		return (p22(xslope, yslope, params) / cos_theta_h_sqr_sqr);
	}
	return 0.0;
}

//---------------------------------------------------------------------------


template <typename Float>
Float microfacet<Float>::p22(Float x, Float y, const params& params) 
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
// Microfacet VNDF
template <typename Float>
Float
microfacet<Float>::vndf(
	const Vector3f& h, const Vector3f& k,
	const microfacet<Float>::params& params
) {
	Float kh = dr::dot(k, h);
    Mask active = true;
    active &= kh > 0.0;
	if (dr::none_or<true>(active)) {
		Float D = ndf(h, params);

		return (kh * D / sigma(k, params));
	}
	return 0.0;
}


/* Tabulated Microfacet API implementation*/

template <typename Float>
Float tabular<Float>::p22_radial(Float r_sqr) 
{
	Float r = dr::sqrt(r_sqr);
	Float u = dr::sqrt(2.0 * dr::atan(r) / M_PI);
 //  return spline::eval<Float>(m_p22, "edge", u);
   MI_IMPORT_CORE_TYPES();
	int edge   = (int)m_p22.size();
    using FloatStorage                = DynamicBuffer<Float>;
	FloatStorage s_data;

  std::unique_ptr<Float[]> data_out(new Float[edge]);
	Float *data_out_ptr = data_out.get();
	Float *points_ptr = m_p22.data();
	 for (uint32_t k = 0; k < edge; ++k)
        *data_out_ptr++ = *points_ptr++ ;
	s_data = dr::load<FloatStorage>(data_out.get(), (edge));
    
    UInt32 intpart=UInt32 ((u * edge)-u);
    Float frac=u * edge -u - dr::floor((u * edge)-u); 

   
    UInt32 i1 = dr::clamp(intpart, 0,edge-1);
	UInt32 i2 = dr::clamp(intpart + 1,0, edge-1);
   

	const Float &p1 = dr::gather<Float>(s_data, i1);
	const Float &p2 = dr::gather<Float>(s_data, i2);
   // Throw("Invalid file structure: %f", u);
	return dr::lerp(p1, p2, frac);
    

	
}
/*
template <typename Float>
float tabular<Float>::sigma_std_radial(float cos_theta_k) const
{
	float u = (float)2.0 * acos(cos_theta_k) / (float)M_PI;
	return spline::eval(m_sigma, spline::uwrap_edge, u);
}
template <typename Float>
float tabular<Float>::cdf_radial(float r) const
{
	float u = atan(r) * (float)2.0 / (float)M_PI;
	if (u < (float_t)0.0) u = (float_t)0.0;
	return spline::eval(m_cdf, spline::uwrap_edge, sqrt(u));
}
template <typename Float>
float tabular<Float>::qf_radial(float u) const
{
	DJB_ASSERT(u > (float)0.0 && u < (float)1.0);
	float qf = spline::eval(m_qf, spline::uwrap_edge, u);
	return tan(qf * (float)M_PI / (float)2.0);
}
*/



/*template <typename Float>
tabular<Float>::tabular(const brdf<Float>& brdf, uint32_t res, bool shadow)
:radial<Float>()
{
	DJB_ASSERT(res > 2 && "Invalid Resolution");
	//set_shadow(shadow);

	// allocate memory
	m_p22.reserve(res);
//	m_sigma.reserve(res);
//	m_cdf.reserve(res);
//	m_qf.reserve(res);

	// eval
	//compute_p22_smith(brdf, res);
	//normalize_p22();
	//compute_sigma();
	//compute_fresnel(brdf, res);

	// sample
//	compute_cdf();
//	compute_qf();
}*/

// Normalize the Slope PDF
/*
template <typename Float>
void tabular<Float>::normalize_p22()
{
	const int ntheta = 128;
	const float dphi = 2.0 * M_PI;
	const float dtheta = M_PI / (float)ntheta;
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
	Throw("djb_verbose: Slope PDF norm. constant = %f", nint);
#endif
}

/**
 * Compute the Smith Masking Term
 *
 * This requires the computation of an integral.
 * The resolution at which the integral is computed is hard coded to 
 * yield sufficient precision for most distributions.
 */ 
/*
template <typename Float>
void tabular<Float>::compute_sigma()
{
	const int ntheta = 90;
	const int nphi   = 180;
	float dtheta = M_PI / (float)ntheta;
	float dphi   = 2.0 * M_PI / (float)nphi;
	int cnt = m_p22.size() - 1;

	for (int i = 0; i < cnt; ++i) {
		float tmp = (float)i / (float)cnt; // in [0, 1)
		float theta_k = tmp * 0.5 * M_PI; // in [0, pi/2)
		float cos_theta_k = cos(theta_k);
		float sin_theta_k = sin(theta_k);
		float nint = 0.0;

		for (int j2 = 0; j2 < nphi; ++j2) {
			float u_j = (float)j2 / (float)nphi; // in [0, 1)
			float phi_h = u_j * 2.0 * M_PI;          // in [0, 2pi)
			for (int j1 = 0; j1 < ntheta; ++j1) {
				float u_i = (float)j1 / (float)ntheta; // in [0, 1)
				float theta_h = u_i * u_i * M_PI * 0.5; // in [0, sqrt(pi/2))
				float sin_theta_h = sin(theta_h);
				float kh = sin_theta_k * sin_theta_h * cos(phi_h)
				           + cos_theta_k * cos(theta_h);

				nint+= max((float)0.0, kh)
				     * ndf(vec3(theta_h, phi_h))
				     * u_i * sin_theta_h;
			}
		}
		nint*= dtheta * dphi;
		m_sigma.push_back(max((float)cos_theta_k, nint));
	}
	m_sigma.push_back(m_sigma.back());

#ifndef NVERBOSE
	Throw("djb_verbose: Projected area term ready");
#endif
}*/

// Fresnel computation 
/*
template <typename Float>
void tabular<Float>::compute_fresnel(const brdf<Float>& brdf, int res)
{
	std::vector<vec3> fresnel(res);
	int cnt = res - 1;

	// compute average ratio between input and Torrance Sparrow equation
	for (int i = 0; i < cnt; ++i) {
		const float_t phi_d = M_PI * 0.5; // this can NOT be tweaked
		const float_t phi_h = 0.0;        // this can be tweaked (no impact on MERL)
		float_t tmp = (float_t)i / (float_t)cnt;
		float_t theta_d = tmp * M_PI * 0.5; // linear parameterization in theta_d
		vec3 f = vec3(0); // Fresnel value at R, G, B wavelengths
		int count[3] = {0, 0, 0};
		float_t theta_h = 0.0;

		for (int j = 0; theta_h < M_PI * 0.5 - theta_d; ++j) {
			float_t tmp1 = (float_t)j / (float_t)cnt;
			theta_h = tmp1 * tmp1 * M_PI * 0.5; // in [0, pi/2)

			if (theta_h > M_PI * 0.5) continue;

			vec3 dir_h = vec3(theta_h, phi_h);
			vec3 dir_d = vec3(theta_d, phi_d);
			vec3 dir_i, dir_o;
			hd_to_io(dir_h, dir_d, &dir_i, &dir_o);

			dir_i = vec3(0, 0, 1); //XXX: hack to reproduce my EGSR fits
			vec3 fr1 = brdf.eval(dir_i, dir_o);
			vec3 fr2 = eval(dir_i, dir_o);

			if (fr2.x > /* epsilon *//*1e-4) {
				float_t ratio = fr1.x / fr2.x;
				f.x+= ratio;
				++count[0];
			}
			if (fr2.y > /* epsilon *//*1e-4) {
				float_t ratio = fr1.y / fr2.y;
				f.y+= ratio;
				++count[1];
			}
			if (fr2.z > /* epsilon *//*1e-4) {
				float_t ratio = fr1.z / fr2.z;
				f.z+= ratio;
				++count[2];
			}
		}

		// compute average
		fresnel[i].x = count[0] == 0 ? 1.0 : min((float_t)1.0, f.x / (float_t)count[0]);
		fresnel[i].y = count[1] == 0 ? 1.0 : min((float_t)1.0, f.y / (float_t)count[1]);
		fresnel[i].z = count[2] == 0 ? 1.0 : min((float_t)1.0, f.z / (float_t)count[2]);
	}
	// copy last value
	fresnel[res - 1] = fresnel[res - 2];
	set_fresnel(fresnel::spline(fresnel));
#ifndef NVERBOSE
	Throw("djb_verbose: Fresnel function ready");
#endif
}
*/
template <typename Float>
void tabular<Float>::compute_p22_smith(const brdf<Float>& brdf, int res)
{   MI_IMPORT_CORE_TYPES();
	int cnt = res - 1;
	float dtheta = sqrt(M_PI * 0.5) / (float)cnt;
	matrix<Float> km(cnt);

	for (int i = 0; i < cnt; ++i) {
		float tmp = (float)i / (float)cnt;
		float theta = tmp * sqrt(M_PI * 0.5);
		float theta_o = theta * theta;
		float cos_theta_o = cos(theta_o);
		float tan_theta_o = tan(theta_o);
        Float ss=dr::sin(theta_o);
        Vector3f xyz(ss*dr::cos(0.0),ss*dr::sin(0.0),dr::cos(theta_o));
		Vector3f fr = brdf.eval(xyz, xyz);
		Float fr_i =intensity(fr);
		Float kji_tmp = (dtheta * dr::pow(cos_theta_o, 6.0)) * (8.0 * fr_i);

		for (int j = 0; j < cnt; ++j) {
			const float dphi_h = M_PI / 180.0;
			float tmp = (float)j / (float)cnt;
			float theta = tmp * sqrt(M_PI * 0.5);
			float theta_h = theta * theta;
			float cos_theta_h = cos(theta_h);
			float tan_theta_h = tan(theta_h);
			float tan_product = tan_theta_h * tan_theta_o;
			float nint = 0.0;

			for (float phi_h = 0.0; phi_h < 2.0 * M_PI; phi_h+= dphi_h)
				nint+= std::max((float)1, tan_product * (float)cos(phi_h));
			nint*= dphi_h;

			km(j, i) = theta * kji_tmp * nint * tan_theta_h
			         / (cos_theta_h * cos_theta_h);
		}
	}

	// compute slope pdf
	const std::vector<Float> v = km.eigenvector(4);
	for (int i = 0; i < (int)v.size(); ++i)
		m_p22.push_back(1e-2 * v[i]);
	m_p22.push_back(0);
}
/*
// Compute the CDF (stored for debugging)
template <typename Float>
void tabular<Float>::compute_cdf()
{
	int cnt = (int)m_p22.size() - 1;
	float dtheta = M_PI / (float)cnt;
	float nint = 0.0;

	m_cdf.resize(0);
	for (int i = 0; i < cnt; ++i) {
		float u = (float)i / (float)cnt;
		float theta_h = u * u * M_PI * 0.5;
		float cos_theta_h = cos(theta_h);
		float r_h = tan(theta_h);
		float p22_r = p22_radial(r_h * r_h);

		nint+= (u * r_h * p22_r) / (cos_theta_h * cos_theta_h);
		m_cdf.push_back(nint * dtheta * /* normalize *//*(2.0 * M_PI));
	}
	m_cdf.push_back(1);

#ifndef NVERBOSE
	Throw("djb_verbose: Slope CDF ready");
#endif
}
*/
// -------------------------------------------------------------------------------------------------
// Compute the Quantile Function
/*
template <typename Float>
void tabular<Float>::compute_qf()
{
	int cnt = (int)m_p22.size() - 1;
	int res = cnt * 8; // resolution of inversion
	int j = 0;

	m_qf.resize(0);
	m_qf.push_back(0);
	for (int i = 1; i < cnt; ++i) {
		float cdf = (float)i / (float)cnt;

		for (; j < res; ++j) {
			float u = (float)j / (float)res;
			float theta_h = u * M_PI * 0.5;
			float qf = cdf_radial(tan(theta_h));

			// lerp lookup
			if (qf >= cdf) {
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
#endif
}*/

/* Constructors */
template <typename Float>
microfacet<Float>::params::params(Float a1, Float a2, Float phi_a)
{
	set_ellipse(a1, a2, phi_a);
	set_location(0, 0);
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
	return microfacet<Float>::params::elliptic(a, a, 0);
}

template <typename Float>
typename microfacet<Float>::params
microfacet<Float>::params::elliptic(Float a1, Float a2, Float phi_a)
{
   
	return microfacet<Float>::params(a1, a2, phi_a);
}
template <typename Float>
void microfacet<Float>::params::set_ellipse(Float a1, Float a2, Float phi_a)
{
	DJB_ASSERT(a1 > 0.0 && a2 > 0.0 && "Invalid ellipse radii");
	m_a1 = a1;
	m_a2 = a2;
	m_phi_a = phi_a;
	//ellipse_to_pdfparams<Float>(a1, a2, phi_a, &m_ax, &m_ay, &m_rho);
    Float cos_phi_a = dr::cos(phi_a);
	Float sin_phi_a = dr::sin(phi_a);
	Float cos_2phi_a = 2.0 * cos_phi_a * cos_phi_a - 1.0f; // cos(2x) = 2cos(x)^2 - 1
	Float a1_sqr = a1 * a1;
	Float a2_sqr = a2 * a2;
	Float tmp1 = a1_sqr + a2_sqr;
	Float tmp2 = a1_sqr - a2_sqr;

	Float m_ax = dr::sqrt(0.5 * (tmp1 + tmp2 * cos_2phi_a));
	Float m_ay = dr::sqrt(0.5 * (tmp1 - tmp2 * cos_2phi_a));
	Float m_rho = (a2_sqr - a1_sqr) * cos_phi_a * sin_phi_a / ((m_ax) * (m_ay));



	m_sqrt_one_minus_rho_sqr = dr::sqrt(1.0 - m_rho * m_rho);
}
template <typename Float>
void microfacet<Float>::params::set_location(Float tx_n, Float ty_n)
{
	m_tx_n = tx_n;
	m_ty_n = ty_n;
	m_n = dr::normalize(Vector3f(-tx_n, -ty_n, 1));
}


template <typename Float>
typename microfacet<Float>::params tabular<Float>::fit_ggx_parameters(const tabular<Float>& tab)
{
	const int ntheta = 128;
	float dtheta = M_PI / (float)ntheta;
	Float nint = 0.0;
	Float alpha=0.0;
/*
	for (int i = 0; i < ntheta; ++i) {
		float u = (float)i / (float)ntheta; // in [0,1)
		float theta_h = u * u * M_PI * 0.5;
		float cos_theta_h = cos(theta_h);
		float r_h = tan(theta_h);
		float r_h_sqr = r_h * r_h;
		//Float p22_r = tab.p22_radial(r_h_sqr);
        Float p22_r = 0.0f;

		nint=nint+ (u * r_h_sqr * p22_r) / (cos_theta_h * cos_theta_h);
	}
	nint=nint* dtheta * /* int_0^2pi fabs(cos(phi)) dphi *//*4.0;
	alpha = nint;*/
/*
#ifndef NVERBOSE
	Throw("Invalid file structure: %f", alpha);
#endif*/
  // return typename microfacet<Float>::params::params(alpha, alpha, alpha);
	return microfacet<Float>::params::isotropic(alpha);
}
template <typename Float>
Float tabular<Float>::fit_ggx_parameter(const tabular<Float> &tab)
{  MI_IMPORT_CORE_TYPES();
	const int ntheta = 128;
	double dtheta = M_PI / (double)ntheta;
	double nint = 0.0;
	Float alpha=0.0;
  //  float p22_r=0.0;
    
	for (int i = 0; i < ntheta; ++i) {
		double u = i / ntheta; // in [0,1)
		double theta_h = u * u * M_PI * 0.5;
		double cos_theta_h = cos(theta_h);
		double r_h = tan(theta_h);
		double r_h_sqr = r_h * r_h;
	//	Float p22_r = tab.p22_radial(r_h_sqr); //this causes the error
        double p22_r=0.0;
       

		nint=nint+ (u * r_h_sqr * p22_r) / dr::maximum((cos_theta_h * cos_theta_h),DJB_EPSILON);
	}
    
	nint=nint* dtheta *4.0;// int_0^2pi fabs(cos(phi)) dphi 
//	alpha = nint;


   return alpha;
	//return microfacet<Float>::params::isotropic(alpha);
}


template <typename Float>
Float ggx<Float>::sigma_std_radial(Float cos_theta_k)
{
	return ((1.0 + cos_theta_k) / 2.0);
}












namespace spline {


template <typename Float>
Float eval(const std::vector<Float>& points, std::string wrap, float u)
{   MI_IMPORT_CORE_TYPES();
	int edge   = (int)points.size();
    using FloatStorage                = DynamicBuffer<Float>;
	FloatStorage s_data;

  std::unique_ptr<Float[]> data_out(new Float[edge]);
	Float *data_out_ptr = data_out.get();
	Float *points_ptr = points.data();
	 for (uint32_t k = 0; k < edge; ++k)
        *data_out_ptr++ = *points_ptr++ ;
	s_data = dr::load<FloatStorage>(data_out.get(), (edge));
    
    UInt32 intpart=UInt32 ((u * edge)-u);
    Float frac=u * edge -u - dr::floor((u * edge)-u); 

    edge=(UInt32) edge;
    UInt32 i1 = dr::clamp(intpart, 0,edge-1);
	UInt32 i2 = dr::clamp(intpart + 1,0, edge-1);
    if (wrap == "edge") {
	UInt32 i1 = dr::clamp(intpart, 0,edge-1);
	UInt32 i2 = dr::clamp(intpart + 1,0, edge-1);
   }
   else if (wrap == "repeat") {

    while (i1 >= edge) i1-= edge;
	while (i1 < 0    ) i1+= edge;
    while (i2 >= edge) i2-= edge;
    while (i2 < 0    ) i2+= edge;
   }

	const Float &p1 = dr::gather<Float>(s_data, i1);
	const Float &p2 = dr::gather<Float>(s_data, i2);

	return dr::lerp(p1, p2, frac);
}
}


NAMESPACE_END(mitsuba)

#endif // DJ_BRDF_IMPLEMENTATION
