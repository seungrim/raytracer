////////////////////////////////////////////////////////////////////////////////////////////////
// RayTracing Surface Classes
//
// Mar 25, 2008 Seung Rim
//
// Surfaces are basic elements for raytracing
// A surface has a set of shapes with find_cross_point function and reflect_refract function 
// find_cross_point() - find the nearest cross point between a surface and a ray
// reflect_refract() - find reflected and refracted rays
// classes implemented here
//	rt_surface - abstract class
//		rt_fresnel_surface:rt_surface - surface follows fresnel's laws
//			rt_absorb_surface:rt_fresnel_surface - absorbs fraction of rays and follow fresnel's laws
//			rt_multilayer_surface::rt_fresnel_surface - surface with thin film multilayer structures
//				multilayers are calculated mlcalc.cpp based on transfer matrix method
//		rt_mirror_surface:rt_surface - perfect mirror to reflect all rays
//			rt_diffuse_surface:rt_mirror_surface - reflect rays isotropically
//		rt_escape_surface:rt_surface - all rays are escaped through this surface i.e. rays disappear
// 
#pragma once

#ifndef _RT_SURFACE_H_
#define _RT_SURFACE_H_

#include <vector>
#include <complex>
#include "rtvector.h"
#include "rtshape.h"
#include "rtmultilayer.h"
#include "rtnklibrary.h"

using namespace std;

class rt_surface {
public:
	// constructors and destructor
	rt_surface() { }
	rt_surface(const rt_surface& s) { set(s.shapes); }
	rt_surface(const vector<rt_shape*>& s) { set(s); }
	~rt_surface(void) { shapes.clear(); }
	void delete_shapes(void) { for(int j=0; j<(int)shapes.size(); j++) { delete shapes[j]; shapes[j]=(rt_shape*)0; } shapes.clear(); }
	// getters and setters
	void set(const vector<rt_shape*>& s) { shapes=s; }
	virtual void set_wavelength(const double lambda) { wavelength=lambda; }
	virtual void set_material_library(rt_nk_library* plib) { pmaterial_library=plib; }
	unsigned int add_shape(rt_shape* s) { shapes.push_back(s); return shapes.size(); }
	void erase_shape(rt_shape* s);
	bool find_cross_point(const rt_vector& from, const rt_vector& direction, rt_vector& to, rt_vector& normal) const;
	// pure virtual functions
	static const int s_pol=0, p_pol=1;
	static const int reflected=1, refracted=2, reflected_refracted=(reflected|refracted), none=0;
	virtual int reflect_refract(const rt_vector& dir_in,const rt_vector& normal, const double intensity_in, const int pol, const string& material_in,
		rt_vector& dir_r, double& intensity_r, rt_vector& dir_t, double& intensity_t, string& material_to) const=0;
protected:
	double wavelength;
	vector<rt_shape*> shapes;
	rt_nk_library* pmaterial_library;
};

class rt_fresnel_surface: public rt_surface {
public:
	// constructors and destructor
	rt_fresnel_surface() {}
	rt_fresnel_surface(const rt_fresnel_surface& s) { rt_surface::set(s.shapes); set(s.materials[0],s.materials[1]); }
	rt_fresnel_surface(const string& m1, const string& m2) { set(m1,m2); }
	~rt_fresnel_surface(void) {}
	// getters and setters
	void set(const string& m1, const string& m2) { 
		materials[0]=m1; materials[1]=m2;
	}
	virtual int reflect_refract(const rt_vector& dir_in,const rt_vector& normal, const double intensity_in, const int pol, const string& material_in,
		rt_vector& dir_r, double& intensity_r, rt_vector& dir_t, double& intensity_t, string& material_to) const;
protected:
	string materials[2];
};

class rt_mirror_surface: public rt_surface {
public:
	rt_mirror_surface() {}
	rt_mirror_surface(const rt_mirror_surface& s) { *this=s; }

	virtual int reflect_refract(const rt_vector& dir_in,const rt_vector& normal, const double intensity_in, const int pol, const string& material_in,
		rt_vector& dir_r, double& intensity_r, rt_vector& dir_t, double& intensity_t, string& material_to) const;
};

class rt_diffuse_surface: public rt_mirror_surface {
public:
	rt_diffuse_surface() {}
	rt_diffuse_surface(const rt_diffuse_surface& s) { *this=s; }

	virtual int reflect_refract(const rt_vector& dir_in,const rt_vector& normal, const double intensity_in, const int pol, const string& material_in,
		rt_vector& dir_r, double& intensity_r, rt_vector& dir_t, double& intensity_t, string& material_to) const;
};

class rt_escape_surface: public rt_surface {
public:
	rt_escape_surface() {}
	rt_escape_surface(const rt_escape_surface& s) { *this=s; }

	int reflect_refract(const rt_vector& dir_in,const rt_vector& normal, const double intensity_in, const int pol, const string& material_in,
		rt_vector& dir_r, double& intensity_r, rt_vector& dir_t, double& intensity_t, string& material_to) const
	{ dir_r.set(0,0,0); intensity_r=0; dir_t.set(0,0,0); intensity_t=0; material_to=material_in; return none; }
};

class rt_absorb_surface: public rt_fresnel_surface {
public:
	rt_absorb_surface(const double fabs, const string& m1, const string& m2):rt_fresnel_surface(m1,m2) {f_abs=fabs;}
	rt_absorb_surface(const rt_absorb_surface& s) { *this=s; }

	virtual int reflect_refract(const rt_vector& dir_in,const rt_vector& normal, const double intensity_in, const int pol, const string& material_in,
		rt_vector& dir_r, double& intensity_r, rt_vector& dir_t, double& intensity_t, string& material_to) const
	{ return rt_fresnel_surface::reflect_refract(dir_in,normal,intensity_in*(1-f_abs),pol,material_in,dir_r,intensity_r,dir_t,intensity_t,material_to);}
private:
	double f_abs;
};

class rt_multilayer_surface: public rt_fresnel_surface {
public:
	rt_multilayer_surface(rt_multilayer* player, const string& m1, const string& m2):rt_fresnel_surface(m1,m2) { 
		set_multilayer(player); 
	}
	rt_multilayer_surface(const rt_multilayer_surface& s) { *this=s; }
	~rt_multilayer_surface(void) { delete pmultilayer; }

	// getters and setters
	void set_multilayer(rt_multilayer* player) { pmultilayer=player; }
	void set_material_library(rt_nk_library* plib) { rt_fresnel_surface::set_material_library(plib); pmultilayer->set_material_library(plib); }
	void set_wavelength(const double lambda) { rt_surface::set_wavelength(lambda); pmultilayer->set_wavelength(lambda); }

	int reflect_refract(const rt_vector& dir_in,const rt_vector& normal, const double intensity_in, const int pol, const string& material_in,
		rt_vector& dir_r, double& intensity_r, rt_vector& dir_t, double& intensity_t, string& material_to) const;

private:
	rt_multilayer* pmultilayer;
};

#endif // _RT_SURFACE_H_

