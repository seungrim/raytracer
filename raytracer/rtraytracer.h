//////////////////////////////////////////////////////////////////////////
// RayTracing RayTracer class
//
// Mar 25, 2008 Seung Rim
//
// rt_raytracer has surfaces where ray traverses.
// After trace or trace_bunch called, absorption and currents are recorded and
// they can be read before calling trace or trace_bunch again.
//
#pragma once

#ifndef _RT_RAYTRACER_H_
#define _RT_RAYTRACER_H_

#include <iostream>
#include <map>
#include <vector>
#include "rtvector.h"
#include "rtsurface.h"
#include "rtnklibrary.h"

using namespace std;

class rt_trace_param {
public:
	rt_trace_param(rt_vector f, rt_vector d, string m, int p, double i) {
		from=f; dir_in=d; material_in=m; pol=p; intensity_in=i;
	}	
	rt_trace_param(const rt_trace_param& t) {
		from=t.from; dir_in=t.dir_in; material_in=t.material_in;
		pol=t.pol; intensity_in=t.intensity_in;
	}
	
public:
	rt_vector from;
	rt_vector dir_in;
	string material_in;
	int pol;
	double intensity_in;
};

class rt_ray_field {
public:
	rt_ray_field(int p,rt_vector f, double int_i, rt_vector t, double int_f, double abs_v, double a, double abs_f) {
		pol=p; from=f; intensity_i=int_i; to=t; intensity_f=int_f; abs_vol=abs_v; alpha=a; abs_film=abs_f;
	}
	rt_ray_field(const rt_ray_field& t) {
		pol=t.pol;
		from=t.from; intensity_i=t.intensity_i; to=t.to; intensity_f=t.intensity_f; 
		abs_vol=t.abs_vol; alpha=t.alpha; abs_film=t.abs_film;
	}
	friend ostream& operator <<(ostream& stream,const rt_ray_field& obj);
public:
	int pol; // 0:s- 1:p-
	rt_vector from, to;
	double intensity_i, intensity_f, abs_vol, abs_film, alpha;
};

class rt_ray_data {
public:
	rt_ray_data() { ray_que.clear(); }
	~rt_ray_data() { clear(); }
	void clear() {
		while(!ray_que.empty()) {
			rt_ray_field *pfield=ray_que.back();
			delete pfield;
			ray_que.pop_back();
		}
	}
	void push(int p,rt_vector f, double int_i, rt_vector t, double int_f, double abs_v, double a, double abs_f) {
		rt_ray_field *pfield=new rt_ray_field(p,f,int_i,t,int_f,abs_v,a,abs_f);
		ray_que.push_back(pfield);
	}
	friend ostream& operator <<(ostream& stream, const rt_ray_data& obj);

private:
	vector<rt_ray_field*> ray_que;
};

class rt_raytracer
{
public:
	typedef map<rt_surface*, double* > surface_map;
	typedef pair<rt_surface*, double* > surface_pair;
	
	rt_raytracer(double min_int=0,double min_sep=0,double lambda=600.0);
	rt_raytracer(const rt_raytracer& t) { surfaces=t.surfaces; intensity_min=t.intensity_min; separation_min=t.separation_min; 
		record_trace=t.record_trace; record_flag=t.record_flag; precstream=t.precstream; wavelength=t.wavelength;
		b_trace_counter=t.b_trace_counter; n_trace_counter=t.n_trace_counter; n_trace_max=t.n_trace_max; }
	~rt_raytracer(void);

	bool add_surface(rt_surface* s);
	bool delete_surfaces(void);
	const static int s_pol=0, p_pol=1;
	bool trace(const rt_vector& from, const rt_vector& dir_in, const string& material_in, const int pol, const double intensity_in,
		double &aqe,double& aqesurf, double& aqevol);
	bool trace_bunch(const rt_vector& from_topleft, const rt_vector& dir1, const rt_vector& dir2, const int n1, const int n2,
		const rt_vector& dir_in, const string& material_in, const int pol,const double intensity_in,double& aqe,double& aqesurf, double& aqevol);
	void set_seperation_min(const double s) { separation_min=s; }
	void set_wavelength(const double l);
	void set_material_library(rt_nk_library* plib);

	void record_begin() { record_trace=true; }
	void record_end(void) { record_trace=false; }
	void record_ray_data(ostream* pos) { *pos << ray_data; }
	void discard_ray_data(void) { ray_data.clear(); }
	
	void set_trace_counter(int max_counter) { b_trace_counter=true; n_trace_max=max_counter; }
	void clear_trace_counter(void) { b_trace_counter=false; }
	void clear_trace_que(void) {
		while(trace_que.size()>0) {
			rt_trace_param* pparam=(rt_trace_param*)trace_que.back();
			trace_que.pop_back();
			if(pparam!=NULL) delete pparam;
		}
		trace_que.clear();
	}

private:
	bool _trace(const rt_vector& from, const rt_vector& dir_in, const string& material_in, const int pol, const double intensity_in);

	bool clear_efficiencies(void);
	bool scale_efficiencies(const double scale);
	bool get_efficiencies(rt_surface* s, double& abs, double& absvol);
	int get_total_efficiencies(double& abs, double& abssurf, double& absvol);
	bool add_efficiencies(rt_surface* s, const unsigned int index_abs, const double abs);

	const static unsigned int i_abs=0, i_vol=1, n_info=2; // index of efficiencies
	surface_map surfaces;
	double intensity_min, separation_min;
	double wavelength;

	bool record_trace;
	unsigned int record_flag;	
	ostream* precstream;

	rt_nk_library* pmaterial_library;

	bool b_trace_counter;
	int n_trace_counter, n_trace_max; // counting number of recursion // for debugging purpose currently

	vector<rt_trace_param*> trace_que;

	rt_ray_data ray_data;
};

#endif //_RT_RAYTRACER_H_
