#include <limits>
#include <fstream>

#include "rtraytracer.h"
#include "rtsurface.h"
#include "rtvector.h"

using namespace std;

ostream& operator<<(ostream& stream,const rt_ray_field& obj)
{
	stream << obj.pol << ' ';
	stream << obj.from.x << ' ' << obj.from.y << ' ' << obj.from.z << ' ' << obj.intensity_i << ' ';
	stream << obj.to.x << ' ' << obj.to.y << ' ' << obj.to.z << ' ' << obj.intensity_f << ' ';
	stream << obj.abs_vol << ' ' << obj.alpha << ' ' << obj.abs_film << ' ';

	return stream;
}

ostream& operator<<(ostream& stream,const rt_ray_data& obj)
{
	for(int i=0; i<obj.ray_que.size(); i++) {
		stream << *(obj.ray_que[i]) << endl;
	}
	return stream;
}

rt_raytracer::rt_raytracer(double min_int, double min_sep, double lambda)
{
	intensity_min=min_int;
	separation_min=min_sep;
	wavelength=lambda;
	record_trace=false;	
	b_trace_counter=false;
}

rt_raytracer::~rt_raytracer(void)
{
	surface_map::iterator pIter;
	for(pIter=surfaces.begin(); pIter!=surfaces.end(); pIter++) {
		delete(pIter->second);
	}
	surfaces.clear();
}

bool rt_raytracer::add_surface(rt_surface* s)
{
	surface_map::const_iterator iter=surfaces.find(s);
	if(iter==surfaces.end()) { // if not found
		double* qe=new (nothrow) double[n_info];
		if(qe==NULL) {
			return false;
		}
		qe[0]=0.0; qe[1]=0.0;
		pair< surface_map::iterator, bool > pr=surfaces.insert(surface_pair(s,qe));
		return pr.second;
	} 
	return false;
}

bool rt_raytracer::delete_surfaces(void)
{
	surface_map::iterator pIter;
	for(pIter=surfaces.begin(); pIter!=surfaces.end(); pIter++) {
		delete pIter->first;
		delete pIter->second;
	}
	surfaces.clear();
	return true;
}

bool rt_raytracer::add_efficiencies(rt_surface* s, const unsigned int index_abs, const double abs)
{
	surface_map::const_iterator iter=surfaces.find(s);
	if(iter!=surfaces.end()) { // if found
		iter->second[index_abs]+=abs; // update absorption
		return true;
	} 
	return false;
}

bool rt_raytracer::clear_efficiencies(void)
{
	surface_map::iterator pIter;
	for(pIter=surfaces.begin(); pIter!=surfaces.end(); pIter++) {
		pIter->second[i_abs]=0;
		pIter->second[i_vol]=0;
	}
	return true;
}

bool rt_raytracer::scale_efficiencies(const double scale)
{
	surface_map::iterator pIter;
	for(pIter=surfaces.begin(); pIter!=surfaces.end(); pIter++) {
		pIter->second[i_abs] /= scale;
		pIter->second[i_vol] /= scale;
	}
	return true;
}

bool rt_raytracer::get_efficiencies(rt_surface* s, double& abs, double& absvol)
{
	surface_map::const_iterator iter=surfaces.find(s);
	if(iter!=surfaces.end()) { // if found
		abs=iter->second[i_abs]; // update absorption
		absvol=iter->second[i_vol];
		return true;
	} 
	return false;
}

int rt_raytracer::get_total_efficiencies(double& aqe, double& abssurf, double& absvol)
{
	aqe=0; abssurf=0; absvol=0;
	int nactsurface=0;
	surface_map::iterator pIter;
	for(pIter=surfaces.begin(); pIter!=surfaces.end(); pIter++) {
		rt_surface* psurface=pIter->first;
		double abs_surf, abs_vol;
		if(get_efficiencies(psurface,abs_surf,abs_vol)) {
			aqe+=abs_surf+abs_vol;
			abssurf+=abs_surf;
			absvol+=abs_vol;
			nactsurface++;
		}
	}
	return nactsurface;
}

void rt_raytracer::set_wavelength(const double lambda)
{
	wavelength=lambda;
	surface_map::iterator pIter;
	for(pIter=surfaces.begin(); pIter!=surfaces.end(); pIter++) {
		rt_surface* psurface=pIter->first;
		psurface->set_wavelength(lambda);
	}
}

void rt_raytracer::set_material_library(rt_nk_library* plib)
{
	pmaterial_library=plib;
	surface_map::iterator pIter;
	for(pIter=surfaces.begin(); pIter!=surfaces.end(); pIter++) {
		rt_surface* psurface=pIter->first;
		psurface->set_material_library(plib);
	}
}

bool rt_raytracer::trace(const rt_vector& from, const rt_vector& dir_in, const string& material_in, const int pol, const double intensity_in,
						 double &aqe, double& aqesurf, double& aqevol)
{
	n_trace_counter=0;
	clear_efficiencies();
	trace_que.push_back(new rt_trace_param(from,dir_in,material_in,pol,intensity_in));
	while(trace_que.size()>0) {
		if(b_trace_counter && (int)trace_que.size()>n_trace_max) {
			clear_trace_que();
			break;
		}
		rt_trace_param* pparam=(rt_trace_param*)trace_que.back();
		trace_que.pop_back();
		_trace(pparam->from,pparam->dir_in,pparam->material_in,pparam->pol,pparam->intensity_in);
		delete pparam;
	}
	get_total_efficiencies(aqe,aqesurf,aqevol);
	return true;
}


bool rt_raytracer::trace_bunch(const rt_vector& from_topleft, const rt_vector& dir1, const rt_vector& dir2, const int n1, const int n2,
		const rt_vector& dir_in, const string& material_in, const int pol,const double intensity_in, double& aqe, double& aqesurf, double& aqevol)
{
	n_trace_counter=0;
	clear_efficiencies();
	rt_vector from=from_topleft;
	for(int i1=0; i1<n1; i1++) {
		for(int i2=0; i2<n2; i2++) {
			from=from_topleft+i1*dir1+i2*dir2;
			trace_que.push_back(new rt_trace_param(from,dir_in,material_in,pol,intensity_in));
			while(trace_que.size()>0) {
				if(b_trace_counter && (int)trace_que.size()>n_trace_max) {
					clear_trace_que();
					break;
				}
				rt_trace_param* pparam=(rt_trace_param*)trace_que.back();
				trace_que.pop_back();
				_trace(pparam->from,pparam->dir_in,pparam->material_in,pparam->pol,pparam->intensity_in);
				delete pparam;
			}
		}
	}
	scale_efficiencies((double)n1*n2);
	get_total_efficiencies(aqe,aqesurf,aqevol);

	return true;
}

bool rt_raytracer::_trace(const rt_vector& from, const rt_vector& dir_in, const string& material_in, const int pol, const double intensity_in)
{
	surface_map::iterator pIter, pMin;
	bool found=false;
	double distance, distance_min=numeric_limits<double>::max();
	rt_vector to_min, normal_min;
	rt_surface* surf_min;

	// find cross points
	for(pIter=surfaces.begin(); pIter!=surfaces.end(); pIter++) {
		rt_vector to, normal;
		rt_surface* psurface=pIter->first;
		if(psurface->find_cross_point(from, dir_in, to, normal)) {
			distance=rt_vector::norm(from-to);
			if(distance<distance_min && distance>separation_min) {
				pMin=pIter;
				surf_min=psurface;
				distance_min=distance;
				to_min=to;
				normal_min=normal;
				found=true;
			}
		}
	}
	// find reflected and refracted rays
	bool generate_ray=false;
	if(found) {
		rt_vector dir_r, dir_t;
		double intensity_r, intensity_t;
		string material_to;
		unsigned int flag=0;
		double k, abs_vol=0.0;
		double intensity_to=intensity_in;
		complex<double> n_in;
		pmaterial_library->find(material_in.data(), wavelength, n_in);
		// absorption by volume
		if((k=n_in.imag())>0) {
			const double pi=3.141592;
			abs_vol=intensity_in*(1-exp(-4*pi*k/wavelength*rt_vector::norm(to_min-from))); // Beer Lambert's law
			intensity_to=intensity_in-abs_vol;
			add_efficiencies(surf_min,i_vol,abs_vol); // change to tell if it's volume absorption or surface absorption
		}
		flag=surf_min->reflect_refract(dir_in,normal_min,intensity_to,pol,material_in,dir_r,intensity_r,dir_t,intensity_t,material_to);
		// absorption by thin films
		double abs=intensity_to-intensity_r-intensity_t;
		if(abs>0) add_efficiencies(surf_min,i_abs,abs); // need to implement absorption between surfaces
		if(record_trace) {
			ray_data.push(pol,from,intensity_in,to_min,intensity_to,abs_vol,4*3.141592*k/wavelength,(abs>0)?abs:0);
		}
		if(flag) {
			if(intensity_r>intensity_min) {
				trace_que.push_back(new rt_trace_param(to_min,dir_r,material_in,pol,intensity_r));
				generate_ray=true;
			}
			if(intensity_t>intensity_min) {
				trace_que.push_back(new rt_trace_param(to_min,dir_t,material_to,pol,intensity_t));
				generate_ray=true;
			}
		}
	}

	return (found && generate_ray);
}
