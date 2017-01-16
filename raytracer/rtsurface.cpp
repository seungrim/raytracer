#include <limits>
#include <cstdlib>
#include "rtsurface.h"
#include "rtvector.h"

using namespace std;

// find_cross_point
// find the nearest t out of all shapes
bool rt_surface::find_cross_point(const rt_vector& from, const rt_vector& direction, rt_vector& to, rt_vector& normal) const {
	double dist_min=numeric_limits<float>::max(), dist;
	double epsilon=numeric_limits<float>::epsilon();
	int degeneracy=0;
	for(int j=0; j<(int)shapes.size(); j++) {
		rt_vector toj, normalj;
		if(shapes[j]->find_cross_point(from,direction,toj,normalj)) {
			dist=rt_vector::norm(from-toj);
			if(dist>epsilon&&dist<dist_min) {
				dist_min=dist;
				to=toj;
				normal=normalj;
				degeneracy=1;
			} else if(dist==dist_min) {
				normal+=normalj;
				degeneracy++;
			}
		}
	}
	if(degeneracy>=1) {
		normal=normal/rt_vector::norm(normal);
		return true;
	}
	return false;
}

void rt_surface::erase_shape(rt_shape* s)
{
	vector<rt_shape*>::iterator iter;
	for(iter=shapes.begin(); iter!=shapes.end(); iter++) {
		if(s==(rt_shape*)*iter) {
			shapes.erase(iter);
		}
	}
}

int rt_fresnel_surface::reflect_refract(const rt_vector& dir_in,const rt_vector& normal, const double intensity_in, const int pol, const string& material_in,
		rt_vector& dir_r, double& intensity_r, rt_vector& dir_t, double& intensity_t, string& material_to) const
{
	complex<double> n_in, n_t;
	pmaterial_library->find(material_in.data(),wavelength,n_in);

	// fresnel refraction and reflection here
	complex<double> ref[2];
	pmaterial_library->find(materials[0].data(), wavelength, ref[0]);
	pmaterial_library->find(materials[1].data(), wavelength, ref[1]);

	if(n_in==ref[0]) {
		n_t=ref[1];
		material_to=materials[1];
	} else { 
		n_t=ref[0];
		material_to=materials[0];
	}
	// direction reflected
	double cos_i=-rt_vector::dot_product(dir_in,normal);
	dir_r=dir_in+normal*cos_i*2.0;
	dir_r=dir_r/rt_vector::norm(dir_r);
	intensity_r=intensity_in;
	// direction refracted
	double nr=n_in.real()/n_t.real();
	double cos_t, sin_t_sq;
	sin_t_sq=nr*nr*(1-cos_i*cos_i);
	if(sin_t_sq>1) { // TIR
		dir_t.set(0,0,0);
	} else {
		cos_t=sqrt(1-sin_t_sq);
		dir_t=dir_in*nr+normal*(nr*cos_i-cos_t);
	}
	// intensity reflected
	complex<double> qi, qt, r, t;
	qi=n_in*rt_vector::dot_product(dir_r,normal);
	qt=-n_t*rt_vector::dot_product(dir_t,normal);
	if(pol==s_pol) { // s-polarization
		r=(qi-qt)/(qi+qt);
	} else { // p-polarization
		complex<double> denom=n_t*n_t*qi+n_in*n_in*qt;
		if(abs(denom)==0) { // brewster angle
			r=0;
		} else {
			r=(n_t*n_t*qi-n_in*n_in*qt)/denom;
		}
	}
	intensity_r=abs(r)*abs(r)*intensity_in;
	// intensity refracted
	if(sin_t_sq<=1) {
		if(pol==s_pol) {
			t=2.0*qi/(qi+qt);
		} else {
			t=2.0*n_in*n_t*qi/(n_t*n_t*qi+n_in*n_in*qt);
		}
		intensity_t=abs(t)*abs(t)*qt.real()/qi.real()*intensity_in;
	} else { // TIR
		intensity_t=0;
	}

	return reflected_refracted;
}

int rt_mirror_surface::reflect_refract(const rt_vector& dir_in,const rt_vector& normal, const double intensity_in, const int pol, const string& material_in,
		rt_vector& dir_r, double& intensity_r, rt_vector& dir_t, double& intensity_t, string& material_to) const
{
	material_to=material_in;
	// mirror reflection here
	double c1=-rt_vector::dot_product(dir_in,normal);
	dir_r=dir_in+normal*c1*2.0;
	dir_r=dir_r/rt_vector::norm(dir_r);
	intensity_r=intensity_in;
	dir_t.set(0,0,0); intensity_t=0; // no refraction

	return reflected;
}

int rt_diffuse_surface::reflect_refract(const rt_vector& dir_in,const rt_vector& normal, const double intensity_in, const int pol, const string& material_in,
		rt_vector& dir_r, double& intensity_r, rt_vector& dir_t, double& intensity_t, string& material_to) const
{
	// diffuse reflection here
	rt_vector ex(1,0,0), ey(0,1,0), ez(0,0,1);
	rt_vector e1, e2;
	if(rt_vector::norm(rt_vector::cross_product(normal,ex))) {
		e1=ex;
	} else if(rt_vector::norm(rt_vector::cross_product(normal,ey))) {
		e1=ey;
	} else if(rt_vector::norm(rt_vector::cross_product(normal,ez))) {
		e1=ez;
	}
	e1=e1-rt_vector::dot_product(e1,normal)*normal;
	e1=e1/rt_vector::norm(e1);
	e2=rt_vector::cross_product(normal,e1);
	e2=e2/rt_vector::norm(e2);
	const double pi=3.141592;
	double phi = rand()/(double)RAND_MAX * 2.0 * pi;
	double theta = rand()/(double)RAND_MAX * pi/2.0;
	dir_r=cos(theta)*cos(phi)*e1+cos(theta)*sin(phi)*e2+sin(theta)*normal;
	intensity_r=intensity_in;
	dir_t.set(0,0,0); intensity_t=0; 
	return reflected;
}

int rt_multilayer_surface::reflect_refract(const rt_vector& dir_in,const rt_vector& normal, const double intensity_in, const int pol, const string& material_in,
		rt_vector& dir_r, double& intensity_r, rt_vector& dir_t, double& intensity_t, string& material_to) const
{
	int return_flag=0;
	
	rt_fresnel_surface::reflect_refract(dir_in,normal,intensity_in,pol,material_in,dir_r,intensity_r,dir_t,intensity_t,material_to);

	double angle_in=acos(-rt_vector::dot_product(dir_in,normal)); // dir_in and normal are assumed to be unit vectors
	double aqe=0.0, eqe=0.0, trans=0.0;
	if(material_in.compare(pmultilayer->get_material_front())==0) {
		pmultilayer->calc_efficiencies(angle_in,pol,aqe,trans);
	} else { // reverse order
		pmultilayer->calc_efficiencies_reverse(angle_in,pol,aqe,trans);
	}
	intensity_r=intensity_in*(1-aqe-trans);
	intensity_t=intensity_in*trans;
		
	if(intensity_r>0) return_flag |= reflected;
	if(intensity_t>0) return_flag |= refracted;

	return return_flag;
}
