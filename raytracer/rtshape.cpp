#include <limits>

#include "rtshape.h"
#include "fminbnd.h"

using namespace std;

bool rt_2d_shape::find_cross_point(const rt_vector& from, const rt_vector& direction, rt_vector& to, rt_vector& n) const {
	double epsilon=numeric_limits<float>::epsilon();
	double t_den=rt_vector::dot_product(direction,normal);
	if(t_den!=0) {
		double t;
		t=rt_vector::dot_product(get_one_vector()-from,normal)/t_den;
		to=from+t*direction;
		if(t>epsilon&&inside(to)) {
			if(t_den>0) n=-normal; else n=normal;
			if(rt_vector::norm(n)==0) {
				return false;
			}
			n=n/rt_vector::norm(n);
			return true;
		}
	}
	return false;
}

bool rt_triangle::inside(const rt_vector& v) const {
	double sa,sb,sc;
	sa=rt_vector::dot_product(rt_vector::cross_product(vertex[1]-vertex[0],v-vertex[0]),normal);
	sb=rt_vector::dot_product(rt_vector::cross_product(vertex[2]-vertex[1],v-vertex[1]),normal);
	sc=rt_vector::dot_product(rt_vector::cross_product(vertex[0]-vertex[2],v-vertex[2]),normal);
	if(sa>=0 && sb>=0 && sc>=0) {
		return true;
	}
	return false;
}

bool rt_triangles::find_cross_point(const rt_vector& from, const rt_vector& direction, rt_vector& to, rt_vector& n) const
{
	double dist_min=numeric_limits<float>::max(), epsilon=numeric_limits<float>::epsilon();
	double dist;
	int degeneracy=0;
	for(int j=0; j<(int)triangles.size(); j++) {
		rt_vector toj, normalj;
		if(triangles[j]->find_cross_point(from,direction,toj,normalj)) {
			dist=rt_vector::norm(from-toj);
			if(dist>epsilon&&dist<dist_min) {
				dist_min=dist;
				to=toj;
				n=normalj;
				degeneracy=1;
			} else if(dist==dist_min) {
				n=n+normalj;
				degeneracy++;
			}
		}
	}
	if(degeneracy>=1) {
		if(rt_vector::norm(n)==0) {
			return false;
		}
		n=n/rt_vector::norm(n);
		return true;
	}
	return false;
}

double distance_point_to_curve(const double t, const void* fopts1, const void* fopts2) {
	rt_vector *from=(rt_vector*)fopts1;
	rt_curve *curve=(rt_curve*)fopts2;
	return rt_vector::norm(*from-curve->get_point_at(t));
}

double rt_curve::find_nearest_point(const rt_vector& from, rt_vector& on_curve, double &t_curve) const {
	double fval;
	t_curve=fminbnd(tmin,tmax,fval,distance_point_to_curve,(const void*)&from,(const void*)this);
	on_curve=get_point_at(t_curve);
	return rt_vector::norm(from-on_curve);
}

double distance_ray_to_curve(const double t, const void* fopts1, const void* fopts2) {
	rt_vector **vs=(rt_vector**)fopts1;
	rt_vector from=*vs[0];
	rt_vector direction=*vs[1];
	rt_tube *tube=(rt_tube*)fopts2;
	rt_vector on_curve;
	double t_curve;

	return abs(tube->get_curve().find_nearest_point(from+t*direction,on_curve,t_curve));
}

double distance_ray_to_tube(const double t, const void* fopts1, const void* fopts2) {
	rt_vector **vs=(rt_vector**)fopts1;
	rt_vector from=*vs[0];
	rt_vector direction=*vs[1];
	rt_tube *tube=(rt_tube*)fopts2;
	rt_vector on_curve;
	double t_curve;

	return abs(tube->get_curve().find_nearest_point(from+t*direction,on_curve,t_curve)-tube->get_radius());
}

bool rt_tube::find_cross_point(const rt_vector& from, const rt_vector& direction, rt_vector& to, rt_vector& normal) const {
	double tmin,tmax; get_curve().gett(tmin,tmax); tmax=1e2*(tmax-tmin);
	double epsilon=numeric_limits<float>::epsilon();
	epsilon=(epsilon>tol)?epsilon:tol;
	rt_vector* opts[]={(rt_vector*)&from,(rt_vector*)&direction};
	double fval=numeric_limits<float>::max();
	tmax=fminbnd(epsilon,tmax,fval,distance_ray_to_curve,(const void*)opts,(const void*)this);
	if(fval>get_radius()) {
		return false;
	}
	double t=fminbnd(epsilon,tmax,fval,distance_ray_to_tube,(const void*)opts,(const void*)this);
	if(t>epsilon && abs(fval)<tolfun) {
		to=from+t*direction;
		rt_vector on_curve; 
		double u_curve;
		get_curve().find_nearest_point(to,on_curve,u_curve);
		rt_vector n=to-on_curve;
		if(rt_vector::dot_product(direction,n)>0) n=-n;
		if(rt_vector::norm(n)==0) {
			return false;
		}
		normal=n/rt_vector::norm(n);
		return true;
	}
	return false;
}

rt_textile::~rt_textile(void)
{
	if(tubes.size()>0) clear_tubes();
	if(curves.size()>0) clear_curves();
}

void rt_textile::set(const double mag, const double per_curve, const double per_x, const double per_y, const int nx, const int ny, 
					 const rt_vector& o, const double r, const double t, const double tf)
{ 
	const double PI=3.141592;

	if(tubes.size()>0) clear_tubes();
	if(curves.size()>0) clear_curves();
	magnitude=mag; period_curve=per_curve; period_x=per_x; period_y=per_y; n_x=nx; n_y=ny; org=o;radius=r; tol=t; tolfun=tf;
	rt_curve* pcurve; 
	rt_tube* ptube;
	rt_vector org_curve=org;
	double sign=-1;
	for(int x=0; x<nx; x++) {
		org_curve=org+rt_vector(x*per_x,0,0);
		pcurve=new rt_sin_curve(magnitude,period_curve,sign*PI/2,org_curve,rt_vector(0,1,0),0,per_y*(ny-1));
		sign*=-1;
		curves.push_back(pcurve);
		ptube=new rt_tube(pcurve,radius,tol,tolfun);
		tubes.push_back(ptube);
	}
	sign=1;
	for(int y=0; y<nx; y++) {
		org_curve=org+rt_vector(0,y*per_y,0);
		pcurve=new rt_sin_curve(magnitude,period_curve,sign*PI/2,org_curve,rt_vector(1,0,0),0,per_x*(nx-1));
		sign*=-1;
		curves.push_back(pcurve);
		ptube=new rt_tube(pcurve,radius,tol,tolfun);
		tubes.push_back(ptube);
	}
}

void rt_textile::clear_tubes(void) {
	vector<rt_tube*>::iterator iter;
	for(iter=tubes.begin(); iter!=tubes.end(); iter++) {
		delete *iter;
	}
	tubes.clear();
}

void rt_textile::clear_curves(void) {
	vector<rt_curve*>::iterator iter;
	for(iter=curves.begin(); iter!=curves.end(); iter++) {
		delete *iter;
	}
	curves.clear();
}

bool rt_textile::find_cross_point(const rt_vector& from, const rt_vector& direction, rt_vector& to, rt_vector& n) const
{
	double epsilon=numeric_limits<float>::epsilon();
	epsilon=(epsilon>tol)?epsilon:tol;
	double dist_min=numeric_limits<double>::max(), dist;
	bool found=false;
	for(int j=0; j<(int)tubes.size(); j++) {
		rt_vector toj, normalj;
		if(tubes[j]->find_cross_point(from,direction,toj,normalj)) {
			dist=rt_vector::norm(from-toj);
			if(dist>tol&&dist<dist_min) {
				dist_min=dist;
				to=toj;
				n=normalj;
				found=true;
			}
		}
	}
	if(found) {
		if(rt_vector::norm(n)==0) {
			return false;
		}
		n=n/rt_vector::norm(n);
		return true;
	}
	return false;
}

bool rt_parabola::find_cross_point(const rt_vector& from, const rt_vector& direction, rt_vector& to, rt_vector& n) const
{
	double epsilon=numeric_limits<float>::epsilon();
	rt_vector x0=from-center;
	double x0sq=rt_vector::dot_product(x0,x0);
	double x0u=rt_vector::dot_product(x0,direction);
	double x0n=rt_vector::dot_product(x0,normal_at_center);
	double un=rt_vector::dot_product(direction,normal_at_center);

	double c2=height*(1-un*un);
	double c1=2*height*(x0u-x0n*un)-un;
	double c0=height*(x0sq-x0n*x0n)-x0n;

	bool found=false;
	double t; 
	if(c2==0) {
		if(c1!=0) {
			t=-c0/c1;
			if(t>epsilon) {
				found=true;
			}
		}
	} else {
		double det=c1*c1-4*c2*c0;
		if(det>0) {
			double t0=(-c1+sqrt(det))/(2*c2);
			double t1=(-c1-sqrt(det))/(2*c2);
			if(t0>epsilon && t1>epsilon) { 
				t=(t0<t1)?t0:t1;
				found=true;
			} else if(t1>epsilon) { // t1 is only solution
				t=t1;
				found=true;
			} else if(t0>epsilon) { // t0 is only solution
				t=t0;
				found=true;
			}
		}
	}
	if(found) {
		to=from+t*direction;
		double d=rt_vector::norm(to-rt_vector::dot_product(to,normal_at_center)*normal_at_center);
		if(d<=d_max) {
			x0=to-center;
			x0n=rt_vector::dot_product(x0,normal_at_center);
			rt_vector n0=normal_at_center;
			rt_vector y=x0-x0n*n0;
			rt_vector grad=n0+2*height*rt_vector::dot_product(n0,y)*n0-2*height*y;
			if(rt_vector::norm(grad)==0) {
				return false;
			}
			n=grad/rt_vector::norm(grad);
			double ug=rt_vector::dot_product(direction,n);
			if(ug>0) n=-n;

			return true;
		}
	}

	return false;
}

bool rt_semi_ellipsoid::find_cross_point(const rt_vector& from, const rt_vector& direction, rt_vector& to, rt_vector& n) const 
{
	double epsilon=numeric_limits<float>::epsilon();
	rt_vector x=from-center;
	double x0=x.x, y0=x.y, z0=x.z;
	double ux=direction.x, uy=direction.y, uz=direction.z;
	double c2=ux*ux/(a*a)+uy*uy/(b*b)+uz*uz/(c*c);
	double c1=2*x0*ux/(a*a)+2*y0*uy/(b*b)+2*z0*uz/(c*c);
	double c0=x0*x0/(a*a)+y0*y0/(b*b)+z0*z0/(c*c)-1;

	bool found=false;
	double t; 
	rt_vector to0, to1;
	if(c2==0) {
		if(c1!=0) {			
			t=-c0/c1;
			to0=from+t*direction;
			if(t>epsilon && to0.z<=center.z) {
				to=to0;
				found=true;
			}
		}
	} else {
		double det=c1*c1-4*c2*c0;
		if(det>0) {
			double t0=(-c1+sqrt(det))/(2*c2);
			double t1=(-c1-sqrt(det))/(2*c2);
			to0=from+t0*direction;
			to1=from+t1*direction;
			// to be a solution t>0 and to.z<=center.z 
			if((t0>epsilon && to0.z<=center.z) && (t1>epsilon && to1.z<=center.z)) { 
				t=(t0<t1)?t0:t1;
				to=(t0<t1)?to0:to1;
				found=true;
			} else if(t1>epsilon && to1.z<=center.z) { // t1 is only solution
				t=t1;
				to=to1;
				found=true;
			} else if(t0>epsilon && to0.z<=center.z) { // t0 is only solution
				t=t0;
				to=to0;
				found=true;
			}
		}
	}
	if(found) {
		rt_vector d=to-center;
		if(abs(d.x)<=a && abs(d.y)<=b && abs(d.z)<=c) {
			rt_vector grad(2*d.x/(a*a),2*d.y/(b*b),2*d.z/(c*c));
			if(rt_vector::norm(grad)==0) {
				return false;
			}
			n=grad/rt_vector::norm(grad);
			double ug=rt_vector::dot_product(direction,n);
			if(ug>0) n=-n;
			return true;
		}
	}

	return false;
}
