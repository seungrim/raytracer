//////////////////////////////////////////////////////////////////////////
// RayTracing Shape Classes
//
// Mar 25, 2008 Seung Rim
//
// Shapes consist of a surface.
// A shape class has a find_cross_point function which finds out a cross point between ray and a shape.
// classes implemented here
//	rt_shape - pure abstract class
//		rt_2d_shape:rt_shape - abstract base class for 2d planar shapes
//			rt_triangle, rt_rectangle, rt_circle :rt_2d_shape
//		rt_triangles:rt_shape - set of triangles useful for shapes with triangular meshes
//			rt_pyramid - pyramidal structure consist of four triangles
//		rt_tube - represents a tube with 1-dim curve and radius from the curve
//					rt_curve (no base class) - represents a curve
//						rt_sin_curve:rt_curve - a sinusoidal curve

#pragma once

#ifndef _RT_SHAPE_H_
#define _RT_SHAPE_H_

#include <vector>
#include "rtvector.h"

using namespace std;

class rt_shape {
public:
	virtual bool find_cross_point(const rt_vector& from, const rt_vector& direction, rt_vector& to, rt_vector& normal) const = 0;
};

class rt_2d_shape: public rt_shape {
public:
	// constructors and destructor
	rt_2d_shape(void) {};
	rt_2d_shape(const rt_2d_shape& s) { set(s.normal); }
	rt_2d_shape(const rt_vector& v) { set(v); }
	// getters and setters
	void set(const rt_vector& v) { normal=v; }
	rt_vector& get_normal(void) {return normal; }
	// rt_shape function
	bool find_cross_point(const rt_vector& from, const rt_vector& direction, rt_vector& to, rt_vector& n) const;

	virtual bool inside(const rt_vector& v) const=0;
	virtual rt_vector get_one_vector(void) const=0;
protected:
	rt_vector normal;
};

class rt_triangle: public rt_2d_shape {
public:
	// constructors and destructor
	rt_triangle(void) {}
	rt_triangle(const rt_triangle& tri) { set(tri.vertex[0],tri.vertex[1],tri.vertex[2]); }
	rt_triangle(const rt_vector& a, const rt_vector& b, const rt_vector& c) { set(a,b,c); }
	~rt_triangle(void) {}
	
	// getters and setters
	void set(const rt_vector& a, const rt_vector& b, const rt_vector& c) {
		vertex[0]=a; vertex[1]=b; vertex[2]=c; set_normal(); 
	}
	
	// rd_2d_shape function
	bool inside(const rt_vector& v) const;
	rt_vector get_one_vector(void) const { return vertex[0]; }

private:	
	void set_normal(void) { rt_vector n=rt_vector::cross_product(vertex[1]-vertex[0],vertex[2]-vertex[0]); rt_2d_shape::set(n/rt_vector::norm(n)); }
	rt_vector vertex[3];
};

class rt_rectangle: public rt_2d_shape {
public:
	// constructors and destructor
	rt_rectangle(void) {}
	rt_rectangle(const rt_rectangle &rc) { *this=rc; }
	rt_rectangle(const rt_vector& a, const rt_vector& edge1, const rt_vector& edge2) { set(a,a+edge1,a+edge2,a+edge1+edge2); }
	~rt_rectangle(void) {}

	// getters and setters
	void set(const rt_vector& a, const rt_vector& b, const rt_vector& c, const rt_vector& d) { 
		triangles[0]=rt_triangle(a,b,c); triangles[1]=rt_triangle(b,c,d); 
		rt_2d_shape::set(triangles[0].get_normal());
	}

	// rd_2d_shape function
	bool inside(const rt_vector& v) const { return (triangles[0].inside(v) || triangles[1].inside(v)); }
	rt_vector get_one_vector(void) const { return triangles[0].get_one_vector(); };
private:
	rt_triangle triangles[2];
};

class rt_circle: public rt_2d_shape {
public:
	// constructors and destructor
	rt_circle(void) {}
	rt_circle(const rt_circle& c) { set(c.center,c.normal,c.radius); }
	rt_circle(const rt_vector& c, const rt_vector& n, const double r) { set(c,n,r); }
	~rt_circle(void) {}
	
	// getters and setters
	void set(const rt_vector& c, const rt_vector& n, const double r) { center=c; normal=n; radius=r; }
	
	// rd_2d_shape function
	bool inside(const rt_vector& v) const { 
		return (rt_vector::dot_product(v-center,normal)==0 && rt_vector::norm(v-center)<=radius); 
	}
	
protected:
	// rd_2d_shape function
	rt_vector get_one_vector(void) const { return center; };
private:
	rt_vector center;
	double radius;
};

class rt_triangles: public rt_shape {
public:
	rt_triangles(void) {}
	rt_triangles(const int ntri,rt_vector* a,rt_vector* b,rt_vector* c) { set(ntri,a,b,c); }
	~rt_triangles(void) {
		for(int j=0; j<(int)triangles.size(); j++) {
			delete triangles[j]; triangles[j]=(rt_triangle*)0;
		}
		triangles.clear();
	}
	void set(const int ntri,rt_vector* a,rt_vector* b,rt_vector* c) {
		for(int j=0; j<ntri; j++) {
			rt_triangle* ptri=new (nothrow) rt_triangle(a[j],b[j],c[j]);
			if(ptri) {
				triangles.push_back(ptri);
			}
		}
	}
	void add(const rt_vector& a,const rt_vector& b, const rt_vector& c) {
		rt_triangle* ptri=new (nothrow) rt_triangle(a,b,c);
		if(ptri) {
			triangles.push_back(ptri);
		}
		triangles.push_back(ptri);
	}
	bool find_cross_point(const rt_vector& from, const rt_vector& direction, rt_vector& to, rt_vector& n) const;
private:
	vector<rt_triangle*> triangles;
};

class rt_pyramid: public rt_triangles {
public:
	rt_pyramid(void) {}
	rt_pyramid(const rt_vector& a, const rt_vector& b, const rt_vector& c, const rt_vector& d, const rt_vector &tip)		
	{
		rt_vector x[]={a,b,c,d};
		rt_vector y[]={b,c,d,a};
		rt_vector z[]={tip,tip,tip,tip};
		rt_triangles::set(4,x,y,z);
	}
};

class rt_curve {
public:
	// constructors and destructor
	rt_curve(void) {}
	rt_curve(const double t_min, const double t_max) { set(t_min,t_max); }
	~rt_curve(void) {}
	
	// getters and setters
	void set(const double t_min, const double t_max) { tmin=t_min; tmax=t_max; }
	void gett(double &t_min, double &t_max) { t_min=tmin; t_max=tmax; }
	
	// public functions
	double find_nearest_point(const rt_vector& from, rt_vector& on_curve, double &t_curve) const;
	virtual rt_vector get_point_at(const double t) const = 0;
	
private:
	double tmin, tmax;
};

class rt_sin_curve: public rt_curve {
public:
	// constructors and destructor
	rt_sin_curve(void) {}
	rt_sin_curve(const double mag, const double per, const double pha, const rt_vector& x0, const rt_vector& u, 
				const double t_min, const double t_max):rt_curve(t_min,t_max) {
		set(mag,per,pha,x0,u); 
	}
	~rt_sin_curve(void) {}
	// getters and setters
	void set(const double mag, const double per, const double pha, const rt_vector& x0, const rt_vector& u) {
		magnitude=mag; period=per; phase=pha; x_initial=x0; x_direction=u;
	}
	// rt_curve functions
	rt_vector get_point_at(const double t) const {
		const double PI=3.141592;
		rt_vector x=x_initial+t*x_direction;
		x.z=magnitude*sin(2*PI/period*t+phase)+x.z;
		return x;
	}
private:
	double magnitude, period, phase;
	rt_vector x_initial, x_direction;
};

class rt_tube: public rt_shape {
public:
	// constructors and destructor
	rt_tube(const double t=0, const double tf=0) {tol=t; tolfun=tf; curve=(rt_curve*)0; }
	rt_tube(const rt_tube& t) { curve=(rt_curve*)0; set(t.curve,t.radius,t.tol,t.tolfun); }
	rt_tube(rt_curve* c, const double r,const double t=0,const double tf=0) { curve=(rt_curve*)0; set(c,r,t,tf); }
	~rt_tube(void) { if(curve) delete curve; }
	
	// getters and setters
	void set(rt_curve* c, const double r, const double t, const double tf) { if(curve) delete curve; curve=c; radius=r; tol=t; tolfun=tf; }
	double get_radius(void) const { return radius; }
	rt_curve& get_curve(void) const { return *curve; }
	
	// rt_shape function
	bool find_cross_point(const rt_vector& from, const rt_vector& direction, rt_vector& to, rt_vector& n) const;
	
private:
	rt_curve* curve;
	double radius;
	double tol, tolfun;
};

class rt_textile: public rt_shape {
public:
	rt_textile(const double t=0, const double tf=0) {tol=t;tolfun=tf;}
	rt_textile(const double mag, const double per_curve, const double per_x, const double per_y, const int nx, const int ny, 
		const rt_vector& o,const double r,const double t=0, const double tf=0)
	{ set(mag,per_curve,per_x,per_y,nx,ny,o,r,t,tf); }
	rt_textile(const rt_textile& t) { set(t.magnitude,t.period_curve,t.period_x,t.period_y,t.n_x,t.n_y,t.org,t.radius); }
	~rt_textile(void);
	void set(const double mag, const double per_curve, const double per_x, const double per_y, const int nx, const int ny, 
		const rt_vector& o, const double r,const double t=0,const double tf=0);
	bool find_cross_point(const rt_vector& from, const rt_vector& direction, rt_vector& to, rt_vector& n) const;

private:
	void clear_tubes(void);
	void clear_curves(void);

	double magnitude, period_curve, period_x, period_y, radius;
	int n_x, n_y;
	rt_vector org;
	vector<rt_tube*> tubes;
	vector<rt_curve*> curves;
	double tol, tolfun;
};

// parabolic curves
// (x-c,n)=h*|(x-c)-(x-c,n)n|^2
class rt_parabola: public rt_shape {
public:
	rt_parabola(void) {}
	rt_parabola(const rt_parabola& p) {set(p.center,p.normal_at_center,p.height,p.d_max);}
	rt_parabola(const rt_vector& c, const rt_vector& n,const double h,const double tmax) {set(c,n,h,tmax);}
	~rt_parabola(void) {}

	void set(const rt_vector& c, const rt_vector& n,const double h,const double tmax) 
	{center=c; normal_at_center=n;height=h;d_max=tmax;}
	bool find_cross_point(const rt_vector& from, const rt_vector& direction, rt_vector& to, rt_vector& n) const;

private:
	rt_vector center, normal_at_center;
	double height, d_max;
};

// ellipsoid
// x^2/a^2+y^2/b^2+z^2/c^2=1
class rt_semi_ellipsoid: public rt_shape {
public:
	rt_semi_ellipsoid(void) {}
	rt_semi_ellipsoid(const rt_vector& o, const double a, const double b, const double c) {set(o,a,b,c);}
	rt_semi_ellipsoid(const rt_semi_ellipsoid& e) {set(e.center,e.a,e.b,e.c);}
	~rt_semi_ellipsoid(void) {}

	void set(const rt_vector& o, const double ai, const double bi, const double ci) { center=o; a=ai; b=bi; c=ci; }
	bool find_cross_point(const rt_vector& from, const rt_vector& direction, rt_vector& to, rt_vector& n) const;

private:
	rt_vector center;
	double a,b,c;
};

#endif //_RT_SHAPE_H_
