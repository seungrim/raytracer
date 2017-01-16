#pragma once

#ifndef _RT_VECTOR_H_
#define _RT_VECTOR_H_

#include <cmath>
#include <iostream>
using namespace std;

class rt_vector
{
public:
	// constructors and destructor
	rt_vector(void) {}
	rt_vector(const rt_vector& v) { set(v.x,v.y,v.z); }
	rt_vector(double xi, double yi, double zi) { set(xi,yi,zi); }
	~rt_vector(void) {}
	
	// getters and setters
	void get(double *xo, double *yo, double *zo) const { *xo=x; *yo=y; *zo=z; }
	void set(const double xi, const double yi, const double zi) { x=xi; y=yi; z=zi; }
	
	// operators
	rt_vector operator+(const rt_vector& v) const { return rt_vector(x+v.x,y+v.y,z+v.z); }
	friend rt_vector& operator+=(rt_vector& v1, const rt_vector& v2) { return v1=v1+v2; }
	rt_vector operator-(const rt_vector& v) const { return rt_vector(x-v.x,y-v.y,z-v.z); }
	rt_vector operator-(void) const { return *this*(-1); }	
	rt_vector operator*(const double s) const { return rt_vector(s*x,s*y,s*z); }
	friend rt_vector operator*(const double s, const rt_vector& v) { return v*s; }
	rt_vector operator/(const double s) const { return *this*(1/s); }
	bool operator==(const rt_vector& v) const { return (x==v.x && y==v.y && z==v.z); }
	bool operator!=(const rt_vector& v) const { return !(x==v.x && y==v.y && z==v.z); }
	rt_vector& operator=(const rt_vector& v) { x=v.x; y=v.y; z=v.z; return *this; }
	friend ostream& operator<<(ostream& os, const rt_vector& v) {
		os << v.x << ' ' << v.y << ' ' << v.z;
		return os;
	}
	
	// static functions
	static double norm(const rt_vector& v) { return sqrt(v.x*v.x+v.y*v.y+v.z*v.z); }
	static double dot_product(const rt_vector& v1, const rt_vector& v2) { 
		return v1.x*v2.x+v1.y*v2.y+v1.z*v2.z; 
	}
	static rt_vector cross_product(const rt_vector& v1, const rt_vector& v2) { 
		return rt_vector(v1.y*v2.z-v1.z*v2.y,v1.z*v2.x-v1.x*v2.z,v1.x*v2.y-v1.y*v2.x); 
	}
	
	// member variables
	double x,y,z;
};

#endif // _RT_VECTOR_H_

