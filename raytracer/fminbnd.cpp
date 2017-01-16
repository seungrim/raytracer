////////////////////////////////////////////
// fminbnd
// minimize func given ax<x<bx
// reference: matlab R2007b
// seung rim at mar 16, 2008

#include <cmath>
#include <limits>
#include <iostream>

using namespace std;

inline int sign(double x) { if(x>0) return 1; else if(x<0) return -1; else return 0; }
inline double max(double a, double b) { return (a>b)?a:b; }

double fminbnd(const double ax, const double bx, double &fret, double func(const double, const void*, const void*), 
			   const void* fopts1, const void* fopts2)
{
	const double tol=numeric_limits<float>::epsilon();
	double tol1, tol2;
	const int maxiter=1024;
	double a,b,c,d,e,v,w,x,xf,xm;
	double p,q,r;
	double fx,fu,fv,fw;
	double si;
	int funccount=0,iter=0,gs;	
	double seps;
	
	// Compute the start point
	seps=sqrt(numeric_limits<double>::epsilon());
	c = 0.5*(3.0 - sqrt(5.0));
	a = ax; b = bx;
	v = a + c*(b-a);
	w = v; xf = v;
	d = 0.0; e = 0.0;
	x= xf; fx = func(x,fopts1,fopts2);
	funccount = funccount + 1;

	fv = fx; fw = fx;
	xm = 0.5*(a+b);
	tol1=seps*abs(xf)+tol/3.0;
	tol2=2.0*tol1;

	// Main loop
	while ( (abs(xf-xm) > (tol2 - 0.5*(b-a))) && (iter<maxiter) ) {
	    gs = 1;
	    // Is a parabolic fit possible
	    if (abs(e) > tol1) {
	        // Yes, so fit parabola
	        gs = 0;
	        r = (xf-w)*(fx-fv);
	        q = (xf-v)*(fx-fw);
	        p = (xf-v)*q-(xf-w)*r;
	        q = 2.0*(q-r);
	        if (q > 0.0)  p = -p;
	        q = abs(q);
	        r = e;  e = d;

	        // Is the parabola acceptable
	        if ( (abs(p)<abs(0.5*q*r)) && (p>q*(a-xf)) && (p<q*(b-xf)) ) {
	            // Yes, parabolic interpolation step
	            d = p/q;
	            x = xf+d;
	            
	            // f must not be evaluated too close to ax or bx
	            if(((x-a) < tol2) || ((b-x) < tol2)) {
	                si = sign(xm-xf) + (((xm-xf) == 0)?1:0);
	                d = tol1*si;
	            }
	        } else {
	            // Not acceptable, must do a golden section step
	            gs=1;
	        }
	    }
	    if (gs>0) {
	        // A golden-section step is required
	        if (xf >= xm) e = a-xf; else e = b-xf;
	        d = c*e;
	    }

	    // The function must not be evaluated too close to xf
	    si = sign(d) + ((d == 0)?1:0);
	    x = xf + si * max( abs(d), tol1 );
	    fu = func(x,fopts1,fopts2);
	    funccount = funccount + 1;
	    iter = iter + 1;

	    // Update a, b, v, w, x, xm, tol1, tol2
	    if (fu <= fx) {
	        if (x >= xf) a = xf; else b = xf; 
	        v = w; fv = fw;
	        w = xf; fw = fx;
	        xf = x; fx = fu;
	    } else { // fu > fx
	        if (x < xf) a = x; else b = x;
	        if ( (fu <= fw) || (w == xf) ) {
	            v = w; fv = fw;
	            w = x; fw = fu;
	        } else if ( (fu <= fv) || (v == xf) || (v == w) ) {
	            v = x; fv = fu;
	        }
	    }
	    xm = 0.5*(a+b);
		tol1=seps*abs(xf)+tol/3.0; tol2=2.0*tol1;
	}
	
	fret=fx;
	
	return xf;
}
/*
// test of fminbnd function
double ftest(double x, const void* f1, const void* f2)
{
	return x*x;
}

int main(void) 
{
	double fval,xmin;
	xmin=fminbnd(-1,3,fval,ftest,(void*)0,(void*)0);
	cout << "xmin=" << xmin << "fval=" << fval << endl;
	return 0;
}	
*/

