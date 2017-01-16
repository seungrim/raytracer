//////////////////////////////////
// multilayer calculation
//
// compiler information:
// 	This codes use C++ Standard Template Library (STL).
// 	Visual C++ in Windows and g++ in linux are tested to compile this codes.
// 	For gcc, you need to include cstdlib to include STL.
//
// created by Seung Rim at Mar 12, 2008
// Apr 02, 2008 correct to obtain coefficients at the last layer 
// Apr 04, 2008 bhjflags added to calculate BHJ using exponential carrier collection
// May 01, 2009 calc_qe takes care of absorption only
// Jun 30, 2009 q=sqrt(n^2-n0.real()^2*sin^2) -->.real() added

#include "mlcalc.h"

#include <cmath>
#include <vector>
#include <complex>
#include <float.h>

using namespace std;

#define PI 3.141592
template<class T> inline const T SQR(const T a) {return a*a;}

int calc_qe(int nlayers, double *d, complex<double> *n, double lambda, double phi0, int pol, double *trans, double *aqe)
{
	double *T=new (nothrow) double[nlayers];
	double *ro=new (nothrow) double[nlayers];
	double *delta=new (nothrow) double[nlayers];
	complex<double> *tp=new (nothrow) complex<double>[nlayers];
	double *aqe_layer=new (nothrow) double[nlayers];
	double eta_a=0;

	if(T==NULL || ro==NULL || delta==NULL || tp==NULL || aqe_layer==NULL) {
		return -1;
	}

	calc_absorption_qe(nlayers, d, n, n[0], lambda, phi0, pol, T, ro, delta, tp, aqe_layer);

	for(int j=0; j<nlayers-1; j++) {
		eta_a=eta_a+aqe_layer[j];
	}
	*aqe=eta_a;
	*trans=T[nlayers-1];

	delete [] T; delete [] ro; delete [] delta; delete [] tp; delete [] aqe_layer;;

	return 0;
}

// Main optical interference effects calculation here 
/* optical profile calculation for pv cell structure
   Output params: tp, rpp and q
   E(x) can be obtained using the formula
    E(x)=tp(exp(i*xi*x)+rpp*exp(-i*xi*(2*d-x)))|E0+|;
	where xi=2*PI*q/lambda
*/
int calc_absorption_qe(int nlayers, double *d, complex<double> *n, complex<double> n_front, double lambda, double phi0, int pol,
		   double *T, double *ro, double *delta, complex<double> *tp, double *aqe)
{
	int layers=nlayers;
	int m=layers-2; /* index of layer at end-1 */

	vector< complex<double> > xi(layers);
	vector< complex<double> > r(layers);
	vector< complex<double> > t(layers);
	vector< complex<double> > rpp(layers);
	vector< complex<double> > q(layers);

	int j;
	for(j=0;j<layers;j++) {
		q[j]=sqrt(SQR(n[j])-SQR(n_front.real())*SQR(sin(phi0))); // n_front to n_front.real() at 06/30/2009
		xi[j]=2*PI/lambda*q[j];
	}
	if(pol==1) { // p-polarization
		for(j=0; j<layers-1; j++) {
			if(SQR(n[j+1])*q[j]+SQR(n[j])*q[j+1]==(double)0) { // brewster angle
				if(SQR(n[j+1])*q[j]==SQR(n[j])*q[j+1]) { // total transmission
					r[j]=0;
					t[j]=1;
				} else { // total reflection
					r[j]=1;
					t[j]=0;
				}
			} else {
				r[j]=-(SQR(n[j+1])*q[j]-SQR(n[j])*q[j+1])/(SQR(n[j+1])*q[j]+SQR(n[j])*q[j+1]); // debugged
				t[j]=((double)2*n[j]*n[j+1]*q[j])/(SQR(n[j+1])*q[j]+SQR(n[j])*q[j+1]);
			}
		}
	} else { // s-polarization
		for(j=0; j<layers-1; j++) {
			if(q[j]+q[j+1]==(double)0) {
				if(q[j]==q[j+1]) {
					r[j]=0;
					t[j]=1;
				} else {
					r[j]=1;
					t[j]=0;
				}
			} else {
				r[j]=(q[j]-q[j+1])/(q[j]+q[j+1]);
				t[j]=2.0*q[j]/(q[j]+q[j+1]);
			}
		}
	}

	vector< complex<double> > rtemp(layers); /* Sm(2,1)/Sm(1,1) */
	vector< complex<double> > ttemp(layers); /* 1/Sp(1,1) */
	vector< complex<double> > rmtemp(layers); /* -Sp(1,2)/Sp(1,1) */

	vector< vector< complex<double> > > S(2,vector< complex<double> >(2,0)), Stemp(2,vector< complex<double> >(2,0));	
	vector< vector< complex<double> > > IL(2,vector< complex<double> >(2,0)), LI(2,vector< complex<double> >(2,0));
	vector< complex<double> > exp_ixid(layers), exp_mixid(layers);
	for(j=1; j<layers; j++) {		
		exp_ixid[j]=exp(complex<double>(0,1)*xi[j]*d[j]);
		exp_mixid[j]=1.0/exp_ixid[j];
	}
	/* Set rtemp at end to zero */
	rtemp[m+1]=0;
	/* Initialize Sm, S=Im(m+1)*/	
	S[0][0]=S[1][1]=1.0/t[m]; 
	S[0][1]=S[1][0]=r[m]/t[m];
	rtemp[m]=S[1][0]/S[0][0];
	for(j=m-1; j>=0; j--) {
		IL[0][0]=exp_mixid[j+1]/t[j];
		IL[0][1]=r[j]*exp_ixid[j+1]/t[j];
		IL[1][0]=r[j]*exp_mixid[j+1]/t[j];
		IL[1][1]=exp_ixid[j+1]/t[j];
		Stemp=S;
		S[0][0]=IL[0][0]*Stemp[0][0]+IL[0][1]*Stemp[1][0];
		S[0][1]=IL[0][0]*Stemp[0][1]+IL[0][1]*Stemp[1][1];
		S[1][0]=IL[1][0]*Stemp[0][0]+IL[1][1]*Stemp[1][0];
		S[1][1]=IL[1][0]*Stemp[0][1]+IL[1][1]*Stemp[1][1];
		rtemp[j]=S[1][0]/S[0][0];
	}
	/* Set ttemp and rmtemp at front to zero */
	ttemp[0]=1;	rmtemp[0]=0;
	/* Initialize Sp, S=I01 */
	S[0][0]=S[1][1]=1.0/t[0];
	S[0][1]=S[1][0]=r[0]/t[0];
	ttemp[1]=1.0/S[0][0];
	rmtemp[1]=-S[0][1]/S[0][0];
	for(j=1; j<=m; j++) {
		LI[0][0]=exp_mixid[j]/t[j];
		LI[0][1]=r[j]*exp_mixid[j]/t[j];
		LI[1][0]=r[j]*exp_ixid[j]/t[j];
		LI[1][1]=exp_ixid[j]/t[j];
		Stemp=S;
		S[0][0]=Stemp[0][0]*LI[0][0]+Stemp[0][1]*LI[1][0];
		S[0][1]=Stemp[0][0]*LI[0][1]+Stemp[0][1]*LI[1][1];
		S[1][0]=Stemp[1][0]*LI[0][0]+Stemp[1][1]*LI[1][0];
		S[1][1]=Stemp[1][0]*LI[0][1]+Stemp[1][1]*LI[1][1];
		ttemp[j+1]=1.0/S[0][0];
		rmtemp[j+1]=-S[0][1]/S[0][0];
	}
	/* calculate tp and rpp */
	tp[0]=1; rpp[0]=rtemp[0];
	for(j=1; j<layers; j++) {
		tp[j]=ttemp[j]/(1.0-rmtemp[j]*rtemp[j]*SQR(exp_ixid[j]));
		rpp[j]=rtemp[j];
	}	
	// glass/air interface
	double f;
	complex<double> q_front=n_front*cos(phi0);
	if(n_front!=n[0]) {
		complex<double> q0=q[0], n0=n[0], r0=rpp[0];
		complex<double> qf=q[m+1], tf=tp[m+1];
		complex<double> r_inc, t_inc;
		double R_inc, R_int;
		if(pol==1) { // p-polarization
			r_inc=-(SQR(n0)*q_front-SQR(n_front)*q0)/(SQR(n0)*q_front+SQR(n_front)*q0); // debugged
			t_inc=2.0*n_front*n0*q_front/(n0*n0*q_front+n_front*n_front*q0);
		} else if(pol==0) { // s-polarization
			r_inc=(q_front-q0)/(q_front+q0);
			t_inc=2.0*q_front/(q_front+q0);
		}
		R_inc=SQR(abs(r_inc)); // T_inc=q0.real()/q_front.real()*SQR(abs(t_inc));
		R_int=SQR(abs(r0));	// T_int=qf.real()/q_front.real()*SQR(abs(tf));
		f=(1-R_inc)/(1-R_inc*R_int);
		// R_tot=1-(R_inc+R_int-R_inc*R_int)/(1-R_inc*R_int)
		// f=(1-R_tot)/(1-R_int)=(1-R_inc)/(1-R_inc*R_int)
		// T_tot=T_inc*T_tot/(1-R_inc*R_int);
	} else {
		f=1;
	}
	// tp and rpp to ro and delta
	for(j=0; j<layers; j++) {
		ro[j]=abs(rpp[j]);
		delta[j]=arg(rpp[j]);
	}
	for(j=0; j<layers; j++) {
		tp[j]=sqrt(f)*tp[j];
		T[j]=q[j].real()/q[0].real()*SQR(abs(tp[j]));
	}
	// aqe calculation
	for(j=0; j<layers; j++) {
		double eta=q[j].real();
		double alpha=4*PI*q[j].imag()/lambda;
		double ead=exp(-alpha*d[j]);
		if(eta==0) { // refractive index is purely imaginary
			aqe[j]=0;	// no absorption
		} else {
			aqe[j]=T[j]*((1.0+SQR(ro[j])*ead)*(1-ead)-2*ro[j]*ead*lambda*alpha/(4*PI*eta)*
					(sin(delta[j])-sin(4*PI*eta*d[j]/lambda+delta[j])));
		}
	}

	return 0;
}

// calculate external quantum efficiency for bilayer 1dim thin film pv
int calc_external_qe(int nlayers, double *vd, complex<double> *vn, complex<double> n_front, double lambda, double phi0, 
					 double *vT, double *vro, double *vdelta, double *vLd, double *vsi, double *vsf,
					 double *qe)
{
	complex<double> q;

	int layers=nlayers;
	int j;
	for(j=0; j<layers; j++) {
		q=sqrt(SQR(vn[j])-SQR(n_front)*SQR(sin(phi0)));
		
		double eta=q.real();
		double alpha=4*PI*q.imag()/lambda;
		double d=vd[j], T=vT[j], ro=vro[j], delta=vdelta[j], ld=vLd[j];
		double si=vsi[j], sf=vsf[j];
		double qei=0.0, qef=0.0;

		if(ld==0 || (si==0 && sf==0) || d==0) { // no eqe can be obtained
			qe[j]=0;
		} else {
			double A,B,C1,C2;
			double beta=1/ld;
			double ead=exp(-alpha*d), ebd=exp(beta*d);
			double k=4*PI*eta/lambda;

			C1=SQR(ro*ead);
			C2=(SQR(beta)-SQR(alpha))*2*ro*ead/(SQR(beta)+SQR(k));
			
			if(si==0 && sf==0) {
				qei=qef=0;
			} else if(si<0 && sf<0) {
				A=((ebd-ead)+C1*(ebd-1/ead)+C2*(ebd*cos(k*d+delta)-cos(delta)))
					/(1/ebd-ebd);
				B=-((1/ebd-ead)+C1*(1/ebd-1/ead)+C2*(1/ebd*cos(k*d+delta)-cos(delta)))
					/(1/ebd-ebd);
			} else if(si<0) {
				A=(-((alpha-sf)*ead-C1*(alpha+sf)/ead+C2*(-sf*cos(delta)-k*sin(delta)))
					-ebd*(beta+sf)*(1+C1+C2*cos(k*d+delta)))
					/((beta+sf)*ebd+(beta-sf)/ebd);
				B=(((alpha-sf)*ead-C1*(alpha+sf)/ead+C2*(-sf*cos(delta)-k*sin(delta)))
					-1/ebd*(beta-sf)*(1+C1+C2*cos(k*d+delta)))
					/((beta+sf)*ebd+(beta-sf)/ebd);
			} else if(sf<0) {
				A=((beta-si)*(ead+C1/ead+C2*cos(delta))+ebd
					*((alpha+si)-C1*(alpha-si)+C2*(si*cos(k*d+delta)-k*sin(k*d+delta))))
					/(-(beta+si)*ebd-(beta-si)/ebd);
				B=((beta+si)*(ead+C1/ead+C2*cos(delta))-1/ebd
					*((alpha+si)-C1*(alpha-si)+C2*(si*cos(k*d+delta)-k*sin(k*d+delta))))
					/(-(beta+si)*ebd-(beta-si)/ebd);
			} else {
				A=((beta-si)*((alpha-sf)*ead-C1*(alpha+sf)/ead+C2*(-sf*cos(delta)-k*sin(delta)))
					-ebd*(beta+sf)*((alpha+si)-C1*(alpha-si)+C2*(si*cos(k*d+delta)-k*sin(k*d+delta))))
					/((beta+sf)*(beta+si)*ebd-(beta-sf)*(beta-si)/ebd);
				B=((beta+si)*((alpha-sf)*ead-C1*(alpha+sf)/ead+C2*(-sf*cos(delta)-k*sin(delta)))
					-1/ebd*(beta-sf)*((alpha+si)-C1*(alpha-si)+C2*(si*cos(k*d+delta)-k*sin(k*d+delta))))
					/((beta+sf)*(beta+si)*ebd-(beta-sf)*(beta-si)/ebd);
			}

			if(si==0||si==-2) {
				qei=0;
			} else {
				qei=alpha*T/(SQR(beta)-SQR(alpha))*(-beta*A+beta*B-alpha+alpha*C1+k*C2*sin(k*d+delta));
			}
			if(sf==0||sf==-2) {
				qef=0;
			} else {
				qef=alpha*T/(SQR(beta)-SQR(alpha))*(beta*A/ebd-beta*B*ebd+alpha*ead
					-alpha*C1/ead-k*C2*sin(delta));
			}
			qe[j]=qei+qef;
		} // if no eqe condition
	} // for layers

	return 0;
}

