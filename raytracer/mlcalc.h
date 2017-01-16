#ifndef _ML_CALC_H_
#define _ML_CALC_H_

//////////////////////////////////////////////////////////////////////////
// mlcalc
// multilayer structure calculator for thin film photovoltaic cells
//
// created by Seung Rim at Mar 12, 2008
// See details in cpp codes
#include <complex>
using namespace std;

// Optical interference calculation 
int calc_absorption_qe(int nlayers, double *d, complex<double> *n, complex<double> n_front, double lambda, double phi0,
			   int pol, double *T, double *ro, double *delta, complex<double> *tp, double *aqe);

// calculate external quantum efficiency for planar junction pv with information from optical calculation
int calc_external_qe(int nlayers, double *vd, complex<double> *vn, complex<double> n_front, double lambda, double phi0, 
					 double *vT, double *vro, double *vdelta, double *vLd, double *vsi, double *vsf,
					 double *qe);

// calculate qe of planar junction cells
// this function just calls calc_absorption_qe and calc_external_qe
// this handles absorption only
int calc_qe(int nlayers, double *d, complex<double> *n, double lambda, double phi0, int pol,double *trans, double *aqe);

#endif // _ML_CALC_H_

