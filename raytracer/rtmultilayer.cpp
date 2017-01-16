#include "rtmultilayer.h"
#include "mlcalc.h"

bool rt_multilayer::calc_efficiencies(const double angle_in, const int pol, double& aqe, double &trans)
{
	double transmitted=0, absorbed=0;

	complex<double> *refractive_index=new (nothrow) complex<double>[nlayer];
	if(refractive_index==NULL) {
		return false;
	}
	for(int j=0; j<nlayer; j++) {
		pmaterial_library->find(materials[j].data(),wavelength,refractive_index[j]);
	}
	calc_qe(nlayer,thickness,refractive_index,wavelength,angle_in,pol,&transmitted,&absorbed);
	aqe=absorbed; trans=transmitted;

	delete [] refractive_index;

	return true;
}

bool rt_multilayer::calc_efficiencies_reverse(const double angle_in, const int pol, double& aqe, double &trans)
{
	double *d=new (nothrow) double[nlayer];
	complex<double> *n=new (nothrow) complex<double>[nlayer];
	double transmitted=0, absorbed=0;
	
	if(d==NULL || n==NULL) {
		return false;
	}

	for(int j=0; j<nlayer; j++) {
		d[j]=thickness[nlayer-j-1];
		pmaterial_library->find(materials[nlayer-j-1].data(),wavelength,n[j]);
	}

	calc_qe(nlayer,d,n,wavelength,angle_in,pol,&transmitted,&absorbed);
	aqe=absorbed; trans=transmitted;

	delete [] d; delete [] n;

	return true;
}
