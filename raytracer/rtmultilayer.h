#pragma once

#ifndef _RT_MULTILAYER_H_
#define _RT_MULTILAYER_H_

#include <complex>
#include "rtnklibrary.h"

using namespace std;

class rt_multilayer
{
public:
	rt_multilayer(void) { 
		thickness=(double*)0; materials=(string*)0;
		wavelength=(double)0.0;
	}
	rt_multilayer(const int nl, double *t, string* m, double lambda=0) {
		thickness=(double*)0; materials=(string*)0; 
		wavelength=(double)0.0;
		set(nl,t,m); set_wavelength(lambda);
	}
	rt_multilayer(const rt_multilayer& m) { 
		thickness=(double*)0; materials=(string*)0; 
		wavelength=(double)0.0; 
		set(m.nlayer,m.thickness,m.materials); 
		set_wavelength(m.wavelength);
	}
	~rt_multilayer(void) { delete thickness; delete materials; }

	void set(const int nl, double *t, string* m) {
		nlayer=nl;
		if(thickness) delete [] thickness; thickness=new double[nlayer]; 
		if(materials) delete [] materials; materials=new string[nlayer];
		for(int j=0; j<nlayer; j++) {
			thickness[j]=t[j]; materials[j]=m[j];
		}
	}
	void set_wavelength(const double lambda) { wavelength=lambda; }
	void set_material_library(rt_nk_library* plib) { pmaterial_library=plib; }
	double get_wavelength(void) const { return wavelength; }
	int get_nlayer(void) const { return nlayer; }
	string& get_material_front(void) const { return materials[0]; }

	static const int s_pol=0, p_pol=1;
	bool calc_efficiencies(const double angle_in, const int pol, double& aqe, double &trans);
	bool calc_efficiencies_reverse(const double angle_in, const int pol, double& aqe, double &trans);

private:
	int nlayer;
	double *thickness;
	string *materials;
	double wavelength;

	rt_nk_library *pmaterial_library;
};

#endif // _RT_MULTILAER_H_

