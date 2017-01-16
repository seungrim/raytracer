///////////////////////////////////////////////
// fminbnd
// minimize 1-dimensional function func
// 
// created by Seung Rim at Mar 14, 2008
//
// references: MATLAB R2007b and Numerical Recipes

#ifndef _FMINBND_H_
#define _FMINBND_H_

double fminbnd(const double ax, const double bx, double &fret, double func(const double, const void*, const void*), 
			   const void* fopts1, const void* fopts2);

#endif // _FMINBND_H_

