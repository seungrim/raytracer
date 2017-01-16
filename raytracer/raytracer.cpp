/*********************************************************************************************
*
*	raytracer.cpp 
*	05/01/2009 by Seung Rim
*
*	main program to call rtraytracer package
*
*********************************************************************************************/

#include <limits>
#include <fstream>
#include <cstring>
#include <cstdlib>

#include "rtmultilayer.h"
#include "rtnklibrary.h"
#include "rtraytracer.h"

struct init_params {
	string nk_filename;
	string materiali;
	double xi,xstep,xf, yi,ystep,yf, zi;
	double xu,yu,zu;
	double lambdai,lambdastep,lambdaf;
} gparam;

void rt_assert(const bool condition) {}

rt_raytracer& parse_inputfile(rt_raytracer& tracer, const char* filename)
{
	// read input file
	ifstream is_infile;
	string str_temp, str_line;
	is_infile.open(filename);
	rt_assert(is_infile.good());
	gparam.nk_filename="nk.asc";

	while(getline(is_infile, str_line)) {
		// blank line or comment line
		if(str_line.length()==0 || str_line[0]=='%') {
			continue;
		}

		istringstream iss_inline(str_line);
		// instruction
		iss_inline >> str_temp;
		if(str_temp.compare("NKFILE")==0) {
			string str_nkfile;
			iss_inline >> str_nkfile;
			gparam.nk_filename=str_nkfile;
			// check if nkfile exists
			fstream is_nkfile;
			is_nkfile.open(gparam.nk_filename.data());
			rt_assert(is_nkfile.good());
			is_nkfile.close();
		} else if(str_temp.compare("RAYIN")==0) {
			// xi:xstep:xf yi:ystep:yf zi
			double xi,xstep,xf,yi,ystep,yf,zi;
			string str_x;
			iss_inline >> str_x;
			size_t found=str_x.find(':');
			if(found!=string::npos) {
				xi=atof(str_x.substr(0,found).data());
				size_t found2=str_x.find(':',found+1);
				xstep=atof(str_x.substr(found+1,found2).data());
				xf=atof(str_x.substr(found2+1,str_x.size()).data());
			} else {
				xi=xf=atof(str_x.data());
				xstep=10;
			}
			string str_y;
			iss_inline >> str_y;
			found=str_y.find(':');
			if(found!=string::npos) {
				yi=atof(str_y.substr(0,found).data());
				size_t found2=str_y.find(':',found+1);
				ystep=atof(str_y.substr(found+1,found2).data());
				yf=atof(str_y.substr(found2+1,str_y.size()).data());
			} else {
				yi=yf=atof(str_y.data());
				ystep=10;
			}
			iss_inline >> zi;
			gparam.xi=xi; gparam.xstep=xstep; gparam.xf=xf;
			gparam.yi=yi; gparam.ystep=ystep; gparam.yf=yf;
			gparam.zi=zi;
		} else if(str_temp.compare("RAYDIR")==0) {
			// xu yu zu
			double xu,yu,zu;
			iss_inline >> xu >> yu >> zu;
			gparam.xu=xu; gparam.yu=yu; gparam.zu=zu;
		} else if(str_temp.compare("RAYMATERIAL")==0) {
			// material
			string materiali;
			iss_inline >> materiali;
			gparam.materiali=materiali;
		} else if(str_temp.compare("WAVELENGTH")==0) {
			double lambdai,lambdastep,lambdaf;
			// lambdai:lambdastep:lambdaf
			string str_lambda;
			iss_inline >> str_lambda;
			size_t found=str_lambda.find(':');
			if(found!=string::npos) {
				lambdai=atof(str_lambda.substr(0,found).data());
				size_t found2=str_lambda.find(':',found+1);
				lambdastep=atof(str_lambda.substr(found+1,found2).data());
				lambdaf=atof(str_lambda.substr(found2+1,str_lambda.size()).data());
			} else {
				lambdai=lambdaf=atof(str_lambda.data());
				lambdastep=10;
			}
			gparam.lambdai=lambdai; gparam.lambdastep=lambdastep; gparam.lambdaf=lambdaf;
		} else if(str_temp.compare("SURFACE")==0) {
			// shape surfacetype
			// read shape
			rt_shape* p_shape=(rt_shape*)0;
			string str_temp2;
			iss_inline >> str_temp2;
			if(str_temp2.compare("MESH")==0) {
				// meshfilename
				string str_mesh, str_surfacetype;
				iss_inline >> str_mesh;
				// read meshfile
				ifstream is_meshfile;
				is_meshfile.open(str_mesh.data());
				rt_assert(is_meshfile.good());
				rt_triangles *pshape_mesh=new (nothrow) rt_triangles;
				rt_assert(!(pshape_mesh==NULL));
				string str_meshline;
				while(getline(is_meshfile,str_meshline)) {
					if(str_meshline.length()==0 || str_meshline[0]=='%') {
						continue;
					}
					istringstream iss_meshline(str_meshline);
					rt_vector a,b,c;
					double x,y,z;
					iss_meshline >> x >> y >> z;
					a.set(x,y,z);
					iss_meshline >> x >> y >> z;
					b.set(x,y,z);
					iss_meshline >> x >> y >> z;
					c.set(x,y,z);
					pshape_mesh->add(a,b,c);
				}
				is_meshfile.close();
				p_shape=(rt_shape*)pshape_mesh;
			} else {
				// unknown shape
			}			
			// read surfacetype
			rt_surface *psurface=(rt_surface*)0;
			string str_surfacetype;
			iss_inline >> str_surfacetype;
			if(str_surfacetype.compare("FRESNEL")==0) {
				// material1 material2
				string m1, m2;
				iss_inline >> m1;
				iss_inline >> m2;
				psurface=new (nothrow) rt_fresnel_surface(m1,m2);
			} else if(str_surfacetype.compare("MULTILAYER")==0) {
				string str_mltemp;
				vector<string> materials;
				vector<double> thicknesses;
				bool read_front=false, read_end=false;
				while(read_end==false && getline(is_infile,str_mltemp)) {
					if(str_mltemp.length()==0 || str_mltemp[0]=='%') {
						continue;
					}
					istringstream is_multilayer(str_mltemp);
					string str_material, str_thickness;
					is_multilayer >> str_material >> str_thickness;
					materials.push_back(str_material);
					if(str_thickness.length()==0) {
						if(read_front) read_end=true;
						else read_front=true;
						thicknesses.push_back(0.0);
					} else {
						thicknesses.push_back(atof(str_thickness.data()));
					}
				}
				int nlayer=materials.size();
				rt_multilayer* pmultilayer=new (nothrow) rt_multilayer(nlayer,(double*)&thicknesses[0],(string*)&materials[0]);
				psurface=(rt_surface*) new (nothrow) rt_multilayer_surface(pmultilayer,materials.front(),materials.back());
			} else if(str_surfacetype.compare("MIRROR")==0) {
				psurface=(rt_surface*) new (nothrow) rt_mirror_surface;
			} else if(str_surfacetype.compare("ABSORBER")==0) {
				double f_abs=0.0;
				string m1,m2;
				iss_inline >> f_abs >> m1 >> m2;
				psurface=(rt_surface*) new (nothrow) rt_absorb_surface(f_abs,m1,m2);
			} else if(str_surfacetype.compare("ESCAPE")==0) {
				psurface=(rt_surface*) new (nothrow) rt_escape_surface;
			} else {
				// unknown surface
			}
			psurface->add_shape(p_shape);
			tracer.add_surface(psurface);
		} else {
			// unknown instructions
		}
	}
	is_infile.close();

	return tracer;
}

int main(int argc, char* argv[])
{
	string option, param; int nopt=1;
	double lambda=600.0, lambda_start=300, lambda_end=900, lambda_step=100;
	double angle=0.0;
	double min_intensity=5e-4; // 0.1% 
	bool record=false, qeout=false;

	char *in_filename, *out_filename, *rec_filename;
	if(argc<=1) {
		cout << "raytracer -i infile [-q qefile] [-r rayfile] [-intmin minimumintensity]" << endl;
		cout << "-q qefile contains absorption efficiencies at each wavelength." << endl;
		cout << "-r rayfile contains information about rays traced." << endl;
		return 0;
	}
	while(nopt<argc) {
		option.assign(argv[nopt++]);
		if(option.compare("-i")==0) { // input filename
			in_filename=argv[nopt++];
		} else if(option.compare("-q")==0) { // output filename
			qeout=true;
			out_filename=argv[nopt++];
		} else if(option.compare("-r")==0) { // recording options
			record=true;
			rec_filename=argv[nopt++];
		} else if(option.compare("-intmin")==0) {
			min_intensity=atof(argv[nopt++]);
		}
	}

	// minimum seperation to resolve points
	double min_separation=numeric_limits<float>::epsilon();
	// minimum intensity to discard rays
	rt_raytracer tracer(min_intensity,min_separation);
	
	parse_inputfile(tracer, in_filename);
		
	// material libraries to store refractive indicies
	rt_nk_library materials;
	materials.read(gparam.nk_filename.data());
	tracer.set_material_library(&materials);

	double aqes,aqep,aqessurf,aqepsurf,aqesvol,aqepvol;
	int nstepx=(int)((gparam.xf-gparam.xi)/gparam.xstep+1), nstepy=(int)((gparam.yf-gparam.yi)/gparam.ystep+1);

	ofstream os;
	if(qeout) {
		os.open(out_filename);
		if(os.good()) {
			os << "% " << out_filename << " contains wavelength(nm)/absorption(s- and p-pol max:1)/absorption(by volume)/absorption(s-pol max:0.5)/absorption(p-pol max:0.5)." << endl;
			os << "% " << nstepx*nstepy << " rays are traced." << endl;
		} else {
			cout << out_filename << " can't be created." << endl;
		}
	}

	rt_vector x0=rt_vector(gparam.xi,gparam.yi,gparam.zi);
	rt_vector xd=rt_vector(gparam.xstep,0,0), yd=rt_vector(0,gparam.ystep,0);
	rt_vector dir=rt_vector(gparam.xu,gparam.yu,gparam.zu);
	string n0=gparam.materiali;

	for(lambda=gparam.lambdai;lambda<=gparam.lambdaf;lambda+=gparam.lambdastep) {
		tracer.set_wavelength(lambda);
		aqes=0; aqep=0;
		aqessurf=0; aqesvol=0;
		aqepsurf=0; aqepvol=0;

		if(record) tracer.record_begin();
		tracer.trace_bunch(x0,xd,yd,nstepx,nstepy,dir,n0,rt_raytracer::s_pol,0.5,aqes,aqessurf,aqesvol);
		tracer.trace_bunch(x0,xd,yd,nstepx,nstepy,dir,n0,rt_raytracer::p_pol,0.5,aqep,aqepsurf,aqepvol);
		if(record) tracer.record_end();

		if(qeout && os.good()) {
			os << lambda << " " << aqes+aqep << " " << aqesvol+aqepvol << " " << aqes << " " << aqep << endl;
		}
		if(record) {
			ofstream os_record;
			string str_recfile(rec_filename);
			if(lambda_start<=lambda_end) { // always use suffix
				size_t nperiod=str_recfile.find_last_of('.');
				stringstream ss_reclambda;
				ss_reclambda << '_' << (int)lambda;
				str_recfile.insert(nperiod,ss_reclambda.str());
			}
			os_record.open(str_recfile.data());
			if(os_record.good()) {
				os_record << "% " << rec_filename << " contains information about rays traced." << endl;
				os_record << "% " << nstepx*nstepy << " rays are traced." << endl;
				os_record << "% " << lambda <<"nm wavelength is used in the tracing." << endl;
				os_record << "% Each row contains numbers as following: " << endl;
				os_record << "% pol,x_start,y_start,z_start,intensity_start," << 
						"x_end,y_end,z_end,intensity_end,absorption_in_volume," << 
						"alpha_in_volume,absorption_in_surface" << endl;
			} else {
				cout << rec_filename << " can't be created." << endl;
			}
			tracer.record_ray_data(&os_record);
			tracer.discard_ray_data();
			os_record.close();
		}
	}

	os.close();

	return 0;
}
