#pragma once

#ifndef _RT_NK_LIBRARY_H_
#define _RT_NK_LIBRARY_H_

#include <vector>
#include <complex>
#include <map>
#include <fstream>

using namespace std;

class rt_nk_library
{
public:
	typedef map< double, complex<double> > material;
	typedef map< string, material* > material_map;
	typedef pair< double, complex<double> > material_pair;
	typedef pair< string, material* > material_map_pair;

	rt_nk_library(void) {}
	rt_nk_library(const rt_nk_library& l) { *this=l; }
	~rt_nk_library(void) {
		material_map::iterator iter;
		for(iter=materials.begin(); iter!=materials.end(); iter++) {
			if(iter->second) delete iter->second;
		}
		materials.clear();
	}

	bool find(const char* name, const double lambda, complex<double>& n) {
		material_map::iterator iter_map=materials.find(string(name));
		if(iter_map!=materials.end()) {
			material::const_iterator iter=iter_map->second->find(lambda);
			if(iter!=iter_map->second->end()) {
				n=iter->second;
				return true;
			}
		}
		return false;
	}

	bool read(const char* filename) {
		ifstream is(filename);
		if(!is.good()) {
			return false;
		}
		string str_line, name;
		bool b_name_stored=false;
		while(b_name_stored || getline(is,str_line)) { // read material's name
			if(str_line.length()==0 || str_line[0]=='%') {
				continue;
			}
			name=str_line;
			// erase control characters
			size_t nerase;
			nerase=name.find_first_of("\t\r\b\v\n");
			while(nerase!=string::npos) {
				name.erase(nerase,1);
				nerase=name.find_first_of("\t\r\b\v\n");
			}
			// check if the materials is in library
			material_map::iterator iter_map=materials.find(name);
			if(iter_map==materials.end()) { // material is not in this library
				material *pmaterial=new (nothrow) material;
				if(pmaterial==NULL) {
					is.close();
					return false;
				}
				// read n and k data
				while(getline(is,str_line)) {
					// found material's name
					size_t found;
					found=str_line.find_first_not_of(' ');
					if(str_line.length()>0 && isdigit(str_line[found])==0) {
						b_name_stored=true;
						break;
					}
					if(str_line.length()==0 || str_line[0]=='%') {
						continue;
					}
					istringstream iss_line(str_line);
					double lambda,n,k;
					iss_line >> lambda >> n >> k;
					pmaterial->insert(material_pair(lambda,complex<double>(n,k)));
					b_name_stored=false;
				}
				materials.insert(material_map_pair(string(name),pmaterial));
			}
		}
		is.close();
		return true;
	}

private:
	material_map materials;
};

#endif // _RT_NK_LIBRARY_H_
