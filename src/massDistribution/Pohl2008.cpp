#include "crpropa/massDistribution/Pohl2008.h"

#include <fstream>
#include <sstream>

namespace crpropa {


void Pohl08::loadGridHI() {
	std::ifstream fin("/rest/CRPropa3/share/crpropa/Pohl_HI.txt");
	if(!fin) {
		std::stringstream ss;
		ss << "load Pohl Grid: Pohl_HI.txt not found";
		throw std::runtime_error(ss.str());
	}
	//skp header lines
	while (fin.peek() == '#')
		fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	
	for(int ix = 0; ix<400;ix++) {
		for(int iy = 0; iy<400;iy++) {
			for(int iz = 0;iz<30;iz++) {
			
				float &b = HIdensity.get(ix,iy,iz);
				fin >>b;
				
				if(fin.eof()) 
					throw std::runtime_error("Reduced Pohl_HI: File too short");
			}
		}
	}
	fin.close();
}

void Pohl08::loadGridH2() {
	std::ifstream fin("/rest/CRPropa3/share/crpropa/Pohl_H2.txt");
	if(!fin) {
		std::stringstream ss;
		ss << "load Pohl Grid: Pohl_H2.txt not found";
		throw std::runtime_error(ss.str());
	}
	//skp header lines
	while (fin.peek() == '#')
		fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	
	for(int ix = 0; ix<300;ix++) {
		for(int iy = 0; iy<300;iy++) {
			for(int iz = 0;iz<10;iz++) {
			
				float &b = H2density.get(ix,iy,iz);
				fin >>b;
				
				if(fin.eof()) 
					throw std::runtime_error("Reduced Pohl_H2: File too short");
			}
		}
	}
	fin.close();
}

double Pohl08::getH2Density(const Vector3d &position) const {
	
	double n = 0; //density in ccm
	Vector3d pos = position;
	
	if(fabs(pos.x)>15*kpc)	//boundary of Grid not repeat periodicly
	{
		return 0;
	}
	if(fabs(pos.y)>15*kpc)
	{
		return 0;
	}
	if(fabs(pos.z)>500*pc)
	{
		return 0;
	}
	
	n= H2density.interpolate(pos);
	
	// check if density is NAN
	// return 0 instead
	bool NaN = std::isnan(n);
	if(NaN == true){
		KISS_LOG_WARNING
			<< "\nDensity with 'nan' occured:\n"
			<< "postion = " << position << "\n"
			<< "density-model: Pohl 2008 \n"
			<< "density-type: H2 (molecular)\n"
			<< "density is set to 0. \n";
			return 0;
	}
	
	return n/cm;
}


double Pohl08::getHIDensity(const Vector3d &position) const {
	
	double n = 0;	// density in ccm
	Vector3d pos = position;
	
	if(fabs(pos.x)>20*kpc)	//boundary of Grid not repeat periodicly
	{
		return 0;
	}
	if(fabs(pos.y)>20*kpc)
	{
		return 0;
	}
	if(fabs(pos.z)>1500*pc)
	{
		return 0;
	}
	
	n= HIdensity.interpolate(pos);
	
	// check if density is NAN
	// return 0 instead
	bool NaN = std::isnan(n);
	if(NaN == true){
		KISS_LOG_WARNING
			<< "\nDensity with 'nan' occured:\n"
			<< "postion = " << position << "\n"
			<< "density-model: Pohl 2008 \n"
			<< "density-type: HI (atomic)\n"
			<< "density is set to 0. \n";
			return 0;
	}
	
	return n/cm;
	
}
double Pohl08::getDensity(const Vector3d &position) const {
	double n=0;
	if(isforHI)
		n+=Pohl08::getHIDensity(position);
	if(isforH2)
		n+=Pohl08::getH2Density(position);
		
	//check if any density is activ and give warning if not
	bool anyDensityActive = isforHI||isforH2;

	if(anyDensityActive == false){
		KISS_LOG_WARNING
			<< "\n tryed to get density although all density-types are deaktivated \n"
			<< "density-module: Ferriere\n"
			<< "returned 0 density\n"
			<< "please use constant Density with 0 \n";
	}	
	
	return n;
}


bool Pohl08::getisforHI() {
	return isforHI;
}

bool Pohl08::getisforHII() {
	return isforHII;
}

bool Pohl08::getisforH2() {
	return isforH2;
}

void Pohl08::setisforHI(bool HI) {
	isforHI=HI;
}

void Pohl08::setisforH2(bool H2) {
	isforH2=H2;
}

} //namespace
