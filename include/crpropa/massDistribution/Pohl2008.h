#ifndef CRPROPA_POHL2008_H
#define CRPROPA_POHL2008_H

#include "crpropa/Units.h"
#include "crpropa/Vector3.h"
#include "crpropa/massDistribution/Density.h"
#include "crpropa/Grid.h"
#include "crpropa/GridTools.h"

#include <string>

#include "kiss/logger.h"


namespace crpropa {

/**
	@class Pohl08 
	@brief Grid for H2 density based on arxiv:0712.4264 Fits avalibal at http://www.app.physik.uni-potsdam.de/gas.html 
	Grid for HI density based on private communication with Martin Pohl.
*/

class Pohl08: public Density {
	ScalarGrid H2density = ScalarGrid(Vector3d(-15*kpc,-15*kpc,-500*pc),300,300,10,100*pc);
	ScalarGrid HIdensity = ScalarGrid(Vector3d(-20*kpc,-20*kpc,-1500*pc),400,400,30,100*pc);	

	bool isforHI = true;
	bool isforHII = false;
	bool isforH2 = true; 


public:

	
	void loadGridHI();
	void loadGridH2();
	double getDensity(const Vector3d &position) const;
	double getH2Density(const Vector3d &position) const;
	double getHIDensity(const Vector3d &position) const;
	
	bool getisforHI();
	bool getisforHII();
	bool getisforH2();
	
	void setisforHI(bool HI);
	void setisforH2(bool H2);

	
};


} //namespace crpropa

#endif //CRPROPA_POHL2008_H
