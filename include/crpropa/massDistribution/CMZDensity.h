#ifndef CRPROPA_CMZDENSITY_H
#define CRPROPA_CMZDENSITY_H

#include "crpropa/massDistribution/Density.h"
#include "crpropa/Vector3.h"
#include "crpropa/Units.h"

#include <math.h>

namespace crpropa {
    /*
     ..
     */
    
    class CMZDensity: public Density {
        
    public:
        CMZDensity();
        Vector3d CMZTrafo(const Vector3d &position) const;
        double getDensity(const Vector3d &position) const;
        double getH2Density(const Vector3d &position) const;
 
        
    };
    
} //namespace

#endif //CRPROPA_CMZDENSITY_H

