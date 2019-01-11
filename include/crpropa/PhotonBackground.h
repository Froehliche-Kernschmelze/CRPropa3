#ifndef CRPROPA_PHOTONBACKGROUND_H
#define CRPROPA_PHOTONBACKGROUND_H

#include <crpropa/Vector3.h>
#include <crpropa/Grid.h>

#include <iostream>  // cout, endl
#include <cmath>  // sqrt, pow
#include <string>
#include <fstream>  // write to file
#include <algorithm>  // max_element
#include <limits>  // for ::max
#include <vector>

namespace crpropa {

/**
 * \addtogroup EnergyLosses
 * @{
 */
// Photon fields
// The default IRB model is that of Kneiske et al. 2004
// The slots PF1 to PF8 may be used for custom photon fields
enum PhotonField {
	CMB,
	IRB,  // same as IRB_Kneiske04
	IRB_Kneiske04,
	IRB_Stecker05,
	IRB_Franceschini08,
	IRB_Finke10,
	IRB_Dominguez11,
	IRB_Gilmore12,
	IRB_Stecker16_upper,
	IRB_Stecker16_lower,
	URB_Protheroe96,
	PF1, PF2, PF3, PF4,  // customizable
	PF5, PF6, PF7, PF8,  // field slots
};

// Returns overall comoving scaling factor
double photonFieldScaling(PhotonField photonField, double z);

// Returns a string representation of the field
std::string photonFieldName(PhotonField photonField);


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~
// custom photon field methods
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~

/** 
 @class Photon_Field
 @brief Handler class for photon fields. Provides the sample_eps method.

 sample_eps draws a photon from a given photon background. This method 
 and all methods it depends on have been taken from the SOPHIA code.
 */
class Photon_Field {
 public:
    /** Constructor for photon field data
     @param fieldPath  path/to/photonField.txt
     */
    explicit Photon_Field(std::string fieldPath);

    /* Empty constructor to ease initialization in some modules
     */
    Photon_Field();

    /** Draws a photon from the photon background
     @param onProton  Primary particle type. true=proton, false=neutron
     @param E_in      Energy of the primary in GeV
     @param z_in      Redshift of primary
     */
    double sample_eps(bool onProton, double E_in, double z_in) const;

 private:
    void init(std::string fieldPath);
        std::vector<double> energy;
        std::vector< std::vector<double> > dn_deps;
        std::vector<double> redshift;
    double get_photonDensity(double eps, int z_pos) const;
    double gaussInt(std::string type, double lowerLimit, double upperLimit, bool onProton, double E_in, int z_pos) const;
    double functs(double s, bool onProton) const;
    double prob_eps(double eps, bool onProton, double E_in, int z_pos) const;
    double crossection(double eps, bool onProton) const;
        double Pl(double eps, double eps_threshold, double eps_max, double weight) const;
        double Ef(double eps, double eps_threshold, double threshold) const;
        double breitwigner(double xsection, double width, double m_resonance, double eps, bool onProton) const;
};

/** @}*/
} // namespace crpropa

#endif // CRPROPA_PHOTONBACKGROUND_H
