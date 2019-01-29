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
 and all methods it depends on have been inspired by the SOPHIA code.
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
     @param onProton  true=proton, false=neutron
     @param E_in      energy of primary
     @param z_in      redshift of primary
     */
    double sample_eps(bool onProton, double E_in, double z_in) const;

    /** Returns the photon field density in 1/(cm³eV).
        Multiply by h*nu for physical photon density.
     @param eps   photon energy in eV
     @param z_in  redshift
    */
    double get_photonDensity(double eps, double z) const;

    /** Returns the crossection of p-gamma interaction in µbarn
     @param eps       photon energy in eV
     @param onProton  true=proton, false=neutron
    */
    double crossection(double eps, bool onProton) const;

 private:
    double prob_eps(double eps, bool onProton, double E_in, double z_in) const;
    double Pl(double x, double xth, double xmax, double alpha) const;
    double Ef(double x, double th, double w) const;
    double breitwigner(double sigma_0, double Gamma, double DMM, double epsPrime, bool onProton) const;
    double gaussInt(double lowerLimit, double upperLimit, bool onProton) const;
    double functs(double s, bool onProton) const;
    void init(std::string fieldPath);
        std::vector<double> energy;
        std::vector<double> redshift;
        std::vector<double> photonDensity;
};

/** @}*/
} // namespace crpropa

#endif // CRPROPA_PHOTONBACKGROUND_H
