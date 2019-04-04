#ifndef CRPROPA_HADRONICINTERACTION_H
#define CRPROPA_HADRONICINTERACTION_H

#include "crpropa/Module.h"

#include <vector>

namespace crpropa {


class HadronicInteraction: public Module {

private:
	bool haveElectrons;
	bool havePhotons;
	bool haveNeutrinos;
	
public:
	HadronicInteraction(bool electrons = true, bool photons = true, bool neutrinos = true);
	void setHaveElectrons(bool b);
	void setHavePhotons(bool b);
	void setHaveNeutrinos(bool b);
    double getDensity(const Vector3d &position) const;
    //crpropa::Vector3d MagField(const Vector3d &position) const;
	void process(Candidate *candidate) const;
	//double distribution_pi(double energy, double x) const;
	double distribution_e(double energy, double x) const;
	double distribution_my1(double energy, double x) const; 
	double distribution_gamma(double energy, double x) const; 
	//double number_pi(double energy) const;
	double number_e(double energy) const;
	double number_my1(double energy) const;
	double number_gamma(double energy) const;
	double distribution_Carceller(double energy, double x, double jcap, double a0, double b0) const;
	double distribution_Carceller_g(double energy, double x, double jcap, double a0, double b0) const;
	double CrossSection_Carceller(double energy) const;
	double CrossSection_Kelner(double energy) const;
	double CrossSection_Galprop(double energy) const;
	//double CrossSection_Kafexhiu(double energy) const;
    crpropa::Vector3d Position(double height, double radius) const;
	//~ double counter(Candidate *candidate, double density) const;
	//~ double counterPion(Candidate *candidate, double density) const;
};

} // namespace crpropa

#endif // CRPROPA_HADRONICINTERACTION_H
