#ifndef CRPROPA_PHOTOPIONPRODUCTIONMOD_H
#define CRPROPA_PHOTOPIONPRODUCTIONMOD_H

#include "crpropa/Module.h"
#include "crpropa/PhotonBackground.h"

#include <vector>

namespace crpropa {
/**
 * \addtogroup EnergyLosses
 * @{
 */

/**
 @class PhotoPionProductionMod
 @brief Photo-pion interactions of nuclei with background photons.
 */
class PhotoPionProductionMod: public Module {
protected:
	PhotonField photonField;
	std::vector<double> tabLorentz; ///< Lorentz factor of nucleus
	std::vector<double> tabRedshifts;  ///< redshifts (optional for haveRedshiftDependence)
	std::vector<double> tabProtonRate; ///< interaction rate in [1/m] for protons
	std::vector<double> tabNeutronRate; ///< interaction rate in [1/m] for neutrons
	double limit; ///< fraction of mean free path to limit the next step
	bool havePhotons;
	bool haveNeutrinos;
	bool haveElectrons;
	bool haveAntiNucleons;
	bool haveRedshiftDependence;

public:
	PhotoPionProductionMod(
		PhotonField photonField = CMB,
		bool photons = false,
		bool neutrinos = false,
		bool electrons = false,
		bool antiNucleons = false,
		double limit = 0.1,
		bool haveRedshiftDependence = false);
	void setPhotonField(PhotonField photonField);
	void setHavePhotons(bool b);
	void setHaveNeutrinos(bool b);
	void setHaveElectrons(bool b);
	void setHaveAntiNucleons(bool b);
	void setHaveRedshiftDependence(bool b);
	void setLimit(double limit);
	void initRate(std::string filename);
	double nucleonMFP(double gamma, double z, bool onProton) const;
	double nucleiModification(int A, int X) const;
	void process(Candidate *candidate) const;
	void performInteraction(Candidate *candidate, bool onProton) const;

	/**
	 Calculates the loss length E dx/dE in [m].
	 This is not used in the simulation.
	 @param	id		PDG particle id
	 @param gamma	Lorentz factor of particle
	 @param z		redshift
	 */
	double lossLength(int id, double gamma, double z = 0);

	/**
	 Direct SOPHIA interface.
	 This is not used in the simulation.
	 Returns a vector of length 2x the amount of produced particles;
	 the first half contains their IDs, the second their energy.
	 @param nature       primary proton or neutron
	 @param Ein          Energy of interacting nucleon
	 @param eps          Energy of scattering photon
	*/
	std::vector<double> sophiaEvent(bool onProton, double Ein, double eps) const;
};
/** @}*/

} // namespace crpropa

#endif // CRPROPA_PHOTOPIONPRODUCTIONMOD_H