#ifndef CRPROPA_HADRONICINTERACTION_H
#define CRPROPA_HADRONICINTERACTION_H

#include "crpropa/Module.h"
#include "crpropa/Vector3.h"

namespace crpropa {
/**
 * \addtogroup EnergyLosses
 * @{
 */

/**
 @class HadronicInteraction
 @brief interactions of nuclei with background nucleons (Hydrogen only).
 */
class HadronicInteraction: public Module {
protected:
	double massDensity;

public:
	HadronicInteraction(double massDensity = 0.);
	void process(Candidate *candidate) const;
	void performInteraction(Candidate *candidate) const;

	double xSectionKelner06(double ePrimary) const;
	double pionSpectrum(double x, double ePrimary) const;
	double etaSpectrum(double x, double ePrimary) const;
	double samplePionEnergy(double ePrimary) const;
	double sampleEtaEnergy(double ePrimary) const;
};

} // namespace crpropa

#endif // CRPROPA_HADRONICINTERACTION_H
