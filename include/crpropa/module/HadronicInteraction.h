#ifndef CRPROPA_HADRONICINTERACTION_H
#define CRPROPA_HADRONICINTERACTION_H

#include "crpropa/Module.h"
#include "crpropa/Vector3.h"
#include "crpropa/Random.h"

namespace crpropa {
/**
 * \addtogroup EnergyLosses
 * @{
 */

/**
 @class HadronicInteraction
 @brief interactions of nuclei with background protons. Work based on Kelner et al. 2006
 */
class HadronicInteraction: public Module {
protected:
	double massDensity;
	double limit;

public:
	HadronicInteraction(
		double massDensity = 0.,
		double limit = 0.1);
	void setMassDensity(double massDensity);
	void setLimit(double limit);
	void process(Candidate *candidate) const;
	void performInteraction(Candidate *candidate) const;

	double xSectionKelner06(double ePrimary) const;
	double spectrumPion(double x, double ePrimary) const;
	double spectrumPhoton(double x, double ePrimary) const;

	template<typename SpectrumFunction>
	double sampleParticleEnergy(SpectrumFunction&& spectrumFunction, double ePrimary) const {
		const double xMin = 1. / 1000.;
		const double xMax = 1.;
		const double stepSize = 1. / 1000.;

		double Fmax = 0.;
		double x = xMin;
		while (x < xMax) {
			const double F = spectrumFunction(x);
			if (F > Fmax)
				Fmax = F;
			x += stepSize;
		}

		Random &random = Random::instance();
		double F = 0.;
		do {
			x = std::pow(10, -3 * random.rand());
			F = spectrumFunction(x);
		} while (F < random.rand() * Fmax);
		return x * ePrimary;
	}
};

} // namespace crpropa

#endif // CRPROPA_HADRONICINTERACTION_H
