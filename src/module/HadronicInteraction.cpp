#include "crpropa/module/HadronicInteraction.h"
#include "crpropa/Common.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"
#include "crpropa/Random.h"

#include <fstream>
#include <limits>
#include <cmath>
#include <stdexcept>

namespace crpropa {


HadronicInteraction::HadronicInteraction(double massDensity, double limit) {
	this->massDensity = massDensity;
	this->limit = limit;
	setDescription("HadronicInteraction");
}

void HadronicInteraction::setMassDensity(double dens) {
	this->massDensity = massDensity;
}

void HadronicInteraction::setLimit(double lim) {
	this->limit = lim;
}

double HadronicInteraction::spectrumPion(double x, double ePrimary) const {
	const double E0pi = 139.571 * MeV;
	if (x < E0pi / ePrimary)
		return 0.;
	const double L = std::log(ePrimary / TeV);
	const double a = 3.67 + 0.83 * L + 0.075 * L * L;
	const double Bpi = a + 0.25;
	const double alpha = 0.98 / std::sqrt(a);
	const double r = 2.6 / std::sqrt(a);
	const double xa = std::pow(x, alpha);
	const double term1 = 4. * alpha * Bpi * std::pow(x, alpha - 1.);
	const double term2 = std::pow((1. - xa) / (1. + r * xa * (1. - xa)), 4);
	const double term3 = 1. / (1. - xa) + r * (1. - 2. * xa) / (1. + r * xa * (1. - xa));
	const double term4 = std::sqrt(1. - E0pi  / (x * ePrimary));
	const double Fpi = term1 * term2 * term3 * term4;
	return Fpi;
}

double HadronicInteraction::spectrumPhoton(double x, double ePrimary) const {
	const double L = std::log(ePrimary / TeV);
	const double B = 1.30 + 0.14 * L + 0.011 * L * L;
	const double beta = 1. / (1.79 + 0.11 * L + 0.008 * L * L);
	const double k = 0.801 + 0.049 * L + 0.014 * L * L;
	const double xb = std::pow(x, beta);
	const double term1 = B * std::log(x) / x;
	const double term2 = (1. - xb) / (1 + k * xb * (1 - xb));
	const double term3 = 1 / std::log(x) - 4 * beta * xb / (1 - xb);
	const double term4 = 4 * k * beta * xb * (1 - 2 * xb) / (1 + k * xb * (1 - xb));
	const double Fgamma = term1 * std::pow(term2, 4) * (term3 - term4);
	return Fgamma;
}

double HadronicInteraction::xSectionKelner06(double ePrimary) const {
	const double L = std::log(ePrimary / TeV);
	const double A = 1 - std::pow(1.22 * 1e-3 * TeV / ePrimary, 4);
	return (34.3 + 1.88 * L + 0.25 * L * L) * A * A * 1e-31;
}

void HadronicInteraction::process(Candidate *candidate) const {
	if (this->massDensity <= 0.)
		return;

	const int id = candidate->current.getId();
	if (!isNucleus(id))
		return;

	const double ePrimary = candidate->current.getEnergy();
	if (ePrimary < 1. * GeV)
		return;

	const double xSection = xSectionKelner06(ePrimary);
	const double meanFreePath = 1. / (xSection * this->massDensity);
	const double stepLength = candidate->getCurrentStep();
	const double interactionProbability = stepLength / meanFreePath;

	Random &random = Random::instance();
	if (random.rand() > interactionProbability)
		candidate->limitNextStep(this->limit * meanFreePath);
		return;

	performInteraction(candidate);
}

void HadronicInteraction::performInteraction(Candidate *candidate) const {
	std::vector<int> outPartID;
	std::vector<double> outPartE;
	double eAvailable = candidate->current.getEnergy();

	Random &random = Random::instance();
	const double startPiont = random.rand();
	bool doPhoton = (startPiont >= 0. && startPiont < 0.5);
	bool doPiCharged = (startPiont >= 0.5 && startPiont <= 1.);
	bool donePhoton = false;
	bool donePiCharged = false;

	do {
		if (doPhoton && !donePhoton) {
			// const double xMin = 1 * GeV / eAvailable;
			const double xMin = std::min(1 * GeV / eAvailable, 1e-3);
			const double xMax = 1.;
			int nPhoton = std::round(gaussInt([this, eAvailable](double x) { return this->spectrumPhoton(x, eAvailable); }, xMin, xMax));
			for (int i = 0; i < nPhoton; ++i) {
				const double ePhoton = sampleParticleEnergy([this, eAvailable](double x) { return this->spectrumPhoton(x, eAvailable); }, eAvailable);
				if (eAvailable >= ePhoton) {
					outPartID.push_back(22);
					outPartE.push_back(ePhoton);
					eAvailable -= ePhoton;
				} else {
					break;
				}
			}
			doPiCharged = true;
			donePhoton = true;
		}

		if (doPiCharged && !donePiCharged) {
			const double xMin = 0.;
			const double xMax = 1.;
			int nPiCharged = 2 * std::round(gaussInt([this, eAvailable](double x) { return this->spectrumPion(x, eAvailable); }, xMin, xMax));
			for (int i = 0; i < nPiCharged / 2; ++i) {
				const double ePiPlus = sampleParticleEnergy([this, eAvailable](double x) { return this->spectrumPion(x, eAvailable); }, eAvailable);
				const double ePiMinus = sampleParticleEnergy([this, eAvailable](double x) { return this->spectrumPion(x, eAvailable); }, eAvailable);
				if (eAvailable >= ePiPlus + ePiMinus) {
					outPartID.push_back(211);
					outPartE.push_back(ePiPlus);
					eAvailable -= ePiPlus;

					outPartID.push_back(-211);
					outPartE.push_back(ePiMinus);
					eAvailable -= ePiMinus;
				} else {
					break;
				}
			}
			doPhoton = true;
			donePiCharged = true;
		}
	} while (!donePhoton || !donePiCharged);

	candidate->current.setEnergy(eAvailable);
	// std::cout << "ePrimary = " << eAvailable / GeV << std::endl;
	const Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
	for (int i = 0; i < outPartID.size(); ++i) {
		candidate->addSecondary(outPartID[i], outPartE[i], pos);
		std::cout << outPartID[i] << "\t" << outPartE[i] / GeV << std::endl;
	}
}

} // namespace CRPropa
