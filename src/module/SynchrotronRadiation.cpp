#include "crpropa/module/SynchrotronRadiation.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"

#include <fstream>
#include <limits>
#include <stdexcept>

namespace crpropa {

SynchrotronRadiation::SynchrotronRadiation(ref_ptr<MagneticField> field, bool havePhotons, std::string tag, double limit) {
	Brms = 0.;
	setField(field);
	initSpectrum();
	this->spaceTimeGrid = ScalarGrid4d();
	setDescription("SynchrotronRadiation_isotropicConstant");
	this->havePhotons = havePhotons;
	this->limit = limit;
	this->tag = tag;
	secondaryThresholdLower = 1e7 * eV;
	secondaryThresholdUpper = std::numeric_limits<double>::max();
}

SynchrotronRadiation::SynchrotronRadiation(ref_ptr<MagneticField> field, ScalarGrid4d spaceTimeGrid, bool havePhotons, std::string tag, double limit) {
	Brms = 0.;
	setField(field);
	initSpectrum();
	this->spaceTimeGrid = spaceTimeGrid;
	setDescription("SynchrotronRadiation_spaceTimeDependent");
	this->havePhotons = havePhotons;
	this->limit = limit;
	this->tag = tag;
	secondaryThresholdLower = 1e7 * eV;
	secondaryThresholdUpper = std::numeric_limits<double>::max();
}

SynchrotronRadiation::SynchrotronRadiation(double Brms, bool havePhotons, std::string tag, double limit) {
	this->Brms = Brms;
	initSpectrum();
	this->spaceTimeGrid = ScalarGrid4d();
	setDescription("SynchrotronRadiation_isotropicConstant");
	this->havePhotons = havePhotons;
	this->limit = limit;
	this->tag = tag;
	secondaryThresholdLower = 1e7 * eV;
	secondaryThresholdUpper = std::numeric_limits<double>::max();
}

SynchrotronRadiation::SynchrotronRadiation(double Brms, ScalarGrid4d spaceTimeGrid, bool havePhotons, std::string tag, double limit) {
	this->Brms = Brms;
	initSpectrum();
	this->spaceTimeGrid = spaceTimeGrid;
	setDescription("SynchrotronRadiation_spaceTimeDependent");
	this->havePhotons = havePhotons;
	this->limit = limit;
	this->tag = tag;
	secondaryThresholdLower = 1e7 * eV;
	secondaryThresholdUpper = std::numeric_limits<double>::max();
}

void SynchrotronRadiation::setField(ref_ptr<MagneticField> f) {
	this->field = f;
}

ref_ptr<MagneticField> SynchrotronRadiation::getField() {
	return field;
}

void SynchrotronRadiation::setBrms(double Brms) {
	this->Brms = Brms;
}

double SynchrotronRadiation::getBrms() {
	return Brms;
}

void SynchrotronRadiation::setHavePhotons(bool havePhotons) {
	this->havePhotons = havePhotons;
}

bool SynchrotronRadiation::getHavePhotons() {
	return havePhotons;
}

void SynchrotronRadiation::setLimit(double limit) {
	this->limit = limit;
}

double SynchrotronRadiation::getLimit() {
	return limit;
}

void SynchrotronRadiation::setSecondaryThresholdLower(double threshold) {
	secondaryThresholdLower = threshold;
}

double SynchrotronRadiation::getSecondaryThresholdLower() const {
	return secondaryThresholdLower;
}

void SynchrotronRadiation::setSecondaryThresholdUpper(double threshold) {
	secondaryThresholdUpper = threshold;
}

double SynchrotronRadiation::getSecondaryThresholdUpper() const {
	return secondaryThresholdUpper;
}

void SynchrotronRadiation::initSpectrum() {
	std::string filename = getDataPath("Synchrotron/spectrum.txt");
	std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error(
				"SynchrotronRadiation: could not open file " + filename);

	// clear previously loaded interaction rates
	tabx.clear();
	tabCDF.clear();

	while (infile.good()) {
		if (infile.peek() != '#') {
			double a, b;
			infile >> a >> b;
			if (infile) {
				tabx.push_back(pow(10, a));
				tabCDF.push_back(b);
			}
		}
		infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
	}
	infile.close();
}

void SynchrotronRadiation::process(Candidate *candidate) const {
	double charge = fabs(candidate->current.getCharge());
	if (charge == 0)
		return; // only charged particles

	// geometric scaling
	Vector3d pos = candidate->current.getPosition();
	const double time = candidate->getTrajectoryLength() / c_light;
	double scaling = 1;
    const std::string description = getDescription();
    if (description == "SynchrotronRadiation_isotropicConstant") {
        // do nothing, just check for correct initialization
    } else if (description == "SynchrotronRadiation_spaceTimeDependent") {
        scaling *= spaceTimeGrid.interpolate(pos, time);
    } else {
        throw std::runtime_error("SynchrotronRadiation: invalid description string");
    }
    if (scaling == 0.)
        return;

	// calculate gyroradius, evaluated at the current position
	double z = candidate->getRedshift();
	double B;
	if (field.valid()) {
		Vector3d Bvec = field->getField(pos, z);
		B = Bvec.cross(candidate->current.getDirection()).getR();
	} else {
		B = sqrt(2. / 3) * Brms; // average perpendicular field component
	}
	B *= pow(1 + z, 2); // cosmological scaling
	B *= scaling;
	double Rg = candidate->current.getMomentum().getR() / charge / B;

	// calculate energy loss
	double lf = candidate->current.getLorentzFactor();
	double dEdx = 1. / 6 / M_PI / epsilon0 * pow(lf * lf - 1, 2) * pow(eplus / Rg, 2); // Jackson p. 770 (14.31)
	double step = candidate->getCurrentStep() / (1 + z); // step size in local frame
	double dE = step * dEdx;

	// apply energy loss and limit next step
	double E = candidate->current.getEnergy();
	candidate->current.setEnergy(E - dE);
	candidate->limitNextStep(limit * E / dEdx);

	// optionally add secondary photons
	if (not(havePhotons))
		return;

	// check if photons with energies > 14 * Ecrit are possible
	double Ecrit = 3. / 4 * h_planck / M_PI * c_light * pow(lf, 3) / Rg;
	if (14 * Ecrit < secondaryThresholdLower)
		return;

	// draw photons up to the total energy loss
	Random &random = Random::instance();
	while (dE > 0) {
		// draw random value between 0 and maximum of corresponding cdf
		// choose bin of s where cdf(x) = cdf_rand -> x_rand
		size_t i = random.randBin(tabCDF); // draw random bin (upper bin boundary returned)
		double binWidth = (tabx[i] - tabx[i-1]);
		double x = tabx[i-1] + random.rand() * binWidth; // draw random x uniformly distributed in bin
		double Egamma = x * Ecrit;

		// if the remaining energy is not sufficient check for random accepting
		if (Egamma > dE)
			if (random.rand() > (dE / Egamma))
				break; // not accepted

		// create synchrotron photon and repeat with remaining energy
		dE -= Egamma;
		// create only photons with energies above threshold
		if ((Egamma > secondaryThresholdLower) && (Egamma < secondaryThresholdUpper)) {
			Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
			candidate->addSecondary(22, Egamma, pos, tag);
		}
	}
}

// std::string SynchrotronRadiation::getDescription() const {
// 	std::stringstream s;
// 	s << "Synchrotron radiation";
// 	if (field.valid())
// 		s << " for specified magnetic field";
// 	else
// 		s << " for Brms = " << Brms / nG << " nG";
// 	if (havePhotons)
// 		s << ", synchrotron photons E > " << secondaryThresholdLower / eV << " eV";
// 	else
// 		s << ", no synchrotron photons";
// 	return s.str();
// }

} // namespace crpropa
