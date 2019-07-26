#include "crpropa/module/PhotoPionProductionMod.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"
#include "crpropa/Random.h"

#include <kiss/convert.h>
#include <kiss/logger.h>
#include "sophia.h"

#include <limits>
#include <cmath>
#include <sstream>
#include <fstream>
#include <stdexcept>

namespace crpropa {

PhotoPionProductionMod::PhotoPionProductionMod(PhotonField field, bool photons, bool neutrinos, bool electrons, bool antiNucleons, double l, bool redshift) {
	havePhotons = photons;
	haveNeutrinos = neutrinos;
	haveElectrons = electrons;
	haveAntiNucleons = antiNucleons;
	haveRedshiftDependence = redshift;
	limit = l;
	setPhotonField(field);
}

void PhotoPionProductionMod::setPhotonField(PhotonField field) {
	photonField = field;
	if (haveRedshiftDependence) {
		std::cout << "PhotoPionProductionMod: tabulated redshift dependence not needed for CMB, switching off" << std::endl;
		haveRedshiftDependence = false;
	}
	std::string fname = photonFieldName(field);
	setDescription("PhotoPionProductionMod: " + fname);
	if (haveRedshiftDependence)
		initRate(getDataPath("PhotoPionProduction/rate_" + fname.replace(0, 3, "IRBz") + ".txt"));
	else
		initRate(getDataPath("PhotoPionProduction/rate_" + fname + ".txt"));
}

void PhotoPionProductionMod::setHavePhotons(bool b) {
	havePhotons = b;
}

void PhotoPionProductionMod::setHaveElectrons(bool b) {
	haveElectrons = b;
}

void PhotoPionProductionMod::setHaveNeutrinos(bool b) {
	haveNeutrinos = b;
}

void PhotoPionProductionMod::setHaveAntiNucleons(bool b) {
	haveAntiNucleons = b;
}

void PhotoPionProductionMod::setHaveRedshiftDependence(bool b) {
	haveRedshiftDependence = b;
	setPhotonField(photonField);
}

void PhotoPionProductionMod::setLimit(double l) {
	limit = l;
}

void PhotoPionProductionMod::initRate(std::string filename) {
	// clear previously loaded tables
	tabLorentz.clear();
	tabRedshifts.clear();
	tabProtonRate.clear();
	tabNeutronRate.clear();

	std::ifstream infile(filename.c_str());
	if (!infile.good())
		throw std::runtime_error("PhotoPionProductionMod: could not open file " + filename);

	if (haveRedshiftDependence) {
		double zOld = -1, aOld = -1;
		while (infile.good()) {
			if (infile.peek() == '#') {
				infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
				continue;
			}
			double z, a, b, c;
			infile >> z >> a >> b >> c;
			if (!infile)
				break;
			if (z > zOld) {
				tabRedshifts.push_back(z);
				zOld = z;
			}
			if (a > aOld) {
				tabLorentz.push_back(pow(10, a));
				aOld = a;
			}
			tabProtonRate.push_back(b / Mpc);
			tabNeutronRate.push_back(c / Mpc);
		}
	} else {
		while (infile.good()) {
			if (infile.peek() == '#') {
				infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
				continue;
			}
			double a, b, c;
			infile >> a >> b >> c;
			if (!infile)
				break;
			tabLorentz.push_back(pow(10, a));
			tabProtonRate.push_back(b / Mpc);
			tabNeutronRate.push_back(c / Mpc);
		}
	}

	infile.close();
}

double PhotoPionProductionMod::nucleonMFP(double gamma, double z, bool onProton) const {
	const std::vector<double> &tabRate = (onProton)? tabProtonRate : tabNeutronRate;

	// scale nucleus energy instead of background photon energy
	gamma *= (1 + z);
	if (gamma < tabLorentz.front() or (gamma > tabLorentz.back()))
		return std::numeric_limits<double>::max();

	double rate;
	if (haveRedshiftDependence)
		rate = interpolate2d(z, gamma, tabRedshifts, tabLorentz, tabRate);
	else
		rate = interpolate(gamma, tabLorentz, tabRate) * photonFieldScaling(photonField, z);

	// cosmological scaling
	rate *= pow(1 + z, 2);

	return 1. / rate;
}

double PhotoPionProductionMod::nucleiModification(int A, int X) const {
	if (A == 1)
		return 1.;
	if (A <= 8)
		return 0.85 * pow(X, 2. / 3.);
	return 0.85 * X;
}

void PhotoPionProductionMod::process(Candidate *candidate) const {
	double step = candidate->getCurrentStep();
	double z = candidate->getRedshift();
	// the loop is processed at least once for limiting the next step
	do {
		// check if nucleus
		int id = candidate->current.getId();
		if (!isNucleus(id))
			return;

		// find interaction with minimum random distance
		Random &random = Random::instance();
		double randDistance = std::numeric_limits<double>::max();
		double meanFreePath;
		double totalRate = 0;
		bool onProton = true; // interacting particle: proton or neutron

		int A = massNumber(id);
		int Z = chargeNumber(id);
		int N = A - Z;
		double gamma = candidate->current.getLorentzFactor();

		// check for interaction on protons
		if (Z > 0) {
			meanFreePath = nucleonMFP(gamma, z, true) / nucleiModification(A, Z);
			randDistance = -log(random.rand()) * meanFreePath;
			totalRate += 1. / meanFreePath;
		}
		// check for interaction on neutrons
		if (N > 0) {
			meanFreePath = nucleonMFP(gamma, z, false) / nucleiModification(A, N);
			totalRate += 1. / meanFreePath;
			double d = -log(random.rand()) * meanFreePath;
			if (d < randDistance) {
				randDistance = d;
				onProton = false;
			}
		}

		// check if interaction does not happen
		if (step < randDistance) {
			if (totalRate > 0.)
				candidate->limitNextStep(limit / totalRate);
			return;
		}

		// interact and repeat with remaining step
		performInteraction(candidate, onProton);
		step -= randDistance;
	} while (step > 0);
}

void PhotoPionProductionMod::performInteraction(Candidate *candidate, bool onProton) const {
	int id = candidate->current.getId();
	int A = massNumber(id);
	int Z = chargeNumber(id);
	double E = candidate->current.getEnergy();
	double EpA = E / A;
	double z = candidate->getRedshift();

	// SOPHIA simulates interactions only for protons / neutrons
	// for anti-protons / neutrons assume charge symmetry and change all
	// interaction products from particle <--> anti-particle
	int sign = (id > 0) ? 1 : -1;

	// arguments for SOPHIA
	int nature = 1 - int(onProton); // interacting particle: 0 for proton, 1 for neutron
	double Ein = EpA / GeV; // energy of in-going nucleon in GeV
	double momentaList[5][2000]; // momentum list, what are the five components?
	int particleList[2000]; // particle id list
	int nParticles; // number of outgoing particles
	double maxRedshift = 100; // IR photon density is zero above this redshift
	int dummy1; // not needed
	double dummy2[2]; // not needed
	int background = (photonField == CMB) ? 1 : 2; // photon background: 1 for CMB, 2 for Kneiske IRB

	// check if below SOPHIA's energy threshold
	double E_threshold = (photonField == CMB) ? 3.72e18 * eV : 5.83e15 * eV;
	if (EpA * (1 + z) < E_threshold)
		return;

#pragma omp critical
	{
		sophiaevent_(nature, Ein, momentaList, particleList, nParticles, z, background, maxRedshift, dummy1, dummy2, dummy2);
	}

	Random &random = Random::instance();
	Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
	std::vector<int> pnType;  // filled with either 13 (proton) or 14 (neutron)
	std::vector<double> pnEnergy;  // corresponding energies of proton or neutron
	for (int i = 0; i < nParticles; i++) { // loop over out-going particles
		double Eout = momentaList[3][i] * GeV; // only the energy is used; could be changed for more detail
		int pType = particleList[i];
		switch (pType) {
		case 13: // proton
		case 14: // neutron
			// proton and neutron data is taken to determine primary particle in a later step
			pnType.push_back(pType);
			pnEnergy.push_back(Eout);
			break;
		case -13: // anti-proton
		case -14: // anti-neutron
			if (haveAntiNucleons)
				try
				{
					candidate->addSecondary(-sign * nucleusId(1, 14 + pType), Eout, pos);
				}
				catch (std::runtime_error &e)
				{
					KISS_LOG_ERROR<< "Something went wrong in the PhotoPionProductionMod (anti-nucleon production)\n" << "Something went wrong in the PhotoPionProductionMod\n"<< "Please report this error on https://github.com/CRPropa/CRPropa3/issues including your simulation setup and the following random seed:\n" << Random::instance().getSeed_base64();
					throw;
				}
			break;
		case 1: // photon
			if (havePhotons)
				candidate->addSecondary(22, Eout, pos);
			break;
		case 2: // positron
			if (haveElectrons)
				candidate->addSecondary(sign * -11, Eout, pos);
			break;
		case 3: // electron
			if (haveElectrons)
				candidate->addSecondary(sign * 11, Eout, pos);
			break;
		case 15: // nu_e
			if (haveNeutrinos)
				candidate->addSecondary(sign * 12, Eout, pos);
			break;
		case 16: // antinu_e
			if (haveNeutrinos)
				candidate->addSecondary(sign * -12, Eout, pos);
			break;
		case 17: // nu_muon
			if (haveNeutrinos)
				candidate->addSecondary(sign * 14, Eout, pos);
			break;
		case 18: // antinu_muon
			if (haveNeutrinos)
				candidate->addSecondary(sign * -14, Eout, pos);
			break;
		default:
			throw std::runtime_error("PhotoPionProductionMod: unexpected particle " + kiss::str(pType));
		}
	}
	double maxEnergy = *std::max_element(pnEnergy.begin(), pnEnergy.end());  // criterion for being declared primary
	for (int i = 0; i < pnEnergy.size(); ++i) {
		if (pnEnergy[i] == maxEnergy) {  // nucleon is primary particle
			if (A == 1) {
				// single interacting nucleon
				candidate->current.setEnergy(pnEnergy[i]);
				try
				{
					candidate->current.setId(sign * nucleusId(1, 14 - pnType[i]));
				}
				catch (std::runtime_error &e)
				{
					KISS_LOG_ERROR<< "Something went wrong in the PhotoPionProductionMod (primary particle, A==1)\n" << "Please report this error on https://github.com/CRPropa/CRPropa3/issues including your simulation setup and the following random seed:\n" << Random::instance().getSeed_base64();
					throw;
				}
			} else {
				// interacting nucleon is part of nucleus: it is emitted from the nucleus
				candidate->current.setEnergy(E - EpA);
				try
				{
					candidate->current.setId(sign * nucleusId(A - 1, Z - int(onProton)));
					candidate->addSecondary(sign * nucleusId(1, 14 - pnType[i]), pnEnergy[i], pos);
				}
				catch (std::runtime_error &e)
				{
					KISS_LOG_ERROR<< "Something went wrong in the PhotoPionProductionMod (primary particle, A!=1)\n" << "Please report this error on https://github.com/CRPropa/CRPropa3/issues including your simulation setup and the following random seed:\n" << Random::instance().getSeed_base64();
					throw;
				}
			}
		} else {  // nucleon is secondary proton or neutron
			candidate->addSecondary(sign * nucleusId(1, 14 - pnType[i]), pnEnergy[i], pos);
		}
	}
}

double PhotoPionProductionMod::lossLength(int id, double gamma, double z) {
	int A = massNumber(id);
	int Z = chargeNumber(id);
	int N = A - Z;

	double lossRate = 0;
	if (Z > 0)
		lossRate += 1 / nucleonMFP(gamma, z, true) * nucleiModification(A, Z);
	if (N > 0)
		lossRate += 1 / nucleonMFP(gamma, z, false) * nucleiModification(A, N);

	// approximate the relative energy loss
	// - nucleons keep the fraction of mass to delta-resonance mass
	// - nuclei lose the energy 1/A the interacting nucleon is carrying
	double relativeEnergyLoss = (A == 1) ? 1 - 938. / 1232. : 1. / A;
	lossRate *= relativeEnergyLoss;

	// scaling factor: interaction rate --> energy loss rate
	lossRate *= (1 + z);

	return 1. / lossRate;
}

std::vector<double> PhotoPionProductionMod::sophiaEvent(bool onProton,  // proton or neutron
                                                    	double Ein,     // primary nucleon's energy
                                                    	double eps      // target photon's energy
                                                    	) const {
	int nature = 1 - int(onProton); // interacting particle: 0 for proton, 1 for neutron
	Ein /= GeV;
	eps /= GeV;
	double momentaList[5][2000]; // momentum list
	int particleList[2000]; // particle id list
	int nParticles; // number of outgoing particles

	sophiaeventmod_(nature, Ein, eps, momentaList, particleList, nParticles);

	std::vector<double> output;
	for (int i = 0; i < nParticles; ++i) {
		int id = 0;
		int pType = particleList[i];
		switch (pType) {
			case 13:  // proton
			case 14:  // neutron
				id = nucleusId(1, 14 - pType);
				break;
			case -13:  // anti-proton
			case -14:  // anti-neutron
				id = -nucleusId(1, 14 - pType);
				break;
			case 1:  // photon
				id = 22;
				break;
			case 2:  // positron
				id = -11;
				break;
			case 3:  // electron
				id = 11;
				break;
			case 15:  // nu_e
				id = 12;
				break;
			case 16:  // anti-nu_e
				id = -12;
				break;
			case 17:  // nu_mu
				id = 14;
				break;
			case 18:  // antu-nu_mu
				id = -14;
				break;
			default:
				throw std::runtime_error("PhotoPionProduction: unexpected particle " + kiss::str(pType));
		}
		output.push_back(id);
	}
	for (int i = 0; i < nParticles; ++i) {
		double Eout = momentaList[3][i] * GeV; // only the energy is used; could be changed for more detail
		output.push_back(Eout);
	}
	return output;
}

} // namespace crpropa