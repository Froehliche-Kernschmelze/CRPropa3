#include "crpropa/module/HadronicInteraction.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"
#include "crpropa/Random.h"

#include <fstream>
#include <limits>
#include <cmath>
#include <stdexcept>

namespace crpropa {


HadronicInteraction::HadronicInteraction(double massDensity) {
    this->massDensity = massDensity;
    setDescription("HadronicInteraction");
}

double HadronicInteraction::pionSpectrum(double x, double ePrimary) const {
    const double E0pi = 139. * MeV;
    if (x < E0pi / ePrimary) return 0.;
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

double HadronicInteraction::etaSpectrum(double x, double ePrimary) const {
    const double E0eta = 547. * MeV;
    if (x < E0eta / ePrimary) return 0.;
    const double term1 = 0.55 + 0.028 * std::log(x);
    const double term2 = std::sqrt(1. - E0eta  / (x * ePrimary));
    const double term3 = pionSpectrum(x, ePrimary);
    const double Feta = term1 * term2 * term3;
    return Feta;
}

int HadronicInteraction::samplePionNumber(double ePrimary) const {
    const double xMin = 1. / 1000.;
    const double xMax = 1.;
    const double stepSize = 1. / 1000.;
    double x = xMin;
    double y = 0.;
    int stepsDone = 0;
    do {
        y += pionSpectrum(x, ePrimary);
        x += stepSize;
        stepsDone++;
    } while (x < xMax);
    return round(y / stepsDone * (x - xMin));
}

int HadronicInteraction::sampleEtaNumber(double ePrimary) const {
    const double xMin = 1. / 1000.;
    const double xMax = 1.;
    const double stepSize = 1. / 1000.;
    double x = xMin;
    double y = 0.;
    int stepsDone = 0;
    do {
        y += etaSpectrum(x, ePrimary);
        x += stepSize;
        stepsDone++;
    } while (x < xMax);
    return round(y / stepsDone * (x - xMin));
}

double HadronicInteraction::samplePionEnergy(double ePrimary) const {
    double Fmax = 0.;
    const double xMin = 1. / 1000.;
    const double xMax = 1.;
    const double stepSize = 1. / 1000.;

    double x = xMin;
    while (x < xMax) {
        double F = pionSpectrum(x, ePrimary);
        if (F > Fmax)
            Fmax = F;
        x += stepSize;
    }

    Random &random = Random::instance();
    double F = 0.;
    do {
        x = std::pow(10, -3 * random.rand());
        F = pionSpectrum(x, ePrimary);
    } while(F < random.rand() * Fmax);
    return x * ePrimary;
}

double HadronicInteraction::sampleEtaEnergy(double ePrimary) const {
    double Fmax = 0.;
    const double xMin = 1. / 1000.;
    const double xMax = 1.;
    const double stepSize = 1. / 1000.;

    double x = xMin;
    while (x < xMax) {
        double F = etaSpectrum(x, ePrimary);
        if (F > Fmax)
            Fmax = F;
        x += stepSize;
    }

    Random &random = Random::instance();
    double F = 0.;
    do {
        x = std::pow(10, -3 * random.rand());
        F = etaSpectrum(x, ePrimary);
    } while(F < random.rand() * Fmax);
    return x * ePrimary;
}

double HadronicInteraction::xSectionKelner06(double ePrimary) const {
    const double L = std::log(ePrimary / TeV);
    const double A = 1 - std::pow(1.22 * 1e-3 * TeV / ePrimary, 4);
    return (34.3 + 1.88 * L + 0.25 * L * L) * A * A * 1e-31;
}

void HadronicInteraction::process(Candidate *candidate) const {
    const double stepLength = candidate->getCurrentStep();
    const double ePrimary = candidate->current.getEnergy();
    const double xSection = xSectionKelner06(ePrimary);
    const double interactionProbability = xSection * stepLength * massDensity;

    const double limit = 1 / interactionProbability * 0.1;
    if (stepLength > limit)
        candidate->limitNextStep(limit);

    if (ePrimary < 1. * GeV)
        return;

    Random &random = Random::instance();
    if (random.rand() > interactionProbability)
        return;

    performInteraction(candidate);
}

void HadronicInteraction::performInteraction(Candidate *candidate) const {
    std::vector<int> outPartID;
    std::vector<double> outPartE;
    double eAvailable = candidate->current.getEnergy();

    Random &random = Random::instance();
    const double startPiont = random.rand();
    bool doPi0 = (startPiont >= 0. && startPiont < 0.25);
    bool doPiCharged = (startPiont >= 0.25 && startPiont < 0.75);
    bool doEta = (startPiont >= 0.75 && startPiont <= 1.);
    bool donePi0 = false;
    bool donePiCharged = false;
    bool doneEta = false;

    do {
        if (doPi0 && !donePi0) {
            int nPi0 = samplePionNumber(eAvailable);
            // std::cout << "nPi0 = " << nPi0 << std::endl;
            for (int i = 0; i < nPi0; ++i) {
                const double ePi0 = samplePionEnergy(eAvailable);
                // std::cout << "ePi0 = " << ePi0 / GeV << std::endl;
                if (eAvailable >= ePi0) {
                    outPartID.push_back(111);
                    outPartE.push_back(ePi0);
                    eAvailable -= ePi0;
                } else {
                    break;
                }
            }
            doPiCharged = true;
            doEta = true;
            donePi0 = true;
        }

        if (doPiCharged && !donePiCharged) {
            int nPiCharged = 2 * samplePionNumber(eAvailable);
            // std::cout << "nPi+/- = " << nPiCharged << std::endl;
            for (int i = 0; i < nPiCharged / 2; ++i) {
                const double ePiPlus = samplePionEnergy(eAvailable);
                const double ePiMinus = samplePionEnergy(eAvailable);
                // std::cout << "ePi+ = " << ePiPlus / GeV << std::endl;    
                // std::cout << "ePi- = " << ePiMinus / GeV << std::endl;
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
            doPi0 = true;
            doEta = true;
            donePiCharged = true;
        }

        if (doEta && !doneEta) {
            int nEta = sampleEtaNumber(eAvailable);
            // std::cout << "nEta = " << nEta << std::endl;
            for (int i = 0; i < nEta; ++i) {
                const double eEta = sampleEtaEnergy(eAvailable);
                // std::cout << "eEta = " << eEta / GeV << std::endl;    
                if (eAvailable >= eEta) {
                    outPartID.push_back(221);
                    outPartE.push_back(eEta);
                    eAvailable -= eEta;
                } else {
                    break;
                }
            }
            doPi0 = true;
            doPiCharged = true;
            doneEta = true;
        }
    } while (!donePi0 || !donePiCharged || !doneEta);

    candidate->current.setEnergy(eAvailable);
    std::cout << "ePrimary = " << eAvailable / GeV << std::endl;
    const Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
    for (int i = 0; i < outPartID.size(); ++i) {
        candidate->addSecondary(outPartID[i], outPartE[i], pos);
        std::cout << outPartID[i] << " " << outPartE[i] / GeV << std::endl;
    }
}

} // namespace CRPropa
