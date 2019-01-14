#include "crpropa/PhotonBackground.h"
#include "crpropa/Common.h"
#include "crpropa/Random.h"
#include "crpropa/Units.h"

#include <vector>
#include <fstream>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <sstream>

namespace crpropa {

// Class to handle global evolution of IRB models (cf. CRPropa3-data/calc_scaling.py)
struct PhotonFieldScaling {
	bool initialized;
	std::string name;
	std::vector<double> tab_z;
	std::vector<double> tab_s;

	PhotonFieldScaling(std::string filename) {
		name = filename;
		initialized = false;
	}

	void init() {
		std::string path = getDataPath("Scaling/scaling_" + name + ".txt");
		std::ifstream infile(path.c_str());

		if (!infile.good())
			throw std::runtime_error(
					"crpropa: could not open file scaling_" + name);

		double z, s;
		while (infile.good()) {
			if (infile.peek() != '#') {
				infile >> z >> s;
				tab_z.push_back(z);
				tab_s.push_back(s);
			}
			infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		}

		infile.close();
		initialized = true;
	}

	double scalingFactor(double z) {
		if (!initialized) 
#pragma omp critical(init)
			init();
		if (z > tab_z.back())
			return 0;  // zero photon background beyond maximum tabulated value
		return interpolate(z, tab_z, tab_s);
	}
};

static PhotonFieldScaling scalingKneiske04("IRB_Kneiske04");
static PhotonFieldScaling scalingStecker05("IRB_Stecker05");
static PhotonFieldScaling scalingFranceschini08("IRB_Franceschini08");
static PhotonFieldScaling scalingFinke10("IRB_Finke10");
static PhotonFieldScaling scalingDominguez11("IRB_Dominguez11");
static PhotonFieldScaling scalingGilmore12("IRB_Gilmore12");
static PhotonFieldScaling scalingStecker16_upper("IRB_Stecker16_upper");
static PhotonFieldScaling scalingStecker16_lower("IRB_Stecker16_lower");

double photonFieldScaling(PhotonField photonField, double z) {
    switch (photonField) {
    case CMB: // constant comoving photon number density
    case PF1: case PF2: case PF3: case PF4:
    case PF5: case PF6: case PF7: case PF8:
        return 1;
    case IRB:
    case IRB_Kneiske04:
        return scalingKneiske04.scalingFactor(z);
    case IRB_Stecker05:
        return scalingStecker05.scalingFactor(z);
    case IRB_Franceschini08:
        return scalingFranceschini08.scalingFactor(z);
    case IRB_Finke10:
        return scalingFinke10.scalingFactor(z);
    case IRB_Dominguez11:
        return scalingDominguez11.scalingFactor(z);
    case IRB_Gilmore12:
        return scalingGilmore12.scalingFactor(z);
    case IRB_Stecker16_upper:
        return scalingStecker16_upper.scalingFactor(z);
    case IRB_Stecker16_lower:
        return scalingStecker16_lower.scalingFactor(z);
    case URB_Protheroe96:
        if (z < 0.8) { return 1; }
        if (z < 6) { return pow((1 + 0.8) / (1 + z), 4); }
        else { return 0; }
    default:
        throw std::runtime_error("PhotonField: unknown photon background");
    }
}

std::string photonFieldName(PhotonField photonField) {
    switch (photonField) {
        case CMB: return "CMB";
        case PF1: return "PF1";
        case PF2: return "PF2";
        case PF3: return "PF3";
        case PF4: return "PF4";
        case PF5: return "PF5";
        case PF6: return "PF6";
        case PF7: return "PF7";
        case PF8: return "PF8";
        case IRB:
        case IRB_Kneiske04: return "IRB_Kneiske04";
        case IRB_Stecker05: return "IRB_Stecker05";
        case IRB_Franceschini08: return "IRB_Franceschini08";
        case IRB_Finke10: return "IRB_Finke10";
        case IRB_Dominguez11: return "IRB_Dominguez11";
        case IRB_Gilmore12: return "IRB_Gilmore12";
        case IRB_Stecker16_upper: return "IRB_Stecker16_upper";
        case IRB_Stecker16_lower: return "IRB_Stecker16_lower";
        case URB_Protheroe96: return "URB_Protheroe96";
        default:
            throw std::runtime_error("PhotonField: unknown photon background");
    }
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// custom photon field methods related to SAMPLING
// These methods are taken from SOPHIA.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Photon_Field::Photon_Field(std::string fieldPath) {
    init(fieldPath);
}


Photon_Field::Photon_Field() {
    // empty constructor for initialization in some modules
}


double Photon_Field::sample_eps(double z_in) const {
/*
    - input: particle type 0=neutron, 1=proton, its energy [GeV], its redshift
    - output: photon energy [eV] of random photon of photon field
    - samples distribution of n(epsilon)/epsilon^2
*/ 
    const double z_min = redshift[0];
    const double z_max = redshift[redshift.size()-1];
    if ( (z_in < z_min) || (z_in > z_max) )
        return 0.;

    // find closest redshift position
    int z_pos;
    double smallestDiff = z_max;
    for (int i = 0; i < redshift.size(); ++i) {
        double diff = std::abs(z_in-redshift[i]);
        if (diff < smallestDiff) {
            smallestDiff = diff;
            z_pos = i;
        }
    }

    const double epsMin = energy[0];
    const double epsMax = energy[energy.size()-1];

    // assume probability propto density
    double sumDens = 0.;
    double maxDens = 0.;
    int maxDensPos;
    for (int i = 0; i < energy.size()-1; ++i) {
        double dens = get_photonDensity(energy[i], z_pos);
        sumDens += dens;
        if (dens > maxDens) {
            maxDens = dens;
            maxDensPos = i;
        }
    }
    const double pMax = get_photonDensity(energy[maxDensPos], z_pos)/sumDens;

    // sample eps randomly between epsMin...epsMax
    Random &random = Random::instance();
    double eps, peps;
    do {
        eps = epsMin+random.rand()*(epsMax-epsMin);
        peps = get_photonDensity(eps, z_pos)/sumDens;
    } while (random.rand()*pMax > peps);
    return eps;
}


void Photon_Field::init(std::string filename) {
    std::ifstream infile(filename.c_str());
    if (!infile.good()) {
        throw std::runtime_error("PhotoPionProduction @ Photon_Field::init : could not open file " + filename);
    }
    std::string line;
    int i = 0;
    while ( std::getline(infile, line) ) {
        if (line.find('#') == 0)
            continue;
        std::istringstream ss(line);
        std::vector<double> vec;
        double n;
        while (ss >> n)
            vec.push_back(n);
        if (i == 0) {
            energy = vec;
            i++;
            continue;
        }
        if (i == 1) {
            redshift = vec;
            i++;
            continue;
        }
        dn_deps.push_back(vec);
    }
    infile.close();
}


double Photon_Field::get_photonDensity(double eps, int z_pos) const {
/*
    - input: photon energy [eV], redshift
    - output: dn_deps(e,z) [#/(eV cm^3)] from input file
    - called by: sample_eps
*/
    // find closest dn_deps
    double smallestDiff = energy[energy.size()-1];
    int closestPos;
    for (int i = 0; i < energy.size(); ++i) {
        double diff = std::abs(eps-energy[i]);
        if (diff < smallestDiff) {
            smallestDiff = diff;
            closestPos = i;
        }
    }
    // linear interpolation of energy
    double realDiff = eps-energy[closestPos];
    double rho;
    if (realDiff >= 0.) {
        rho = realDiff / (energy[closestPos+1] - energy[closestPos]);
        return (1.-rho)*dn_deps[closestPos][z_pos]
               + rho*dn_deps[closestPos+1][z_pos];
    } else {
        rho = 1. - (std::abs(realDiff)/energy[closestPos-1]);
        return (1.-rho)*dn_deps[closestPos-1][z_pos]
               + rho*dn_deps[closestPos][z_pos];
    }
}

} // namespace crpropa
