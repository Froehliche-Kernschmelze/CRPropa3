#include "crpropa/module/HadronicInteraction.h"
#include "crpropa/massDistribution/CMZDensity.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"
#include "crpropa/Random.h"


#include <fstream>
#include <limits>
#include <cmath>
#include <stdexcept>

namespace crpropa {
    
    
    HadronicInteraction::HadronicInteraction(bool electrons, bool photons, bool neutrinos) {
        haveElectrons = electrons;
		havePhotons = photons;
		haveNeutrinos = neutrinos;
        setDescription("HadronicInteraction");
    }
    

void HadronicInteraction::setHaveElectrons(bool b) {
	haveElectrons = b;
}

void HadronicInteraction::setHavePhotons(bool b) {
	havePhotons = b;
}

void HadronicInteraction::setHaveNeutrinos(bool b) {
	haveNeutrinos = b;
}



    double HadronicInteraction::getDensity(const Vector3d &position) const {
        CMZDensity obj;
        double density=obj.getDensity(position);
        return density;
    }
    Vector3d HadronicInteraction::Position(double height, double radius) const{
        Random &random = Random::instance();
        int i=0;
        Vector3d pos(0,0,0);
        double phi=random.rand()*2*M_PI;
        int j=0;
        do{
            double r = random.rand()* radius;
            double yr = random.rand();
            double Fr=exp(-r*r/(2*4.2*4.2*kpc*kpc));
            
            if (yr < Fr){
                pos=Vector3d(cos(phi)*r, sin(phi)*r, 0);
                j++;
            }
        }while(j==0);
        do{
            double z = random.rand()* height;
            double yz = random.rand();
            double Fz=exp(-z/(10*pc));
            if (yz < Fz){
                double a = random.rand();
                if (a <= 0.5){
                    z=-z;
                }
                pos= pos + Vector3d(0,0,z);
                j++;
                
            }
        }while(j==1);
        return pos;
    }

    
    //Distribution function (energy) for electrons, electron neutrinso and (second) muon neutrinos based on Kelner 2006
    
    double HadronicInteraction::distribution_e(double energy, double x) const{
        
        double L=log(energy / TeV);
        double Be= 1/(69.5+2.65*L+0.3*pow(L,2.));
        double betae=1/pow((0.201+0.062*L+0.00042*pow(L,2.)), 0.25);
        double ke=(0.279 + 0.141 *L + 0.0172* pow(L, 2.))/(0.3+ pow((2.3+L), 2.));
        double F=Be*pow((1+ke*pow(log(x),2.)), 3.) /(x*(1+0.3/pow(x, betae)))*(pow(-log(x), 5.));
        
        return F;
    }
    
    //Number of electrons, electron neutrinos and (second) muons neutrinos produced in a given interaction  based on Kelner 2006
    double HadronicInteraction::number_e(double energy) const{
        double x=1/100000.;
        double i=1/100000.;
        double y=0;
        double j=0;
        do{
            y=y+distribution_e(energy,x);
            x=x+i;
            j++;
        }while(x < 1);
        return y/j*(x-1/1000.);
    }
    
    //Distribution function (energy) for (first) muon neutrino based on Kelner 2006
    double HadronicInteraction::distribution_my1(double energy, double x) const{
        double L=log(energy / TeV);
        double Bm= 1.75+0.204*L+0.01 * pow(L,2.);
        double betam=1/(1.67+0.111*L+0.0038*pow(L,2.));
        double km=1.07-0.086*L+0.002*pow(L,2.);
        x=x/0.427;
        double aa=(1-pow(x,betam))/(1+km*pow(x, betam)*(1-pow(x,betam)));
        double A=Bm*log(x)/x*pow(aa, 4.);
        double B=1/log(x)-4*betam*pow(x,betam)/(1- pow(x,betam))-4*km*betam*pow(x, betam)*(1-2*pow(x,betam))/(1+km*pow(x,betam)*(1-pow(x,betam)));
        double F=A*B;
        
        return F;
    }
    
    
    //Number of (first) muon neutrinos produced in a given interaction  based on Kelner 2006
    
    double HadronicInteraction::number_my1(double energy) const{
        double x=1/(100000.);
        double i=1/100000.;
        double y=0;
        double j=0;
        do{
            y=y+distribution_my1(energy,x);
            x=x+i;
            j++;
        }while(x < 0.427);
        return y/j*(x-1/1000.);
    }
    
    //Distribution function (energy) for gamma rays based on Kelner 2006
    double HadronicInteraction::distribution_gamma(double energy, double x) const{
        double L=log(energy / TeV);
        double Bg=1.3+0.14*L+0.011*L*L;
        double betag=1/(1.79+0.11*L+0.008*L*L);
        double kg=1/(0.801+0.049 *L+0.014 *L*L);
        double A=Bg*log(x)/x;
        double B=(1-pow(x, betag))/(1+kg*pow(x, betag)*(1-pow(x, betag)));
        double C=1/log(x)-4*betag*pow(x, betag)/(1-pow(x, betag))-4*kg*betag*pow(x, betag)*(1-2*pow(x, betag))/(1+kg*pow(x, betag)*(1-pow(x, betag)));
        double F=A*pow(B, 4.)*C;
        
        return F;
    }
    
    //Number of gamma rays produced in a given interaction  based on Kelner 2006
    double HadronicInteraction::number_gamma(double energy) const{
        double x=1/100000.;
        double i=1/100000.;
        double y=0;
        double j=0;
        do{
            y=y+distribution_gamma(energy,x);
            x=x+i;
            j++;
        }while(x < 1);
        return y/j*(x-1/1000.);
    }
    
    //Energy distribution for lepton secondaries of pp interactions based on Carceller 2017
    double HadronicInteraction::distribution_Carceller(double energy, double x, double jcap, double a0, double b0) const{
        double a = a0 * (1+ 0.073*log(energy/ PeV)+ 0.0070*log(energy/PeV)*log(energy/PeV));
        double b = b0 * (1+0.020*log(energy/PeV)+0.0018*log(energy/PeV)*log(energy/PeV));
        double A = a * pow((1-jcap*x), 3.)/x;
        double B = exp(-b*pow(jcap*x, 0.43))/pow(1+pow(0.1*GeV/(x*energy), 0.5), 2.);
        double F=A*B;
        
        return F;
    }
    
    
    //Energy distribution for gamma photons based on Carceller 2017
    double HadronicInteraction::distribution_Carceller_g(double energy, double x, double jcap, double a0, double b0) const{
        double a = a0 * (1+ 0.073*log(energy/ PeV)+ 0.0070*log(energy/PeV)*log(energy/PeV));
        double b = b0 * (1+0.02*log(energy/PeV)+0.0018*log(energy/PeV)*log(energy/PeV));
        double A = a * pow((1-jcap*x), 3.)/x;
        double B = exp(-b*pow(jcap*x, 0.43))/pow(1+pow(0.2*GeV/(x*energy), 0.5), 2.);
        double F=A*B;
        
        return F;
    }
    
    //Cross Section of inelastic pp interaction based on Tan & Ng 1983 (Used in Galprop)
    double HadronicInteraction::CrossSection_Galprop(double energy) const{
        double cs_inel;
        double U = log(energy/ GeV * 1/200);
        if (U >= 0 and energy >= 3 * GeV){
            cs_inel=(32.2 * (1+0.0273*U))*1e-31+32.2*0.01*pow(U,2.)*1e-31;
        }
        if (U < 0 and energy >= 3 * GeV){
            cs_inel=(32.2 * (1+0.0273*U))*1e-31;
        }
        if (energy <= 0.3 * GeV){
            cs_inel = 0;
        }
        return cs_inel;
    }
    
    //Cross Section of inelastic pp interaction based on Kelner 2006
    double HadronicInteraction::CrossSection_Kelner(double energy) const{
        double L=log(energy / TeV);
        double A=1-pow(1.22*1e-3*TeV/energy, 4.);
        double cs_inel=(34.3 + 1.88*L+0.25 *L*L)*A*A*1e-31;
        return cs_inel;
    }
    
    //Cross Section of inelastic pp interaction based on Carceller 2017
    double HadronicInteraction::CrossSection_Carceller(double energy) const{
        double cs_inel=17.7*pow(energy/GeV, 0.082)*1e-31;
        return cs_inel;
    }
    
    
    
    //Process Function:
    
    void HadronicInteraction::process(Candidate *candidate) const {
        
        //Get candidate's properties
        double step = candidate->getCurrentStep();
        double energy = candidate->current.getEnergy();
        double id = candidate->current.getId();
        Vector3d pos= candidate->current.getPosition();
        
        //No interaction for secondaries
        if (id != 1000010010){
            return;
        }
        
        
        double cs_inel=0;
 
        
        if (id == 1000010010) {
            cs_inel=CrossSection_Kelner(energy);
        }
        
        
        Random &random = Random::instance();
       
        CMZDensity obj;
        double density=obj.getDensity(pos);
        //double density=pow(10.,15.);

        //Probability of interaction
        double p_pp=cs_inel*density*step;
        double ra = random.rand();
        
        // limit next step to mean free path
        double limit = 1 / p_pp*0.1;
        
        if (step > limit) {
            candidate->limitNextStep(limit);
        }
        
        //Interaction?
        if (ra > p_pp or energy < 1*GeV){
            return;
        }

        //Initialize energies of secondaries
        double Eout=0;
        double Emo=0;
        double Ee=0;
        double Ene=0;
        double Emt=0;
        double Eg=0;
        double Etot=0;
        double Econ=0;
        
        double gamma=0;
        double j=0;

        // Establish number of secondaries
        
            
   //Gammas NG
        double FG=number_gamma(energy);
        double NG=round(FG);
        
        //First myon neutrino Nmy1
        double Fmy1=number_my1(energy);
        double Nmy1=round(Fmy1);
        
        
        //Electron NE
        double FE=number_e(energy);
        double NE=round(FE);
        
        //Electron Neutrino NEN
        double NEN=NE;
    
        
        //Second Myon Neutrino Nmy2
        double Nmy2=NE;
    
        
        //Total number of secondaries in the interaction
        double N_tot=NE+NG+Nmy1+NEN+Nmy2;
        

        //Initialize for stopping criteria
        double i=1;
        double iG=1;
        double imy1=1;
        double ie=1;
        double ien=1;
        double imy2=1;
        double test;
        double threshold=0.0001;
        double end=-1;
        
        //Jump to random starting point:
        double sp=random.rand();
        if (sp <= 0.2)
            goto label1;
        if (sp <= 0.4)
            goto label2;
        if (sp <= 0.6)
            goto label3;
        if (sp <= 0.8)
            goto label4;
        if (sp <= 1.0)
            goto label5;

        do{
            //~ Gamma rays
            label1:
            test=iG;

            //Check if all gamma rays were created
            if (iG <= NG){
				if(end==-1)
                {
                //pick gamma ray's energy
                do{
                double x=threshold+random.rand()*(1-Eout/energy-threshold);
                Eout=x*energy;
                double E=distribution_gamma(energy, x);
                double Emax=distribution_gamma(energy, threshold);
                double y=random.rand()*Emax;
                    
                    if (y < E and (Etot+Eout)<energy )
                    {
                        candidate->addSecondary(22, Eout, pos);
                        Eg=Eg+Eout;
                        Etot=Etot+Eout;
                        i++;
                        iG++;
                        if(Etot/energy>=(1-threshold)){
                            end=i;
                        }
                    }
                    
                }while(test == iG);
			}else {
                    Eout=(energy-Etot)/(N_tot-end);
                    candidate->addSecondary(22, Eout, pos);
                    Eg=Eg+Eout;
                    i++;
                    iG++;
                }
			
            }
            label2:

            //~ First myon neutrino 14
            test = imy1;
			
            if (imy1 <= Nmy1){
				if(end==-1){
                do{
                
                double x=threshold+random.rand()*(0.427-threshold);
                Eout=x*energy;
                double E=distribution_my1(energy, x);
                double Emax=distribution_my1(energy, threshold);
                double y=random.rand()*Emax;
                    
                    if (y < E and (Etot+Eout)<energy)
                    {
                        candidate->addSecondary(14, Eout, pos);
                        Emo=Emo+Eout;
                        Etot=Etot+Eout;
                        i++;
                        imy1++;
                        if(Etot/energy>=(1-threshold)){
                            end=i;
                        }
                    }
                    
                }while(test == imy1);
            }else {
                    Eout=(energy-Etot)/(N_tot-end);
                    candidate->addSecondary(14, Eout, pos);
                    Emo=Emo+Eout;
                    i++;
                    imy1++;
                }
		}
			label3:
            //~ Electron 11
            test = ie;
			
            if (ie <= NE){
				
				if (end==-1){
                do{
                double x=threshold+random.rand()*(1-Eout/energy-threshold);
                Eout=x*energy;
                double E=distribution_e(energy, x);
                double Emax=distribution_e(energy, threshold);
                double y=random.rand()*Emax;
                
                    if (y < E and (Etot+Eout)<energy )
                    {
                        candidate->addSecondary(11, Eout, pos);
                        Etot=Etot+Eout;
                        Ee=Ee+Eout;
                        i++;
                        ie++;
                        if(Etot/energy>=(1-threshold)){
                            end=i;
                        }
                    }
                    
                }while(test == ie);
            }else {
                    Eout=(energy-Etot)/(N_tot-end);
                    candidate->addSecondary(11, Eout, pos);
                    i++;
                    ie++;
                    Ee=Ee+Eout;
                }
		}
            
            //~ Electron neutrino 12
            label4:
            test=ien;

            if (ien <= NEN ){
				if (end==-1){
                do{
                double x=threshold+random.rand()*(1-Eout/energy-threshold);
                Eout= x*energy;
                double E=distribution_e(energy, x);
                double Emax=distribution_e(energy, threshold);
                double y=random.rand()*Emax;
                    
                    if (y < E and (Etot+Eout)<energy)
                    {
                        candidate->addSecondary(12, Eout, pos);
                        Ene=Ene+Eout;
                        Etot=Etot+Eout;
                        i++;
                        ien++;
                        if(Etot/energy>=(1-threshold)){
                            end=i;
                        }
                    }
                    
                }while(ien==test);
            }else {
                    Eout=(energy-Etot)/(N_tot-end);
                    candidate->addSecondary(12, Eout, pos);
                    i++;
                    ien++;
                    Ene=Ene+Eout;
                }
            }

            //~ Second myon neutrino 14
			label5:
            test=imy2;
            if(imy2 <= Nmy2){
				if (end==-1){
                do{
                    double x=threshold+random.rand()*(1-Eout/energy-threshold);
                    Eout= x*energy;
                    double E=distribution_e(energy, x);
                    double Emax=distribution_e(energy, threshold);
                    double y=random.rand()*Emax;

                        if (y < E and (Etot+Eout)<energy)
                        {
                            candidate->addSecondary(14, Eout, pos);
                            Emt=Emt+Eout;
                            Etot=Etot+Eout;
                            i++;
                            imy2++;
                            if(Etot/energy>=(1-threshold)){
                            end=i;
                        }
                        }
                    
                }while(imy2==test);
            }else {
                    Eout=(energy-Etot)/(N_tot-end);
                    candidate->addSecondary(14, Eout, pos);
                    i++;
                    imy2++;
                    Emt=Emt+Eout;
                }
            }
        }while (i <= N_tot);
    
		if(end != -1){
			std::cout<<0<<std::endl;
		}
			
      //Reduce primary's energy
        candidate->current.setEnergy(energy-(Ene+Emt+Ee+Emo+Eg));
        
        return;
    }
    
} //~ namespace CRPropa



