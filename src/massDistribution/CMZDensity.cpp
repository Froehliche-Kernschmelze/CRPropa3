#include "crpropa/massDistribution/CMZDensity.h"

namespace crpropa {
    
    CMZDensity::CMZDensity(){
		}
    
    Vector3d CMZDensity::CMZTrafo(const Vector3d &position) const {
        
        double xC = -50*pc;        //offset
        double yC = 50*pc;
        double ThettaC = 70/180*M_PI; //radiant
        
        Vector3d pos;
        pos.x = (position.x - xC)*cos(ThettaC) + (position.y -yC)*sin(ThettaC);
        pos.y = -(position.x -xC)*sin(ThettaC) + (position.y - yC)*cos(ThettaC);
        pos.z = position.z;
        return pos;
    }
    double CMZDensity::getDensity(const Vector3d &position) const {
        
        return getH2Density(position);
    }
    
    double CMZDensity::getH2Density(const Vector3d &position) const {
        
        //VALUES FOR n PER M^3
        double n=0;
        double pi=3.1415926535;
        double x=position.x;
        double y=position.y;
        double z=position.z;
        double cmm=1000000;
        
        ////SgrC//"The Galactic Center Molecular Cloud Survey" Kaufmann et al. 2016
        double cloud_mid_y=std::sin(359.45*pi/180)*8.5*kpc;
        double cloud_mid_x=0;
        double cloud_mid_z=std::sin(-0.11*pi/180)*8.5*kpc;
        double r=1.7*pc;
        

        if (pow((cloud_mid_z-z)*(cloud_mid_z-z)+(cloud_mid_x-x)*(cloud_mid_x-x)+(cloud_mid_y-y)*(cloud_mid_y-y), 0.5)<r){
                n=2*1.8*10000;
            
            return n*cmm;
        }
//        //20kms //already included in the inner 10pc
//
//        cloud_mid_x=0;
//        cloud_mid_y=std::sin(395.87*pi/180)*8.5*kpc;
//        cloud_mid_z=std::sin(-0.08*pi/180)*8.5*kpc;
//        r=5.1*pc;
//
//        if (pow((cloud_mid_z-z)*(cloud_mid_z-z)+(cloud_mid_x-x)*(cloud_mid_x-x)+(cloud_mid_y-y)*(cloud_mid_y-y), 0.5)<r){
//                n=2*0.9*1000;
//            return n;
//        }
//
//        //50kms// already included in the inner 10pc
//
//        cloud_mid_x=0;
//        cloud_mid_y=std::sin(395.98*pi/180)*8.5*kpc;
//        cloud_mid_z=std::sin(-0.07*pi/180)*8.5*kpc;
//        r=2.7*pc;
//
//        if (pow((cloud_mid_z-z)*(cloud_mid_z-z)+(cloud_mid_x-x)*(cloud_mid_x-x)+(cloud_mid_y-y)*(cloud_mid_y-y), 0.5)<r){
//            n=2*1.2*1000;
//            return n;
//        }
        
        //G0.253+0.016// from "MAGNETIC FIELDS IN HIGH-MASS INFRARED DARK CLOUDS" Pillai et al. 2015
        cloud_mid_x=0;
        cloud_mid_y=std::sin(0.25*pi/180)*8.5*kpc;
        cloud_mid_z=std::sin(-0.02*pi/180)*8.5*kpc;
        r=2.9*pc;
        if (pow((cloud_mid_z-z)*(cloud_mid_z-z)+(cloud_mid_x-x)*(cloud_mid_x-x)+(cloud_mid_y-y)*(cloud_mid_y-y), 0.5)<r){
            n=2*8*10000;
            return n*cmm;
        }
        
        //Sgr B1-off//"The Galactic Center Molecular Cloud Survey" Kaufmann et al. 2016
        
        cloud_mid_x=0;
        cloud_mid_y=std::sin(0.48*pi/180)*8.5*kpc;
        cloud_mid_z=std::sin(-0.00*pi/180)*8.5*kpc;
        r=3.6*pc;
        
        if (pow((cloud_mid_z-z)*(cloud_mid_z-z)+(cloud_mid_x-x)*(cloud_mid_x-x)+(cloud_mid_y-y)*(cloud_mid_y-y), 0.5)<r){
            n=2*1.1*10000;
            return n*cmm;
        }
        
        //Dust Ridge B//"Star formation in a high-pressure environment: an SMA view of the Galactic Centre dust ridge" Walker et al 2018 and the references therein
        cloud_mid_x=0;
        cloud_mid_y=std::sin(0.34*pi/180)*8.5*kpc;
        cloud_mid_z=std::sin(0.055*pi/180)*8.5*kpc;
        r=1.9*pc;
        if (pow((cloud_mid_z-z)*(cloud_mid_z-z)+(cloud_mid_x-x)*(cloud_mid_x-x)+(cloud_mid_y-y)*(cloud_mid_y-y), 0.5)<r){
            n=2*0.7*10000;
            return n*cmm;
       }
        //Dust Ridge C//"Star formation in a high-pressure environment: an SMA view of the Galactic Centre dust ridge" Walker et al 2018 and the references therein
        cloud_mid_x=0;
        cloud_mid_y=std::sin(0.38*pi/180)*8.5*kpc; //-16.3 pc
        cloud_mid_z=std::sin(0.05*pi/180)*8.5*kpc;
        r=1.9*pc;
        if (pow((cloud_mid_z-z)*(cloud_mid_z-z)+(cloud_mid_x-x)*(cloud_mid_x-x)+(cloud_mid_y-y)*(cloud_mid_y-y), 0.5)<r){
            n=2*0.9*10000;
            return n*cmm;
        }
        //Dust Ridge D//"Star formation in a high-pressure environment: an SMA view of the Galactic Centre dust ridge" Walker et al 2018 and the references therein
        cloud_mid_x=0;
        cloud_mid_y=std::sin(0.412*pi/180)*8.5*kpc; //-16.3 pc
        cloud_mid_z=std::sin(0.05*pi/180)*8.5*kpc;
        r=3.2*pc;
        if (pow((cloud_mid_z-z)*(cloud_mid_z-z)+(cloud_mid_x-x)*(cloud_mid_x-x)+(cloud_mid_y-y)*(cloud_mid_y-y), 0.5)<r){
            n=2*0.8*10000;
            return n*cmm;
        }
        //Dust Ridge E//"Star formation in a high-pressure environment: an SMA view of the Galactic Centre dust ridge" Walker et al 2018 and the references therein
        cloud_mid_x=0;
        cloud_mid_y=std::sin(0.478*pi/180)*8.5*kpc; //-16.3 pc
        cloud_mid_z=std::sin(0.005*pi/180)*8.5*kpc;
        r=2.4*pc;
        if (pow((cloud_mid_z-z)*(cloud_mid_z-z)+(cloud_mid_x-x)*(cloud_mid_x-x)+(cloud_mid_y-y)*(cloud_mid_y-y), 0.5)<r){
            n=2*2.8*10000;
            return n*cmm;
        }
        //Dust Ridge F//"Star formation in a high-pressure environment: an SMA view of the Galactic Centre dust ridge" Walker et al 2018 and the references therein
        cloud_mid_x=0;
        cloud_mid_y=std::sin(0.496*pi/180)*8.5*kpc; //-16.3 pc
        cloud_mid_z=std::sin(0.02*pi/180)*8.5*kpc;
        r=2.0*pc;
        if (pow((cloud_mid_z-z)*(cloud_mid_z-z)+(cloud_mid_x-x)*(cloud_mid_x-x)+(cloud_mid_y-y)*(cloud_mid_y-y), 0.5)<r){
            n=2*3.2*10000;
            return n*cmm;
        }
        //Sgr D//HII region//#from LaRosa et al (2000) "A WIDE-FIELD 90 CENTIMETER VLA IMAGE OF THE GALACTIC CENTER REGION"-> density of SgrC,Sgr D column density is 3/5 of the SgrC column Density--> from Lis et al (1991) "Millimeter continuum observations of Galactic center giant molecular cloud cores" Radius is calculated from this relation
        cloud_mid_x=0;
        cloud_mid_y=std::sin(1.12*pi/180)*8.5*kpc; //-16.3 pc
        cloud_mid_z=std::sin(-0.07*pi/180)*8.5*kpc;
        r=1.8*pc;
        if (pow((cloud_mid_z-z)*(cloud_mid_z-z)+(cloud_mid_x-x)*(cloud_mid_x-x)+(cloud_mid_y-y)*(cloud_mid_y-y), 0.5)<r){
            n=2*1.08*10000;
            return n*cmm;
        }
        //Sgr B2//"High angular resolution submillimeter bservations of Sagittarius B2" Goldsmith et al 1990
       cloud_mid_x=0;
       cloud_mid_y=std::sin(0.66*pi/180)*8.5*kpc;
        cloud_mid_z=std::sin(-0.04*pi/180)*8.5*kpc;
        double r1=1.25*pc;
        double r2=22.5*pc;
       double RR=pow((cloud_mid_z-z)*(cloud_mid_z-z)+(cloud_mid_x-x)*(cloud_mid_x-x)+(cloud_mid_y-y)*(cloud_mid_y-y), 0.5);
        if (RR<r1){
           n=5.5*10000 + 2.2*1000;
           return n*cmm;
           }
        if  (r2>=RR>=r1) {
            n = 5.5*10000 *(r1/RR)*(r1/RR) +2.2*1000;
            return n*cmm;
           }
       
        
        double mSun=1.988e30; // in kg
        double m_p=1.672e-27; //proton mas in kg
        //pc=3.08567758e18##pc in cm
        
//Distribution from "Interstellar gas within 10pc of SgrA*" by Katia Ferriere 2012
        x=position.x/pc;
        y=position.y/pc;
        z=position.z/pc;
    
//Central Cavity
        double a=2.9/2;
        double b=a;
        double c=2.1/2;
    

//Circumnuclear ring
        r1=1.2;
        r2=3.0;
        double h1=0.4;
        double h2=1.0;
 
        double m=(h2-h1) / (r2-r1);

//SgrA East
        double cen1=-2.0;
        double cen2= 1.2;
        double cen3=-1.5;
        double a1=9.0/2;
        double b1=9.0/2;
        double c1=6.7/2;
        double x1=x-cen1;
        double y1=y-cen2;
        double z1=z-cen3;

//Radio Halo
        double rr1=0.0;
        double rr2=9.0;

//SC M-0.13-0.08
        double cen4=8.0;
        double cen5=-11.0;
        double cen6=-5.0;
        double x2=x-cen4;
        double y2=y-cen5;
        double z2=z-cen6;
        double a2=7.5/2;
        double b2=15./2;
        double c2=7.5/2;

//EC M-0.02-0.07
        double cen7=-3.0;
        double cen8=7.0;
        double cen9=-4.5;
        double x3=x-cen7;
        double y3=y-cen8;
        double z3=z-cen9;
        double rrr1=0.0;
        double rrr2=4.5;

//MR Molecular ridge
        //Beginning of the bridge
        double cen10=-3.0;
        double cen11=2.5;
        double cen12=-4.5;
//end of the bridge, here the y component is important as in this case the z component is irrelevant since z(EC)=z(SC), and x component depends on y.
        double cen102=4.8;
        double cen112=-6.0;
        double cen122=-5.0;
        double x4=x-cen10; //new coordinates
        double y4=y-cen11; //new coordinates
        double z4=z-cen12; //new coordinates

// SS Southern Streamer
//beginning part of the bridge with CNR
        double cen13=0.0;
        double cen14=-3.0;
        double cen15=0.0;
        
        double cen132=7.0;
        double cen142=-6.0;
        double cen152=-5.0;
        
        double x5=x-cen13;
        double y5=y-cen14;
        double z5=z-cen15;

//Western Streamer and Northern ridge

        double cen16=-2.0;
        double cen17=1.2;
        double cen18=-1.5;
        double a5=9.0/2;
        double b5=9.0/2;
        double c5=6.7/2;
        //????????
        x5=x-cen16;
        y5=y-cen17;
        z5=z-cen18;
        double V, M;

//Central Cavity, Form=Ellipsoid, Centered at SgrA*
        if (x*x /(a*a) + y*y /(b*b) + z*z /(c*c) <=1.)
        {
            V=4/3 *pi*a*b*c *pc*pc*pc;
            M=362.12*mSun;
            n+=M/V *1/m_p;
        }
        
//Circumnuclear ring, Form= Trapetoidal ring, Centered at SgrA*
        if (r1*r1 <= (x*x +y*y) and (x*x +y*y) <= r2*r2 and (-m/2*(pow(x*x+y*y, 0.5) -r1) -h1/2) <= z and z <= (m/2*(pow(x*x+y*y, 0.5)-r1) +h1/2))
        {
            
            V=pi/2 *(r2*r2 -r1*r1)*(h1+h2)*pc*pc*pc;
            M=2.013e5 *mSun;
            n+= M/V *1/m_p;
        }
        
//SgrA East, Form=Ellipsoid, Centered at (-2.0,1.2,-1.5)
        if ( x1*x1 /(a1*a1) + y1*y1 /(b1*b1) + (z1*z1) /(c1*c1) <1.)
        {
            V=4/3 *pi*a1*b1*c1*pc*pc*pc;
            M=19*mSun;
            n+= 3.*cmm; //M/V *1/m_p;
        }
        
//Radio Halo, Form=Sphere, Centered at (-2.0,1.2,-1.5)
        if (rr1*rr1<= x1*x1 +y1*y1 +z1*z1  and x1*x1 +y1*y1 +z1*z1 <=rr2*rr2){
            V=4/3 *pi*rr2*rr2*rr2 *pc*pc*pc;
            M= 1.3*10000 *mSun;
            n+= M/V *1/m_p;
        }
        
//SC M-0.13-0.08, Form=Ellipsoid, Centered at (4-12,-11,-5)
        if (x2*x2 /(a2*a2) + y2*y2 /(b2*b2) + z2*z2 /(c2*c2) <=1.){
            V=4/3*pi * a2*b2*c2*pc*pc*pc;
            M= 2.2e5 *mSun;
            n+= M/V *1/m_p;
        }
//EC M-0.02-0.07, Form=Sphere, Centered at (-3,7,-4.5)
        if (rrr1*rrr1<=x3*x3 +y3*y3 +z3*z3 and x3*x3 +y3*y3 +z3*z3 <=rrr2*rrr2){
            V=4/3*pi*rrr2*rrr2*rrr2 *pc*pc*pc;
            M= 1.9e5 *mSun;
            n+= M/V *1/m_p;
        }
//MR Molecular ridge, Form=Curved Cylinder of radius 1pc and length of 12 pc# curved cylinder which connects the southern part of EC with north eastern part of SC
        if (0+x4 <=pow(y4*y4 + z4*z4, 0.5) and pow(y4*y4 + z4*z4, 0.5) <=1.+x4 and cen112-3.0<=y4 and y4 <=0. and -1.<=z4 and z4<=1.){ //the cylinder changes just with x and y is symmetric around z--see coordinates# first if condition: because of the circle of radius 1 which is shifted with x4# second if condition: due to the variation scale of y from begining of EC to the end of SC# third if condition: due to the maximum value which z can take if y is zero (as circle has the radius 1).
            V=pi*1* 12*pc*pc*pc;
            M= 6.6e4*mSun;
            n+= 3.e4*cmm; //M/V *1/m_p;
        }
//SS Southern Streamer## curved cylinder of radius 1 pc and length of 10.5pc which connects northern end part of SC to southeastern part of CNR ##see comments above in MR Molecular ridge
        if (0+x5<=pow(y5*y5 + z5*z5, 0.5) and pow(y5*y5 + z5*z5, 0.5)<=1.+x5 and cen152 <= z and z <=cen15 and cen142<=y and y<=cen14) { //the cylinder changes with x,y and z-> has no symmetry #first of condition: circle of radius 1 which is shifted by x5# second if condition: variation scale of z# third if condition: variation scale of y
            V=pi*1* 10.5*pc*pc*pc;
            M= 2.2e4*mSun;
            n+= 3.e4*cmm; //M/V *1/m_p;
            }
        
            double theta=50.;
//WS Western Streamer #curved cylinder of radius 0.5 pc and length of 8 pc#along western body of SgrA East
        if (x5>=0 and std::abs(y5/b5)<1.){//because we are interested in the western part otherwise we would obtain two times the same value which are axissymmetric
            theta=std::asin(y5/b5);
//phi=np.arcsin(x1/a1 *1/np.cos(theta))
            if (1<=x5*x5 /(a5*a5) + y5*y5 /(b5*b5) + z5*z5 /(c5*c5) and x5*x5 /(a5*a5) + y5*y5 /(b5*b5) + z5*z5 /(c5*c5) <=1.5 and -0.5<=z5 and z5 <=0.5  and -5/12*pi<=theta<=2/12 *pi){  //##first if condition: ellipsoid ring with 1pc diameter# second if condition: cylinder changing in direction of x but not in z direction and cylinder surface is also given in z direction ->diameter of 1pc# third if condition: we need only a part of the ellipsoid ring which is also a "curved cylinder" ->length of 8pc
                V=pi *0.5*0.5 *8*pc*pc*pc;
                M= 4.5e3 *mSun;
                n+= 3.e4*cmm; //M/V *1/m_p;
            }}
        
//NR Northern ridge #curved cylinder of radius 0.5 pc and length of 4pc#along norther body of SgrA East
        if (y5/b5>=0 and std::abs(x5/a5)<1. ){ //}#because we are interested in the western part otherwise we would obtain two times the same value which are axissymmetric
            theta=std::asin(x5/a5);
        
//phi=np.arcsin(x1/a1 *1/np.cos(theta))
            if (1.<=x5*x5 /(a5*a5) + y5*y5 /(b5*b5) + z5*z5 /(c5*c5) and x5*x5 /(a5*a5) + y5*y5 /(b5*b5) + z5*z5 /(c5*c5) <=1.5 and -0.5<=z5 and z5 <=0.5  and -pi*2./12.<=theta and theta <= pi *1.5/12.){ //first if condition: ellipsoid ring with 1pc diameter, second if condition: cylinder changing in direction of x and not in z direction and cylinder surface is also in z direction ->diameter of 1pc, third if condition: we need only a part of the ellipsoid ring which is also a "curved cylinder" ->length of 4pc
                V=pi*0.5*0.5 *4*pc*pc*pc;
                M= 2.2e3 *mSun;
                n+= 3.e4*cmm; //M/V *1/m_p;
            }}
        else{// intercloud medium in CMZ from Ferriere 2007
            double x = position.x;
            double y = position.y;
            double z = position.z;
            n+=2*128.4*std::exp(-pow(((std::sqrt(x*x+(2.5*y)*(2.5*y))-125.*pc)*1./(137.*pc)),4.)) *std::exp(-pow((z/(18.*pc)),2.))*cmm;//##molecular Density H2#actually 150 instead of 128.4# 128.4 substract the considered MC and the inner 10 pc part
            n+=7.5*std::exp(-pow((std::sqrt(x*x+(2.5*y)*(2.5*y))-125.*pc)*1./(137.*pc), 4.)) *std::exp(-pow((z/(54.*pc)),2.))*cmm;// Density H#7.2 insteadt of 8.8 (see above)
        }
        return n;
    }
    
    

    
} //END NAMESPACE CRPROPA

