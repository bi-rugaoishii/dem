#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>
#include <cmath>
#include <string>
#define PI 3.14159265359
#include "particle.H"
#include "misc.H"



particle::particle(){ //default contructor
    density_ = 1.;
    radius_ = .001;
    volume_ = pow(radius_,3)*PI*4./3.;
    mass_ = volume_*density_;
    moi_ = 2./5.*mass_*pow(radius_,2.); //for sphere
    moiInv_ = 1./moi_;
}

particle::particle(double density, double radius){
    density_ = density;
    radius_ = radius;
    volume_ = pow(radius_,3)*PI*4./3.;
    mass_ = volume_*density_;
    moi_ = 2./5.*mass_*pow(radius_,2.);
    moiInv_ = 1./moi_;
}

