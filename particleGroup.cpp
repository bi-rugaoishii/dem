#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>
#include <cmath>
#include <string>
#define PI 3.14159265359
#define MAX_HISTORY_NEIGH 50
#include "particleGroup.H"



particleGroup::particleGroup(){ //default contructor
}

particleGroup::particleGroup(double density, double radius,int numParticles){
    int dim=3;
    density_.resize(numParticles,1) ;
    density_.fill(density);
    radius_.resize(numParticles,1) ;
    radius_.fill(radius);
    radiusSq_.resize(numParticles,1) ;
    radiusSq_.fill(radius*radius);
    volume_.resize(numParticles,1) ;
    volume_.fill(pow(radius,3)*PI*4./3.);
//    mass_.resize(numParticles,1);
    mass_ = volume_.array()*density_.array();
 //   massInv_.resize(numParticles,1);
    massInv_ = mass_.cwiseInverse();
    moi_ = 2./5.*mass_*pow(radius,2.);
    moiInv_ = moi_.cwiseInverse();
    pos.resize(dim,numParticles);
    pos.setZero();
    vel.resize(dim,numParticles);
    vel.setZero();
    tmpDist.resize(dim,numParticles);
    tmpDist.setZero();
    acc.resize(dim,numParticles);
    acc.setZero();
    w.resize(dim,numParticles);
    w.setZero();
    aw.resize(dim,numParticles);
    aw.setZero();

    posRef.resize(dim,numParticles);
    posRef.setZero();

    collisionHistory.resize(numParticles,MAX_HISTORY_NEIGH);
    collisionHistory.setZero();

    collisionHistoryNew.resize(numParticles,MAX_HISTORY_NEIGH);
    collisionHistoryNew.setZero();


    deltatHistory.resize(dim,numParticles*MAX_HISTORY_NEIGH);
    deltatHistory.setZero();

    deltatHistoryNew.resize(dim,numParticles*MAX_HISTORY_NEIGH);
    deltatHistoryNew.setZero();

    collisionHistoryWall.resize(numParticles,MAX_HISTORY_NEIGH);
    collisionHistoryWall.setZero();
    
    collisionHistoryNewWall.resize(numParticles,MAX_HISTORY_NEIGH);
    collisionHistoryNewWall.setZero();

    deltatHistoryWall.resize(dim,numParticles*MAX_HISTORY_NEIGH);
    deltatHistoryWall.setZero();

    deltatHistoryNewWall.resize(dim,numParticles*MAX_HISTORY_NEIGH);
    deltatHistoryNewWall.setZero();

    distTravelled.resize(numParticles,1);
    distTravelled.setZero();
    neighborList_.resize(numParticles,MAX_HISTORY_NEIGH);
    neighborList_.setConstant(-1);
    numParticles_=numParticles;

    numCollisionHistory.resize(numParticles,1);
    numCollisionHistory.setZero();

    numCollisionHistoryWall.resize(numParticles,1);
    numCollisionHistoryWall.setZero();

    numCollisionHistoryNew.resize(numParticles,1);
    numCollisionHistoryNew.setZero();

    numCollisionHistoryNewWall.resize(numParticles,1);
    numCollisionHistoryNewWall.setZero();

}

