#ifndef _readFiles_H_
#define _readFiles_H_
#include "misc.H"
#define PI 3.14159265359
#include "triangles.H"
#include <fstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <limits>
#include <boost/format.hpp>
#include "picojson.h"

class readFiles{
    public:
        readFiles();
        void readSetting();

        triangles readStl(std::string filename); 
        std::vector<double> normal;
        std::vector<double> particles;
        double density, radius,cor,k,mu;
        double dtFactor,startTime,endTime,outputTiming,wallAmp,refreshFreq;
        double boxXmin,boxXmax,boxYmin,boxYmax,boxZmin,boxZmax;
        double shuffleXmin,shuffleXmax,shuffleYmin,shuffleYmax,shuffleZmin,shuffleZmax;
        picojson::object particleTypes;

        int numPart,seed, numWalls, numParticleTypes;
};

#endif
