#ifndef _fileOutput_H_
#define _fileOutput_H_

#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include "particle.H"
#include "triangles.H"
#include "particleGroup.H"

class fileOutput{
    public:
        fileOutput();
        void writeHeaders(std::ofstream &output);
        void createResultDirectory();
        void writeHeadersRestart(std::ofstream &output);
        void writeHeadersWalls(std::ofstream &output);
        void writeResults(const double currentTime,const particleGroup &particle, std::ofstream &output);
        void writeResults(const double currentTime,const triangles &wall, std::ofstream &output);
        void writeVtk(const double currentTime,const particleGroup &particle, std::ofstream &output);
        void writeResultsRestart(const double currentTime,const int timeStepNum, std::ofstream &output);

};

#endif
