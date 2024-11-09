#include <iostream>
#include <fstream>
#include "fileOutput.H"

fileOutput::fileOutput(){};
void fileOutput::createResultDirectory(){
    std::string dir="results";
    if (mkdir(dir.c_str(),0777)==0){
        std::cout << "results folder created!!" <<std::endl<<std::endl;
    }else{
        std::cout << "results folder already exists" <<std::endl<<std::endl;
    }
}
void fileOutput::writeHeaders(std::ofstream &output){
    output << "#time(s) x y z vx vy vz ax ay az wx wy wz" << std::endl;
}
void fileOutput::writeHeadersWalls(std::ofstream &output){
    output << "#time(s) x y z normX normY normZ" << std::endl;
}

void fileOutput::writeHeadersRestart(std::ofstream &output){
    output << "#time(s) timestepNum" << std::endl;
}

void fileOutput::writeVtk(const double currentTime,const particleGroup &particle, std::ofstream &output){
    output << "# vtk DataFile Version 3.0" << std::endl; 
    output << "results.vtk" << std::endl; 
    output << "ASCII" << std::endl; 
    output << "DATASET UNSTRUCTURED_GRID" << std::endl; 
    output <<  std::endl; 
    output << "POINTS " << particle.numPart()<< " float" << std::endl; 
    for(int pInd=0; pInd<particle.numPart(); pInd++){
        output << particle.getPos()(0,pInd)<<" "<< particle.getPos()(1,pInd) << " " <<particle.getPos()(2,pInd) << std::endl;
    }
    output <<  std::endl; 

    output << "CELL_TYPES " << particle.numPart() << std::endl; 
    for(int pInd=0; pInd<particle.numPart(); pInd++){
        output << "1" << std::endl;
    }
    output <<  std::endl; 

    output << "POINT_DATA " << particle.numPart() << std::endl; 
    output << "SCALARS radius float" << std::endl; 
    output << "LOOKUP_TABLE default" << std::endl; 
    for(int pInd=0; pInd<particle.numPart(); pInd++){
        output << particle.radius()(pInd)<< std::endl;
    }
}

void fileOutput::writeResults(const double currentTime,const particleGroup &particle, std::ofstream &output){
    for(int pInd=0; pInd<particle.numPart(); pInd++){
        output << currentTime << " " << particle.getPos()(0,pInd)<<" "<< particle.getPos()(1,pInd) << " " <<particle.getPos()(2,pInd)
            << " "  << particle.getVel()(0,pInd) << " " << particle.getVel()(1,pInd) << " "<< particle.getVel()(2,pInd) << " "
            <<particle.getAcc()(0,pInd) << " "<<particle.getAcc()(1,pInd) << " " << particle.getAcc()(2,pInd) <<" "
            <<particle.getW()(0,pInd) << " "<<particle.getW()(1,pInd) << " " << particle.getW()(2,pInd) 
            << std::endl; 

    }
}


void fileOutput::writeResultsRestart(const double currentTime,const int timeStepNum, std::ofstream &output){
    output << currentTime << " " << timeStepNum<<std::endl; 
}



