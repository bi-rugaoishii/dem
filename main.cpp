#include <iostream>
#include <chrono>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>
#include <string>
#include "particle.H"
#include "readFiles.H"
#include "fileOutput.H"
#include "demCalc.H"
#include "triangles.H"
#include "boundingBox.H"
#include "vector3d.H"
#include "misc.H"
#define PI 3.14159265359
#include "particleGroup.H"
#include "timeAdvTest.H"
#include <Eigen/Dense>
#include <fstream>

int main(){
    auto start = std::chrono::steady_clock::now();
    readFiles readFile;

    readFile.readSetting();
    
    int seed = readFile.seed;
    double dt = readFile.dt;
    double startTime = readFile.startTime;
    double wallAmp = readFile.wallAmp;
    double endTime = readFile.endTime;
    double outputTiming = readFile.outputTiming;
    double density;
    double radius;
    double maxRad = 0.;


    double e = readFile.particleTypes["particle1"].get<picojson::object>()["CoR"].get<double>();
    //std::cout << "hello" << std::endl;
    double k = readFile.particleTypes["particle1"].get<picojson::object>()["k"].get<double>();
    double mu = readFile.particleTypes["particle1"].get<picojson::object>()["mu"].get<double>();

    int numPart=readFile.numPart;
    int numParticleTypes=readFile.particleTypes.size();
    std::vector<particle> particleTypes(numParticleTypes);


    int i = 0;
    for (picojson::object::iterator p=readFile.particleTypes.begin(); p!= readFile.particleTypes.end(); p++){
        density = p->second.get<picojson::object>()["density"].get<double>();
        radius = p->second.get<picojson::object>()["radius"].get<double>();
        particleTypes[i] = particle(density, radius);
        i++;
        /* find maximum Radius */
        if (maxRad < radius){
            maxRad = radius;
        }
    }

    std::cout << "maxRad is " <<  maxRad << std::endl;

    particleGroup particleGroup1(particleTypes[0].density(),particleTypes[0].radius(),numPart);

    std::cout << "radius = " << particleGroup1.radius()(0) << std::endl;
    std::cout << "density = " << particleGroup1.density()(0) << std::endl;
    std::cout << "mass = " << particleGroup1.mass()(0) << std::endl;
    std::cout << "moi = " << particleGroup1.moi()(0) << std::endl;

    //calculate eta from coeff of res
    double eta = -2.*log(e)*sqrt(particleGroup1.mass()(0)*k/(pow(PI,2.)+pow(log(e),2.0)));


    std::cout << "boxXmin boxXmax boxYmin boxYmax boxZmin boxZmax" << std::endl;
    std::cout << readFile.boxXmin << " " << readFile.boxXmax << " " << readFile.boxYmin << " "<< readFile.boxYmax 
        << " " << readFile.boxZmin << " " << readFile.boxZmax << std::endl;
    std::cout <<std::endl;


    boundingBox box(readFile.boxXmin,readFile.boxXmax,readFile.boxYmin,readFile.boxYmax,readFile.boxZmin,readFile.boxZmax);
    box.showRange();
    box.split(maxRad);
    box.showSplit();

    std::cout << std::endl;
    triangles walls;
    walls = readFile.readStl("geometry/box.stl");
    walls.getEdgeInfo();
    box.trianglesIntersect(walls);
    std::cout << "Created " << walls.numTriangles << " walls." << std::endl;
    std::cout << std::endl;
    std::cout << "Normals are " << walls.fNormal.transpose() <<  std::endl;
    std::cout << std::endl;
    
    std::cout << "eNormals are " << walls.eNormal.transpose() <<  std::endl;

    std::cout << "NumParticleTypes is " << readFile.numParticleTypes << std::endl;

    demCalc demCalc1(box,particleGroup1,walls);
    demCalc1.setRefreshFreq(readFile.refreshFreq);
    std::cout << "numPart = " << demCalc1.numPart() << std::endl;

    demCalc1.shuffle(seed,readFile.shuffleXmin,readFile.shuffleXmax,readFile.shuffleYmin,readFile.shuffleYmax,readFile.shuffleZmin,readFile.shuffleZmax);

    /*
    demCalc1.getWriteParticle().pos.col(0)= Eigen::Vector3d(0.0,0.0,0.0);
    demCalc1.getWriteParticle().pos.col(1)= Eigen::Vector3d(0.1,0.0,0.1);
    demCalc1.getWriteParticle().pos.col(2)= Eigen::Vector3d(-0.1,0.0,-0.1);
    */

    /*
    for (int i = 0; i<demCalc1.numPart(); i++){
        std::cout << demCalc1.getParticle().getPos()(0,i) << " " << demCalc1.getParticle().getPos()(1,i) << " " << demCalc1.getParticle().getPos()(2,i) << std::endl;
    }
    */

    std::cout << "numPart seed dt startTime endTime outputTiming" << std::endl;
    std::cout << readFile.numPart << " " << readFile.seed << " " << readFile.dt << " " 
        << readFile.startTime << " " << readFile.endTime << " " << readFile.outputTiming << std::endl;


    //demCalc1.getParticleInCellList();
    // demCalc1.showNumPartInCell();
    demCalc1.initLinkedList();
    std::cout << "initial refreshLinkedList" << std::endl;
    demCalc1.refreshLinkedList();
    std::cout << "refreshLinkedListDone" << std::endl;
    //demCalc1.printCellIniPart();
    //demCalc1.printLinkedList();
    //demCalc1.getDist();




    demCalc1.gravity3D(0.,-9.81,0.); //sets gravity

    demCalc1.setParameters(k,eta,mu);

    fileOutput output;
    output.createResultDirectory();

    timeAdvTest timeAdv(dt,startTime,endTime,outputTiming,wallAmp,start);

    // time advancement// 

    timeAdv.gogogo(demCalc1);

    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    std::cout << "Calculation done "<< std::endl;
    std::cout << "Execution time = ";
    std::cout << std::chrono::duration<double, std::milli>(diff).count()*1e-3 << " s"<< std::endl;

    return 0;
}
