#ifndef _timeAdvTest_H_
#define _timeAdvTest_H_
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>
#include <cmath>
#include <string>
#include <chrono>
#define PI 3.14159265359
#include "particle.H"
#include "triangles.H"
#include "fileOutput.H"
#include "demCalc.H"
#include "vector3d.H"
#include "boundingBox.H"
#include "misc.H"


class timeAdvTest{
    public:
        timeAdvTest();

       // timeAdvTest(forceCalculation &forceCalc,demCalc &demCalc1, double dt, double startTime, double endTime,double outputTiming,double wallAmp, fileOutput &output);
        timeAdvTest(double dt, double startTime, double endTime,double outputTiming,double wallAmp, std::chrono::steady_clock::time_point startCalcTime);


        inline bool isOutOfBound(const particle  &particle1, const boundingBox &box){
            if (particle1.getPos()[0]<box.minX()){
                return true;
            }else if(particle1.getPos()[0]>box.maxX()){
                return true;
            }else if(particle1.getPos()[1]<box.minY()){
                return true;
            }else if(particle1.getPos()[1]>box.maxY()){
                return true;
            }else if(particle1.getPos()[2]>box.minZ()){
                return true;
            }else if(particle1.getPos()[2]>box.maxZ()){
                return true;
            }else{
                return false;
            }
        }

        void nextStep(demCalc& demCalc1);

        void gogogo(demCalc& demCalc1);

    private:
        fileOutput output_;

        int outputTimeStepNum_,timeStepNum_;
        double dt_,startTime_,endTime_,currentTime_,dtsq_,outputTiming_,wallAmp_;
       std::chrono::steady_clock::time_point startCalcTime_;

};

#endif
