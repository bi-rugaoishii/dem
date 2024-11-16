#include "timeAdvTest.H"

timeAdvTest::timeAdvTest(){}

/*
timeAdvTest::timeAdvTest(forceCalculation &forceCalc,demCalc &demCalc1, double dt, double startTime, double endTime,double outputTiming,double wallAmp, fileOutput &output){
    forceCalc_ = forceCalc;
    demCalc_ = demCalc1;
    dt_ = dt;
    dtsq_ = pow(dt,2.);
    startTime_ = startTime;
    endTime_ = endTime;
    currentTime_ = 0.;
    output_ = output;
    timeStepNum_ = 0;
    outputTiming_ = outputTiming;
    wallAmp_ = wallAmp;

    outputTimeStepNum_ = static_cast<int>(outputTiming_/dt_+1e-10); //1e-10 to prevent round off error
    std::cout << "outputTimeStepNum = " << outputTimeStepNum_ << std::endl;
}
*/

timeAdvTest::timeAdvTest(double dt, double startTime, double endTime,double outputTiming, double wallAmp, std::chrono::steady_clock::time_point startCalcTime){
    dt_ = dt;
    dtsq_ = pow(dt,2.);
    startTime_ = startTime;
    endTime_ = endTime;
    currentTime_ = 0.;
    timeStepNum_ = 0;
    outputTiming_ = outputTiming;
    wallAmp_ = wallAmp;
    startCalcTime_ = startCalcTime;

    outputTimeStepNum_ = static_cast<int>(outputTiming_/dt_+1e-10); //1e-10 to prevent round off error
    std::cout << "outputTimeStepNum = " << outputTimeStepNum_ << std::endl;
}


void timeAdvTest::nextStep(demCalc& demCalc1){
    int i;
    int numPart=demCalc1.numPart();
    std::string filename;
    std::ofstream fp;
    int checkMod = 0;

    demCalc1.getAcc(dt_);

    //position integration
    demCalc1.getWriteParticle().vel = demCalc1.getParticle().vel + dt_*demCalc1.getParticle().acc;
    demCalc1.getWriteParticle().pos = demCalc1.getParticle().pos+demCalc1.getParticle().vel*dt_;
    for (int k=0; k<numPart; k++){
        demCalc1.getWriteParticle().distTravelled(k)=(demCalc1.getParticle().pos.col(k) - demCalc1.getParticle().posRef.col(k)).norm();
    }
    demCalc1.getWriteParticle().w = demCalc1.getWriteParticle().w + dt_*demCalc1.getParticle().aw;
    currentTime_ += dt_;
    timeStepNum_ += 1;

    //demCalc1.getWriteWall().move(currentTime_,wallAmp_,dt_);

    checkMod = timeStepNum_%outputTimeStepNum_;
    if (checkMod == 0){

        std::cout << "Writing result at t = " << currentTime_ << std::endl;

        filename = "results/"+ std::to_string(currentTime_) + ".txt";
        fp.open(filename);
        output_.writeHeaders(fp);
        output_.writeResults(currentTime_,demCalc1.getParticle(),fp);
        fp.close();

        /*
        filename = "results/"+ std::to_string(currentTime_) + "Walls.txt";
        fp.open(filename);
        output_.writeHeadersWalls(fp);
        output_.writeResults(currentTime_,demCalc1.getWall(),fp);
        fp.close();
        */

        /*
           filename = "results/restart.txt";
           fp.open(filename);
           output_.writeHeadersRestart(fp);
           output_.writeResultsRestart(currentTime_,timeStepNum_,fp);
           fp.close();
           */

        filename = "results/results"+ std::to_string(timeStepNum_) + ".vtk";
        fp.open(filename);
        output_.writeVtk(currentTime_,demCalc1.getParticle(),fp);
        fp.close();

        std::cout << "Result written." <<std::endl;

    }

    double maxDistTravelled = demCalc1.getParticle().getDistTravelled().maxCoeff();
    if (demCalc1.getRefreshThreshold()<maxDistTravelled){
        std::cout << "refreshing the linkedList" << std::endl;
        demCalc1.refreshLinkedList();
    }

    checkMod = (timeStepNum_)%1000;
    if (checkMod == 0){
        int checkParticleInd = 0;
        if (numPart<checkParticleInd){
            std::cout << "Check particle index larger than number of particles!!!!!!!!" << std::endl;
        }else{
            std::cout << currentTime_ << " " << endTime_ <<" "<< demCalc1.getParticle().getPos()(0,checkParticleInd) << " " << demCalc1.getParticle().getPos()(1,checkParticleInd) << " " << demCalc1.getParticle().getPos()(2,checkParticleInd) << std::endl; 
            auto end = std::chrono::steady_clock::now();
            auto diff = end - startCalcTime_;
            std::cout << "Execution time = ";
            std::cout << std::chrono::duration<double, std::milli>(diff).count()*1e-3 << " s"<< std::endl;
        }
    }
}

void timeAdvTest::gogogo(demCalc& demCalc1){
    while (currentTime_ < endTime_ ){
        this->nextStep(demCalc1);
    }
}


