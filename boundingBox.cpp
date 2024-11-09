#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>
#include <cmath>
#include <string>
#define PI 3.14159265359
#include "boundingBox.H"
#include "misc.H"


boundingBox::boundingBox(){};
boundingBox::boundingBox(double minX, double maxX, double minY, double maxY,double minZ, double maxZ){
    if (minX>=maxX || minY>=maxY || minZ>=maxZ){
        std::cout << "min is greater than max!\n" << std::endl;
        abort();
    }else{
        minX_=minX;
        maxX_=maxX;
        minY_=minY;
        maxY_=maxY;
        minZ_=minZ;
        maxZ_=maxZ;
        rangeX_=maxX_-minX_;
        rangeY_=maxY_-minY_;
        rangeZ_=maxZ_-minZ_;
        isSplit_ = false;
    }
} 

void boundingBox::showRange(){
    std::cout << "Range is x = [" << minX_ <<","<< maxX_ <<"] y = ["<<minY_ <<","<< maxY_ <<"]"<<std::endl;
    std::cout << "RangeSize is x = " << rangeX_<< " y = "<< rangeY_ <<std::endl;
}

void boundingBox::showSplit(){
    if (isSplit_ == false){
        std::cout << "not split yet!" << std::endl;
    }else{
        std::cout << "sizeSplitX = " << sizeSplitX_ <<" "<<"sizeSplitY = "<< sizeSplitY_ << " sizeSplitZ = " << sizeSplitZ_ << std::endl;

    }
}

void boundingBox::split(double maxRad){
    double maxDiam = maxRad*2.;
    if (maxDiam>=rangeX_ ||maxDiam>=rangeY_ || maxDiam>=rangeZ_){
        std::cout <<"radius is too big!" <<std::endl;
    }else{
        dx_ = maxDiam*2.;
        dy_ = maxDiam*2.;
        dz_ = maxDiam*2.;

        sizeSplitX_=ceil(rangeX_/dx_) + 2;
        sizeSplitY_=ceil(rangeY_/dy_) + 2;
        sizeSplitZ_=ceil(rangeZ_/dz_) + 2;
        indBox_.resize(sizeSplitX_*sizeSplitY_*sizeSplitZ_);
        isSplit_ = true;
    }
}





