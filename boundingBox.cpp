#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <random>
#include <cmath>
#include <string>
#include <boost/format.hpp>
#define PI 3.14159265359
#define NUM_MAX_TRIANGLE 100 //max number of possible triangle in a grid
#define NUM_VERTEX 3
#include "boundingBox.H"
#include <Eigen/Dense>
#include "misc.H"
#include "triangles.H"


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

   this->num_max_triangle=NUM_MAX_TRIANGLE;
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
        dx_ = maxDiam*4.;
        dy_ = maxDiam*4.;
        dz_ = maxDiam*4.;

        sizeSplitX_=ceil(rangeX_/dx_) + 2;
        sizeSplitY_=ceil(rangeY_/dy_) + 2;
        sizeSplitZ_=ceil(rangeZ_/dz_) + 2;
        indBox_.resize(sizeSplitX_*sizeSplitY_*sizeSplitZ_);
        triangleList.resize(sizeSplitX_*sizeSplitY_*sizeSplitZ_*(NUM_MAX_TRIANGLE),-1); 
        isSplit_ = true;
    }
}

void boundingBox::trianglesIntersect(triangles &walls){
    //iterate over the grid
    double halfdx=dx_/2.;
    //axis vector
    Eigen::Vector3d x0(1.0,0.0,0.0);
    Eigen::Vector3d y0(0.0,1.0,0.0);
    Eigen::Vector3d z0(0.0,0.0,1.0);
    Eigen::MatrixXd listAxis(3,13);

    for (int triangleInd=0; triangleInd<walls.numTriangles; triangleInd++){
        //+1 for the ghost cell
     //   std::cout << boost::format("%f %f %f %f %f %f")%walls.minX[triangleInd]%walls.minY[triangleInd]%walls.minZ[triangleInd]%walls.maxX[triangleInd]%walls.maxY[triangleInd]%walls.maxZ[triangleInd] <<std::endl;
        int minX=static_cast<int>((walls.minX[triangleInd]-this->minX())/this->dx())+1;
        int minY=static_cast<int>((walls.minY[triangleInd]-this->minY())/this->dx())+1;
        int minZ=static_cast<int>((walls.minZ[triangleInd]-this->minZ())/this->dx())+1;
        int maxX=static_cast<int>((walls.maxX[triangleInd]-this->minX())/this->dx())+1;
        int maxY=static_cast<int>((walls.maxY[triangleInd]-this->minY())/this->dx())+1;
        int maxZ=static_cast<int>((walls.maxZ[triangleInd]-this->minZ())/this->dx())+1;
      //  std::cout << boost::format("%d %d %d %d %d %d")%minX%minY%minZ%maxX%maxY%maxZ <<std::endl;

        for(int z=minZ; z<=maxZ; z++){
            for(int y=minY; y<=maxY; y++){
                for(int x=minX; x<=maxX; x++){
                    bool hasCollided = true;
                    //check if triangle intersects with the grid
                    //shifted by 1 as minXs doesn't include the ghost cell
                    Eigen::Vector3d cellCenterPos(dx_*(x-1)+halfdx+this->minX(),dy_*(y-1)+halfdx+this->minY(),dz_*(z-1)+halfdx+this->minZ());

                    //std::cout <<"centerPos " << cellCenterPos.transpose() <<std::endl;
                    //shift vertices to cell center coordinate
                    Eigen::Vector3d v0 = walls.f.col(triangleInd*NUM_VERTEX+0)-cellCenterPos; 
                    Eigen::Vector3d v1 = walls.f.col(triangleInd*NUM_VERTEX+1)-cellCenterPos; 
                    Eigen::Vector3d v2 = walls.f.col(triangleInd*NUM_VERTEX+2)-cellCenterPos; 
                    //std::cout <<"v0 before shift " << walls.f.col(triangleInd*NUM_VERTEX+0).transpose() <<std::endl;
                    //std::cout <<"v0 " << v0.transpose() <<std::endl;
                    //std::cout << std::endl;

                    const Eigen::Vector3d e0 = walls.e.col(triangleInd*NUM_VERTEX+0); 
                    const Eigen::Vector3d e1 = walls.e.col(triangleInd*NUM_VERTEX+1); 
                    const Eigen::Vector3d e2 = walls.e.col(triangleInd*NUM_VERTEX+2); 

                    //create 13 axis to be checked
                    listAxis.col(0)=x0.cross(e0); 
                    listAxis.col(0).normalize();
                    listAxis.col(1)=y0.cross(e0); 
                    listAxis.col(1).normalize();
                    listAxis.col(2)=z0.cross(e0); 
                    listAxis.col(2).normalize();
                    listAxis.col(3)=x0.cross(e1); 
                    listAxis.col(3).normalize();
                    listAxis.col(4)=y0.cross(e1); 
                    listAxis.col(4).normalize();
                    listAxis.col(5)=z0.cross(e1); 
                    listAxis.col(5).normalize();
                    listAxis.col(6)=x0.cross(e2); 
                    listAxis.col(6).normalize();
                    listAxis.col(7)=y0.cross(e2); 
                    listAxis.col(7).normalize();
                    listAxis.col(8)=z0.cross(e2); 
                    listAxis.col(8).normalize();
                    listAxis.col(9)=x0; 
                    listAxis.col(9).normalize();
                    listAxis.col(10)=y0; 
                    listAxis.col(10).normalize();
                    listAxis.col(11)=z0; 
                    listAxis.col(11).normalize();
                    listAxis.col(12)=walls.fNormal.col(triangleInd); 
                    listAxis.col(12).normalize();

                    //seperating axis test
                    for (int j=0; j<13; j++){ //13 is the number of axis to be checked
                        double p0=v0.dot(listAxis.col(j));
                        double p1=v1.dot(listAxis.col(j));
                        double p2=v2.dot(listAxis.col(j));

                        double r = halfdx*(std::abs(x0.dot(listAxis.col(j)))+std::abs(y0.dot(listAxis.col(j)))+std::abs(z0.dot(listAxis.col(j))))+dx_;//add dx at the end to make the condition little bit bigger
                        //std::cout <<"r = " << r << std::endl;
                        if (std::max({p0,p1,p2})<-r || std::min({p0,p1,p2})>r){
                            hasCollided = false;
                            break;
                        }
                    }
                    if (hasCollided){
                        //std::cout <<"collided!"<<std::endl;
                        int k=0;
                        int indCell=index(z,y,x,sizeSplitX_,sizeSplitY_);
                        while(triangleList[indCell*this->num_max_triangle+k]!=-1){
                            k+=1;
                        }
                        if (k<this->num_max_triangle){
                            triangleList[indCell*this->num_max_triangle+k]=triangleInd;
                        }else{
                            std::cerr <<"TRIANGLE OVER FLOW IN CELL!!"<< indCell << std::endl;
                            std::abort();
                        }
                    }
                }
            }
        }
    }
}




