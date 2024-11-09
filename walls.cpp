#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>
#include <cmath>
#include <string>
#define PI 3.14159265359
#include "misc.H"

class walls{
    public:
        walls();
        walls(double startPointX, double startPointY, double endPointX, double endPointY);
        inline std::vector<double> &normal(){
            return normal_;
        }

        inline const std::vector<double> &startPoint() const{
            return startPoint_;
        }
        
        inline const std::vector<double> &endPoint() const{
            return endPoint_;
        }

        inline std::vector<double> &writeStartPoint(){
            return startPoint_;
        }
        
        inline std::vector<double> &writeEndPoint(){
            return endPoint_;
        }

        inline void recalcC(){
            c_ = -(normal_[0]*startPoint_[0] + normal_[1]*startPoint_[1]);
        }

    private:
        double c_;
        std::vector<double> startPoint_;
        std::vector<double> endPoint_;
        std::vector<double> normal_;
};


walls::walls(){}

walls::walls(double startPointX, double startPointY, double endPointX, double endPointY){
    startPoint_=std::vector<double>{startPointX,startPointY};
    endPoint_=std::vector<double>{endPointX,endPointY};
    std::vector<double> tmpVec(2,0.);
    normal_ = std::vector<double>(2,0.);

    tmpVec[0] = endPoint_ [0] - startPoint_[0];
    tmpVec[1] = endPoint_ [1] - startPoint_[1];

    if (tmpVec[0] != 0){
        normal_[1] = pow(tmpVec[1]/tmpVec[0],2.)+1.;
        normal_[1] = sqrt(1./(normal_[1]));
        normal_[0] = -tmpVec[1]*normal_[1]/tmpVec[0];

    }else if (tmpVec[1]!=0){
        normal_[0] = pow(tmpVec[0]/tmpVec[1],2.)+1.;
        normal_[0] = sqrt(1./(normal_[0]));
        normal_[1] = -tmpVec[0]*normal_[0]/tmpVec[1];

    }

    c_ = -(normal_[0]*startPoint_[0] + normal_[1]*startPoint_[1]);
};
