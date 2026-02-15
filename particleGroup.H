#ifndef _particleGroup_H_
#define _particleGroup_H_

#include <vector>
#include <Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::MatrixXi;

class particleGroup{
    public:
        MatrixXd pos; //position
        MatrixXd posRef; //reference position recorded for refreshing the neighbor list
        MatrixXd vel; //velocity
        MatrixXd acc; //acceleration
        MatrixXd distTravelled;

        MatrixXd w; //angular velocity
        MatrixXd aw; //angular acceleration
       
        MatrixXi collisionHistory;
        MatrixXi collisionHistoryNew;
        MatrixXi numCollisionHistory;
        MatrixXi numCollisionHistoryWall;
        MatrixXi collisionFlag;
        MatrixXi collisionFlagWall;
        MatrixXd deltatHistory;
        MatrixXd deltatHistoryNew;

        MatrixXi collisionHistoryWall;
        MatrixXi collisionHistoryNewWall;
        MatrixXi numCollisionHistoryNew;
        MatrixXi numCollisionHistoryNewWall;
        MatrixXd deltatHistoryWall;
        MatrixXd deltatHistoryNewWall;

        MatrixXd tmpDist;

        particleGroup();
        particleGroup(double density, double radius, int numParticles);

        inline const MatrixXd &getDeltatHistoryWall() const{
            return deltatHistoryWall;
        }

        inline const MatrixXd &getDeltatHistory() const{
            return deltatHistory;
        }


        /*
        inline void printPos() const{
            std::cout << this->pos[0] << " " << this->pos[1] << " " << this->pos[2] << std::endl;
        }

        inline void printVel() const{
            std::cout << this->vel[0] << " " << this->vel[1] << " " << this->vel[2] << std::endl;
        }

        inline void printAcc() const{
            std::cout << this->acc[0] << " " << this->acc[1] << " " << this->acc[2] << std::endl;
        }
        */

        inline const MatrixXd  &getPos() const{
            return pos;
        }


        inline const MatrixXd  &getVel() const{
            return vel;
        }


        inline const MatrixXd  &getAcc() const{
            return acc;
        }

        inline const MatrixXi  &getNumCollisionHistory() const{
            return numCollisionHistory;
        }

        inline const MatrixXi  &getNumCollisionHistoryWall() const{
            return numCollisionHistoryWall;
        }

        inline const MatrixXd  &getW() const{
            return w;
        }


        inline const MatrixXd  &getAw() const{
            return aw;
        }



        inline const MatrixXd  &radius() const{
            return radius_;
        }

        inline const MatrixXd  &radiusSq() const{
            return radiusSq_;
        }
        inline const MatrixXd &density() const{
            return density_;
        }

        inline const MatrixXd &mass() const{
            return mass_;
        }
        
        inline const MatrixXd &massInv() const{
            return massInv_;
        }
        inline const MatrixXd &moi() const{
            return moi_;
        }

        inline const MatrixXd &moiInv() const{
            return moiInv_;
        }

        inline const MatrixXd &getTmpDist() const{
            return tmpDist;
        }


        inline const int &numPart() const{
            return numParticles_;
        }

        inline MatrixXi &writeNeighborList(){
            return neighborList_;
        }

        inline const MatrixXi &getNeighborList() const{
            return neighborList_;
        }

        inline const MatrixXd &getDistTravelled() const{
            return distTravelled;
        }


    private:
        MatrixXd density_,radius_,radiusSq_,volume_,mass_,massInv_,moi_,moiInv_; //moi = moment of inertia
        MatrixXi neighborList_;
        int numParticles_;
};

#endif
