#ifndef _particle_H_
#define _particle_H_

#include <vector>
#include "vector3d.H"

class particle{
    public:
        vector3d pos; //position
        vector3d vel; //velocity
        vector3d acc; //acceleration

        vector3d w; //angular velocity
        vector3d aw; //angular acceleration
        std::vector<int> collisionHistory;
        std::vector<vector3d> deltatHistory;

        std::vector<int> collisionHistoryWall;
        std::vector<vector3d> deltatHistoryWall;


        particle();
        particle(double density, double radius);

        inline const std::vector<vector3d> &getDeltatHistoryWall() const{
            return deltatHistoryWall;
        }

        inline const std::vector<vector3d> &getDeltatHistory() const{
            return deltatHistory;
        }

        inline const int collisionHistorySize() const{
            return collisionHistory.size();
        }

        inline const int collisionHistoryWallSize() const{
            return collisionHistoryWall.size();
        }

        inline void printPos() const{
            std::cout << this->pos[0] << " " << this->pos[1] << " " << this->pos[2] << std::endl;
        }

        inline void printVel() const{
            std::cout << this->vel[0] << " " << this->vel[1] << " " << this->vel[2] << std::endl;
        }

        inline void printAcc() const{
            std::cout << this->acc[0] << " " << this->acc[1] << " " << this->acc[2] << std::endl;
        }

        inline const vector3d  &getPos() const{
            return pos;
        }


        inline const vector3d  &getVel() const{
            return vel;
        }


        inline const vector3d  &getAcc() const{
            return acc;
        }


        inline const vector3d  &getW() const{
            return w;
        }


        inline const vector3d  &getAw() const{
            return aw;
        }



        inline const double  &radius() const{
            return radius_;
        }

        inline const double &density() const{
            return density_;
        }

        inline const double &mass() const{
            return mass_;
        }
        
        inline const double &moi() const{
            return moi_;
        }

        inline const double &moiInv() const{
            return moiInv_;
        }



    private:
        double density_,radius_,volume_,mass_,moi_,moiInv_; //moi = moment of inertia
};

#endif
