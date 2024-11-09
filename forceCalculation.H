#ifndef _forceCalculation_H_
#define _forceCalculation_H_

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>
#include <cmath>
#include <string>
#define PI 3.14159265359
#include "particle.H"
#include "demCalc.H"
#include "planeBoundary.H"
#include "vector3d.H"
#include <Eigen/Dense>
#include "misc.H"

class forceCalculation{
    public:
        forceCalculation();
        inline  forceCalculation(double k, double eta, double mu){
            eta_ = eta;
            k_ = k;
            mu_ = mu;
        }; 

        inline void gravity3D(double x,double y, double z){

            gravSource_(0) = x;
            gravSource_(1) = y;
            gravSource_(2) = z;
        }

        /*
           inline double getDist(const particle &particle1, const particle &particle2){
           double dist;
           dist = distVec(particle1.getPos(),particle2.getPos());
           return dist;
           }
           */


        inline double getDist(const particleGroup &particle, int indPart, const planeBoundary &wall, int wallInd){
            double dist;
            dist = fabs(particle.getPos().col(indPart).dot(wall.getNormal().col(wallInd)) + wall.d()(wallInd));

            return dist;
        }



        /*
        //gets vector from particle2 to particle1
        inline vector3d getNorm(const particle &particle1, const particle &particle2, double dist){
        vector3d normal;
        normal = 1./dist*(particle1.getPos()-particle2.getPos());
        return normal;
        }
        */

        Eigen::Vector3d getNorm(const particleGroup &particle,int indPart, const planeBoundary &wall, int indWall);

        /*
           inline vector3d getDelta(const vector3d &normal, const particle &particle1, const particle &particle2,const double &dist){
           vector3d delta;
           double delta_mag;

           delta_mag = dist - (particle1.radius()+particle2.radius());
           delta = delta_mag * normal;
           return delta;
           }
           */

        inline Eigen::Vector3d getDelta(const Eigen::Vector3d &normal, const particleGroup &particle,int indPart, const planeBoundary &wall,int indWall,const double &dist){
            Eigen::Vector3d delta;
            double delta_mag;

            delta_mag = dist - (particle.radius()(indPart));
            delta = delta_mag * normal;

            return delta;
        }


        /*
           inline  vector3d getvn(const vector3d &normal, const particle &particle1, const particle &particle2){
           vector3d vn;
           double vn_mag;

           vn_mag = (particle1.getVel() - particle2.getVel()) & normal;
           vn = vn_mag*normal;

           return vn;
           }
           */

        inline Eigen::Vector3d getvn(const Eigen::Vector3d &normal, const particleGroup &particle1,int indPart, const planeBoundary &wall, int indWall){
            Eigen::Vector3d vn;
            double vn_mag;

            vn_mag = (particle1.getVel().col(indPart)-wall.getVel().col(indWall)).dot(normal);
            vn = vn_mag * normal;

            return vn;
        }

        /*
           inline  vector3d getvt(const vector3d &vn, const vector3d &normal, const particle &particle1, const particle &particle2){
           vector3d vt;
           vector3d rot;
           vector3d vr; //relative velocity

           vr = particle1.getVel() - particle2.getVel();

           rot = particle1.radius()*particle1.getW() + particle2.radius()*particle2.getW();
           rot = rot*normal ; //cross product


           vt = vr - vn + rot;
           return vt;
           }

           inline  vector3d getvt(const vector3d &vn,const vector3d &normal, const particle &particle1, const planeBoundary &wall){
           vector3d vt;
           vector3d rot;;

           vector3d vr; //relative velocity

           vr = particle1.getVel() - wall.getVel();

           rot = particle1.radius()*particle1.getW();
           rot = rot*normal ; //cross product


           vt = vr - vn + rot;

           return vt;
           }
           */

        inline Eigen::Vector3d getNormalForce(const Eigen::Vector3d &delta_n,const  Eigen::Vector3d &vn){
            Eigen::Vector3d resVec = delta_n; //initialize

            resVec = -k_*delta_n-eta_*vn;
            return resVec;
        }

        inline vector3d getTangentialForce(const vector3d &delta_t, const vector3d &vt){
            vector3d resVec = delta_t; //initialize

            resVec = -k_*delta_t - eta_*vt;
            return resVec;
        }

        inline vector3d getMoment(const vector3d &normal, const particle &particle1, const vector3d &ft){
            vector3d resMom;
            resMom = particle1.radius()*(normal*ft); //normal*ft is cross product

            return resMom;
        }


        void getAcc(demCalc &demCalc1, const double &dt);


    private:
        Eigen::Vector3d gravSource_;
        double k_;
        double eta_;
        double mu_;
};

#endif
