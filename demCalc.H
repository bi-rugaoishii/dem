#ifndef _demCalc_H_
#define _demCalc_H_

#include <vector>
#include <random>
#include <iostream>
#define PI 3.14159265359
#include "misc.H"
#include "boundingBox.H"
#include "triangles.H"
#include "particle.H"
#include "particleGroup.H"
#include <Eigen/Dense>


class demCalc{
    public:
        demCalc();
        demCalc(boundingBox &box,particleGroup &particles);

        demCalc(boundingBox box,particleGroup particles, triangles boundaries);

        /*
        inline double getDist(const particleGroup &particle, int indPart, const triangles &wall, int wallInd){
            double dist;
            const Eigen::Vector3d &particle1 = particle.getPos().col(indPart);
            dist = fabs(particle1.dot(wall.getNormal().col(wallInd)) + wall.d()(wallInd));

            return dist;
        }
        */
        

        inline double getDistSq(const particleGroup &particle, int indPart1,int indPart2){
            double distSq;
            const Eigen::Vector3d &particle1 = particle.getPos().col(indPart1);
            const Eigen::Vector3d &particle2 = particle.getPos().col(indPart2);
            distSq = (particle1-particle2).squaredNorm();

            return distSq;
        }

        Eigen::Vector3d getNorm(const particleGroup &particle,int indPart, const triangles &wall, int indWall);

        inline Eigen::Vector3d getNorm(const particleGroup &particle,int indPart1,int indPart2,double dist){
            Eigen::Vector3d normal(0.,0.,0.);
            const Eigen::Vector3d &particle1 = particle.getPos().col(indPart1);
            const Eigen::Vector3d &particle2 = particle.getPos().col(indPart2);

            normal = (particle1-particle2)/dist;
            return normal;
        }

        inline Eigen::Vector3d getDelta(const Eigen::Vector3d &normal, const particleGroup &particle,int indPart1, int indPart2 ,double dist){
            Eigen::Vector3d delta;
            double delta_mag;

            delta_mag = dist - (particle.radius()(indPart1)+particle.radius()(indPart2));
            delta = delta_mag * normal;

            return delta;
        }

        inline Eigen::Vector3d getDelta(const Eigen::Vector3d &normal, const particleGroup &particle,int indPart, const triangles &wall,int indWall,const double &dist){
            Eigen::Vector3d delta;
            double delta_mag;

            delta_mag = dist - (particle.radius()(indPart));
            delta = delta_mag * normal;

            return delta;
        }

        inline Eigen::Vector3d getNormalForce(const Eigen::Vector3d &delta_n,const  Eigen::Vector3d &vn){
            Eigen::Vector3d resVec = delta_n; //initialize

            resVec = -k_*delta_n-eta_*vn;
            return resVec;
        }

        inline Eigen::Vector3d getTangentialForce(const Eigen::Vector3d &delta_t, const Eigen::Vector3d &vt){
            Eigen::Vector3d resVec = delta_t; //initialize

            resVec = -k_*delta_t - eta_*vt;
            return resVec;
        }

        inline Eigen::Vector3d getMoment(const Eigen::Vector3d &normal, const particleGroup &particle, int p1Ind, const Eigen::Vector3d &ft){
            Eigen::Vector3d resMom;
            resMom = particle.radius()(p1Ind)*(normal.cross(ft)); //normal*ft is cross product

            return resMom;
        }


        void getAcc(double dt);

        inline Eigen::Vector3d getvn(const Eigen::Vector3d &normal, const particleGroup &particle,int indPart1,int indPart2){
            Eigen::Vector3d vn;
            double vn_mag;

            vn_mag = (particle.getVel().col(indPart1)-particle.getVel().col(indPart2)).dot(normal);
            vn = vn_mag * normal;

            return vn;
        }

        inline Eigen::Vector3d getvn(const Eigen::Vector3d &normal, const particleGroup &particle1,int indPart, const triangles &wall, int indWall){
            Eigen::Vector3d vn;
            double vn_mag;

            vn_mag = (particle1.getVel().col(indPart)-wall.v.col(indWall)).dot(normal);
            vn = vn_mag * normal;

            return vn;
        }

        inline  Eigen::Vector3d getvt(const Eigen::Vector3d &vn, const Eigen::Vector3d &normal, const particleGroup &particle, int p1Ind, int p2Ind){
            Eigen::Vector3d vt;
            Eigen::Vector3d rot;
            Eigen::Vector3d vr; //relative velocity

            vr = particle.getVel().col(p1Ind) - particle.getVel().col(p2Ind);

            rot = particle.radius()(p1Ind)*particle.getW().col(p1Ind) + particle.radius()(p2Ind)*particle.getW().col(p2Ind);
            rot = rot.cross(normal) ; //cross product
            vt = vr - vn + rot;
            return vt;
        }

        inline  Eigen::Vector3d getvt(const Eigen::Vector3d &vn,const Eigen::Vector3d &normal, const particleGroup &particle, int p1Ind, const triangles &wall, int wallInd){
            Eigen::Vector3d vt;
            Eigen::Vector3d rot;

            Eigen::Vector3d vr; //relative velocity

            vr = particle.getVel().col(p1Ind) - wall.v.col(wallInd);

            rot = particle.radius()(p1Ind)*particle.getW().col(p1Ind);
            rot = rot.cross(normal) ; //cross product


            vt = vr - vn + rot;

            return vt;
        }

        inline int numPart(){
            return particles_.numPart();
        }

        void showNumPartInCell();

        /* temporarly just count numbers in cell */
        void getParticleInCellList();
        void getParticleInCellListTest();
        void fillInStructuredCells();

        void shuffle(int seed, double minX_in, double maxX_in, double minY_in, double maxY_in, double minZ_in, double maxZ_in);

        void getDistList();

        void refreshLinkedList();
        void printCellIniPart();
        void printLinkedList();

        void getDist();
        inline double distSq(const Eigen::Vector3d &part1, const Eigen::Vector3d &part2){
            return (part1-part2).squaredNorm();
        }

        inline std::vector<int>& getLinkedList(){
            return linkedList_;
        }

        inline std::vector<int>& getCellIniPart(){
            return cellIniPart_;
        }


        inline const particleGroup &getParticle(){
            return particles_;
        }

        inline const boundingBox &getBoundingBox() const{
            return box_;
        }

        inline const triangles &getWall(){
            return boundaries_;
        }

        inline triangles &getWriteWall(){
            return boundaries_;
        }

        inline int getNumTriangles(){
            return boundaries_.numTriangles;
        }

        inline particleGroup &getWriteParticle(){
            return particles_;
        }

        inline void initLinkedList(){
            linkedList_ = std::vector<int>(this->numPart(),-1);
            std::cout << "Linked list initialized" << std::endl;
        }

        inline void setParameters(double k,double eta, double mu){
            k_ = k;
            eta_ = eta;
            mu_ = mu;
        }

        inline void setRefreshFreq(double freq){
            refreshFreq_ = freq;
            refreshRadius_ = this->refreshFreq_*this->minDiam_;
            refreshThreshold_ = 0.5*this->refreshRadius_;
        }

        inline const double getRefreshThreshold() const{
            return refreshThreshold_;
        }

        inline void gravity3D(double x,double y, double z){
            gravSource_(0) = x;
            gravSource_(1) = y;
            gravSource_(2) = z;
        }

    private:
        boundingBox box_;
        triangles boundaries_;
        particleGroup particles_;
        std::vector<int> numPartInCell_;
        std::vector<int> cellIniPart_;
        std::vector<int> linkedList_;
        Eigen::Vector3d gravSource_;
        double dt_;
        double k_;
        double eta_;
        double mu_;
        double boxDx_,boxDy_,boxDz_;
        double boxDxInv_,boxDyInv_,boxDzInv_;
        double boxMinX_,boxMinY_,boxMinZ_;
        double minDiam_;
        double maxRad_;
        double refreshFreq_;
        double refreshThreshold_;
        double refreshRadius_;

};



#endif
