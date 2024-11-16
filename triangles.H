#ifndef _triangles_H_
#define _triangles_H_

#include <vector>
#include <iostream>
#include <Eigen/Dense>

class triangles{
    public:
        triangles();
        triangles(int numTriangles);
        int dim = 3;
        int numTriangles = 0;
        Eigen::Matrix3Xd e; //edges
        Eigen::Matrix3Xd f; //faces which consists of verticies coordinates
        int numCollVorE = 0; //number which collided with vertex or edge. used later to avoid double counting

        int numVertex = 3;
        
        Eigen::MatrixXd lengthSqInv;
        Eigen::MatrixXd d; //face's d, when considered as a plane ax+by+cz+d=0
        Eigen::Matrix3Xd v;
        Eigen::Matrix3Xd fNormal;
        Eigen::Matrix3Xd eNormal; //normal of edge out of the triangle and is perpendicular to fNormal
        Eigen::Matrix3Xd collVorEHist;

        void getEdgeInfo();
        Eigen::Vector3d getClosestPt(const Eigen::Vector3d &p, int wallInd);
        inline bool isOutThePlane(const Eigen::Vector3d &v0, const Eigen::Vector3d &normal,const Eigen::Vector3d &p) const;
        inline bool isBetweenZeroToOne(const double a) const;
        inline double projectEdge(const Eigen::Vector3d &v0, const Eigen::Vector3d &e0, double lengthSqInv, const Eigen::Vector3d &p) const;
};


inline bool triangles::isBetweenZeroToOne(const double a) const{
    if (a<0. || a>1.){
        return false;
    }else{
        return true;
    }
}

inline bool triangles::isOutThePlane(const Eigen::Vector3d &v0, const Eigen::Vector3d &normal,const Eigen::Vector3d &p) const{
    Eigen::Vector3d v0p = p-v0;

    if (v0p.dot(normal)<0){
        return false;
    }else{
        return true;
    }
}

inline double triangles::projectEdge(const Eigen::Vector3d &v0, const Eigen::Vector3d &e0, double lengthSqInv, const Eigen::Vector3d &p) const{
    return (p-v0).dot(e0)*lengthSqInv;
}

#endif
