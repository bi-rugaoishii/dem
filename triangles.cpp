#include "triangles.H"

triangles::triangles(){
}


triangles::triangles(int const numTriangles){
    this->numTriangles=numTriangles;
    lengthSqInv.resize(numTriangles*numVertex,1);
    lengthSqInv.setZero();
    fNormal.resize(this->dim,numTriangles);
    fNormal.setZero();
    eNormal.resize(this->dim,numTriangles*numVertex);
    eNormal.setZero();
    f.resize(this->dim,numTriangles*numVertex);
    f.setZero();
    e.resize(this->dim,numTriangles*numVertex);
    e.setZero();
    v.resize(this->dim,numTriangles);
    v.setZero();
    d.resize(numTriangles,1);
    d.setZero();
    minX.resize(numTriangles,0.);
    minY.resize(numTriangles,0.);
    minZ.resize(numTriangles,0.);
    maxX.resize(numTriangles,0.);
    maxY.resize(numTriangles,0.);
    maxZ.resize(numTriangles,0.);

    int numHistory = 10;
    collVorEHist.resize(3,numHistory);
    collVorEHist.setZero();
}

Eigen::Vector3d triangles::getClosestPt(const Eigen::Vector3d &p, int wallInd){
    Eigen::Vector3d result(0.,0.,0.);
    
   const Eigen::Vector3d &v0 = this->f.col(wallInd*this->numVertex+0); 
   const Eigen::Vector3d &v1 = this->f.col(wallInd*this->numVertex+1); 
   const Eigen::Vector3d &v2 = this->f.col(wallInd*this->numVertex+2); 
   const Eigen::Vector3d &e0 = this->e.col(wallInd*this->numVertex+0); 
   const Eigen::Vector3d &e1 = this->e.col(wallInd*this->numVertex+1); 
   const Eigen::Vector3d &e2 = this->e.col(wallInd*this->numVertex+2); 
   const Eigen::Vector3d &eNormal0 = this->eNormal.col(wallInd*this->numVertex+0); 
   const Eigen::Vector3d &eNormal1 = this->eNormal.col(wallInd*this->numVertex+1); 
   const Eigen::Vector3d &eNormal2 = this->eNormal.col(wallInd*this->numVertex+2); 
   const Eigen::Vector3d &normalP = this->fNormal.col(wallInd); 
   const double &lengthSqInv0 = this->lengthSqInv(wallInd*this->numVertex+0);
   const double &lengthSqInv1 = this->lengthSqInv(wallInd*this->numVertex+1);
   const double &lengthSqInv2 = this->lengthSqInv(wallInd*this->numVertex+2);
   const double &d = this->d(wallInd);




   double pro01 = this->projectEdge(v0,e0,lengthSqInv0,p);
   double pro12 = this->projectEdge(v1,e1,lengthSqInv1,p);




   if (pro01 > 1. && pro12 < 0.){
       return v1;
   }

   double pro20 = this->projectEdge(v2,e2,lengthSqInv2,p);
   if (pro12 > 1. && pro20 < 0.){
       return v2;
   }

   if (pro20 > 1. && pro01 < 0.){
       return v0;
   }

   if (this->isBetweenZeroToOne(pro01) && this->isOutThePlane(v0,eNormal0,p)){
       return v0+pro01*e0;
   }

   if (this->isBetweenZeroToOne(pro12) && this->isOutThePlane(v1,eNormal1,p)){
       return v1+pro12*e1;
   }

   if (this->isBetweenZeroToOne(pro20) && this->isOutThePlane(v2,eNormal2,p)){
       return v2+pro20*e2;
   }

   if ((p-v0).dot(normalP)>0){
       return p-fabs(p.dot(normalP)+d)*normalP;
   }else{
       return p+fabs(p.dot(normalP)+d)*normalP;
   }

}

void triangles::getEdgeInfo(){
    for (int triInd=0; triInd< this->numTriangles; triInd++){
        int ind0 = triInd*3+0;
        int ind1 = triInd*3+1;
        int ind2 = triInd*3+2;

        
        //make bounding box of the triangle
        this->minX[triInd] = std::min({this->f.col(ind0)(0), this->f.col(ind1)(0),this->f.col(ind2)(0)});  
        this->minY[triInd] = std::min({this->f.col(ind0)(1), this->f.col(ind1)(1),this->f.col(ind2)(1)});  
        this->minZ[triInd] = std::min({this->f.col(ind0)(2), this->f.col(ind1)(2),this->f.col(ind2)(2)});  
        this->maxX[triInd] = std::max({this->f.col(ind0)(0), this->f.col(ind1)(0),this->f.col(ind2)(0)});  
        this->maxY[triInd] = std::max({this->f.col(ind0)(1), this->f.col(ind1)(1),this->f.col(ind2)(1)});  
        this->maxZ[triInd] = std::max({this->f.col(ind0)(2), this->f.col(ind1)(2),this->f.col(ind2)(2)});  

        this->e.col(ind0) = this->f.col(ind1) - this->f.col(ind0);  

        this->e.col(ind1) = this->f.col(ind2) - this->f.col(ind1);  

        this->e.col(ind2) = this->f.col(ind0) - this->f.col(ind2);  

        /*
           this->fNormal.col(triInd) = this->e.col(ind0).cross(this->e.col(ind1));
           this->fNormal.col(triInd).normalize(); 
           */

        this->eNormal.col(ind0) = this->e.col(ind0).cross(this->fNormal.col(triInd));
        this->eNormal.col(ind0).normalize();

        this->eNormal.col(ind1) = this->e.col(ind1).cross(this->fNormal.col(triInd));
        this->eNormal.col(ind1).normalize();

        this->eNormal.col(ind2) = this->e.col(ind2).cross(this->fNormal.col(triInd));
        this->eNormal.col(ind2).normalize();

        this->lengthSqInv(ind0) = 1./this->e.col(ind0).squaredNorm();
        this->lengthSqInv(ind1) = 1./this->e.col(ind1).squaredNorm();
        this->lengthSqInv(ind2) = 1./this->e.col(ind2).squaredNorm();
        this->d(triInd) = -1.*this->f.col(ind0).dot(this->fNormal.col(triInd));
    }
}



