#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>
#include <cmath>
#include <string>
#define PI 3.14159265359
#define smallNum 1.e-16
#include "particle.H"
#include "misc.H"
#include "forceCalculation.H"
#include "vector3d.H"
#include <Eigen/Dense>

forceCalculation::forceCalculation(){};

Eigen::Vector3d forceCalculation::getNorm(const particleGroup &particle,int indPart, const planeBoundary &wall, int indWall){
    Eigen::Vector3d resultNormal;
    double innerProd;
    Eigen::Vector3d v_tmp;
    double sign =1. ;

    v_tmp = particle.getPos().col(indPart) - wall.getPos().col(indWall);
    innerProd = v_tmp.dot(wall.getNormal().col(indWall));

    if (innerProd < 0){
        sign = -1.;
    }

    resultNormal = sign*wall.getNormal().col(indWall);

    return resultNormal;
}



void forceCalculation::getAcc(demCalc &demCalc1,const double &dt){
    int numPart = demCalc1.numPart();
    int x,y,z,ind,indHistory;
    int wallNum = demCalc1.getNumWalls();
    int wallInd;
    int tmpx, tmpy, tmpz;
    bool hasCollidedBefore;
    double vtsq,vtNorm,deltatOldNorm,tCoeff,ftNorm,fnNorm;
    double dist;
    vector3d mom;
    boundingBox const &box = demCalc1.getBoundingBox();
    Eigen::Vector3d normal,delta,vn,fn;
    vector3d tangential,deltat,vt,ft;
    vector3d tangentialUnit;
    std::vector<vector3d> tmpDeltat;
    std::vector<vector3d> tmpDeltatWall;
    vector3d zeroVec(0.);
//    std::cout << "Calculating the force and the acceleration " << std::endl;
    demCalc1.getWriteParticle().acc.setZero();

    particleGroup const &particle = demCalc1.getParticle();
    planeBoundary const &wall = demCalc1.getWall();
    for (int i=0; i<numPart; i++){
        /*
        //initialize history related arrays
        std::vector<int> historyNew;
        std::vector<vector3d> deltatHistoryNew;

        std::vector<int> historyNewWall;
        std::vector<vector3d> deltatHistoryNewWall;
        
        //initialize acceleration with the source
        vector3d zeroVec(0.);
        demCalc1.getWriteParticle(i).acc = zeroVec;
        demCalc1.getWriteParticle(i).aw= zeroVec;

        particle const &particle1 = demCalc1.getParticle(i);
        tmpx = (int)((particle1.pos[0] - box.minX())/box.dx())+1;
        tmpy = (int)((particle1.pos[1] - box.minY())/box.dy())+1;
        tmpz = (int)((particle1.pos[2] - box.minZ())/box.dz())+1;
        for (z= tmpz-1; z<tmpz+2; z++ ){
            for (y=tmpy-1; y<tmpy+2; y++){
                for (x=tmpx-1; x<tmpx+2; x++){
                    ind = demCalc1.getCellIniPart()[index(z,y,x,box.sizeSplitX(),box.sizeSplitY())];
                    //get particles in neighboring cells
                    while (ind!=-1){
                        if (i!=ind){
                            particle const &particle2 = demCalc1.getParticle(ind);
                            dist = this->getDist(particle1,particle2);
                            // check if collided
                            if (dist < particle1.radius()+particle2.radius()){
                                //std::cout << i << " and "<< ind <<" collided with "<< dist << "!" << std::endl;

                                //normal force
                                normal = this->getNorm(particle1,particle2,dist);
                                delta = this->getDelta(normal,particle1,particle2,dist);
                                vn = this->getvn(normal,particle1,particle2);
                                fn = this->getNormalForce(delta,vn);
                                demCalc1.getWriteParticle(i).acc= demCalc1.getWriteParticle(i).acc + fn;



                                //tangential Force Calculation
                                vt = this->getvt(vn,normal,particle1,particle2);

                                //write in the history and get tangential displacement
                                //get tangential displacement

                                vtNorm = normVec(vt);
                                tmpDeltat = particle1.getDeltatHistory();
                                hasCollidedBefore = false;
                                indHistory = -1;
                                for (int j=0; j<particle1.collisionHistorySize(); j++){
                                    if (particle1.collisionHistory[j] == ind){
                                        indHistory=j;
                                        deltatOldNorm = normVec(tmpDeltat[indHistory]);
                                        tangentialUnit = vt/(vtNorm+smallNum);
                                        deltat = tangentialUnit*deltatOldNorm;
                                        deltat = deltat+vt*dt;
                                        deltatHistoryNew.emplace_back(deltat+tmpDeltat[indHistory]);
                                        historyNew.emplace_back(ind);
                                        hasCollidedBefore = true;
                                        break;
                                    }
                                }

                                if (hasCollidedBefore == false){
                                    deltat = vt*dt;
                                    deltatHistoryNew.emplace_back(deltat);
                                    historyNew.emplace_back(ind);
                                }


                                ft = this->getTangentialForce(deltat,vt);
                                ftNorm = normVec(ft);
                                fnNorm = normVec(fn);

                                if (ftNorm > mu_*fnNorm){
                                    if (hasCollidedBefore == true){
                                        deltatHistoryNew.pop_back();
                                        deltatHistoryNew.emplace_back(tmpDeltat[indHistory]);
                                        ft = -mu_*fnNorm*tangentialUnit;
                                    }else{
                                        deltatHistoryNew.pop_back();
                                        deltatHistoryNew.emplace_back(zeroVec);
                                        ft = zeroVec; //need check
                                    }
                                }



                                //will be divided by mass later
                                demCalc1.getWriteParticle(i).acc = demCalc1.getWriteParticle(i).acc + ft;
                                

                                
                                   std::cout << "deltat = " << deltat[0] << " " << deltat[1] << std::endl;
                                   std::cout << "ft = " << ft[0] << " " << ft[1] << std::endl;
                                   std::cout << "vt = " << vt[0] << " " << vt[1] << std::endl;
                                   */
                                   

                                //add rotation

                                /*
                                mom = getMoment(normal,particle1,ft);
                                //will be divided by moment of intertia later
                                demCalc1.getWriteParticle(i).aw = demCalc1.getWriteParticle(i).aw + mom;

                                //  std::cout << "mom = " << mom << std::endl;


                            }// end if (dist < particle1.radius()+particle2.radius())

                        }
                        ind=demCalc1.getLinkedList()[ind];
                    }
                }
            }
        }
*/


        //check for wall collisions

        for (wallInd=0; wallInd<wallNum; wallInd++){
            dist = this->getDist(particle,i,wall,wallInd);
            if (dist < particle.radius()(i)){
                std::cout << i << " and "<< wallInd <<" collided with "<< dist << "!" << std::endl;
                //normal force
                normal = this->getNorm(particle,i,wall,wallInd);
                delta = this->getDelta(normal,particle,i,wall,wallInd,dist);
                vn = this->getvn(normal,particle,i,wall,wallInd);
                fn = this->getNormalForce(delta,vn);

                demCalc1.getWriteParticle().acc.col(i) = demCalc1.getWriteParticle().acc.col(i) + fn;


            }
        }




        
        //add source term

        demCalc1.getWriteParticle().acc.col(i) = demCalc1.getParticle().getAcc().col(i) + gravSource_;

        //std::cout << "acc = " << demCalc1.getParticle(i).acc[0] <<" " << demCalc1.getParticle(i).acc[1]<< " " << demCalc1.getParticle(i).acc[2] << std::endl;

    }
}

