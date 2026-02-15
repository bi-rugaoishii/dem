#include "demCalc.H"
#include "particle.H"
#include "boundingBox.H"
#include "triangles.H"
#include "particleGroup.H"
#include <omp.h>
#include <vector>
#include <iostream>
#define smallNum 1.e-16
#define MAX_HISTORY_NEIGH 50
#define MAX_WALL_EorV_HISTORY 10

demCalc::demCalc(){}
demCalc::demCalc(boundingBox &box,particleGroup &particles){
    particles_= particles;
    box_ = box;
    std::vector<int> a(box_.sizeSplitY()*box_.sizeSplitX()*box_.sizeSplitZ(),0);
    std::vector<int> b(box_.sizeSplitY()*box_.sizeSplitX()*box_.sizeSplitZ(),-1);
    numPartInCell_ = a;
    cellIniPart_ = b;
}

demCalc::demCalc(boundingBox box,particleGroup particles, triangles boundaries){
    boundaries_ = boundaries;
    particles_= particles;
    box_ = box;
    std::vector<int> a(box_.sizeSplitY()*box_.sizeSplitX()*box_.sizeSplitZ(),0);
    std::vector<int> b(box_.sizeSplitY()*box_.sizeSplitX()*box_.sizeSplitZ(),-1);
    numPartInCell_ = a;
    cellIniPart_ = b;
    boxDx_ = box_.dx();
    boxDy_ = box_.dy();
    boxDz_ = box_.dz();
    boxDxInv_ = 1./boxDx_;
    boxDyInv_ = 1./boxDy_;
    boxDzInv_ = 1./boxDz_;
    boxMinX_ = box_.minX();
    boxMinY_ = box_.minY();
    boxMinZ_ = box_.minZ();
    minDiam_ = particles_.radius().minCoeff()*2.;
    maxRad_ = particles_.radius().maxCoeff()*2.;
    refreshFreq_ = 1.0;
    refreshRadius_ = refreshFreq_*minDiam_;
    refreshThreshold_ = 0.5*refreshFreq_*minDiam_;
}



void demCalc::getAcc(double dt){
    int numPart = this->numPart();
    int x,y,z,ind,indHistory;
    int wallNum = this->getNumTriangles();
    int wallInd;
    double dist;
    bool hasCollidedBefore,hasCollidedVorE;
    double vtNorm,deltatOldNorm,ftNorm,fnNorm;
    Eigen::Vector3d mom(0.,0.,0.);
    boundingBox const &box = this->getBoundingBox();
    Eigen::Vector3d normal,delta,vn,fn;
    Eigen::Vector3d closestPt;
    Eigen::Vector3d const zeroVec(0.,0.,0.);
    Eigen::Vector3d deltat(0.,0.,0.),vt(0.,0.,0.),ft(0.,0.,0);
    Eigen::Vector3d tangentialUnit(0.,0.,0.);
    this->particles_.acc.setZero();
    this->particles_.aw.setZero();

    particleGroup const &particle = this->getParticle();
    triangles const &walls = this->getWall();
    const Eigen::MatrixXd &tmpDeltat = particles_.getDeltatHistory();
    #pragma omp parallel
    {
        #pragma omp for private(ind,dist,hasCollidedBefore,indHistory,normal,delta,vn,fn,vt,vtNorm,deltatOldNorm,deltat,ft,ftNorm,fnNorm,mom,tangentialUnit,wallInd,closestPt,hasCollidedVorE)
        for (int i=0; i<numPart; i++){
            //initialize history related arrays




            int const p1Ind = i;
            int nInd = 0;
            int newIndHistory=0;

            ind = particle.getNeighborList()(p1Ind, nInd);
            //get particles in neighboring cells
            while (ind!=-1){

                int p2Ind = ind;

                double distSq = this->getDistSq(particles_,p1Ind,p2Ind);
                // check if collided
                double sumRadius = particles_.radius()(p1Ind)+particles_.radius()(p2Ind);
                double sumRadiusSq = sumRadius * sumRadius;


                if (distSq < sumRadiusSq){

                    dist = sqrt(distSq);
                    //std::cout << i << " and "<< ind <<" collided with "<< dist << "!" << std::endl;
                    //normal force
                    normal = this->getNorm(particles_,p1Ind,p2Ind,dist);
                    delta = this->getDelta(normal,particles_,p1Ind,p2Ind,dist);
                    //check overlap
                    if (delta.norm()/particles_.radius()(p1Ind) > 0.05){
                        std::cout << "particle overlap over 5 %!!!!!!" << std::endl;
                    }

                    //std::cout << "normal = " << delta(0) << " " << delta(1) <<" "<<delta(2) << std::endl;
                    vn = this->getvn(normal,particles_,p1Ind,p2Ind);
                    // std::cout << "vn = " << vn(0) << " " << vn(1) <<" "<<vn(2) << std::endl;
                    fn = this->getNormalForce(delta,vn,this->eta_);
                    particles_.acc.col(p1Ind)= particles_.acc.col(p1Ind) + fn;

                    //std::cout << "fn = " << fn(0) << " " << fn(1) <<" "<<fn(2) << std::endl;


                    //tangential Force Calculation
                    vt = this->getvt(vn,normal,particles_,p1Ind,p2Ind);

                    //write in the history and get tangential displacement
                    //get tangential displacement

                    vtNorm = vt.norm();
                    hasCollidedBefore = false;
                    indHistory = -1;

                    for (int j=0; j<particles_.getNumCollisionHistory()(p1Ind); j++){
                        if (particles_.collisionHistory(p1Ind,j) == p2Ind){
                            indHistory=j;
                            deltatOldNorm = tmpDeltat.col(p1Ind*MAX_HISTORY_NEIGH+indHistory).norm();
                            tangentialUnit = vt/(vtNorm+smallNum);
                            deltat = tangentialUnit*deltatOldNorm;
                            deltat = deltat+vt*dt;

                            particles_.deltatHistoryNew.col(p1Ind*MAX_HISTORY_NEIGH+newIndHistory)=deltat;
                            particles_.collisionHistoryNew(p1Ind*MAX_HISTORY_NEIGH+newIndHistory) = p2Ind;
                            newIndHistory += 1;
                            hasCollidedBefore = true;
                            break;
                        }
                    }

                    if (hasCollidedBefore == false){
                        deltat = vt*dt;
                        tangentialUnit = vt/(vtNorm+smallNum);
                        particles_.deltatHistoryNew.col(p1Ind*MAX_HISTORY_NEIGH+newIndHistory)=deltat;
                        particles_.collisionHistoryNew(p1Ind*MAX_HISTORY_NEIGH+newIndHistory) = p2Ind;
                        newIndHistory += 1;
                    }


                    ft = this->getTangentialForce(deltat,vt,this->eta_);
                    ftNorm = ft.norm();
                    fnNorm = fn.norm();

                    if (ftNorm > mu_*fnNorm){
                        if (hasCollidedBefore == true){
                            particles_.deltatHistoryNew.col(p1Ind*MAX_HISTORY_NEIGH+newIndHistory-1)=tmpDeltat.col(p1Ind*MAX_HISTORY_NEIGH+indHistory);
                            ft = -mu_*fnNorm*tangentialUnit;
                        }else{
                            particles_.deltatHistoryNew.col(p1Ind*MAX_HISTORY_NEIGH+newIndHistory-1)=zeroVec;
                            ft = -mu_*fnNorm*tangentialUnit;
                        }
                    }



                    //will be divided by mass later
                    particles_.acc.col(p1Ind) = particles_.acc.col(p1Ind) + ft;



                    //add rotation
                    mom = getMoment(normal,particles_,p1Ind,ft);
                    //will be divided by moment of intertia later
                    particles_.aw.col(p1Ind) = particles_.aw.col(p1Ind) + mom;



                    //  std::cout << "mom = " << mom << std::endl;


                }// end if (dist < particle1.radius()+particle2.radius())

                nInd +=1;
                ind = particle.getNeighborList()(p1Ind, nInd);
            }//done particle collisions

            particles_.numCollisionHistory(p1Ind) = newIndHistory;



            //check for wall collisions

            newIndHistory = 0;
            const Eigen::MatrixXd &tmpDeltatWall = particles_.getDeltatHistoryWall();
            //this->getWriteWall().numCollVorE=0;
            int numCollVorE=0;
            Eigen::Matrix3Xd collVorEHist(3,MAX_WALL_EorV_HISTORY); //3 is dimension. Memory allocated here to avoid parallization conflict. Maybe include it as particleGroupMember?
            std::vector<int> histCollWall(box_.num_max_triangle,-1); //Memory allocated here to avoid parallization conflict. Maybe include it as particleGroupMember?
            collVorEHist.setZero();
            int countHitWall=0;

            //get list of neighbor cells
            int tmpx = (int)((particles_.getPos()(0,i) - boxMinX_)*boxDxInv_)+1;
            int tmpy = (int)((particles_.getPos()(1,i) - boxMinY_)*boxDyInv_)+1;
            int tmpz = (int)((particles_.getPos()(2,i) - boxMinZ_)*boxDzInv_)+1;
            for (int z= tmpz-1; z<tmpz+2; z++ ){
                for (int y=tmpy-1; y<tmpy+2; y++){
                    for (int x=tmpx-1; x<tmpx+2; x++){
                        int cellInd = index(z,y,x,box_.sizeSplitX(),box_.sizeSplitY());
                        int wallCountInCell = 0;
                        //get walls in neighborlist
                        int wallInd = box_.triangleList[cellInd*box_.num_max_triangle+wallCountInCell];
                        //std::cout << "wallInd is " << wallInd<< std::endl;
                            while(wallInd!=-1){
                                bool hasCollidedWithTheWallBefore=false;
                                for (int wallCounter =0; wallCounter<countHitWall; wallCounter++){
                                    if(histCollWall[wallCounter]==wallInd){
                                        hasCollidedWithTheWallBefore=true;
                                        break;
                                           // std::cout << "was hit before!" << std::endl;
                                    }
                                }
                                if(hasCollidedWithTheWallBefore == true){
                                    //std::cout <<"continue "<< std::endl;
                                    wallCountInCell+=1;
                                    wallInd = box_.triangleList[cellInd*box_.num_max_triangle+wallCountInCell]  ;//check next wall in cell
                                    //std::cout << "wallInd refreshed as " << wallInd<< std::endl;
                                    continue;
                                }
                                const Eigen::Vector3d &particle1Pos = particles_.getPos().col(p1Ind);
                                closestPt = this->getWriteWall().getClosestPt(particle1Pos,wallInd);
                                //std::cout << hasCollidedVorE << std::endl;


                                normal = particle1Pos-closestPt; //will be normalized later
                                double distSq = (normal).squaredNorm();
                                if (distSq < particle.radiusSq()(p1Ind)){
                                    //was Hit
                                    histCollWall[countHitWall]=wallInd;
                                    countHitWall +=1;


                                    //check if the closestPt is double counting
                                    bool hasCollidedBeforeVorE = false;
                                    for (int k = 0; k<numCollVorE; k++){
                                        double distSqVorE = (closestPt - collVorEHist.col(k)).squaredNorm();
                                        //                    std::cout << "distSqVorE was " << distSqVorE << std::endl;
                                        if(distSqVorE <smallNum){
                                            std::cout << "had repetition!" << std::endl;
                                            hasCollidedBeforeVorE = true;
                                            continue;
                                        }
                                    }
                                    if(hasCollidedBeforeVorE == true){
                                        //std::cout <<"continue "<< std::endl;
                                        wallCountInCell+=1;
                                        wallInd = box_.triangleList[cellInd*box_.num_max_triangle+wallCountInCell]  ;
                                        //std::cout << "wallInd refreshed as " << wallInd<< std::endl;
                                        continue;
                                    }
                                    collVorEHist.col(numCollVorE) = closestPt;
                                    numCollVorE += 1;
                                    //this->getWriteWall().numCollVorE += 1;
                                    /*
                                       std::cout << "wallNum is " << wallInd << std::endl;
                                       std::cout << "particlePos is " << particle1Pos.transpose() << std::endl;
                                       std::cout << "closestPt is " << closestPt.transpose() << std::endl;
                                       std::cout <<  std::endl;
                                       */

                                    dist = sqrt(distSq);
                                    //normal force
                                    normal.normalize();
                                    delta = this->getDelta(normal,particle,p1Ind,walls,wallInd,dist);
                                    vn = this->getvn(normal,particle,p1Ind,walls,wallInd);
                                    fn = this->getNormalForce(delta,vn,this->etaWall_);

                                    //will be divided by mass later
                                    particles_.acc.col(p1Ind) = particles_.acc.col(p1Ind) + fn;

                                    //tangential Force Calculation
                                    vt = this->getvt(vn,normal,particles_,p1Ind,walls, wallInd);

                                    //write in the history and get tangential displacement
                                    //get tangential displacement


                                    vtNorm = vt.norm();
                                    hasCollidedBefore = false;
                                    indHistory = -1;
                                    for (int j=0; j<particles_.numCollisionHistoryWall(p1Ind); j++){
                                        if (particles_.collisionHistoryWall(p1Ind,j) == wallInd){
                                            indHistory=j;
                                            deltatOldNorm = tmpDeltatWall.col(p1Ind*MAX_HISTORY_NEIGH+indHistory).norm();
                                            tangentialUnit = vt/(vtNorm+smallNum);
                                            deltat = tangentialUnit*deltatOldNorm;
                                            deltat = deltat+vt*dt;

                                            particles_.deltatHistoryNewWall.col(p1Ind*MAX_HISTORY_NEIGH+newIndHistory)=deltat;
                                            particles_.collisionHistoryNewWall(p1Ind*MAX_HISTORY_NEIGH+newIndHistory) = wallInd;
                                            newIndHistory += 1;
                                            hasCollidedBefore = true;
                                            break;
                                        }
                                    }

                                    if (hasCollidedBefore == false){
                                        deltat = vt*dt;
                                        tangentialUnit = vt/(vtNorm+smallNum);
                                        particles_.deltatHistoryNewWall.col(p1Ind*MAX_HISTORY_NEIGH+newIndHistory)=deltat;
                                        particles_.collisionHistoryNewWall(p1Ind*MAX_HISTORY_NEIGH+newIndHistory) = wallInd;
                                        newIndHistory += 1;
                                    }


                                    ft = this->getTangentialForce(deltat,vt,this->etaWall_);


                                    ftNorm = ft.norm();
                                    fnNorm = fn.norm();


                                    if (ftNorm > mu_*fnNorm){
                                        if (hasCollidedBefore == true){
                                            particles_.deltatHistoryNewWall.col(p1Ind*MAX_HISTORY_NEIGH+newIndHistory-1)=tmpDeltatWall.col(p1Ind*MAX_HISTORY_NEIGH+indHistory);
                                            ft = -mu_*fnNorm*tangentialUnit;
                                        }else{
                                            particles_.deltatHistoryNewWall.col(p1Ind*MAX_HISTORY_NEIGH+newIndHistory-1)=zeroVec;
                                            ft = -mu_*fnNorm*tangentialUnit;
                                        }
                                    }



                                    //will be divided by mass later
                                    particles_.acc.col(p1Ind) = particles_.acc.col(p1Ind) + ft;



                                    //add rotation
                                    mom = getMoment(normal,particles_,p1Ind,ft);
                                    //will be divided by moment of intertia later
                                    particles_.aw.col(p1Ind) = particles_.aw.col(p1Ind) + mom;



                                }
                                wallCountInCell+=1;
                                wallInd = box_.triangleList[cellInd*box_.num_max_triangle+wallCountInCell];// check next wall in cell
                                //std::cout << "wallInd refreshed as " << wallInd<< std::endl;
                            }
                    }
                }
            }

            particles_.numCollisionHistoryNewWall(p1Ind) = newIndHistory;

            //update the history





            //std::cout << "acc = " << demCalc1.getParticle(i).acc[0] <<" " << demCalc1.getParticle(i).acc[1]<< " " << demCalc1.getParticle(i).acc[2] << std::endl;
        }
    }
    particles_.deltatHistory = particles_.deltatHistoryNew;
    particles_.collisionHistory = particles_.collisionHistoryNew;
    particles_.deltatHistoryWall = particles_.deltatHistoryNewWall;
    particles_.collisionHistoryWall = particles_.collisionHistoryNewWall;
    //add source term and divide by mass inverse

    for (int p1Ind=0; p1Ind<numPart; p1Ind++){
        particles_.acc.col(p1Ind) = particles_.getAcc().col(p1Ind)*particles_.massInv()(p1Ind) + gravSource_;
        particles_.aw.col(p1Ind) = particles_.aw.col(p1Ind)*particles_.moiInv()(p1Ind);
    }
}

void demCalc::showNumPartInCell(){
    for (int k=0; k<box_.sizeSplitZ(); k++){
        for (int i=0; i<box_.sizeSplitY(); i++){
            for (int j=0; j<box_.sizeSplitX(); j++){
                std::cout << numPartInCell_[index(box_.sizeSplitZ()-1-k,box_.sizeSplitY()-1-i,j,box_.sizeSplitX(),box_.sizeSplitY())] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}

void demCalc::fillInStructuredCells(){
    int numPart = this->numPart();
    int partCounter = 0;
    double particleRad = particles_.radius()(0);
    double diam = particles_.radius()(0)*2.0;

    for (int k=1; k<box_.sizeSplitZ()-1; k++){
        for (int i=1; i<box_.sizeSplitY()-1; i++){
            for (int j=1; j<box_.sizeSplitX()-1; j++){
                particles_.pos(0,partCounter)=(double)(j-1)*diam+box_.minX()+particleRad;
                particles_.pos(1,partCounter)=(double)(i-1)*diam+box_.minY()+particleRad;
                particles_.pos(2,partCounter)=(double)(k-1)*diam+box_.minZ()+particleRad;
                partCounter++;
                if (partCounter>=numPart){
                    break;
                }
            }
            if (partCounter>=numPart){
                break;
            }

        }
    }
    std::cout << "Filled the box with particles!" << std::endl;

}

/* temporarly just count numbers in cell */
void demCalc::getParticleInCellList(){
    int numPart = this->numPart();
    int tmpx,tmpy,tmpz;
    for (int i=0; i<numPart; i++){
        tmpx = (int)((particles_.getPos()(0,i) - box_.minX())/box_.dx())+1;
        tmpy = (int)((particles_.getPos()(1,i) - box_.minY())/box_.dy())+1;
        tmpz = (int)((particles_.getPos()(2,i) - box_.minZ())/box_.dz())+1;
        numPartInCell_[index(tmpz,tmpy,tmpx,box_.sizeSplitX(),box_.sizeSplitY())]+=1;
    }
    std::cout << "counted the number of particles in a cell " << std::endl;
}

void demCalc::getParticleInCellListTest(){
    int ind;
    for (int k=0; k<box_.sizeSplitZ(); k++){
        for (int i=0; i<box_.sizeSplitY(); i++){
            for (int j=0; j<box_.sizeSplitX(); j++){
                ind = index(k,i,j,box_.sizeSplitX(),box_.sizeSplitY());
                numPartInCell_[ind]=ind;
            }
        }
    }
    std::cout << "Got the index of cell array" << std::endl;
}


void demCalc::shuffle(int seed,double minX_in,double maxX_in, double minY_in, double maxY_in, double minZ_in, double maxZ_in){
    double distSq,tmpPosx,tmpPosy,tmpPosz;
    int numPart = this->numPart();
    int placedNum = 0;
    double minX,maxX,minY,maxY, minZ, maxZ; //possible position without collision with the box walls
    std::mt19937 mt(seed);
    std::uniform_real_distribution<double> rand01(0.,1.);

    std::cout<<"Shuffling"<<std::endl;
    // for initial particle
    minX = minX_in+particles_.radius()(0);
    minY = minY_in+particles_.radius()(0);
    minZ = minZ_in+particles_.radius()(0);
    maxX = maxX_in-particles_.radius()(0);
    maxY = maxY_in-particles_.radius()(0);
    maxZ = maxZ_in-particles_.radius()(0);
    tmpPosx = rand01(mt)*(maxX-minX)+minX;
    tmpPosy = rand01(mt)*(maxY-minY)+minY;
    tmpPosz = rand01(mt)*(maxZ-minZ)+minZ;
    particles_.pos(0,0)=tmpPosx;
    particles_.pos(1,0)=tmpPosy;
    particles_.pos(2,0)=tmpPosz;
    placedNum = 1;
    std::cout<<"Placed the initial"<<std::endl;

    for (int i=1; i<numPart; i++){
        //set the range so that a particle doesn't hit the wall
        minX = minX_in+particles_.radius()(i);
        minY = minY_in+particles_.radius()(i);
        minZ = minZ_in+particles_.radius()(i);
        maxX = maxX_in-particles_.radius()(i);
        maxY = maxY_in-particles_.radius()(i);
        maxZ = maxZ_in-particles_.radius()(i);

        distSq = -1.;
        //check the distance between other particles
        tmpPosx = rand01(mt)*(maxX-minX)+minX;
        tmpPosy = rand01(mt)*(maxY-minY)+minY;
        tmpPosz = rand01(mt)*(maxZ-minZ)+minZ;

        for (int j=0; j<placedNum; j++){
            Eigen::Vector3d tmpPos(tmpPosx,tmpPosy,tmpPosz);
            if (i==j){
                continue;
            }else{
                distSq = (tmpPos-particles_.pos.col(j)).squaredNorm();
                //    std::cout << distSq << std::endl;
                if (distSq<pow(particles_.radius()(i)+particles_.radius()(j),2)){
                    tmpPosx = rand01(mt)*(maxX-minX)+minX;
                    tmpPosy = rand01(mt)*(maxY-minY)+minY;
                    tmpPosz = rand01(mt)*(maxZ-minZ)+minZ;
                    j=-1; //reroll and restart the loop if the particles collide
                }
            }
        }

        particles_.pos(0,i)=tmpPosx;
        particles_.pos(1,i)=tmpPosy;
        particles_.pos(2,i)=tmpPosz;
        placedNum++;
        std::cout<<"placed"<<std::endl;
    }
} //end void shuffle() 

void demCalc::getDistList(){
    int numPart = this->numPart();
    double dist;
    for (int i=0; i<numPart; i++){
        for (int j=i; j<numPart; j++){
            if (i != j){
                dist = sqrt(pow(particles_.getPos()(0,i)-particles_.getPos()(0,j),2)+pow(particles_.getPos()(1,i)-particles_.getPos()(1,j),2)+pow(particles_.getPos()(2,i)-particles_.getPos()(2,j),2));
                if (dist < particles_.radius()(i)*2){
                    std::cout << "particle "<<i << " and " << j<< " is too close!"<< std::endl;
                    std::cout << "Dist was " << dist << std::endl;
                }
            }
        }
    }
}


void demCalc::refreshLinkedList(){
    int numPart = this->numPart();
    int ind, tmpx,tmpy,tmpz;
    int indLast;

    //std::cout <<"Refreshing the linked List" << std::endl;

    //initialize
    //
    for (int i=0; i<cellIniPart_.size(); i++){
        cellIniPart_[i] = -1;
    }

    for (int i=0; i<linkedList_.size(); i++){
        linkedList_[i] = -1;
    }

    for (int i=0; i<numPart; i++){
        tmpx = (int)((particles_.getPos()(0,i) - boxMinX_)*boxDxInv_)+1;
        tmpy = (int)((particles_.getPos()(1,i) - boxMinY_)*boxDyInv_)+1;
        tmpz = (int)((particles_.getPos()(2,i) - boxMinZ_)*boxDzInv_)+1;
        //std::cout << "tmpx = " << tmpx << " tmpy = " << tmpy << " tmpz = " << tmpz <<std::endl;
        // std::cout << "posx = " << particles_[i].getPos()[0] << " posy = " << particles_[i].getPos()[1] << " posz = " << particles_[i].getPos()[2] << std::endl;
        ind = index(tmpz,tmpy,tmpx,box_.sizeSplitX(),box_.sizeSplitY());
        if (cellIniPart_[ind] == -1){
            cellIniPart_[ind] = i;
        }else{
            /*
            indLast = cellIniPart_[ind];
            while (linkedList_[indLast]!=-1){
                indLast = linkedList_[indLast];
            }
            linkedList_[indLast] = i;
            */

            linkedList_[i] = cellIniPart_[ind];
            cellIniPart_[ind]=i;
        }
    }

    //refreshing the neighborlist//

    //initialize//

    particles_.writeNeighborList().setConstant(-1);

    for (int i=0; i<numPart; i++){
        int p1Ind = i;
        int neighborCount = 0;

        tmpx = (int)((particles_.getPos()(0,p1Ind) - boxMinX_)*boxDxInv_)+1;
        tmpy = (int)((particles_.getPos()(1,p1Ind) - boxMinY_)*boxDyInv_)+1;
        tmpz = (int)((particles_.getPos()(2,p1Ind) - boxMinZ_)*boxDzInv_)+1;
        for (int z= tmpz-1; z<tmpz+2; z++ ){
            for (int y=tmpy-1; y<tmpy+2; y++){
                for (int x=tmpx-1; x<tmpx+2; x++){
                    ind = this->getCellIniPart()[index(z,y,x,box_.sizeSplitX(),box_.sizeSplitY())];
                    //get particles in neighboring cells
                    while (ind!=-1){
                        if (i!=ind){
                            int p2Ind = ind;
                            double distSq = this->getDistSq(particles_,p1Ind,p2Ind);
                            // check if in range
                            double sumRadius = particles_.radius()(p1Ind)+particles_.radius()(p2Ind)+refreshRadius_;
                            double sumRadiusSq = sumRadius * sumRadius;

                            if (distSq < sumRadiusSq){
                                //std::cout << "Dist between " << p1Ind << " and " << p2Ind<< " is " << sqrt(distSq) << std::endl;
                                particles_.writeNeighborList()(p1Ind,neighborCount) = p2Ind;
                                neighborCount +=1;
                            }
                        }
                        ind=this->getLinkedList()[ind];
                    }
                }
            }
        }
    }

    particles_.posRef = particles_.getPos();

}

void demCalc::printCellIniPart(){
    int numPart = this->numPart();
    int ind;
    std::cout << "Printing the initial particle List" << std::endl;
    for (int k=0; k<box_.sizeSplitZ(); k++){
        for (int i=0; i<box_.sizeSplitY(); i++){
            for (int j=0; j<box_.sizeSplitX(); j++){
                std::cout << cellIniPart_[index(box_.sizeSplitZ()-1-k,box_.sizeSplitY()-1-i,j,box_.sizeSplitX(),box_.sizeSplitY())] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}

void demCalc::printLinkedList(){
    int numPart = this->numPart();
    int ind;
    std::cout << "Printing the linked List" << std::endl;
    for (int i=0; i<numPart; i++){
        std::cout << "["<<i<<","<<linkedList_[i] << "] " ;
    }
    std::cout << std::endl;
}

void demCalc::getDist(){
    int numPart = this->numPart();
    int i,j;
    int x,y,z;
    int tmpx,tmpy,tmpz;
    int ind;
    double dist;
    std::cout << "Getting the distance" << std::endl;
    for (int i=0; i<numPart; i++){
        tmpx = (int)((particles_.getPos()(i,0) - box_.minX())/box_.dx())+1;
        tmpy = (int)((particles_.getPos()(i,1) - box_.minY())/box_.dy())+1;
        tmpz = (int)((particles_.getPos()(i,2) - box_.minZ())/box_.dz())+1;
        for (z=tmpz-1; z<tmpz+2; z++){
            for (y=tmpy-1; y<tmpy+2; y++){
                for (x=tmpx-1; x<tmpx+2; x++){
                    ind = cellIniPart_[index(z,y,x,box_.sizeSplitX(),box_.sizeSplitY())];
                    //get particles in neighboring cells
                    while (ind!=-1){
                        if (i!=ind){
                            dist = (particles_.pos.col(i)-particles_.pos.col(ind)).norm();
                            std::cout << "For " << i <<"," <<ind<< " distance is "<< dist<<"."<<std::endl;
                        }
                        ind=linkedList_[ind];
                    }
                }
            }
        }
    }
}


