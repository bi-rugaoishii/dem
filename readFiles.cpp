#include "readFiles.H"


readFiles::readFiles(){
    normal = std::vector<double>(2,0.);
    particles = std::vector<double>(1);
};

triangles readFiles::readStl(std::string filename){
    int dim = 3;
    std::ifstream ifs(filename);

    if (!ifs){
        std::cout << "geometry doesn't exists" << std::endl;
        std::abort();
    }

    int numTriangles = 0;

    std::string str;
    //get number of triangles
    while (ifs >> str){
        if (str == "facet"){
            numTriangles +=1;
        }
    }

    ifs.close();

    std::cout << "number of triangles are " << numTriangles << std::endl;
    triangles walls(numTriangles);

    ifs.open(filename);


    int countTriangles = 0;
    int numVertex = 0;

    while (ifs >> str){
        if (str == "normal"){
            ifs >> walls.fNormal(0,countTriangles) >> walls.fNormal(1,countTriangles) >> walls.fNormal(2,countTriangles) ; 
            countTriangles += 1;
            std::cout << countTriangles<< std::endl;
        }
        if (str == "vertex"){
            int indTriangle = countTriangles-1;
            ifs >> walls.f(0,indTriangle*walls.numVertex+numVertex) >> walls.f(1,indTriangle*walls.numVertex+numVertex) >> walls.f(2,indTriangle*walls.numVertex+numVertex) ; 
            numVertex += 1;
            if (numVertex==walls.numVertex){
                numVertex = 0;
            }
        }
    }

    ifs.close();
    return walls;
}
void readFiles::readSetting(){
    std::ifstream ifs("demSettings.txt");

    
    
    if (!ifs){
        std::cout << " demSettings.txt does not exist!!!" << std::endl;
        std::cout << " demSettings.txt does not exist!!!" << std::endl;
        std::cout << " demSettings.txt does not exist!!!" << std::endl;
        std::cout << " demSettings.txt does not exist!!!" << std::endl;
    }
    
    std::string str;
    double tmpNum;
    int numColumns;

    while (ifs >> str){
        
        if (str == "boundingBox"){
            std::cout << "str was boundingBox" << std::endl;
            ifs.ignore(std::numeric_limits<std::streamsize>::max(),'\n'); //go to the next line
            ifs.ignore(std::numeric_limits<std::streamsize>::max(),'\n'); //ignore {
            ifs.ignore(std::numeric_limits<std::streamsize>::max(),'\n'); //ignore #
            ifs >> boxXmin >> boxXmax  >> boxYmin  >> boxYmax >> boxZmin >> boxZmax;
            ifs.ignore(std::numeric_limits<std::streamsize>::max(),'\n'); //go to the next line
            ifs.ignore(std::numeric_limits<std::streamsize>::max(),'\n'); //ignore }
        }else if(str == "shuffle"){
            std::cout << "str was shuffle" << std::endl;
            ifs.ignore(std::numeric_limits<std::streamsize>::max(),'\n'); //go to the next line
            ifs.ignore(std::numeric_limits<std::streamsize>::max(),'\n'); //ignore {
            ifs.ignore(std::numeric_limits<std::streamsize>::max(),'\n'); //ignore #
            ifs >> shuffleXmin >> shuffleXmax  >> shuffleYmin  >> shuffleYmax >> shuffleZmin >> shuffleZmax;
            ifs.ignore(std::numeric_limits<std::streamsize>::max(),'\n'); //go to the next line
            ifs.ignore(std::numeric_limits<std::streamsize>::max(),'\n'); //ignore }
        }else if (str == "others"){
            ifs.ignore(std::numeric_limits<std::streamsize>::max(),'\n'); //go to the next line
            ifs.ignore(std::numeric_limits<std::streamsize>::max(),'\n'); //ignore {
            ifs.ignore(std::numeric_limits<std::streamsize>::max(),'\n'); //ignore #
            ifs >> numPart >> seed >> dt >> startTime >> endTime >> outputTiming >> wallAmp >> refreshFreq;
            std::cout << "numparticle is " << numPart << std::endl;
            ifs.ignore(std::numeric_limits<std::streamsize>::max(),'\n'); //go to the next line
            ifs.ignore(std::numeric_limits<std::streamsize>::max(),'\n'); //ignore }
        }else if (str == "particleTypes"){
            numColumns = 5;
            ifs >> numParticleTypes;
            std::cout << "numparticleTypes is " << numParticleTypes << std::endl;
            particles = std::vector<double>(numWalls*numColumns);
            ifs.ignore(std::numeric_limits<std::streamsize>::max(),'\n'); //go to the next line
            ifs.ignore(std::numeric_limits<std::streamsize>::max(),'\n'); //ignore {
            ifs.ignore(std::numeric_limits<std::streamsize>::max(),'\n'); //ignore #

            for (int i=0 ; i<numParticleTypes; i++){
                for (int j=0; j<numColumns; j++){
                    ifs >> particles[index(i,j,numColumns)];
                }
            }

            ifs.ignore(std::numeric_limits<std::streamsize>::max(),'\n'); //go to the next line
            ifs.ignore(std::numeric_limits<std::streamsize>::max(),'\n'); //ignore }

        }
    }
    
    

    ifs.close();
    std::cout << "Reading done"<< std::endl;
}

