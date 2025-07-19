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
    std::ifstream ifs("demSettings.json");

    
    
    if (!ifs){
        std::cout << " demSettings.json does not exist!!!" << std::endl;
        std::cout << " demSettings.json does not exist!!!" << std::endl;
        std::cout << " demSettings.json does not exist!!!" << std::endl;
        std::cout << " demSettings.json does not exist!!!" << std::endl;
        std::abort();
    }
    
    picojson::value settings;
    std::string str;
    ifs >> settings;
    ifs.close(); 
    
    double tmpNum;
    int numColumns;

    //gather all the objects
    picojson::object settingObjects = settings.get<picojson::object>();

    picojson::object boundingBox = settingObjects["boundingBox"].get<picojson::object>();


    boxXmin = boundingBox["Xmin"].get<double>();
    boxYmin = boundingBox["Ymin"].get<double>();
    boxZmin = boundingBox["Zmin"].get<double>();
    boxXmax = boundingBox["Xmax"].get<double>();
    boxYmax = boundingBox["Ymax"].get<double>();
    boxZmax = boundingBox["Zmax"].get<double>();

    picojson::object others = settingObjects["others"].get<picojson::object>();


    numPart = others["numPart"].get<double>(); 
    seed = static_cast<int>(others["seed"].get<double>()); 
    dt = others["dt"].get<double>(); 
    startTime = others["startTime"].get<double>(); 
    endTime = others["endTime"].get<double>(); 
    outputTiming = others["outputTiming"].get<double>(); 
    wallAmp = others["wallAmp"].get<double>(); 
    refreshFreq = others["refreshFreq"].get<double>(); 

    picojson::object shuffle = settingObjects["shuffle"].get<picojson::object>();
    shuffleXmin = shuffle["Xmin"].get<double>();
    shuffleYmin = shuffle["Ymin"].get<double>();
    shuffleZmin = shuffle["Zmin"].get<double>();
    shuffleXmax = shuffle["Xmax"].get<double>();
    shuffleYmax = shuffle["Ymax"].get<double>();
    shuffleZmax = shuffle["Zmax"].get<double>();

    std::cout << boost::format("numparticle is %d")%numPart << std::endl;

   particleTypes = settingObjects["particleTypes"].get<picojson::object>();

    ifs.close();
    std::cout << "Reading done"<< std::endl;
}

