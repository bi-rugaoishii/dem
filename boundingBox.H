#ifndef _boundingBox_H_
#define _boundingBox_H_

class boundingBox{
    public:
        boundingBox();
        boundingBox(double minX, double maxX, double minY, double maxY, double minZ, double maxZ);

        inline const double &minX() const{
            return minX_;
        }

        inline const double &minY() const{
            return minY_;
        }

        inline const double &minZ() const{
            return minZ_;
        }

        inline const double &maxX() const{
            return maxX_;
        }

        inline const double &maxY() const{
            return maxY_;
        }

        inline const double &maxZ() const{
            return maxZ_;
        }

        inline const int &sizeSplitX() const{
            return sizeSplitX_;
        }

        inline const int &sizeSplitY() const{
            return sizeSplitY_;
        }

        inline const int &sizeSplitZ() const{
            return sizeSplitZ_;
        }

        inline const double &dx() const{
            return dx_;
        }
        inline const double &dy() const{
            return dy_;
        }

        inline const double &dz() const{
            return dz_;
        }

        void showRange();

        void showSplit();

        void split(double maxRad);





    private:
        double minX_,maxX_,minY_,maxY_,minZ_,maxZ_;
        double dx_,dy_,dz_;
        double rangeX_,rangeY_,rangeZ_;  
        bool isSplit_;
        int sizeSplitX_,sizeSplitY_,sizeSplitZ_;
        std::vector<int> indBox_; 

};


#endif
