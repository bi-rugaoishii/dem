#ifndef _vector3d_H_
#define _vector3d_H_

#include <vector>

class vector3d{
    public:
        vector3d(){
            this->v = std::vector<double>(3);
            this->v[0] = 0.;
            this->v[1] = 0.;
            this->v[2] = 0.;
        }
        vector3d(double a){
            this->v = std::vector<double>(3);
            this->v[0] = a;
            this->v[1] = a;
            this->v[2] = a;
        }

        vector3d(double a,double b, double c){
            this->v = std::vector<double>(3);
            this->v[0] = a;
            this->v[1] = b;
            this->v[2] = c;
        }

        std::vector<double> v;

        inline double& operator[](int i){
            return v[i];
        }

        inline const double& operator[](int i) const{
            return v[i];
        }

};

//dot product
inline double operator&(const vector3d& a, const vector3d& b){
    double result(a[0]*b[0]+a[1]*b[1]+a[2]*b[2]);
    return result;
}


//scalar product
inline vector3d operator*(const double& scalar, const vector3d& b){
    vector3d result(scalar*b[0],scalar*b[1],scalar*b[2]);
    return result;
}

inline vector3d operator*(const vector3d& b, const double& scalar){
    vector3d result(scalar*b[0],scalar*b[1],scalar*b[2]);
    return result;
}

//scalar division

inline vector3d operator/(const vector3d& b, const double& scalar){
    vector3d result(b[0]/scalar,b[1]/scalar,b[2]/scalar);
    return result;
}

//cross product
inline vector3d operator*(const vector3d& a, const vector3d& b){
    vector3d result(a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]);
    return result;
}

inline vector3d operator+(const vector3d& a, const vector3d& b){
    vector3d result(a[0]+b[0],a[1]+b[1],a[2]+b[2]);
    return result;
}

inline vector3d operator-(const vector3d& a, const vector3d& b){
    vector3d result(a[0]-b[0],a[1]-b[1],a[2]-b[2]);
    return result;
}


#endif
