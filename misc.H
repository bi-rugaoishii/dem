#ifndef _misc_H_
#define _misc_H_

#include <Eigen/Dense>

#include <cmath>

inline int index(const int &z, const int &y, const int &x, const int &numx, const int &numy){
    return z*numy*numx+y*numx+x;
}

inline int index(const int &y, const int &x, const int &numx){
    return y*numx+x;
}



#endif
