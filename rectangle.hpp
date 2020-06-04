#ifndef RECTANGLE
#define RECTANGLE

#include <math.h>
#include <eigen3/Eigen/Core>
#include <vector>
#include <iostream>
float Rectangle2D(float* params, 
            Eigen::VectorXf & position, 
            Eigen::VectorXf & exitpoint,
            Eigen::VectorXf & normal);

#endif