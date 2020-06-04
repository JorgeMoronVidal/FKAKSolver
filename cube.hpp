#ifndef CUBE
#define CUBE

#include <math.h>
#include <eigen3/Eigen/Core>
#include <vector>
#include <iostream>
float Cube(float* params, 
            Eigen::VectorXf & position, 
            Eigen::VectorXf & exitpoint,
            Eigen::VectorXf & normal);

#endif