#include <iostream>
#include <stdio.h> 
#include <eigen3/Eigen/Core>
#include "cube.hpp"

int main(void){
    float params[4] = {1.0f, 0.0f, 0.0f, 0.0f}, dist;
    Eigen::VectorXf X0, N, E_P;
    X0.resize(3); N.resize(3); E_P.resize(3);
    X0[0] = 1.5f; X0[1] = -1.5; X0[2] = -1.5f;
    dist = Cube(params, X0, E_P, N);
    printf(" Point %f %f %f\n Dist %f\n ExitPoint %f %f %f\n Normal %f %f %f\n", X0(0), X0(1), X0(2), dist, E_P(0), E_P(1), E_P(2), N(0), N(1), N(2));
    return 0;
}