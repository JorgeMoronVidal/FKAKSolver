#include "rbfunction.hpp"

RBFunction::RBFunction(void){
    function = Default_RBF;
}

void RBFunction::Init(pRBF input){
    function = input;
}


float RBFunction::Value (Eigen::VectorXf x, 
                          Eigen::VectorXf x_i,
                          float c2){
    return function(x, x_i, c2);
}

float Default_RBF(Eigen::VectorXf x, 
                  Eigen::VectorXf x_i,
                  float c2){
    return 0.0f;
}