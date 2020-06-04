#include "equation.hpp"
#include <math.h>
#include <iostream>

float Equation_c(Eigen::VectorXf X, float t){
    return 0;
}
float Equation_f(Eigen::VectorXf X, float t){
    return 1.0f;
}
float Equation_g(Eigen::VectorXf X, float t){
    float aux = 0.0f;
    for(int i = 0; i < X.size(); i++){
        aux += X(i);
    }
    return aux;
}
Eigen::VectorXf Equation_b(Eigen::VectorXf X, float t){
    return 0.0f*X;
}
Eigen::VectorXf Equation_F(Eigen::VectorXf X, float t){
    Eigen::VectorXf aux = X;
    for (int i = 0 ; i < X.size(); i++){
        aux(i) = (-2.0f/3)*X(i)+1;
    }
    return -aux;
}
Eigen::MatrixXf Equation_sigma(Eigen::VectorXf X, float t){
    Eigen::MatrixXf sigma(X.size(),X.size());
    for (int i = 0; i < X.size(); i++){
        sigma(i,i) = 1.0f;
    }
    return sigma;
}
float Equation_u(Eigen::VectorXf X, float t){
    return ((1.0f/3)*(1-pow(X.norm(),2.0)))+Equation_g(X,t);
}
bool Stopping(Eigen::VectorXf position){

    return true;
}