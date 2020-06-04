#include "equation.hpp"
#include <math.h>
#include <iostream>

float Equation_c(Eigen::VectorXf X, float t){
    return 0.0f;
}
float Equation_f(Eigen::VectorXf X, float t){
    return 1.0f;
}
float Equation_g(Eigen::VectorXf X, float t){
    return 0.0f;
}
Eigen::VectorXf Equation_b(Eigen::VectorXf X, float t){
    return X*0.0f;
}
Eigen::VectorXf Equation_F(Eigen::VectorXf X, float t){
    Eigen::VectorXf F(2);
    float aux = 0.0f;
    for(int i = 0.0f; i < 10.0f; i ++){
        for(int j = 0.0f; j < 10.0f; j++){
            aux += -pow(-1.0f,(float)i+j)*0.5*M_PI*sin((2.0f*i+1.0f)*0.5*M_PI*X(0))*cos((2.0f*j+1.0f)*0.5*M_PI*X(1))*
            (1.0f/((2.0f*j+1.0f)*((2.0f*i+1)*(2.0f*i+1)+(2.0f*j+1)*(2.0f*j+1))));
        }
    }

    F(0) = aux;
    aux = 0.0f;

    for(int i = 0.0f; i < 10.0f; i ++){
        for(int j = 0.0f; j < 10.0f; j++){
            aux += -pow(-1.0f,(float)i+j)*cos((2.0f*i+1.0f)*0.5*M_PI*X(0))*0.5*M_PI*sin((2.0f*j+1.0f)*0.5*M_PI*X(1))*
            (1.0f/((2.0f*i+1.0f)*((2.0f*i+1)*(2.0f*i+1)+(2.0f*j+1)*(2.0f*j+1))));
        }
    }

    F(1) = aux;

    return -F *128.0f/pow(M_PI,4.0);
}
Eigen::MatrixXf Equation_sigma(Eigen::VectorXf X, float t){
    Eigen::MatrixXf sigma(2,2);
    sigma << 1.0f ,0.0f,
            0.0f, 1.0f;
    return sigma;
}
float Equation_u(Eigen::VectorXf X, float t){
    float aux = 0.0f;
    for(int i = 0.0f; i < 50.0f; i ++){
        for(int j = 0.0f; j < 50.0f; j++){
            aux += pow(-1.0f,(float)i+j)*cos((2.0f*i+1.0f)*0.5*M_PI*X(0))*cos((2.0f*j+1.0f)*0.5*M_PI*X(1))*
            (1.0f/((2.0f*i+1.0f)*(2.0f*j+1.0f)*((2.0f*i+1)*(2.0f*i+1)+(2.0f*j+1)*(2.0f*j+1))));
        }
    }

    return aux*128.0f/pow((double)M_PI,4.0);
}

float Equation_RBF(Eigen::VectorXf x , Eigen::VectorXf xj, float c2){
    float r2 = pow(x(0)-xj(0),2) + pow(x(1)-xj(1),2);
    return sqrt(r2 + c2);
}

float Equation_d2dx2RBF(Eigen::VectorXf x , Eigen::VectorXf xj, float c2){
    float r2 = pow(x(0)-xj(0),2) + pow(x(1)-xj(1),2);
    return (pow(x(1)-xj(1),2)+c2)/pow(r2+c2, 1.5);
}

float Equation_d2dy2RBF(Eigen::VectorXf x , Eigen::VectorXf xj, float c2){
    float r2 = pow(x(0)-xj(0),2) + pow(x(1)-xj(1),2);
    return (pow(x(0)-xj(0),2)+c2)/pow(r2+c2, 1.5);
}

float Equation_uRBF(Eigen::VectorXf x , Eigen::VectorXf xj, float c2){
    return Equation_d2dx2RBF(x,xj,c2) + Equation_d2dy2RBF(x,xj,c2) + 2.0;
}

bool Stopping(Eigen::VectorXf position){

    return true;
}