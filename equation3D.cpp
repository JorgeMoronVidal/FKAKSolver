#include "edepaco.hpp"
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
float Equation_p(Eigen::VectorXf X, float t){
    return 0.0f;
}
Eigen::VectorXf Equation_b(Eigen::VectorXf X, float t){
    return X*0.0f;
}
Eigen::VectorXf Equation_F(Eigen::VectorXf X, float t){
    Eigen::VectorXf F(3);
    float aijk;
    for(int i = 1; i < 15; i+= 2){
        for(int j = 1; j < 15; j += 2){
            for(int k = 1; k < 15; k += 2){

                aijk = (pow(-1.0,(i+j+k+1)*0.5f)/(i*j*k*(i*i+j*j+k*k)))*(1.0f-exp(-M_PI*M_PI*(i*i+j*j+k*k)*(t)/8.0f));
                F[0] += -aijk*i*sin(i*M_PI*X(0)*0.5)*cos(j*M_PI*X(1)*0.5)*cos(k*M_PI*X(2)*0.5);
                F[1] += -aijk*j*cos(i*M_PI*X(0)*0.5)*sin(j*M_PI*X(1)*0.5)*cos(k*M_PI*X(2)*0.5);
                F[2] += -aijk*k*cos(i*M_PI*X(0)*0.5)*cos(j*M_PI*X(1)*0.5)*sin(k*M_PI*X(2)*0.5);

            }
        }
    }
    return -F*256.0f/pow(M_PI,4.0);
}
Eigen::MatrixXf Equation_sigma(Eigen::VectorXf X, float t){
    return Eigen::Matrix<float,3,3>::Identity();
}
float Equation_u(Eigen::VectorXf X, float t){
    float aux = 0.0f,aijk;
    for(int i = 1; i < 50; i+= 2){
        for(int j = 1; j < 50; j += 2){
            for(int k = 1; k < 50; k += 2){

                aijk = (pow(-1.0,(i+j+k+1)*0.5f)/(i*j*k*(i*i+j*j+k*k)))*(1.0f-exp(-M_PI*M_PI*(i*i+j*j+k*k)*(t)/8.0f));
                aux += aijk*cos(i*M_PI*X(0)*0.5)*cos(j*M_PI*X(1)*0.5)*cos(k*M_PI*X(2)*0.5);

            }
        }
    }

    return aux*512.0f/pow(M_PI,5.0);
}
float Equation_Psi(Eigen::VectorXf X, Eigen::VectorXf normal, float t){
    return 0.0f;
}
float Equation_Varphi(Eigen::VectorXf X, Eigen::VectorXf normal, float t){
    return 0.0f;
}
bool Stopping(Eigen::VectorXf position){

    return true;
}