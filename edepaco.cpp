#include "edepaco.hpp"

float Equation_c(Eigen::VectorXf X, float t){
    return -pow(X.norm(),2)* exp(-X(0));
}
float Equation_f(Eigen::VectorXf X, float t){
    return 2.0f * sin(5.0f*t + 7.0f*X(0) + 9.0f*X(1)) + 130.0f*
    cos(5.0f*t + 7.0f*X(0) + 9.0f*X(1)) -Equation_c(X,t) * Equation_u(X,t);
}
float Equation_p(Eigen::VectorXf X, float t){
    return 2.1f + cos(7.0f * X(0) + 9.0f * X(1));
}
float Equation_g(Eigen::VectorXf X, float t){
    return 2.1f + cos(5.0f * t + 7.0f * X(0) + 9.0f * X(1));
}
float Equation_Psi(Eigen::VectorXf X, Eigen::VectorXf normal, float t){
    return -(normal(0) * 7.0f + normal(1) * 9.0f) * sin(5.0f * t + 7.0f * X(0) + 9.0f * X(1)) 
            + Equation_u(X,t);
}
float Equation_Varphi(Eigen::VectorXf X, Eigen::VectorXf normal, float t){
    return -1.0f;
}
bool Stopping(Eigen::VectorXf X){
    if(X(1) >0.0f) return true;
    return false;
}
Eigen::VectorXf Equation_b(Eigen::VectorXf X, float t){
    Eigen::VectorXf b(X.size());
    b << 1.0f, 0.0f;
    return b;
}
Eigen::VectorXf Equation_F(Eigen::VectorXf X, float t){
    Eigen::VectorXf F(X.size());
    F << 7.0f, 9.0f;
    F = F*sin(5.0 *t + 7.0f * X(0) + 9.0f * X(1) );
    return F;
}
Eigen::MatrixXf Equation_sigma(Eigen::VectorXf X, float t){
    Eigen::MatrixXf sigma(X.size(), X.size());
    sigma << 1.41421356237f, 0.0f,
             0.0f, 1.41421356237f;
    return sigma;
}
float Equation_u(Eigen::VectorXf X, float t){
    return 2.1f + cos(5.0f * t + 7.0f * X(0) +
    9.0f * X(1));
}
