#include <eigen3/Eigen/Core>

float Equation_c(Eigen::VectorXf X, float t);
float Equation_f(Eigen::VectorXf X, float t);
float Equation_g(Eigen::VectorXf X, float t);
Eigen::VectorXf Equation_b(Eigen::VectorXf X, float t);
Eigen::VectorXf Equation_F(Eigen::VectorXf X, float t);
Eigen::MatrixXf Equation_sigma(Eigen::VectorXf X, float t);
float Equation_u(Eigen::VectorXf X, float t);
float Equation_RBF(Eigen::VectorXf x , Eigen::VectorXf xj, float c2);
float Equation_d2dx2RBF(Eigen::VectorXf x , Eigen::VectorXf xj, float c2);
float Equation_d2dy2RBF(Eigen::VectorXf x , Eigen::VectorXf xj, float c2);
float Equation_uRBF(Eigen::VectorXf x , Eigen::VectorXf xj, float c2);
bool Stopping(Eigen::VectorXf position);