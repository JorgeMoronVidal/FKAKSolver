#include "boundary.hpp"
Boundary::Boundary(void){

    analytic = true;

}

void Boundary::_init_(pfbound fbound, pfstop fstop){

    Distance_Analytic = fbound;
    stop = fstop;
    analytic = true;
}

void Boundary::_init_(int dim, std::string lupbound, pfstop fstop){
    Distance_Numeric.Init(dim, lupbound);
    stop = fstop;
    analytic = false;
}

float Boundary::Dist(float* params, 
            Eigen::VectorXf & position, 
            Eigen::VectorXf & exitpoint,
            Eigen::VectorXf & normal){
    
    if(analytic){

        return Distance_Analytic(params,position,exitpoint,normal);
    }

    return Distance_Numeric.Eval(params, position, exitpoint, normal);
}
