#include "sphere.hpp"
#include "iostream"
float Sphere(float* params, 
            Eigen::VectorXf & position, 
            Eigen::VectorXf & exitpoint,
            Eigen::VectorXf & normal){
    
    //Give the distance to a sphere of radious parameters[0] centered in 
    //(parameters[1],...,parameters[N])
    float r = 0.0f; 

    for(int i = 0; i < position.size(); i++){
        
        r += pow(position(i)-params[i+1],2.0);

    }

    r = sqrt(r);
    normal = position/r;
    for(int i = 0; i < position.size(); i++){
        exitpoint(i) = normal(i) * params[0];
    }
    
    if(r < params[0]){
        return r - params[0];
    }

    return r - params[0];
    
}
