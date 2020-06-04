#include "cube.hpp"

float Cube(float* params, 
            Eigen::VectorXf & position, 
            Eigen::VectorXf & exitpoint,
            Eigen::VectorXf & normal){
    
    float dist = 0.0f, distance[6];
    int index = 0;
    bool out = false;
    exitpoint = position;
    normal.Zero(3);
    for(int i = 0; i < 3; i++){
        distance[i] = (position(i) - params[i+1]) - params[0];
        if(distance[i] >= 0.0f){
            dist += distance[i] * distance[i];
            exitpoint(i) = params[1+i] + params[0];
            normal(i) = 1.0f;
            out = true;
        }
    }
    for(int i = 0; i < 3; i++){
        distance[3+i] =(params[1+i] - position(i)) - params[0];
        if(distance[3+i] >= 0.0f){
            dist += distance[3+i] * distance[3+i];
            exitpoint(i) = params[1+i] - params[0];
            normal(i) = -1.0f;
            out = true;
        }
    }

    if(! out){
        dist = -params[0] -0.1;
        for(int i = 0; i < 6; i++){
            if(distance[i] > dist){
                dist = distance[i];
                index = i;
            }
        }
        if(index < 3){
            normal(index) = 1.0f;
            exitpoint(index) = params[1+index] + params[0];
        } else {
            normal(index-3) = -1.0f;
            exitpoint(index-3) = params[1+index-3] - params[0];
        }
        return dist;
    } else {
        normal = normal/normal.norm();
        return sqrt(dist);
    }

}