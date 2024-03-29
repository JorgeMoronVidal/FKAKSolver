#include "rectangle.hpp"

float Rectangle2D(float* params, 
            Eigen::VectorXf & position, 
            Eigen::VectorXf & exitpoint,
            Eigen::VectorXf & normal){

    float center[2] = {0.5f* (params[0]+ params[2]), 0.5f* (params[1]+ params[3])};
    float halfside[2] = {params[2]-center[0], params[3] - center[1]};
    float distances[2] = {fabs(position(0) - center[0]), fabs(position(1) - center[1])};
    float r;
    exitpoint = position;
    if(position(0) > center[0]){

        if(position(1) > center[1]){
            //First quadrant
            if((distances[0] <= halfside[0]) && (distances[1] <= halfside[1])){
                //Is inside
                if(distances[0] -halfside[0] > distances[1] - halfside[1]){
                    //It's closer to the right side
                    exitpoint(0) = halfside[0] + center[0];
                    normal(0) = 1.0f;
                    normal(1) = 0.0f;
                    return distances[0] - halfside[0];

                } else {
                    //It's closer to the top side
                    exitpoint(1) = halfside[1] + center[1];
                    normal(0) = 0.0f;
                    normal(1) = 1.0f;
                    return distances[1] - halfside[1];

                }
            }  else {
                //Is outside

                if(distances[0] > halfside[0]){

                    if(distances[1] > halfside[1]){
                        //It's closer to the first quadrant corner
                        
                        r = sqrt((distances[0]-halfside[0])*(distances[0]-halfside[0])+ 
                        (distances[1]-halfside[1])*(distances[1]-halfside[1]));
                        exitpoint(0) = halfside[0] + center[0];
                        exitpoint(1) = halfside[1] + center[1];
                        normal(0) = 0.70710678118f;
                        normal(1) = 0.70710678118f;
                        return r;

                    } else {
                        //It's closer to the right side
                        exitpoint(0) = halfside[0] + center[0];
                        normal(0) = 1.0f;
                        normal(1) = 0.0f;
                        return  distances[0] - halfside[0];

                    }

                } else {
                        //It's closer to the bottom side
                        exitpoint(1) = halfside[1] + center[1];
                        normal(0) = 0.0f;
                        normal(1) = 1.0f;
                        return distances[1] - halfside[1];
                    }
                }

            } else {
                //Fourth quadrant
               if((distances[0] <= halfside[0]) && (distances[1] <= halfside[1])){
                //Is inside
                if((distances[0] -halfside[0] > distances[1] - halfside[1])){
                    //It's closer to the right side
                    exitpoint(0) = halfside[0] + center[0];
                    normal(0) = 1.0f;
                    normal(1) = 0.0f;
                    return distances[0] - halfside[0];

                } else {
                    //It's closer to the bottom side
                    exitpoint(1) = - halfside[1] + center[1];
                    normal(0) = 0.0f;
                    normal(1) = -1.0f;
                    return distances[1] - halfside[1];

                }
            }  else {
                //Is outside
                if(distances[0] > halfside[0]){

                    if(distances[1] > halfside[1]){
                        //It's closer to the first quadrant corner
                        r = sqrt((distances[0]-halfside[0])*(distances[0]-halfside[0])+ 
                        (distances[1]-halfside[1])*(distances[1]-halfside[1]));
                        exitpoint(0) = halfside[0] + center[0];
                        exitpoint(1) = -halfside[1] + center[1];
                        normal(0) = 0.70710678118f;
                        normal(1) = -0.70710678118f;
                        return r;

                    } else {
                        //It's closer to the right side
                        exitpoint(0) = halfside[0] + center[0];
                        normal(0) = 1.0f;
                        normal(1) = 0.0f;
                        return  distances[0] - halfside[0];

                    }

                } else {
                        //It's closer to the bottom side
                        exitpoint(1) = -halfside[1] + center[1];
                        normal(0) = 0.0f;
                        normal(1) = -1.0f;
                        return distances[1] - halfside[1];
                    }
                }
            }
        } else {
            if(position(1) > center[1]){

            //Second quadrant
                if((distances[0] <= halfside[0]) && (distances[1] <= halfside[1])){
                //Is inside
                    if((distances[0] -halfside[0] > distances[1] - halfside[1])){
                        //It's closer to the left side
                        exitpoint(0) = -halfside[0] + center[0];
                        normal(0) = -1.0f;
                        normal(1) = 0.0f;
                        return distances[0] - halfside[0];

                    } else {
                        //It's closer to the top side
                        exitpoint(1) = halfside[1] + center[1];
                        normal(0) = 0.0f;
                        normal(1) = 1.0f;
                        return distances[1] - halfside[1];

                }
            }  else {
                //Is outside
                if(distances[0] > halfside[0]){

                    if(distances[1] > halfside[1]){
                        //It's closer to the second quadrant corner
                        r = sqrt((distances[0]-halfside[0])*(distances[0]-halfside[0])+ 
                        (distances[1]-halfside[1])*(distances[1]-halfside[1]));
                        exitpoint(0) = -halfside[0] + center[0];
                        exitpoint(1) = halfside[1] + center[1];
                        normal(0) = -0.70710678118f;
                        normal(1) = 0.70710678118f;
                        return r;

                    } else {
                        //It's closer to the left side
                        exitpoint(0) = -halfside[0] + center[0];
                        normal(0) = -1.0f;
                        normal(1) = 0.0f;
                        return  distances[0] - halfside[0];

                    }

                } else {
                        //It's closer to the top side
                        exitpoint(1) = halfside[1] + center[1];
                        normal(0) = 0.0f;
                        normal(1) = 1.0f;
                        return distances[1] - halfside[1];
                    }
                }
        } else {
            //Third quadrant
                if((distances[0] <= halfside[0]) && (distances[1] <= halfside[1])){

                    //Is inside
                    if(distances[0] -halfside[0] > distances[1] - halfside[1]){
                        //It's closer to the left side
                        exitpoint(0) = -halfside[0] + center[0];
                        normal(0) = -1.0f;
                        normal(1) = 0.0f;
                        return distances[0] - halfside[0];

                    } else {
                    //It's closer to the bottom side
                        exitpoint(1) = -halfside[1] + center[1];
                        normal(0) = 0.0f;
                        normal(1) = -1.0f;
                        return distances[1] - halfside[1];

                    }
            }  else {

                //Is outside
                if(distances[0] > halfside[0]){

                    if(distances[1] > halfside[1]){
                        //It's closer to the first quadrant corner
                        r = sqrt((distances[0]-halfside[0])*(distances[0]-halfside[0])+ 
                        (distances[1]-halfside[1])*(distances[1]-halfside[1]));
                        exitpoint(0) = -halfside[0] + center[0];
                        exitpoint(1) = -halfside[1] + center[1];
                        normal(0) = -0.70710678118f;
                        normal(1) = -0.70710678118f;
                        return r;

                    } else {
                        //It's closer to the right side

                        exitpoint(0) = -halfside[0] + center[0];
                        normal(0) = -1.0f;
                        normal(1) = 0.0f;
                        return  distances[0] - halfside[0];

                    }

                } else {
                        //It's closer to the bottom side
                        exitpoint(1) = -halfside[1] + center[1];
                        normal(0) = 0.0f;
                        normal(1) = -1.0f;
                        return distances[1] - halfside[1];
                    }
                }
            }
        }
    }