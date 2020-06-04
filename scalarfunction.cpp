#include"scalarfunction.hpp"

ScalarFunction::ScalarFunction(void){
    function = Default_Scalar;
    analytic = true;
}

void ScalarFunction::Init(pfscalar input){
    function = input;
    analytic = true;
}

void ScalarFunction::Init(unsigned int dim, 
                          std::string input){
    lookuptable.Init(dim, input);
    analytic = false;
}

float ScalarFunction::Value(Eigen::VectorXf position, 
                            float t){
    if(analytic){
        return function(position, t);
    }

    return lookuptable.Eval(position);
}

float Default_Scalar(Eigen::VectorXf X, float t){
    return 0.0f;
}
