#include"scalarfunctionN.hpp"

ScalarFunctionN::ScalarFunctionN(void){
    function = Default_ScalarN;
    analytic = true;
}

void ScalarFunctionN::Init(pfscalarN input){
    function = input;
    analytic = true;
}

void ScalarFunctionN::Init(unsigned int dim, 
                          std::string input){
    lookuptable.Init(dim, input);
    analytic = false;
}

float ScalarFunctionN::Value(Eigen::VectorXf position, 
                            Eigen::VectorXf normal, 
                            float t){
    if(analytic){
        return function(position, normal, t);
    }

    return lookuptable.Eval(position);
}

float Default_ScalarN(Eigen::VectorXf X, Eigen::VectorXf N, float t){
    return 0.0f;
}
