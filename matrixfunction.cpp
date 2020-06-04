#include"matrixfunction.hpp"
MatrixFunction::MatrixFunction(void){

    function = Default_Matrix;
    analytic = true;

}

void MatrixFunction::Init(pfmatrix input){

    function = input;
    analytic = true;

}

void MatrixFunction::Init(int dim,  std::string input){
    std::string aux;
    lookuptable.resize(dim);
    for(int i = 0; i < dim; i++){

        lookuptable[i].resize(dim);

        for(int j = 0; j < dim; j++){

            aux = input + "/" + std::to_string(i) +std::to_string(j) +"_component";
            lookuptable[i][j].Init(dim, aux);

        }

    }

    analytic = false;
}

Eigen::MatrixXf MatrixFunction::Value(Eigen::VectorXf position, 
                                      float t){
    if(analytic){
        return function(position, t);
    }

    Eigen::MatrixXf out;

    out.resize(position.size(), position.size());

    for(int i = 0; i < position.size(); i++ ){

        for(int j = 0; j < position.size(); j++){

            out(i,j) = lookuptable[i][j].Eval(position);

        }
    }
    return out;
}

Eigen::MatrixXf Default_Matrix(Eigen::VectorXf position, float t){
    
    return Eigen::MatrixXf::Zero(position.size(), position.size());
}
