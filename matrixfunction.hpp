#ifndef MATRIXFUNCTION
#define MATRIXFUNCTION

#include <iostream>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>
#include "lookuptable.hpp"

/*This class takes care of the functions which image is a matrix in BVP*/
class MatrixFunction{

    typedef Eigen::MatrixXf (*pfmatrix)(Eigen::VectorXf , float);

    private:

        //stores the function if it is analytically defined
        pfmatrix function;

        //stores the function if it is stored in a LUT
        std::vector<std::vector<LookUpTable>> lookuptable;

        //True if the function is analytic, false if it is not
        bool analytic;

    public:
    
        //Initializes the class
        MatrixFunction(void);

        /*
        Initialization of the object with the function
        it is suposed to perform.
        input has to be of the kind:
        Eigen::MatrixXf (*input)(Eigen::VectorXf X, float t);
        */
        void Init(pfmatrix input);

        /*
        Initialization of the object with the lut where
        The values of the function are stored
        dim has to be an unsigned integer
        input has to be a std::string
        */
        void Init(int dim, 
                  std::string input);

        /*
        Returns the value of the function in X with
        normal vector N.
        Inputs: std::EigenvectorXf, float
        */
        Eigen::MatrixXf Value(Eigen::VectorXf position, 
                    float t);
};  

//Default function which always returns a matrix full of 0.0f
Eigen::MatrixXf Default_Matrix(Eigen::VectorXf position, float t);

#endif 