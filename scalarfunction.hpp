#ifndef SCALARFUNCTION
#define SCALARFUNCTION

#include <iostream>
#include <string>
#include <eigen3/Eigen/Core>
#include "lookuptable.hpp"

/*This class takes care of the scalar functions in BVP*/
class ScalarFunction{

    typedef float (*pfscalar)(Eigen::VectorXf , float );

    private:

        //stores the function if it is analytically defined
        pfscalar function;

        //stores the function if it is stored in a LUT
        LookUpTable lookuptable;

        //True if the function is analytic, false if it is not
        bool analytic;

    public:

        //Initializes the class
        ScalarFunction(void);

        /*
        Initialization of the object with the function
        it is suposed to perform.
        input has to be of the kind:
        float (*input)(Eigen::VectorXf X, Eigen::VectorXf N);
        */
        void Init(pfscalar input);

        /*
        Initialization of the object with the lut where
        The values of the function are stored
        dim has to be an unsigned integer
        input has to be a std::string
        */
        void Init(unsigned int dim, 
                  std::string input);

        /*
        Returns the value of the function in X with
        normal vector N.
        Inputs: std::EigenvectorXf, std::EigenvectorXf
        */
        float Value(Eigen::VectorXf position, 
                    float t);
};  
//Default function which always returns  0.0f
float Default_Scalar(Eigen::VectorXf X, float t);
#endif 

