#ifndef SCALARFUNCTIONN
#define SCALARFUNCTIONN

#include <iostream>
#include <string>
#include <eigen3/Eigen/Core>
#include "lookuptable.hpp"

/*This class takes care of the scalar functions in BVP*/
class ScalarFunctionN{

    typedef float (*pfscalarN)(Eigen::VectorXf , Eigen::VectorXf, float );

    private:

        //stores the function if it is analytically defined
        pfscalarN function;

        //stores the function if it is stored in a LUT
        LookUpTable lookuptable;

        //True if the function is analytic, false if it is not
        bool analytic;

    public:

        //Initializes the class
        ScalarFunctionN(void);

        /*
        Initialization of the object with the function
        it is suposed to perform.
        input has to be of the kind:
        float (*input)(Eigen::VectorXf X, Eigen::VectorXf N);
        */
        void Init(pfscalarN input);

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
                    Eigen::VectorXf normal,
                    float t);
};  
//Default function which always returns  0.0f
float Default_ScalarN(Eigen::VectorXf X, Eigen::VectorXf N, float t);
#endif 

