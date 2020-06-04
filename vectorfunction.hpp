#ifndef VECTORFUNCTION
#define VECTORFUNCTION

#include <iostream>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>
#include "lookuptable.hpp"

/*This class takes care of the functions which image is a vector in BVP*/
class VectorFunction{

    typedef Eigen::VectorXf (*pfvector)(Eigen::VectorXf,float);

    private:

        //stores the function if it is analytically defined
        pfvector function;

        //stores the function if it is stored in a LUT
        std::vector<LookUpTable> lookuptable;

        //True if the function is analytic, false if it is not
        bool analytic;

    public:

        //Initializes the class
        VectorFunction(void);

        /*
        Initialization of the object with the function
        it is suposed to perform.
        input has to be of the kind:
        Eigen::VectorXf (*input)(Eigen::VectorXf X, Eigen::VectorXf N);
        */
        void Init(pfvector input);

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
        Inputs: std::EigenvectorXf, std::EigenvectorXf
        */
        Eigen::VectorXf Value(Eigen::VectorXf position,float t);
};  

//Default function which always returns a vector full of 0.0f
Eigen::VectorXf Default_Vector(Eigen::VectorXf position, float t);

#endif 