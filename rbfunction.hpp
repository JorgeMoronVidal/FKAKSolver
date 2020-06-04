#ifndef RBP
#define RBP
#include <eigen3/Eigen/Core>
/*This class takes care of the RBP functions of the meshless interpolator*/
class RBFunction{

	typedef float (*pRBF)(Eigen::VectorXf , Eigen::VectorXf, float c2);

    private:

        //stores the function if it is analytically defined
        pRBF function;
        
    public:

        //Initializes the class
        RBFunction(void);

        /*
        Initialization of the object with the function
        it is suposed to perform.
        input has to be of the kind:
        float (*input)(Eigen::VectorXf X, Eigen::VectorXf N);
        */
        void Init(pRBF input);

        /*
        Returns the value of the function in X with
        normal vector N.
        Inputs: std::EigenvectorXf, std::EigenvectorXf
        */
        float Value(Eigen::VectorXf position, 
                    Eigen::VectorXf normal,
                    float c2);
};  

//Default function which always returns  0.0f
float Default_RBF(Eigen::VectorXf x, 
                  Eigen::VectorXf x_i,
                  float c2);
#endif