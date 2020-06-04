#ifndef BoundaryValueProblem
#define BoundaryValueProblem

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>
#include "lookuptable.hpp"
#include "scalarfunction.hpp"
#include "scalarfunctionN.hpp"
#include "vectorfunction.hpp"
#include "matrixfunction.hpp"
#include "rbfunction.hpp"
#include "boundary.hpp"
typedef float (*pfscalar)(Eigen::VectorXf, float);
 typedef float (*pfscalarN)(Eigen::VectorXf , Eigen::VectorXf, float);
typedef Eigen::VectorXf (*pfvector)(Eigen::VectorXf,float);
typedef Eigen::MatrixXf (*pfmatrix)(Eigen::VectorXf , float);
typedef float (*pfbound)(float* params, 
                             Eigen::VectorXf & position, 
                             Eigen::VectorXf & exitpoint,
                             Eigen::VectorXf & normal);
typedef bool (*pfstop)(Eigen::VectorXf);
typedef float (*pRBF)(Eigen::VectorXf , Eigen::VectorXf, float c);
/*
 BVP stores the functions that can be used to compute the solution of 
 such kind of problems using the Feynman-Kak formula.
*/
class BVP
{  

    private:

        std::map<std::string, bool> control;

    public:

        //Functions whose image is scalar s.t. c(X), f(X), u(X), varphi(X)
        ScalarFunction f,c, u, g, p;
        ScalarFunctionN psi, varphi;
        //Functions whose image is a vector s.t b(X), F(x), g(X), psi(X)
        VectorFunction F, mu, b;

        //sigma is the /sigma matrix
        MatrixFunction sigma;
        
        //The boundary is stored in boundary
        Boundary boundary;

        //RBP function for the meshless method
        RBFunction rbf;
        //All the entries of the maps are set as Default. bvplat remains empty
        BVP(void);

        //The object is properly initialized given the functions or the directorys where the look up tables are stored.
        void BVP_init(int dim,
                    std::map<std::string, pfscalar> map_fscalar,
                    std::map<std::string, pfscalarN> map_fscalarN,
                    std::map<std::string, pfvector> map_fvector,
                    std::map<std::string, pfmatrix> map_fmatrix,
                    std::map<std::string, std::string> map_lut);

        void BVP_init(int dim,
                    std::map<std::string, pfscalar> map_fscalar,
                    std::map<std::string, pfscalarN> map_fscalarN,
                    std::map<std::string, pfvector> map_fvector,
                    std::map<std::string, pfmatrix> map_fmatrix,
                    std::map<std::string, std::string> map_lut,
                    pRBF rbfunc);
        
        /*Initialization of the surface variable with analytic boundary
        -boundary is a function s.t. (*float)(float* params, Eigen::VectorXf & position, 
                             Eigen::VectorXf & exitpoint, Eigen::VectorXf & normal)
        -stopf is a function s.t. (*pfstop)(Eigen::MatrixXf)
        */
        void Boundary_init(pfbound boundary, pfstop stopf);
        
        /*Initialization of the surface variable with LUT boundary
        -dim is the dimension of the problem
        -boundary is the directory where the lup of the distance is storedAh per
        -stopf is a function s.t. (*pfstop)(Eigen::MatrixXf)
        */
        void Boundary_init(int dim, std::string boundary , pfstop stopf);

};




#endif