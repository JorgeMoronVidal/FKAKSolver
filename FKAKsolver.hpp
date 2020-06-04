#ifndef FKAKSOLVER
#define FKAKSOLVER

#include <stdlib.h>
#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "BVP.hpp"
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// https://gcc.gnu.org/onlinedocs/cpp/Variadic-Macros.html
// variadic macro to print to both file and stdout
#define PRINTF2(fp, ...) {printf(__VA_ARGS__);fprintf(fp,__VA_ARGS__);}

class FKAKSolver
{
private:
    /*Boundary value problem's equations'*/
    BVP bvp;
    /*Flux of the class flow
     -0 if class wasn't properly initializated
     -1 if class was properly initializated*/
    int flux;
    unsigned int N_rngcalls, N_trayectories;
    /*Params: Parameters of the surface
      -t:Time (For parabolic equations).
      -h: Time discretization.
      -sqrth: square root of h.
      -Y: FKAK Y.
      -Z: FKAK Z.
      -ji_t auxiliary time in order to deal with reflecint BC's.
      -xi: Control variates variable.*/
    float *params,t, h, sqrth, Y, Z, ji_t, xi;
    /*Stores some quantities in order to compute variances and averages*/
    double sums[7]; 
    /*-N: Normal vector to the boundary.
      -E_P: Exit point.
      -X: Position.
      -Variance reduction mu function.
      -increment: Random increment in position.*/
    Eigen::VectorXf N, E_P, X, mu, increment;
    /*Sigma matrix*/
    Eigen::MatrixXf sigma;
    /*RNG type*/
    const gsl_rng_type * T;
    /*RNG*/
    gsl_rng * rng;
    /*Outcome stores the status of the particle*/
    enum outcome {in, stop, reflect, time_out};
    outcome status;
    /*Updates increment with random numbers*/
    void increment_update(void);
    /*Reset of FKAK */
    void FKAK_reset(Eigen::VectorXf X0);
    void FKAK_reset(Eigen::VectorXf X0, float T_start);
    /*Updates statistical quantities of Milstein algorithm*/
    void Update_Stat(float sol_0, float xi);
    /*One plain step of the Milstein's algorithm*/
    void Step(void);
    /*One step of the Euler's discretization with Variance Reduction*/
    void VR_Step(void);
    /*One step of the Euler's discretization with Control Variates*/
    void CV_Step(void);
    /*One step of the Euler's discretization with Variance Reduction and 
    Control Variates*/
    void VR_CV_Step(void);
    /*Returns true if the particle is inside the boundary false if it is not.*/
    bool Inside_Elliptic(Eigen::VectorXf & position, 
                         Eigen::VectorXf & normal, 
                         Eigen::VectorXf & exitpoint);
    bool Inside_Elliptic(Eigen::VectorXf & position, 
                         Eigen::VectorXf & normal, 
                         Eigen::VectorXf & exitpoint,
                         Eigen::MatrixXf & sigma);
    bool Inside_Parabolic(Eigen::VectorXf & position, 
                         Eigen::VectorXf & normal, 
                         Eigen::VectorXf & exitpoint,
                         float & time);
    bool Inside_Parabolic(Eigen::VectorXf & position, 
                         Eigen::VectorXf & normal, 
                         Eigen::VectorXf & exitpoint,
                         float & time,
                         float & sqrdt,
                         Eigen::MatrixXf & sigma);
    /*Multilevel related variables and functions*/
    /*-l level
      -option 1 -- sub-sampling, no offset
              2 -- nothing
              3 -- sub-sampling, offset
      -M is 
    */
    uint16_t l, option, M, D, L;
    uint32_t dNl[21];
    Eigen::VectorXf Xc, incrementc, Nc, E_Pc;
    float hc, Yc, Zc, tc, ji_tc, sqrthc;
    /*Updates coarsed increment with random numbers*/
    void incrementc_update(void);
    /*Step of FKAK formula for corased set of variables*/
    void Stepc(void);
    /*Reset of FKAK */
    void FKAK_resetc(Eigen::VectorXf X0);
    void FKAK_resetc(Eigen::VectorXf X0, float T_start);
    /*Computes one multilevel trayectory of level l for elliptic equations*/
    void Solve_l(Eigen::VectorXf X0);
    /*Computes one multilevel trayectory of level l for parabolic equations*/
    void Solve_l(Eigen::VectorXf X0, float T_start);
    /*Linear regression routine*/
    void regression(int N_sample, float *x, float *y, float &a, float &b);

public:
    float mean, var, std, covar, var_xi, pearson_c, Cl[21];
    uint32_t Nl[21];
    /*Empty initizialization*/
    FKAKSolver(void);
    /*Proper initialization with default random seed equal to 0
    -boundary_value_problem is a BVP object which stores all problem's equations
    -surface_parameters stores the parameters for the boundary construction
    -discretization stores the time discretization for the Stochastic process */
    FKAKSolver(BVP boundary_value_problem, 
      std::vector<float> boundary_parameters,
      float discretization,
      unsigned int seed);
    /*Changes class problem and parameters 
    -boundary_value_problem is a BVP object which stores all problem's equations
    -surface_parameters stores the parameters for the boundary construction
    -discretization stores the time discretization for the Stochastic process */
    void Init(BVP boundary_value_problem, 
        std::vector<float> surface_parameters,
        float discretization,
        unsigned int seed);
    /*Solves an elliptic equation on inital point X0 using FKAK formula 
    approximated by Euler-Maruyama's integrator */
    float Solve_EM(Eigen::VectorXf X0, float tolerance);
    float Solve_GM(Eigen::VectorXf X0, float tolerance);
    /*Solves a parabolic equation on inital point X0 using FKAK formula 
    approximated by Euler-Maruyama's integrator */
    float Solve_EM(Eigen::VectorXf X0, float T_start, float tolerance);
    float Solve_GM(Eigen::VectorXf X0, float T_start, float tolerance);
    /*Solves an elliptic equation on inital point X0 using Multilevel 
    algorithm*/
    float Solve_Multilevel(Eigen::VectorXf X0, int Lmin, int Lmax, int N0, float eps,
                           float alpha_0,float beta_0,float gamma_0);
    float Solve_Multilevel(Eigen::VectorXf X0, float T_start, int Lmin, int Lmax, int N0, float eps,
                           float alpha_0,float beta_0,float gamma_0);
    /*Multilevel test*/
    void mlmc_test(Eigen::VectorXf X0, uint32_t N_test,int L_test, int N0, float *Eps, int Lmin, int Lmax, FILE *fp);
    void mlmc_test(Eigen::VectorXf X0, float T_start, uint32_t N_test,int L_test, int N0, float *Eps, int Lmin, int Lmax, FILE *fp);
    void Elliptic_test(Eigen::VectorXf X0, float tol, FILE *fp);
    void Parabolic_test(Eigen::VectorXf X0, float T_start, float tol, FILE *fp);
    ~FKAKSolver();
};
#endif