#include "FKAKsolver.hpp"
#include "edepaco.hpp"
#include "cube.hpp"
int main(void){
    Eigen::VectorXf X0;
    X0.resize(3);
    int Lmin = 2, Lmax = 10, L_test = 6, N0 = 200;
    uint32_t N_test = 50000;
    float Eps[5] = {0.0005f, 0.001f, 0.002f, 0.005f, 0.01f};
    X0[0] = 0.0f; X0[1] = 0.0f; X0[2] = 0.0f;
    std::vector<float> params;
    float rerr = 0.4f, T_start = 1.0f;
    params.push_back(1.0f); params.push_back(0.0f); params.push_back(0.0f); params.push_back(0.0f);
    BVP bvp;
    std::map<std::string, pfscalar> scalar_init;
    std::map<std::string, pfscalarN> scalarN_init;
    std::map<std::string, pfvector> vector_init;
    std::map<std::string, pfmatrix> matrix_init;
    std::map<std::string, std::string> string_init;
    //BVP initialization 
    scalar_init["f"] = Equation_f;
    scalar_init["c"] = Equation_c;
    scalar_init["u"] = Equation_u;
    scalar_init["g"] = Equation_g;
    scalar_init["p"] = Equation_p;
    vector_init["b"] = Equation_b;
    matrix_init["sigma"] = Equation_sigma;
    bvp.Boundary_init(Cube, Stopping);
    bvp.BVP_init(3,scalar_init, scalarN_init, vector_init, matrix_init,string_init);
    FKAKSolver Solver(bvp, params, 0.1f, 0);
    FILE *fp;
    fp = fopen("Output/mlmc_test.txt", "w");
    Solver.mlmc_test(X0, T_start, N_test, L_test, N0, Eps, Lmin, Lmax, fp);
    fclose(fp);
}
