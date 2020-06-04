#include "FKAKsolver.hpp"
#include "edepaco.hpp"
#include "sphere.hpp"
int main(){
    Eigen::VectorXf X0;
    X0.resize(2);
    X0[0] = 0.3f; X0[1] = 0.6f;
    std::vector<float> params;
    float rerr = 0.4f;
    params.push_back(1.0f); params.push_back(0.0f); params.push_back(0.0f);
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
    vector_init["F"] = Equation_F;
    scalarN_init["psi"]  = Equation_Psi;
    scalarN_init["varphi"] = Equation_Varphi;
    matrix_init["sigma"] = Equation_sigma;
    bvp.Boundary_init(Sphere, Stopping);
    bvp.BVP_init(2,scalar_init, scalarN_init, vector_init, matrix_init,string_init);
    FILE *fp;
    fp = fopen("Output/parabolic_testGM.txt", "w");
    FKAKSolver Solver[10];
    time_t now = time(NULL);
    char *date = ctime(&now);
    int len = strlen(date);
    date[len-1] = ' ';
    PRINTF2(fp,"\n");
    PRINTF2(fp,"**********************************************************\n");
    PRINTF2(fp,"*** FKAKsolver class file version 0.1    produced by J. Moron ***\n");
    PRINTF2(fp,"*** C++ parabolic_test on %s         ***\n",date);
    PRINTF2(fp,"**********************************************************\n");
    PRINTF2(fp,"\n");
    PRINTF2(fp,"**********************************************************\n");
    PRINTF2(fp,"******************** Convergence tests *******************\n");
    PRINTF2(fp,"*******************Gobet-Menozi Integrator*******************\n");
    PRINTF2(fp,"\n      h  average  variance     err    ");
    PRINTF2(fp,"   rerr        cost \n--------------------------");
    PRINTF2(fp,"-------------------------------------------------------------\n");

    float h = 0.01;
    #pragma omp parallel for    
    for(unsigned int i = 0; i < 10; i++){
        Solver[i].Init(bvp, params,(float) h*pow(2,-1.0*i),i);
        Solver[i].Parabolic_test(X0, 1.0f, rerr, fp);
    }
    fclose(fp);
}