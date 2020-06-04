#include "FKAKsolver.hpp"

FKAKSolver::FKAKSolver()
{
    //Default values
    h = 0.001;
    sqrth = sqrt(h);
    t = 0;
    flux = 0;
    T = gsl_rng_default; //Mt1997
    rng = gsl_rng_alloc(T);
}

FKAKSolver::FKAKSolver(BVP boundary_value_problem, 
    std::vector<float> surface_parameters,
    float discretization,
    unsigned int seed)
{
    Init(boundary_value_problem, surface_parameters, discretization, seed);
    
}


void FKAKSolver::Init(BVP boundary_value_problem, 
    std::vector<float> surface_parameters,
    float discretization,
    unsigned int seed)
{  
    h = discretization;
    sqrth = sqrt(h);
    params = new float[surface_parameters.size()];
    for(unsigned int i = 0; i < surface_parameters.size(); i++) params[i] = surface_parameters[i];
    bvp = boundary_value_problem;
    t = 0;
    flux = 1;
    option = 2;
    N_rngcalls = 0;
    N_trayectories = 0;
    T = gsl_rng_default; //Mt1997
    rng = gsl_rng_alloc(T);
    gsl_rng_set(rng, seed);
}

inline void FKAKSolver::increment_update(void){
    for(int j = 0; j < increment.size(); j++)
    {
        increment(j) = (float) gsl_ran_gaussian_ziggurat(rng,1)*sqrth;
    }
}

inline void FKAKSolver::incrementc_update(void){
    for(int j = 0; j < incrementc.size(); j++)
    {
        incrementc(j) = (float) gsl_ran_gaussian_ziggurat(rng,1)*sqrthc;
    }
}

void FKAKSolver::FKAK_reset(Eigen::VectorXf X0)
{
    xi = 0.0f;
    X = X0;
    Y = 1.0f;
    Z = 0.0f;
    ji_t = 0.0f;
}

void FKAKSolver::FKAK_resetc(Eigen::VectorXf X0){
    xi = 0.0f;
    Xc = X0;
    Yc = 1.0f;
    Zc = 0.0f;
    ji_tc = 0.0f;
}

void FKAKSolver::FKAK_reset(Eigen::VectorXf X0, float T_start)
{
    FKAK_reset(X0);
    ji_t = 0.0f;
    t = T_start;
}

void FKAKSolver::FKAK_resetc(Eigen::VectorXf X0, float T_start){
    FKAK_resetc(X0);
    ji_tc = 0.0f;
    t = T_start;
}

void FKAKSolver::Update_Stat(float sol_0, float xi)
{

    float sol = sol_0 + xi;
    sums[0] += sol;
    sums[1] += sol*sol;
    sums[2] += sol_0;
    sums[3] += sol_0*sol_0;
    sums[4] += xi;
    sums[5] += xi*xi;
    sums[6] += xi*sol_0;
    N_trayectories ++;
    float aux = 1.0f/N_trayectories;
    mean = sums[0] * aux;
    var = sums[1]*aux - mean*mean;
    std = sqrt(aux*var);

}

inline void FKAKSolver::Step(void){
    
    sigma = bvp.sigma.Value(X,t);
    Z += bvp.f.Value(X,t)*Y*h + bvp.psi.Value(X,N,t)*Y*ji_t;
    Y += bvp.c.Value(X,t)*Y*h + bvp.varphi.Value(X,N,t) * Y * ji_t;
    X += bvp.b.Value(X,t)*h + sigma*increment;
    t -= h;
}

inline void FKAKSolver::Stepc(void)
{
    
    sigma = bvp.sigma.Value(Xc,tc);
    Zc += bvp.f.Value(Xc,tc)*Y*hc + bvp.psi.Value(X,N,t)*Yc*ji_tc;
    Yc += bvp.c.Value(Xc,tc)*Y*hc + bvp.varphi.Value(X,N,t) * Yc * ji_tc;
    Xc += bvp.b.Value(Xc,tc)*h + sigma*incrementc;
    tc -= hc;
}

inline void FKAKSolver::VR_Step(void)
{   
    sigma = bvp.sigma.Value(X,t);
    mu = bvp.mu.Value(X,t);
    Z += bvp.f.Value(X,t)*Y*h + bvp.psi.Value(X,N,t)*Y*ji_t;
    Y += bvp.c.Value(X,t)*Y*h + Y*bvp.mu.Value(X,t).transpose()*increment + bvp.varphi.Value(X,N,t) * Y * ji_t;
    X += (bvp.b.Value(X,t) - sigma*mu)*h + sigma*increment;
    t -= h;
}

inline void FKAKSolver::CV_Step(void)
{   
    sigma = bvp.sigma.Value(X,t);
    xi+= Y*(sigma.transpose()*bvp.F.Value(X,t)).dot(increment);
    Z += bvp.f.Value(X,t)*Y*h + bvp.psi.Value(X,N,t)*Y*ji_t;
    Y += bvp.c.Value(X,t)*Y*h + bvp.varphi.Value(X,N,t) * Y * ji_t;
    X += bvp.b.Value(X,t)*h + sigma*increment;
    t -= h;
}

inline void FKAKSolver::VR_CV_Step(void)
{   
    sigma = bvp.sigma.Value(X,t);
    mu = bvp.mu.Value(X,t);
    xi+= Y*(sigma.transpose()*bvp.F.Value(X,t)).dot(increment);
    Z += bvp.f.Value(X,t)*Y*h + bvp.psi.Value(X,N,t)*Y*ji_t;
    Y += bvp.c.Value(X,t)*Y*h + Y*mu.transpose()*increment + bvp.varphi.Value(X,N,t) * Y * ji_t;
    X += (bvp.b.Value(X,t) - sigma*mu)*h + sigma*increment - N*ji_t;
    t += -h;
}


float FKAKSolver::Solve_EM(Eigen::VectorXf X0, float tolerance)
{   
    for(unsigned int i = 0; i < 7; i++) sums[i] = 0.0;
    increment.resize(X0.size());
    X.resize(X0.size());
    E_P.resize(X0.size());
    N.resize(X0.size());
    float sol_0, sol_a = bvp.u.Value(X0,t);
    //Sample test
    for(unsigned int i = 0; i < 100; i++)
    {
        FKAK_reset(X0);

        do{
            increment_update();
            N_rngcalls += X0.size();
            VR_CV_Step();
        }while(Inside_Elliptic(X,N,E_P));
        sol_0 = Y*bvp.g.Value(E_P,t)+Z;
        Update_Stat(sol_0,xi);

    }
    //Trayectories are computed till the MSE we want to obtain is achieved
    do{
        FKAK_reset(X0);
        do{
            increment_update();
            N_rngcalls += X0.size();
            VR_CV_Step();
        }while(Inside_Elliptic(X,N,E_P));
        sol_0 = Y*bvp.g.Value(E_P,t)+Z;
        Update_Stat(sol_0,xi);
    }while(2.0f*std > tolerance*fabs(mean-sol_a));

    double aux = 1.0f/N_trayectories;
    covar =  (sums[6] - sums[2]*sums[4]*aux)*aux;
    var_xi = sums[5]*aux - sums[4]*sums[4]*aux*aux;
    pearson_c = covar/sqrt((sums[3]*aux - (sums[0]-sums[4])*
    (sums[0]-sums[4])*aux*aux)*var_xi);

    return mean;
}

float FKAKSolver::Solve_EM(Eigen::VectorXf X0, float T_start, float tolerance)
{   
    for(unsigned int i = 0; i < 7; i++) sums[i] = 0.0;
    increment.resize(X0.size());
    X.resize(X0.size());
    E_P.resize(X0.size());
    N.resize(X0.size());
    double aux;
    float sol_0, sol_a = bvp.u.Value(X0,T_start);
    //Sample test
    for(unsigned int i = 0; i < 100; i++)
    {
        FKAK_reset(X0, T_start);

        do{
            increment_update();
            N_rngcalls += X0.size();
            VR_CV_Step();
        }while(Inside_Parabolic(X,N,E_P,t));

        if(status == time_out){
          sol_0 = Y*bvp.p.Value(X,0.0f)+Z;
        } else {
          sol_0 = Y*bvp.g.Value(E_P,t)+Z;
        }
        Update_Stat(sol_0,xi);
    }
    //Trayectories are computed till the MSE we want to obtain is achieved
    do{
        FKAK_reset(X0, T_start);
        do{
            increment_update();
            N_rngcalls += X0.size();
            VR_CV_Step();
        }while(Inside_Parabolic(X,N,E_P,t));

        if(status == time_out){
          sol_0 = Y*bvp.p.Value(X,0.0f)+Z;
        } else {
          sol_0 = Y*bvp.g.Value(E_P,t)+Z;
        }
        Update_Stat(sol_0,xi);
    }while(2.0f*std > tolerance*fabs(mean-sol_a));
    aux = 1.0f/N_trayectories;
    covar =  (sums[6] - sums[2]*sums[4]*aux)*aux;
    var_xi = sums[5]*aux - sums[4]*sums[4]*aux*aux;
    pearson_c = covar/sqrt((sums[3]*aux - (sums[0]-sums[4])*
    (sums[0]-sums[4])*aux*aux)*var_xi);

    return mean;
}

float FKAKSolver::Solve_GM(Eigen::VectorXf X0, float tolerance)
{   
    for(unsigned int i = 0; i < 7; i++) sums[i] = 0.0;
    increment.resize(X0.size());
    X.resize(X0.size());
    E_P.resize(X0.size());
    N.resize(X0.size());
    float sol_0, sol_a = bvp.u.Value(X0,t);
    //Sample test
    for(unsigned int i = 0; i < 100; i++)
    {
        FKAK_reset(X0);

        do{
            increment_update();
            N_rngcalls += X0.size();
            VR_CV_Step();
        }while(Inside_Elliptic(X,N,E_P,sigma));
        sol_0 = Y*bvp.g.Value(E_P,t)+Z;
        Update_Stat(sol_0,xi);

    }
    //Trayectories are computed till the MSE we want to obtain is achieved
    do{
        FKAK_reset(X0);
        do{
            increment_update();
            N_rngcalls += X0.size();
            VR_CV_Step();
        }while(Inside_Elliptic(X,N,E_P,sigma));
        sol_0 = Y*bvp.g.Value(E_P,t)+Z;
        Update_Stat(sol_0,xi);
    }while(2.0f*std > tolerance*fabs(mean-sol_a));

    double aux = 1.0f/N_trayectories;
    covar =  (sums[6] - sums[2]*sums[4]*aux)*aux;
    var_xi = sums[5]*aux - sums[4]*sums[4]*aux*aux;
    pearson_c = covar/sqrt((sums[3]*aux - (sums[0]-sums[4])*
    (sums[0]-sums[4])*aux*aux)*var_xi);

    return mean;
}

float FKAKSolver::Solve_GM(Eigen::VectorXf X0, float T_start, float tolerance)
{   
    for(unsigned int i = 0; i < 7; i++) sums[i] = 0.0;
    increment.resize(X0.size());
    X.resize(X0.size());
    E_P.resize(X0.size());
    N.resize(X0.size());
    double aux;
    float sol_0, sol_a = bvp.u.Value(X0,T_start);
    //Sample test
    for(unsigned int i = 0; i < 100; i++)
    {
        FKAK_reset(X0, T_start);

        do{
            increment_update();
            N_rngcalls += X0.size();
            Step();
        }while(Inside_Parabolic(X,N,E_P,t,sqrth,sigma));

        if(status == time_out){
          sol_0 = Y*bvp.p.Value(X,0.0f)+Z;
        } else {
          sol_0 = Y*bvp.g.Value(E_P,t)+Z;
        }
        Update_Stat(sol_0,xi);
    }
    //Trayectories are computed till the MSE we want to obtain is achieved
    do{
        FKAK_reset(X0, T_start);
        do{
            increment_update();
            N_rngcalls += X0.size();
            VR_CV_Step();
        }while(Inside_Parabolic(X,N,E_P,t));

        if(status == time_out){
          sol_0 = Y*bvp.p.Value(X,0.0f)+Z;
        } else {
          sol_0 = Y*bvp.g.Value(E_P,t)+Z;
        }
        Update_Stat(sol_0,xi);
    }while(2.0f*std > tolerance*fabs(mean-sol_a));
    aux = 1.0f/N_trayectories;
    covar =  (sums[6] - sums[2]*sums[4]*aux)*aux;
    var_xi = sums[5]*aux - sums[4]*sums[4]*aux*aux;
    pearson_c = covar/sqrt((sums[3]*aux - (sums[0]-sums[4])*
    (sums[0]-sums[4])*aux*aux)*var_xi);

    return mean;
}

bool FKAKSolver::Inside_Elliptic(Eigen::VectorXf & position, 
                                 Eigen::VectorXf & normal, 
                                 Eigen::VectorXf & exitpoint){
  /*True if the proyection of the previous point on the boundary is
  associated with stopping BC's false otherwise*/
  bool stoppingbc = bvp.boundary.stop(exitpoint);
  //Distance to the boundary 
  float dist = bvp.boundary.Dist(params, position, exitpoint,normal);
  if( dist < 0.0f){

    status = in;

  }else{

    if(stoppingbc){

      status = stop;

  } else {

      status = reflect;

      }

  }
    switch(status){
      
      case stop:
        ji_t = 0.0f;
        position = exitpoint;
        return false;
        break;

      case reflect:
        ji_t = dist;
        position = exitpoint;
        return true;
        break;

      default:
        ji_t = 0.0f;
        return true;
        break;
    }
  
}

bool FKAKSolver::Inside_Parabolic(Eigen::VectorXf & position, 
                                  Eigen::VectorXf & normal, 
                                  Eigen::VectorXf & exitpoint,
                                  float &time){
    bool stoppingbc = bvp.boundary.stop(exitpoint);
    float dist = bvp.boundary.Dist(params, position, exitpoint, normal);
    if( dist < 0.0f){

        if (time > 0.0f) {

            status = in;

        } else {

            status = time_out;

        }
    }

    else{
        
        if (time <= 0.0f) {

            status = time_out;

        } else {

            if(stoppingbc){

                status = stop;

            } else {

                status = reflect;

            }
        }
    }

    switch(status){
      case stop:
        ji_t = 0.0f;
        position = exitpoint;
        return false;
        break;

      case reflect:
        ji_t = dist;
        position = exitpoint;
        return true;
        break;

      case time_out:
        ji_t = 0.0f;
        time = 0.0f;
        return false;
        break;

      default:
        ji_t = 0.0f;
        return true;
        break;
    }
  
}

bool FKAKSolver::Inside_Elliptic(Eigen::VectorXf & position, 
                                 Eigen::VectorXf & normal, 
                                 Eigen::VectorXf & exitpoint,
                                 Eigen::MatrixXf & sigma){

    bool stoppingbc = bvp.boundary.stop(exitpoint);
    float dist = bvp.boundary.Dist(params, position, exitpoint, normal);
    if(stoppingbc){
      if( dist < -0.5826*(normal.transpose()*sigma).norm()*sqrth){

        status = in;

      }else{

        status = stop;
      }

    } else {

      if( dist < 0.0f){

        status = in;

      }else{

        status = reflect;
      }

    }
    switch(status){
      
      case stop:
        ji_t = 0.0f;
        position = exitpoint;
        return false;
        break;

      case reflect:
        ji_t = dist;
        position = exitpoint;
        return true;
        break;

      default:
        ji_t = 0.0f;
        return true;
        break;
    }
  
}

bool FKAKSolver::Inside_Parabolic(Eigen::VectorXf & position, 
                         Eigen::VectorXf & normal, 
                         Eigen::VectorXf & exitpoint,
                         float & time,
                         float & sqrdt,
                         Eigen::MatrixXf & sigma){
    bool stoppingbc = bvp.boundary.stop(exitpoint);
    float dist = bvp.boundary.Dist(params, position, exitpoint, normal);
    if(stoppingbc){
      if( dist < -0.5826*(normal.transpose()*sigma).norm()*sqrdt){

          if (time > 0.0f) {

              status = in;

          } else {

              status = time_out;

          }
      } else {

        if (time > 0.0f) {

              status = stop;

          } else {

              position = exitpoint;
              status = time_out;

          }

      }

    }else{
        
      if( dist <= 0){

          if (time > 0.0f) {

              status = in;

          } else {
              status = time_out;

          }
      } else {

        if (time > 0.0f) {
          
              status = reflect;

          } else {

              position = exitpoint;
              status = time_out;

          }

      }
    }

    switch(status){
      case stop:
        ji_t = 0.0f;
        position = exitpoint;
        return false;
        break;

      case reflect:
        ji_t = dist;
        position = exitpoint;
        return true;
        break;

      case time_out:
        ji_t = 0.0f;
        time = 0.0f;
        return false;
        break;

      default:
        ji_t = 0.0f;
        return true;
        break;
    }
  
}

void FKAKSolver::Elliptic_test(Eigen::VectorXf X0, float tol, FILE *fp){
   // current date/time based on current system
  float sol_n = 0.0f, sol_a = 0.0f;
  sol_n = Solve_EM(X0, tol);
  //sol_n = Solve_GM(X0, tol);
  sol_a = bvp.u.Value(X0,t);
  PRINTF2(fp,"%f  %f  %f  %f  %f  %u\n", h, mean, var, fabs(sol_n- sol_a), fabs((sol_n-sol_a)/sol_a), N_rngcalls);

}

void FKAKSolver::Parabolic_test(Eigen::VectorXf X0, float T_start, float tol, FILE *fp){
   // current date/time based on current system
  float sol_n = 0.0f, sol_a = 0.0f;
  //sol_n = Solve_EM(X0, T_start, tol);
  sol_n = Solve_GM(X0, T_start, tol);
  sol_a = bvp.u.Value(X0,T_start);
  PRINTF2(fp,"%f  %f  %f  %f  %f  %u\n", h, mean, var, fabs(sol_n- sol_a), fabs((sol_n-sol_a)/sol_a), N_rngcalls);

}

float FKAKSolver::Solve_Multilevel(Eigen::VectorXf X0, int Lmin, int Lmax, 
                int N0, float eps,float alpha_0,float beta_0,float gamma_0) 
{

  double  suml[3][21];
  float  ml[21], Vl[21], NlCl[21], x[21], y[21],
         alpha, beta, gamma, sum, theta;
  int    L, converged;

  int    diag = 1;  // diagnostics, set to 0 for none 

  //
  // check input parameters
  //

  if (Lmin<2) {
    fprintf(stderr,"error: needs Lmin >= 2 \n");
    exit(1);
  }
  if (Lmax<Lmin) {
    fprintf(stderr,"error: needs Lmax >= Lmin \n");
    exit(1);
  }

  if (N0<=0 || eps<=0.0f) {
    fprintf(stderr,"error: needs N>0, eps>0 \n");
    exit(1);
  }

  //
  // initialisation
  //

  alpha = fmax(0.0f,alpha_0);
  beta  = fmax(0.0f,beta_0);
  gamma = fmax(0.0f,gamma_0);
  theta = 0.25f;             // MSE split between bias^2 and variance

  L = Lmin;
  converged = 0;

  for(int l=0; l<=Lmax; l++) {
    Nl[l]   = 0;
    Cl[l]   = powf(2.0f,(float)l*gamma);
    NlCl[l] = 0.0f;

    for(int n=0; n<3; n++) suml[n][l] = 0.0;
  }

  for(int l=0; l<=Lmin; l++) dNl[l] = N0;

  //
  // main loop
  //

  while (!converged) {

    //
    // update sample sums
    //

    for (l=0; l<=L; l++) {
      if (diag) printf(" %d ",dNl[l]);

      if (dNl[l]>0) {
        Solve_l(X0);
        suml[0][l] += (float) dNl[l];
        suml[1][l] += sums[1];
        suml[2][l] += sums[2];
        NlCl[l]    += sums[0];  // sum total cost
      }
    }
    if (diag) printf(" \n");

    //
    // compute absolute average, variance and cost,
    // correct for possible under-sampling,
    // and set optimal number of new samples
    //

    sum = 0.0f;

    for (int l=0; l<=L; l++) {
      ml[l] = fabs(suml[1][l]/suml[0][l]);
      Vl[l] = fmaxf(suml[2][l]/suml[0][l] - ml[l]*ml[l], 0.0f);
      if (gamma_0 <= 0.0f) Cl[l] = NlCl[l] / suml[0][l];

      if (l>1) {
        ml[l] = fmaxf(ml[l],  0.5f*ml[l-1]/powf(2.0f,alpha));
        Vl[l] = fmaxf(Vl[l],  0.5f*Vl[l-1]/powf(2.0f,beta));
      }

      sum += sqrtf(Vl[l]*Cl[l]);
    }

    for (int l=0; l<=L; l++) {
      dNl[l] = ceilf( fmaxf( 0.0f, 
                       sqrtf(Vl[l]/Cl[l])*sum/((1.0f-theta)*eps*eps)
                     - suml[0][l] ) );
    }
 
    //
    // use linear regression to estimate alpha, beta, gamma if not given
    //

    if (alpha_0 <= 0.0f) {
      for (int l=1; l<=L; l++) {
        x[l-1] = l;
        y[l-1] = - log2f(ml[l]);
      }
      regression(L,x,y,alpha,sum);
      if (diag) printf(" alpha = %f \n",alpha);
    }

    if (beta_0 <= 0.0f) {
      for (int l=1; l<=L; l++) {
        x[l-1] = l;
        y[l-1] = - log2f(Vl[l]);
      }
      regression(L,x,y,beta,sum);
      if (diag) printf(" beta = %f \n",beta);
    }

     if (gamma_0 <= 0.0f) {
      for (int l=1; l<=L; l++) {
        x[l-1] = l;
        y[l-1] = log2f(Cl[l]);
      }
      regression(L,x,y,gamma,sum);
      if (diag) printf(" gamma = %f \n",gamma);
    }

    //
    // if (almost) converged, estimate remaining error and decide 
    // whether a new level is required
    //

    sum = 0.0;
      for (int l=0; l<=L; l++)
        sum += fmaxf(0.0f, (float)dNl[l]-0.01f*suml[0][l]);

    if (sum==0) {
      if (diag) printf(" achieved variance target \n");

      converged = 1;
      float rem = ml[L] / (powf(2.0f,alpha)-1.0f);

      if (rem > sqrtf(theta)*eps) {
        if (L==Lmax)
          printf("*** failed to achieve weak convergence *** \n");
        else {
          converged = 0;
          L++;
          Vl[L] = Vl[L-1]/powf(2.0f,beta);
          Cl[L] = Cl[L-1]*powf(2.0f,gamma);

          if (diag) printf(" L = %d \n",L);

          sum = 0.0f;
          for (int l=0; l<=L; l++) sum += sqrtf(Vl[l]*Cl[l]);
          for (int l=0; l<=L; l++)
            dNl[l] = ceilf( fmaxf( 0.0f, 
                            sqrtf(Vl[l]/Cl[l])*sum/((1.0f-theta)*eps*eps)
                          - suml[0][l] ) );
        }
      }
    }
  }

  //
  // finally, evaluate multilevel estimator and set outputs
  //

  float P = 0.0f;
  for (int l=0; l<=L; l++) {
    P    += suml[1][l]/suml[0][l];
    Nl[l] = suml[0][l];
    Cl[l] = NlCl[l] / Nl[l];
  }

  return P;
}

float FKAKSolver::Solve_Multilevel(Eigen::VectorXf X0, float T_start, int Lmin, int Lmax, 
                int N0, float eps,float alpha_0,float beta_0,float gamma_0) 
{

  double  suml[3][21];
  float  ml[21], Vl[21], NlCl[21], x[21], y[21],
         alpha, beta, gamma, sum, theta;
  int    converged;

  int    diag = 1;  // diagnostics, set to 0 for none 

  //
  // check input parameters
  //

  if (Lmin<2) {
    fprintf(stderr,"error: needs Lmin >= 2 \n");
    exit(1);
  }
  if (Lmax<Lmin) {
    fprintf(stderr,"error: needs Lmax >= Lmin \n");
    exit(1);
  }

  if (N0<=0 || eps<=0.0f) {
    fprintf(stderr,"error: needs N>0, eps>0 \n");
    exit(1);
  }

  //
  // initialisation
  //

  alpha = fmax(0.0f,alpha_0);
  beta  = fmax(0.0f,beta_0);
  gamma = fmax(0.0f,gamma_0);
  theta = 0.25f;             // MSE split between bias^2 and variance

  L = Lmin;
  converged = 0;

  for(int l=0; l<=Lmax; l++) {
    Nl[l]   = 0;
    Cl[l]   = powf(2.0f,(float)l*gamma);
    NlCl[l] = 0.0f;

    for(int n=0; n<3; n++) suml[n][l] = 0.0;
  }

  for(int l=0; l<=Lmin; l++) dNl[l] = N0;

  //
  // main loop
  //

  while (!converged) {

    //
    // update sample sums
    //

    for (l=0; l<=L; l++) {
      if (diag) printf(" %d ",dNl[l]);

      if (dNl[l]>0) {
        Solve_l(X0, T_start);
        suml[0][l] += (float) dNl[l];
        suml[1][l] += sums[1];
        suml[2][l] += sums[2];
        NlCl[l]    += sums[0];  // sum total cost
      }
    }
    if (diag) printf(" \n");

    //
    // compute absolute average, variance and cost,
    // correct for possible under-sampling,
    // and set optimal number of new samples
    //

    sum = 0.0f;

    for (int l=0; l<=L; l++) {
      ml[l] = fabs(suml[1][l]/suml[0][l]);
      Vl[l] = fmaxf(suml[2][l]/suml[0][l] - ml[l]*ml[l], 0.0f);
      if (gamma_0 <= 0.0f) Cl[l] = NlCl[l] / suml[0][l];

      if (l>1) {
        ml[l] = fmaxf(ml[l],  0.5f*ml[l-1]/powf(2.0f,alpha));
        Vl[l] = fmaxf(Vl[l],  0.5f*Vl[l-1]/powf(2.0f,beta));
      }

      sum += sqrtf(Vl[l]*Cl[l]);
    }

    for (int l=0; l<=L; l++) {
      dNl[l] = ceilf( fmaxf( 0.0f, 
                       sqrtf(Vl[l]/Cl[l])*sum/((1.0f-theta)*eps*eps)
                     - suml[0][l] ) );
    }
 
    //
    // use linear regression to estimate alpha, beta, gamma if not given
    //

    if (alpha_0 <= 0.0f) {
      for (int l=1; l<=L; l++) {
        x[l-1] = l;
        y[l-1] = - log2f(ml[l]);
      }
      regression(L,x,y,alpha,sum);
      if (diag) printf(" alpha = %f \n",alpha);
    }

    if (beta_0 <= 0.0f) {
      for (int l=1; l<=L; l++) {
        x[l-1] = l;
        y[l-1] = - log2f(Vl[l]);
      }
      regression(L,x,y,beta,sum);
      if (diag) printf(" beta = %f \n",beta);
    }

     if (gamma_0 <= 0.0f) {
      for (int l=1; l<=L; l++) {
        x[l-1] = l;
        y[l-1] = log2f(Cl[l]);
      }
      regression(L,x,y,gamma,sum);
      if (diag) printf(" gamma = %f \n",gamma);
    }

    //
    // if (almost) converged, estimate remaining error and decide 
    // whether a new level is required
    //

    sum = 0.0;
      for (int l=0; l<=L; l++)
        sum += fmaxf(0.0f, (float)dNl[l]-0.01f*suml[0][l]);

    if (sum==0) {
      if (diag) printf(" achieved variance target \n");

      converged = 1;
      float rem = ml[L] / (powf(2.0f,alpha)-1.0f);

      if (rem > sqrtf(theta)*eps) {
        if (L==Lmax)
          printf("*** failed to achieve weak convergence *** \n");
        else {
          converged = 0;
          L++;
          Vl[L] = Vl[L-1]/powf(2.0f,beta);
          Cl[L] = Cl[L-1]*powf(2.0f,gamma);

          if (diag) printf(" L = %d \n",L);

          sum = 0.0f;
          for (int l=0; l<=L; l++) sum += sqrtf(Vl[l]*Cl[l]);
          for (int l=0; l<=L; l++)
            dNl[l] = ceilf( fmaxf( 0.0f, 
                            sqrtf(Vl[l]/Cl[l])*sum/((1.0f-theta)*eps*eps)
                          - suml[0][l] ) );
        }
      }
    }
  }

  //
  // finally, evaluate multilevel estimator and set outputs
  //

  float P = 0.0f;
  for (int l=0; l<=L; l++) {
    P    += suml[1][l]/suml[0][l];
    Nl[l] = suml[0][l];
    Cl[l] = NlCl[l] / Nl[l];
  }

  return P;
}

void FKAKSolver::regression(int N_sample, float *x, float *y, float &a, float &b){

  float sum0=0.0f, sum1=0.0f, sum2=0.0f, sumy0=0.0f, sumy1=0.0f;

  for (int i=0; i<N_sample; i++) {
    sum0  += 1.0f;
    sum1  += x[i];
    sum2  += x[i]*x[i];

    sumy0 += y[i];
    sumy1 += y[i]*x[i];
  }

  a = (sum0*sumy1 - sum1*sumy0) / (sum0*sum2 - sum1*sum1);
  b = (sum2*sumy0 - sum1*sumy1) / (sum0*sum2 - sum1*sum1);
}

void FKAKSolver::Solve_l(Eigen::VectorXf X0){
  int16_t split, in_f, in_c, set_f, set_c;
  uint32_t counter = 0;
  float  T, h0, Pf, Pc, dP;
  double t_split;
  Eigen::VectorXf X_split;
    increment.resize(X0.size());
    incrementc.resize(X0.size());
    N.resize(X0.size());
    E_P.resize(X0.size());
    M = 4;
    D = X0.size();
  // printf(" **** l, Nl = %d, %d **** \n",l,N);

    h0 = 0.1;
    T  = 1.0;

    h = h0 / powf((float)M,(float) l);
    sqrth = sqrt(h);
    hc = h0 / powf((float)M,(float) (l-1));
    sqrthc = sqrt(hc);
  for (unsigned int k=0; k<7; k++) sums[k] = 0.0;
  if (option==1 || option==3) {
    //  split = 4 * (1<<l);
    split = 1<<l;
  }
  else if (option==2) {
    split = 1;
  }
  for (unsigned int np = 0; np < dNl[l]; np++) {

    X = X0;
    Xc = X0;
    t = 0.0;
    tc   = 0.0;
    in_f = 1;
    in_c = 1;
    set_f = 0;
    set_c = 0;

    Pf = 0.0;
    Pc = 0.0;

    // level 0
    FKAK_reset(X0);
    if (l==0) {
      do {
        counter += D;
        increment_update();

        Step();

        if(bvp.boundary.Dist(params,X,E_P,N) >= 0.0f) in_f = 0;
      } while (in_f);
      std::cout << bvp.g.Value(E_P,t) << "\n";
      Pf = Y*bvp.g.Value(E_P,t)+Z;
      std::cout << Pf << std::endl;
    }

    // level l>0

    else {
      FKAK_reset(X0);
      FKAK_resetc(X0);
      do {

        for (int d=0; d<D; d++) incrementc(d) = 0.0f;

        for (int m=0; m<M; m++) {
          counter += D;
          increment_update();
          incrementc += increment;
          Step();

          if(bvp.boundary.Dist(params,X,E_P,N) >= 0.0f) in_f = 0;

          if ( (!in_f) && (!set_f) ) {
            set_f = 1;
            Pf = Y*bvp.g.Value(E_P,t)+Z;
          }
      	}

        Stepc();

        if(bvp.boundary.Dist(params,Xc,E_P,N) >= 0.0f) in_c = 0;

        if ( (!in_c) && (!set_c) ) {
          set_c = 1;
          Pc    = Yc*bvp.g.Value(E_P,t)+Zc;
        }

      } while (in_f && in_c);

    // split continuation paths 

      if (in_f) {
        t_split = t;
        X_split = X;

      for (int s=0; s<split; s++) {
	  // reset state at split
          in_f = 1;
          t = t_split;
          FKAK_reset(X_split);

	  // continue until exit
          do {
            counter += D;
            increment_update();
            Step();
            if(bvp.boundary.Dist(params,X,E_P,N) >= 0.0f) in_f = 0;

          } while (in_f);

          Pf += (Y*bvp.g.Value(E_P,t)+Z)/ ((float) split);
        }
      }

      if (in_c) {
        t_split = tc;
        X_split = Xc;

        for (int s=0; s<split; s++) {
	  // reset state at split
          in_c = 1;
          tc = t_split;
          FKAK_resetc(X_split);

	  // continue until exit
          do {
            counter += D;
            incrementc_update();

            Stepc();

            if(bvp.boundary.Dist(params,Xc,E_P,N) >= 0.0f) in_c = 0;

          } while (in_c);

          Pc += (Yc*bvp.g.Value(E_P,t)+Zc) / ((float) split);
        }
      }
    }

    dP = Pf-Pc;

    sums[0] += counter; // add number of RNG calls
    sums[1] += dP;
    sums[2] += dP*dP;
    sums[3] += dP*dP*dP;
    sums[4] += dP*dP*dP*dP;
    sums[5] += Pf;
    sums[6] += Pf*Pf;
  }
}

void FKAKSolver::Solve_l(Eigen::VectorXf X0, float T_start){
  int16_t split, in_f, in_c, set_f, set_c;
  uint32_t counter = 0;
  float  h0, Pf, Pc, dP;
  double t_split, Y_split, Z_split;
  Eigen::VectorXf X_split;
    increment.resize(X0.size());
    incrementc.resize(X0.size());
    N.resize(X0.size());
    E_P.resize(X0.size());
    Nc.resize(X0.size());
    E_Pc.resize(X0.size());
    M = 4;
    D = X0.size();
  // printf(" **** l, Nl = %d, %d **** \n",l,N);

    h0 = 0.1;
    h = h0 / powf((float)M,(float) l);
    sqrth = sqrt(h);
    hc = h0 / powf((float)M,(float) (l-1));
    sqrthc = sqrt(hc);
  for (unsigned int k=0; k<7; k++) sums[k] = 0.0;
  if (option==1 || option==3) {
    //  split = 4 * (1<<l);
    split = 1<<l;
    //split = 0;
  }
  else if (option==2) {
    split = 1;
  }
  for (unsigned int np = 0; np < dNl[l]; np++) {

    X = X0;
    Xc = X0;
    t = T_start;
    tc   = T_start;
    in_f = 1;
    in_c = 1;
    set_f = 0;
    set_c = 0;

    Pf = 0.0;
    Pc = 0.0;

    // level 0
    FKAK_reset(X0, T_start);
    if (l==0) {
      do {
        increment_update();
        counter += D;
        Step();
        if(!Inside_Parabolic(X,N,E_P,t,sqrth,sigma)) in_f = 0;
      } while (in_f);
      if(status == time_out){
          Pf = Y*bvp.p.Value(X,0.0f)+Z;
      } else {
          Pf = Y*bvp.g.Value(E_P,t)+Z;
      }
    }

    // level l>0

    else {

      FKAK_reset(X0, T_start);
      FKAK_resetc(X0, T_start);
      do {

        for (int d=0; d<D; d++) incrementc(d) = 0.0f;

        for (int m=0; m<M; m++){
          counter += D;
          increment_update();
          incrementc += increment;
          Step();
          if(!Inside_Parabolic(X,N,E_P,t,sqrth,sigma)) in_f = 0;
          if ( (!in_f) && (!set_f) ) {
            set_f = 1;
            if(status == time_out){
                Pf = Y*bvp.p.Value(X,0.0f)+Z;
            } else {
                Pf = Y*bvp.g.Value(E_P,t)+Z;
            }
          }
      	}
        Stepc();
        if(!Inside_Parabolic(Xc,Nc,E_Pc,tc,sqrthc,sigma)) in_c = 0;
        if ( (!in_c) && (!set_c) ) {
          set_c = 1;
          if(status == time_out){
                Pc = Yc*bvp.p.Value(Xc,0.0f)+Zc;
          } else {
                Pc = Yc*bvp.g.Value(E_Pc,tc)+Zc;
          }
        }

      } while (in_f && in_c);

    // split continuation paths 

      if (in_f) {
        t_split = t;
        X_split = X;
        Y_split = Y;
        Z_split = Z;
      for (int s=0; s<split; s++) {
	  // reset state at split
          in_f = 1;
          t = t_split;
          FKAK_reset(X_split, t_split);
          Y = Y_split;
          Z = Z_split;
	  // continue until exit
          do {
            counter += D;
            increment_update();
            Step();
            if(!Inside_Parabolic(X,N,E_P,t,sqrth,sigma)) in_f = 0;

          } while (in_f);
          if(status == time_out){
                Pf += (Y*bvp.p.Value(X,0.0f)+Z)/((float) split);
            } else {
                Pf += (Y*bvp.g.Value(E_P,t)+Z)/((float) split);
          }
        }
      }

      if (in_c) {
        t_split = tc;
        X_split = Xc;
        Y_split = Yc;
        Z_split = Zc;

        for (int s=0; s<split; s++) {
	        // reset state at split
          in_c = 1;
          FKAK_resetc(X_split, t_split);
          Yc = Y_split;
          Zc = Z_split;
	        // continue until exit
          do {
            counter += D;
            incrementc_update();

            Stepc();

            if(!Inside_Parabolic(Xc,Nc,E_Pc,tc,sqrthc,sigma)) in_c = 0;

          } while (in_c);
          if(status == time_out){
                Pc += (Yc*bvp.p.Value(Xc,0.0f)+Zc)/ ((float) split);
          } else {
                Pc += (Yc*bvp.g.Value(E_Pc,tc)+Zc)/ ((float) split);
          }
        }
      }
    }

    dP = Pf-Pc;

    sums[0] += counter; // add number of RNG calls
    sums[1] += dP;
    sums[2] += dP*dP;
    sums[3] += dP*dP*dP;
    sums[4] += dP*dP*dP*dP;
    sums[5] += Pf;
    sums[6] += Pf*Pf;
  }
}

void FKAKSolver::mlmc_test(Eigen::VectorXf X0, uint32_t N_test,int L_test, int N0, float *Eps, int Lmin, int Lmax, FILE *fp) {

//
// first, convergence tests
//

  // current date/time based on current system
  time_t now = time(NULL);
  char *date = ctime(&now);
  int len = strlen(date);
  date[len-1] = ' ';

  PRINTF2(fp,"\n");
  PRINTF2(fp,"**********************************************************\n");
  PRINTF2(fp,"*** MLMC file version 0.9     produced by              ***\n");
  PRINTF2(fp,"*** C++ mlmc_test on %s         ***\n",date);
  PRINTF2(fp,"**********************************************************\n");
  PRINTF2(fp,"\n");
  PRINTF2(fp,"**********************************************************\n");
  PRINTF2(fp,"*** Convergence tests, kurtosis, telescoping sum check ***\n");
  PRINTF2(fp,"*** using N =%7d samples                           ***\n",N_test);
  PRINTF2(fp,"**********************************************************\n");
  PRINTF2(fp,"\n l   ave(Pf-Pc)    ave(Pf)   var(Pf-Pc)    var(Pf)");
  PRINTF2(fp,"    kurtosis     check        cost \n--------------------------");
  PRINTF2(fp,"-------------------------------------------------------------\n");

  double sums[7];
  float *cost = (float *)malloc((L_test+1)*sizeof(float));
  float *del1 = (float *)malloc((L_test+1)*sizeof(float));
  float *del2 = (float *)malloc((L_test+1)*sizeof(float));
  float *var1 = (float *)malloc((L_test+1)*sizeof(float));
  float *var2 = (float *)malloc((L_test+1)*sizeof(float));
  float *chk1 = (float *)malloc((L_test+1)*sizeof(float));
  float *kur1 = (float *)malloc((L_test+1)*sizeof(float));

  for (l=0; l<=L_test; l++) {
    dNl[l] = N_test;
    Solve_l(X0);

    for (int m=0; m<7; m++) sums[m] = sums[m]/N_test;

    if (M>0)
      cost[l] = powf((float)M,(float)l);
    else
      cost[l] = sums[0];
    del1[l] = sums[1];
    del2[l] = sums[5];
    var1[l] = fmax(sums[2]-sums[1]*sums[1], 1e-10);
    var2[l] = fmax(sums[6]-sums[5]*sums[5], 1e-10);

    kur1[l]  = (      sums[4]
                - 4.0*sums[3]*sums[1]
                + 6.0*sums[2]*sums[1]*sums[1]
                - 3.0*sums[1]*sums[1]*sums[1]*sums[1] )
             / (var1[l]*var1[l]);

    if (l==0)
      chk1[l] = 0.0f;
    else
      chk1[l] = sqrtf((float) N_test) * 
                fabsf(  del1[l]  +       del2[l-1]  -       del2[l] )
         / (3.0f*(sqrtf(var1[l]) + sqrtf(var2[l-1]) + sqrtf(var2[l])));

    PRINTF2(fp,"%2d  %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e \n",
    l,del1[l],del2[l],var1[l],var2[l],kur1[l],chk1[l],cost[l]);
  }

//
// print out a warning if kurtosis or consistency check looks bad
//

  if (kur1[L] > 100.0f) {
    PRINTF2(fp,"\n WARNING: kurtosis on finest level = %f \n",kur1[L]);
    PRINTF2(fp," indicates MLMC correction dominated by a few rare paths; \n");
    PRINTF2(fp," for information on the connection to variance of sample variances,\n");
    PRINTF2(fp," see http://mathworld.wolfram.com/SampleVarianceDistribution.html \n");
  }

  float max_chk = 0.0f;
  for (int l=0; l<=L; l++) max_chk = fmaxf(max_chk,chk1[l]);
  if (max_chk > 1.0f) {
    PRINTF2(fp,"\n WARNING: maximum consistency error = %f \n",max_chk);
    PRINTF2(fp," indicates identity E[Pf-Pc] = E[Pf] - E[Pc] not satisfied \n");
  }

//
// use linear regression to estimate alpha, beta, gamma
//

  float alpha, beta, gamma, foo;
  float *x = (float *)malloc(L*sizeof(float));
  float *y = (float *)malloc(L*sizeof(float));

  for (int l=1; l<=L; l++) {
    x[l-1] = l;
    y[l-1] = - log2f(fabsf(del1[l]));
  } 
  regression(L,x,y,alpha,foo);

  for (int l=1; l<=L; l++) {
    x[l-1] = l;
    y[l-1] = - log2f(var1[l]);
  } 
  regression(L,x,y,beta,foo);

  for (int l=1; l<=L; l++) {
    x[l-1] = l;
    y[l-1] = log2f(cost[l]);
  } 
  regression(L,x,y,gamma,foo);

  PRINTF2(fp,"\n******************************************************\n");
  PRINTF2(fp,"*** Linear regression estimates of MLMC parameters ***\n");
  PRINTF2(fp,"******************************************************\n");
  PRINTF2(fp,"\n alpha = %f  (exponent for MLMC weak convergence)\n",alpha);
  PRINTF2(fp," beta  = %f  (exponent for MLMC variance) \n",beta);
  PRINTF2(fp," gamma = %f  (exponent for MLMC cost) \n",gamma);

//
// second, mlmc complexity tests
//

  PRINTF2(fp,"\n");
  PRINTF2(fp,"***************************** \n");
  PRINTF2(fp,"*** MLMC complexity tests *** \n");
  PRINTF2(fp,"***************************** \n\n");
  PRINTF2(fp,"  eps       value   mlmc_cost   std_cost  savings     N_l \n");
  PRINTF2(fp,"--------------------------------------------------------- \n");
 
  int i=0;
  int   *Nl = (int *)malloc((Lmax+1)*sizeof(int));
  float *Cl = (float *)malloc((Lmax+1)*sizeof(float));

  while (Eps[i]>0) {
    float eps = Eps[i++];

    float P = Solve_Multilevel(X0, Lmin, Lmax, N0, eps, alpha, beta, gamma);

    float std_cost = 0.0f, mlmc_cost = 0.0f, theta=0.25f;

    for (int l=0; l<=Lmax; l++) {
      if (Nl[l]>0) {
        // printf(" l, Cl, cost = %d  %f  %f \n",l,Cl[l],cost[l]);
        mlmc_cost += Nl[l]*Cl[l];
        if (l<=L)
          std_cost = var2[l]*Cl[l] / ((1.0f-theta)*eps*eps);
        else
          std_cost = var2[L]*Cl[l] / ((1.0f-theta)*eps*eps);
      }
    }

    PRINTF2(fp,"%.4f  %.4e  %.3e  %.3e  %7.2f ",
	    eps, P, mlmc_cost, std_cost, std_cost/mlmc_cost);
    for (int l=0; Nl[l]>0; l++) PRINTF2(fp,"%9d",Nl[l]);
    PRINTF2(fp,"\n");
  }
  PRINTF2(fp,"\n");
}

void FKAKSolver::mlmc_test(Eigen::VectorXf X0, float T_start, uint32_t N_test, int L_test, int N0, float *Eps, int Lmin, int Lmax, FILE *fp) {

//
// first, convergence tests
//

  // current date/time based on current system
  option = 1;
  time_t now = time(NULL);
  char *date = ctime(&now);
  int len = strlen(date);
  date[len-1] = ' ';

  PRINTF2(fp,"\n");
  PRINTF2(fp,"**********************************************************\n");
  PRINTF2(fp,"*** MLMC file version 0.9     produced by              ***\n");
  PRINTF2(fp,"*** C++ mlmc_test on %s         ***\n",date);
  PRINTF2(fp,"**********************************************************\n");
  PRINTF2(fp,"\n");
  PRINTF2(fp,"**********************************************************\n");
  PRINTF2(fp,"*** Convergence tests, kurtosis, telescoping sum check ***\n");
  PRINTF2(fp,"*** using N =%7d samples                           ***\n",N_test);
  PRINTF2(fp,"**********************************************************\n");
  PRINTF2(fp,"\n l   ave(Pf-Pc)    ave(Pf)   var(Pf-Pc)    var(Pf)");
  PRINTF2(fp,"    kurtosis     check        cost \n--------------------------");
  PRINTF2(fp,"-------------------------------------------------------------\n");
  
  float *cost = (float *)malloc((L_test+1)*sizeof(float));
  float *del1 = (float *)malloc((L_test+1)*sizeof(float));
  float *del2 = (float *)malloc((L_test+1)*sizeof(float));
  float *var1 = (float *)malloc((L_test+1)*sizeof(float));
  float *var2 = (float *)malloc((L_test+1)*sizeof(float));
  float *chk1 = (float *)malloc((L_test+1)*sizeof(float));
  float *kur1 = (float *)malloc((L_test+1)*sizeof(float));

  for (l=0; l<=L_test; l++) {
    dNl[l] = N_test;
    Solve_l(X0, T_start);

    for (int m=0; m<7; m++) sums[m] = sums[m]/N_test;

    if (M>0)
      cost[l] = powf((float)M,(float)l);
    else
      cost[l] = sums[0];
    del1[l] = sums[1];
    del2[l] = sums[5];
    var1[l] = fmax(sums[2]-sums[1]*sums[1], 1e-10);
    var2[l] = fmax(sums[6]-sums[5]*sums[5], 1e-10);

    kur1[l]  = (      sums[4]
                - 4.0*sums[3]*sums[1]
                + 6.0*sums[2]*sums[1]*sums[1]
                - 3.0*sums[1]*sums[1]*sums[1]*sums[1] )
             / (var1[l]*var1[l]);

    if (l==0)
      chk1[l] = 0.0f;
    else
      chk1[l] = sqrtf((float) N_test) * 
                fabsf(  del1[l]  +       del2[l-1]  -       del2[l] )
         / (3.0f*(sqrtf(var1[l]) + sqrtf(var2[l-1]) + sqrtf(var2[l])));

    PRINTF2(fp,"%2d  %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e \n",
    l,del1[l],del2[l],var1[l],var2[l],kur1[l],chk1[l],cost[l]);
    L = l;
  }

//
// print out a warning if kurtosis or consistency check looks bad
//

  if (kur1[L] > 100.0f) {
    PRINTF2(fp,"\n WARNING: kurtosis on finest level = %f \n",kur1[L]);
    PRINTF2(fp," indicates MLMC correction dominated by a few rare paths; \n");
    PRINTF2(fp," for information on the connection to variance of sample variances,\n");
    PRINTF2(fp," see http://mathworld.wolfram.com/SampleVarianceDistribution.html \n");
  }

  float max_chk = 0.0f;
  for (int l=0; l<=L; l++) max_chk = fmaxf(max_chk,chk1[l]);
  if (max_chk > 1.0f) {
    PRINTF2(fp,"\n WARNING: maximum consistency error = %f \n",max_chk);
    PRINTF2(fp," indicates identity E[Pf-Pc] = E[Pf] - E[Pc] not satisfied \n");
  }

//
// use linear regression to estimate alpha, beta, gamma
//

  float alpha, beta, gamma, foo;
  float *x, *y;
  x = new float[L];
  y = new float[L];

  for (int l=1; l<=L; l++) {
    x[l-1] = l;
    y[l-1] = - log2f(fabsf(del1[l]));
  } 
  regression(L,x,y,alpha,foo);

  for (int l=1; l<=L; l++) {
    x[l-1] = l;
    y[l-1] = - log2f(var1[l]);
  } 
  regression(L,x,y,beta,foo);

  for (int l=1; l<=L; l++) {
    x[l-1] = l;
    y[l-1] = log2f(cost[l]);
  } 
  regression(L,x,y,gamma,foo);

  PRINTF2(fp,"\n******************************************************\n");
  PRINTF2(fp,"*** Linear regression estimates of MLMC parameters ***\n");
  PRINTF2(fp,"******************************************************\n");
  PRINTF2(fp,"\n alpha = %f  (exponent for MLMC weak convergence)\n",alpha);
  PRINTF2(fp," beta  = %f  (exponent for MLMC variance) \n",beta);
  PRINTF2(fp," gamma = %f  (exponent for MLMC cost) \n",gamma);

//
// second, mlmc complexity tests
//

  PRINTF2(fp,"\n");
  PRINTF2(fp,"***************************** \n");
  PRINTF2(fp,"*** MLMC complexity tests *** \n");
  PRINTF2(fp,"***************************** \n\n");
  PRINTF2(fp,"  eps       value   mlmc_cost   std_cost  savings     N_l \n");
  PRINTF2(fp,"--------------------------------------------------------- \n");
 
  int i=0;
  int   *Nl = (int *)malloc((Lmax+1)*sizeof(int));
  float *Cl = (float *)malloc((Lmax+1)*sizeof(float));

  while (Eps[i]>0) {
    float eps = Eps[i++];

    float P = Solve_Multilevel(X0, T_start, Lmin, Lmax, N0, eps, alpha, beta, gamma);

    float std_cost = 0.0f, mlmc_cost = 0.0f, theta=0.25f;

    for (int l=0; l<=Lmax; l++) {
      if (Nl[l]>0) {
        // printf(" l, Cl, cost = %d  %f  %f \n",l,Cl[l],cost[l]);
        mlmc_cost += Nl[l]*Cl[l];
        if (l<=L)
          std_cost = var2[l]*Cl[l] / ((1.0f-theta)*eps*eps);
        else
          std_cost = var2[L]*Cl[l] / ((1.0f-theta)*eps*eps);
      }
    }

    PRINTF2(fp,"%.4f  %.4e  %.3e  %.3e  %7.2f ",
	    eps, P, mlmc_cost, std_cost, std_cost/mlmc_cost);
    for (int l=0; Nl[l]>0; l++) PRINTF2(fp,"%9d",Nl[l]);
    PRINTF2(fp,"\n");
  }
  PRINTF2(fp,"\n");
  delete x;
  delete y;
}

FKAKSolver::~FKAKSolver()
{
    delete params;
    N.resize(0);
    E_P.resize(0);
    mu.resize(0);
    X.resize(0);
    sigma.resize(0,0);
}