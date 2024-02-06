/**
 * @file function.cpp
 * @brief This code uses Ising model Universality class to introduce a critical point and an alternative expansion scheme to Extend the QCD based Eqaution of state to finite densities and match lattice QCD at low densities
 *
 *
 * @author Micheal Kahangirwe <mkahangi@central.uh.edu>
 */

/**
 * @brief main function
 *
 * - function to read inputs T-range, muB-range, step size amd mapping parameters  w,  rho  and muBC from inputs/parameterfile
 * The main procedure is:
 *  - Maps Ising co-ordinates  to Alternative expansion scheme then to QCD co-ordinates
 *  - Merges Lattice to Ising
 *  Computes all the Thermodynamic Observables
 *  - Baryon Density, Pressure, Entropy, Baryon number susceptability, Energy, Speed of Sound
 *
 *
 * @return Out puts thermodynamic observables
 */

#include "function.hpp"
#include <iostream>
#include <string>
#include <filesystem>
#include <cstdio>    // Instead of #include <stdio.h>
#include <cstdlib>   // Instead of #include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <cmath>     // Instead of #include <math.h>
#include <unistd.h>
#include <cstring>   // Instead of #include <string.h>
#include <ctime>     // Instead of #include <time.h>
#include <cerrno>    // Instead of #include <err.h>
#include <iomanip>

using namespace std;

/*** Quardratic mapping compared to Linear mapping (BEST Collaboration) 
 * We relate the Quadratic mapping parameters wp, rhop, alpha12p to Linear mapoping pararmeters w,rho, alpha12 = Pi - alphadiff(BEST) ( ***/
Parameters param;

void read_param(char *param_file) 
{
   
    std::cerr << "Reading parameter file.\n";

    double muBCx, alpha12x, ww, rhox, lowT_outx, highT_outx,T_stepx, lowMU_outx, highMU_outx,muB_stepx;

    std::ifstream file(param_file);
    char line[FILENAME_MAX];
    if (!file.is_open())
        std::cerr << "Could not open " << param_file << std::endl;

    while (file.getline(line, sizeof(line)))
    {
        char *c = strchr(line, '#');
        if (c)
            *c = 0;
        if (sscanf(line, "lowT_out=%lf ", &lowT_outx) == 1)
            continue;
        if (sscanf(line, "highT_out=%lf ", &highT_outx) == 1)
            continue;
        if (sscanf(line, "T_step=%lf ", &T_stepx) == 1)
            continue;
        if (sscanf(line, "lowMU_out=%lf ", &lowMU_outx) == 1)
            continue;
        if (sscanf(line, "highMU_out=%lf ", &highMU_outx) == 1)
            continue;
        if (sscanf(line, "muB_step=%lf ", &muB_stepx) == 1)
            continue;
        if (sscanf(line, "muBC=%lf ", &muBCx) == 1)
            continue;
        if (sscanf(line, "alpha12=%lf ", &alpha12x) == 1)
            continue;
        if (sscanf(line, "ww=%lf ", &ww) == 1)
            continue;
        if (sscanf(line, "rho=%lf ", &rhox) == 1)
            continue;

        std::cerr << "Line: <" << line << "> was not understood" << std::endl;
    }

    file.close();

    // Assign values to a struct
    // Parameters param; 
    param.lowT_out = lowT_outx;
    param.highT_out = highT_outx;
    param.T_step = T_stepx;
    param.lowMU_out = lowMU_outx;
    param.highMU_out = highMU_outx;
    param.muB_step = muB_stepx;
    param.muBC = muBCx;
    param.rho = rhox;
    param.w = ww;
    // param.angle1 = (atan((2*Kappa2n(TC(muBCx))*muBCx)/(TC(muBCx)*dTprimedT(TC(muBCx),muBCx))));
    param.angle1 = ((2*Kappa2n(TC(muBCx))*muBCx)/(TC(muBCx)*dTprimedT(TC(muBCx),muBCx)));
    // param.angle12 = param.angle1;
    param.angle12 = PI/180*alpha12x ;//- PI/180*alpha12x;
    param.angle2 = param.angle1 - param.angle12;
    param.angle12p = atan(tan(param.angle1) - tan(param.angle2));
    param.wp = param.w*(1/cos(param.angle1))*sqrt(pow(cos(param.angle1)*cos(param.angle2),2)+ pow(sin(param.angle12),2));
    param.rhop = param.rho*(pow(cos(param.angle1),2)/(sqrt(pow(cos(param.angle1)*cos(param.angle2),2)+ pow(sin(param.angle12),2))));

    // update the struct for global memory access


end_param:
    std::cerr << "Finished importing and setting parameters.\n";
} 


int get_options(int argc, char *argv[]) {
    int o;
    extern int optind;

    while ((o = getopt(argc, argv, "r:R:vpLSh")) != -1) {
        switch (o) {
            
        }
    }


    return 0;
}


// Parameters param; // global memory access using extern
// double muBC = param.muBC;
// double rho = param.rho;
// double w = param.w;
// double angle1 = param.angle1;
// double angle12 = param.angle12;
// double angle2 = param.angle2;

// double angle12p= atan(tan(angle1) - tan(angle2));
// double wp = w*(1/cos(angle1))*sqrt(pow(cos(angle1)*cos(angle2),2)+ pow(sin(angle12),2)); 
// double rhop = rho*(pow(cos(angle1),2)/(sqrt(pow(cos(angle1)*cos(angle2),2)+ pow(sin(angle12),2)))); 


// void print(){
//     double muBC = param.muBC;
//     double rho = param.rho;
//     double w = param.w;
//     double angle1 = param.angle1;
//     double angle12 = param.angle12;
//     double angle2 = param.angle2;
//     cout << "Parameters Chosen New:\n";
//     cout << left << setw(25) << "    muBC" << "= " << muBC << " MeV\n";
//     cout << left << setw(25) << "    TC" << "= " << TC(muBC) << " MeV\n";

//     cout << left << setw(25) << "    w" << "= " << w << '\n';
//     cout << left << setw(25) << "    rho" << "= " << rho << '\n';

// }

//******************************************************************************************************************//
//****************** g(theta) and derivatives *************************************************************************************//
double g(double theta){
    return c0 + c1*(1-pow(theta,2)) + c2*pow(1-pow(theta,2),2) + c3*pow(1-pow(theta,2),3);
}
double dgdth(double theta){
    return -2*theta *(c1 + (-1 + pow(theta,2))*(-2 *c2 + 3*c3 *(-1 + pow(theta,2))));
}
double d2gdth2(double theta){
    return -2*(c1 + c2*(2 - 6*pow(theta,2)) + 3*c3*(1 - 6*pow(theta,2) + 5*pow(theta,4)));
}

//******************************************************************************************************************//
//****************** htilda(theta) and derivatives *******************************************************************************//

double htild(double theta){
    return theta*(1 + aa*pow(theta,2) + bb*pow(theta,4));
}
double dhtilddth(double theta){
    return 1 + 3*aa*pow(theta,2) + 5*bb*pow(theta,4);
}
double d2htilddth2(double theta){
    return 6*aa*theta + 20*bb*pow(theta,3);
}
//******************************************************************************************************************//

//****************** Pressure in (R, theta)  ***********************************************************************//
double Press(double R,double theta){
    return h0*M0*pow(R,2-alpha)*(theta*htild(theta)-g(theta));
}

//****************** 1th derivative  ********************************************************************************//
double dPdR(double R, double theta){
    return (2 - alpha)*h0*M0*(-g(theta) + theta*htild(theta))*pow(R,1 - alpha);
}
double dPdth(double R, double theta){
    return h0*M0*(-dgdth(theta) + theta*dhtilddth(theta) + htild(theta))*pow(R,2 - alpha);
}
//******************************************************************************************************************//
//****************** 2th derivative  ********************************************************************************//
double d2PdR2(double R, double theta){
    return (1 - alpha)*(2 - alpha)*h0*M0*(-g(theta) + theta*htild(theta))*pow(R,-alpha);
}
double d2PdRdth(double R, double theta){
    return (2 - alpha)*h0*M0*(-dgdth(theta) + theta*dhtilddth(theta) + htild(theta))*pow(R,1 - alpha);
}
double d2Pdth2(double R, double theta){
    return h0*M0*(-d2gdth2(theta) + theta*d2htilddth2(theta) + 2*dhtilddth(theta))*pow(R,2 - alpha);
}
//******************************************************************************************************************//

//****************** 1th derivative  ********************************************************************************//
double dRdrr(double R, double theta){
    return dhtilddth(theta)*pow(dhtilddth(theta) + 2*bede*theta*htild(theta) - dhtilddth(theta)*pow(theta,2),-1);
}  
double dRdh(double R, double theta){
    return 2*theta*pow(h0,-1)*pow(R,1 - bede)*pow(dhtilddth(theta) + 2*bede*theta*htild(theta) - dhtilddth(theta)*pow(theta,2),-1);
}
double dthdrr(double R, double theta){
    return bede*htild(theta)*pow(R,-1)*pow(-2*bede*theta*htild(theta) + dhtilddth(theta)*(-1 + pow(theta,2)),-1);
}  
double dthdh(double R, double theta){
    return pow(h0,-1)*pow(R,-(bede))*(-1 + pow(theta,2))*pow(-2*bede*theta*htild(theta) + dhtilddth(theta)*(-1 + pow(theta,2)),-1);
}  
//******************************************************************************************************************//

//****************** 2th derivative  *******************************************************************************//
double d2Rdrr2(double R, double theta){
    return -(dhtilddth(theta)*(-2*theta*dhtilddth(theta)*dthdrr(R,theta) + 
             2*beta*delta*theta*dhtilddth(theta)*dthdrr(R,theta) + 
              2*beta*delta*dthdrr(R,theta)*htild(theta) + 
              d2htilddth2(theta)*dthdrr(R,theta)*(1 - pow(theta,2)))*
            pow(2*beta*delta*theta*htild(theta) + 
          dhtilddth(theta)*(1 - pow(theta,2)),-2)) + 
             d2htilddth2(theta)*dthdrr(R,theta)*
             pow(2*beta*delta*theta*htild(theta) + 
           dhtilddth(theta)*(1 - pow(theta,2)),-1);
}
//******************************************************************************************************************//
double d2Rdh2(double R, double theta){
    return -2*theta*pow(h0,-1)*pow(R,1 - beta*delta)*
          (d2htilddth2(theta)*dthdh(R,theta) - 
          2*theta*dhtilddth(theta)*dthdh(R,theta) + 
            2*beta*delta*theta*dhtilddth(theta)*dthdh(R,theta) + 
             2*beta*delta*dthdh(R,theta)*htild(theta) - 
              d2htilddth2(theta)*dthdh(R,theta)*pow(theta,2))*
           pow(dhtilddth(theta) + 2*beta*delta*theta*htild(theta) - 
             dhtilddth(theta)*pow(theta,2),-2) + 
              2*(1 - beta*delta)*theta*dRdh(R,theta)*pow(h0,-1)*
             pow(R,-(beta*delta))*pow(dhtilddth(theta) + 
              2*beta*delta*theta*htild(theta) - dhtilddth(theta)*pow(theta,2),
          -1) + 2*dthdh(R,theta)*pow(h0,-1)*pow(R,1 - beta*delta)*
             pow(dhtilddth(theta) + 2*beta*delta*theta*htild(theta) - 
          dhtilddth(theta)*pow(theta,2),-1);
} 
//******************************************************************************************************************//
double d2thdrr2(double R, double theta){
    return beta*delta*htild(theta)*pow(R,-1)*
              (-2*theta*dhtilddth(theta)*dthdrr(R,theta) + 
            2*beta*delta*theta*dhtilddth(theta)*dthdrr(R,theta) + 
            2*beta*delta*dthdrr(R,theta)*htild(theta) + 
              d2htilddth2(theta)*dthdrr(R,theta)*(1 - pow(theta,2)))*
             pow(2*beta*delta*theta*htild(theta) + 
              dhtilddth(theta)*(1 - pow(theta,2)),-2) + 
             beta*delta*dRdrr(R,theta)*htild(theta)*pow(R,-2)*
           pow(2*beta*delta*theta*htild(theta) + 
             dhtilddth(theta)*(1 - pow(theta,2)),-1) - 
             beta*delta*dhtilddth(theta)*dthdrr(R,theta)*pow(R,-1)*
           pow(2*beta*delta*theta*htild(theta) + 
         dhtilddth(theta)*(1 - pow(theta,2)),-1);
} 
//******************************************************************************************************************//
double d2thdh2(double R, double theta){
    return   -(pow(h0,-1)*pow(R,-(beta*delta))*
              (-2*theta*dhtilddth(theta)*dthdh(R,theta) + 
              2*beta*delta*theta*dhtilddth(theta)*dthdh(R,theta) + 
               2*beta*delta*dthdh(R,theta)*htild(theta) + 
               d2htilddth2(theta)*dthdh(R,theta)*(1 - pow(theta,2)))*
              (1 - pow(theta,2))*pow(2*beta*delta*theta*htild(theta) + 
            dhtilddth(theta)*(1 - pow(theta,2)),-2)) - 
            2*theta*dthdh(R,theta)*pow(h0,-1)*pow(R,-(beta*delta))*
             pow(2*beta*delta*theta*htild(theta) + 
         dhtilddth(theta)*(1 - pow(theta,2)),-1) - 
         beta*delta*dRdh(R,theta)*pow(h0,-1)*pow(R,-1 - beta*delta)*
            (1 - pow(theta,2))*pow(2*beta*delta*theta*htild(theta) + 
             dhtilddth(theta)*(1 - pow(theta,2)),-1);
} 
//******************************************************************************************************************//
double d2Rdrrdh(double R, double theta){
    return -(dhtilddth(theta)*(-2*theta*dhtilddth(theta)*dthdh(R,theta) + 
              2*beta*delta*theta*dhtilddth(theta)*dthdh(R,theta) + 
                 2*beta*delta*dthdh(R,theta)*htild(theta) + 
              d2htilddth2(theta)*dthdh(R,theta)*(1 - pow(theta,2)))*
             pow(2*beta*delta*theta*htild(theta) + 
               dhtilddth(theta)*(1 - pow(theta,2)),-2)) + 
             d2htilddth2(theta)*dthdh(R,theta)*
             pow(2*beta*delta*theta*htild(theta) + 
              dhtilddth(theta)*(1 - pow(theta,2)),-1);
} 
//******************************************************************************************************************//
double d2thdrrdh(double R, double theta){
    return beta*delta*pow(R,-2)*(htild(theta)*
             (2*beta*delta*(theta*dRdh(R,theta) + R*dthdh(R,theta))*
              htild(theta) - R*d2htilddth2(theta)*dthdh(R,theta)*
               (-1 + pow(theta,2))) + 
             dhtilddth(theta)*(htild(theta)*
              (-2*R*theta*dthdh(R,theta) - dRdh(R,theta)*(-1 + pow(theta,2)))
                + R*dhtilddth(theta)*dthdh(R,theta)*(-1 + pow(theta,2))))*
             pow(-2*beta*delta*theta*htild(theta) + 
              dhtilddth(theta)*(-1 + pow(theta,2)),-2);
} 
//******************************************************************************************************************//

//****************** 1th derivative  ******************************************************************************//
double dPdrr(double R, double theta){
    return dPdR(R,theta)*dRdrr(R,theta) + dPdth(R,theta)*dthdrr(R,theta);
}  
//******************************************************************************************************************//
double dPdh(double R, double theta){
    return dPdR(R,theta)*dRdh(R,theta) + dPdth(R,theta)*dthdh(R,theta);
} 
//****************** 2th derivative  ******************************************************************************//

double d2Pdrr2(double R, double theta){
    return d2Rdrr2(R,theta)*dPdR(R,theta) + d2thdrr2(R,theta)*dPdth(R,theta) + 
   2*d2PdRdth(R,theta)*dRdrr(R,theta)*dthdrr(R,theta) + d2PdR2(R,theta)*pow(dRdrr(R,theta),2) + 
   d2Pdth2(R,theta)*pow(dthdrr(R,theta),2);
}
//******************************************************************************************************************//
double d2Pdh2(double R, double theta){
    return d2Rdh2(R,theta)*dPdR(R,theta) + d2thdh2(R,theta)*dPdth(R,theta) + 
   2*d2PdRdth(R,theta)*dRdh(R,theta)*dthdh(R,theta) + d2PdR2(R,theta)*pow(dRdh(R,theta),2) + 
   d2Pdth2(R,theta)*pow(dthdh(R,theta),2);
}
//******************************************************************************************************************//
double d2Pdrrdh(double R, double theta){
    return d2Rdrrdh(R,theta)*dPdR(R,theta) + d2thdrrdh(R,theta)*dPdth(R,theta) + 
   d2PdR2(R,theta)*dRdh(R,theta)*dRdrr(R,theta) + d2PdRdth(R,theta)*dRdrr(R,theta)*dthdh(R,theta) + 
   d2PdRdth(R,theta)*dRdh(R,theta)*dthdrr(R,theta) + d2Pdth2(R,theta)*dthdh(R,theta)*dthdrr(R,theta);
}
//******************************************************************************************************************//




double hh(double T, double muB) { // h as a function of T and muB
  return -(Tprime(T,muB) - T0) / (param.wp *TC(param.muBC)* dTprimedT(TC(param.muBC),param.muBC) * sin(param.angle12p));
}
double rr(double T, double muB) { // r as a function of T and muB
  return -((pow(muB, 2) - pow(param.muBC, 2)) /
               ( 2*param.wp * param.muBC*TC(param.muBC)) +
           (hh(T, muB)* cos(param.angle12p)))/param.rhop ;
}

double dhhdmuB(double T, double muB) {
  return -1 /  (param.wp *TC(param.muBC)* dTprimedT(TC(param.muBC),param.muBC) * sin(param.angle12p)) * dTprimedmuB(T, muB);
}
double d2hhdTdmuB(double T, double muB) {
   return -1 /  (param.wp *TC(param.muBC)* dTprimedT(TC(param.muBC),param.muBC) * sin(param.angle12p)) * d2TprimedTdmuB(T, muB);
 }
double dhhdT(double T, double muB) {
  return -1 /  (param.wp *TC(param.muBC)* dTprimedT(TC(param.muBC),param.muBC) * sin(param.angle12p)) * dTprimedT(T, muB);
}
double d2hhdmuB2(double T, double muB) {
  return -1 / (param.wp *TC(param.muBC)* dTprimedT(TC(param.muBC),param.muBC) * sin(param.angle12p)) * d2TprimedmuB2(T, muB);
}
double drrdmuB(double T, double muB) {
  return -((2 * muB) / ( 2*param.wp * param.muBC*TC(param.muBC)) +
           dhhdmuB(T, muB)* cos(param.angle12p))/param.rhop ;
}
double d2rrdTdmuB(double T, double muB) {
  return - d2hhdTdmuB(T, muB)* cos(param.angle12p)/param.rhop ;
}
double drrdT(double T, double muB) {
  return -(dhhdT(T, muB)* cos(param.angle12p))/param.rhop ;
}
double d2rrdmuB2(double T, double muB) {
  return -(2 / ( 2*param.wp * param.muBC*TC(param.muBC)) +
           d2hhdmuB2(T, muB) * cos(param.angle12p))/param.rhop ;
}



//******************************************************************************************************************//

//****************** Mapping   (R,theta) to (T, muB)
//************************************************************************//


//Generating cordinates for the mapping (R,theta) to (T, muB)

double RRxx(double T, double muB) { return get<0>(FindRoot(rr(T, muB), hh(T,muB))); }
double Thetax(double T, double muB) { return get<1>(FindRoot(rr(T,muB), hh(T, muB))); }


//******************************************************************************************************************//
//******************************************************************************************************************//

//First derivative with respect to muB
double dPdmuB(double T,double muB, double R, double theta){
    return dhhdmuB(T,muB)*dPdh(R,theta) + dPdrr(R,theta)*drrdmuB(T,muB);
}

//second derivative with respect to muB
double d2PdmuB2(double T,double muB, double R, double theta){
    return d2Pdh2(R,theta)*pow(dhhdmuB(T,muB),2) + d2hhdmuB2(T,muB)*dPdh(R,theta) + 
   d2rrdmuB2(T,muB)*dPdrr(R,theta) + 2*d2Pdrrdh(R,theta)*dhhdmuB(T,muB)*drrdmuB(T,muB) + 
   d2Pdrr2(R,theta)*pow(drrdmuB(T,muB),2);
}
double d2PdTdmuB(double T,double muB, double R, double theta){
  return d2Pdh2(R,theta)*dhhdT(T,muB)*dhhdmuB(T,muB) + d2Pdrrdh(R,theta)*drrdT(T,muB)*dhhdmuB(T,muB) + dPdh(R,theta)*d2hhdTdmuB(T,muB) +d2Pdrrdh(R,theta)*dhhdT(T,muB)*drrdmuB(T,muB)+ d2Pdrr2(R,theta)*drrdT(T,muB)*drrdmuB(T,muB) + dPdrr(R,theta)*d2rrdTdmuB(T,muB) ;
  
}

// Analytical derivative
double Pcrit(double T, double muB){
    return Press(RRxx( T,muB),Thetax( T,muB));
}
double BardIsingA(double T, double muB){
    return T*dPdmuB(T,muB,RRxx( T,muB),Thetax( T,muB));
}
double dBardIsingAdT(double T, double muB){
    return BardIsingA(T,muB)/T + T*d2PdTdmuB(T,muB,RRxx( T,muB),Thetax( T,muB));
}
double dBardIsingdmuBA(double T, double muB){
    return T*T*d2PdmuB2(T,muB,RRxx( T,muB),Thetax( T,muB));
}
double d3BardIsingdmuB3N(double T, double muB){
  double diffx=0.01;
    return T*T*(-dBardIsingdmuBA(T,muB - 2*diffx) + 16*dBardIsingdmuBA(T,muB - diffx) -30*dBardIsingdmuBA(T,muB) + 16*dBardIsingdmuBA(T,muB + diffx) - dBardIsingdmuBA(T,muB + 2*diffx))/(12*diffx*diffx);
}
double ddBardIsingdmuBdT(double T, double muB){
    return (dBardIsingdmuBA(T- 2*diff,muB) - 8*dBardIsingdmuBA(T- diff,muB) + 8*dBardIsingdmuBA(T+ diff,muB) - dBardIsingdmuBA(T + 2*diff,muB))/(12*diff);
}
double dd3BardIsingdmuB3dT(double T, double muB){
    return (d3BardIsingdmuB3N(T- 2*diff,muB) - 8*d3BardIsingdmuB3N(T- diff,muB) + 8*d3BardIsingdmuB3N(T+ diff,muB) - d3BardIsingdmuB3N(T + 2*diff,muB))/(12*diff);
}
// Numerical derivatives
//******************************************************************************************************************//
//******************************************************************************************************************//
double PresIsingN(double T, double muB){
    return (Press(RRxx( T,  muB),Thetax(T, muB))) ;
}
double BardIsingN(double T, double muB){
  return T*(-PresIsingN(T,muB - 3*diff) + 9*PresIsingN(T,muB - 2*diff) - 45*PresIsingN(T,muB - diff) + 45*PresIsingN(T,muB + diff) - 9*PresIsingN(T,muB + 2*diff) + PresIsingN(T,muB + 3*diff))/(12*diff);
}
double dBardIsingdmuBN(double T, double muB){
    return T*(BardIsingN(T,muB - 2*diff) - 8*BardIsingN(T,muB - diff) + 8*BardIsingN(T,muB + diff) - BardIsingN(T,muB + 2*diff))/(12*diff);
}



 //Defining Tprime Full and Tprime Crit
//******************************************************************************************************************//
// Tprime Critical
double Tcrit(double T, double muB){
   return (muB == 0) ? pow(dchi2dT(T0),-1)*dBardIsingdmuBA(T, 0) : 1/dchi2dT(T0)*BardIsingA(T, muB)/(muB/T);
}
double TyTcrit(double T, double muB){
  // return  1/dchi2dT(T0)*(dBardIsingdmuBA(T,0) + (1/6)*d3BardIsingdmuB3N(T,0)*pow(muB/T,2));
  return pow(dchi2dT(T0),-1)*(dBardIsingdmuBA(T,0) + (d3BardIsingdmuB3N(T,0)*pow(muB/T,2))/6.);
}
// Tprime Full, Its contains both Critical and non- Critical contributions
double Tfull(double T, double muB){
    return Tprime(T,muB) + Tcrit(T, muB) - TyTcrit(T,muB);
}

double dTfulldmuB(double T, double muB){
  return 2*Kappa2n(T)*muB + 1/dchi2dT(T0)*((1/(muB/T) )*dBardIsingdmuBA(T,muB) - BardIsingA(T, muB)/(pow(muB/T,2))+(1/3)*d3BardIsingdmuB3N(T,0)*pow(muB/T,1));
}



double dTfulldTn(double T, double muB){
  return (Tfull(T - 2*diff,muB) - 8*Tfull(T - 1*diff,muB) + 8*Tfull(T + 1*diff, muB) - Tfull(T + 2*diff,muB))/(12*diff);
}
double dTfulldT(double T, double muB) {
    return (muB==0) ? 1 + 1/ (dchi0dT(T0)) * (1/T*dBardIsingdmuBA(T,0) + ddBardIsingdmuBdT(T, 0) ) : dTprimedT(T,muB) +  1/ (dchi0dT(T0)) * ((T/muB)*dBardIsingAdT(T, muB)  + BardIsingA(T, muB) / muB + ddBardIsingdmuBdT(T, 0) + (1 / 6) * dd3BardIsingdmuB3dT(T, 0) * pow(muB / T, 2) - (1 / 3) * d3BardIsingdmuB3N(T, 0) * pow(muB, 2) / pow(T, 3));
}


//******************************************************************************************************************//
////////////////////////////////////// THERMODYNAMIC OBSERVABLES/////////////////////////////////
//******************************************************************************************************************//

// Full baryon Density
//******************************************************************************************************************//


double nBfull(double T, double muB){
  return muB/T*chi2(Tfull(T,muB))*pow(T,3);
  
}
// Baryon Number susceptability
//******************************************************************************************************************//
double chi2full(double T, double muB){
  return muB == 0 ? chi2(T)*pow(T,2) : (chi2(Tfull(T,muB)) + (muB/T)*dchi2dT(Tfull(T,muB))*dTfulldmuB(T,muB))*pow(T,2); 
}

double dnBfulldT(double T, double muB){
  return 2*nBfull(T,muB)/T + (muB/T)*dchi2dT(Tfull(T,muB))*dTfulldT(T,muB)*pow(T,3); 
}
double dnBfulldTn(double T, double muB){
   
  return (nBfull(T- 2*diff,muB) - 8*nBfull(T- diff,muB ) + 8*nBfull(T+ diff,muB ) - nBfull(T+ 2*diff,muB ))/(12*diff);
   }


//******************************************************************************************************************//
//******************************************************************************************************************//
// Pressure
//******************************************************************************************************************//




double Pressfull(double T, double muB)
{
  
    double integral = 0.0;

    // Gauss Quadrature points and weights
    std::vector<double> gaussPoints = {-0.90618, -0.538469, 0.0, 0.538469, 0.90618};
    std::vector<double> gaussWeights = {0.236927, 0.478629, 0.568889, 0.478629, 0.236927};

    double y, step =  0.1; // Integration step size
    //double y, step =  (fabs(T - TC(muB)) < 5 )? 1 : 1 ;

    std::vector<double> values;
    values.push_back(nBfull(T, 0)); // Store initial value

    int gaussIndex = 0;
    for (y = step; y <= muB; y += step)
    {
        values.push_back(nBfull(T, y));
        integral += (values.back() + values[values.size() - 2]) * gaussWeights[gaussIndex] / 2.0;
        gaussIndex = (gaussIndex + 1) % gaussPoints.size();
    }
    return  pow(T, 4) * chi0(T) + integral; 
}

double Entropyn(double T, double muB){
  
  return (Pressfull(T- 2*diff,muB) - 8*Pressfull(T- diff,muB ) + 8*Pressfull(T+ diff,muB ) - Pressfull(T+ 2*diff,muB ))/(12*diff);
   }



double Entropy(double T, double muB)
{
    double integral = 0.0;

    // Gauss Quadrature points and weights
    std::vector<double> gaussPoints = {-0.90618, -0.538469, 0.0, 0.538469, 0.90618};
    std::vector<double> gaussWeights = {0.236927, 0.478629, 0.568889, 0.478629, 0.236927};

    //double y, step = 1; // Integration step size
    double y, step =  0.1; // Integration step size

    std::vector<double> values;
    values.push_back(0); // Store initial value

    int gaussIndex = 0;
    for (y = step; y <= muB; y += step)
    {
        values.push_back(dnBfulldT(T, y));
        integral += (values.back() + values[values.size() - 2]) * gaussWeights[gaussIndex] / 2.0;
        gaussIndex = (gaussIndex + 1) % gaussPoints.size();
    }

    return 4*pow(T,3)*chi0(T) + pow(T,4)*dchi0dT(T) + integral;
   
}


double dEntrpydT(double T, double muB){
 
  return (Entropy(T- 2*diff,muB) - 8*Entropy(T- diff,muB) + 8*Entropy(T+ diff,muB) - Entropy(T+ 2*diff,muB))/(12*diff);
}

double Energy(double T, double muB){
  return  - Pressfull(T,muB) + T*Entropy(T,muB) + muB*nBfull(T,muB);
}

// }
double Cv(double T, double muB){
  return T*(dEntrpydT(T,muB)*chi2full(T,muB) - pow(dnBfulldT(T,muB),2))/chi2full(T,muB) ;

}
double Cs2(double T, double muB){
  return (pow(nBfull(T,muB),2)*dEntrpydT(T,muB) - 2*Entropy(T,muB)*nBfull(T,muB)*dnBfulldT(T,muB) + pow(Entropy(T,muB),2)*chi2full(T,muB))/(((Energy(T,muB) + Pressfull(T,muB)))*(dEntrpydT(T,muB)*chi2full(T,muB) - pow(dnBfulldT(T,muB),2) ));
}


