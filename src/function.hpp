#ifndef FUNCTION_HPP
#define FUNCTION_HPP



#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <tuple>
#include <vector>
#include <numeric>
#include <algorithm>

#define aa -0.76201
#define bb 0.00804
#define bede 1.5648
#define beta 0.326
#define delta 4.80
#define alpha 0.11
#define M0 0.605
#define h0 0.394
#define c0 0.0424369
#define c1 0.321329
#define c2 -1.20375
#define c3 -0.00126032
#define PI 3.14159265358979323846

#define diff 0.01

using namespace std;

#define T0 158.0




void read_param(char *param_file);
int get_options(int argc, char *argv[]);



double chi0(double T);

double dchi0dT(double T);

double d2chi0dT2(double T);

double chi2(double T);

double chi4(double T);

double dchi2dT(double T);

double Kappa2n(double T) ;

double dKappa2ndT(double T);

double Tprime(double T, double muB);

double dTprimedmuB(double T, double muB);

double d2TprimedTdmuB(double T, double muB);

double dTprimedT(double T, double muB);

double d2TprimedTdmuB(double T, double muB);

double d2TprimedmuB2(double T, double muB);


double TC (double muB);

double muBCC(double T);


struct Parameters {
    double muBC;
    double rho;
    double w;
    double angle1;
    double angle12;
    double angle2;
    double angle12p;
    double wp;
    double rhop;
    double lowT_out;
    double highT_out;
    double T_step;
    double lowMU_out;
    double highMU_out;
    double muB_step;
};

// extern Parameters param;


// #define angle2  (angle1 - angle12)

// #define angle12p atan(tan(angle1) - tan(angle2));
// #define wp  w*(1/cos(angle1))*sqrt(pow(cos(angle1)*cos(angle2),2)+ pow(sin(angle12p),2)); 
// #define rhop  rho*(pow(cos(angle1),2)/(sqrt(pow(cos(angle1)*cos(angle2),2)+ pow(sin(angle12p),2)))); 

const double EPS = 1e-5; // tolerance




/** Function Declarations **/



/*** Variable  Declarations ***/




double g(double theta);

double dgdth(double theta);

double d2gdth2(double theta);

double htild(double theta);

double dhtilddth(double theta);

double d2htilddth2(double theta);

double Press(double R,double theta);

double dPdR(double R, double theta);

double dPdth(double R, double theta);

double d2PdR2(double R, double theta);

double d2PdRdth(double R, double theta);

double d2Pdth2(double R, double theta);

double dRdrr(double R, double theta);

double dRdh(double R, double theta);

double dthdrr(double R, double theta);

double dthdh(double R, double theta);

double d2Rdrr2(double R, double theta);

double d2Rdh2(double R, double theta);

double d2thdrr2(double R, double theta);

double d2thdh2(double R, double theta);

double d2Rdrrdh(double R, double theta);

double d2thdrrdh(double R, double theta);

double dPdrr(double R, double theta);

double dPdh(double R, double theta);

double d2Pdrr2(double R, double theta);

double d2Pdh2(double R, double theta);

double d2Pdrrdh(double R, double theta);


double hh(double T, double muB);

double rr(double T, double muB);

double dhhdmuB(double T, double muB) ;

double d2hhdTdmuB(double T, double muB);

double dhhdT(double T, double muB);

double d2hhdmuB2(double T, double muB);

double drrdmuB(double T, double muB);

double d2rrdTdmuB(double T, double muB);

double drrdT(double T, double muB);

double d2rrdmuB2(double T, double muB);

double h(double R, double theta);

double r(double R, double theta);

tuple<double, double> FindRoot(double T, double muB);

double RRxx(double T, double muB) ;

double Thetax(double T, double muB) ;

double dPdmuB(double T,double muB, double R, double theta);

double d2PdmuB2(double T,double muB, double R, double theta);

double d2PdTdmuB(double T,double muB, double R, double theta);

double Pcrit(double T, double muB);

double BardIsingA(double T, double muB);

double dBardIsingAdT(double T, double muB);

double dBardIsingdmuBA(double T, double muB);

double d3BardIsingdmuB3N(double T, double muB);

double ddBardIsingdmuBdT(double T, double muB);

double dd3BardIsingdmuB3dT(double T, double muB);

double PresIsingN(double T, double muB);

double BardIsingN(double T, double muB);

double dBardIsingdmuBN(double T, double muB);

double Tcrit(double T, double muB);

double TyTcrit(double T, double muB);

double Tfull(double T, double muB);

double dTfulldmuB(double T, double muB);

double dTfulldTn(double T, double muB);

double dTfulldT(double T, double muB);

double nBfull(double T, double muB);

double chi2full(double T, double muB);

double dnBfulldT(double T, double muB);

double dnBfulldTn(double T, double muB);

double Pressfull(double T, double muB);

double Entropyn(double T, double muB);

double Entropy(double T, double muB);

double dEntrpydT(double T,double muB);

double Energy(double T, double muB);

double Cv(double T, double muB);

double Cs2(double T, double muB);

double radiansToDegrees(double radians);

void storevector(vector<double>&vec,int row, double rowstep, double (*function)(double)) ;

void storeMatrix(vector<vector<double>>& matrix, int row, double rowstep, int col, double colstep, double (*function)(double, double));

vector<vector<double>> integrateMatrix(const vector<vector<double>>& matrix, int row, double rowstep, int col, double colstep, bool direction);


vector<vector<double>> deriv_matrix(const vector<vector<double>>& matrix, int rowSize, double rowStep, int colSize, double colStep, int derivativeDirection);

#endif
