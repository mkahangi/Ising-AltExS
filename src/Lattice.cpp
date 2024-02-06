#include "function.hpp"
//Defining the functions for the Lattice and the derivatives (Taylor coefficient)


// extern Parameters param; // global memory access using extern
// double muBC = param.muBC;
// double rho = param.rho;
// double w = param.w;
// double angle1 = param.angle1;
// double angle12 = param.angle12;
// double angle2 = param.angle2;


double chi0(double T) {
  double aa00 = 7.53891;
  double aa10 = -6.18858;
  double aa20 = -5.37961;
  double aa30 = 7.08750;
  double aa40 = -0.977970;
  double aa50 = 0.0302636;
  double bb00 = 2.24530;
  double bb10 = -6.02568;
  double bb20 = 15.3737;
  double bb30 = -19.6331;
  double bb40 = 10.2400;
  double bb50 = 0.799479;
  return (aa00 + aa10 * pow(154 / T, 1) + aa20 * pow(154 / T, 2) +
          aa30 * pow(154 / T, 3) + aa40 * pow(154 / T, 4) +
          aa50 * pow(154 / T, 5)) /
         (bb00 + bb10 * pow(154 / T, 1) + bb20 * pow(154 / T, 2) +
          bb30 * pow(154 / T, 3) + bb40 * pow(154 / T, 4) +
          bb50 * pow(154 / T, 5));
}
double dchi0dT(double T) {
  double aa00 = 7.53891;
  double aa10 = -6.18858;
  double aa20 = -5.37961;
  double aa30 = 7.08750;
  double aa40 = -0.977970;
  double aa50 = 0.0302636;
  double bb00 = 2.24530;
  double bb10 = -6.02568;
  double bb20 = 15.3737;
  double bb30 = -19.6331;
  double bb40 = 10.2400;
  double bb50 = 0.799479;
return -((((-433085465120*bb50)/pow(T,6) - (2249794624*bb40)/pow(T,5) - 
          (10956792*bb30)/pow(T,4) - (47432*bb20)/pow(T,3) - 
          (154*bb10)/pow(T,2))*
        (aa00 + (86617093024*aa50)/pow(T,5) + 
          (562448656*aa40)/pow(T,4) + (3652264*aa30)/pow(T,3) + 
          (23716*aa20)/pow(T,2) + (154*aa10)/T))/
      pow(bb00 + (86617093024*bb50)/pow(T,5) + 
        (562448656*bb40)/pow(T,4) + (3652264*bb30)/pow(T,3) + 
        (23716*bb20)/pow(T,2) + (154*bb10)/T,2)) + 
   ((-433085465120*aa50)/pow(T,6) - (2249794624*aa40)/pow(T,5) - 
      (10956792*aa30)/pow(T,4) - (47432*aa20)/pow(T,3) - 
      (154*aa10)/pow(T,2))/
    (bb00 + (86617093024*bb50)/pow(T,5) + 
      (562448656*bb40)/pow(T,4) + (3652264*bb30)/pow(T,3) + 
      (23716*bb20)/pow(T,2) + (154*bb10)/T);
}

double d2chi0dT2(double T) {
  double aa00 = 7.53891;
  double aa10 = -6.18858;
  double aa20 = -5.37961;
  double aa30 = 7.08750;
  double aa40 = -0.977970;
  double aa50 = 0.0302636;
  double bb00 = 2.24530;
  double bb10 = -6.02568;
  double bb20 = 15.3737;
  double bb30 = -19.6331;
  double bb40 = 10.2400;
  double bb50 = 0.799479;
return ((2*pow((-433085465120*bb50)/pow(T,6) - (2249794624*bb40)/pow(T,5) - (10956792*bb30)/pow(T,4) - (47432*bb20)/pow(T,3) - (154*bb10)/pow(T,2),2))/
       pow(bb00 + (86617093024*bb50)/pow(T,5) + (562448656*bb40)/pow(T,4) + (3652264*bb30)/pow(T,3) + (23716*bb20)/pow(T,2) + (154*bb10)/T,3) - 
      ((2598512790720*bb50)/pow(T,7) + (11248973120*bb40)/pow(T,6) + (43827168*bb30)/pow(T,5) + (142296*bb20)/pow(T,4) + (308*bb10)/pow(T,3))/
       pow(bb00 + (86617093024*bb50)/pow(T,5) + (562448656*bb40)/pow(T,4) + (3652264*bb30)/pow(T,3) + (23716*bb20)/pow(T,2) + (154*bb10)/T,2))*
    (aa00 + (86617093024*aa50)/pow(T,5) + (562448656*aa40)/pow(T,4) + (3652264*aa30)/pow(T,3) + (23716*aa20)/pow(T,2) + (154*aa10)/T) - 
   (2*((-433085465120*aa50)/pow(T,6) - (2249794624*aa40)/pow(T,5) - (10956792*aa30)/pow(T,4) - (47432*aa20)/pow(T,3) - (154*aa10)/pow(T,2))*
      ((-433085465120*bb50)/pow(T,6) - (2249794624*bb40)/pow(T,5) - (10956792*bb30)/pow(T,4) - (47432*bb20)/pow(T,3) - (154*bb10)/pow(T,2)))/
    pow(bb00 + (86617093024*bb50)/pow(T,5) + (562448656*bb40)/pow(T,4) + (3652264*bb30)/pow(T,3) + (23716*bb20)/pow(T,2) + (154*bb10)/T,2) + 
   ((2598512790720*aa50)/pow(T,7) + (11248973120*aa40)/pow(T,6) + (43827168*aa30)/pow(T,5) + (142296*aa20)/pow(T,4) + (308*aa10)/pow(T,3))/
    (bb00 + (86617093024*bb50)/pow(T,5) + (562448656*bb40)/pow(T,4) + (3652264*bb30)/pow(T,3) + (23716*bb20)/pow(T,2) + (154*bb10)/T);
}


double chi2(double T){
  double d1 = 0.7327203479284046;
  double d2 = 11.186232699501845;
  double d3 = 0.31757988857623776;
  double d4 = 0.19950943628460033;
  double d5 = 0.6886799470187275;
  double mp = 938.2720813;

  return pow((2*mp/200)/(PI*(T/200)),3/2)*exp(-(mp/200)/(T/200))/(1+pow((T/200)/d1,d2)) + d3*exp(-pow(d4,2)/pow((T/200),2)-pow(d5,4)/pow((T/200),4))/(1+ pow((T/200)/d1,-d2));
}






double chi4(double T) {
  double a0 = 0.0009270171006903285;
  double a1 = -0.0013607301311995574;
  double a2 = -0.0005190343146423107;
  double a3 = 0.0015368880694823762;
  double a4 = -0.0006434870283073376;
  double a5 = 0.000062480577701191;
  double b0 = -0.032898097117572744;
  double b1 = 0.5274850805193245;
  double b2 = -1.7320609548220323;
  double b3 = 2.376461108582519;
  double b4 = -1.5011718623101546;
  double b5 = 0.362232372641797;
  return (a0 + a1 * pow(154 / T, 1) + a2 * pow(154 / T, 2) +
          a3 * pow(154 / T, 3) + a4 * pow(154 / T, 4) + a5 * pow(154 / T, 5)) /
         (b0 + b1 * pow(154 / T, 1) + b2 * pow(154 / T, 2) +
          b3 * pow(154 / T, 3) + b4 * pow(154 / T, 4) + b5 * pow(154 / T, 5));
}
// double dchi2dT(double T) {
//   double h1 = -0.318451210804933;
//   double h2 = 0.4927954599817167;
//   double f3 = 0.14931911904005374;
//   double f4 = 6.6991991785479295;
//   double f5 = -5.110080379554422;
//   return exp(-h1/(T/200) - h2/pow(T/200,2))*f3*(pow(1/cosh(f4*(T/200)+f5),2)*f4/200) + (h1/(pow(T,2)/200) + h2*80000/pow(T,3))*exp(-h1/(T/200) - h2/pow(T/200,2))*f3*(1+tanh(f4*(T/200) + f5));
// }

double dchi2dT(double T){
  double d1 = 0.7327203479284046;
  double d2 = 11.186232699501845;
  double d3 = 0.31757988857623776;
  double d4 = 0.19950943628460033;
  double d5 = 0.6886799470187275;
  double mp = 938.2720813;

  return exp(-(pow(T,-4)*(40000*(40000*pow(d5,4) + pow(d4,2)*pow(T,2)) + mp*pow(T,3))))*
   pow(PI,-1.5)*pow(T,-5)*(d3*exp(mp*pow(T,-1))*pow(PI,1.5)*pow(T*pow(d1,-1),d2)*
      (d2*pow(200,d2)*pow(T,4) + 6400000000*pow(d5,4)*(pow(200,d2) + pow(T*pow(d1,-1),d2)) + 
        80000*pow(d4,2)*pow(T,2)*(pow(200,d2) + pow(T*pow(d1,-1),d2))) - 
     pow(2,0.5 + 3*d2)*pow(25,d2)*exp(40000*pow(T,-4)*(40000*pow(d5,4) + pow(d4,2)*pow(T,2)))*pow(T,3)*(2*d2*T*pow(T*pow(d1,-1),d2) - 2*mp*(pow(200,d2) + pow(T*pow(d1,-1),d2)) + 
        3*T*(pow(200,d2) + pow(T*pow(d1,-1),d2)))*pow(mp*pow(T,-1),1.5))*
   pow(pow(200,d2) + pow(T*pow(d1,-1),d2),-2);
}


//******************************************************************************************************************//

//Defining the functions for the Lattice and the derivatives (Alternative expansion coefficient)
double Kappa2n(double T){
   
    double a1x = 0.48292140725481037847765989016973055;
    double a2x = -0.6408168549933056781061051325390577;
    double a3x = 0.519862524527988047270578759013881583;
    double a4x = -0.77583169896916358584790610556541513;
    double a5x = 0.532593649712493080242282347568790899;
    double b0x = 11.89642700608466971492008198159742386;
    double b1x = -9.624855783606008220736908683766264107;
    double b2x = -7.5071562888344658307286070034536059506;
    double b3x = 10.7732408388688673612104462968093193582;
   
      return  ( a1x/pow(200/T,1) + a2x/pow(200/T,2) + a3x/pow(200/T,3) + a4x/pow(200/T,4) + a5x/pow(200/T,5))/(b0x + b1x/pow(200/T,1) + b2x/pow(200/T,2) + b3x/pow(200/T,3) );
}
double dKappa2ndT(double T){
    double a1x = 0.48292140725481037847765989016973055;
    double a2x = -0.6408168549933056781061051325390577;
    double a3x = 0.519862524527988047270578759013881583;
    double a4x = -0.77583169896916358584790610556541513;
    double a5x = 0.532593649712493080242282347568790899;
    double b0x = 11.89642700608466971492008198159742386;
    double b1x = -9.624855783606008220736908683766264107;
    double b2x = -7.5071562888344658307286070034536059506;
    double b3x = 10.7732408388688673612104462968093193582;
   
  return -((b1x/200. + (b2x*T)/20000. + (3*b3x*pow(T,2))/8.e6)*
      ((a1x*T)/200. + (a2x*pow(T,2))/40000. + (a3x*pow(T,3))/8.e6 + 
        (a4x*pow(T,4))/1.6e9 + (a5x*pow(T,5))/3.2e11)*
      pow(b0x + (b1x*T)/200. + (b2x*pow(T,2))/40000. + (b3x*pow(T,3))/8.e6,-2)) + 
   (a1x/200. + (a2x*T)/20000. + (3*a3x*pow(T,2))/8.e6 + (a4x*pow(T,3))/4.e8 + 
      (a5x*pow(T,4))/6.4e10)*pow(b0x + (b1x*T)/200. + (b2x*pow(T,2))/40000. + 
      (b3x*pow(T,3))/8.e6,-1);
}



double Tprime(double T, double muB) { // Tprime from Lattice
  return T * (1 + Kappa2n(T) * pow(muB / T, 2));
}
double dTprimedmuB(double T, double muB) {
  return 2 * Kappa2n(T) * pow(muB / T, 1);
}
double dTprimedT(double T, double muB) {
  return 1 + dKappa2ndT(T)*pow(muB,2)/T -   Kappa2n(T) * pow(muB / T, 2);
}
double d2TprimedTdmuB(double T, double muB) {
  return 2*dKappa2ndT(T)*pow(muB/T,1) -  2*Kappa2n(T) * muB*pow(1 / T, 2);
}
double d2TprimedmuB2(double T, double muB) {
  return 2 * Kappa2n(T) / T;
}


double TC(double muB) {
    double TCm = 5.0; // Initial guess for TC
    // Broyden's method iteration
    double dif = 1.0; // Difference between consecutive iterations
    while (std::abs(dif) > EPS) {
        double Tprime_TCm = Tprime(TCm, muB);
        double J = (Tprime(TCm + EPS, muB) - Tprime_TCm) / EPS; // Numerical approximation of the Jacobian
        
        // Update TC using Broyden's method
        TCm -= (Tprime_TCm - T0) / J;
        
        dif = Tprime(TCm, muB) - T0;
    }
    return TCm;
}
double muBCC(double T) {
    double muBm = 400.0; // Initial guess for TC
    // Broyden's method iteration
    double dif = 1.0; // Difference between consecutive iterations
    while (std::abs(dif) > EPS) {
        double Tprime_muBm = Tprime(T, muBm);
        double J = (Tprime(T , muBm + EPS) - Tprime_muBm) / EPS; // Numerical approximation of the Jacobian
        
        // Update TC using Broyden's method
        muBm -= (Tprime_muBm - T0) / J;
        
        dif = Tprime(T, muBm) - T0;
    }
    return muBm;
}