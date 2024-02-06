#include "function.hpp"

// extern Parameters param; // global memory access using extern
// double muBC = param.muBC;
// double rho = param.rho;
// double w = param.w;
// double angle1 = param.angle1;
// double angle12 = param.angle12;
// double angle2 = param.angle2;

//****************** RootFinding  (R,theta) to (r,h)
//************************************************************************//
double h(double R, double theta) {
    return h0 * std::pow(R, bede) * theta * (1 + aa * std::pow(theta, 2) + bb * std::pow(theta, 4));
}

double r(double R, double theta) {
    return R * (1 - std::pow(theta, 2));
}

std::tuple<double, double> FindRoot(double xr, double xh) {

  // Initial conditions for R and theta
  // Use previous solution as the initial guess.

  double theta = ((xh > 0.0) ? 1.0 : -1.0);
  double R =  200000*sqrt(xr * xr + xh * xh);

  double h_R_theta = h(R, theta);
  double r_R_theta = r(R, theta);

  // Secant method for global convergence
  double R_prev = R;
  double theta_prev = theta;

  while (sqrt((h_R_theta - xh) * (h_R_theta - xh) +
              (r_R_theta - xr) * (r_R_theta - xr)) > EPS) {

    double h_R_theta_R = (h(R + EPS, theta) - h(R_prev, theta)) / (R + EPS - R_prev);
    double h_R_theta_theta =
        (h(R, theta + EPS) - h(R, theta_prev)) / (theta + EPS - theta_prev);
    double r_R_theta_R = (r(R + EPS, theta) - r(R_prev, theta)) / (R + EPS - R_prev);
    double r_R_theta_theta =
        (r(R, theta + EPS) - r(R, theta_prev)) / (theta + EPS - theta_prev);

    double J11 = h_R_theta_R;
    double J12 = h_R_theta_theta;
    double J21 = r_R_theta_R;
    double J22 = r_R_theta_theta;

    double det_J = J11 * J22 - J12 * J21;
    double inv_J11 = J22 / det_J;
    double inv_J12 = -J12 / det_J;
    double inv_J21 = -J21 / det_J;
    double inv_J22 = J11 / det_J;

    double F1 = h_R_theta - xh;
    double F2 = r_R_theta - xr;
    R_prev = R;
    theta_prev = theta;
    R -= inv_J11 * F1 + inv_J12 * F2;
    theta -= inv_J21 * F1 + inv_J22 * F2;

    h_R_theta = h(R, theta);
    r_R_theta = r(R, theta);
  }
      
  return std::make_tuple( R, theta);
}




//******************************************************************************************************************//
