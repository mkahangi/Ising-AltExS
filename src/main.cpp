/**
 * @file main.cpp
 * @brief This code uses Ising model Universality class to introduce a critical point and an alternative expansion scheme to Extend the QCD based Eqaution of state to finite densities and match lattice QCD at low densities
 *
 *
 * @author Micheal Kahangirwe <mkahangi@central.uh.edu>
 */

/**
 * @brief main function
 *
 * Takes inputs from w,  rho  and muBC from the USER
 * The main procedure is:
 *  - Maps Ising co-ordinates  to Alternative expansion scheme then to QCD co-ordinates
 *  - Merges Lattice to Ising
 *  Computes all the Thermodynamic Observables
 *  - Baryon Density, Pressure, Entropy, Baryon number susceptability, Energy, Speed of Sound
 *
 *
 * @return Out puts thermodynamic observables
 */
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




#include "function.hpp"

namespace fs = std::filesystem;

// declare global variables
extern Parameters param; 

// double inter(double T, double muB) {
//     return boost::math::quadrature::trapezoidal(
//         [&](double x) { return pow(T,4)*chi2(T) + nBfull(T, x); }, 0.0, muB);
// }

int main(int argc, char *argv[])
{
    char buff[FILENAME_MAX];

	get_options(argc,argv);
	
	extern int optind;


 	/* Assign the name of the main folder where the program lives and 
	 * the files we wish to import are located. */
	getcwd(buff,FILENAME_MAX);	/* store current working dir in buf */

	read_param(argv[optind]);		/* read parameter file */

   

   
clock_t tStart = clock(); // Time

    

    // double muBC = param.muBC;
    // double rho = param.rho;
    // double w = param.w;
    // double angle1 = param.angle1;
    // double angle12 = param.angle12;
    // double angle2 = param.angle2;

    // Define matrix dimensions
    int muB_min = param.lowMU_out;
    int muB_max = param.highMU_out + param.T_step;
    // double muB_step = 1;
    int NmuB = round((muB_max - muB_min) / param.muB_step );

    int  T_min = param.lowMU_out-param.T_step;
    int  T_max = param.highT_out + param.T_step;
    // double T_step = 1;
    int NdT = round((T_max - T_min) / param.T_step );

    int NdT_out = round((param.highT_out - param.lowT_out) / param.T_step );
    int NmuB_out = round((param.highMU_out - param.lowMU_out) / param.muB_step);

  


   

    // Define the Parameter choice
   
    cout << "Grid Range:\n";
    cout << left << setw(25) << "    T_range " << "=  [" << param.lowT_out << " MeV, " <<param.highT_out <<  " MeV] \n";
    cout << left << setw(25) << "    T_step" << "=   " << param.T_step << " MeV\n";
    cout << left << setw(25) << "    muB_range " << "=  [" << param.lowMU_out << " MeV, " <<param.highMU_out <<  " MeV] \n";
    cout << left << setw(25) << "    muB_step" << "=   " << param.muB_step << " MeV\n";

    cout << "Parameters Chosen:\n";
    cout << left << setw(25) << "    muBC" << "= " << param.muBC << " MeV\n";
    cout << left << setw(25) << "    TC" << "= " << TC(param.muBC) << " MeV\n";
    cout << left << setw(25) << "    alpha1" << "= " << radiansToDegrees(param.angle1) << " degrees\n";
    cout << left << setw(25) << "    alpha2" << "= " << radiansToDegrees(param.angle2) << " degrees\n";
    cout << left << setw(25) << "    alpha12" << "= " << radiansToDegrees(param.angle12) << " degrees\n";
    cout << left << setw(25) << "    w" << "= " << param.w << '\n';
    cout << left << setw(25) << "    rho" << "= " << param.rho << '\n';

    




    // Create an empty matrix
    vector<double> Press0_vec(NdT, 0.0);
    vector<double> dPress0dT_vec(NdT, 0.0);
    vector<double> chi20_vec(NdT, 0.0);
    vector<vector<double>> nB_mat(NdT, vector<double>(NmuB, 0.0));
    vector<vector<double>> chi2_mat(NdT, vector<double>(NmuB, 0.0));
    vector<vector<double>> press_mat(NdT, vector<double>(NmuB, 0.0));

    //Thermodynamic Observables
    vector<vector<double>> Pressure(NdT, vector<double>(NmuB, 0.0));
    vector<vector<double>> dPressdT(NdT, vector<double>(NmuB, 0.0));
    vector<vector<double>> d2PressdT2(NdT, vector<double>(NmuB, 0.0));
    vector<vector<double>> dnBdT(NdT, vector<double>(NmuB, 0.0));
    vector<vector<double>> dnBdT_mat(NdT, vector<double>(NmuB, 0.0));
    vector<vector<double>> Entro_mat(NdT, vector<double>(NmuB, 0.0));
    vector<vector<double>> Entropy_mat(NdT, vector<double>(NmuB, 0.0));
    vector<vector<double>> d2nBdTdmuB(NdT, vector<double>(NmuB, 0.0));
    vector<vector<double>> dnBdmuB(NdT, vector<double>(NmuB, 0.0));
    vector<vector<double>> Energy(NdT, vector<double>(NmuB, 0.0));
    vector<vector<double>> Cv(NdT, vector<double>(NmuB, 0.0));
    vector<vector<double>> Cs2(NdT, vector<double>(NmuB, 0.0));
    


    //Creating floders to store the output
    std::cout << std::setprecision(2) << std::fixed;

    const string prefix = "FilePar_" + to_string(param.muBC) + "_" + to_string(TC(param.muBC)) + "_" + to_string(param.w) + "_" + to_string(param.rho);
    fs::path folderPath = "../output/" + prefix;

    if (fs::exists(folderPath)) {
        cout << "Folder '" << folderPath << "' already exists in the output folder.\n\n";
        cout<< "Please delete the file in the output folder and try again.\n\n";
        return 1;
    }

        fs::create_directories(folderPath);
        cout << "Folder " << folderPath << " created successfully.\n\n";

    // Store a sample function in the matrix
    cout << "***************** /Creating Baryon density matrix and storing in memory/ *****************\n\n";

    
    storevector(Press0_vec, NdT, param.T_step,chi0);
    storevector(dPress0dT_vec, NdT, param.T_step,dchi0dT);
    storevector(chi20_vec, NdT, param.T_step, chi2);
    storeMatrix(nB_mat, NdT, param.T_step, NmuB, param.muB_step, nBfull);
    // storeMatrix(dnBdT_mat, NdT, param.T_step, NmuB, param.muB_step, dnBfulldT);
    cout << "***************** /Creating chi2 density Matrix and storing in memory/ *****************\n\n";
    storeMatrix(chi2_mat, NdT, param.T_step, NmuB, param.muB_step, chi2full);

    cout << "***************** /Calculating  Pressure density by integrating BaryonDensity Matrix / *****************\n\n";
    press_mat = integrateMatrix(nB_mat, NdT, param.T_step, NmuB, param.muB_step,true);
    

    for (int j = 0; j < NmuB; j+=1) {
        for (int i = 0; i < NdT; i+=1) {
            double Tval = static_cast<double>(i) * param.T_step;
            Pressure[i][j] = press_mat[i][j] +  Press0_vec[i]*pow(Tval,4);
            

        }
    }
    //Calculate the derivatives of both Temperature and Baryon density
    dPressdT = deriv_matrix(Pressure, NdT, param.T_step, NmuB, param.muB_step, 0);
    d2PressdT2 = deriv_matrix(dPressdT, NdT, param.T_step, NmuB, param.muB_step,0);
    dnBdT = deriv_matrix(nB_mat, NdT, param.T_step, NmuB, param.muB_step,0);
    Entro_mat = integrateMatrix(dnBdT, NdT, param.T_step, NmuB, param.muB_step,true);
    dnBdmuB = deriv_matrix(nB_mat, NdT, param.T_step, NmuB, param.muB_step, 1);
    d2nBdTdmuB = deriv_matrix(dnBdT,  NdT, param.T_step, NmuB, param.muB_step, 1);

     for (int j = 0; j < NmuB; j+=1) {
        for (int i = 0; i < NdT; i+=1) {
            double Tval = static_cast<double>(i) * param.T_step;
            Entropy_mat[i][j] = dPress0dT_vec[i]*pow(Tval,4) + 4*pow(Tval,3)*Press0_vec[i]+  Entro_mat[i][j];

        }
    }
    
    // Print the matrix
    
    // Files created for thermodynamic observables
    //File name correspond to parameters chosen
    const string FileTest =  "Test_" + prefix + ".dat";
    const string FileBarDen =  "BarDensity_" + prefix + ".dat";
    const string FileChi2Den =  "chi2density_" + prefix + ".dat";
    const string FilePressBarDen =  "Pressdensity_" + prefix + ".dat";
    const string FileEntropyBarDen =  "Entropysdensity_" + prefix + ".dat";
    const string FileEnergyBarDen =  "Energysdensity_" + prefix + ".dat";
    const string FileCvBarDen =  "CvDensity_" + prefix + ".dat";
    const string FileSpeedofSoundBarDen =  "SpeedofSound_" + prefix + ".dat";



        cout << "***************** /Outputting all thermodynamic Variable/ *****************\n\n";
        ofstream Test(folderPath / FileTest);
        
        ofstream Bardensity(folderPath / FileBarDen);     // Relative path
        ofstream chi2density(folderPath / FileChi2Den);     // Relative path
        ofstream Pressdensity(folderPath / FilePressBarDen);     // Relative path
        ofstream Entropysdensity(folderPath / FileEntropyBarDen);     // Relative path
        ofstream Energysdensity(folderPath / FileEnergyBarDen);     // Relative path
        ofstream Cvdensity(folderPath / FileCvBarDen);     // Relative path
        ofstream SpeedofSound(folderPath / FileSpeedofSoundBarDen);     // Relative path

        for (double muB = 0; muB < param.highMU_out; muB += param.muB_step) {
            for (double T = 0; T < param.highT_out; T += param.T_step) {
            if (fabs(T - TC(param.muBC)) < 3) {
                param.T_step = 0.05;
            } else {
                param.T_step = 0.50;
            }
            
            Test << muB << " " << T << " " << BardIsingA(T,muB) << endl;
        }
    }

    // for (double T = 0; T < param.highT_out; T += param.T_step) {
    //      Test  << T <<  " " << dBardIsingdmuBA(T,500) << " " << d3BardIsingdmuB3N(T,500) <<" " << TyTcrit(T,500) << endl;
    // }
    

    for (int j = 0; j < NmuB_out; j+=1) {
        for (int i = 0; i < NdT_out; i+=1) {
            double Tval = static_cast<double>(i) * param.T_step;
            double muBval = static_cast<double>(j) * param.muB_step;

            Energy[i][j] = - Pressure[i][j] + Tval*dPressdT[i][j] + muBval*nB_mat[i][j];  

            Cv[i][j] = Tval*(d2PressdT2[i][j]*chi2_mat[i][j]- pow(dnBdT[i][j],2))/chi2_mat[i][j]; 
            
            Cs2[i][j] = (pow(nB_mat[i][j],2)*d2PressdT2[i][j] - 2*dPressdT[i][j]*nB_mat[i][j]*dnBdT[i][j] + pow(dPressdT[i][j],2) *chi2_mat[i][j])/((Energy[i][j] + Pressure[i][j])*(d2PressdT2[i][j]*chi2_mat[i][j] - pow(dnBdT[i][j],2)));


            Bardensity << muBval << " " << Tval << " " << nB_mat[i][j]/pow(Tval,3)  << endl;
            chi2density << muBval << " " << Tval << " " <<chi2_mat[i][j]/pow(Tval,2) << endl;
            Pressdensity << muBval << " " << Tval <<" " << Pressure[i][j]/pow(Tval,4)  << endl;
            Entropysdensity << muBval << " " << Tval << " " << dPressdT[i][j]/pow(Tval,3) <<" "<<Entropy_mat[i][j]/pow(Tval,3)<< endl;
            Energysdensity << muBval << " " << Tval << " " << Energy[i][j]/pow(Tval,4)  << endl;
            Cvdensity << muBval << " " << Tval << " " << Cv[i][j]/pow(Tval,3)  << endl;
            SpeedofSound << muBval << " " << Tval << " " << Cs2[i][j]  << endl;


        }
    }

    Test.close();
    Bardensity.close();
    Pressdensity.close();
    Entropysdensity.close();
    Energysdensity.close();
    Cvdensity.close();
    SpeedofSound.close();
    cout << "***************** //// Successfully Completed ////  *****************\n\n";

   printf("Time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC); 
    return 0;
}









