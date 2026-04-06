#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <iomanip> // For setting precision
#include <cmath>
#include <vector>
#include <time.h>
using namespace std;

#include "../random.h"
#include "../delayed_functions.h"

int main(int argc, char* argv[])
{
    INIT_RANDOM;

    if (argc != 2) {
        cerr << "Usage: program.x <params.dat>" << endl;
        return 1;
    }

    // Parameters
    double alpha, b0, d, sigma, tau;
    double dt, tmax;
    int samples;

    // control parameters in the loop
    double alpha_min, alpha_max, d_alpha;
    double tau_min, tau_max, d_tau;


    // read parameters from file
    if (!readParameters(argv[1], alpha, b0, d, sigma, tau, 
                        dt, tmax, samples, 
                        alpha_min, alpha_max, d_alpha, 
                        tau_min, tau_max, d_tau)) {
        return 1; // Exit if the parameters could not be read
    }



    // output file
    stringstream filename;
    filename << "output_" << argv[1];

    ofstream outfile(filename.str());
    if (!outfile) {
        cerr << "Error creating file: " << filename.str() << endl;
        return 1; // Exit if file creation fails
        }

    
    // loop in parameters

    for(alpha=alpha_min; alpha<=alpha_max; alpha+=d_alpha)
    {
        for(tau=tau_min; tau<=tau_max; tau+=d_tau)
        {

            vector<double> traj_x;
            vector<double> traj_env;

            for (int i = 0; i < samples; ++i) {
                integrate_nonlineardynamics(alpha, b0, d, sigma, tau, 1.0, 0.5, dt, tmax, 1.0, traj_x, traj_env);
            }

            double xmean, xstd;
            calculate_mean_and_std(traj_x, xmean, xstd);

            outfile << alpha << " " << " " << tau << " " << xmean << " " << xstd << endl;

            traj_x.clear();
            traj_env.clear();
            traj_x.shrink_to_fit();
            traj_env.shrink_to_fit();
        
            
        }
        outfile << endl;
    }

    outfile.close();


    FREE_RANDOM;
    return 0;
}
