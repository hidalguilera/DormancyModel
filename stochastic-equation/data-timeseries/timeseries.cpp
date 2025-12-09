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
    double sigma_min, sigma_max, d_sigma;


    // read parameters from file
    if (!readParameters(argv[1], alpha, b0, d, sigma, tau, 
                        dt, tmax, samples, 
                        alpha_min, alpha_max, d_alpha, 
                        sigma_min, sigma_max, d_sigma)) {
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


        // Print headers to outfile
        outfile << "# alpha sigma xmean xstd" << endl;
        outfile << "# x(t) for histogram" << endl;


    vector<double> traj_x;
    vector<double> traj_env;
    double dt_measure = dt;

    for (int s = 0; s<samples; s++)
    {
        integrate_nonlineardynamics(alpha, b0, d, sigma, tau, 1.0, 0.5, dt, tmax, dt_measure, traj_x, traj_env);

        // print sample
        outfile << "# sample " << s << endl;

        double t = 0.;
        for(int k = 0; k < traj_x.size(); ++k)
        {
            outfile << t << " " << traj_x[k] << endl;
            t += dt_measure;   
        }

        outfile << endl << endl;
    }


    traj_x.clear();
    traj_env.clear();
    traj_x.shrink_to_fit();
    traj_env.shrink_to_fit();



    FREE_RANDOM;
    return 0;
}
