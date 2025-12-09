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

        vector<double> alpha_values = {50.};
//        for (double alpha = alpha_min; alpha <= alpha_max; alpha += d_alpha) {
//            alpha_values.push_back(alpha);
//        }

        vector<double> N_values = {10., 30., 100., 300., 1000., 3000., 10000.};

        // Print headers to outfile
        outfile << "# b0 alpha N Tmean Tstd" << endl;

    for (double alpha : alpha_values)
    {
        for (double N : N_values)
        {
           
            cerr << alpha << " " << N << endl;

            vector<double> T;

            for (int i = 0; i < samples; ++i)
                T.push_back(calculate_extinction_time(alpha, b0, d, sigma, tau, 1.0, 0.5, dt, tmax, 1./N));

            double Tmean, Tstd;
            calculate_mean_and_std(T, Tmean, Tstd);

            outfile << b0 << " " << alpha << " " << N << " " << Tmean << " " << Tstd;
            for (double t : T) // print also all the times
                outfile << " " << t;
            outfile << endl;
            outfile.flush();

            T.clear();
            T.shrink_to_fit();
        }
    }



    FREE_RANDOM;
    return 0;
}
