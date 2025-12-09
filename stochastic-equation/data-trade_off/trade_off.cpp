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

#ifdef __linux__ // Linux system
    #include "/home/jhidalgo/Dropbox/random.h"
#elif __APPLE__ // macOS system
    #include "/Users/jhidalgo/Dropbox/random.h"
#elif _WIN32 // Windows system
    #include "C:\\path\\to\\random.h"
#else
    #error "Unsupported operating system. Please define the correct path for random.h."
#endif

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

        vector<double> alpha_values;
        for (double alpha = alpha_min; alpha <= alpha_max; alpha += d_alpha) {
            alpha_values.push_back(alpha);
        }

        vector<double> sigma_values;
        for (double sigma = sigma_min; sigma <= sigma_max; sigma += d_sigma) {
            sigma_values.push_back(sigma);
        }

        // Print headers to outfile
        outfile << "alpha b0 d sigma tau Gmean Gstd Tmean Tstd" << endl;
    double N = 1000.;

    for (double alpha : alpha_values)
    {
        for (double sigma : sigma_values)
        {
//            cerr << alpha << " " << sigma << endl;

// calculate Gmean
	    double Gmean, Gstd;

            if(alpha<d_alpha)
                Gmean = calculate_meanG_nodelay(b0, d, sigma, tau, dt, tmax, samples, &Gstd);
            else
                Gmean = calculate_meanG_delay(alpha, b0, d, sigma, tau, dt, tmax, true, samples, &Gstd);     



// calculate Tmean            
            vector<double> T;

            for (int i = 0; i < samples; ++i)
                T.push_back(calculate_extinction_time(alpha, b0, d, sigma, tau, 1.0, 0.5, dt, tmax, 1./N));

                double Tmean, Tstd;
                calculate_mean_and_std(T, Tmean, Tstd);



                T.clear();
                T.shrink_to_fit();


	outfile << alpha << " " << b0 << " " << d << " " << sigma << " " << tau << " " << Gmean << " "  << Gstd << " " << Tmean  << " " << Tstd << endl;


        }
    }



    FREE_RANDOM;
    return 0;
}
