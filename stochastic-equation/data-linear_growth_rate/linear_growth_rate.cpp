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

    
    // loop in parameters
    outfile << "alpha" << " " << "b0" << " " << "d" << " " <<  "sigma" << " " << "tau" << " " << "Gmean" << " " << "Gstd" << endl;


    for(alpha=alpha_min; alpha<=alpha_max; alpha+=d_alpha)
    {
        for(sigma=sigma_min; sigma<=sigma_max; sigma+=d_sigma)
        {
            double Gmean, Gstd;
            if(alpha<d_alpha)
                Gmean = calculate_meanG_nodelay(b0, d, sigma, tau, dt, tmax, samples, &Gstd);
            else
                Gmean = calculate_meanG_delay(alpha, b0, d, sigma, tau, dt, tmax, true, samples, &Gstd);

            cerr << alpha << " " << sigma << endl;
            outfile << alpha << " " << b0 << " " << d << " " <<  sigma << " " << tau << " " << Gmean << " " << Gstd << endl;
            
        }
        outfile << endl;
    }

    outfile.close();


    FREE_RANDOM;
    return 0;
}
