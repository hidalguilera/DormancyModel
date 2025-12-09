
double integrate_lineardynamics_delay(double alpha, double b0, double d, double sigma, double tau, double x0, double dt, double tmax, bool log_integration);
double integrate_lineardynamics_nodelay(double b0, double d, double sigma, double tau, double x0, double dt, double tmax);



void integrate_nonlineardynamics(double alpha, double b0, double d, double sigma, double tau, double K, double x0, double dt, double tmax, double traj_every_dt, std::vector<double>& traj_x, std::vector<double>& traj_env);
double calculate_extinction_time(double alpha, double b0, double d, double sigma, double tau, double K, double x0, double dt, double tmax, double x_abs);

double calculate_G_delay(double alpha, double b0, double d, double sigma, double tau, double dt, double tmax, bool log_integration);
double calculate_G_nodelay(double b0, double d, double sigma, double tau, double dt, double tmax);

double calculate_meanG_delay(double alpha, double b0, double d, double sigma, double tau, double dt, double tmax, bool log_integration, int samples, double *Gstd);
double calculate_meanG_nodelay(double b0, double d, double sigma, double tau, double dt, double tmax, int samples, double *Gstd);


bool readParameters(
    const string& filename, 
    double& alpha, double& b0, double& d, double& sigma, double& tau, 
    double& dt, double& tmax, int& samples, 
    double& x_min, double& x_max, double& d_x, 
    double& y_min, double& y_max, double& d_y)

{
    std::ifstream infile(filename);
    if (!infile) {
        cerr << "Error opening file of parameters" << filename << endl;
        return false;
    }

    string line;
    while (std::getline(infile, line)) {
        string key;
        double value;
        std::stringstream ss(line);

        if (std::getline(ss, key, '=') && ss >> value) {
            key.erase(key.find_last_not_of(" \t") + 1); // Trim trailing whitespace
            if (key == "alpha") alpha = value;
            else if (key == "b0") b0 = value;
            else if (key == "d") d = value;
            else if (key == "sigma") sigma = value;
            else if (key == "tau") tau = value;
            else if (key == "dt") dt = value;
            else if (key == "tmax") tmax = value;
            else if (key == "samples") samples = static_cast<int>(value);
            else if (key == "x_min") x_min = value;
            else if (key == "x_max") x_max = value;
            else if (key == "d_x") d_x = value;
            else if (key == "y_min") y_min = value;
            else if (key == "y_max") y_max = value;
            else if (key == "d_y") d_y = value;
        }
    }


/*
    cout << "alpha = " << alpha << endl;
    cout << "b0 = " << b0 << endl;
    cout << "d = " << d << endl;
    cout << "sigma = " << sigma << endl;
    cout << "tau = " << tau << endl;
    cout << "dt = " << dt << endl;
    cout << "tmax = " << tmax << endl;
    cout << "samples = " << samples << endl;
    cout << "x_min = " << x_min << endl;
    cout << "x_max = " << x_max << endl;
    cout << "d_x = " << d_x << endl;
    cout << "y_min = " << y_min << endl;
    cout << "y_max = " << y_max << endl;
    cout << "d_y = " << d_y << endl;
    cout << endl;
*/

    infile.close();
    return true;
}

void calculate_mean_and_std(const std::vector<double>& data, double& mean, double& stddev) {
    double sum = 0.0;
    double sum2 = 0.0;
    int n = data.size();

    for (double value : data) {
        sum += value;
        sum2 += value * value;
    }

    mean = sum / n;
    stddev = sqrt((sum2 / n) - (mean * mean));
}

double integrate_lineardynamics_delay(double alpha, double b0, double d, double sigma, double tau, double x0, double dt, double tmax, bool log_integration)
{
    // check alpha is multiple of dt
    // double check_alpha = fmod(alpha, dt);

    double x = x0;
    if(log_integration)
    {
        double y = log(x);
        double env;
        if(RANDOM < 0.5)
            env = +1.;
        else env = -1.;

        // past history
        int n_past = 2*int(alpha/dt);
        double *y_past = (double *)malloc(n_past*sizeof(double));
        double *env_past = (double *)malloc(n_past*sizeof(double));

        for(int i=0; i<n_past; i++)
        {
            y_past[i] = y;
            env_past[i] = env;
        }

        int index_t_alpha = 0;
        double t = 0.;
        
        double t_next_change_env = t+RANDOM_EXP(tau);
        
        while(t<tmax)
        {

            // compute next environments
            double env_next = env;
            double env_prenext = env;
            if((t+dt)>=t_next_change_env)
            {
                env_next = -env_next;
                if((t+dt/2.)>=t_next_change_env)
                    env_prenext = -env_prenext;

                t_next_change_env += RANDOM_EXP(tau);
            }

            // RK4 - compute intermediate points
            double k1 = (b0+sigma*env_past[index_t_alpha])*exp(y_past[index_t_alpha]-y)-d;
            double k2 = (b0+sigma*env_past[(index_t_alpha+1)%n_past])*exp(y_past[(index_t_alpha+1)%n_past]-(y+dt*k1/2.))-d;
            double k3 = (b0+sigma*env_past[(index_t_alpha+1)%n_past])*exp(y_past[(index_t_alpha+1)%n_past]-(y+dt*k2/2.))-d;
            double k4 = (b0+sigma*env_past[(index_t_alpha+2)%n_past])*exp(y_past[(index_t_alpha+2)%n_past]-(y+dt*k3))-d;

            // RK4
            double y_next = y + dt/6.*(k1+2.*k2+2.*k3+k4);

            // update past, interpolating at the intermediate step x(t+dt/2)
            y_past[index_t_alpha] = y;
            y_past[index_t_alpha+1] = (y_next+y)/2.;
            
            env_past[index_t_alpha] = env;
            env_past[index_t_alpha+1] = env_prenext;
            
            index_t_alpha = (index_t_alpha+2)%n_past;


            y = y_next;
            env = env_next;

            t += dt;
        }

            free(y_past);
            free(env_past);

            return exp(y);
    }
    else
    {
        double env;
        if(RANDOM < 0.5)
            env = +1.;
        else env = -1.;
        
        // past // TODO: INCLUDE ALPHA=0 IN MY SIMULATIONS!
        int n_past = 2*int(alpha/dt);
        double *x_past = (double *)malloc(n_past*sizeof(double));
        double *env_past = (double *)malloc(n_past*sizeof(double));

        for(int i=0; i<n_past; i++)
        {
            x_past[i] = x;
            env_past[i] = env;
        }

        int index_t_alpha = 0;
        double t = 0.;
        double t_next_change_env = t+RANDOM_EXP(tau);
        
        while(t<tmax)
        {

            // compute next environments
            double env_next = env;
            double env_prenext = env;
            if((t+dt)>=t_next_change_env)
            {
                env_next = -env_next;
                if((t+dt/2.)>=t_next_change_env)
                    env_prenext = -env_prenext;

                t_next_change_env += RANDOM_EXP(tau);
            }

        


            // RK4 - compute intermediate points
            double k1 = (b0+sigma*env_past[index_t_alpha])*x_past[index_t_alpha]-d*x;
            double k2 = (b0+sigma*env_past[(index_t_alpha+1)%n_past])*x_past[(index_t_alpha+1)%n_past]-d*(x+dt*k1/2.);
            double k3 = (b0+sigma*env_past[(index_t_alpha+1)%n_past])*x_past[(index_t_alpha+1)%n_past]-d*(x+dt*k2/2.);
            double k4 = (b0+sigma*env_past[(index_t_alpha+2)%n_past])*x_past[(index_t_alpha+2)%n_past]-d*(x+dt*k3);

            // RK4
            double x_next = x + dt*(k1+2.*k2+2.*k3+k4)/6.;

            // update past, interpolating at the intermediate step x(t+dt/2)
            x_past[index_t_alpha] = x;
            x_past[index_t_alpha+1] = (x_next+x)/2.;

            env_past[index_t_alpha] = env;
            env_past[index_t_alpha+1] = env_prenext;

            index_t_alpha = (index_t_alpha+2)%n_past;

            x = x_next;
            env = env_next;

            t += dt;
 
        }

        free(x_past);
        free(env_past);
        
        return x;
    }

}

void integrate_nonlineardynamics(double alpha, double b0, double d, double sigma, double tau, double K, double x0, double dt, double tmax, double traj_every_dt, std::vector<double>& traj_x, std::vector<double>& traj_env)
{
    // check alpha is multiple of dt
    // double check_alpha = fmod(alpha, dt);


	// statistics
	double tmin = tmax/2;

	double t_next_traj = tmin;


    // initial conditions
	double x = x0;
	double env;
	if(RANDOM < 0.5)
	    env = +1.;
	else env = -1.;

    double t = 0.;
    double t_next_change_env = t+RANDOM_EXP(tau);

    // choose between non-delayed and delayed integration
	if(alpha < dt)
	{
	
	
		while(t<tmax)
		{

		    // compute next environments
		    double env_next = env;
		    double env_prenext = env;
		    if((t+dt)>=t_next_change_env)
		    {
			env_next = -env_next;
			if((t+dt/2.)>=t_next_change_env)
			    env_prenext = -env_prenext;

			t_next_change_env += RANDOM_EXP(tau);
		    }

		    // RK4 - compute intermediate points
		    double k1 = (b0+sigma*env)*x*(1.-x/K)-d*x;
		    double x1 = x+dt*k1/2.;
		    double k2 = (b0+sigma*env_prenext)*x1*(1.-x1/K)-d*x1;
		    double x2 = x+dt*k2/2.;
		    double k3 = (b0+sigma*env_prenext)*x2*(1.-x2/K)-d*x2;
		    double x3 = x+dt*k3;
		    double k4 = (b0+sigma*env_next)*x3*(1.-x3/K)-d*x3;

		    // RK4
		    double x_next = x + dt*(k1+2.*k2+2.*k3+k4)/6.;

		    x = x_next;
		    env = env_next;

		    t += dt;
		    
		    if(t>=t_next_traj)
		    {
                traj_x.push_back(x);
                traj_env.push_back(env);

	    		t_next_traj += traj_every_dt;
		    }

		}
	}
	
	else
	{
		int n_past = 2*int(alpha/dt);
		double *x_past = (double *)malloc(n_past*sizeof(double));
		double *env_past = (double *)malloc(n_past*sizeof(double));

		for(int i=0; i<n_past; i++)
		{
		    x_past[i] = x;
		    env_past[i] = env;
		}

		int index_t_alpha = 0;

		while(t<tmax)
		{

		    // compute next environments
		    double env_next = env;
		    double env_prenext = env;
		    if((t+dt)>=t_next_change_env)
		    {
			env_next = -env_next;
			if((t+dt/2.)>=t_next_change_env)
			    env_prenext = -env_prenext;

			t_next_change_env += RANDOM_EXP(tau);
		    }




		    // RK4 - compute intermediate points
		    double k1 = (b0+sigma*env_past[index_t_alpha])*x_past[index_t_alpha]*(1.-x/K)-d*x;
		    double x1 = x+dt*k1/2.;
		    double k2 = (b0+sigma*env_past[(index_t_alpha+1)%n_past])*x_past[(index_t_alpha+1)%n_past]*(1.-x1/K)-d*x1;
		    double x2 = x+dt*k2/2.;
		    double k3 = (b0+sigma*env_past[(index_t_alpha+1)%n_past])*x_past[(index_t_alpha+1)%n_past]*(1.-x2/K)-d*x2;
		    double x3 = x+dt*k3;
		    double k4 = (b0+sigma*env_past[(index_t_alpha+2)%n_past])*x_past[(index_t_alpha+2)%n_past]*(1.-x3/K)-d*x3;

		    // RK4
		    double x_next = x + dt*(k1+2.*k2+2.*k3+k4)/6.;

		    // update past, interpolating at the intermediate step x(t+dt/2)
		    x_past[index_t_alpha] = x;
		    x_past[index_t_alpha+1] = (x_next+x)/2.;

		    env_past[index_t_alpha] = env;
		    env_past[index_t_alpha+1] = env_prenext;

		    index_t_alpha = (index_t_alpha+2)%n_past;

		    x = x_next;
		    env = env_next;

		    t += dt;
		    
		    if(t>=t_next_traj)
		    {
			traj_x.push_back(x);
            traj_env.push_back(env);

			t_next_traj += traj_every_dt;
		    }



		}

		free(x_past);
		free(env_past);
	}

    return;
}



void integrate_nonlineardynamics_twopopulations(double alpha[2], double b0[2], double d[2], double sigma[2], double tau, double K, double x0[2], double xabs, double dt, double tmax, double traj_every_dt, std::vector<double> traj_x[2], std::vector<double>& traj_env)
{
    // check alpha is multiple of dt
    // double check_alpha = fmod(alpha, dt);

	double t_next_traj = 0;


	double b01 = b0[0];
	double b02 = b0[1];
	double d1 = d[0];
	double d2 = d[1];
	double sigma1 = sigma[0];
	double sigma2 = sigma[1];
	double alpha1 = alpha[0];
	if(alpha1 < dt) alpha1 = dt;
	double alpha2 = alpha[1];
	if(alpha2 < dt) alpha2 = dt;

	// choose longer posible delay to allocate past
	double alpha_max = std::max(alpha1, alpha2);
	int n_past = 2*int(alpha_max/dt);


   // initial conditions
	double x1 = x0[0];
	double x2 = x0[1];
	double xtot = x1+x2;
	
	double env;
	if(RANDOM < 0.5)
	    env = +1.;
	else env = -1.;

    double t = 0.;
    double t_next_change_env = t+RANDOM_EXP(tau);


	double *x1_past = (double *)malloc(n_past*sizeof(double));
	double *x2_past = (double *)malloc(n_past*sizeof(double));
	double *env_past = (double *)malloc(n_past*sizeof(double));
	// INTRODUCE TWO ENVIROMENTS
	for(int i=0; i<n_past; i++)
	{
		x1_past[i] = x1;
		x2_past[i] = x2;
		env_past[i] = env;
	}

	int index_t_alpha_max = 0;
	int index_t_alpha1 = 2*int((alpha_max-alpha1)/dt);
	int index_t_alpha2 = 2*int((alpha_max-alpha2)/dt);

	while(t<tmax)
	{

		// compute next environments
		double env_next = env;
		double env_prenext = env;
		if((t+dt)>=t_next_change_env)
		{
		env_next = -env_next;
		if((t+dt/2.)>=t_next_change_env)
			env_prenext = -env_prenext;

		t_next_change_env += RANDOM_EXP(tau);
		}

		// RK4 - compute intermediate points
		double k1_1 = (b01+sigma1*env_past[index_t_alpha1])*x1_past[index_t_alpha1]*(1.-xtot/K)-d1*x1;
		double x1_1 = x1+dt*k1_1/2.;
		double k2_1 = (b02+sigma2*env_past[index_t_alpha2])*x2_past[index_t_alpha2]*(1.-xtot/K)-d2*x2;
		double x2_1 = x2+dt*k2_1/2.;
		double xtot_1 = x1_1+x2_1;

		double k1_2 = (b01+sigma1*env_past[(index_t_alpha1+1)%n_past])*x1_past[(index_t_alpha1+1)%n_past]*(1.-xtot_1/K)-d1*x1_1;
		double x1_2 = x1+dt*k1_2/2.;
		double k2_2 = (b02+sigma2*env_past[(index_t_alpha2+1)%n_past])*x2_past[(index_t_alpha2+1)%n_past]*(1.-xtot_1/K)-d2*x2_1;
		double x2_2 = x2+dt*k2_2/2.;
		double xtot_2 = x1_2+x2_2;

		double k1_3 = (b01+sigma1*env_past[(index_t_alpha1+1)%n_past])*x1_past[(index_t_alpha1+1)%n_past]*(1.-xtot_2/K)-d1*x1_2;
		double x1_3 = x1+dt*k1_3;
		double k2_3 = (b02+sigma2*env_past[(index_t_alpha2+1)%n_past])*x2_past[(index_t_alpha2+1)%n_past]*(1.-xtot_2/K)-d2*x2_2;
		double x2_3 = x2+dt*k2_3;
		double xtot_3 = x1_3+x2_3;

		double k1_4 = (b01+sigma1*env_past[(index_t_alpha1+2)%n_past])*x1_past[(index_t_alpha1+2)%n_past]*(1.-xtot_3/K)-d1*x1_3;
		double k2_4 = (b02+sigma2*env_past[(index_t_alpha2+2)%n_past])*x2_past[(index_t_alpha2+2)%n_past]*(1.-xtot_3/K)-d2*x2_3;

		// RK4
		double x1_next = x1 + dt*(k1_1+2.*k1_2+2.*k1_3+k1_4)/6.;
		double x2_next = x2 + dt*(k2_1+2.*k2_2+2.*k2_3+k2_4)/6.;	

		if(x1_next<xabs) x1_next = 0;
		if(x2_next<xabs) x2_next = 0;
		

		// update past, interpolating at the intermediate step x(t+dt/2)
		x1_past[index_t_alpha_max] = x1;
		x1_past[index_t_alpha_max+1] = (x1_next+x1)/2.;

		x2_past[index_t_alpha_max] = x2;
		x2_past[index_t_alpha_max+1] = (x2_next+x2)/2.;

		env_past[index_t_alpha_max] = env;
		env_past[index_t_alpha_max+1] = env_prenext;

		index_t_alpha1 = (index_t_alpha1+2)%n_past;
		index_t_alpha2 = (index_t_alpha2+2)%n_past;
		index_t_alpha_max = (index_t_alpha_max+2)%n_past;

		x1 = x1_next;
		x2 = x2_next;

		env = env_next;

		t += dt;
		
		if(t>=t_next_traj)
		{
		traj_x[0].push_back(x1);
		traj_x[1].push_back(x2);

		traj_env.push_back(env);

		t_next_traj += traj_every_dt;
		}



	}

	free(x1_past);
	free(x2_past);
	free(env_past);

    return;
}



double calculate_extinction_time(double alpha, double b0, double d, double sigma, double tau, double K, double x0, double dt, double tmax, double x_abs)
{
    // check alpha is multiple of dt
    // double check_alpha = fmod(alpha, dt);


    // initial conditions
	double x = x0;
	double env;
	if(RANDOM < 0.5)
	    env = +1.;
	else env = -1.;

    double t = 0.;
    double t_next_change_env = t+RANDOM_EXP(tau);

    // choose between non-delayed and delayed integration
	if(alpha < dt)
	{

		while(t<tmax && x>x_abs)
		{

		    // compute next environments
		    double env_next = env;
		    double env_prenext = env;
		    if((t+dt)>=t_next_change_env)
		    {
			env_next = -env_next;
			if((t+dt/2.)>=t_next_change_env)
			    env_prenext = -env_prenext;

			t_next_change_env += RANDOM_EXP(tau);
		    }

		    // RK4 - compute intermediate points
		    double k1 = (b0+sigma*env)*x*(1.-x/K)-d*x;
		    double x1 = x+dt*k1/2.;
		    double k2 = (b0+sigma*env_prenext)*x1*(1.-x1/K)-d*x1;
		    double x2 = x+dt*k2/2.;
		    double k3 = (b0+sigma*env_prenext)*x2*(1.-x2/K)-d*x2;
		    double x3 = x+dt*k3;
		    double k4 = (b0+sigma*env_next)*x3*(1.-x3/K)-d*x3;

		    // RK4
		    double x_next = x + dt*(k1+2.*k2+2.*k3+k4)/6.;

		    x = x_next;
		    env = env_next;

		    t += dt;
		    
	    }
	}
	
	else
	{
		int n_past = 2*int(alpha/dt);
		double *x_past = (double *)malloc(n_past*sizeof(double));
		double *env_past = (double *)malloc(n_past*sizeof(double));

		for(int i=0; i<n_past; i++)
		{
		    x_past[i] = x;
		    env_past[i] = env;
		}

		int index_t_alpha = 0;

		while(t<tmax && x>x_abs)
		{

		    // compute next environments
		    double env_next = env;
		    double env_prenext = env;
		    if((t+dt)>=t_next_change_env)
		    {
			env_next = -env_next;
			if((t+dt/2.)>=t_next_change_env)
			    env_prenext = -env_prenext;

			t_next_change_env += RANDOM_EXP(tau);
		    }




		    // RK4 - compute intermediate points
		    double k1 = (b0+sigma*env_past[index_t_alpha])*x_past[index_t_alpha]*(1.-x/K)-d*x;
		    double x1 = x+dt*k1/2.;
		    double k2 = (b0+sigma*env_past[(index_t_alpha+1)%n_past])*x_past[(index_t_alpha+1)%n_past]*(1.-x1/K)-d*x1;
		    double x2 = x+dt*k2/2.;
		    double k3 = (b0+sigma*env_past[(index_t_alpha+1)%n_past])*x_past[(index_t_alpha+1)%n_past]*(1.-x2/K)-d*x2;
		    double x3 = x+dt*k3;
		    double k4 = (b0+sigma*env_past[(index_t_alpha+2)%n_past])*x_past[(index_t_alpha+2)%n_past]*(1.-x3/K)-d*x3;

		    // RK4
		    double x_next = x + dt*(k1+2.*k2+2.*k3+k4)/6.;

		    // update past, interpolating at the intermediate step x(t+dt/2)
		    x_past[index_t_alpha] = x;
		    x_past[index_t_alpha+1] = (x_next+x)/2.;

		    env_past[index_t_alpha] = env;
		    env_past[index_t_alpha+1] = env_prenext;

		    index_t_alpha = (index_t_alpha+2)%n_past;

		    x = x_next;
		    env = env_next;

		    t += dt;
		    
		}

		free(x_past);
		free(env_past);
	}

    return t;
}


double integrate_lineardynamics_nodelay(double b0, double d, double sigma, double tau, double x0, double dt, double tmax)
{
    // check alpha is multiple of dt
    // double check_alpha = fmod(alpha, dt);

    double y = log(x0);
    double env;
    if(RANDOM < 0.5)
        env = +1.;
    else env = -1.;

    double t = 0.;
    double t_next_change_env = t+RANDOM_EXP(tau);
    
    while(t<tmax)
    {

        // compute next environments
	double env_next = env;
	double env_prenext = env;
        if((t+dt)>=t_next_change_env)
        {
            env_next = -env_next;
            if((t+dt/2.)>=t_next_change_env)
            	env_prenext = -env_prenext;

	    t_next_change_env += RANDOM_EXP(tau);
        }



        double y_next = y + dt*(b0+sigma*env-d);

        y = y_next;
        env = env_next;

        t += dt;
    }
        return exp(y);
}


double calculate_G_delay(double alpha, double b0, double d, double sigma, double tau, double dt, double tmax, bool log_integration)
{
    return log(integrate_lineardynamics_delay(alpha, b0, d, sigma, tau, 1., dt, tmax, log_integration))/tmax;
}

double calculate_G_nodelay(double b0, double d, double sigma, double tau, double dt, double tmax)
{
    return log(integrate_lineardynamics_nodelay(b0, d, sigma, tau, 1., dt, tmax))/tmax;
}

double calculate_meanG_delay(double alpha, double b0, double d, double sigma, double tau, double dt, double tmax, bool log_integration, int samples, double *Gstd)
{
    double Gsum = 0.;
    double Gsum2 = 0.;
    
    double G;
    
    for(int k=0; k<samples; k++)
    {
        G = calculate_G_delay(alpha, b0, d, sigma, tau, dt, tmax, log_integration);
        Gsum += G;
        Gsum2 += G*G;
    }
    
    Gsum /= samples;
    Gsum2 /= samples;
    

    *Gstd = sqrt(Gsum2-Gsum*Gsum);
    
    return Gsum;
}


double calculate_meanG_nodelay(double b0, double d, double sigma, double tau, double dt, double tmax, int samples, double *Gstd)
{
    double Gsum = 0.;
    double Gsum2 = 0.;
    
    double G;
    
    for(int k=0; k<samples; k++)
    {
        G = calculate_G_nodelay(b0, d, sigma, tau, dt, tmax);
        Gsum += G;
        Gsum2 += G*G;
    }
    
    Gsum /= samples;
    Gsum2 /= samples;
    

    *Gstd = sqrt(Gsum2-Gsum*Gsum);
    
    return Gsum;
}
