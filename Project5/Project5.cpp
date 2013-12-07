// Project5.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "armadillo"
#include <random>
#include "SolarSystem.h"
#include "Solvers.h"
#include "CelestialBodyInitializer.h"

using namespace arma;
using namespace std;

// fit radial distribution to curve (using least squares)
// returns a vector: [ r0, n0, r0/N^(-1/3), n0/N^2 ]
vec radialDistFitLSq(mat radialDist, int N, int nNR)
{
	// these factors are educated guesses as to where we will find the best match
	const double nFact = 2.0;
	const double rFact = 1.0;

	// dependence of n0 and r0 on N taken from http://arxiv.org/abs/1011.0614
	double nScale = pow((double)N, 2.0);
	double rScale = pow((double)N, -(1.0 / 3.0));

	// range to search
	double minN0 = 0.0;
	double minR0 = 0.0;
	double maxN0 = nFact * nScale;
	double maxR0 = rFact * rScale;

	// step size for n0 and r0 trials, respectively
	double deltaN = (maxN0 - minN0) / (nNR - 1.0);
	double deltaR = (maxR0 - minR0) / (nNR - 1.0);

	// return this:
	vec ret = vec(4);

	// determine number of bins in histogram
	int boxes = radialDist.n_rows;

	double sum; // sum of deviations squared
	double minsum; // minimal sum (found so far)

	// loop over n0 values
	for (int n = 0; n < nNR; n++)
	{
		double n0 = minN0 + n * deltaN;

		// loop over r0 values
		for (int r = 0; r < nNR; r++)
		{
			double r0 = minR0 + r * deltaR;

			sum = 0.0; // reset

			// for every bin, compare value to formula
			// add square of deviation to sum
			for (int i = 0; i < boxes; i++)
			{
				double formula = n0 / (1.0 + pow((radialDist(i, 0) / r0), 4.0));
				sum += pow(formula - radialDist(i, 2), 2.0);
			}

			if ((n == 0) && (r == 0)) // first attempt: nothing to compare to
			{
				minsum = sum;
				ret(0) = r0;
				ret(1) = n0;
			}
			else
			{
				if (sum < minsum) // better match than any previous
				{
					minsum = sum;
					ret(0) = r0;
					ret(1) = n0;
				}
			}
		}
	}

	ret(2) = ret(0) / rScale; // r0/N^(-1/3)
	ret(3) = ret(1) / nScale; // n0/N^2

	return ret;
}

// (for benchmarking): convert 2D polar coordinates to cartesian
vec toCartesian2D(double r, double theta)
{
	vec p = vec(2);
	p[0] = r*cos(theta); // x
	p[1] = r*sin(theta); // y
	return p;
}

// (for benchmarking): returns an angle in radians that is orthogonal to theta
double orthogonal2D(double theta, bool clockwise = false)
{
	if (clockwise)
	{
		return (theta - 0.5 * cPI);
	}
	else // counterclockwise
	{
		return (theta + 0.5 * cPI);
	}
}

// (for benchmarking): initialize 2-body problem (in 2D)
void initial2D(CelestialBody* cb, double d, double v, double theta = -1.0)
{
	if (theta == -1.0) // no angle given
	{
		// create random angle
		long* seed = new long(-time(NULL));
		theta = (2 * cPI), GaussPDF::ran2(seed);
		delete seed; // farming in a star cluster is risky business, anyway
	}

	*(cb->position) = toCartesian2D(d, theta);
	*(cb->velocity) = toCartesian2D(v, orthogonal2D(theta)); // orthogonal (counterclockwise) vs. position
}

// (for benchmarking): returns an angle in radians that is orthogonal to theta
double orthogonal3D(double theta, bool clockwise = false)
{
	return (theta + 0.5 * cPI);
}

vec toCartesian3D(double r, double theta, double psi)
{
	vec p = vec(3);
	p[0] = r*sin(theta)*cos(psi); // x
	p[1] = r*sin(theta)*sin(psi); // y
	p[2] = r*cos(theta);

	return p;
}

void initial3D(CelestialBody* cb, double d, double v)
{
	// determine position(r0) - independent of system dimension
	int dim = 3;
	long* seedTheta = new long(-time(NULL));
	long* seedPsi = new long(-time(NULL)); // long time, no C?

	double theta;
	double psi;

	vec x = vec(dim);
	theta = (2 * cPI), GaussPDF::ran2(seedTheta);
	psi = (2 * cPI), GaussPDF::ran2(seedPsi);

	delete seedTheta;
	delete seedPsi;

	*(cb->position) = toCartesian3D(d, theta, psi); // scale to correct radius
	*(cb->velocity) = toCartesian3D(v, orthogonal3D(theta), psi); // orthogonal (counterclockwise) vs. position
}



// Entry point for console application
int _tmain(int argc, _TCHAR* argv[])
{
	#pragma region Flags and settings

	// dimensions
	const int DIM = 3;

	// initialization & time steps (common)
	const int N = 500; // number of celestial bodies
	const double EPSILON = (0.1 / sqrt(5.0)); // correction to Newton in ly to avoid infinite forces at close range
	const double R0 = 20.0; // initial radius in ly
	const double AVG_M = 2.0; // solar masses
	const double STD_M = 1.0; // solar masses
	const double CRUNCH_TIMES = 4.0; // # of crunch times to simulate for
	const int N_STEPS = 1000; // number of steps total
	const int N_PLOT = 100; // number of steps to plot (must be <= N_STEPS)
	const double STEP = CRUNCH_TIMES / ((double)N_STEPS - 1.0); // step size (in crunch times)
	const int PLOT_EVERY = N_STEPS / N_PLOT; // plot every ...th step
	const double CURVEFIT_STDDEV = 1.0; // max r limit for curve fitting in standard deviations
	const double AVG_BIN = 20.0; // avg. number of particles in each bin (curve fitting)

	// initialization & time steps (run many with different n/epsilon, same total mass)
	const int N_SIMS = 4; // number of simulations to run (set to 1 to run just once)
	const int N_END = 2000; // max N for last sim (ignored if N_SIMS == 1)
	const double EPSILON_END = 0.15; // max epsilon for last sim (ignored if N_SIMS == 1)
	const double TOTAL_M = AVG_M * N; // total mass (to be kept constant)
	const double STD_FACTOR = STD_M / AVG_M; // scale std. dev. to average

	// physical constants
	const double LY = 9.4607e15; // m
	const double MYR = cYr * 1.0e6; // s
	const double M_SUN = 1.9891e30; // kg
	const double G_YLS = cG * M_SUN * pow(cYr, 2.0) / pow(LY, 3.0); // G in years, ly, solar masses

	// flags
	const bool EPSILON_LOOP = false; // vary epsilon instead of n
	const bool USE_LEAPFROG = true; // use Leapfrog method
	const bool USE_RK4 = false; // use Runge-Kutta method
	const bool USE_EULER = false; // use Euler-Cromer method
	const bool DEBUG = false; // for debugging only
	const bool BENCHMARK = false; // To test against the project 3 code
	#pragma endregion

	if (!BENCHMARK)
	{
#pragma region Initialization

		int deltaN = 0;
		double deltaEpsilon = 0.0;

		int nSims = N_SIMS; // avoid const error

		if (!EPSILON_LOOP) // loop over n
		{
			if (nSims > 1) // avoid division by 0
			{
				deltaN = (N_END - N) / (nSims - 1);
			}
		}
		else // loop over epsilon
		{
			if (nSims > 1) // avoid division by 0
			{
				deltaEpsilon = (EPSILON_END - EPSILON) / (nSims - 1);
			}
		}

		// run one or more simulations
		for (int isim = 0; isim < N_SIMS; isim++)
		{
			ostringstream fname; // for dynamic file names

			int nParticles;
			double eps;

			// create gravity (we will update it later)
			Gravity g = Gravity(0.0, 0.0);

			if (!EPSILON_LOOP) // loop over n
			{
				// determine number of particles
				nParticles = N + isim * deltaN;

				// set epsilon (see report for explanation)
				eps = sqrt((double)N / (double)nParticles) * EPSILON;
				g.setEpsilon(eps);
			}
			else // loop over epsilon
			{
				// set number of particles
				nParticles = N;

				// set epsilon (since we use epsilon squared, we make the intervals linear for that)
				eps = EPSILON + isim * deltaEpsilon;
				g.setEpsilon(eps);
			}

			cout << endl << "--- SIMULATION " << (isim + 1) << " OF " << N_SIMS << " ---" << endl;
			cout << "N = " << nParticles << endl;
			cout << "epsilon = " << eps << endl << endl;

			// create system
			SolarSystem* system = new SolarSystem(DIM, N_STEPS, PLOT_EVERY, &g);

			// conserve the total mass & scale std.dev to average
			double avgMass = TOTAL_M / N;
			double stdMass = STD_FACTOR * avgMass;

			// add N randomly initialized celestial bodies (also sets G)
			CelestialBodyInitializer::initialize(system, nParticles, avgMass, stdMass, R0);

			// id of simulation (for file names)
			ostringstream id = ostringstream();
			id << isim;

			// call this only when initialization is 100% complete!
			Solvers solv = Solvers(system, id.str(),  USE_RK4, USE_LEAPFROG, USE_EULER);

			if (DEBUG)
			{
				for (int i = 0; i < system->n(); i++)
				{
					CelestialBody* cb = system->body(i);
					cout << "mass = " << cb->mass << endl;
					cout << "position = " << *(cb->position) << endl << endl;
				}

				cout << "T_CRUNCH = " << CelestialBodyInitializer::tCrunch(R0, system, G_YLS) << endl;
				cout << "G = " << CelestialBodyInitializer::G(R0, system) << endl;
				cout << "G_YLS = " << G_YLS << endl;
				cout << endl;
			}

#pragma endregion

#pragma region Solve and plot

			// average total energy before simulation starts
			double EtotBefore = system->EpAvg(false) + system->EkAvg(false);

			// Plot energies before simulation (all should be bound)
			cout << "E_k before: " << system->EkAvg(false) << endl;
			cout << "E_p before: " << system->EpAvg(false) << endl;
			cout << "E_tot before: " << EtotBefore << endl;

			if (DEBUG)
			{
				cout << "E_k before (bound): " << system->EkAvg(true) << endl;
				cout << "E_p before (bound): " << system->EpAvg(true) << endl;
				cout << "E_tot before (bound): " << (system->EpAvg(true) + system->EkAvg(true)) << endl;

				cout << "Center of mass (all): " << system->centerOfMass(false) << endl;
				cout << "Center of mass (bound): " << system->centerOfMass(true) << endl;

				double maxR = R0;
				mat radial = system->radialDistribution(maxR, AVG_BIN, true);
				fname = ostringstream();
				fname << "radial_before_" << isim << ".dat";
				radial.save(fname.str(), raw_ascii);

				cout << radial << endl << endl;

				cout << "n(r0) = " << system->numDens(R0) << endl; // nice to compare this to n(r) in the bins

				cout << "(bound) avg = " << system->avgDistCoM(true) << endl;
				cout << "(bound) stdDev = " << system->stdDevDistCoM(true) << endl;
			}

			cout << endl << "Running simulation";

			// this is where the magic happens :)
			vector<SolarSystem*>* systems = solv.Solve(STEP);

			// for each algorithm used
			for (int i = 0; i < systems->size(); i++)
			{
				string alg = "";
				switch (i) // name of algorithm - must match what's going on in Solvers
				{
				case 0:
					alg = "leapfrog";
					break;
				case 1:
					alg = "rk4";
					break;
				case 2:
					alg = "euler";
					break;
				}

				// get system reference
				system = systems->at(i);

				// average total/bound energy after simulation
				double Etot = system->EpAvg(false) + system->EkAvg(false);
				double EtotBound = system->EpAvg(true) + system->EkAvg(true);

				// average kinetic/potential energy for bound particles
				double EkBound = system->EkAvg(true);
				double EpBound = system->EpAvg(true);

				// distance to bound center of mass
				double avgComBound = system->avgDistCoM(true);
				double stdComBound = system->stdDevDistCoM(true);

				cout << "E_k after (bound): " << EkBound << endl;
				cout << "E_p after (bound): " << EpBound << endl;
				cout << "E_tot after: " << Etot << endl;
				cout << "E_tot after (bound): " << EtotBound << endl;

				if (DEBUG)
				{
					cout << "E_k after: " << system->EkAvg(false) << endl;
					cout << "E_p after: " << system->EpAvg(false) << endl;

					cout << "Center of mass (all): " << system->centerOfMass(false) << endl;
					cout << "Center of mass (bound): " << system->centerOfMass(true) << endl;
				}

				// get radial distribution & save to file
				// set max radius to average + x standard deviations (where we choose x)
				double maxR = system->avgDistCoM(true) + CURVEFIT_STDDEV * system->stdDevDistCoM(true);
				mat radial = system->radialDistribution(maxR, AVG_BIN, true);
				fname = ostringstream();
				fname << "radial_after_" << isim << "_" << alg << ".dat";
				radial.save(fname.str(), raw_ascii);

				if (DEBUG)
				{
					cout << radial << endl << endl;
				}

				// save the number of bound particles per time step for plotting
				fname = ostringstream();
				fname << "nbound_" << isim << "_" << alg << ".dat";
				system->nBoundPlot().save(fname.str(), raw_ascii);

				// curve fitting, save results to file
				vec testFit = radialDistFitLSq(radial, system->n(), 1000);
				fname = ostringstream();
				fname << "curveFit_" << isim << "_" << alg << ".dat";
				testFit.save(fname.str(), raw_ascii);

				if (DEBUG)
				{
					cout << "r0 = " << testFit(0) << endl;
					cout << "n0 = " << testFit(1) << endl;
					cout << "r0 / N^(-1/3) = " << testFit(2) << endl;
					cout << "n0 / N^2 = " << testFit(3) << endl << endl;

					cout << "(bound) avg = " << system->avgDistCoM(true) << endl;
					cout << "(bound) stdDev = " << system->stdDevDistCoM(true) << endl;
					cout << "(all) avg = " << system->avgDistCoM(false) << endl;
					cout << "(all) stdDev = " << system->stdDevDistCoM(false) << endl;
				}

				// calculate "classic" potential energy instead for bound particles
				// (i.e. ignoring epsilon in the potential)
				g.setEpsilon(0.0);
				system->calculate(); // re-compute potential energies
				double EpBoundClassic = system->EpAvg(true);

				// save misc. data about system to file
				vec sysdata = vec(11);
				// parameters
				sysdata(0) = nParticles;
				sysdata(1) = eps;
				// time
				sysdata(2) = solv.totalTime;
				// energy conservation
				sysdata(3) = EtotBefore;
				sysdata(4) = Etot;
				sysdata(5) = EtotBound;
				// virial theorem data
				sysdata(6) = EkBound;
				sysdata(7) = EpBound;
				sysdata(8) = EpBoundClassic;
				// distance to bound center of mass
				sysdata(9) = avgComBound;
				sysdata(10) = stdComBound;
				fname = ostringstream();
				fname << "sysdata_" << isim << "_" << alg << ".dat";
				sysdata.save(fname.str(), raw_ascii);
			}

			// free up resources
			while (systems->size() > 0)
			{
				system = systems->at(systems->size() - 1); // get last system
				systems->pop_back(); // pop it off the stack
				// NOTE: Do not delete system! That is handled when solv expires!
			}
			delete systems; // finally, delete the vector
		}
	#pragma endregion
	}
	#pragma region Benchmark

	if (BENCHMARK) // re-create Project 3 for comparison - can we reproduce the results?
	{
		//getchar(); // pause

		// Time steps
		const int STEP_BM = 24 * 60 * 60; // step length (s)
		const int N_STEPS_BM = 300 * 365; // number of steps total
		const int PLOT_EVERY_BM = 1; // plot every ...th step
		// create gravity
		Gravity g_Bench = Gravity(cG, 0);
		printf("Entering the Benchmark part \n");
		// create system
		SolarSystem* system_BM = new SolarSystem(3, N_STEPS_BM, PLOT_EVERY_BM, &g_Bench);
		const double M_SUN_2 = 2e30;
		const double M_EARTH = 6e24;
		const double D_SUN = 0.0 * cAU;
		const double D_EARTH = 1.0 * cAU;
		const double V_SUN = 0.0;
		const double V_EARTH = 29.8e3;
		//system->grav()->setG(cG);
		CelestialBody* sun = new CelestialBody("Sun", M_SUN_2, system_BM, true);
		CelestialBody* earth = new CelestialBody("Earth", M_EARTH, system_BM, false);
		initial3D(earth, D_EARTH, V_EARTH);

		double E_k_init = system_BM->EkAvg(false);
		double E_p_init = system_BM->EpAvg(false);
		double E_tot_init = system_BM->EpAvg(false) + system_BM->EkAvg(false);
		

		// call this only when initialization is 100% complete!
		Solvers solv = Solvers(system_BM, "BM", USE_RK4, USE_LEAPFROG, USE_EULER);

		// this is where the magic happens :)
		vector<SolarSystem*>* systemsBM = solv.Solve(STEP_BM);
		printf("Aaand ... We're done \n");

		cout << "E_k before : " << E_k_init << endl;
		cout << "E_p before : " << E_p_init << endl;
		cout << "E_tot before : " << E_tot_init << endl;

		cout << "E_tot after: " << (system_BM->EpAvg(false) + system_BM->EkAvg(false)) << endl;		
		cout << "E_k after: " << system_BM->EkAvg(false) << endl;
		cout << "E_p after: " << system_BM->EpAvg(false) << endl;
	}

	cout << endl << "Press ENTER to exit...";
	getchar(); // pause

	#pragma endregion 
	return 0;
}