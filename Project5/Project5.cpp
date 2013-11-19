// Project5.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "armadillo"

#include "SolarSystem.h"
#include "Solvers.h"
#include "CelestialBodyInitializer.h"

using namespace arma;
using namespace std;

// fit radial distribution to curve (using least squares)
vec radialDistFitLSq(mat radialDist, int nNR, double minR0, double maxR0, double minN0, double maxN0)
{
	vec ret = vec(2); // vec(0) = r0, vec(1) = n0

	int boxes = radialDist.n_rows;
	double sum;
	double minsum;

	double deltaN = (maxN0 - minN0) / (nNR - 1.0);
	double deltaR = (maxR0 - minR0) / (nNR - 1.0);

	for (int n = 0; n < nNR; n++)
	{
		double n0 = minN0 + n * deltaN;

		for (int r = 0; r < nNR; r++)
		{
			double r0 = minR0 + r * deltaR;

			sum = 0.0;
			for (int i = 0; i < boxes; i++)
			{
				double formula = n0 / (1.0 + pow(radialDist(i, 0) / r0, 4.0));
				sum += pow(formula - radialDist(i, 1), 2.0);
			}

			if ((n == 0) && (r == 0))
			{
				minsum = sum;
				ret(0) = r0;
				ret(1) = n0;
			}
			else
			{
				if (sum < minsum)
				{
					minsum = sum;
					ret(0) = r0;
					ret(1) = n0;
				}
			}
		}
	}

	return ret;
}

// Entry point for console application
int _tmain(int argc, _TCHAR* argv[])
{
	#pragma region Flags and settings

	// dimensions
	const int DIM = 3;

	// randomization seeds (must be negative)
	long IDUM = -1337;
	long IDUM2 = -1;

	// initialization
	const int N = 100; // number of celestial bodies
	const double R0 = 20.0; // initial radius in ly
	const double AVG_M = 10.0; // solar masses
	const double STD_M = 1.0; // solar masses

	// physical constants
	const double LY = 9.4607e15; // m
	const double MYR = cYr * 1.0e6; // s
	const double M_SUN = 1.9891e30; // kg
	const double G_YLS = cG * M_SUN * pow(cYr, 2.0) / pow(LY, 3.0); // G in years, ly, solar masses
	const double EPSILON = 0.0; // correction to Newton in ly to avoid infinite forces at close range

	// time steps
	const int N_STEPS = 1000; // number of steps total
	const int N_PLOT = 1000; // number of steps to plot (must be <= N_STEPS)
	const double CRUNCH_TIMES = 10.0; // # of crunch times to simulate for
	const double STEP = CRUNCH_TIMES / ((double)N_STEPS - 1.0); // step size (in crunch times)
	const int PLOT_EVERY = N_STEPS / N_PLOT; // plot every ...th step

	// flags
	const bool USE_LEAPFROG = true; // use Leapfrog method
	const bool USE_RK4 = false; // use Runge-Kutta method
	const bool USE_EULER = false; // use Euler-Cromer method
	const bool DEBUG = true; // for debugging only

	#pragma endregion

	#pragma region Initialization

	if (DEBUG)
	{
		cout << "T_CRUNCH = " << CelestialBodyInitializer::tCrunch(R0, AVG_M, N, G_YLS) << endl;
		cout << "G = " << CelestialBodyInitializer::G(R0, AVG_M, N) << endl;
		cout << "G_YLS = " << G_YLS << endl;
		cout << endl;
	}

	// create gravity (we will update G later)
	Gravity g = Gravity(0.0, EPSILON);

	// create system
	SolarSystem* system = new SolarSystem(DIM, N_STEPS, PLOT_EVERY, &g);

	// add N randomly initialized celestial bodies
	CelestialBodyInitializer::initialize(system, N, AVG_M, STD_M, R0, &IDUM, &IDUM2);

	// call this only when initialization is 100% complete!
	Solvers solv = Solvers(system, USE_RK4, USE_LEAPFROG, USE_EULER);

	if (DEBUG)
	{
		for (int i = 0; i < system->n(); i++)
		{
			CelestialBody* cb = system->body(i);
			cout << "mass = " << cb->mass << endl;
			cout << "position = " << *(cb->position) << endl << endl;
		}
	}

	#pragma endregion

	#pragma region Solve and plot

	cout << "E_k before: " << system->EkAvg(false) << endl;
	cout << "E_k before (bound): " << system->EkAvg(true) << endl;
	cout << "E_p before: " << system->EpAvg(false) << endl;
	cout << "E_p before (bound): " << system->EpAvg(true) << endl;
	cout << "E_tot before: " << (system->EpAvg(false) + system->EkAvg(false)) << endl;
	cout << "E_tot before (bound): " << (system->EpAvg(true) + system->EkAvg(true)) << endl;

	if (DEBUG)
	{
		cout << "Center of mass: " << system->centerOfMass() << endl;

		mat radial = system->radialDistribution(R0, 100, true);
		radial.save("test.dat", raw_ascii);

		vec testFit = radialDistFitLSq(radial, 100, 0.0, 100.0, 0.0, 100.0);
		cout << "r0 = " << testFit(0) << endl;
		cout << "n0 = " << testFit(1) << endl;

		cout << "avg = " << system->avgDistCoM(true) << endl;
		cout << "stdDev = " << system->stdDevDistCoM(true) << endl;
	}

	cout << endl << "Running simulation";

	// this is where the magic happens :)
	system = solv.Solve(STEP, PLOT_EVERY);

	if (system != nullptr)
	{
		cout << "E_k after: " << system->EkAvg(false) << endl;
		cout << "E_k after (bound): " << system->EkAvg(true) << endl;
		cout << "E_p after: " << system->EpAvg(false) << endl;
		cout << "E_p after (bound): " << system->EpAvg(true) << endl;
		cout << "E_tot after: " << (system->EpAvg(false) + system->EkAvg(false)) << endl;
		cout << "E_tot after (bound): " << (system->EpAvg(true) + system->EkAvg(true)) << endl;

		if (DEBUG)
		{
			cout << "Center of mass: " << system->centerOfMass() << endl;

			mat radial = system->radialDistribution(R0, 100, true);
			radial.save("test.dat", raw_ascii);

			system->nBoundPlot().save("test2.dat", raw_ascii);

			vec testFit = radialDistFitLSq(radial, 100, 0.0, 100.0, 0.0, 100.0);
			cout << "r0 = " << testFit(0) << endl;
			cout << "n0 = " << testFit(1) << endl;

			cout << "avg = " << system->avgDistCoM(true) << endl;
			cout << "stdDev = " << system->stdDevDistCoM(true) << endl;
		}
	}

	cout << endl << "Press ENTER to exit...";
	getchar(); // pause

	#pragma endregion

	return 0;
}