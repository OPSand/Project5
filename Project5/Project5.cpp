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
vec radialDistFitLSq(mat radialDist, double minR0, double maxR0, double minN0, double maxN0, int nNR)
{
	vec ret = vec(2); // vec(0) = r0, vec(1) = n0

	int boxes = radialDist.n_rows;
	double sum;
	double minsum;

	double deltaN = (maxN0 - minN0) / (nNR - 1.0);
	double deltaR = (maxR0 - minR0) / (nNR - 1.0);

	// debug
	cout << radialDist << endl << endl;

	// radialDist gives N(r)
	// we want n(r) = N(r)/V(r)
	double rPrev = 0.0;
	for (int i = 0; i < boxes; i++)
	{
		double r = radialDist(i, 0);
		double delta_r = (r - rPrev);
		double Vmin = CelestialBodyInitializer::volume(r - delta_r, 3);
		double Vmax = CelestialBodyInitializer::volume(r + delta_r, 3);
		double V = (Vmax - Vmin);

		assert(V > 0.0);
		
		radialDist(i, 1) = radialDist(i, 1) / V;

		rPrev = r + delta_r;
	}

	// debug
	cout << radialDist << endl << endl;

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


vec toCartesian2D(double r, double theta)
{
	vec p = vec(2);
	p[0] = r*cos(theta); // x
	p[1] = r*sin(theta); // y
	return p;
}
// returns an angle in radians that is orthogonal to theta
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
// return a random floating point number in the interval [minIncl, maxExcl).
double randDbl(double minIncl, double maxExcl, minstd_rand* eng)
{
	uniform_real<double> rd(0.0, (2 * cPI));
	return rd(*eng);
}
void initial2D(CelestialBody* cb, double d, double v, minstd_rand* eng, double theta = -1.0)
{
	if (theta == -1.0) // no angle given
	{
		// create random angle
		theta = randDbl(0.0, (2 * cPI), eng);
	}
	*(cb->position) = toCartesian2D(d, theta);
	*(cb->velocity) = toCartesian2D(v, orthogonal2D(theta)); // counterclockwise
}


// Entry point for console application
int _tmain(int argc, _TCHAR* argv[])
{
	#pragma region Flags and settings

	// dimensions
	const int DIM = 2;

	// initialization
	const int N = 2; // number of celestial bodies
	const double R0 = 20.0; // initial radius in ly
	const double AVG_M = 10.0; // solar masses
	const double STD_M = 1.0; // solar masses

	// physical constants
	const double LY = 9.4607e15; // m
	const double MYR = cYr * 1.0e6; // s
	const double M_SUN = 1.9891e30; // kg
	const double G_YLS = cG * M_SUN * pow(cYr, 2.0) / pow(LY, 3.0); // G in years, ly, solar masses
	const double EPSILON = 0;//sqrt(0.0225); // correction to Newton in ly to avoid infinite forces at close range

	// time steps
	const int N_STEPS = 1000; // number of steps total
	const int N_PLOT = 1000; // number of steps to plot (must be <= N_STEPS)
	const double CRUNCH_TIMES = 5; // # of crunch times to simulate for
	const double STEP = CRUNCH_TIMES / ((double)N_STEPS - 1.0); // step size (in crunch times)
	const int PLOT_EVERY = N_STEPS / N_PLOT; // plot every ...th step
	
	// flags
	const bool USE_LEAPFROG = false; // use Leapfrog method
	const bool USE_RK4 = true; // use Runge-Kutta method
	const bool USE_EULER = false; // use Euler-Cromer method
	const bool DEBUG = true; // for debugging only
	const bool BENCHMARK = true; // To test against the project 3 code
	#pragma endregion

	#pragma region Initialization

	// create gravity (we will update G later)
	Gravity g = Gravity(0.0, EPSILON);

	// create system
	SolarSystem* system = new SolarSystem(DIM, N_STEPS, PLOT_EVERY, &g);

	// add N randomly initialized celestial bodies
	CelestialBodyInitializer::initialize(system, N, AVG_M, STD_M, R0);

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

		cout << "T_CRUNCH = " << CelestialBodyInitializer::tCrunch(R0, system->avgMass(), system->n(), G_YLS, DIM) << endl;
		cout << "G = " << CelestialBodyInitializer::G(R0, system->avgMass(), system->n(), DIM) << endl;
		cout << "G_YLS = " << G_YLS << endl;
		cout << endl;
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

		mat radial = system->radialDistribution(R0, 10, true);
		radial.save("radial_before.dat", raw_ascii);

		vec testFit = radialDistFitLSq(radial, 0.0, 100.0, 0.0, 1.0, 100);
		cout << "r0 = " << testFit(0) << endl;
		cout << "n0 = " << testFit(1) << endl;

		cout << "(bound) avg = " << system->avgDistCoM(true) << endl;
		cout << "(bound) stdDev = " << system->stdDevDistCoM(true) << endl;
	}

	cout << endl << "Running simulation";

	// this is where the magic happens :)
	system = solv.Solve(STEP);

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

			mat radial = system->radialDistribution(R0, 10, true);
			radial.save("radial_after.dat", raw_ascii);

			system->nBoundPlot().save("nbound.dat", raw_ascii);

			vec testFit = radialDistFitLSq(radial, 0.0, 100.0, 0.0, 1.0, 100);
			cout << "r0 = " << testFit(0) << endl;
			cout << "n0 = " << testFit(1) << endl;

			cout << "(bound) avg = " << system->avgDistCoM(true) << endl;
			cout << "(bound) stdDev = " << system->stdDevDistCoM(true) << endl;
			cout << "(all) avg = " << system->avgDistCoM(false) << endl;
			cout << "(all) stdDev = " << system->stdDevDistCoM(false) << endl;
		}
	}

	cout << endl << "Press ENTER to exit...";
	getchar(); // pause

	#pragma endregion
	
	#pragma region Benchmark

	if (BENCHMARK)
	{
		// Time steps
		const int N_STEPS = 300 * 365; // number of steps total
		const int N_PLOT = 300 * 365; // number of steps to plot (must be <= N_STEPS)
		const double STEP = 24 * 60 * 60;
		const int PLOT_EVERY = N_STEPS / N_PLOT; // plot every ...th step

		// create gravity (we will update G later)
		Gravity g_Bench = Gravity(cG, EPSILON);
		printf("Entering the Benchmark part \n");
		// create system
		SolarSystem* system_BM = new SolarSystem(DIM, N_STEPS, PLOT_EVERY, &g_Bench);
		const double M_SUN_2 = 2e30;
		const double M_EARTH = 6e24;
		const double D_SUN = 0.0 * cAU;
		const double D_EARTH = 1.0 * cAU;
		const double V_SUN = 0.0;
		const double V_EARTH = 29.8e3;
		//system->grav()->setG(cG);
		CelestialBody* sun = new CelestialBody("Sun", M_SUN_2, system_BM, true);
		CelestialBody* earth = new CelestialBody("Earth", M_EARTH, system_BM, false);
		// intitialze random number engine
		minstd_rand eng;
		int randomSeed = (int)clock();
		eng.seed(randomSeed);
		initial2D(earth, D_EARTH, V_EARTH, &eng);


		// call this only when initialization is 100% complete!
		Solvers solv = Solvers(system_BM, USE_RK4, USE_LEAPFROG, USE_EULER);
		printf("Aaand ... We're done \n");

		// this is where the magic happens :)
		system_BM = solv.Solve(STEP);
		getchar();
	}

	#pragma endregion 
	return 0;
}