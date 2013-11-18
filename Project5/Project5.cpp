// Project5.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "armadillo"

#include "SolarSystem.h"
#include "Solvers.h"
#include "CelestialBodyInitializer.h"

using namespace arma;
using namespace std;

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
	const int N = 10; // number of celestial bodies
	const double R0 = 20.0; // initial radius in ly
	const double AVG_M = 10.0; // solar masses
	const double STD_M = 1.0; // solar masses

	// physical constants
	const double LY = 9.4607e15; // m
	const double MYR = cYr * 1.0e6; // s
	const double M_SUN = 1.9891e30; // kg
	const double G_YLS = cG * M_SUN * pow(cYr, 2.0) / pow(LY, 3.0); // G in years, ly, solar masses
	const double EPSILON = 0.0; // correction to Newton in ly to avoid infinite forces at close range

	// calculate time scale & G in correct units
	const double V0 = (4.0 / 3.0)* cPI *pow(R0, 3.0); // initial volume [ly^3]
	const double RHO0 = N * AVG_M / V0; // initial mass density [solar masses / ly^3]
	const double T_CRUNCH = sqrt(3.0 * cPI / (32.0 * G_YLS * RHO0)); // crunch time [years]
	const double G = G_YLS * pow(T_CRUNCH, 2.0); // G in t_crunch, ly, solar masses

	// time steps
	const int N_STEPS = 1000; // number of steps total
	const int N_PLOT = 1000; // number of steps to plot (must be <= N_STEPS)
	const double CRUNCH_TIMES = 5.0; // # of crunch times to simulate for
	const double STEP = CRUNCH_TIMES / ((double)N_STEPS - 1.0); // step size (in crunch times)
	const int PLOT_EVERY = N_STEPS / N_PLOT; // plot every ...th step

	// flags
	const bool USE_LEAPFROG = true; // use Leapfrog method
	const bool USE_RK4 = false; // use Runge-Kutta method
	const bool USE_EULER = false; // use Euler-Cromer method
	const bool DEBUG = false; // for debugging only

	#pragma endregion

	#pragma region Initialization

	if (DEBUG)
	{
		cout << "T_CRUNCH = " << T_CRUNCH << endl;
		cout << "G = " << G << endl;
		cout << "G_YLS = " << G_YLS << endl;
	}

	// create gravity
	Gravity g = Gravity(G, EPSILON);

	// create system
	SolarSystem system = SolarSystem(DIM, N_STEPS, PLOT_EVERY, &g);

	// add N randomly initialized celestial bodies
	CelestialBodyInitializer::initialize(&system, N, AVG_M, STD_M, R0, &IDUM, &IDUM2);

	// call this only when initialization is 100% complete
	Solvers solv = Solvers(&system, USE_RK4, USE_LEAPFROG, USE_EULER);

	if (DEBUG)
	{
		for (int i = 0; i < system.n(); i++)
		{
			CelestialBody* cb = system.body(i);
			cout << "mass = " << cb->mass << endl;
			cout << "position = " << *(cb->position) << endl << endl;
		}
	}

	#pragma endregion

	#pragma region Solve and plot

	cout << "E_k after: " << system.EkAvg(false) << endl;
	cout << "E_k after (bound): " << system.EkAvg(true) << endl;
	cout << "E_p after: " << system.EpAvg(false) << endl;
	cout << "E_p after (bound): " << system.EpAvg(true) << endl;

	// this is where the magic happens :)
	system = *(solv.Solve(STEP, PLOT_EVERY));

	if (&system != nullptr)
	{
		cout << "E_k after: " << system.EkAvg(false) << endl;
		cout << "E_k after (bound): " << system.EkAvg(true) << endl;
		cout << "E_p after: " << system.EpAvg(false) << endl;
		cout << "E_p after (bound): " << system.EpAvg(true) << endl;
	}

	getchar(); // pause

	#pragma endregion

	return 0;
}