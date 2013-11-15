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

	// initialization
	const int N = 20; // number of celestial bodies
	const double R0 = 20.0; // ly
	const double AVG_M = 10.0; // solar masses
	const double STD_M = 1.0;

	// time steps
	const int STEP = 1; 
	const int N_STEPS = 1000; // number of steps total
	const int PLOT_EVERY = 1; // plot every ...th step
	const int N_PLOT = (N_STEPS / PLOT_EVERY); // how many steps we actually plot

	// flags
	const bool USE_LEAPFROG = false; // use Leapfrog method
	const bool USE_RK4 = true; // use Runge-Kutta method
	const bool USE_EULER = true; // use Euler-Cromer method

	#pragma endregion

	#pragma region Initialization

	// create system
	SolarSystem system = SolarSystem(DIM, N_STEPS, N_PLOT);

	// add N randomly initialized celestial bodies
	CelestialBodyInitializer::initialize(&system, N, AVG_M, STD_M, R0);

	// call this only when initialization is 100% complete
	Solvers solv = Solvers(&system, USE_RK4, USE_LEAPFROG, USE_EULER);

	#pragma endregion

	#pragma region Solve and plot

	solv.Solve(STEP, PLOT_EVERY); // this is where the magic happens :)

	getchar(); // pause

	#pragma endregion

	return 0;
}