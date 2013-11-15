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

	// randomization seed (must be negative)
	const long IDUM = -1337;

	// physical constants
	const double LY = 9.4607e15; // m
	const double MYR = cYr * 1.0e6; // s
	const double M_SUN = 1.9891e30; // kg
	const double G = cG * M_SUN * pow(MYR, 2.0) / pow(LY, 3.0); // G in Myr, ly, solar masses
	const double EPSILON = 0.0; // correction to Newton in ly to avoid infinite forces at close range

	// initialization
	const int N = 100; // number of celestial bodies
	const double R0 = 20.0; // ly
	const double AVG_M = 10.0; // solar masses
	const double STD_M = 1.0; // solar masses

	// time steps
	const double STEP = 1.0; // step size (Myr)
	const int N_STEPS = 1000; // number of steps total
	const int PLOT_EVERY = 1; // plot every ...th step

	// flags
	const bool USE_LEAPFROG = true; // use Leapfrog method
	const bool USE_RK4 = true; // use Runge-Kutta method
	const bool USE_EULER = true; // use Euler-Cromer method

	#pragma endregion

	#pragma region Initialization

	// create gravity
	Gravity g = Gravity(G, EPSILON);

	// create system
	SolarSystem system = SolarSystem(DIM, N_STEPS, PLOT_EVERY, &g);

	// add N randomly initialized celestial bodies
	CelestialBodyInitializer::initialize(&system, N, AVG_M, STD_M, R0, IDUM);

	// call this only when initialization is 100% complete
	Solvers solv = Solvers(&system, USE_RK4, USE_LEAPFROG, USE_EULER);

	#pragma endregion

	#pragma region Solve and plot

	// this is where the magic happens :)
	solv.Solve(STEP, PLOT_EVERY);

	getchar(); // pause

	#pragma endregion

	return 0;
}