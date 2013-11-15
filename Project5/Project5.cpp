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
	const double EPSILON = 0.15; // correction to Newton in ly to avoid infinite forces at close range

	// calculate time scale & G in correct units
	const double V0 = (4.0 / 3.0)* cPI *pow(R0, 3.0); // initial volume
	const double RHO0 = N * AVG_M / V0; // initial mass density
	const double T_CRUNCH = sqrt(3.0 * cPI / (32.0 * G_YLS * RHO0)); // crunch time in years
	const double G = pow(cPI, 2.0) * pow(R0, 3.0) / (8 * pow(T_CRUNCH, 2.0) * AVG_M * N); // G in t_crunch, ly, solar masses

	// time steps
	const double CRUNCH_TIMES = 5.0; // # of crunch times to simulate for
	const int N_STEPS = 1000; // number of steps total
	const double STEP = 1.0; // step size (in crunch times)
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