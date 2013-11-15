// Project5.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include"armadillo"
#include<list>

using namespace arma;
using namespace std;

// Entry point for console application
int _tmain(int argc, _TCHAR* argv[])
{
	#pragma region Flags and settings

	// dimensions
	const int DIM = 2;

	// number of celestial bodies
	const int N = 2;

	// time steps
	const int N_STEPS = 1000; // number of steps total
	const int PLOT_EVERY = 1; // plot every ...th step
	const int N_PLOT = (N_STEPS / PLOT_EVERY); // how many steps we actually plot

	// flags
	const bool USE_LEAPFROG = false; // use Leapfrog method
	const bool USE_RK4 = true; // use Runge-Kutta method
	const bool DEBUG = false; // use for debugging only

	#pragma endregion

	getchar(); // pause

	return 0;
}