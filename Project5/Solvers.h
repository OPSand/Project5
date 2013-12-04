#pragma once

#include "SolarSystem.h"
#include "CelestialBody.h"

// class to solve the equation of motion for a SolarSystem N-body simulation
// using the Runge-Kutta (4th order), Leapfrog and/or Euler-Cromer algoritms
class Solvers
{
protected:
	// flags
	bool _useRK4;
	bool _useLeapfrog;
	bool _useEuler;

	// references to SolarSystem objects to iterate on
	SolarSystem* _system;
	SolarSystem* _rk4;
	SolarSystem* _leapfrog;
	SolarSystem* _euler;

	// individual step for each of the supported algoritms
	void RK4(double step);
	void Leapfrog(double step);
	void Euler(double step);

public:

	// constructors/destructors
	Solvers(SolarSystem* system, bool useRK4, bool useLeapfrog, bool useEuler);
	~Solvers();

	// call this to solve the equations and save results to files
	// the system will determine the number of steps; all we need is to provide a step length
	// returns the system in the state the Leapfrog algoritm leaves it in (or nullptr if that algorithm is not used)
	SolarSystem* Solve(double step);
};

