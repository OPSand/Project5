#pragma once

#include "SolarSystem.h"
#include "CelestialBody.h"

class Solvers
{
protected:
	bool _useRK4;
	bool _useLeapfrog;
	bool _useEuler;
	SolarSystem* _system;
	SolarSystem* _rk4;
	SolarSystem* _leapfrog;
	SolarSystem* _euler;
	void RK4(int step);
	void Leapfrog(int step);
	void Euler(int step);
public:
	Solvers(SolarSystem* system, bool useRK4, bool useLeapfrog, bool useEuler);
	~Solvers();
	void Solve(int step, int plotEvery);
};

