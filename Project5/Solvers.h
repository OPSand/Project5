#pragma once

#include "SolarSystem.h"
#include "CelestialBody.h"

class Solvers
{
protected:
	bool _useRK4;
	bool _useLeapfrog;
	bool _useEuler;
	const SolarSystem* _system;
	SolarSystem* _rk4;
	SolarSystem* _leapfrog;
	SolarSystem* _euler;
	void RK4(int step, int nSteps, int plotEvery, int nPlot);
	void Leapfrog(int step, int nSteps, int plotEvery, int nPlot);
	void Euler(int step, int nSteps, int plotEvery, int nPlot);
public:
	Solvers(const SolarSystem* system, bool useRK4, bool useLeapfrog, bool useEuler);
	~Solvers();
	void Solve(int step, int nSteps, int plotEvery);
};

