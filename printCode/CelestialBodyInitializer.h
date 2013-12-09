// CelestialBodyInitializer.h

#pragma once

#include "CelestialBody.h"
#include "SolarSystem.h"
#include "GaussPDF.h"

// static class to initialize SolarSystem objects with uniform distribution
static class CelestialBodyInitializer
{
public:
	// call this to initialize system with n objects of average mass avgM (std. dev stdM) within radius r0
	static void initialize(SolarSystem* system, int n, double avgM, double stdM, double r0, double timestep);

	// initialize a single particle (CelestialBody)
	static void randInit(SolarSystem* system, const string& name, double avgM, double stdM, double r0, long* idum, long* idum2);
	
	// calculate G in units of t_crunch (after initialization)
	// and units depending on the mass units of the system and r0
	// (here we will use solar masses and ly, respectively)
	double static G(double r0, SolarSystem* system);

	// calculate epsilon in your length unit of choice (see report)
	double static epsilon(double r0, SolarSystem* system, double timestep);

	// calculate the crunch time of the system in years (after initialization)
	double static tCrunch(double r0, SolarSystem* system, double Gyls);
};

