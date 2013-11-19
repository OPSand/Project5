#pragma once

#include "CelestialBody.h"
#include "SolarSystem.h"
#include "GaussPDF.h"

static class CelestialBodyInitializer
{
public:
	static void randInit(SolarSystem* system, const string& name, double avgM, double stdM, double r0, long* idum, long* idum2);
	static void initialize(SolarSystem* system, int n, double avgM, double stdM, double r0);
	double static volume(double r, int dim);
	double static rho(double r, double avgM, double n, int dim);
	double static G(double r0, double avgM, double n, int dim);
	double static tCrunch(double r0, double avgM, double n, double Gyls, int dim);
};

