#pragma once

#include "CelestialBody.h"
#include "SolarSystem.h"
#include "GaussPDF.h"

static class CelestialBodyInitializer
{
public:
	static void randInit(SolarSystem* system, const string& name, double avgM, double stdM, double r0, long idum);
	static void initialize(SolarSystem* system, int n, double avgM, double stdM, double r0, long idum);
};

