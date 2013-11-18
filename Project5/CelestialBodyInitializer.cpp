#include "stdafx.h"
#include "CelestialBodyInitializer.h"

void CelestialBodyInitializer::randInit(SolarSystem* system, const string& name, double avgM, double stdM, double r0, long *idum, long* idum2)
{
	// determine mass(avgM, stdM)
	double mass = 0.0;

	// avoid non-positive masses
	while (mass <= 0.0)
	{
		double gauss = GaussPDF::gaussian_deviate(idum);
		mass = avgM + stdM * gauss;
	}

	// construct object & add to system
	CelestialBody* cb = new CelestialBody(name, mass, system);

	// determine position(r0) - independent of system dimension
	int dim = system->dim();

	vec r = vec(dim); // dim-dimensional unit vector
	r.fill(1.0); // we need the while loop to run at least once

	// ensure we are inside the din-dimensional unit sphere so the norm of the position is always less than 1
	while (norm(r, dim) > 1.0)
	{
		// for each dimension
		for (int i = 0; i < dim; i++)
		{
			double uniform = GaussPDF::ran2(idum2);
			r(i) = (2.0 * uniform - 1.0); // between -1.0 and 1.0

			// this required no while loop, but do not use - it favorizes central locations (I think)
			// r(i) = (2.0 * uniform - 1.0) * sqrt(1.0 - norm(r, dim));
		}
	}

	*(cb->position) = r0*r; // scale to correct radius
}

void CelestialBodyInitializer::initialize(SolarSystem* system, int n, double avgM, double stdM, double r0, long *idum, long* idum2)
{
	for (int i = 0; i < n; i++)
	{
		// Generate name
		ostringstream name = ostringstream();
		name << "Star " << i;

		// Initialize CB and add to system
		randInit(system, name.str(), avgM, stdM, r0, idum, idum2); // adds new CB to system
	}
}