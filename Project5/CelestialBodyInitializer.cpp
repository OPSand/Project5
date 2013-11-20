#include "stdafx.h"
#include "CelestialBodyInitializer.h"

void CelestialBodyInitializer::randInit(SolarSystem* system, const string& name, double avgM, double stdM, double r0, long* idum, long* idum2)
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

	// Create uniform distribution in the n-ball using the recipe at:
	// https://en.wikipedia.org/wiki/Hypersphere#Uniformly_at_random_from_the_.28n.C2.A0.E2.88.92.C2.A01.29-sphere

	vec x = vec(dim);
	
	// for each dimension
	for (int i = 0; i < dim; i++)
	{
		x(i) = GaussPDF::gaussian_deviate(idum);
	}
	
	double r = norm(x, dim);
	x = (1.0 / r)*x; // x is now a unit vector - its possible directions are now uniformly distributed on the surface!

	r = GaussPDF::ran2(idum2); // uniformly distributed between 0.0 and 1.0
	x = pow(r, (1.0 / dim)) * x; // uniformly distributed in the unit n-ball

	*(cb->position) = r0*x; // scale to correct radius
}

void CelestialBodyInitializer::initialize(SolarSystem* system, int n, double avgM, double stdM, double r0)
{
	// generate random (negative) seeds based on clock
	long* idum = new long(-time(NULL));
	long* idum2 = new long(-time(NULL));

	for (int i = 0; i < n; i++)
	{
		// Generate name
		ostringstream name = ostringstream();
		name << "Star " << i;

		// Initialize CB and add to system
		randInit(system, name.str(), avgM, stdM, r0, idum, idum2); // adds new CB to system
	}

	// finally, update Gravity using average mass (calculated, not avgM)
	// TODO: Make dimensionless? https://en.wikipedia.org/wiki/N-sphere#Closed_forms
	system->grav()->setG(G(r0, system->avgMass(), system->n(), system->dim()));

	// garbage collection
	delete idum;
	delete idum2;
}

// calculate the volume of an n-sphere
// https://en.wikipedia.org/wiki/N-sphere#Closed_forms
double CelestialBodyInitializer::volume(double r, int dim)
{
	double k = (double)dim / 2.0; // dim = 2k
	double unitSphereVolume = (pow(cPI, k) / tgamma(k + 1.0)); // tgamma is the gamma function
	return pow(r, dim)*unitSphereVolume;
}

// calculate mass density in solar masses per cubic light yeats
double CelestialBodyInitializer::rho(double r, double avgM, double n, int dim)
{
	double V = volume(r, dim); // initial volume [ly^3]
	return (n * avgM / V); // mass density [solar masses / ly^3]
}

// calculate G in units of t_crunch, ly, solar masses
double CelestialBodyInitializer::G(double r0, double avgM, double n, int dim)
{
	return (3.0 * cPI / (32.0 * rho(r0, avgM, n, dim))); // G in t_crunch, ly, solar masses
}

// calculate t_crunch in years
double CelestialBodyInitializer::tCrunch(double r0, double avgM, double n, double Gyls, int dim)
{
	return sqrt(3.0 * cPI / (32.0 * Gyls * rho(r0, avgM, n, dim))); // crunch time [years]
}