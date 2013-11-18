#pragma once

#include "stdafx.h"
#include "CelestialBody.h"
#include "Gravity.h"

class CelestialBody; // forward declaration to avoid circular reference
class Gravity; // forward declaration to avoid circular reference

class SolarSystem
{
protected:
	int _dim;
	int _nSteps;
	int _nPlot;
	Gravity* _grav;
	vector<CelestialBody*>* _bodies; // list of celestial bodies in solar system (use pointers to avoid needless copying)
	SolarSystem add(SolarSystem other, bool plus);
	vec centerOfMass();

public:
	SolarSystem(int dim, int nSteps, int plotEvery, Gravity* grav);
	SolarSystem(const SolarSystem& other);
	~SolarSystem(void);
	SolarSystem operator =(const SolarSystem& other); 
	void setForces();
	double Ep(CelestialBody*);
	double EpAvg(bool boundOnly);
	double EkAvg(bool boundOnly);
	int nBound(); // number of bound particles in the system

	// return dimension of system
	int dim(void)
	{
		return _dim;
	}

	// return number of time steos
	int nSteps(void)
	{
		return _nSteps;
	}

	// return number of steps to plot
	int nPlot(void)
	{
		return _nPlot;
	}

	// return number of celestial bodies in system
	const int n(void)
	{
		return _bodies->size();
	}

	// return a celestial body at the index i
	CelestialBody* body(int i);

	// add a new celecstial body to solar system
	void add(CelestialBody* cb);

	// return total momentum of system
	vec totalMomentum(void);

	// plot dimension # i for all elements to "<path>.dat" (rows: time - cols: elements)
	void plotDim(int i, const string& path);

	// plots all element positions if the condition is met
	// returns true if room, false if not
	bool plotCurrentPositions(bool condition);

	// radial distribution of particles
	mat radialDistribution(double maxR, int boxes, bool boundOnly);

	// total mass of the system
	double totalMass();

	// distance to center of mass
	double distCoM(CelestialBody* cb);
};