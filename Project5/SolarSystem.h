#pragma once

#include "stdafx.h"
#include "CelestialBody.h"

class CelestialBody; // forward declaration to avoid circular reference

class SolarSystem
{
protected:
	int _dim;
	int _nSteps;
	int _nPlot;
	vector<CelestialBody*>* _bodies; // list of celestial bodies in solar system (use pointers to avoid needless copying)
	SolarSystem add(SolarSystem other, bool plus);

public:
	SolarSystem(int dim, int nSteps, int nPlot);
	SolarSystem(const SolarSystem& other);
	~SolarSystem(void);
	SolarSystem operator =(const SolarSystem& other); 
	void setForces(void);
	double Ep(CelestialBody*);

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
};