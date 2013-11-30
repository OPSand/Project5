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
	mat* _nBoundPlot;
	int _plotStep; // used for plotting
	vec _com; // total center of mass
	vec _bcom; // center of mass for bound particles

public:
	SolarSystem(int dim, int nSteps, int plotEvery, Gravity* grav);
	SolarSystem(const SolarSystem& other);
	~SolarSystem(void);
	SolarSystem operator =(const SolarSystem& other); 
	void setForces();
	void calculate();
	double EpAvg(bool boundOnly);
	double EkAvg(bool boundOnly);
	int nBound(); // number of bound particles in the system
	double volume(double r);
	double rho(double r);
	double numDens(double r);

	// return dimension of system
	int dim(void)
	{
		return this->_dim;
	}

	// return number of time steos
	int nSteps(void)
	{
		return this->_nSteps;
	}

	// increment the step counter
	void nextStep()
	{
		this->_plotStep++;
	}

	// return number of steps to plot
	int nPlot(void)
	{
		return this->_nPlot;
	}

	// return number of celestial bodies in system
	const int n(void)
	{
		return this->_bodies->size();
	}

	// we plot every ..th step (may differ from the original parameter due to long division)
	int plotEvery()
	{
		return (this->_nSteps / this->_nPlot);
	}

	// which
	int plotStep()
	{
		return this->_plotStep;
	}

	Gravity* grav()
	{
		return this->_grav;
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
	bool plotCurrentStep(bool condition);

	// radial distribution of particles
	mat radialDistribution(double maxR, int boxes, bool boundOnly);

	// # of bound particles per time step
	mat nBoundPlot();

	// total mass of the system
	double totalMass(bool boundOnly);

	// average particle mass
	double avgMass();

	// distance to center of mass
	double distCoM(CelestialBody* cb, bool boundOnly);

	// average distance to center of mass
	double avgDistCoM(bool boundOnly);

	// standard deviation of distance to center of mass
	double stdDevDistCoM(bool boundOnly);

	// coordinates of center of mass
	vec centerOfMass(bool boundOnly);
};