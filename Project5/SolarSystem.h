#pragma once

#include "stdafx.h"
#include "CelestialBody.h"
#include "Gravity.h"

class CelestialBody; // forward declaration to avoid circular reference
class Gravity; // forward declaration to avoid circular reference

// class for running an N-body simulation in any dimension
// (numerical solving of the equations of motion must be handled by another class, such as Solvers)
class SolarSystem
{
protected:
	int _dim; // dimensions
	int _nSteps; // number of steps to simulate
	int _nPlot; // number of steps to plot
	Gravity* _grav; // gravity object (provides forces & potential energies)
	vector<CelestialBody*>* _bodies; // list of celestial bodies in solar system (use pointers to avoid needless copying)
	mat* _nBoundPlot; // used for plotting number of bound particles
	int _plotStep; // used for determining which steps to plot
	vec _com; // total center of mass
	vec _bcom; // center of mass for bound particles

public:

	// constructors, destructors, copying
	SolarSystem(int dim, int nSteps, int plotEvery, Gravity* grav);
	SolarSystem(const SolarSystem& other);
	~SolarSystem(void);
	SolarSystem operator =(const SolarSystem& other); 

	// call this every step to update force vectors on all elements based on current positions
	void setForces();

	// call this every step to calculate potential energy & total/unbound centers of mass
	// (determining which elements are bound and their radial distribution)
	void calculate();

	// average kinetic and potential energies (for the virial theorem)
	double EpAvg(bool boundOnly);
	double EkAvg(bool boundOnly);

	int nBound(); // number of bound particles in the system

	// calculate system parameters within a radius r
	double volume(double r);
	double rho(double r); // mass density
	double numDens(double r); // number density

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

	// which step are we currenty at?
	int plotStep()
	{
		return this->_plotStep;
	}

	// return reference to gravity object
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

	// number of (bound) particles within radius r
	// from the (bound) center of mass
	int nWithin(double r, bool boundOnly);

	// radial distribution of particles within radius maxR
	// avgBin: average number of particles per bin
	mat radialDistribution(double maxR, double avgBin, bool boundOnly);

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

	// average MINIMUM distance between a pair of particles
	// O(n^2), so do not call this every step!
	double avgMinDist();
};