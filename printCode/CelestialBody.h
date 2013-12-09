#pragma once

#include "stdafx.h"
#include "SolarSystem.h"

using namespace arma;
using namespace std;

class SolarSystem; // forward declaration to avoid circular reference

// class that simulates a particle in an N-body simulation
// must belong to a SolarSystem object
class CelestialBody
{
protected: 
	int _dim; // number of dimensions (set by _system)
	SolarSystem* _system; // reference to parent object

public:
	// construction, copying, destruction
	CelestialBody(const string& name, double mass, SolarSystem* system, bool fixed = false);
	CelestialBody(const CelestialBody &cb);
	~CelestialBody(void);
	CelestialBody operator = (const CelestialBody &cb);

	// public parameters
	string name;
	double mass;
	vec* position; // _dim-dimensional vector
	vec* velocity;
	vec* force;
	bool fixed; // flag: never update position if true
	mat* plot; // _dim columns, _system->nPlot() rows
	double Ep; // potential energy, calculated by _system

public:
	// kinetic energy
	double Ek(void);

	// returns the acceleration when the force is set
	vec acc(void)
	{
		return (*(this->force)/this->mass);
	}

	// returns the position of cb relative to this in vector form
	vec position_diff(CelestialBody* cb)
	{
		return (*(cb->position) - *(this->position));
	}

	// returns the distance between this and cb as a scalar
	double dist(CelestialBody* cb)
	{
		return norm(this->position_diff(cb), this->_dim);
	}

	// add current position to plot matrix
	// returns true if room, false if not
	bool plotCurrentPosition(void);

	// plot position matrix to file (after simulation)
	void positionToFile()
	{
		this->plot->save(this->name + ".dat", raw_ascii);
	}

	// true if total energy is negative, false otherwise
	bool isBound();

	// set system reference
	void setSystem(SolarSystem* system)
	{
		this->_system = system;
	}
};