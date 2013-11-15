#pragma once

#include "stdafx.h"
#include "SolarSystem.h"

using namespace arma;
using namespace std;

class SolarSystem; // forward declaration to avoid circular reference

class CelestialBody
{
public:
	CelestialBody(const string& name, double mass, SolarSystem* system, bool fixed = false);
	CelestialBody(const CelestialBody &cb);
	~CelestialBody(void);
	CelestialBody operator = (const CelestialBody &cb);
	void diff(void);
	string name;
	double mass;
	vec* position;
	vec* velocity;
	vec* force;
	bool fixed;
	mat* plot;
	double Ek(void);

protected:
	int _dim;
	int _currentStep; // used for plotting
public:

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

	// add current position to plot matrix (increments _currentStep afterwards)
	// returns true if room, false if not
	bool plotCurrentPosition(void);

	// plot position matrix to file (after simulation)
	void positionToFile()
	{
		this->plot->save(this->name + ".dat", raw_ascii);
	}
};