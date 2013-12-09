#pragma once

#include "CelestialBody.h"

class CelestialBody; // forward declaration to avoid circular reference

// class that simulates the gravity in a SolarSystem
// must be able to return a force and a potential energy
static class Gravity
{
protected:
	double _G; // gravitational constant
	double _epsilon; // smoothing term
public:
	// constructors, copying and destructors
	Gravity(double G, double epsilon);
	Gravity(const Gravity& other);
	~Gravity();
	Gravity operator= (const Gravity& other);

	// return gravitational constant
	double G()
	{
		return this->_G;
	}

	// return smoothing (epsilon)
	double epsilon()
	{
		return this->_epsilon;
	}

	// ability to update G and Epsilon after construction
	void setG(double G);
	void setEpsilon(double epsilon);

	// return gravitational force between two objects
	vec force(CelestialBody* cb_i, CelestialBody* cb_j);

	// return potential energy between two objects
	double potEnergy(CelestialBody* cb_i, CelestialBody* cb_j);
};

