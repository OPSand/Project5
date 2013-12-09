// Gravity.cpp

#include "stdafx.h"
#include "Gravity.h"

// constructor
Gravity::Gravity(double G, double epsilon)
{
	this->_G = G;
	this->_epsilon = epsilon;
}

// copy constructor
Gravity::Gravity(const Gravity& other)
{
	this->_G = other._G;
	this->_epsilon = other._epsilon;
}

Gravity::~Gravity()
{
	// empty
}

// = operator
Gravity Gravity::operator= (const Gravity& other)
{
	if (this != &other) // protect against invalid self-assignment
	{
		this->_G = other._G;
		this->_epsilon = other._epsilon;
	}

	return *this; // to allow chaining of operators
}

// return gravitational force
vec Gravity::force(CelestialBody* cb_i, CelestialBody* cb_j)
{
	double dist = cb_i->dist(cb_j); // distance (absolute value)
	vec r = cb_i->position_diff(cb_j); // gives the force the proper direction

	// Newton's law of gravity (with smoothing factor epsilon)
	vec F = ( this->_G * cb_i->mass * cb_j->mass / ( pow(dist, 3.0) + pow(this->_epsilon, 2.0) ) ) * r;
	return F;
}

// Return potential energy of cb_i with regards to cb_j. NOTE: potential energy is negative
double Gravity::potEnergy(CelestialBody* cb_i, CelestialBody* cb_j)
{
	if (this->_epsilon == 0.0) // no smoothing
	{
		return -(this->_G * cb_i->mass * cb_j->mass / cb_i->dist(cb_j));
	}
	else // smoothed potential
	{
		return ((this->_G * cb_i->mass * cb_j->mass / this->_epsilon)*(atan(cb_i->dist(cb_j) / this->_epsilon) - (cPI / 2.0)));
	}
}

// update the value of G
void Gravity::setG(double G)
{
	this->_G = G;
}

// update the value of epsilon
void Gravity::setEpsilon(double epsilon)
{
	this->_epsilon = epsilon;
}
