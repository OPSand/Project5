#include "stdafx.h"
#include "CelestialBody.h"
#include "SolarSystem.h"

using namespace std;
using namespace arma;

// constructor (system passed by reference so we do not copy it)
CelestialBody::CelestialBody(const string& name, double mass, SolarSystem* system, bool fixed)
{
	// set protected member variables
	this->name = name;
	this->mass = mass;
	this->_dim = system->dim();
	this->_currentStep = 0;

	// initialize vectors with correct dimension
	this->position = new vec(this->_dim);
	this->velocity = new vec(this->_dim);
	this->force = new vec(this->_dim);

	// default values
	this->fixed = fixed;
	this->position->fill(0);
	this->velocity->fill(0);
	this->force->fill(0);

	// plot matrix
	this->plot = new mat(system->nPlot(), this->_dim);
	plot->fill(0.0);

	// add to system
	system->add(this);
}

// copy constructor
CelestialBody::CelestialBody(const CelestialBody &cb)
{
	this->name = cb.name; // unlikely to change
	this->mass = cb.mass;
	this->_dim = cb._dim;
	this->_currentStep = cb._currentStep;

	this->fixed = cb.fixed;

	// These may be changed, so we copy them
	this->position = new vec(*cb.position);
	this->velocity = new vec(*cb.velocity);
	this->force = new vec(*cb.force);
	this->plot = new mat(*cb.plot);

	// NOTE: we do not copy system
}

// destructor
CelestialBody::~CelestialBody(void)
{
	delete this->position;
	delete this->velocity;
	delete this->force;
	delete this->plot;
}

// operator =
CelestialBody CelestialBody::operator = (const CelestialBody &cb)
{
	if( this != &cb ) // protect against invalid self-assignment
	{
		this->name = cb.name; // unlikely to change
		this->mass = cb.mass;
		this->_dim = cb._dim;
		this->_currentStep = cb._currentStep;

		this->fixed = cb.fixed;

		// These may be changed, so we copy them
		this->position = new vec(*cb.position);
		this->velocity = new vec(*cb.velocity);
		this->force = new vec(*cb.force);
		this->plot = new mat(*cb.plot);

		// NOTE: we do not copy system
	}

	return *this; // to allow operator chaining: a = b = c
}

// add current position to plot matrix (increments _currentStep afterwards)
// returns true if room, false if not
bool CelestialBody::plotCurrentPosition()
{
	if( this->_currentStep < this->plot->n_rows )
	{
		for( int j = 0; j < this->plot->n_cols; j++ )
		{
			this->plot->at(this->_currentStep, j) = this->position->at(j);
		}

		this->_currentStep++;
		return true;
	}
	else // no more room in matrix
	{
		return false;
	}
}

// differentiate (changes the current object)
void CelestialBody::diff()
{
	if( ! this->fixed ) // moving body
	{
		this->position = new vec(*this->velocity);
		this->velocity = new vec(this->acc());
	}
	else // a fixed body will never move
	{
		this->position->fill(0.0);
		this->velocity->fill(0.0);
	}
}

// kinetic energy
double CelestialBody::Ek()
{
	double v2 = dot(*(this->velocity), *(this->velocity)); // inner product equals v^2
	return (0.5 * this->mass * v2);
}