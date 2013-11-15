#include "stdafx.h"
#include "SolarSystem.h"

// constructor
SolarSystem::SolarSystem(int dim, int nSteps, int nPlot)
{
	this->_dim = dim;
	this->_nSteps = nSteps;
	this->_nPlot = nPlot;
	this->_bodies = new vector<CelestialBody*>();
}

// copy constructor
SolarSystem::SolarSystem(const SolarSystem& other)
{
	this->_dim = other._dim;
	this->_nSteps = other._nSteps;
	this->_nPlot = other._nPlot;

	// deep copy
	this->_bodies = new vector<CelestialBody*>();
	vector<CelestialBody*> ob = *(other._bodies); // avoid const
	int n = ob.size();
	for( int i = 0; i < n; i++ )
	{
		CelestialBody* cb = new CelestialBody(*(ob.at(i))); // copy
		this->add(cb);
	}
}

// destructor
SolarSystem::~SolarSystem(void)
{
	while( this->n() > 0 ) // counter will be updated automatically
	{
		CelestialBody* cb = this->body(this->n() - 1); // return last element
		this->_bodies->pop_back(); // remove it from vector
		delete cb; // delete it
	}
	delete _bodies; // delete vector
}

// operator =
SolarSystem SolarSystem::operator = (const SolarSystem& other)
{
	if( this != &other ) // protect against invalid self-assignment
	{
		this->_dim = other._dim;
		this->_nSteps = other._nSteps;
		this->_nPlot = other._nPlot;

		// deep copy
		this->_bodies = new vector<CelestialBody*>();
		vector<CelestialBody*> ob = *(other._bodies); // avoid const
		int n = ob.size();
		for( int i = 0; i < n; i++ )
		{
			CelestialBody* cb = new CelestialBody(*(ob.at(i))); // copy
			this->add(cb);
		}
	}

	return *this; // to allow chaining of operators
}

// set force vectors on all elements
void SolarSystem::setForces(void)
{
	int n = this->n();
	for( int i = 0; i < (n - 1); i++ ) // i: 0 -> n-2
	{
		CelestialBody* cb_i = this->body(i);
		
		// reset on first pass (0)
		if( i == 0 )
		{
			cb_i->force->fill(0.0);
		}

		for( int j = (i + 1); j < n; j++ ) // j: i+1 -> n-1
		{
			CelestialBody* cb_j = this->body(j);

			// reset on first pass (1 -> n-1)
			if( i == 0 )
			{
				cb_j->force->fill(0.0);
			}
			
			double dist = cb_i->dist(cb_j); // distance (absolute value)
			vec r = cb_i->position_diff(cb_j); // gives the force the proper direction
			vec F = (cG * cb_i->mass * cb_j->mass / pow(dist, 3.0)) * r; // Newton's law of gravity

			*(cb_i->force) += F; // add force contribution to i
			*(cb_j->force) -= F; // add force contribution to j (Newton's 3rd law)
		}
	}
}

// return a celestial body at the index i
CelestialBody* SolarSystem::body(int i)
{
	return _bodies->at(i);
}

// add a new celestial body to solar system
void SolarSystem::add(CelestialBody* cb)
{
	this->_bodies->push_back(cb); // append to end of vector
}

// return total momentum of system
vec SolarSystem::totalMomentum()
{
	vec mom(this->_dim);
	mom.fill(0.0);

	for(int i = 0; i < this->n(); i++)
	{
		CelestialBody* cb_i = this->body(i);
		if( ! cb_i->fixed ) // fixed bodies never move
		{
			mom += (cb_i->mass * *(cb_i->velocity)); // add p = m*v (non-relativistic)
		}
	}

	return mom;
}

// plots all element positions if the condition is met
// returns true if room, false if not
bool SolarSystem::plotCurrentPositions(bool condition)
{
	bool success = true; // there was room, or the condition wasn't met

	if( condition )
	{
		int n = this->n();

		for( int i = 0; i < n; i++ )
		{
			if( ! this->body(i)->plotCurrentPosition() ) // plot this element
			{
				success = false; // it is sufficient that one element has no room
			}
		}
	}

	return success;
}

// save dimension # i for all elements to "<path>.dat" (rows: time - cols: elements)
void SolarSystem::plotDim(int i, const string& path)
{
	assert(i < this->dim() );

	int n = this->n();
	mat plot(this->nPlot(), n);
	for( int j = 0; j < n; j++ ) // loop through elements
	{
		CelestialBody* cb = this->body(j);
		plot.col(j) = cb->plot->col(i);
	}

	plot.save(path, raw_ascii); // save to file
}

// returns the potential energy of the element cb (if it exists in the system)
double SolarSystem::Ep(CelestialBody* cb)
{
	int n = this->n();

	double potEnergy = 0;

	int i = 0;
	while( this->body(i) != cb )
	{
		i++; // wrong planet, increase index
	}
	
	assert( i < n ); // triggers if cb is not part of the system (then i == n)
	CelestialBody* cb_i = this->body(i);
	
	// for all other celestial bodies
	for( int j = 0; j < n; j++ )
	{
		if( i != j ) // don't add potential energy from self
		{
			CelestialBody* cb_j = this->body(j);
			potEnergy -= (cG * cb_j->mass / cb_i->dist(cb_j)); // potential energy is negative
		}
	}

	return potEnergy;
}