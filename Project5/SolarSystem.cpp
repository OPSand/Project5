#include "stdafx.h"
#include "SolarSystem.h"

// constructor
SolarSystem::SolarSystem(int dim, int nSteps, int plotEvery, Gravity* grav)
{
	this->_dim = dim;
	this->_nSteps = nSteps;
	this->_plotStep = 0;
	this->_nPlot = (nSteps / plotEvery);
	this->_grav = grav;
	this->_bodies = new vector<CelestialBody*>();
	this->_nBoundPlot = new mat(this->_nPlot, 2);
}

// copy constructor
SolarSystem::SolarSystem(const SolarSystem& other)
{
	this->_dim = other._dim;
	this->_nSteps = other._nSteps;
	this->_plotStep = other._plotStep;
	this->_nPlot = other._nPlot;
	this->_grav = other._grav;
	this->_nBoundPlot = other._nBoundPlot;

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

	delete _nBoundPlot;
}

// operator =
SolarSystem SolarSystem::operator = (const SolarSystem& other)
{
	if( this != &other ) // protect against invalid self-assignment
	{
		this->_dim = other._dim;
		this->_nSteps = other._nSteps;
		this->_plotStep = other._plotStep;
		this->_nPlot = other._nPlot;
		this->_grav = other._grav;
		this->_nBoundPlot = other._nBoundPlot;

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
void SolarSystem::setForces()
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
			
			// use Gravity object to determine force between the pair
			vec F = this->_grav->force(cb_i, cb_j);

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
	cb->setSystem(this); // update cb's system reference
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
bool SolarSystem::plotCurrentStep(bool condition)
{
	bool success = true; // there was room, or the condition wasn't met

	if( condition )
	{
		// calculate potential energy
		this->calculateEp();

		int n = this->n();

		for( int i = 0; i < n; i++ )
		{
			if( ! this->body(i)->plotCurrentPosition() ) // try to plot this element
			{
				success = false; // it is sufficient that one element has no room
			}
			
			// plot # of bound particles (if room)
			if (this->_plotStep < this->_nBoundPlot->n_rows)
			{
				this->_nBoundPlot->at(this->_plotStep, 0) = this->_plotStep;
				this->_nBoundPlot->at(this->_plotStep, 1) = this->nBound();
			}
			else
			{
				success = false;
			}
		}

		// increment step counter
		this->_plotStep++;

		// display progress
		cout << ".";
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

// calculates potential energy of all particles in the system
void SolarSystem::calculateEp()
{
	int n = this->n();

	for (int i = 0; i < n; i++)
	{
		CelestialBody* cb_i = this->body(i);
		double potEnergy = 0;

		// for all other celestial bodies
		for (int j = 0; j < n; j++)
		{
			if (i != j) // don't add potential energy from self
			{
				CelestialBody* cb_j = this->body(j);
				potEnergy += this->_grav->potEnergy(cb_i, cb_j);
			}
		}

		cb_i->Ep = potEnergy;
	}
}

// average potential energy
double SolarSystem::EpAvg(bool boundOnly)
{
	double sum = 0.0;
	int n = this->n();

	for (int i = 0; i < this->n(); i++)
	{
		CelestialBody* cb_i = this->body(i);

		if ((!boundOnly) || (cb_i->isBound()))
		{
			sum += cb_i->Ep;
		}
		else
		{
			n -= 1;
		}
	}

	return (sum / this->n());
}

// average kinetic energy
double SolarSystem::EkAvg(bool boundOnly)
{
	double sum = 0.0;
	int n = this->n();

	for (int i = 0; i < this->n(); i++)
	{
		CelestialBody* cb_i = this->body(i);

		if ((!boundOnly) || (cb_i->isBound()))
		{
			sum += cb_i->Ek();
		}
		else
		{
			n -= 1;
		}
	}

	return (sum / this->n());
}

int SolarSystem::nBound()
{
	int n = 0;

	for (int i = 0; i < this->n(); i++)
	{
		CelestialBody* cb_i = this->body(i);

		if (cb_i->isBound())
		{
			n += 1;
		}
	}

	return n;
}

double SolarSystem::totalMass()
{
	double sum = 0.0;

	for (int i = 0; i < this->n(); i++)
	{
		CelestialBody* cb_i = this->body(i);

		sum += cb_i->mass;
	}

	return sum;
}

double SolarSystem::avgMass()
{
	return (this->totalMass() / this->n());
}

vec SolarSystem::centerOfMass()
{
	vec com = vec(this->_dim);
	com.fill(0.0);

	double totalMass = this->totalMass();

	for (int i = 0; i < this->n(); i++)
	{
		CelestialBody* cb_i = this->body(i);

		com += (cb_i->mass / totalMass) * *(cb_i->position);
	}

	return com;
}

double SolarSystem::distCoM(CelestialBody* cb)
{
	return norm(*(cb->position) - this->centerOfMass(), this->_dim);
}

mat SolarSystem::radialDistribution(double maxR, int boxes, bool boundOnly)
{
	double histogramWidth = (maxR / (double)boxes);
	mat plot = mat(boxes, 2);
	plot.fill(0.0);

	for (int i = 0; i < boxes; i++)
	{
		// use the midpoint of the interval
		plot(i, 0) = (i + 0.5)*histogramWidth;
	}

	for (int i = 0; i < this->n(); i++)
	{
		CelestialBody* cb_i = this->body(i);

		if ((!boundOnly) || (cb_i->isBound()))
		{
			// distance to center of mass
			double comDist = this->distCoM(cb_i);

			// box to put it in
			double box = (comDist / histogramWidth);

			// don't plot if outside max radius
			if (box < boxes)
			{
				plot((int) box, 1) += 1.0; // drop decimal part
			}
		}
	}

	return plot;
}

double SolarSystem::avgDistCoM(bool boundOnly)
{
	int n = this->n();
	double sum = 0.0;

	for (int i = 0; i < this->n(); i++)
	{
		CelestialBody* cb_i = this->body(i);

		if ((!boundOnly) || (cb_i->isBound()))
		{
			sum += distCoM(cb_i);
		}
		else
		{
			n -= 1;
		}
	}

	return (sum / n);
}

double SolarSystem::stdDevDistCoM(bool boundOnly)
{
	int n = this->n();
	double sum = 0.0;

	double avg = this->avgDistCoM(boundOnly);

	for (int i = 0; i < this->n(); i++)
	{
		CelestialBody* cb_i = this->body(i);

		if ((!boundOnly) || (cb_i->isBound()))
		{
			sum += pow((distCoM(cb_i) - avg), 2.0); // deviation squared
		}
		else
		{
			n -= 1;
		}
	}

	return sqrt(sum / n);
}

mat SolarSystem::nBoundPlot()
{
	return *(this->_nBoundPlot);
}