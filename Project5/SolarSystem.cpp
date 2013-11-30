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
		// calculate potential energy ++
		this->calculate();

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

// calculates potential energy & total/unbound centers of mass
void SolarSystem::calculate()
{
	int n = this->n();

	// reset potential energies
	for (int i = 0; i < n; i++)
	{
		CelestialBody* cb_i = this->body(i);
		cb_i->Ep = 0.0;
	}

	// this loop MUST complete before bound / unbound is evaluated!
	for (int i = 0; i < n; i++)
	{
		CelestialBody* cb_i = this->body(i);

		// for all other celestial bodies not yet iterated over
		for (int j = (i+1); j < n; j++)
		{
			CelestialBody* cb_j = this->body(j);
			cb_i->Ep += this->_grav->potEnergy(cb_i, cb_j);
			cb_j->Ep += this->_grav->potEnergy(cb_i, cb_j); // add to both
		}
	}

	// now that potential energy is calculated for all particles
	// we can finally start evaluating bound / unbound

	// initialize center of mass vectors
	this->_com = vec(this->dim()); this->_com.fill(0.0);
	this->_bcom = vec(this->dim()); this->_bcom.fill(0.0);

	// total masses (all + bound)
	double totalMass = this->totalMass(false);
	double boundMass = this->totalMass(true);

	for (int i = 0; i < n; i++)
	{
		CelestialBody* cb_i = this->body(i);

		// contribution to center of mass
		this->_com += (cb_i->mass / totalMass) * *(cb_i->position);

		// contribution to bound center of mass
		if (cb_i->isBound())
		{
			this->_bcom += (cb_i->mass / boundMass) * *(cb_i->position);
		}
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

	return (sum / (2.0 *this->n() ));
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

double SolarSystem::totalMass(bool boundOnly)
{
	double sum = 0.0;

	for (int i = 0; i < this->n(); i++)
	{
		CelestialBody* cb_i = this->body(i);

		if ((!boundOnly) || (cb_i->isBound()))
		{
			sum += cb_i->mass;
		}
	}

	return sum;
}

double SolarSystem::avgMass()
{
	return (this->totalMass(false) / this->n());
}

vec SolarSystem::centerOfMass(bool boundOnly)
{
	if (boundOnly)
	{
		return this->_bcom;
	}
	else
	{
		return this->_com;
	}
}

double SolarSystem::distCoM(CelestialBody* cb, bool boundOnly)
{
	return norm(*(cb->position) - this->centerOfMass(boundOnly), this->_dim);
}

// calculate the volume of an n-sphere
// https://en.wikipedia.org/wiki/N-sphere#Closed_forms
double SolarSystem::volume(double r)
{
	double k = ((double)this->dim() / 2.0); // dim = 2k
	double unitSphereVolume = (pow(cPI, k) / tgamma(k + 1.0)); // tgamma is the gamma function
	return pow(r, this->dim())*unitSphereVolume;
}

// mass density within radius r
double SolarSystem::rho(double r)
{
	double V = this->volume(r); // volume [ly^3]
	return (this->totalMass(false) / V); // mass density [solar masses / ly^3]
}

// mass density within radius r
double SolarSystem::numDens(double r)
{
	double V = this->volume(r); // volume [ly^3]
	return (this->n() / V); // mass density [solar masses / ly^3]
}

mat SolarSystem::radialDistribution(double maxR, int boxes, bool boundOnly)
{
	double histogramWidth = (maxR / (double)boxes);
	mat plot = mat(boxes, 3);
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
			double comDist = this->distCoM(cb_i, boundOnly);

			// box to put it in
			double box = (comDist / histogramWidth);

			// don't plot if outside max radius
			if (box < boxes)
			{
				plot((int) box, 1) += 1.0; // drop decimal part
			}
		}
	}

	// second column is N(r)
	// third column: n(r) = N(r)/V(r)

	for (int i = 0; i < boxes; i++)
	{
		double r = plot(i, 0);
		double delta_r = 0.5*histogramWidth;
		double Vmin = this->volume(r - delta_r);
		double Vmax = this->volume(r + delta_r);
		double V = (Vmax - Vmin);

		assert(V > 0.0);

		plot(i, 2) = plot(i, 1) / V;
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
			sum += distCoM(cb_i, boundOnly);
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
			sum += pow((distCoM(cb_i, boundOnly) - avg), 2.0); // deviation squared
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