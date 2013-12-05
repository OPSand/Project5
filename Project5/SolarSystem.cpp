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

// call this every step to update force vectors on all elements based on current positions
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

	// for every body
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
		// calculate potential energy ++ before plotting
		this->calculate();

		int n = this->n();

		// for every element
		for( int i = 0; i < n; i++ )
		{
			if( ! this->body(i)->plotCurrentPosition() ) // try to plot this element
			{
				success = false; // it is sufficient that one element has no room to declare failure
			}
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

		// increment plot step counter
		this->_plotStep++;

		// display progress
		cout << ".";
	}

	return success;
}

// save dimension # i for all elements to "<path>.dat" (rows: time - cols: elements)
void SolarSystem::plotDim(int i, const string& path)
{
	// check that this is a dimension we're actually simulating
	assert(i < this->dim() );

	int n = this->n();
	mat plot(this->nPlot(), n);

	// loop through elements
	for( int j = 0; j < n; j++ )
	{
		CelestialBody* cb = this->body(j);
		plot.col(j) = cb->plot->col(i);
	}

	plot.save(path, raw_ascii); // save to file
}

// call this every step to calculate potential energy & total/unbound centers of mass
// (determining which elements are bound and their radial distribution)
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

	// for all elements
	for (int i = 0; i < n; i++)
	{
		CelestialBody* cb_i = this->body(i);

		// add contribution to center of mass
		this->_com += (cb_i->mass / totalMass) * *(cb_i->position);

		if (cb_i->isBound()) // bound particle
		{
			// add contribution to bound center of mass
			this->_bcom += (cb_i->mass / boundMass) * *(cb_i->position);
		}
	}
}

// average potential energy
double SolarSystem::EpAvg(bool boundOnly)
{
	double sum = 0.0;
	int n = this->n();

	// for all particles
	for (int i = 0; i < this->n(); i++)
	{
		CelestialBody* cb_i = this->body(i);

		// check bound conditions (if they apply)
		if ((!boundOnly) || (cb_i->isBound()))
		{
			// count contribution from this particle
			sum += cb_i->Ep;
		}
		else // unbound (in the case that we want bound only)
		{
			// don't count this particle
			n -= 1;
		}
	}

	// NOTE: divide by 2 so we do not count each vertex twice!
	// (each Ep is determined by PAIRS of bodies)
	return (sum / (2.0 * this->n() ));
}

// average kinetic energy
double SolarSystem::EkAvg(bool boundOnly)
{
	double sum = 0.0;
	int n = this->n();

	// for all particles
	for (int i = 0; i < this->n(); i++)
	{
		CelestialBody* cb_i = this->body(i);

		// check bound conditions (if they apply)
		if ((!boundOnly) || (cb_i->isBound()))
		{
			// count contribution from this particle
			sum += cb_i->Ek();
		}
		else // unbound (in the case that we want bound only)
		{
			// don't count this particle
			n -= 1;
		}
	}

	return (sum / this->n());
}

// return number of bound particles at the current step
int SolarSystem::nBound()
{
	int n = 0; // bound counter

	// for all particles
	for (int i = 0; i < this->n(); i++)
	{
		CelestialBody* cb_i = this->body(i);

		// count if bound
		if (cb_i->isBound())
		{
			n += 1;
		}
	}

	return n;
}

// returns total mass of the system (or its bound particles)
double SolarSystem::totalMass(bool boundOnly)
{
	double sum = 0.0;

	// for all particles
	for (int i = 0; i < this->n(); i++)
	{
		CelestialBody* cb_i = this->body(i);

		// check bound conditions (if they apply)
		if ((!boundOnly) || (cb_i->isBound()))
		{
			// count contribution from this particle
			sum += cb_i->mass;
		}
	}

	return sum;
}

// return average mass of the particles in the system
double SolarSystem::avgMass()
{
	return (this->totalMass(false) / this->n());
}

// return total/bound-particle center of mass
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

// return a particle's distance from total/bound-particle center of mass
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

// number of (bound) particles within radius r
// from the (bound) center of mass
int SolarSystem::nWithin(double r, bool boundOnly)
{
	int sum = 0; // count particles within r

	// for every particle
	for (int i = 0; i < this->n(); i++)
	{
		CelestialBody* cb_i = this->body(i);

		// check bound conditions (if they apply)
		if ((!boundOnly) || (cb_i->isBound()))
		{
			// check distance to the appropriate center of mass
			if (distCoM(cb_i, boundOnly) < r)
			{
				sum++;
			}
		}
	}

	return sum;
}

// radial distribution of particles within radius maxR
// avgBin: average number of particles per bin
mat SolarSystem::radialDistribution(double maxR, double avgBin, bool boundOnly)
{
	// determine number of bins from average and counted number of particles
	// within maxR
	int bins = ((double)this->nWithin(maxR, boundOnly) / avgBin);
	bins++; // we just threw away the decimal part, so add 1 to round up

	// use number of bins to determine width of each bin
	double histogramWidth = (maxR / (double)bins);
	mat plot = mat(bins, 3);
	plot.fill(0.0);

	// set r values
	for (int i = 0; i < bins; i++)
	{
		// use the midpoint of the interval
		plot(i, 0) = (i + 0.5)*histogramWidth;
	}

	// loop through particles and sort into appropriate bins
	for (int i = 0; i < this->n(); i++)
	{
		CelestialBody* cb_i = this->body(i);

		// check bound conditions (if they apply)
		if ((!boundOnly) || (cb_i->isBound()))
		{
			// distance to center of mass
			double comDist = this->distCoM(cb_i, boundOnly);

			// box to put it in
			double bin = (comDist / histogramWidth);

			// don't plot if outside max radius
			if (bin < bins)
			{
				plot((int) bin, 1) += 1.0;
			}
		}
	}

	// second column is N(r), the total number of particles
	// third column is n(r) = N(r)/V(r), the number density

	// for every bin, divide N(r) by the volume V(r)
	for (int i = 0; i < bins; i++)
	{
		double r = plot(i, 0);
		double delta_r = 0.5 * histogramWidth;
		double Vmin = this->volume(r - delta_r);
		double Vmax = this->volume(r + delta_r);
		double V = (Vmax - Vmin);

		assert(V > 0.0);

		plot(i, 2) = plot(i, 1) / V;
	}

	return plot;
}

// average distance to total/bound center of mass
double SolarSystem::avgDistCoM(bool boundOnly)
{
	int n = this->n();
	double sum = 0.0;

	// for every particle
	for (int i = 0; i < this->n(); i++)
	{
		CelestialBody* cb_i = this->body(i);

		// check bound conditions (if they apply)
		if ((!boundOnly) || (cb_i->isBound()))
		{
			// count contribution from this particle
			sum += distCoM(cb_i, boundOnly);
		}
		else // unbound (in the case that we want bound only)
		{
			// don't count this particle
			n -= 1;
		}
	}

	return (sum / n);
}

// standard deviation of distance to total/bound center of mass
double SolarSystem::stdDevDistCoM(bool boundOnly)
{
	int n = this->n();
	double sum = 0.0;

	// we need the average to calculate the variance
	double avg = this->avgDistCoM(boundOnly);

	// for every particle
	for (int i = 0; i < this->n(); i++)
	{
		CelestialBody* cb_i = this->body(i);

		// check bound conditions (if they apply)
		if ((!boundOnly) || (cb_i->isBound()))
		{
			// count contribution from this particle
			sum += pow((distCoM(cb_i, boundOnly) - avg), 2.0); // deviation squared
		}
		else // unbound (in the case that we want bound only)
		{
			// don't count this particle
			n -= 1;
		}
	}

	// sum is now the variance, return std.dev. instead
	return sqrt(sum / n);
}

// returns a matrix that can be used for plotting number of bound particles
// for every step
mat SolarSystem::nBoundPlot()
{
	return *(this->_nBoundPlot);
}

// average MINIMUM distance between a pair of particles
// O(n^2), so do not call this every step!
double SolarSystem::avgMinDist()
{
	double sum = 0.0;
	int n = this->n();

	// for each pair of particles
	for (int i = 0; i < n; i++)
	{
		CelestialBody* cb_i = this->body(i);
		double min = INFINITY;

		for (int j = (i + 1); j < n; j++)
		{
			CelestialBody* cb_j = this->body(j);

			// calculate distance
			double dist = cb_i->dist(cb_j);

			// is this the minimal distance encountered so far
			if (dist < min)
			{
				min = dist;
			}
		}

		sum += min;
	}

	// return average
	return (sum / n);
}