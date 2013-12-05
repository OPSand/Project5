#include "stdafx.h"
#include "Solvers.h"

// constructor (isim = number of sim in a series)
Solvers::Solvers(SolarSystem* system, const string& id, bool useRK4, bool useLeapfrog, bool useEuler)
{
	// save pointer to original system
	this->_system = system;

	// set ID
	this->_id = id;

	// set flags
	this->_useRK4 = useRK4;
	this->_useLeapfrog = useLeapfrog;
	this->_useEuler = useEuler;

	// calculate potential energy ++ before copying
	this->_system->calculate();

	// create system copies as necessary
	if (this->_useRK4)
	{
		this->_rk4 = new SolarSystem(*system); // deep copy
	}

	if (this->_useLeapfrog)
	{
		this->_leapfrog = new SolarSystem(*system); // deep copy
	}

	if (this->_useEuler)
	{
		this->_euler = new SolarSystem(*system); // deep copy
	}
}

// destructor
Solvers::~Solvers()
{
	delete this->_rk4;
	delete this->_leapfrog;
	delete this->_euler;
}

// call this to solve the equations and save results to files
// the system will determine the number of steps; all we need is to provide a step length
// returns the system in the state the Leapfrog algoritm leaves it in (or nullptr if that algorithm is not used)
SolarSystem* Solvers::Solve(double step)
{
	clock_t start, finish;
	start = clock();
	double elapsedTime = 0.0;
	int nSize = this->_system->nSteps();
	vec  tabElapsedTime = vec(nSize);
	for (int i = 0; i < this->_system->nSteps(); i++) // for each time step
	{
		if (this->_useRK4)
		{
			this->_rk4->plotCurrentStep(i % this->_rk4->plotEvery() == 0); // if we want to plot this step, do it
			RK4(step); // perform step
		}

		if (this->_useLeapfrog)
		{
			this->_leapfrog->plotCurrentStep(i % this->_leapfrog->plotEvery() == 0); // if we want to plot this step, do its
			Leapfrog(step); // perform step
		}

		if (this->_useEuler)
		{
			this->_euler->plotCurrentStep(i % this->_euler->plotEvery() == 0); // if we want to plot this step, do it
			Euler(step); // perform step
		}
	}


	// Timing part: End ...
	finish = clock();
	elapsedTime = (double)(finish - start) / CLOCKS_PER_SEC; // To convert this into seconds !

	if (this->_useRK4)
	{
		// plot to file (independent of dimension)
		for (int i = 0; i < this->_rk4->dim(); i++)
		{
			ostringstream fname = ostringstream();
			fname << "sim_" << this->_id << "_pos" << i << "_rk4.dat";
			this->_rk4->plotDim(i, fname.str());
		}
	}

	if (this->_useLeapfrog)
	{
		// plot to file (independent of dimension)
		for (int i = 0; i < this->_leapfrog->dim(); i++)
		{
			ostringstream fname = ostringstream();
			fname << "sim_" << this->_id << "_pos" << i << "_leapfrog.dat";
			this->_leapfrog->plotDim(i, fname.str());
		}
	}

	if (this->_useEuler)
	{
		// plot to file (independent of dimension)
		for (int i = 0; i < this->_euler->dim(); i++)
		{
			ostringstream fname = ostringstream();
			fname << "sim_" << this->_id << "_pos" << i << "_euler.dat";
			this->_euler->plotDim(i, fname.str());
		}
	}

	cout << "DONE!" << endl << endl;

	// And siplay of the timing part 
	printf("Elapsed Time : %lf \n", elapsedTime);

	// Return calculation of the average time spent for a step of leapfrog/RG4
	//avTime(tabElapsedTime);
	// we will use leapfrog to return results by default
	SolarSystem* returnSystem = this->_leapfrog;

	if (returnSystem != nullptr)
	{
		// make sure potential energy ++ is up to date
		returnSystem->calculate();

		return returnSystem;
	}
	else // return nothing
	{
		return nullptr;
	}
}

// a single Runge-Kutta step (parameter = step length)
void Solvers::RK4(double step)
{
	int n = this->_rk4->n();

	// matrices for Runge-Kutta values:
	mat k1_v(this->_rk4->dim(), n); // k1, velocity
	mat k1_p(this->_rk4->dim(), n); // k1, position
	mat k2_v(this->_rk4->dim(), n); // k2, velocity (etc.)
	mat k2_p(this->_rk4->dim(), n);
	mat k3_v(this->_rk4->dim(), n);
	mat k3_p(this->_rk4->dim(), n);
	mat k4_v(this->_rk4->dim(), n);
	mat k4_p(this->_rk4->dim(), n);

	// initial position and velocity for each celestial body
	mat orig_v(this->_rk4->dim(), n);
	mat orig_p(this->_rk4->dim(), n);

	#pragma region k1

	// calculate forces/accelerations based on current postions
	this->_rk4->setForces();

	for (int j = 0; j < n; j++) // for each celestial body
	{
		CelestialBody* cb = this->_rk4->body(j);

		if (!cb->fixed) // a fixed celestial body will never move
		{
			// store initial position and velocity
			orig_v.col(j) = *(cb->velocity);
			orig_p.col(j) = *(cb->position);

			// save values
			k1_v.col(j) = step * cb->acc();
			k1_p.col(j) = step * *(cb->velocity);

			// advance to mid-point after k1
			*(cb->velocity) = orig_v.col(j) + 0.5 * k1_v.col(j);
			*(cb->position) = orig_p.col(j) + 0.5 * k1_p.col(j);
		}
	}

	#pragma endregion

	#pragma region k2

	// calculate forces/accelerations based on current postions
	this->_rk4->setForces();

	for (int j = 0; j < n; j++) // for each celestial body
	{
		CelestialBody* cb = this->_rk4->body(j);

		if (!cb->fixed) // a fixed celestial body will never move
		{
			// save values
			k2_v.col(j) = step * cb->acc();
			k2_p.col(j) = step * *(cb->velocity);

			// switch to new mid-point using k2 instead
			*(cb->velocity) = orig_v.col(j) + 0.5 * k2_v.col(j);
			*(cb->position) = orig_p.col(j) + 0.5 * k2_p.col(j);
		}
	}

	#pragma endregion

	#pragma region k3

	// calculate forces/accelerations based on current postions
	this->_rk4->setForces();

	for (int j = 0; j < n; j++) // for each celestial body
	{
		CelestialBody* cb = this->_rk4->body(j);

		if (!cb->fixed) // a fixed celestial body will never move
		{
		// save values
			k3_v.col(j) = step * cb->acc();
			k3_p.col(j) = step * *(cb->velocity);

			// switch to end-point
			*(cb->velocity) = orig_v.col(j) + k3_v.col(j);
			*(cb->position) = orig_p.col(j) + k3_p.col(j);
		}
	}

	#pragma endregion

	#pragma region k4

	// calculate forces/accelerations based on current postions
	this->_rk4->setForces();

	for (int j = 0; j < n; j++) // for each celestial body
	{
		CelestialBody* cb = this->_rk4->body(j);

		if (!cb->fixed) // a fixed celestial body will never move
		{
			// save values
			k4_v.col(j) = step * cb->acc();
			k4_p.col(j) = step * *(cb->velocity);

			// finally, update position and velocity
			*(cb->velocity) = orig_v.col(j) + (1.0 / 6.0)*(k1_v.col(j) + 2.0 * k2_v.col(j) + 2.0 * k3_v.col(j) + k4_v.col(j));
			*(cb->position) = orig_p.col(j) + (1.0 / 6.0)*(k1_p.col(j) + 2.0 * k2_p.col(j) + 2.0 * k3_p.col(j) + k4_p.col(j));
		}
	}

	#pragma endregion
}

// a single Leapfrog step (parameter = step length)
void Solvers::Leapfrog(double step)
{
	int n = this->_leapfrog->n(); // number of celestial bodies

	double halfStep = (step / 2.0);

	// calculate forces/accelerations based on current postions
	this->_leapfrog->setForces();

	for (int j = 0; j < n; j++) // for each celestial body
	{
		CelestialBody* cb = this->_leapfrog->body(j);

		if (!cb->fixed) // a fixed celestial body will never move
		{
			// Leapfrog algorithm, step 1
			*(cb->velocity) += halfStep * cb->acc();

			// Leapfrog algorithm, step 2
			*(cb->position) += step * *(cb->velocity);
		}
	}

	// calculate the forces using the new positions
	this->_leapfrog->setForces();

	for (int j = 0; j < n; j++) // for each celestial body
	{
		CelestialBody* cb = this->_leapfrog->body(j);

		if (!cb->fixed) // a fixed celestial body will never move
		{
			// Leapfrog algorithm, step 3
			*(cb->velocity) += halfStep * cb->acc();
		}
	}
}

// a single Euler-Cromer step (parameter = step length)
void Solvers::Euler(double step)
{
	int n = this->_euler->n(); // number of celestial bodies

	// calculate forces/accelerations based on current postions
	this->_euler->setForces();

	for (int j = 0; j < n; j++) // for each celestial body
	{
		CelestialBody* cb = this->_euler->body(j);

		if (!cb->fixed) // a fixed celestial body will never move
		{
			// acc -> velocity (Euler-Cromer, for testing only)
			*(cb->velocity) += step * cb->acc();

			// velocity -> position (Euler-Cromer, for testing only)
			*(cb->position) += step * *(cb->velocity);
		}
	}
}

void Solvers::avTime(vec elapsedTime)
{
	double av = 0.0;
	for (int i = 0; i < elapsedTime.n_cols; i++)
	{
		av += elapsedTime(i);
	}
	av = av / elapsedTime.n_cols;
	printf("Average time spent for a step here : %lf \n", av);
}
