#include "stdafx.h"
#include "Solvers.h"

// constructor (isim = number of sim in a series)
Solvers::Solvers(SolarSystem* system, bool useRK4, bool useLeapfrog, bool useEuler)
{
	// save pointer to original system
	this->_system = system;

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
		this->_rk4->name += "_rk4";
	}

	if (this->_useLeapfrog)
	{
		this->_leapfrog = new SolarSystem(*system); // deep copy
		this->_leapfrog->name += "_leapfrog";
	}

	if (this->_useEuler)
	{
		this->_euler = new SolarSystem(*system); // deep copy
		this->_euler->name += "_euler";
	}
}

// destructor
Solvers::~Solvers()
{
	if (this->_rk4 != nullptr)
	{
		delete this->_rk4;
	}

	if (this->_leapfrog != nullptr)
	{
		delete this->_leapfrog;
	}

	if (this->_euler != nullptr)
	{
		delete this->_euler;
	}
}

// call this to solve the equations and save results to files
// the system will determine the number of steps; all we need is to provide a step length
// returns the system in the state the Leapfrog algoritm leaves it in (or nullptr if that algorithm is not used)
vector<SolarSystem*>* Solvers::Solve(double step)
{
	clock_t start, finish; // timers
	int n = this->_system->nSteps();

	if (this->_useRK4)
	{
		start = clock(); // start timer

		for (int i = 0; i < n; i++) // for each time step
		{
			this->_rk4->plotCurrentStep(i, step); // if we want to plot this step, do it
			RK4(step); // perform step
		}

		// Timing part: End ...
		finish = clock();
		this->rk4Time = (double)(finish - start) / CLOCKS_PER_SEC; // To convert this into seconds !
	}

	if (this->_useLeapfrog)
	{
		start = clock(); // start timer

		for (int i = 0; i < n; i++) // for each time step
		{
			this->_leapfrog->plotCurrentStep(i, step); // if we want to plot this step, do its
			Leapfrog(step); // perform step
		}

		// Timing part: End ...
		finish = clock();
		this->leapfrogTime = (double)(finish - start) / CLOCKS_PER_SEC; // To convert this into seconds !
	}

	if (this->_useEuler)
	{
		start = clock(); // start timer

		for (int i = 0; i < n; i++) // for each time step
		{
			this->_euler->plotCurrentStep(i, step); // if we want to plot this step, do it
			Euler(step); // perform step
		}

		// Timing part: End ...
		finish = clock();
		this->eulerTime = (double)(finish - start) / CLOCKS_PER_SEC; // To convert this into seconds !
	}

	// put systems we want to return in here
	vector<SolarSystem*>* returnSystems = new vector<SolarSystem*>();

	if (this->_useLeapfrog)
	{
		// plot to file (independent of dimension)
		for (int i = 0; i < this->_leapfrog->dim(); i++)
		{
			ostringstream fname = ostringstream();
			fname << "sim_" << this->_leapfrog->name << "_pos" << i << ".dat";
			this->_leapfrog->plotDim(i, fname.str());
		}

		// make sure the system is up to date
		this->_leapfrog->calculate();

		// add to list of systems to return
		returnSystems->push_back(this->_leapfrog);
	}

	if (this->_useRK4)
	{
		// plot to file (independent of dimension)
		for (int i = 0; i < this->_rk4->dim(); i++)
		{
			ostringstream fname = ostringstream();
			fname << "sim_" << this->_rk4->name << "_pos" << i << ".dat";
			this->_rk4->plotDim(i, fname.str());
		}

		// make sure the system is up to date
		this->_rk4->calculate();

		// add to list of systems to return
		returnSystems->push_back(this->_rk4);
	}

	if (this->_useEuler)
	{
		// plot to file (independent of dimension)
		for (int i = 0; i < this->_euler->dim(); i++)
		{
			ostringstream fname = ostringstream();
			fname << "sim_" << this->_euler->name << "_pos" << i << ".dat";
			this->_euler->plotDim(i, fname.str());
		}

		// make sure the system is up to date
		this->_euler->calculate();

		// add to list of systems to return
		returnSystems->push_back(this->_euler);
	}

	cout << "DONE!" << endl << endl;

	this->totalTime = this->leapfrogTime + this->rk4Time + this->eulerTime;

	return returnSystems;
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
