#include "stdafx.h"
#include "Solvers.h"


Solvers::Solvers(const SolarSystem* system, bool useRK4, bool useLeapfrog, bool useEuler)
{
	// save pointer to original system
	this->_system = system;

	// set flags
	this->_useRK4 = useRK4;
	this->_useLeapfrog = useLeapfrog;
	this->_useEuler = useEuler;

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


Solvers::~Solvers()
{
	delete this->_rk4;
	delete this->_leapfrog;
	delete this->_euler;
}


void Solvers::Solve(int step, int nSteps, int plotEvery)
{
	// number of steps to actually plot to file
	int nPlot = (nSteps / plotEvery);

	if (this->_useRK4)
	{
		RK4(step, nSteps, plotEvery, nPlot);
	}

	if (this->_useLeapfrog)
	{
		Leapfrog(step, nSteps, plotEvery, nPlot);
	}

	if (this->_useEuler)
	{
		Euler(step, nSteps, plotEvery, nPlot);
	}
}


void Solvers::RK4(int step, int nSteps, int plotEvery, int nPlot)
{
	cout << "Runge-Kutta..." << endl << endl;

	int n = this->_rk4->n();

	// iterate and plot coordinates
	for (int i = 0; i < nSteps; i++) // for each time step
	{
		this->_rk4->plotCurrentPositions(i % plotEvery == 0); // if we want to plot this step, do it

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

	this->_rk4->plotDim(0, "X_rk4.dat");
	this->_rk4->plotDim(1, "Y_rk4.dat");
	if (this->_rk4->dim() >= 3)
	{
		this->_rk4->plotDim(2, "Z_rk4.dat");
	}

	cout << "Finished plotting " << nPlot << " of " << nSteps << " steps (Runge-Kutta)!" << endl << endl;
}


void Solvers::Leapfrog(int step, int nSteps, int plotEvery, int nPlot)
{
	cout << "Leapfrog..." << endl << endl;

	int n = this->_leapfrog->n(); // number of celestial bodies

	// iterate and plot coordinates
	for (int i = 0; i < nSteps; i++) // for each time step
	{
		this->_leapfrog->plotCurrentPositions(i % plotEvery == 0); // if we want to plot this step, do it

		// calculate forces/accelerations based on current postions
		this->_leapfrog->setForces();

		for (int j = 0; j < n; j++) // for each celestial body
		{
			CelestialBody* cb = this->_leapfrog->body(j);

			if (!cb->fixed) // a fixed celestial body will never move
			{
				// TODO: implement Leapfrog algorithm
			}
		}
	}

	this->_leapfrog->plotDim(0, "X_leapfrog.dat");
	this->_leapfrog->plotDim(1, "Y_leapfrog.dat");
	if (this->_leapfrog->dim() >= 3)
	{
		this->_leapfrog->plotDim(2, "Z_leapfrog.dat");
	}

	cout << "Finished plotting " << nPlot << " of " << nSteps << " steps (Leapfrog)!" << endl << endl;
}


void Solvers::Euler(int step, int nSteps, int plotEvery, int nPlot)
{
	cout << "Euler-Cromer..." << endl << endl;

	int n = this->_euler->n(); // number of celestial bodies

	// iterate and plot coordinates
	for (int i = 0; i < nSteps; i++) // for each time step
	{
		this->_euler->plotCurrentPositions(i % plotEvery == 0); // if we want to plot this step, do it

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

	this->_euler->plotDim(0, "X_euler.dat");
	this->_euler->plotDim(1, "Y_euler.dat");
	if (this->_euler->dim() >= 3)
	{
		this->_euler->plotDim(2, "Z_euler.dat");
	}

	cout << "Finished plotting " << nPlot << " of " << nSteps << " steps (Euler-Cromer)!" << endl << endl;
}