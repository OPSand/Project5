// Project5.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "armadillo"
#include <random>
#include "SolarSystem.h"
#include "Solvers.h"
#include "CelestialBodyInitializer.h"

using namespace arma;
using namespace std;

#pragma region Local Methods
// fit radial distribution to curve (using least squares)
// returns a vector: [ n0, r0, n0/N^2r0, /N^(-1/3) ]
vec radialDistFitLSq(mat radialDist, int N, int nNR)
{
	// these factors are educated guesses as to where we will find the best match
	const double nFact = 2.0;
	const double rFact = 1.0;

	// dependence of n0 and r0 on N taken from http://arxiv.org/abs/1011.0614
	double nScale = pow((double)N, 2.0);
	double rScale = pow((double)N, -(1.0 / 3.0));

	// range to search
	double minN0 = 0.0;
	double minR0 = 0.0;
	double maxN0 = nFact * nScale;
	double maxR0 = rFact * rScale;

	// step size for n0 and r0 trials, respectively
	double deltaN = (maxN0 - minN0) / (nNR - 1.0);
	double deltaR = (maxR0 - minR0) / (nNR - 1.0);

	// return this:
	vec ret = vec(4);

	// determine number of bins in histogram
	int boxes = radialDist.n_rows;

	double sum; // sum of deviations squared
	double minsum; // minimal sum (found so far)

	// loop over n0 values
	for (int n = 0; n < nNR; n++)
	{
		double n0 = minN0 + n * deltaN;

		// loop over r0 values
		for (int r = 0; r < nNR; r++)
		{
			double r0 = minR0 + r * deltaR;

			sum = 0.0; // reset

			// for every bin, compare value to formula
			// add square of deviation to sum
			for (int i = 0; i < boxes; i++)
			{
				double formula = n0 / (1.0 + pow((radialDist(i, 0) / r0), 4.0));
				sum += pow(formula - radialDist(i, 2), 2.0);
			}

			if ((n == 0) && (r == 0)) // first attempt: nothing to compare to
			{
				minsum = sum;
				ret(0) = n0;
				ret(1) = r0;
			}
			else
			{
				if (sum < minsum) // better match than any previous
				{
					minsum = sum;
					ret(0) = n0;
					ret(1) = r0;
				}
			}
		}
	}

	ret(2) = ret(0) / nScale; // n0/N^2
	ret(3) = ret(1) / rScale; // r0/N^(-1/3)

	return ret;
}

// (for benchmarking): convert 2D polar coordinates to cartesian
vec toCartesian2D(double r, double theta)
{
	vec p = vec(2);
	p[0] = r*cos(theta); // x
	p[1] = r*sin(theta); // y
	return p;
}

// (for benchmarking): returns an angle in radians that is orthogonal to theta
double orthogonal2D(double theta, bool clockwise = false)
{
	if (clockwise)
	{
		return (theta - 0.5 * cPI);
	}
	else // counterclockwise
	{
		return (theta + 0.5 * cPI);
	}
}

// (for benchmarking): initialize 2-body problem (in 2D)
void initial2D(CelestialBody* cb, double d, double v, double theta = -1.0)
{
	if (theta == -1.0) // no angle given
	{
		// create random angle
		long* seed = new long(-time(NULL));
		theta = (2 * cPI), GaussPDF::ran2(seed);
		delete seed; // farming in a star cluster is risky business, anyway
	}

	*(cb->position) = toCartesian2D(d, theta);
	*(cb->velocity) = toCartesian2D(v, orthogonal2D(theta)); // orthogonal (counterclockwise) vs. position
}

// (for benchmarking): returns an angle in radians that is orthogonal to theta
double orthogonal3D(double theta, bool clockwise = false)
{
	return (theta + 0.5 * cPI);
}

vec toCartesian3D(double r, double theta, double psi)
{
	vec p = vec(3);
	p[0] = r*sin(theta)*cos(psi); // x
	p[1] = r*sin(theta)*sin(psi); // y
	p[2] = r*cos(theta);

	return p;
}

void initial3D(CelestialBody* cb, double d, double v)
{
	// determine position(r0) - independent of system dimension
	int dim = 3;
	long* seedTheta = new long(-time(NULL));
	long* seedPsi = new long(-time(NULL)); // long time, no C?

	double theta;
	double psi;

	vec x = vec(dim);
	theta = (2 * cPI), GaussPDF::ran2(seedTheta);
	psi = (2 * cPI), GaussPDF::ran2(seedPsi);

	delete seedTheta;
	delete seedPsi;

	*(cb->position) = toCartesian3D(d, theta, psi); // scale to correct radius
	*(cb->velocity) = toCartesian3D(v, orthogonal3D(theta), psi); // orthogonal (counterclockwise) vs. position
}
#pragma endregion

// Entry point for console application
int _tmain(int argc, _TCHAR* argv[])
{
	#pragma region Flags and settings

	// dimensions
	const int DIM = 3;

	// initialization & time steps (common)
	const int N = 100; // number of celestial bodies
	const double R0 = 20.0; // initial radius in ly
	const double TOTAL_M = 1000.0; // solar masses (to be kept constant if we loop over N)
	const double STD_FACTOR = 0.1; // % factor of average
	const double CRUNCH_TIMES = 4.0; // # of crunch times to simulate for
	const double EPSILON = 0.0; // first epsilon value to try (ignored if not looping over epsilon values)
	const int N_STEPS = 1000; // number of steps total
	const int N_PLOT = 100; // number of steps to plot (must be <= N_STEPS)
	const int N_NR = 1000; // number of n0 and r0 values to try (curve fitting)
	const double CURVEFIT_STDDEV = 1.0; // max r limit for curve fitting in standard deviations
	const double AVG_BIN = 20.0; // avg. number of particles in each bin (curve fitting)

	// initialization & time steps (run many with different n/epsilon, same total mass)
	const int N_SIMS = 16; // number of simulations to run (set to 1 to run just once)
	const int N_END = 2500; // max N for last sim (ignored if N_SIMS == 1 or if EPSILON_LOOP == true)
	const double EPSILON_END = 0.15; // max epsilon for last sim (ignored if N_SIMS == 1 or if EPSILON_LOOP == false)

	// constants calculated from other constants
	const double AVG_M = (TOTAL_M / (double)N); // average mass, solar masses
	const double STEP = CRUNCH_TIMES / ((double)N_STEPS - 1.0); // step size (in crunch times)
	const int PLOT_EVERY = N_STEPS / N_PLOT; // plot every ...th step

	// physical constants
	const double LY = 9.4607e15; // m
	const double MYR = cYr * 1.0e6; // s
	const double M_SUN = 1.9891e30; // kg
	const double G_YLS = cG * M_SUN * pow(cYr, 2.0) / pow(LY, 3.0); // G in years, ly, solar masses

	// flags
	const bool EPSILON_LOOP = true; // vary epsilon instead of n
	const bool USE_LEAPFROG = true; // use Leapfrog method
	const bool USE_RK4 = false; // use Runge-Kutta method
	const bool USE_EULER = false; // use Euler-Cromer method
	const bool DEBUG = false; // for debugging only
	const bool BENCHMARK = false; // To test against the project 3 code
	#pragma endregion

	if (!BENCHMARK)
	{
		#pragma region Initialize Series

		int deltaN = 0;
		double deltaEpsilon = 0.0;

		int nSims = N_SIMS; // avoid const error

		if (!EPSILON_LOOP) // loop over n
		{
			if (nSims > 1) // avoid division by 0
			{
				deltaN = (N_END - N) / (nSims - 1);
			}
		}
		else // loop over epsilon
		{
			if (nSims > 1) // avoid division by 0
			{
				deltaEpsilon = (EPSILON_END - EPSILON) / (nSims - 1);
			}
		}

		// store info for matlab plots for the entire series of simulations
		mat* leapfrogPlot = new mat(N_SIMS, 18);
		mat* rk4Plot = new mat(N_SIMS, 18);
		mat* eulerPlot = new mat(N_SIMS, 18);

		#pragma endregion

		// run one or more simulations
		for (int isim = 0; isim < N_SIMS; isim++)
		{
			#pragma region Initialize Simulation

			ostringstream fname; // for dynamic file names

			// system name
			ostringstream isimstr = ostringstream();
			isimstr << isim;

			int nParticles;
			double eps;

			// create gravity (we will update it later)
			Gravity g = Gravity(0.0, 0.0);

			if (!EPSILON_LOOP) // loop over n
			{
				// determine number of particles
				nParticles = N + isim * deltaN;

				// epsilon will be set by the initializer
			}
			else // loop over epsilon
			{
				// set number of particles
				nParticles = N;

				// calculate epsilon value to use
				eps = EPSILON + isim * deltaEpsilon;
			}

			// create system
			SolarSystem* system = new SolarSystem(DIM, N_STEPS, PLOT_EVERY, &g);
			system->name = isimstr.str(); // set name

			// conserve the total mass & scale std.dev to average
			double avgMass = TOTAL_M / N;
			double stdMass = STD_FACTOR * avgMass;

			// add N randomly initialized celestial bodies (also sets G and epsilon)
			CelestialBodyInitializer::initialize(system, nParticles, avgMass, stdMass, R0, STEP);

			// overwrite epsilon if we loop over epsilon values
			if (EPSILON_LOOP)
			{
				g.setEpsilon(eps);
			}
			else // if not, might as well save it for later
			{
				eps = g.epsilon();
			}

			#pragma endregion

			#pragma region Solve and Plot Simulation

			cout << endl << "--- SIMULATION " << (isim + 1) << " OF " << N_SIMS << " ---" << endl;
			cout << "N = " << nParticles << endl;
			cout << "epsilon = " << g.epsilon() << endl << endl;

			// call this only when initialization is 100% complete!
			Solvers solv = Solvers(system,  USE_RK4, USE_LEAPFROG, USE_EULER);

			if (DEBUG)
			{
				for (int i = 0; i < system->n(); i++)
				{
					CelestialBody* cb = system->body(i);
					cout << "mass = " << cb->mass << endl;
					cout << "position = " << *(cb->position) << endl << endl;
				}

				cout << "T_CRUNCH = " << CelestialBodyInitializer::tCrunch(R0, system, G_YLS) << endl;
				cout << "G = " << g.G() << endl;
				cout << "G_YLS = " << G_YLS << endl;
				cout << endl;
			}

			// total energy before simulation starts
			double EtotBefore = system->EpTotal(false) + system->EkTotal(false);

			// Plot energies before simulation (all should be bound)
			cout << "E_k before: " << system->EkTotal(false) << endl;
			cout << "E_p before: " << system->EpTotal(false) << endl;
			cout << "E_tot before: " << EtotBefore << endl;

			if (DEBUG)
			{
				cout << "E_k before (bound): " << system->EkTotal(true) << endl;
				cout << "E_p before (bound): " << system->EpTotal(true) << endl;
				cout << "E_tot before (bound): " << (system->EpTotal(true) + system->EkTotal(true)) << endl;

				cout << "Center of mass (all): " << system->centerOfMass(false) << endl;
				cout << "Center of mass (bound): " << system->centerOfMass(true) << endl;

				double maxR = R0;
				mat radial = system->radialDistribution(maxR, AVG_BIN, true);
				fname = ostringstream();
				fname << "radial_before_" << system->name << ".dat";
				radial.save(fname.str(), raw_ascii);

				cout << radial << endl << endl;

				cout << "n(r0) = " << system->numDens(R0) << endl; // nice to compare this to n(r) in the bins

				cout << "(bound) avg = " << system->avgDistCoM(true) << endl;
				cout << "(bound) stdDev = " << system->stdDevDistCoM(true) << endl;
			}

			cout << endl << "Running simulation";

			// this is where the magic happens :)
			vector<SolarSystem*>* systems = solv.Solve(STEP);

			// for each algorithm used
			for (int i = 0; i < systems->size(); i++)
			{
				// get system reference
				system = systems->at(i);

				cout << " ** " << system->name << " **" << endl << endl;

				// Plot time used and...
				double elapsedTime;
				// ...locate correct matrix for this algorithm
				mat* isimsPlot;

				// determine which algorithm was used from system name
				if (system->name.find("rk4") != string::npos)
				{
					elapsedTime = solv.rk4Time;
					isimsPlot = rk4Plot;
				}
				else if (system->name.find("euler") != string::npos)
				{
					elapsedTime = solv.eulerTime;
					isimsPlot = eulerPlot;
				}
				else // leapfrog
				{
					elapsedTime = solv.leapfrogTime;
					isimsPlot = leapfrogPlot;
				}

				cout << "Elapsed time:" << elapsedTime << endl;

				// number of particles
				double n = system->n();

				// record final percentage of bound particles
				double boundParticles = system->nBound();
				double finalBound = (boundParticles / n);

				// total/bound energy after simulation
				double Etot = system->EpTotal(false) + system->EkTotal(false);
				double EtotBound = system->EpTotal(true) + system->EkTotal(true);

				// energy conservation (relative deltas: close to 0 is good)
				double deltaErel = (Etot - EtotBefore) / abs(EtotBefore);
				double deltaErelBound = ((EtotBound/finalBound) - EtotBefore) / abs(EtotBefore);
				// NOTE: We have to scale EtotBound up to account for the lower number of bound particles remaining

				// total kinetic/potential energy for bound particles
				double EkBound = system->EkTotal(true);
				double EpBound = system->EpTotal(true);

				// total energy lost to particle ejection
				double lostEnergy = (Etot - EtotBound);

				// distance to bound center of mass
				double avgComBound = system->avgDistCoM(true);
				double stdComBound = system->stdDevDistCoM(true);

				cout << "E_k after (bound): " << EkBound << endl;
				cout << "E_p after (bound): " << EpBound << endl;
				cout << "E_tot after: " << Etot << endl;
				cout << "E_tot after (bound): " << EtotBound << endl;

				if (DEBUG)
				{
					cout << "E_k after: " << system->EkTotal(false) << endl;
					cout << "E_p after: " << system->EpTotal(false) << endl;

					cout << "Center of mass (all): " << system->centerOfMass(false) << endl;
					cout << "Center of mass (bound): " << system->centerOfMass(true) << endl;
				}

				// get radial distribution & save to file
				// set max radius to average + x standard deviations (where we choose x)
				double maxR = system->avgDistCoM(true) + CURVEFIT_STDDEV * system->stdDevDistCoM(true);
				mat radial = system->radialDistribution(maxR, AVG_BIN, true);
				fname = ostringstream();
				fname << "radial_after_" << system->name << ".dat";
				radial.save(fname.str(), raw_ascii);

				if (DEBUG)
				{
					cout << radial << endl << endl;
				}

				// save the number of bound particles per time step for plotting
				mat nbp = system->nBoundPlot();
				fname = ostringstream();
				fname << "nbound_" << system->name << ".dat";
				nbp.save(fname.str(), raw_ascii);

				// perform curve fitting
				vec testFit = radialDistFitLSq(radial, system->n(), N_NR);

				if (DEBUG)
				{
					cout << "n0 = " << testFit(0) << endl;
					cout << "r0 = " << testFit(1) << endl;
					cout << "n0 / N^2 = " << testFit(2) << endl << endl;
					cout << "r0 / N^(-1/3) = " << testFit(3) << endl;

					cout << "(bound) avg = " << avgComBound << endl;
					cout << "(bound) stdDev = " << stdComBound << endl;
					cout << "(all) avg = " << system->avgDistCoM(false) << endl;
					cout << "(all) stdDev = " << system->stdDevDistCoM(false) << endl;
				}

				// save info for the series of simulations:
				// parameters
				isimsPlot->at(isim, 0) = nParticles;
				isimsPlot->at(isim, 1) = g.epsilon();
				// time
				isimsPlot->at(isim, 2) = elapsedTime;
				// energy conservation
				isimsPlot->at(isim, 3) = EtotBefore;
				isimsPlot->at(isim, 4) = Etot;
				isimsPlot->at(isim, 5) = EtotBound;
				// relative change in energy
				isimsPlot->at(isim, 6) = deltaErel;
				isimsPlot->at(isim, 7) = deltaErelBound;
				// virial theorem data
				isimsPlot->at(isim, 8) = EkBound;
				isimsPlot->at(isim, 9) = EpBound;
				// NOTE: # 10 is below! (since we reset epsilon, we save that one for last)
				// curve fitting
				isimsPlot->at(isim, 11) = testFit(0); // n0
				isimsPlot->at(isim, 12) = testFit(1); // r0
				isimsPlot->at(isim, 13) = testFit(2); // n0 / N^2
				isimsPlot->at(isim, 14) = testFit(3); // r0 / N^(-1/3)
				// distance to bound center of mass
				isimsPlot->at(isim, 15) = avgComBound;
				isimsPlot->at(isim, 16) = stdComBound;
				// final number of bound particles
				isimsPlot->at(isim, 17) = finalBound;

				// Finally:
				// calculate "classic" potential energy instead for bound particles
				// (i.e. ignoring epsilon in the potential)
				// NOTE: After this point in the loop, use eps and not g.epsilon()!
				g.setEpsilon(0.0);
				system->calculate(); // re-compute potential energies
				double EpBoundClassic = system->EpTotal(true);
				isimsPlot->at(isim, 10) = EpBoundClassic;

				if (DEBUG)
				{
					cout << "isimsPlot:" << *isimsPlot << endl;
					getchar();
					if (leapfrogPlot != nullptr)
					{
						cout << "leapfrogPlot:" << *leapfrogPlot << endl;
					}
					getchar();
				}
			}

			// free up resources
			while (systems->size() > 0)
			{
				system = systems->at(systems->size() - 1); // get last system
				systems->pop_back(); // pop it off the stack
				// NOTE: Do not delete system! That is handled when solv expires!
			}
			delete systems; // finally, delete the vector

			#pragma endregion
		}

		#pragma region Plot Series

		// save data for the series of simulations
		if (USE_LEAPFROG)
		{
			if (DEBUG)
			{
				cout << *leapfrogPlot;
			}
			leapfrogPlot->save("plotLeapfrog.dat", raw_ascii);
		}
		if (USE_RK4)
		{
			rk4Plot->save("plotRK4.dat", raw_ascii);
		}
		if (USE_EULER)
		{
			eulerPlot->save("plotEuler.dat", raw_ascii);
		}

		// free up resources
		delete leapfrogPlot;
		delete rk4Plot;
		delete eulerPlot;

		#pragma endregion
	}

	#pragma region Benchmark
	if (BENCHMARK) // re-create Project 3 for comparison - can we reproduce the results?
	{
		//getchar(); // pause

		// Time steps
		const int STEP_BM = 24 * 60 * 60; // step length (s)
		const int N_STEPS_BM = 300 * 365; // number of steps total
		const int PLOT_EVERY_BM = 1; // plot every ...th step
		// create gravity
		Gravity g_Bench = Gravity(cG, 0);
		printf("Entering the Benchmark part \n");
		// create system
		SolarSystem* system_BM = new SolarSystem(3, N_STEPS_BM, PLOT_EVERY_BM, &g_Bench);
		const double M_SUN_2 = 2e30;
		const double M_EARTH = 6e24;
		const double D_SUN = 0.0 * cAU;
		const double D_EARTH = 1.0 * cAU;
		const double V_SUN = 0.0;
		const double V_EARTH = 29.8e3;
		//system->grav()->setG(cG);
		CelestialBody* sun = new CelestialBody("Sun", M_SUN_2, system_BM, true);
		CelestialBody* earth = new CelestialBody("Earth", M_EARTH, system_BM, false);
		initial3D(earth, D_EARTH, V_EARTH);

		double E_k_init = system_BM->EkTotal(false);
		double E_p_init = system_BM->EpTotal(false);
		double E_tot_init = system_BM->EpTotal(false) + system_BM->EkTotal(false);
		
		system_BM->name = "BM";

		// call this only when initialization is 100% complete!
		Solvers solv = Solvers(system_BM, USE_RK4, USE_LEAPFROG, USE_EULER);

		// this is where the magic happens :)
		vector<SolarSystem*>* systemsBM = solv.Solve(STEP_BM);
		printf("Aaand ... We're done \n");

		cout << "E_k before : " << E_k_init << endl;
		cout << "E_p before : " << E_p_init << endl;
		cout << "E_tot before : " << E_tot_init << endl;

		cout << "E_tot after: " << (system_BM->EpTotal(false) + system_BM->EkTotal(false)) << endl;		
		cout << "E_k after: " << system_BM->EkTotal(false) << endl;
		cout << "E_p after: " << system_BM->EpTotal(false) << endl;
	}

	cout << endl << "Press ENTER to exit...";
	getchar(); // pause

	#pragma endregion 

	return 0;
}