#pragma once

// Downloaded from https://www.uio.no/studier/emner/matnat/fys/FYS3150/h13/gaussiandeviate.cpp and adapted into a class

static class GaussPDF
{
public:
	// ran2 for uniform deviates, initialize with negative seed.
	static double ran2(long *);

	// random numbers with gaussian distribution in units of standard deviations (avg = 0)
	static double gaussian_deviate(long *);
};

