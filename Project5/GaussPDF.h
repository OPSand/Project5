#pragma once

// Downloaded from https://www.uio.no/studier/emner/matnat/fys/FYS3150/h13/gaussiandeviate.cpp and adapted into a class

static class GaussPDF
{
public:
	// ran2 for uniform deviates, initialize with negative seed.
	static double ran2(long *);

	// function for gaussian random numbers
	static double gaussian_deviate(long *);
};

