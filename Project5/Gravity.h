#pragma once

#include "CelestialBody.h"

class CelestialBody; // forward declaration to avoid circular reference

static class Gravity
{
protected:
	double _G;
	double _epsilon;
public:
	Gravity(double G, double epsilon);
	Gravity(const Gravity& other);
	~Gravity();
	Gravity operator= (const Gravity& other);
	vec force(CelestialBody* cb_i, CelestialBody* cb_j);
};

