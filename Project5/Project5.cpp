// Project5.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include"armadillo"
#include<list>

using namespace arma;
using namespace std;

int _tmain(int argc, _TCHAR* argv[])
{
	mat A = randu<mat>(5, 5);
	cout << det(A) << endl << endl;

	mat L, U, P;
	if (lu(L, U, P, A))
	{
		cout << "LU success!" << endl;
		cout << L << endl << endl;
		cout << U << endl << endl;
	}

	vec v = randu<vec>(5);
	vec x = solve(L, v);
	cout << x << endl << endl;

	getchar(); // pause

	return 0;
}