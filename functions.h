#pragma once
#include <math.h>
#include <iomanip>
# define M_PI	3.14159265358979323846

typedef double mytype;
mytype	EPS = 10e-3;
int TEST = 4;
int	NBR = 3;
int	SIZE = 100;

mytype	my_func(mytype x, int test)
{
	if (test == 1)
		return pow(x, 2);
	else if (test == 2)
		return 1 / (1 + pow(x, 2));
	else if (test == 3)
		return 1 / atan(1 + 10 * pow(x, 2));
	else if (test == 4)
		return pow(4 * pow(x, 3)+ 2 * pow(x, 2) - 4 * x + 2, sqrt(2)) + asin(1 / (5 + x - pow(x, 2))) - 5;
}

void	vars(int test, mytype &a, mytype &b)
{
	if (test == 1 || test == 2 || test == 4)
	{
		a = -1; b = 1;
	}
	else if (test == 3)
	{
		a = -3; b = 3;
	}
}

void	print(int test, mytype a, mytype b)
{
	std::cout << "Test " << test << std::endl;
	if (test == 1)
		std::cout << "Your function is y = x^2" << std::endl;
	else if (test == 2)
		std::cout << "Your function is y = 1 / (1 + x^2)" << std::endl;
	else if (test == 3)
		std::cout << "Your function is y = 1 / atan(1 + 10x^2)" << std::endl;
	else if (test == 4)
		std::cout << "Your function is y = (4x^3 + 2x^2 - 4x + 2)^sqrt(2) + asin(1 / (5 + x - x^2)) - 5" << std::endl;
	std::cout << "a = " << a << '\t' << "b = " << b << std::endl;
}

mytype	**uniform_mesh(int test, int n, mytype a, mytype b)
{
	mytype	**mesh = new mytype* [n];

	for (int i = 0; i < n; i++)
		mesh[i] = new mytype[2];
	for (int i = 0; i < n; i++)
	{
		mesh[i][0] = a + i * (b - a) / n;
		mesh[i][1] = my_func(mesh[i][0], test);
	}
	return (mesh);
}

mytype	**cheb_mesh(int test, int n, mytype a, mytype b)
{
	mytype	**mesh = new mytype*[n];

	for (int i = 0; i < n; i++)
		mesh[i] = new mytype[2];
	for (int i = 0; i < n; i++)
	{
		mesh[i][0] = (a + b) / 2 + (b - a) / 2 * cos((2 * i + 1) * M_PI / (2 * (n + 1)));
		mesh[i][1] = my_func(mesh[i][0], test);
	}
	return (mesh);
}

void	print_mesh(mytype **mesh, int n)
{
	std::cout << "\t xi \t\t | \t y(xi) " << std::endl;
	for (int i = 0; i < n; i++)
	{
		std::cout.precision(6);
		std::cout << std::setw(8) << std::fixed << '\t' << mesh[i][0] << "\t | \t" << mesh[i][1] << std::endl;
	}
}

mytype *values(int n, int test, mytype** val)
{
    for (int i = 0; i != n; i++)
    {
        val[i][1] = my_func(val[i][0],test);
    }
    return(val[1]);
}

