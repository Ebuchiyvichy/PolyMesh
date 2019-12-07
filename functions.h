#pragma once
#include <math.h>
#include <iomanip>
# define M_PI	3.14159265358979323846

typedef double mytype;
mytype	EPS = 10e-15;
int TEST = 1;
int	NBR = 100;//количество разбиений
int	SIZE = 500;

mytype	my_func(mytype x, int test)
{
	if (test == 1)
		return pow(x, 2);
	else if (test == 2)
		return 1 / (1 + pow(x, 2));
	else if (test == 3)
		return 1 / atan(1 + 10 * pow(x, 2));
	else if (test == 4)
		return pow(4 * pow(x, 3) + 2 * pow(x, 2) - 4 * x + 2, sqrt(2)) + asin(1 / (5 + x - pow(x, 2))) - 5;
	else if (test == 5)
		return (exp(x));
	else if (test == 6)
		return (1);
	else if (test == 7 || test == 8)
		return sin(M_PI * x);
	else if (test == 9 || test == 10)
	{
		if (x < 0)
			return M_PI * x;
		else
			return my_func(x, 7);
	}
	else if (test == 11) // Runge
		return 1 / (1 + 25 * x * x);
}

void	vars(int test, mytype &a, mytype &b)
{
	if (test == 1 || test == 2 || test == 4 || test == 6 || test == 7 || test == 9 || test == 11)
	{
		a = -1; b = 1;
	}
	else if (test == 3)
	{
		a = -3; b = 3;
	}
	else if (test == 5)
	{
		a = 0; b = 2;
	}
	else if (test == 8 || test == 10)
	{
		a = -1.25; b = 1.25;
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

	for (int i = 0; i <= n; i++)
		mesh[i] = new mytype[2];
	for (int i = 0; i <= n; i++)
	{
		mesh[i][0] = a + i * (b - a) / n;
		mesh[i][1] = my_func(mesh[i][0], test);
	}
	return (mesh);
}

mytype	**cheb_mesh(int test, int n, mytype a, mytype b)
{
	mytype	**mesh = new mytype*[n];

	for (int i = 0; i <= n; i++)
		mesh[i] = new mytype[2];
	for (int i = 0; i <= n; i++)
	{
		mesh[i][0] = (a + b) / 2 + (b - a) / 2 * cos((2 * i + 1) * M_PI / (2 * (n + 1)));
		mesh[i][1] = my_func(mesh[i][0], test);
	}
	return (mesh);
}

void	print_mesh(mytype **mesh, int n)
{
	std::cout << "\t xi \t\t | \t y(xi) " << std::endl;
	for (int i = 0; i <= n; i++)
	{
		std::cout.precision(6);
		std::cout << std::setw(8) << std::fixed << '\t' << mesh[i][0] << "\t | \t" << mesh[i][1] << std::endl;
	}
}