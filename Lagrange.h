#pragma once
#include "functions.h"

mytype	*lagrange_val(mytype **mesh, int n, int test, int size, mytype **x)
{
	mytype *lagr = new mytype[size];
	mytype	C = 1;

	for (int i = 0; i <= size; i++)
		lagr[i] = 0;
	for (int i = 0; i <= size; i++)
	{
		for (int k = 0; k <= n; k++)
		{
			C = 1;
			for (int j = 0; j <= n; j++)
			{
				if (k != j)
					C *= ((x[i][0] - mesh[j][0]) / (mesh[k][0] - mesh[j][0]));
			}
			C *= mesh[k][1];
			lagr[i] += C;
		}
	}

	return lagr;
}

mytype	lagrange_one(mytype **mesh, int n, int test, mytype x)
{
	mytype	C = 1;
	mytype	tmp = 0;

	for (int k = 0; k <= n; k++)
	{
		C = 1;
		for (int j = 0; j <= n; j++)
		{
			if (k != j)
				C *= ((x - mesh[j][0]) / (mesh[k][0] - mesh[j][0]));
		}
		C *= mesh[k][1];
		tmp += C;
	}

	return tmp;
}

mytype	find_max(int TEST, mytype a, mytype b, int n)
{
	mytype	max = pow(my_func(a, TEST), n);
	mytype	**mesh;

	mesh = cheb_mesh(TEST, n, a, b);
	for (int i = 0; i <= n; i++)
	{
		if (max < pow(mesh[i][1], n))
			max = pow(mesh[i][1], n);
	}
	return max;
}

int		find_node(int TEST, mytype eps, mytype a, mytype b)
{
	int		n = 0;
	int		fact = 1;
	mytype	interval = b - a;
	int		power = 2;
	mytype	M = my_func(a, TEST);

	if (my_func(b, TEST) > M)
		M = my_func(b, TEST);
	while (eps < M / fact * interval / power)
	{
		M = find_max(TEST, a, b, n + 1);
		interval *= (b - a);
		fact *= n + 1;
		power *= 2;
		n++;
	}
	return n;
}

mytype	lagr_norm(mytype *lagr, int size, mytype **x)
{
	mytype	max;
	max = fabs(lagr[0] - x[0][1]);
	for (int i = 1; i <= size; i++)
		if (fabs(lagr[i] - x[i][1]) > max)
			max = fabs(lagr[i] - x[i][1]);
	return max;
}

mytype* progon(int DIM, mytype* a, mytype* b, mytype* c, mytype* d)

{
	mytype* alfa = new mytype[DIM];
	mytype* betta = new mytype[DIM];

	alfa[0] = -c[0] / b[0]; betta[0] = d[0] / b[0];

	for (int i = 1; i < DIM; i++) {
		alfa[i] = (-1)*c[i] / (b[i] + a[i] * alfa[i - 1]);
		betta[i] = (-a[i] * betta[i - 1] + d[i]) / (a[i] * alfa[i - 1] + b[i]);
	}

	for (int i = DIM - 2; i > -1; i--) {
		betta[i] += alfa[i] * betta[i + 1];
	}

	delete[] alfa;

	return  betta;
}

std::vector<double> Spline(int nbr, int N, mytype** val, mytype** x) {
	std::vector<double> a(nbr - 1);
	for (size_t i = 0; i < nbr - 1; ++i) {
		a[i] = val[i][1];
	}

	std::vector<double> c(nbr);
	std::vector<double> b(nbr - 1);
	std::vector<double> d(nbr - 1);

	std::vector<double> h(nbr - 1);
	std::vector<double> g(nbr - 1);
	std::vector<double> S(N);

	for (size_t i = 0; i < nbr - 1; ++i) {
		h[i] = val[i + 1][0] - val[i][0];
		g[i] = (val[i + 1][1] - val[i][1]) / h[i];
	}

	std::vector<double> ksi(nbr - 1);
	std::vector<double> eta(nbr - 1);

	ksi[nbr - 2] = -h[nbr - 3] / (2 * (h[nbr - 3] + h[nbr - 2]));
	eta[nbr - 2] = 3 * (g[nbr - 2] - g[nbr - 3]) / (2 * (h[nbr - 3] + h[nbr - 2]));

	for (size_t i = nbr - 3; i > 0; --i) {
		ksi[i] = -h[i - 1] / (2 * (h[i - 1] + h[i]) + h[i] * ksi[i + 1]);
		eta[i] = (3 * (g[i] - g[i - 1]) - h[i] * eta[i + 1]) / (2 * (h[i - 1] + h[i]) + h[i] * ksi[i + 1]);
	}

	c[0] = 0;

	for (size_t i = 1; i < nbr - 1; ++i) {
		c[i] = ksi[i] * c[i - 1] + eta[i];
	}

	c[nbr - 1] = 0;

	for (size_t i = 0; i < nbr - 1; ++i) {
		b[i] = g[i] - (c[i + 1] + 2 * c[i]) * h[i] / 3;
		d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
	}

	for (int j = 0; j != nbr - 1; j++) {
		for (int i = 0; i != N; i++) {
			if (x[i][0] >= val[j][0] && x[i][0] <= val[j + 1][0])
				S[i] = a[j] + b[j] * (x[i][0] - val[j][0]) + c[j] * pow((x[i][0] - val[j][0]), 2) + d[j] * pow((x[i][0] - val[j][0]), 3);
		}
	}
	return(S);
}

mytype	lagr_norm(std::vector<double> lagr, int size, mytype **x)
{
	mytype	max;
	max = fabs(lagr[0] - x[0][1]);
	for (int i = 0; i < size; i++) {
		if (fabs(lagr[i] - x[i][1]) > max) {
			max = fabs(lagr[i] - x[i][1]);
		}
	}
	return max;
}