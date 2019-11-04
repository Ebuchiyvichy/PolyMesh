#pragma once
#include "functions.h"

mytype	*lagrange_val(mytype **mesh, int n, int test, int size, mytype **x)
{
	mytype *lagr = new mytype[size];
	mytype	C = 1;

	for (int i = 0; i < size; i++)
		lagr[i] = 0;
	for (int i = 0; i < size; i++)
	{
		for (int k = 0; k < n; k++)
		{
			C = 1;
			for (int j = 0; j < n; j++)
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

mytype	find_max(int TEST, mytype a, mytype b, int n)
{
	mytype	max = pow(my_func(a, TEST), n);
	mytype	**mesh;

	mesh = cheb_mesh(TEST, n, a, b);
	for (int i = 0; i < n; i++)
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