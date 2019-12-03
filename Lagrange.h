#pragma once
#include "functions.h"

std::vector<mytype> lagrange_val(mytype **mesh, int n, int test, int size, mytype **x)
{
	std::vector<mytype> lagr(size);
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

mytype* progon(int DIM, mytype* a, mytype* b, mytype* c, mytype* d)

{
    mytype* alfa = new mytype[DIM];
    mytype* betta = new mytype[DIM];

    alfa[0] = -c[0] / b[0]; betta[0] = d[0] / b[0];

    for (int i = 1; i < DIM; i++) {
        alfa[i] = -c[i] / (b[i] + a[i] * alfa[i - 1]);
        betta[i] = (-a[i] * betta[i - 1] + d[i]) / (a[i] * alfa[i - 1] + b[i]);
    }

    for (int i = DIM - 2; i > -1; i--) {
        betta[i] += alfa[i] * betta[i + 1];
    }

    delete[] alfa;

    return  betta;
}

mytype *splain(mytype **mesh, int n, int test, int size, mytype **x)
{
    mytype *a = new mytype[n];
    mytype *b = new mytype[n];
    mytype *c = new mytype[n];
    mytype *d = new mytype[n];

    mytype* h = new mytype[n];
    mytype* g = new mytype[n];

    for (int i = 0; i != n; i++) {
        a[i] = mesh[i][1];
        h[i] = mesh[i+1][0] - mesh[i][0];
        g[i] = (mesh[i+1][1]- mesh[i][1])/h[i];
    }
    mytype *d1 = new mytype[n-1];
    mytype *d2 = new mytype[n-1];
    mytype *d3 = new mytype[n-1];
    mytype *rvalue = new mytype[n-1];
    for (int  i = 0; i != n -1; i++)
    {
        d2[i] = 2 * (h[i+1] + h[i]);
        d1[i] = h[i];
        d3[i] = h[i + 1];
        rvalue[i] = 3 * (g[i+1] - g[i]);
    }
    d1[0] = 0; d3[n-2] = 0;
    rvalue = progon(n-1, d1, d2, d3, rvalue);
    c[0] = 0;
    for (int i = 1; i < n; i++)
        c[i] = rvalue[i - 1];

    for (int i = 0; i < n - 1; i++) {
        b[i] = g[i] - (c[i + 1] + 2 * c[i])* h[i] / 3;
        d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
    }

    b[n - 1] = g[n - 1] - 2 * c[n - 1] * h[n - 1] / 3;
    d[n - 1] = -c[n - 1] / (3 * h[n - 1]);


    mytype* S = new mytype[size + 1];

    for (int j = 0, i = 0; j < n; j++)
        for (; (i < size + 1) && (x[i][0] <= mesh[j + 1][0]); i++)
            S[i] = a[j] + b[j] * (x[i][0] - mesh[j][0]) + c[j] * pow((x[i][0] - mesh[j][0]), 2) + d[j] * pow((x[i][0] - mesh[j][0]), 3);
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;

    delete[] h;
    delete[] g;

    delete[] d1;
    delete[] d2;
    delete[] d3;
    delete[] rvalue;

    return S;
}
mytype errors(int n, mytype *in, mytype *tr) {
    mytype max = 0;

    for (int i = 0; i != n; i++)
        if (fabs(tr[i] - in[i]) > max)
            max = tr[i] - in[i];

        return max;
}