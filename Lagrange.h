
#pragma once
#include    "functions.h"
#include    <vector>

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

mytype	lagr_norm(std::vector<double> lagr, int size, mytype **x)
{
    mytype	max;
    max = fabs(lagr[0] - x[0][1]);
    for (int i = 0; i <= size; i++) {
        if (fabs(lagr[i] - x[i][1]) > max) {
            max = fabs(lagr[i] - x[i][1]);
        }
    }
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

//mytype *splain(mytype **mesh, int n, int test, int size, mytype **x)
//{
//    mytype *a = new mytype[n-1];
//    mytype *b = new mytype[n-1];
//    mytype *c = new mytype[n];
//    mytype *d = new mytype[n-1];
//
//    mytype* h = new mytype[n];
//    mytype* g = new mytype[n];
//
//    for (int i = 0; i != n-1; i++) {
//        a[i] = mesh[i][1];
//        h[i] = mesh[i+1][0] - mesh[i][0];
//        g[i] = (mesh[i+1][1]- mesh[i][1])/h[i];
//    }
//    mytype *d1 = new mytype[n-1];
//    mytype *d2 = new mytype[n-1];
//    mytype *d3 = new mytype[n-1];
//    mytype *rvalue = new mytype[n-1];
//    for (int  i = 0; i != n; i++)
//    {
//        d2[i] = 2 * (h[i+1] + h[i]);
//        d1[i] = h[i];
//        d3[i] = h[i + 1];
//        rvalue[i] = 3 * (g[i+1] - g[i]);
//    }
//    d1[0] = 0; d3[n-2] = 0;
//    rvalue = progon(n, d1, d2, d3, rvalue);
//    c[0] = 0;
//    for (int i = 1; i < n; i++)
//        c[i] = rvalue[i - 1];
//    mytype *betta= new mytype[n-1];
//    mytype *alfa= new mytype[n-1];
//
//    betta[n-2] = -h[n-3]/(2*h[n-3]+h[n-2]);
//    alfa[n-2] = 3*(g[n-2]-g[n-3])/(2*(h[n-3]+h[n-2]));
//    for (int i = n-3; i > 0; i--)
//    {
//        betta[i] = -h[i-1]/(2*(h[i-1]+h[i]) +h[i]*betta[i+1]);
//        alfa[i] = 3*((g[i]-g[i-1])-h[i] * alfa[i+1])/(2*(h[i-1]+h[i])+h[i]*betta[i]);
//    }
//    c[0] = 0; c[n-1] = 0;
//    for (int i = 1; i < n-1; ++i)
//        c[i]= betta[i]*c[i-1]+alfa[i];
//
//    for (int i = 0; i < n - 1; ++i) {
//        b[i] = g[i] - (c[i + 1] + 2 * c[i])* h[i] / 3.;
//        d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
//    }
//
//
//    mytype* S = new mytype[size + 1];
//
//    for (int j = 0; j != n; j++) {
//        for (int i = 0; i != SIZE; i++) {
//            if (x[i][0] >= mesh[j][0] && x[i][0] <= mesh[j+1][0])
//             S[i] = a[j] + b[j] * (x[i][0] - mesh[j][0]) + c[j] * pow((x[i][0] - mesh[j][0]), 2) + d[j] * pow((x[i][0] - mesh[j][0]), 3);
//        }
//    }
//    delete[] a;
//    delete[] b;
//    delete[] c;
//    delete[] d;
//
//    delete[] h;
//    delete[] g;
//
//    delete[] d1;
//    delete[] d2;
//    delete[] d3;
//    delete[] rvalue;
//
//    return S;
//}
//mytype* Spline(int n, int N, mytype** val, mytype* x)
//{
//    mytype* a = new mytype[n];
//    mytype* b = new mytype[n];
//    mytype* c = new mytype[n];
//    mytype* d = new mytype[n];
//
//    for (int i = 0; i < n; i++) { a[i] = val[i][1]; }//определяем коэф a
//
//    mytype* h = new mytype[n];
//    mytype* g = new mytype[n];
//
//    mytype* diag3 = new mytype[n - 1];
//    mytype* diag2 = new mytype[n - 1];
//    mytype* diag1 = new mytype[n - 1];
//    mytype* right = new mytype[n - 1];
//
//    for (int i = 0; i < n; i++) { h[i] = val[i + 1][0] - val[i][0];  g[i] = (val[i + 1][1] - val[i][1]) / h[i]; }
//    for (int i = 0; i < n - 1; i++) {
//        diag2[i] = 2 * (h[i + 1] + h[i]);
//        diag1[i] = h[i];
//        diag3[i] = h[i + 1];
//        right[i] = 3 * (g[i + 1] - g[i]);
//    }
//
//    diag1[0] = 0.; diag3[n - 2] = 0.;
//
//    right = progon(n - 1, diag1, diag2, diag3, right);
//
//    c[0] = 0.0;//
//
//
//    for (int i = 1; i < n; i++) { c[i] = right[i - 1]; };
//
//    for (int i = 0; i < n - 1; i++) {
//        b[i] = g[i] - (c[i + 1] + 2 * c[i])* h[i] / 3;
//        d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
//    }
//
//    b[n - 1] = g[n - 1] - 2 * c[n - 1] * h[n - 1] / 3;
//    d[n - 1] = -c[n - 1] / (3 * h[n - 1]);
//
//
//    mytype* S = new mytype[N + 1];
//
//    for (int j = 0, i = 0; j < n; j++)
//        for (; (i < N + 1) && (x[i] <= val[j + 1][0]); i++)
//            S[i] = a[j] + b[j] * (x[i] - val[j][0]) + c[j] * pow((x[i] - val[j][0]), 2) + d[j] * pow((x[i] - val[j][0]), 3);
//
//    delete[] a;
//    delete[] b;
//    delete[] c;
//    delete[] d;
//
//    delete[] h;
//    delete[] g;
//
//    delete[] diag1;
//    delete[] diag2;
//    delete[] diag3;
//    delete[] right;
//
//    return S;
//}
//
std::vector<double> Spline(int nbr, int N, mytype** val, mytype** x) {
    std::vector<double> a(nbr - 1);
    for (size_t i = 0; i < nbr; ++i) {
        a[i] = val[i][1];
    }

     std::vector<double> c(nbr);
     std::vector<double> b(nbr - 1);
     std::vector<double> d(nbr - 1);

    std::vector<double> h(nbr - 1);
    std::vector<double> g(nbr - 1);
    std::vector<double> S(N);

    for (size_t i = 0; i < nbr; ++i) {
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

    for (size_t i = 1; i <= nbr - 1; ++i) {
        c[i] = ksi[i] * c[i - 1] + eta[i];
    }

    c[nbr - 1] = 0;

    for (size_t i = 0; i <= nbr - 1; ++i) {
        b[i] = g[i] - (c[i + 1] + 2 * c[i]) * h[i] / 3;
        d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
    }

    for (int j = 0; j != nbr+1; j++) {
        for (int i = 0; i != N; i++) {
            if (x[i][0] >= val[j][0] && x[i][0] <= val[j+1][0])
             S[i] = a[j] + b[j] * (x[i][0] - val[j][0]) + c[j] * pow((x[i][0] - val[j][0]), 2) + d[j] * pow((x[i][0] - val[j][0]), 3);
        }
    }
    return(S);
}
