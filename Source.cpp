#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include "functions.h"
#include "Lagrange.h"

int	main()
{
    mytype	a, b;
    mytype	**mesh;
    mytype	**val;
    mytype	*lagr;
    std::vector<double> sp;
    std::ofstream fout;

    vars(TEST, a, b);
    print(TEST, a, b);
//    mesh = cheb_mesh(TEST, NBR, a, b);
//    print_mesh(mesh, NBR);
//    val = uniform_mesh(TEST, SIZE, a, b);
//	  val = uniform_mesh(TEST, SIZE, a, b + 0.2); //for question #3
//    lagr = lagrange_val(mesh, NBR, TEST, SIZE, val);
//
//
//    fout.open("mass.txt");
//    fout.trunc;
//    for (int i = 0; i <= SIZE; i++)
//        fout << val[i][0] << '\t' << val[i][1] << '\t' << lagr[i] << std::endl;
//    fout.close();
//    std::cout << "Nbr of nodes in Chebyshev mesh for 0.01 is " << find_node(TEST, 0.0001, a, b) << std::endl;
//    std::cout << "Extrapolation error with Lagrange polynom is " << std::scientific << fabs(my_func(2.2, TEST) - lagrange_one(mesh, NBR, TEST, 2.2)) << std::endl; //for question #3
//    std::cout << "Interpolarion error norme " << std::scientific << lagr_norm(lagr, SIZE, val) << std::endl;
    val = uniform_mesh(TEST, NBR, a, b);
    mesh = uniform_mesh(TEST, SIZE, val[0][0], val[NBR][0]);
    sp = Spline(NBR, SIZE, val, mesh);
    std::cout << "Interpolarion error norme " << std::scientific << lagr_norm(sp, SIZE, mesh) << std::endl;
    fout.open("splain.txt");
    fout.trunc;
    for (int i = 0; i < SIZE; i++)
        fout << mesh[i][0] << '\t' << mesh[i][1] << '\t' << sp[i] << std::endl;
    fout.close();

  //  system("pause");
    return (0);
}