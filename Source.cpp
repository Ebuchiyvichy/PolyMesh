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
	mytype	**mesh_lagr;
	mytype	**mesh_splain;
	mytype	**val_lagr;
	mytype	**val_splain;
	mytype	*lagr;
	std::vector<double> sp;
	std::ofstream fout;
	
	vars(TEST, a, b);
	print(TEST, a, b);
	mesh_lagr = cheb_mesh(TEST, NBR, a, b);
//	print_mesh(mesh_lagr, NBR);
	val_lagr = uniform_mesh(TEST, SIZE, a, b);
//	val_lagr = uniform_mesh(TEST, SIZE, a, b + 0.2); //for question #3
	lagr = lagrange_val(mesh_lagr, NBR, TEST, SIZE, val_lagr);
	std::cout << "Interpolarion error norme Lagrange " << std::scientific << lagr_norm(lagr, SIZE, val_lagr) << std::endl;

	mesh_splain = uniform_mesh(TEST, NBR, a, b);
	val_splain = uniform_mesh(TEST, SIZE, mesh_splain[0][0], mesh_splain[NBR][0]);
	sp = Spline(NBR + 1, SIZE, mesh_splain, val_splain);
	std::cout << "Interpolarion error norme splain " << std::scientific << lagr_norm(sp, SIZE, val_splain) << std::endl;
	fout.open("mass.txt");
	fout.trunc;
	for (int i = 0; i <= SIZE; i++)
		fout << val_lagr[i][0] << '\t' << val_lagr[i][1] << '\t' << lagr[i] << '\t' << sp[i] << std::endl;
	fout.close();
//	std::cout << "Nbr of nodes in Chebyshev mesh for 0.01 is " << find_node(TEST, 0.0001, a, b) << std::endl;
//	std::cout << "Extrapolation error with Lagrange polynom is " << std::scientific << fabs(my_func(2.2, TEST) - lagrange_one(mesh_lagr, NBR, TEST, 2.2)) << std::endl; //for question #3
	delete[] mesh_lagr[0];
	delete[] mesh_lagr[1];
	delete[] mesh_splain[0];
	delete[] mesh_splain[1];
	delete[] val_lagr[0];
	delete[] val_lagr[1];
	delete[] val_splain[0];
	delete[] val_splain[1];

	system("pause");
	return (0);
}