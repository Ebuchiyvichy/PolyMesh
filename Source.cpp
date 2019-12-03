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
	mytype  *sp;
	std::vector<mytype>	lagr;
	std::ofstream fout;
	
	vars(TEST, a, b);
	print(TEST, a, b);
	mesh = cheb_mesh(TEST, NBR+1, a, b);
	print_mesh(mesh, NBR);
	val = uniform_mesh(TEST, SIZE, a, b);
	lagr = lagrange_val(mesh, NBR, TEST, SIZE, val);
	fout.open("mass.txt");
	fout.trunc;
	for (int i = 0; i < SIZE; i++)
		fout << val[i][0] << '\t' << val[i][1] << '\t' << lagr[i] << std::endl;
	fout.close();
	std::cout << "Nbr of nodes in Chebyshev mesh for 0.01 is " << find_node(TEST, 0.01, a, b) << std::endl;
	// построение сплайн интерполяции
	sp = splain(mesh, NBR, TEST, SIZE, val);
    for (int i = 0; i < SIZE; i++)
        fout << val[i][0] << '\t' << val[i][1] << '\t' << sp[i] << std::endl;
    std::cout << "Error in splain \t" << errors(SIZE, sp, values(SIZE, TEST, val));

	//system("pause");
	return (0);
}