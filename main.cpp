#include "tests.h"

int main(int argc, char **argv)
{
	if (argc > 2)
	{
		double p = (double)atoi(argv[2]) / 100.0;

		test3(argv[1], p);
	}
	else cout << "Enter with <ppmfilepath> <energytokeep 0-100>" << endl;
	
	return 0;
}
