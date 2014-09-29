#include "tests.h"

int main(int argc, char **argv)
{
	if (argc > 1)
	{
		double p = (double)atoi(argv[1]) / 100.0;

		test6(p);
	}
	else cout << "Enter with <energytocompress 0-100>" << endl;
	
	return 0;
}
