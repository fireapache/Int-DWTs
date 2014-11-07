#include "tests.h"

int main(int argc, char **argv)
{
	if (argc > 1)
	{
		test7(argv[1]);
	}
	else cout << "Enter with <ppmfilepath>" << endl;

	return 0;
}
