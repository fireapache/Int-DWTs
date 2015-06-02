#include "tests.h"

int main(int argc, char **argv)
{
	if (argc > 1)
	{
		fundamentalTest(atoi(argv[1]));
	}
	else cout << "Enter with a test number" << endl;

	return 0;
}
