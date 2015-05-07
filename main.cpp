#include "tests.h"

int main(int argc, char **argv)
{
	if (argc > 2)
	{
		test8(atoi(argv[1]), atoi(argv[2]));
	}
	else cout << "Enter with a test number" << endl;

	return 0;
}
