#include "tests.h"

int main(int argc, char **argv)
{
	if (argc == 2)
	{
		fundamentalTest(atoi(argv[1]));
	}
	else cout << "Enter <MatrixOrder> <MaxRandomValue> <LevelsOfTransformation>" << endl;

	return 0;
}
