#include "tests.h"

int main(int argc, char **argv)
{
	if (argc == 4)
	{
		test9(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
	}
	else cout << "Enter <MatrixOrder> <MaxRandomValue> <LevelsOfTransformation>" << endl;

	return 0;
}
