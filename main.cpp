#include "tests.h"

int main(int argc, char **argv)
{
	if (argc == 3)
	{
		for (int i = 0; i < 4; ++i)
		{
			test10(argv[1], i, atoi(argv[2]));
		}
	}
	else
	{
		cout << "Enter <ppmFilePath> <levels>" << endl;
		cout << endl;
	}

	return 0;
}
