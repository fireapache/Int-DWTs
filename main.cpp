#include "tests.h"

int main(int argc, char **argv)
{
	if (argc == 2)
	{
		for (int i = 0; i < 4; ++i)
		{
			test11(argv[1], i);
		}
	}
	else
	{
		cout << "Enter <ppmFilePath>" << endl;
		cout << endl;
	}

	return 0;
}
