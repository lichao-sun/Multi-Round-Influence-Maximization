#include <iostream>
#include "mi_command_line.h"


// #define _TEST_

#ifdef _TEST_

#include "tests/tester_entry.h"
int main(int argc, char* argv[])
{
	TestEntry::Main(argc, argv);
}

#else



int main(int argc, char* argv[])
{
	try
	{
        MICommandLine cmd;
		return cmd.Main(argc, argv);
	}
	catch (std::exception e) 
	{
		std::cerr << e.what() << std::endl;
	}
	return -1;
}

#endif


