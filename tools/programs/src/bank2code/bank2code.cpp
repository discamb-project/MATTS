#include "discamb/BasicUtilities/string_utilities.h"

#include <iostream>
#include <fstream>

using namespace discamb;
using namespace std;


int main(int argc, char *argv[])
{

    try 
    {
        if (argc != 2)
            on_error::throwException("expected bank file name as an argument, output will be placed in file 'output'", __FILE__, __LINE__);
		ofstream out("out");
		ifstream in(argv[1]);
		string line;
		while (in.good())
		{
			getline(in, line);
			if (!line.empty())
				if (line[0] == '#')
					continue;
			
			out << "\"" << line << "\",\n";
		}
		in.close();
		out.close();

    }
    catch (exception &e)
    {
        cout << e.what() << endl;
    }
}
