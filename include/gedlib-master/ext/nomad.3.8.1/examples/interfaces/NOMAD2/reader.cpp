/*-------------------------------------------------------------------------------------*/
/*  This program read informations files from NOMAD version 2 and create the           */
/*  parameters files used by NOMAD version 3                                           */
/*-------------------------------------------------------------------------------------*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

// This method fills the array 'coord' with numbers read from the file 'fin'.
void fillArray(ifstream & fin, double * coord, int dimension)
{
	char ch;
	char line[80];
	int i = 0;
  // We go in the loop: the next character is extracted from the file.
	while ((i < dimension) && (fin.get(ch)))
	{
      // '%' marks a comment line.
		if (ch == '%')
	// The line is read.
			fin.getline(line, 80);
		else   // If the end of line is reached we do nothing.
			if (ch == '\n') {}
		else   // If it's not a comment line or an end of line, it's a number.
		{
	    // The character is reinserted in the file.
			fin.putback(ch);
	    // The number is read and put in the 'coord' array.
			fin >> coord[i++];
		}
	}
}

// This method fills the array 'coord' with string read from the file 'fin'.
void fillArrayString(ifstream & fin, string * coord, int dimension)
{
	char ch;
	char line[80];
	int i = 0;
	
  // We go in the loop: the next character is extracted from the file.
	while ((i < dimension) && (fin.get(ch)))
	{
      // '%' marks a comment line.
		if (ch == '%')
	// The line is read.
			fin.getline(line, 80);
		else   // If the end of line is reached we do nothing.
			if (ch == '\n') {}
		else   // If it's not a comment line or an end of line, it's a string.
		{
	    // The character is reinserted in the file.
			fin.putback(ch);
	    // The string is read and put in the 'coord' array.
			fin >> coord[i++];
		}
	}
}

// This method remplace the path from s
string remplace (string s)
{
	char c = ' ';
	char s2[255];
	int 	i = 0,
			j = 0;

	while(c != '\0')
	{
		c = s[i++];
		if (c == '/') j = 0;
		s2[j++] = c;
	}
	s2[j] = '\0';
	return (s2);
}

int main(int argc, char** argv)
{
	if (argc > 1 && argc <= 3)
	{
		// Reading the description file
		ifstream desc(argv[1], ios::in);
		
		if(desc)
		{
			string		mot = "",
					bounds = "",
					start = "",
					results = "",
					blackBox = "",
					input = "",
					truth = "",
					constrains = "",
					surrogate = "",
					caches = "";
						
			char	path[50];
						
			bool	isBounds = false,
					isBB = false,
					isSurrogate = false,
					isCaches = false,
					pollComplete,
					speculativeSearch,
					initialComplete,
					isParam[21];
					
			int 	dim = -1,
					cons = -1,
					coarseningExponent,
					refiningExponent,
					filterNorm,
					seed,
					truthEvals,
					newTruthEvals,
					iterations,
					iterativePoints,
					iterativeSearch,
					initialSearch,
					initialPoints,
					pollDirections,
					displayFactor;
			
			double	pollBasis,
						frameCenterTrigger,
						hmin,
						hmax,
						pollSizeTerm;

			// Initialize the presence of the parameters
			for (int i = 0; i < 21; ++i) isParam[i] = false;

			// Searching some parameters
			while (!desc.eof())
			{
				desc >> mot;
				
				if (mot == "DIMENSION")
				{
					desc >> dim;
				}
				else if (mot == "USE_CACHES")
				{
					desc >> mot;
					if (mot == "1") isCaches = true;
				}
				else if (mot == "CACHE_FILE")
				{
					desc >> caches;
					caches = remplace (caches);
				}
				else if (mot == "GEN_CONS_NB")
				{
					desc >> cons;
				}
				else if (mot == "USE_BOUNDS")
				{
					desc >> mot;
					if (mot == "1") isBounds = true;
				}
				else if (mot == "BOUNDS_FILE")
				{
					desc >> bounds;
				}
				else if (mot == "START_PT_FILE")
				{
					desc >> start;
				}
				else if (mot == "RESULTS_FILE")
				{
					desc >> results;
					results = remplace (results);
				}
				else if (mot == "USE_BLACK_BOXES")
				{
					desc >> mot;
					if (mot == "1") isBB = true;
				}
				else if (mot == "INPUT_FILE")
				{
					desc >> input;
					input = remplace (input);
				}
				else if (mot == "TRUTH_EXE")
				{
					desc >> truth;
					truth = remplace (truth);
				}
				else if (mot == "GEN_CONS_FILE")
				{
					desc >> constrains;
				}
				else if (mot == "USE_SURROGATE")
				{
					desc >> mot;
					if (mot == "1") isSurrogate = true;
				}
				else if (mot == "SURROGATE")
				{
					desc >> surrogate;
					surrogate = remplace (surrogate);
				}
				
				else desc >> mot;
			}

			desc.close();

			// Path of the input file
			if (input != "")
			{
				char c = ' ';
				int 	i = 0,
	 					tmp = 0;
				while (c != '\0')
				{
					c = argv[1][i];
					if (c == '/') tmp = i;
					path[i++] = c;
				}
				path[tmp] = '\0';
			}
			else
			{
				cout<< "Error : you must give an input file."<<endl;
				return 1;
			}
			
			if (cons == -1)
			{
				cout <<"Error : you must give the number of constraints." << endl;
				return 1;
			}

			string * cons_val = new string [cons];

			// Constraints files
			if (constrains != "")
			{
				ifstream cons_file(constrains.c_str(), ios::in);
				fillArrayString (cons_file, cons_val, cons);
				cons_file.close();
				
				for (int i = 0; i < cons ; ++i)
				{
					cons_val[i] = remplace (cons_val[i]);
				}
			}

			// Creation of the Black Box
			if (isBB)
			{
				ofstream bb("bb.cpp", ios::trunc);
				bb << "// This file is generated by the reader program.\n\n";
				bb << "#include <iostream>\n#include <stdlib.h>\n#include <fstream>\nusing namespace std;\n";
				bb << "int main(int argc, char** argv)\n{\n";
				bb << "if (argc == 2){\n";
				bb << "\tdouble entree[";
				bb << dim;
				bb << "];\n\tint i = 0;\n\tifstream file(argv[1], ios::in);\n";
				bb << "\tif(file){\n";
				bb << "\t/* Reading the point */\n\t\twhile ( !file.eof()) file >> entree[i++];\n";
				bb << "\t\tfile.close();\n";
				bb << "\t\tofstream fic(\"input.txt\", ios::trunc);\n";
				bb << "\t\tfic.precision(15);\n";
				bb << "\t\tfor (i = 0; i < ";
				bb << dim;
				bb << "; ++i){\n";
				bb << "\t\tfic << entree[i] << \" \";}\n";
				bb << "\t\tfic.close();";
				bb << "system(\"" << path << truth << "\");\n";
				bb << "\t\tcout << endl;\n";
				for (int i = 0; i < cons; ++i)
				{
					bb << "\t\tsystem(\"" << path << cons_val[i] << "\");\n";
					bb << "\t\tcout << endl;\n";
				}
				bb << "\t\t}\n\t}\n";
				bb << "\telse{\n\t\tcout << \"Error.\" << endl;\n";
				bb << "\t\treturn 1;}\n";

				bb << "\n\treturn 0;\n}"<<endl;
				
				bb.close();
			}
			else
			{
				cout << "Error : you must use a Black Box."<<endl;
				return 1;
			}
			
		
			delete [] cons_val;


			// Reading the parameters file
			ifstream param(argv[2], ios::in);
			
			if (param)
			{
				char ch;
				char line[100];
				
				while (param.get(ch))
				{
					if(ch == '#') param.getline(line, 100);
					else if (ch == '\n') {}
					else 
					{
						param.putback(ch);
						param >> mot;

						// Searching the parameters
						if (mot == "POLL_BASIS")
						{
							param >> pollBasis;
							isParam[0] = true;
						}
						else if (mot == "COARSENING_EXPONENT")
						{
							param >> coarseningExponent;
							isParam[1] = true;
						}
						else if (mot == "REFINING_EXPONENT")
						{
							param >> refiningExponent;
							isParam[2] = true;
						}
						else if (mot == "RANDOM_SEED")
						{
							param >> seed;
							isParam[3] = true;
						}
						else if (mot == "POLL_COMPLETE")
						{
							param >> pollComplete;
							isParam[4] = true;
						}
						else if (mot == "POLL_DIRECTIONS")
						{
							param >> pollDirections;
							isParam[5] = true;
						}
						else if (mot == "INITIAL_SEARCH")
						{
							param >> initialSearch;
							isParam[6] = true;
						}
						else if (mot == "INITIAL_COMPLETE")
						{
							param >> initialComplete;
							isParam[7] = true;
						}
						else if (mot == "INITIAL_POINTS")
						{
							param >> initialPoints;
							isParam[8] = true;
						}
						else if (mot == "ITERATIVE_SEARCH")
						{
							param >> iterativeSearch;
							isParam[9] = true;
						}
						else if (mot == "ITERATIVE_POINTS")
						{
							param >> iterativePoints;
							isParam[10] = true;
						}
						else if (mot == "SPECULATIVE_SEARCH")
						{
							param >> speculativeSearch;
							isParam[11] = true;
						}
						else if (mot == "POLL_SIZE_TERM")
						{
							param >> pollSizeTerm;
							isParam[12] = true;
						}
						else if (mot == "ITERATIONS")
						{
							param >> iterations;
							isParam[13] = true;
						}
						else if (mot == "TRUTH_EVALS")
						{
							param >> truthEvals;
							isParam[14] = true;
						}
						else if (mot == "NEW_TRUTH_EVALS")
						{
							param >> newTruthEvals;
							isParam[15] = true;
						}
						else if (mot == "HMAX")
						{
							param >> hmax;
							isParam[16] = true;
						}
						else if (mot == "HMIN")
						{
							param >> hmin;
							isParam[17] = true;
						}
						else if (mot == "FILTER_NORM")
						{
							param >> filterNorm;
							isParam[18] = true;
						}
						else if (mot == "FRAME_CENTER_TRIGGER")
						{
							param >> frameCenterTrigger;
							isParam[19] = true;
						}
						else if (mot == "DISPLAY_FACTOR")
						{
							param >> displayFactor;
							isParam[20] = true;
						}
					}	
				}

				param.close();
			}
				
				// Writing the output file
				ofstream fic("param.txt", ios::trunc);
				
				if (fic)
				{
					// Parameters from description file
					if (dim != -1)
					{
						fic << "DIMENSION\t\t" << dim  << "\n\n";
					}
					else
					{
						cout << "Error : DIMENSION is missing."<<endl;
						return 1;
					}
					
					fic << "BB_OUTPUT_TYPE\t\tOBJ ";
					for (int i = 0; i < cons; ++i)
					{
						fic << "CSTR ";
					}
					fic << "\n# Output type, see the documentation for details\n\n";
					
					if (isCaches)
					{
						if (caches != "")
						{
							fic << "CACHE_FILE\t\t" << caches << "\n\n";
						}
						else
						{
							cout << "Error : you are using the caches, but CACHE_FILE is not defined."<<endl;
							return 1;
						}
					}

					fic << "BB_EXE\t\tbb.exe\t#'bb.exe' is a C++ executable used for run NOMAD with NOMAD v3.\n\n";

					if (start != "")
					{
						ifstream start_file(start.c_str(), ios::in);
						if (start_file)
						{
							double * start_val = new double [dim];
							fillArray (start_file, start_val, dim);
							start_file.close();
							
							fic << "X0 ( ";
							for (int i = 0; i < dim; ++i)
							  fic << start_val[i] << " ";

							delete [] start_val;

							fic << ")\t# Starting point\n\n";
						}
						else cout << start << "can't be open."<<endl;
					}

					if (isBounds)
					{
						if (bounds != "")
						{
							ifstream bounds_file(bounds.c_str(), ios::in);
							if (bounds_file)
							{
							  double * bounds_val = new double [2*dim];
							  fillArray (bounds_file, bounds_val, 2*dim);
								
							  bounds_file.close();

							  fic << "LOWER_BOUND ( " ;
							  for (int i = 0; i < dim; ++i)
							    fic << bounds_val[i] << " ";
							  fic << ")\n\n";
								
							  fic << "UPPER_BOUND ( " ;
							  for (int i = dim; i<2*dim; ++i)
							    fic << bounds_val[i] << " ";
							  fic << ")\n\n";
							  delete [] bounds_val;
							}
							else cout << bounds << "can't be open."<<endl;
						}
						else cout << "Warning : bounds are declared but not defined"<<endl;
					}
					
					if (results != "")
					{
						fic << "SOLUTION_FILE\t\t" << path << results << "\n\n";
					}

					if (isSurrogate)
					{
						if (surrogate != "")
						{
							fic << "SGTE_EXE\t\t" ;
							fic << path << surrogate <<"\n\n";
						}
						else
						{
							cout<< "Warning : you are using surrogate, but no surrogate function is given."<<endl;
						}
					}	
					// Parameters from parameters file
					if (isParam[0])
					{
						fic << "MESH_UPDATE_BASIS\t\t" << pollBasis << "\n\n";
					}
					if (isParam[1])
					{
						fic << "MESH_COARSENING_EXPONENT\t\t" << coarseningExponent << "\n\n";
					}
					if (isParam[2])
					{
						fic << "MESH_REFINING_EXPONENT\t\t" << refiningExponent << "\n\n";
					}
					if (isParam[3])
					{
						fic << "SEED\t\t" << seed << "\n\n";
					}
					if (isParam[4])
					{
						fic << "OPPORTUNISTIC_EVAL\t\t";
						if (pollComplete) fic << "no\n\n";
						else fic << "yes\n\n";
					}
					if (isParam[5])
					{
						fic << "DIRECTION_TYPE\t\t";
						if (!pollDirections) fic << "GPS 2N";
						else if (pollDirections == 1) fic << "GPS N+1";
						else if (pollDirections == 2) fic << "GPS N+1 STATIC UNIFORM";
						else if (pollDirections == 3) fic << "ORTHO 2N";
						else fic << "LT N+1";
						fic << "\n\n";
					}
					if (isParam[6])
					{
						if (initialSearch == 2)
						{
							fic << "LH_SEARCH\t\t";
							if (isParam[8])
							{
								fic << initialPoints << " ";
							}
							else fic << "0 ";
							
							if (isParam[9] && iterativeSearch == 2 && isParam[10])
							{
								fic << "iterativePoints\n\n";
							}
							else fic << "0\n\n";
						}
					}
					if (isParam[7])
					{
						fic << "OPPORTUNISTIC_LH\t\t";
						if (initialComplete) fic << "no\n\n";
						else fic << "yes\n\n";
					}
					if (isParam[11])
					{
						fic << "SPECULATIVE_SEARCH\t\t";
						if (speculativeSearch) fic << "yes\n\n";
						else fic << "no\n\n";
					}
					if (isParam[12] && pollSizeTerm)
					{
						fic << "MIN_POLL_SIZE\t\t" << pollSizeTerm << "\n\n";
					}
					if (isParam[13] && iterations)
					{
						fic << "MAX_ITERATION\t\t" << iterations << "\n\n";
					}
					if (isParam[14] && truthEvals)
					{
						fic << "MAX_BB_EVAL\t\t" << truthEvals << "\n\n";
					}
					if (isParam[15] && newTruthEvals)
					{
						fic << "MAX_SIM_BB_EVAL\t\t" << newTruthEvals << "\n\n";
					}
					if (isParam[16])
					{
						fic << "H_MAX_0\t\t" << hmax << "\n\n";
					}
					if (isParam[17])
					{
						fic << "H_MIN\t\t" << hmin << "\n\n";
					}
					if (isParam[18])
					{
						fic << "H_NORM\t\t";
						if (!filterNorm) fic << "Linf\n\n";
						else if (filterNorm == 1) fic << "L1\n\n";
						else fic << "L2\n\n";
					}
					if (isParam[19])
					{
						fic << "RHO\t\t" << frameCenterTrigger << "\n\n";
					}
					if (isParam[20])
					{
						fic << "DISPLAY_DEGREE\t\t";
						if (displayFactor < 4) fic << displayFactor << "\n\n";
						else fic << "4\n\n";
					}
					
					fic.close();
				}
		}
		else cout <<"Error, cannot open the description.dat file."<<endl;
	}
	else 
            cout << "usage: ./reader description_file parameters_file (optional)\n";
	
	return 0;
}
