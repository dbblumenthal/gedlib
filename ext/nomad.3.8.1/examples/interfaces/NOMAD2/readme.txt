The reader.cpp file is a program which can help you to use NOMAD version 3
with problems for the version 2.

In this part, we assume that the directory in which this readme.txt file is located is:

   $nomad2_interface

"./HS23" is an example of input directory.

We also assume that the NOMAD v2 black-box executables are compiled.
In the example directory, a 'compile' script in given.

     HOW TO USE READER.CPP:

First, you have to compile it. For example, use the following command:

   cd $nomad2_interface
   g++ -o reader.exe reader.cpp

Then, you can run the executable with one or two argument(s),
like this:

   ./reader.exe $dir/description_file ($dir/parameter_file)

description_file and parameter_file are two files used by NOMAD (version 2) to solve a problem.
You have to put all the required files in the same directory ("dir" in our example) :
the description_file and the parameter file,
the executables files (for constraints and objective function), the input...

Somes files are created when the reader is executed:

* bb.cpp: This is the Black Box used by NOMAD.
You have to compile it with:

   g++ -o bb.exe bb.cpp

Warning: The name of this executable must be "bb.exe",
or you have to change this name in the param.txt file (after "BB_EXE").

* param.txt: This file contains all the parameters used by NOMAD.


     HOW TO USE NOMAD:

You just have to call:

   $NOMAD_HOME/bin/nomad.exe $nomad2_interface/param.txt


WARNING: this program works under LINUX; with WINDOWS, you may have to manually change the executable paths directly into the created source file bb.cpp.
