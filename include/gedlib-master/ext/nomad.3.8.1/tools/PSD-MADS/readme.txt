Procedure for Unix/Linux/mac OSX

1- Set the NOMAD_HOME environment variable
2- Compile with the command make
3- Go to one of the problem directory
4- Compile the black box (g++ -o bb.exe bb.cpp)
5- Run psdmads  (mpirun -np 3 ../../psdmads.exe param.txt) 

Procedure for Windows

1- Open de nomad MicrosoftVisualStudio solution located in %NOMAD_EXAMPLES%\VisualStudio
2- Build the solution
3- In VisualStudio, add to the solution the existing project %NOMAD_EXAMPLES%\VisualStudio\psdmads.vcxproj
4- Build this new project
5- Open a VisualStudio command shell
6- Change the directory to %NOMAD_EXAMPLES%\examples\PSD-MADS\problems\G2_10
7- Build the blackbox executable with the command: cl bb.cpp
8- Start the optimization: mpiexec -n 4 ..\..\psdmads.exe param.txt 20 3