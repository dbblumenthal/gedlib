This example illustrates an application of the NOMAD library in order to graphically represent the outputs of an execution. It must not be considered as a graphical user interface of the code, as it is very simple and incomplete for such a purpose. The objective of this example is simply to illustrate a library usage of NOMAD and some links that can be made with another language such as JAVA.

The example has been developed by Quentin Reynaud from ISIMA. It has been tested on a MAC with the C++ and JAVA compilers included in the XCODE package.

The example is composed of two programs: The first is a JAVA application that runs in background and periodically checks NOMAD outputs in order to draw them. The other program uses the NOMAD library in order to optimize single or bi-objective problems and to communicate outputs to the JAVA application.

Open a console, and from the plot directory, go to GUI/src. Compile the JAVA code with the command 'javac *.java'. This generates .class files that you have to move in plot/GUI/bin. Go to this last directory, and execute the JAVA application in the background with the command 'java Prog'.

To execute the NOMAD code associated with the JAVA graphical interface, go to the plot directory and compile with 'make' (NOMAD must be installed and the $NOMAD_HOME environment variable must be defined).
This generates the executable file 'nomad_plot.exe', and you can run it on the two test problems located in the plot/problems directory.

For examples, type './nomad_plot.exe ./problems/01/param_single_obj.txt' to execute on a single-objective problem. Try also the bi-objective version of the same problem with the parameters file param_bi_obj, or 02 instead of 01 for another bi-objective execution. Black-box executables must have been compiled beforehand with command 'g++ -o bb.exe bb.cpp'.