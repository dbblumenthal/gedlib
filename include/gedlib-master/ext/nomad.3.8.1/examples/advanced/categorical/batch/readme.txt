This example illustrates how to run the batch mode of NOMAD
on a problem with categorical variables.

In addition to the usual blackbox executable, the user must provide
an executable for the neighborhoods.
This executable takes a point as input and displays a list of neighbors.
The number of variables must remain constant.

For example:

  neighbors.exe x.txt

  with x.txt containing '0 100 1 100'

  displays:

  '2 100 1 100
   0 100 2 100'

  (two neighbors).

To run the example:

1. compile blackbox executable : g++ -o bb.exe bb.cpp -O3
2. compile neighbors executable: g++ -o neighbors.exe neighbors.cpp -O3
3. run nomad: nomad param.txt

See more details in the NOMAD user guide.

NOTE: another neighbor executable is given as a python script.
Replace NEIGHBORS_EXE neighbors.exe
with    NEIGHBORS_EXE  "$python neighbors.py"
in the parameters file to use it.
