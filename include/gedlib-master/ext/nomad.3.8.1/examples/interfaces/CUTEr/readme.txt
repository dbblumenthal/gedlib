
This example illustrates how to use NOMAD with a CUTEr test problem,
given in a SIF files.

You need to have CUTEr installed, including a SIF decoder.
Everything is available at http://hsl.rl.ac.uk/cuter-www/.

Once the CUTEr installation is complete, define these environment variables:

(bash)
export CUTER=/home/user_name/CUTEr
export MYCUTER=$CUTER/CUTEr.***** (depends on specified options during install)
export SIFDEC=/home/user_name/sifdec
export MYSIFDEC=$SIFDEC/SifDec.****  (depends on specified options during install)

Then use the sifdecode program on the problem file:
$MYSIFDEC/bin/sifdecode PROBLEM.SIF


Compile with the C wrapper 'bb.c' by running the script 'compile'
(it uses gfortran and gcc, but you can use other C and Fortran compilers).

The black-box executable 'bb.exe' should be created.

You can test it with the command 'bb.exe x0.txt', and run nomad with the command
'nomad.exe param.txt'.
