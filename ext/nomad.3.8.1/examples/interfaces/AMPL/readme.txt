Instructions to use AMPL models with NOMAD:

1. Install the AMPL interface library.

   This interface can be downloaded at
   http://www.netlib.org/ampl/solvers
   (for a simple installation and the generation of dynamic libraries).

   The example given in the NOMAD package uses the dynamic library version.

   Complete instructions on the AMPL interface are available at:
   http://ampl.com/REFS/HOOKING/index.html

2. Put the AMPL files in directory MODELDIR;
   Suppose that the name of your model is MODELDIR/modelname.

3. Create a .nl file from your model:

   Open AMPL, load your model with the command 'model modelname;'
   (and possibly your data file), and create the nl file with the
   command 'write gmodelname;'
   (use the 'g' character to generate an ascii nl file).

4. Create x0 and bound files from the AMPL model.

   Warning: check the variable indexing that is not necessarily
   the same than in the model file.

5. Edit the file 'bb.c' and define the appropriate string for MODEL_NAME.

6. Compile bb.c with the command 'make'; this creates 'bb.exe'
   (be sure that your makefile uses the right path for your AMPL interface
    library).

   You can test the bb.exe program by defining DISPLAY in bb.c

7. Write your NOMAD parameters file; pay attention to
   the number of variables, the constraints, starting point
   and the bounds.

8. Run NOMAD.


Many thanks to Dominique Orban, Quentin Reynaud, and Anthony Guillou.
