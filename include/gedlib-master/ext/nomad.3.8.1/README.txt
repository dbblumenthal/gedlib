################################################################################
#                                                                              #
#                                    README                                    #
#                                                                              #
#------------------------------------------------------------------------------#
#  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct search -             #
#          version 3.8.1                                                       #
#                                                                              #
#  NOMAD - version 3.8.1 has been created by                                   #
#                 Charles Audet        - Ecole Polytechnique de Montreal       #
#                 Sebastien Le Digabel - Ecole Polytechnique de Montreal       #
#                 Christophe Tribes    - Ecole Polytechnique de Montreal       #
#                                                                              #
#  The copyright of NOMAD - version 3.8.1 is owned by                          #
#                 Sebastien Le Digabel - Ecole Polytechnique de Montreal       #
#                 Christophe Tribes    - Ecole Polytechnique de Montreal       #
#                                                                              #
#  NOMAD v3 has been funded by AFOSR, Exxon Mobil, Hydro Québec, Rio Tinto     #
#  and IVADO.                                                                  #
#                                                                              #
#  NOMAD v3 is a new version of NOMAD v1 and v2. NOMAD v1 and v2 were created  #
#  and developed by Mark Abramson, Charles Audet, Gilles Couture, and John E.  #
#  Dennis Jr., and were funded by AFOSR and Exxon Mobil.                       #
#                                                                              #
#  Contact information:                                                        #
#    Ecole Polytechnique de Montreal - GERAD                                   #
#    C.P. 6079, Succ. Centre-ville, Montreal (Quebec) H3C 3A7 Canada           #
#    e-mail: nomad@gerad.ca                                                    #
#    phone : 1-514-340-6053 #6928                                              #
#    fax   : 1-514-340-5665                                                    #
#                                                                              #
#  This program is free software: you can redistribute it and/or modify it     #
#  under the terms of the GNU Lesser General Public License as published by    #
#  the Free Software Foundation, either version 3 of the License, or (at your  #
#  option) any later version.                                                  #
#                                                                              #
#  This program is distributed in the hope that it will be useful, but WITHOUT #
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       #
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License #
#  for more details.                                                           #
#                                                                              #
#  You should have received a copy of the GNU Lesser General Public License    #
#  along with this program. If not, see <http://www.gnu.org/licenses/>.        #
#                                                                              #
#  You can find information on the NOMAD software at www.gerad.ca/nomad        #
#------------------------------------------------------------------------------#


DESCRIPTION :

NOMAD is a C++ implementation of the Mesh Adaptive Direct Search (MADS)
algorithm, designed for constrained optimization of black-box functions.


WEB PAGE :

https://www.gerad.ca/nomad/


FURTHER INSTRUCTIONS :

Please visit the web page for futher instruction on the following:

* Downloading, configuring, compiling, and installing NOMAD
* Using NOMAD and setting the parameters
* Reports on NOMAD
* How to report bugs and make enhancement requests
* And more...


BATCH OR LIBRARY MODE :

NOMAD is designed to be used in two different modes : batch and library.
The batch mode is intended for a basic and simple usage of the MADS method,
while the library mode allows more flexibility.
For example, in batch mode, users must define their separate black-box program,
that will be called with system calls by NOMAD.
In library mode, users can define their black-box function as C++ code
that will be directly called by NOMAD, without system calls and temporary files.


EXECUTABLES :

On Windows and Mac OS X platforms, NOMAD batch mode executable is located in
directory $NOMAD_HOME/bin or %NOMAD_HOME%\bin.
In order to avoid compiling the code, you can simply use this executable.

COMPILATION :

On Linux, Unix, and Mac OS X, NOMAD can be compiled using the makefile located
in $NOMAD_HOME/bin. Alternatively, using the script install.sh located in
directory $NOMAD_HOME/install will install all NOMAD binaries.

On Windows, NOMAD can be compiled using the Microsoft Visual Studio projet
located in the %NOMAD_EXAMPLES%\VisualStudio directory.


HOW TO EXECUTE NOMAD :

For informations about the execution of NOMAD, please read the user guide :

$NOMAD_HOME/doc/user_guide.pdf or %NOMAD_HOME%\doc\user_guide.pdf

or

https://www.gerad.ca/NOMAD/Downloads/user_guide.pdf
