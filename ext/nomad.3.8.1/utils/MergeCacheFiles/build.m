%% GERAD NOMAD Build for Matlab

% This file will help you compile NOMAD for use with MATLAB.

% To recompile you will need to get / do the following:

% 1) Get NOMAD
% NOMAD is available from https://www.gerad.ca/nomad
% Complete the download form then download the latest version. Define the
% $NOMAD_HOME environment variable.

% 2) Start Matlab and go into $NOMAD_HOME/examples/interfaces/Matlab_Mex
% The NOMAD MEX Interface is a simple MEX interface written to use NOMAD.

% 3) Compile the MEX File
% The code below will build the NOMAD MEX file and set the Matlab path.

clear nomad

% Default library names
nameLibNomad = 'libnomad.so';

% Default update LDFLAGS
updateLDFLAGS= '';

% Current directory
cdir = cd;

% Check and set nomad_home and create variables for path
clear nomad_home nomad_src;
nomad_home = getenv('NOMAD_HOME');


if ( length(nomad_home)<1 )
% Get a default directory for NOMAD_HOME
cd ..
cd ..
nomad_home = cd;

if ( ~ exist(nomad_home,'dir') )
error('The default NOMAD_HOME directory does not exist. Please provide a correct value for the NOMAD_HOME variables with the command setenv(''NOMAD_HOME'',ARG1)');
end

setenv('NOMAD_HOME',nomad_home);
warning('The NOMAD_HOME variable is set with its default value %s. The default can be replaced by using the command setenv(''NOMAD_HOME'',ARG1)! before runner the GERAD_NOMAD_build command.',nomad_home);
cd(cdir);
end
nomad_src=[nomad_home filesep 'src' filesep];
nomad_lib=[nomad_home filesep 'lib' filesep];
nomad_bin=[nomad_home filesep 'bin' filesep];
nomad_sgtelib_src=[nomad_home filesep 'ext' filesep 'sgtelib' filesep 'src'];

switch(computer)
case 'PCWIN'
%libdir = ' -Lwin32\';
%nameLibNomad = 'nomad.dll';
error('The installation script for Windows is not yet available.');
case 'PCWIN64'
%libdir = ' -Lwin64\';
%nameLibNomad = 'nomad.dll';
error('The installation script for Windows is not yet available.');
case 'GLNX86'
libdir = 'glnx86/';
updateLDFLAGS = 'LDFLAGS=''$LDFLAGS -Wl,-rpath,''''$ORIGIN/../lib/'''' '' ';
case 'GLNXA64'
libdir = 'glnxa64/';
updateLDFLAGS = 'LDFLAGS=''$LDFLAGS -Wl,-rpath,''''$ORIGIN/../lib/'''' '' ';
case 'MACI64'
libdir = 'maci64/';
end

%Get NOMAD Libraries
post = [' -I.  -I' nomad_src ' -I' nomad_sgtelib_src ' -lm -lut -lnomad -L ' nomad_lib ' -output mergeM'];

%Compile & Move
pre = ['mex -v -largeArrayDims ' updateLDFLAGS ' mergeM.cpp ' ];

try
    eval([pre post])
    if ( strcmp(computer,'MACI64') )
        install_cmd = ['install_name_tool -change libnomad.so ' nomad_lib 'libnomad.so mergeM.mexmaci64' ]; 
        system(install_cmd);
    end
    clear nomad_home nomad_src nomad_sgtelib_src cdir post pre libdir;
    fprintf('Done!\n');
catch ME
	clear nomad_home nomad_src nomad_sgtelib_src cdir post pre libdir;
    error('opti:nomad','Error Compiling NOMAD!\n%s',ME.message);
end
fprintf('------------------------------------------------\n');
