%% GERAD NOMAD Build for Matlab

% This file will help you compile NOMAD for use with MATLAB. 

% To recompile you will need to get / do the following:

% 1) Get NOMAD
% NOMAD is available from https://www.gerad.ca/nomad
% Complete the download form then download the latest version. Define the
% $NOMAD_HOME environment variable.

% 2) Start Matlab and go into $NOMAD_HOME/examples/interfaces/Matlab_Mex
% The NOMAD MEX Interface is a simple MEX interface written to use NOMAD.

% 3) Compile the MEX File by executing this file in Matlab
%
%
% The code below will build the NOMAD MEX file and set the Matlab path. 

clear nomad

% Current directory
cdir = cd;

% Check and set nomad_home and create variables for path
clear nomad_home nomad_src nomad_src_sgtelib;

% Default values
nameLibNomad = '';
updateLDFLAGS= '';
install_name_tool='';

if ( strcmp(computer,'PCWIN64') == 1 || strcmp(computer,'PCWIN32') == 1 )

    nomad_home = getenv('NOMAD_EXAMPLES');

    if ( length(nomad_home) > 1)
        nomad_home = [ nomad_home '\VisualStudio' ];
        warning('The nomad_home variable for Matlab is set to %s (that is NOMAD_EXAMPLES\\VisualStudio).The default can be replaced by using the command setenv(''NOMAD_EXAMPLES'',ARG1)! before running the GERAD_NOMAD_build command.',nomad_home);
        if ( ~isempty( find(isspace(nomad_home),1) ) )
            error('The compilation of Nomad for Matlab must be performed in the NOMAD_EXAMPLES directory. The NOMAD_EXAMPLES directory should not contain empty space. Please consider moving the NOMAD_EXAMPLES directory and reset the NOMAD_EXAMPLES environment variable accordingly.');
        end
        
    else
        cd ..
        cd .. 
        cd ..
        cd 'VisualStudio';    
        nomad_home = cd; 
    
        if ( ~ exist(nomad_home,'dir') )
            error('The default NOMAD_EXAMPLES\VisualStudio directory does not exist. Please make sure that a VisualStudio directory exists.');
        end
        
        cd(cdir);
    end
      
    nomad_src=[nomad_home filesep 'src' filesep];
    nomad_src_sgtelib=[nomad_home filesep 'src_sgtelib' filesep];
    nomad_bin=[nomad_home filesep 'bin' filesep];
    nomad_lib='';

    %Compile & Move (Windows) ---> recompile Nomad and sgtelib
    post = [' -I.  -I' nomad_src ' -I' nomad_src_sgtelib ' -lut -output ' nomad_bin filesep 'nomad.' mexext];
    pre = ['mex -v -largeArrayDims nomadmex.cpp ' nomad_src 'Parameters.cpp ' nomad_src 'Barrier.cpp ' nomad_src 'Cache.cpp '...
    nomad_src 'Cache_File_Point.cpp ' nomad_src 'Cache_Point.cpp ' nomad_src 'Cache_Search.cpp ' nomad_src 'Clock.cpp '...
    nomad_src 'Direction.cpp ' nomad_src 'Directions.cpp ' nomad_src 'Display.cpp '...
    nomad_src 'Double.cpp ' nomad_src 'Eval_Point.cpp ' nomad_src 'Evaluator.cpp ' nomad_src 'Evaluator_Control.cpp ' nomad_src 'Exception.cpp '...
    nomad_src 'Extended_Poll.cpp ' nomad_src 'GMesh.cpp ' nomad_src 'L_Curve.cpp ' nomad_src 'LH_Search.cpp ' nomad_src 'OrthogonalMesh.cpp ' nomad_src 'Mads.cpp ' nomad_src 'Model_Sorted_Point.cpp '...
    nomad_src 'Model_Stats.cpp ' nomad_src 'Multi_Obj_Evaluator.cpp ' nomad_src 'Parameter_Entries.cpp '...
    nomad_src 'Parameter_Entry.cpp ' nomad_src 'Pareto_Front.cpp ' nomad_src 'Pareto_Point.cpp ' nomad_src 'Phase_One_Evaluator.cpp '...
    nomad_src 'Phase_One_Search.cpp ' nomad_src 'Point.cpp ' nomad_src 'Priority_Eval_Point.cpp ' nomad_src 'Quad_Model.cpp '...
    nomad_src 'Sgtelib_Model_Evaluator.cpp ' nomad_src 'Sgtelib_Model_Search.cpp ' nomad_src 'Sgtelib_Model_Manager.cpp ' ...
    nomad_src 'Quad_Model_Evaluator.cpp ' nomad_src 'Quad_Model_Search.cpp ' nomad_src 'Random_Pickup.cpp ' nomad_src 'RNG.cpp '...
    nomad_src 'Signature.cpp ' nomad_src 'Slave.cpp ' nomad_src 'SMesh.cpp ' nomad_src 'Speculative_Search.cpp ' nomad_src 'Stats.cpp ' nomad_src 'utils.cpp '...
    nomad_src 'Variable_Group.cpp ' nomad_src 'VNS_Search.cpp ' nomad_src 'XMesh.cpp ' ...
    nomad_src_sgtelib 'Kernel.cpp ' nomad_src_sgtelib 'Surrogate_Ensemble.cpp ' nomad_src_sgtelib 'Surrogate_LOWESS.cpp	'...
    nomad_src_sgtelib 'Surrogate_Parameters.cpp	' nomad_src_sgtelib 'TrainingSet.cpp ' nomad_src_sgtelib 'Matrix.cpp '...
    nomad_src_sgtelib 'Surrogate_Factory.cpp ' nomad_src_sgtelib 'Surrogate_PRS.cpp ' nomad_src_sgtelib 'Surrogate_RBF.cpp '...
    nomad_src_sgtelib 'sgtelib.cpp ' nomad_src_sgtelib 'Surrogate.cpp ' nomad_src_sgtelib 'Surrogate_KS.cpp ' nomad_src_sgtelib 'Surrogate_PRS_CAT.cpp '...
    nomad_src_sgtelib 'Surrogate_Utils.cpp ' nomad_src_sgtelib 'sgtelib_help.cpp '  nomad_src_sgtelib 'Surrogate_CN.cpp ' ...
    nomad_src_sgtelib 'Surrogate_Kriging.cpp '  nomad_src_sgtelib 'Surrogate_PRS_EDGE.cpp ' nomad_src_sgtelib 'Tests.cpp ' ];

else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LINUX AND OSX  ---> use dynamic libraries
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Default library names
    nameLibNomad = 'libnomad.so';

    % Default update LDFLAGS (linux only)
    updateLDFLAGS= '';
    % Post compilation tool for path to library (osx only)
    install_name_tool='';

    nomad_home = getenv('NOMAD_HOME');

    if ( length(nomad_home) < 1 )
        % Get a default directory for NOMAD_HOME
        cd ..
        cd .. 
        cd ..
        nomad_home = cd; 
        if ( ~ exist(nomad_home,'dir') )
            error('The default NOMAD_HOME directory does not exist. Please provide a correct value for the NOMAD_HOME variables with the command setenv(''NOMAD_HOME'',ARG1)');
        end
        warning('The NOMAD_HOME variable for Matlab is set with its default value %s. The default can be replaced by using the command setenv(''NOMAD_HOME'',ARG1)! before running the GERAD_NOMAD_build command.',nomad_home);
        cd(cdir);
    else
        if ( ~isempty( find(isspace(nomad_home),1) ) )
            error('The compilation of Nomad for Matlab uses the sources located in the NOMAD_HOME directory. The NOMAD_HOME directory should not contain empty space. Please consider moving the NOMAD_HOME directory and reset the NOMAD_HOME environment variable accordingly.');
        end
    end
      
    nomad_src=[nomad_home filesep 'src' filesep];
    sgtelib_src=[nomad_home filesep 'ext' filesep 'sgtelib' filesep 'src'];
    nomad_lib=[nomad_home filesep 'lib' filesep];
    nomad_bin=[nomad_home filesep 'bin' filesep]; 
    
    switch(computer)
        case 'GLNX86'
            updateLDFLAGS = 'LDFLAGS=''$LDFLAGS -Wl,-rpath,''''$ORIGIN/../lib/'''' -Wl,-rpath-link,''''../lib/'''' '' ';
        case 'GLNXA64'
            updateLDFLAGS = 'LDFLAGS=''$LDFLAGS -Wl,-rpath,''''$ORIGIN/../lib/'''' -Wl,-rpath-link,''''../lib/'''' '' ';
        case 'MACI64'
            install_name_tool=['install_name_tool -change ' nameLibNomad ' @loader_path/../lib/' nameLibNomad ' ' nomad_bin filesep 'nomad.' mexext];
    end
   
    %Compile & Move (Default) --> use shared object library
    post = [' -I.  -I' nomad_src ' -I' sgtelib_src ' -lut -lnomad -L' nomad_lib ' -output ' nomad_bin filesep 'nomad.' mexext ];
    pre =[ 'mex -v -largeArrayDims nomadmex.cpp ' updateLDFLAGS ];
    
    if ( ~ exist([nomad_lib filesep nameLibNomad],'file') )
        error('The Nomad library file %s is not available. Please perform Nomad project compilation before proceeding.',nameLibNomad);      
    end
    
end
    

fprintf('\n------------------------------------------------\n');
fprintf('NOMAD MEX FILE BUILD --- GERAD VERSION \n\n');

%CD to Source Directory
cd 'Source';

try

    if ( ~ exist([nomad_lib filesep nameLibNomad],'file') )
        error('The Nomad library file %s is not available. Please perform Nomad project compilation before proceeding.',nameLibNomad);      
    end

    eval([pre post])
    
    if ( strcmp(computer,'MACI64') == 1 )
        system(install_name_tool);
    end
    
    cd(cdir);
    fprintf('Compilation done!\n');
    fprintf('\n----------------------------------------------------------------------------------------------\n');
    fprintf(' To be able to use the nomad functions, you may need to modify the Matlab path \n');
    qstring = 'To be able to use the nomad functions, you may need to modify the Matlab path. Do you want to update the Matlab path?';
    choice = questdlg(qstring,'Set path','Yes','No','Yes');
    if ( strcmp(choice,'Yes') )
        addpath([ cdir filesep 'Functions']);
        addpath(nomad_bin);
        fprintf('  ---> The Matlab path has been modified but not saved.\n');
    end
    clear nomad_home nomad_lib nomad_bin nomad_src cdir post pre updateLDFLAGS qstring choice install_name_tool nameLibNomad;
catch ME
    cd(cdir);
	clear nomad_home nomad_lib nomad_bin nomad_src cdir post pre updateLDFLAGS qstring choice install_name_tool nameLibNomad;
    error('Error Compiling NOMAD!\n%s',ME.message);
end
