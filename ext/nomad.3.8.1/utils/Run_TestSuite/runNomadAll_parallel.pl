#!/usr/bin/perl
# use strict;
use warnings;
 
use Data::Dumper;
 
use threads;
use threads::shared;
use Thread::Semaphore;
 
 
if (! exists $ENV{NOMAD_HOME})  {
    print "Please provide an environment variable NOMAD_HOME. \n";
   exit 1;
}

print "NOMAD_HOME = $ENV{NOMAD_HOME} \n";

#if (! exists $ENV{LIBAMPLDIR}) {
#	print "!!!!!!!!!!!!!!!!   The LIBAMPLDIR environment variable is not set -> ampl example may not work! \n";
#} 

#if (! exists $ENV{CUTER} | ! exists $ENV{MYCUTER} | ! exists $ENV{SIFDEC} | ! exists $ENV{MYSIFDEC} ) {
#	print "!!!!!!!!!!!!!!!! The environment variables for CUTEr are not set -> CUTEr example may not work! \n";
#}


my $MPIcommand:shared;
open(aComRes,"uname |") or die "Can't run command: $!\n";
if ( grep(/MINGW/,<aComRes>) ){
	$MPIcommand="mpiexec  -localonly";    # version mingw
} else
{
	$MPIcommand="mpirun -n";    # version OSX and linux
}

my $keySearch=" Erreur | Error | error | error: | Erreur | erreur | Exception | NOMAD::Exception | Failed | Arrêt | Stop";

my $nombre_de_jobs_en_parallele:shared;
if ( ! exists $ARGV[0]) {
$nombre_de_jobs_en_parallele = 1;    # 1 but still 3 processes minimun required for mpi 
} else {
$nombre_de_jobs_en_parallele = $ARGV[0];
}


my $nomadEXE="nomad";
my $nomadMPI_EXE="nomad.MPI";

my $maxLinesForLogProblem:shared;
$maxLinesForLogProblem = 20;
 
my @list = (
	["cd $ENV{NOMAD_HOME}/examples/basic/batch/single_obj ; if [ ! -e bb.exe  ] ; then g++ -o bb.exe  bb.cpp 2>&1 ; fi  ; $ENV{NOMAD_HOME}/bin/$nomadEXE param.txt 2>&1"],
	["cd $ENV{NOMAD_HOME}/examples/basic/batch/single_obj_parallel ; if [ ! -e bb.exe  ] ; then g++ -o bb.exe  bb.cpp 2>&1 ; fi  ; $ENV{NOMAD_HOME}/bin/$nomadEXE param.txt 2>&1"],
	["cd $ENV{NOMAD_HOME}/examples/basic/batch/bi_obj ; if [ ! -e bb.exe  ] ; then g++ -o bb.exe  bb.cpp 2>&1 ; fi ; $ENV{NOMAD_HOME}/bin/$nomadEXE param.txt 2>&1"],
	["cd $ENV{NOMAD_HOME}/examples/basic/batch/single_obj  ; if [ ! -e bb.exe  ] ; then g++ -o bb.exe  bb.cpp 2>&1 ; fi ; $MPIcommand 3 $ENV{NOMAD_HOME}/bin/$nomadMPI_EXE param.txt 2>&1"],
	["cd $ENV{NOMAD_HOME}/examples/basic/batch/bi_obj ; if [ ! -e bb.exe  ] ; then g++ -o bb.exe  bb.cpp 2>&1 ; fi  ; $MPIcommand 3 $ENV{NOMAD_HOME}/bin/$nomadMPI_EXE param.txt 2>&1"] ,
	["cd $ENV{NOMAD_HOME}/examples/basic/library/single_obj ; if [ -e basic_lib.exe  ] ; then rm -f basic_lib.exe  ; fi  ; make clean 2>&1 ; make 2>&1; ./basic_lib.exe  2>&1"],
	["cd $ENV{NOMAD_HOME}/examples/basic/library/bi_obj ; if [ -e basic_lib.exe  ] ; then rm -f basic_lib.exe  ; fi  ; make clean 2>&1 ; make 2>&1 ; ./basic_lib.exe  2>&1"],
	["cd $ENV{NOMAD_HOME}/examples/basic/library/single_obj ; if [ -e basic_lib_MPI.exe  ] ; then rm -f basic_lib_MPI.exe   ; fi  ; make clean 2>&1 ; make mpi 2>&1; $MPIcommand 3 ./basic_lib_MPI.exe  2>&1"],
	["cd $ENV{NOMAD_HOME}/examples/basic/library/bi_obj ; if [ -e basic_lib_MPI.exe  ] ; then rm -f basic_lib_MPI.exe  ; fi  ; make clean 2>&1 ; make mpi 2>&1; $MPIcommand 3 ./basic_lib_MPI.exe  2>&1"],	
	["cd $ENV{NOMAD_HOME}/examples/basic/library/single_obj_parallel ; if [ -e basic_lib.exe  ] ; then rm -f basic_lib.exe  ; fi  ; make clean 2>&1 ; make 2>&1; ./basic_lib.exe  2>&1"],
    ["cd $ENV{NOMAD_HOME}/examples/advanced/categorical/batch ; if [ ! -e bb.exe  -o ! -e neighbors.exe  ] ; then g++ -o bb.exe  bb.cpp 2>&1 ; g++ -o neighbors.exe  neighbors.cpp ; fi ; echo ; $ENV{NOMAD_HOME}/bin/$nomadEXE param.txt 2>&1"],
    ["cd $ENV{NOMAD_HOME}/tools/PSD-MADS ; if [ -e psdmads.exe  ] ; then rm -f *.exe  *.o ; fi ; make 2>&1 ; cd problems/G2_20 ; g++ -o bb.exe  bb.cpp 2>&1; $MPIcommand 6 ../../psdmads.exe  param.txt 50 5 2>&1"],
	["cd $ENV{NOMAD_HOME}/examples/advanced/categorical/batch ; if [ ! -e bb.exe  -o ! -e neighbors.exe  ] ; then g++ -o bb.exe  bb.cpp 2>&1 ; g++ -o neighbors.exe  neighbors.cpp ; fi ; echo ; $MPIcommand 3 $ENV{NOMAD_HOME}/bin/$nomadMPI_EXE param.txt 2>&1"],
	["cd $ENV{NOMAD_HOME}/examples/advanced/categorical/single_obj ; if [ -e categorical.exe  ] ; then rm -f categorical.exe  ; echo ; fi ; make clean 2>&1 ; make 2>&1; ./categorical.exe  2>&1 "],
	["cd $ENV{NOMAD_HOME}/examples/advanced/categorical/single_obj ; if [ -e categorical_MPI.exe  ] ; then rm -f categorical_MPI.exe  ; echo ; fi ; make clean 2>&1 ; make mpi 2>&1; $MPIcommand 3 ./categorical_MPI.exe  2>&1 "],
	["cd $ENV{NOMAD_HOME}/examples/advanced/categorical/bi_obj ; if [ -e categorical.exe  ] ; then rm -f categorical.exe  ; echo ; fi ; make clean 2>&1 ; make 2>&1; ./categorical.exe  2>&1 "],
	["cd $ENV{NOMAD_HOME}/examples/advanced/categorical/bi_obj ; if [ -e categorical_MPI.exe  ] ; then rm -f categorical_MPI.exe  ; echo ; fi ; make clean 2>&1 ; make mpi 2>&1; $MPIcommand 3 ./categorical_MPI.exe  2>&1 "],
	["cd $ENV{NOMAD_HOME}/examples/advanced/restart ; if [ -e restart.exe  ] ; then rm -f restart.exe  ; echo ; fi ; make clean 2>&1 ; make 2>&1; ./restart.exe  2>&1 "],
	["cd $ENV{NOMAD_HOME}/examples/advanced/restart ; if [ -e restart_MPI.exe  ] ; then rm -f restart_MPI.exe  ; echo ; fi ; make clean 2>&1 ; make mpi 2>&1; $MPIcommand 3 ./restart_MPI.exe  2>&1 "],
	["cd $ENV{NOMAD_HOME}/examples/advanced/multi_start ; if [ -e multi.exe  ] ; then rm -f multi.exe  ; echo ; fi ; make clean 2>&1 ; make 2>&1; ./multi.exe  param.txt 4 5 2>&1 "],
	["cd $ENV{NOMAD_HOME}/examples/advanced/user_search ; if [ -e user_search.exe  ] ; then rm -f user_search.exe  ; echo ; fi ; make clean 2>&1 ; make 2>&1; ./user_search.exe  2>&1 "],
	["cd $ENV{NOMAD_HOME}/examples/advanced/user_search ; if [ -e user_search_MPI.exe  ] ; then rm -f user_search_MPI.exe  ; echo ; fi ; make clean 2>&1 ; make mpi 2>&1; $MPIcommand 3 ./user_search_MPI.exe  2>&1 "],
# 	["cd $ENV{NOMAD_HOME}/examples/interfaces/AMPL ; if [ ! -e bb.exe  ] ; then make 2>&1 ; echo ; fi  ; $ENV{NOMAD_HOME}/bin/$nomadEXE param.txt 2>&1 "],
# 	["cd $ENV{NOMAD_HOME}/examples/interfaces/AMPL ; if [ ! -e bb.exe  ] ; then make 2>&1 ; echo ; fi  ; $MPIcommand 3 $ENV{NOMAD_HOME}/bin/$nomadMPI_EXE param.txt 2>&1 "],
# 	["cd $ENV{NOMAD_HOME}/examples/interfaces/CUTEr ; if [ ! -e bb  ] ; then $ENV{MYSIFDEC}/bin/sifdecode PROBLEM.SIF ; ./compile 2>&1 ; fi  ; $ENV{NOMAD_HOME}/bin/$nomadEXE parameters.txt 2>&1 "],
# 	["cd $ENV{NOMAD_HOME}/examples/interfaces/CUTEr ; sleep 10 ; if [ ! -e bb  ] ; then ./compile 2>&1 ; fi ; $MPIcommand 3 $ENV{NOMAD_HOME}/bin/$nomadMPI_EXE parameters.txt 2>&1 "],
# 	["cd $ENV{NOMAD_HOME}/examples/interfaces/DLL/single_obj ; if [ -e nomad_for_dll  ] ; then rm -f nomad_for_dll  ; fi ; echo ; make 2>&1; ./nomad_for_dll  parameters.txt 2>&1"],
# 	["cd $ENV{NOMAD_HOME}/examples/interfaces/DLL/bi_obj ; if [ -e nomad_for_dll  ] ; then rm -f nomad_for_dll  ; fi ; echo ; make clean 2>&1 ; make 2>&1; ./nomad_for_dll  parameters.txt 2>&1"],
#   ["cd $ENV{NOMAD_HOME}/examples/interfaces/FORTRAN/example1 ; if [ -e test.exe  ] ; then rm -f test.exe  ; fi ; make clean 2>&1 ; make 2>&1; ./test.exe  2>&1"],
# 	["cd $ENV{NOMAD_HOME}/examples/interfaces/FORTRAN/example2 ; if [ -e test.exe  ] ; then rm -f test.exe  ; fi ; make clean 2>&1 ; make 2>&1; ./test.exe  2>&1"],
	["cd $ENV{NOMAD_HOME}/tools/COOP-MADS ; if [ -e coopmads.exe  ] ; then rm -f *.exe  *.o ; fi ; make 2>&1 ; cd problems/G2_10 ; g++ -o bb.exe  bb.cpp 2>&1; $MPIcommand 3 ../../coopmads.exe  param.txt 2>&1"],
	["cd $ENV{NOMAD_HOME}/tools/COOP-MADS ; sleep 10 ; if [ ! -e coopmads.exe  ] ; then  make 2>&1 ; fi ; cd problems/RHEOLOGY ; g++ -o bb.exe  bb.cpp 2>&1; $MPIcommand 3 ../../coopmads.exe  param.txt 2>&1"],
	["cd $ENV{NOMAD_HOME}/tools/PSD-MADS ; if [ -e psdmads.exe  ] ; then rm -f *.exe  *.o ; fi ; make 2>&1 ; cd problems/G2_20 ; g++ -o bb.exe  bb.cpp 2>&1; $MPIcommand 6 ../../psdmads.exe  param.txt 50 5 2>&1"]
	); 
 
my @NOMADcompilations = (
    ["cd $ENV{NOMAD_HOME}/src ; make del "],
    ["cd $ENV{NOMAD_HOME}/src ; make all -j $nombre_de_jobs_en_parallele 2>&1"],
    ["cd $ENV{NOMAD_HOME}/src ; make mpi -j $nombre_de_jobs_en_parallele 2>&1"]
        );

 
# le sémaphore permet de réguler le nombre de jobs en simultané, bloquant le processus principal tant qu'il n'y a pas de place de libre
my $semaphoreProblems = Thread::Semaphore->new($nombre_de_jobs_en_parallele);

 
# quelques variables partagées
my $cpt:shared;
my $failed:shared;

# Initialisation de compteur et de tests
$cpt = 0;
$failed = 0;
 
my $started = 0;
 
################################################################ 
# Sub pour résoudre les problèmes en parallèle
################################################################
sub RunProblem($$$$$){
	my $aRef = shift;
	my $index = shift;
	my $cpt_ref = shift;
	my $failed_ref = shift;
	my $sema_ref = shift;
 
    open(WRITE_LOG,"> log$index.txt"); 	
 	for my $command (@{$aRef}){

 	    my $fail=1;
 	    
	 	my @Problem=split(/;/,$command);
	 	my $nmax=$#Problem;
	 	my @Path=split(/\s+/,$Problem[0]);

 	    print WRITE_LOG "Path to problem: $Path[1] \n Command: $Problem[$nmax] ; Managed as process $index \n";
 	    open(LOG,"$command |") or die "Can't run program: $!\n";
		if ($? != 0) {
	        print "Failed execution: command $! du process $index \n";
			$$failed_ref++;
			$fail=0;
    	}
    	else {
 		   	@lines = <LOG>;
 		   	foreach $line (@lines){
				if (my @foo = grep(/$keySearch/, $line) ) {
					print "????????? Problem executing command of process $index:\n     ----->   @foo\n";
 	      			print WRITE_LOG "     !!!!!!!!!!!!!!!   @foo \n ";
 	      			$$failed_ref++;
					$fail=0;
					last;
 	      		}
 	      	}
 	      	my $nL=$maxLinesForLogProblem;
 	      	if ( $#lines < $maxLinesForLogProblem){
 	      	    $nL=$#lines;
 	      	}
	 	    for ($i = $#lines - $nL ; $i != $#lines+1; $i++) {
    				print WRITE_LOG "    ----->   $lines[$i] ";
			}	
 	    }
 	    if ($fail==1){
 	    	print WRITE_LOG "++++++++++++ OK! \n";
 	    }
 	    close (LOG);
 	}
    close (WRITE_LOG);

	# incrémente le nombre de jobs finis
	$$cpt_ref++;
 
	# on a une place de libre. Ne pas oublier de libérer le sémaphore même en cas d'erreur
	$$sema_ref->up();
 
	return;
}

#############################
# Sub pour les compilations
#############################
sub CompileNOMAD($$){
	my $aRef = shift;
	my $index = shift;
	my $failed = 0; 
	
 	open(WRITE_LOG,"> logNOMADCompile$index.txt"); 	
 	for my $command (@{$aRef}){
 	
	 	my @Problem=split(/;/,$command);
	 	my @Path=split(/\s+/,$Problem[0]);

 	    print WRITE_LOG "Path: $Path[1] \n Command: $Problem[1] ; Managed as process $index \n";
 	    open(LOG,"$command |") or die "Can't run program: $!\n";
        my $errSearch="Erreur|Error|error|erreur|commande introuvable|compilation failed";
		if ($? != 0) {
	        print "Failed compilation of NOMAD: command $! process $index \n";
			$failed=1;
    	}
    	else {
			while (<LOG>){
				if (my @foo = grep(/$errSearch/, $_) ) {
					print "??????? Problem encountered when compiling NOMAD in process $index:\n     ----->     @foo\n";
 	      			print WRITE_LOG "     ----->   @foo \n ";
					$failed=1;
 	      			last;   			
 	      		}
 	      		if (my @foo = grep(/warning/, $_) ) {
 	      			print WRITE_LOG "     ----->   @foo \n ";
 	      		}
 	      	}
 	    }
 	    if ($failed==0){
 	    	print WRITE_LOG "++++++++++ OK! \n";
 	    }
 	    close (LOG);
 	}
  	close (WRITE_LOG);
	return $failed;
}

#####################################
##### Debut du programme principal
#####################################

# nettoie les fichiers de log
print "########################################################\n";
print "Cleaning old log files\n";
system ("rm -f log*.txt");
print "########################################################\n";

# Fait clean
print "NOMAD clean started";
my $thrCleanNOMAD = threads->create("CompileNOMAD",($NOMADcompilations[0],1));
$failedCleanNOMAD = $thrCleanNOMAD->join();

if ($failedCleanNOMAD!=0){
    print "\n NOMAD clean failed. Stopping here! \n";
    exit 0;
}
print " ---> success. \n";
print "########################################################\n";


# démarre la compilation de nomad
print "NOMAD compilation (not mpi) started";
my $thrNOMAD = threads->create("CompileNOMAD",($NOMADcompilations[1],2));
my $failedCompileNOMAD = $thrNOMAD->join();

if ($failedCompileNOMAD!=0){
    print "\n NOMAD compilation (not mpi) failed. Stopping here! \n";
    exit 0;
}
print " ---> success. \n";
print "########################################################\n";

# démarre la compilation de nomad_mpi
print "NOMAD compilation (mpi) started";
my $thrNOMAD_MPI = threads->create("CompileNOMAD",($NOMADcompilations[2],3));
my $failedCompileNOMAD_MPI = $thrNOMAD_MPI->join();

if ($failedCompileNOMAD_MPI!=0){
	print "\n NOMAD compilation (mpi) failed. Stopping here! \n";
    exit 0;
}
print " ---> success. \n";
print "########################################################\n\n";

print "########################################################\n";
print "Starting parallel executions of problems \n";
print "########################################################\n";


# démarre tous les jobs 
while ($started < scalar @list ){
	my $aRefToAListOfCommands = $list[$started];
	
	# incrémente le compteur
	$started++;
 
	# avons nous une place de libre ?
	$semaphoreProblems->down();
 
	# si le sémaphore est a 0, le processus principal va se bloquer en attendant une nouvelle place
	print "Creating job $started\n";
	my $thr = threads->create("RunProblem", (
		$aRefToAListOfCommands,
		$started,
		\$cpt,
		\$failed,
		\$semaphoreProblems
		)
	);
	# détache le job du thread principal, rend la main au thread principal
	$thr->detach();
 
	# si on veut attendre la fin du job pour redonner la main, on utilise
	# $thr->join();
}
 
print "########################################################\n";
 
# attend les derniers jobs
my $cpt_Prev = $cpt;
print "\n $cpt jobs completed for $started jobs started, patience!\n";
# Disable buffering 
$| = 1;
while ($cpt < $started){
    print ".";
	if ( $cpt > $cpt_Prev) {
		print "\n $cpt jobs completed for $started jobs started, patience!\n";
		}
	$cpt_Prev = $cpt;	
	sleep(3);
}
print "\n $cpt jobs started, $failed jobs failed, ".scalar @list." jobs to be done\n";

print "########################################################\n";

if ($failed !=0) {
	print "----->  Check the readme file for the failed problem(s)!!! \n"; 
}

print "All logs are combined in a single log file\n";
my $i=0;
open (LOGALL, '>', "logAll.txt") or die("Not able to write in logAll.txt, $!");
while ($i < scalar @list ){
	# incrémente le compteur
	$i++;
	open (LOGI, '<', "log$i.txt") or die("Cannot read file log$i.txt, $!\n");	
	print LOGALL "*********************************************************\n";
  	while (my $Ligne = <LOGI> ) {
    	print LOGALL $Ligne;
  	}
  	close(LOGI);
}

print "The End.\n";
