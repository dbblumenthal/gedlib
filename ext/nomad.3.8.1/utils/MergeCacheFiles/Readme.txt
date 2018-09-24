Cache files from different runs of the same problem can be merged and reused in a later optimization.

Please note that binary cache files must be from the same OS. Please note that the provided cacheX.bin files have been created under OSX.

1- Using the OS (linux and OSX) command line 

- The following sources codes allow to merge two binary files (cache1.bin and cache2.bin) into a new cache file called cache3.bin.

- On the command line execute the make command

- Execute  ./merge.exe 

2- Using Matlab environment

- If not already done, setup the mex Matlab External Interface with the command mex -setup at the command prompt. Please refer to Matlab documentation for details.

- Change the Current Folder to %NOMAD_HOME%/utils/MergeCacheFiles

- Run the command build.m

- Run the command mergeM(‘cache1.bin’,’cache2.bin’,’cache3.bin’)
