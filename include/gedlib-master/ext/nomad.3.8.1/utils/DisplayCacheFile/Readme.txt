A cache file from a different run of the same problem can be reused in a later optimization or to continue an optimization.
The code can  display a binary cache file whose name is provided as first argument.

Please note that binary cache files must be from the same OS. The provided cache.bin file has been created under OSX.

1- Using the OS (linux and OSX) command line 

- On the command line execute the make command

- Execute  ./displayCache.exe cache.bin 

2- Using Matlab environment

- If not already done, setup the mex Matlab External Interface with the command mex -setup at the command prompt. Please refer to Matlab documentation for details.

- Change the Current Folder to %NOMAD_HOME%/utils/DisplayCacheFile

- Run the command build.m

- Run the command displayM(‘cache.bin’)
