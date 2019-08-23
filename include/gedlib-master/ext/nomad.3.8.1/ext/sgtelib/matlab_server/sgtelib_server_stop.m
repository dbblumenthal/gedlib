function sgtelib_server_stop
disp('Kill sgtelib.exe');
!touch flag_quit
%!killName sgtelib.exe
