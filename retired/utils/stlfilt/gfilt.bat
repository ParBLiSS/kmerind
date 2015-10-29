@echo off
rem gfilt.bat: Compile with g++ (2.95.x or 3.x) using STL Error Decryption

rem ********************************************************************
rem Configure COMPILER, COPTS, STLFILT and DOPTS as appropriate for your system
rem (Note: If the env. variable GFILTOPTS is defined, it is appended to DOPTS):

rem Name of your g++ compiler command:
set COMPILER=c++2

rem Compiler options:
set COPTS=-fmessage-length=0

rem Full pathname of the Perl script:
set STLFILT=d:\src\cl\gcc\gSTLFilt.pl

rem Decryptor options (to the Perl script):
set DOPTS=/hdr:m /cand:s /iter:L %GFILTOPTS%

rem ********************************************************************

if not "%1" == "" goto compile

echo Usage: gfilt [options] source-file
echo Compiles source-file, filters errors through STL Error Decryptor.
echo.
goto done

:compile
rem
rem ********************************************************************
rem Note: you can change the /hdr and /cand options for varying detail.
rem See the comments in gSTLFilt.pl for an explanation of the options.
rem ********************************************************************

%COMPILER% %COPTS% %*  2>&1 | perl -W %STLFILT% %DOPTS%

rem Alternative compile line, with pagination of filtered output:
rem %COMPILER% %COPTS% %*  2>&1 | perl %STLFILT% %DOPTS% | more

:done
