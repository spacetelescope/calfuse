@echo off
rem
rem  - - - - - - - -
rem   R E P . B A T
rem  - - - - - - - -
rem
rem  Update one module in the SLALIB library
rem
rem  Command:
rem
rem     REP module
rem
rem  File SLALIB.BAK is deleted.
rem
rem  P.T.Wallace   Starlink   5 April 1992
rem
fl/c /FPi %1.for
lib slalib -+%1;
del %1.obj
del slalib.bak
echo:
