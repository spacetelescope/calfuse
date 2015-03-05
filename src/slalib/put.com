$!
$!  - - - -
$!   P U T
$!  - - - -
$!
$!  Update one SLALIB routine from Fortran source
$!
$!  DCL command is @PUT file
$!
$!  The default directory must be the one containing the
$!  Fortran source, with the libraries in [.RELEASE].
$!
$!  P T Wallace   Starlink   22 January 1993
$!
$!  Save supplied file name and strip recognized extensions
$     FILE=''P1''
$     P1=P1-".FOR"-".VAX"
$!
$!  No action required for TEST program
      IF P1.EQS."TEST" THEN $ GOTO DONE
$!
$!  If platform-specific module, make .FOR file ...
$     IF F$SEARCH("''P1'.VAX").NES."" THEN $ COPY 'P1'.VAX *.FOR
$!
$!  Update the source library
$     LIBR/REPL/TEXT [.RELEASE]SLALIB.TLB 'P1'.FOR
$!
$!  Compile, update object library, delete object
$     FORTRAN/NOLIST 'P1'.FOR
$     LIBR/REPL [.RELEASE]SLALIB.OLB 'P1'.OBJ
$     DELETE 'P1'.OBJ;*
$!
$!  If module just updated was platform-specific, delete the .FOR version
$     IF F$SEARCH("''P1'.VAX",1).NES."" THEN $ DELETE 'P1'.FOR;*
$!
$!  Finished
$DONE:
$     EXIT
