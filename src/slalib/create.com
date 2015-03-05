$!
$!  - - - - - - -
$!   C R E A T E
$!  - - - - - - -
$!
$!  Create SLALIB releases from source - VAX and Unix
$!
$!  The command is @CREATE
$!
$!  The default directory must be the one containing the source
$!  modules.  The VAX release will be inserted into [.RELEASE]
$!  and the Unix release into [.UNIX].  Prior contents of any
$!  pre-existing [.RELEASE] and [.UNIX] directories will be lost.
$!
$!  P.T.Wallace   Starlink   22 January 1998
$!
$!----------------------------------------------------------------------
$!
$!  VMS
$!
$!  Create an empty [.RELEASE] directory
$     IF F$SEARCH("RELEASE.DIR").EQS."" THEN $CREATE/DIR [.RELEASE]
$     IF F$SEARCH("[.RELEASE]*.*").NES."" THEN $DELETE [.RELEASE]*.*;*
$!
$!  Copy the document and news item
$     COPY SUN67.TEX [.RELEASE]*.*
$     COPY SLA.NEWS [.RELEASE]*.*
$!
$!  Initialise the text and object libraries
$     LIBR/CREATE/TEXT [.RELEASE]SLALIB.TLB
$     PURGE [.RELEASE]SLALIB.TLB
$     LIBR/CREATE [.RELEASE]SLALIB.OLB
$     PURGE [.RELEASE]SLALIB.OLB
$!
$!  Update the libraries
$UPOBJ:
$      FILE = F$SEARCH("*.FOR")
$      IF FILE .EQS. "" THEN GOTO UPOBJX
$      NAME = F$PARSE(FILE,,,"NAME")
$      @PUT 'NAME'
$     GOTO UPOBJ
$UPOBJX:
$     @PUT GRESID.VAX
$     @PUT RANDOM.VAX
$     @PUT WAIT.VAX
$!
$!  Compress
$     LIBR/COMP/DATA=REDUCE/TEXT/OUTPUT=[.RELEASE]SLALIB.TLB -
                                        [.RELEASE]SLALIB.TLB
$     PURGE [.RELEASE]SLALIB.TLB
$     LIBR/COMP/OUTPUT=[.RELEASE]SLALIB.OLB [.RELEASE]SLALIB.OLB
$     PURGE [.RELEASE]SLALIB.OLB
$!
$!----------------------------------------------------------------------
$!
$!  UNIX
$!
$!  Create an empty [.UNIX] directory
$     IF F$SEARCH("UNIX.DIR").EQS."" THEN $CREATE/DIR [.UNIX]
$     IF F$SEARCH("[.UNIX]*.*").NES."" THEN $DELETE [.UNIX]*.*;*
$!
$!  Copy the platform-independent Fortran source
$FLOOP:
$     FILE = F$SEARCH("*.FOR")
$     IF FILE.EQS."" THEN $GOTO FLOOPX
$     FILE=F$EXTRACT(0,F$LOCATE(";",FILE),FILE)
$     NAME = F$PARSE(FILE,,,"NAME")
$     IF NAME .EQS. "GRESID" THEN $GOTO FLOOP
$     IF NAME .EQS. "RANDOM" THEN $GOTO FLOOP
$     IF NAME .EQS. "WAIT" THEN $GOTO FLOOP
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT 'FILE' [.UNIX]'NAME'.F
$     GOTO FLOOP
$FLOOPX:
$!
$!  Copy the platform specific source code
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT GRESID.VAX [.UNIX]GRESID.VAX
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT GRESID.PCM [.UNIX]GRESID.PCM
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT GRESID.CNVX [.UNIX]GRESID.CNVX
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT GRESID.MIPS [.UNIX]GRESID.MIPS
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT GRESID.SUN4 [.UNIX]GRESID.SUN4
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT GRESID.LNX [.UNIX]GRESID.LNX
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT GRESID.DEC [.UNIX]GRESID.DEC
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT RANDOM.VAX [.UNIX]RANDOM.VAX
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT RANDOM.PCM [.UNIX]RANDOM.PCM
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT RANDOM.CNVX [.UNIX]RANDOM.CNVX
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT RANDOM.MIPS [.UNIX]RANDOM.MIPS
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT RANDOM.SUN4 [.UNIX]RANDOM.SUN4
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT RANDOM.LNX [.UNIX]RANDOM.LNX
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT RANDOM.DEC [.UNIX]RANDOM.DEC
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT RTL_RANDOM.C [.UNIX]RTL_RANDOM.C
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT WAIT.VAX [.UNIX]WAIT.VAX
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT WAIT.PCM [.UNIX]WAIT.PCM
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT WAIT.CNVX [.UNIX]WAIT.CNVX
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT WAIT.MIPS [.UNIX]WAIT.MIPS
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT WAIT.SUN4 [.UNIX]WAIT.SUN4
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT WAIT.LNX [.UNIX]WAIT.LNX
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT WAIT.DEC [.UNIX]WAIT.DEC
$!
$!  Copy the miscellaneous files
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT MK [.UNIX]MK
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT MAKEFILE [.UNIX]MAKEFILE
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT SLA_LINK [.UNIX]SLA_LINK
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT SLA_LINK_ADAM [.UNIX]SLA_LINK_ADAM
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT SUN67.TEX [.UNIX]SUN67.TEX
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT READ.ME [.UNIX]READ.ME
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT SLA.NEWS [.UNIX]SLA.NEWS
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT PC.BAT [.UNIX]PC.BAT
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT REP.BAT [.UNIX]REP.BAT
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT CREATE.COM [.UNIX]CREATE.COM
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT PUT.COM [.UNIX]PUT.COM
$     CONVERT/FDL=SYS$SYSTEM:UCX$CONVERT VAX_TO_UNIX.USH [.UNIX]VAX_TO_UNIX.
$!                                                   ^^^
$!                                    The USH suffix prevents accidental
$!                                    execution of the script when the
$!                                    working directory is .RELEASE instead
$!                                    of .UNIX
$!
$!  Explain what to do next
$     COPY SYS$INPUT SYS$OUTPUT

   To complete building the Unix release, please do the following:

   1) Login to the Unix machine.

   2) Locate the NFS-served directory corresponding to subdirectory
      [.UNIX] of the current default directory.

   3) For efficiency, and to avoid possible problems involving case
      sensitivity in filenames, copy all the files in that directory
      to a scratch directory on the Unix machine.

   4) Type "vax_to_unix" to archive all the source files.

$!
$!----------------------------------------------------------------------
$!
$!  Wrap up
$     PURGE [...]
$     SET FILE/TRUNCATE [...]*.*/EXCLUDE=CREATE.COM
$     EXIT
