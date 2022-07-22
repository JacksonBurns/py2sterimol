====================================

STERIMOL DOCUMENTATION IS AVAILABLE FROM CCL AT:
http://www.ccl.net/cca/documents/STERIMOL_Documentation/

====================================    

Original documentation changed in few places by Jan Labanowski on Feb 28, 2009

                            STERIMOL for the PC
 
               Calculation of Verloop's STERIMOL Parameters



            Based on the program STERIMOL by W. Hoogenstraaten
                       Duphar B.V., the Netherlands

                Converted for the IBM PC by Stephen Bowlus
                       Sandoz Crop Protection Corp.
                           Palo Alto, California


     The enclosed distribution contains source code and executable files 
for STERIMOL, the program developed for estimation of the Verloop steric 
parameters for QSAR.  The original program was written by Verloop and Hoogen-
straaten at Duphar B.V.  This went through several revisions (see comments
in the file STERIMOL.SRC).  This PC version is based on the program obtained
from E. L. Plummer at FMC Ag Chem.

            The FORTRAN code for the PC has been somewhat cleaned up, to
remove many references and explanations of the generations of changes that
were made, except where these were deemed important to understanding the
current version of the program.  If the history of the source code is of
importance or interest, the file STERIMOL.SRC contains the unconverted VAX
version of the program.  Rather than clutter up the code with (more) extra-
neous remarks, notes on the conversion are appended to the documentation.  
 


                        INSTALLATION and OPERATION


            The distribution contains the following files:

      sterimol.zip -- ZIP archive of the STERIMOL package. Unpack with:
                         unzip sterimol.zip
                      or use MS Windoz equivalent.
      sterimol.tar.gz -- Tar ball of STERIMOL package. Unpack with:
                      gunzip sterimol.tar.gz; tar xvf sterimol.tar

      README.txt   -- This file
      convert.doc  -- original documentation
      sterimol.f   -- minute changes to original sterimol.for to make program
                      compile under GNU g77 under Linux
      sterimol.src -- VAX source code, from which this program was written.
      sterimol.for -- Source code for MicroSoft FORTRAN, v. 3.31 or 5.0.
      STERIMOL.EXE -- Executable image compiled with g77 under Fedora Core 6
                      (will most likely not run on any other version of Linux).
      sterimol.exe -- Executable image, compiled under FORTRAN 5.0 for the
                      80286 or 80386 with math coprocessor.
      sterimol.52i -- Executable image, compiled under FORTRAN 5.0 for the
                      80286 or 80386 using coprocessor emulation (coprocessor
                      optional).
      sterimol.50i -- Executable image, compiled under FORTRAN 5.0 for the
                      8086 using coprocessor emulation.
      sterimol.331 -- Executable image, compiled under FORTRAN 3.31, using
                      the default math library (coprocessor not used).
      example.inp  -- Sample input deck
      example.out  -- Sample results file.

To compile program under Linux
==============================
To compile it under Linux, you need to have the g77 installed. If
   which g77
does not show it, you need to install it.
On Fedora you would do the following as root:
   yum list all | grep g77
It will list you some packages. like:
compat-gcc-XX-g77.YY ....
Install all thiese pacjages as:
yum install compat-gcc-XX-g77.YY
Note... It will install many more packages that are needed by g77.

If you have g77, compile it as:  
      g77 -o STERIMOL.EXE sterimol.f
This should produce STERIMOL.EXE executable. You can place this
executable (STERIMOL.EXE) somewhere in your Unix PATH (e.g., in your
$HOME/bin directory or if you have a root access then in the 
/usr/local/bin or /usr/bin). If you do not have the executable in your PATH
environment variable, you need to use the full path to this file to run it,
say, if you compiled the program in /home/johny/sterimol, you will need
to run it as:
     /home/johny/sterimol/STERIMOL.EXE < file.inp > file.out

To compile program under MS Windoz
===============================
You are well advised to install Cygwin on your Windoz.
Cygwin is just a Linux-like environment that works within MS-Windoz.
More at: http://www.cygwin.com/
It is free. The g77 will not install by default. You need to specifically
select it during installation (you will be presented a page to
add additional software, click on [+] by Development and then
click on gcc and g77). For details Google for: "Windows Cygwin g77"

----------------------------------------------------------------------
To install the program on a hard disk, copy the files STERIMOL.EXE
(or other executable file listed above if you have these platforms listed,
renaming if necessary) and EXAMPLE.INP to an appropriate directory.

   The input file, example.inp is edited with any ASCII text editor
for the compound and radicals (fragments) of interest.  The rules for formula
notation, written by Hoogenstraaten, are given in the following section.  For
the remainder of the file, note the following, taken from the program
listing:

C  *****   INPUT  CARDS   *****                                     
C  THE FIRST CARD CONTAINS THE  33  USED SYMBOLS; STARTING WITH 1 IN COLUMN 
C  1.  SYMBOL NR  27 IS A SPACE.                        
C      1234567890ABCDEFHINOPRSTXZ (&*),=                      


      þ  This line no longer appears in the input file.  This data declara-
         tion is now in the body of the program.  This line should be not be
         included in any test file based on Hoogenstraaten's examples in the
         following section of the documentation.


C    NEXT FOLLOW (IF WANTED) THE FORMULAE FOR THE 'RADICALS', OF THE       
C  FOLLOWING GENERAL TYPE:                              
C       ZKK=FORMULA*                                                
C  WHERE  KK  IS A TWO-DIGIT INTEGER NUMBER (01 - 99).  A RADICAL MAY BE ANY
C  SYMBOL STRING WHICH IS PART OF A VALID FORMULA.  IF NECESSARY, IT MAY BE
C  CONTINUED LIKE A FORMULA CARD.  CLOSE OR REPLACE THE RADICAL-CARD SET BY
C  A BLANK CARD.                            
C     FOR EACH MOLECULE THE FOLLOWING CARDS ARE GIVEN:                    
C  (1)  THE INTEGER  IPR  (I4):  THE OUTPUT INDEX (1 OR 2);               
C  (2)  THE FORMULA.  IF ONE CARD IS NOT SUFFICIENT, END THIS CARD WITH A 
C       &  SIGN AND CONTINUE ON THE NEXT CARD.          
C  (3)  NUMBER OF THE TORSION ANGLES (I4).  IF NONE ARE TO BE GIVEN, REPLACE
C       THIS CARD BY A BLANK ONE.                       
C  (4)  THE TORSION ANGLES (1X,F7.2); DECIMAL POINTS AT COLUMNS  6, 14, 22,
C       ....) .                                   
C    AFTER THE SETS OF MOLECULE-INPUT CARDS:  ONE BLANK CARD.             

      þ  Only one molecule may be processed in each input file.  Input files
         may be chained in a DOS .BAT file, if a number of molecules or   
         conformers are being considered.


1) The command under Linux (if you compiled STERIMOL with g77 as described 
   above):
  ./STERIMOL.EXE < sterimol.inp

  Under MS-DOS it would be (assuming that the original sterimol.exe still
   work on your version of MS Windoz:
  C:\> sterimol.exe < sterimol.inp
 
or, correspondingly

2) ./STERIMOL.EXE < sterimol.inp > STERIMOL.OUT   # Linux
    C:\> sterimol < sterimol.inp > STERIMOL.OUT 
 
In the first command (1), all output will be displayed on the screen. 
With the second command (2), all output will be written to the file:
STERIMOL.OUT, without screen echo.  Since the screen scrolls very fast
(too fast to read much of the output), the second command is preferred.
Alternatively, the output can be directed to an attached printer
(> PRN: under MS-DOS,).  
 
For the FORTRAN v. 3.31 version of the program, rather than use
DOS redirection, the program will prompt you for Unit 5 and Unit 6.  File
names for input and output, respectively, should be given at the prompts.

The input and output files may have any legal DOS name (or Linux name).
The extension is not required.  Files can be named for the molecule being
examined, for example, COMPND1.INP and COMPND1.OUT, or JUNKIN and JUNKOUT.


C  *****   OUTPUT   *****                                           
C       FOR  IPR = 1 :                                        
C  THE RADICAL FORMULAE (IF PRESENT);                               
C  THE MOLECULE FORMULA AS IT WAS INPUT;                      
C  SERIAL NUMBER, NUMBER OF ATOMS, OUTPUT INDEX;                    
C  FULLY WRITTEN FORMULA, WITH 'TRANSLATION' OF THE RADICALS;             
C  LIST OF THE TORSION ANGLES;                                      
C  TABLE OF THE VANDERWAALS RADII AND COORDINATES OF THE ATOMS;           
C  LIST OF THE STERIC PARAMETERS: L; B(1) - B(4) AND B(5) (MAXIMUM WIDTH);
C  A WARNING IF NEIGHBOURING ATOMS HAVE OVERLAPPING VANDERWAALS SPHERES.
C       FOR  IPR = 2  IN ADDITION:                            
C  AN  X, Y  AND A  Y, Z  PROJECTION OF THE MOLECULE;               
C  A LIST OF THE ATOMS, ORDERED BY INCREASING  X  COORDINATES, SERVING TO  
C  RAPIDLY IDENTIFY THE PLOTTED ATOMS.                  

      
      þ  The output type parameter (single digit above the structure nota- 
         tion) is set to 2 for a "connect-the-dots" structure representation,
         or to 1 for normal, numerical output.  The printed structures may
         be meaningless, due to line wrapping.  Success is more likely with
         a printer set to 132 columns.


 
                             CONVERSION NOTES
                             ================

      If you are really interested in this type of bumpf, the following des-
cription will probably make more sense if reviewed together with listings of
STERIMOL.SRC (the original VAX code) and STERIMOL.FOR (the Microsoft FORTRAN
code).  While some of the corrections are serious, others were made only to
stop complaints from the compiler:


             The BLOCK DATA segment was altered to correct a statement out
of order.  The COMMON/CHAR/ declaration was moved to the top of the block.
The COMMON block CHAR was redefined throughout the program to ICH(33).  The
DATA statement in which the characters are defined was altered so the block
definition and declared data are in agreement.

            The COMMON block VAR was made into a blank common block, to
accommodate different definitions of this block.
 
            The variables P and KR were removed where they were referenced
but not used.  (N.B.  There is a surviving constant P!)
 
            I/O unit assignments were changed.  Unit 5 is the standard input
(keyboard), redirected to the input file in the DOS command line when the
program is started.  Unit 6 is the standard output (screen) and the output
file is created by DOS redirection.  (Note that MS Fortran 3.31 and 5.0 han-
dle these differently!)  References to other I/O units in the program (mostly
to get the VAX to write to the terminal), and unit assignments to input and
output files were deleted. 

            All references or calls to EXIT were removed.  Statement 9999 was
changed from "CALL EXIT" to "CONTINUE".  Conditional branches to CALL EXIT
are now directed to this continuation.  

            The DO loops associated with statement 130 were renumbered (now
131) and the inner loop now terminating with 132 was revised, to stop com-
plaints of using a label across blocks.


