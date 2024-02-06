!===============================================================================
                           MODULE moduleGlobal

!This module defines all the variables globally(i.e for the entire program with
!all the other modules)
!===============================================================================
IMPLICIT NONE

!Defining double precision:
INTEGER,PARAMETER   :: DP=KIND(1.D0)

!Variable used while opening files:
INTEGER             :: error 
CHARACTER(LEN=50)   :: fileName 

!Initializing input file units:
INTEGER, PARAMETER  :: eang = 10, sym = 11, mineang = 12, cnt = 15

!Initializing output file units:
INTEGER, PARAMETER  :: outpt = 13, ggdata = 14, enrgy = 16, &
                       potstereo = 17

!Total number of symmetry matrices in a cube
INTEGER, PARAMETER :: totSymm = 24

!Total elements of each symmetry matrix:
INTEGER, PARAMETER :: totalSymmElements = 9

!Total sets of euler angle in the input file
INTEGER :: totalEulerAngles

!Total sets of Potential Minima euler angle in the input file
INTEGER :: totalMinEulerAngles

!Lattice dimension
INTEGER :: xSide, ySide, zSide
REAL(DP) :: width

!Radius bounding the gaussian minima
REAL(DP) :: radius

!Mobility,Scaling factors
REAL(DP) :: mu , E0  , A, delT
REAL(DP) :: sdNoise

!Depth and SD of the inverse Gaussians
REAL(DP), ALLOCATABLE, DIMENSION(:) :: depth, SD

!Matrix dimension
INTEGER, PARAMETER :: dim = 3

!Symmetry & Lattice options
INTEGER :: symm, potOption

!Value of constants
REAL(DP) :: pi, factor

!===============================================================================
!Public interface:
PUBLIC :: SetInputs             !Called at tata.f90
PUBLIC :: SetConstants          !moduleInput - Subroutine ReadminEulerAngles
PUBLIC :: SetInverseGaussian    
PUBLIC :: FileError
PUBLIC :: Fatal
!===============================================================================
                                 CONTAINS
!===============================================================================


!===============================================================================
                   SUBROUTINE SetInputs
!This subroutine is for giving the initial Inputs
!===============================================================================
IMPLICIT NONE

!Local
CHARACTER(LEN=1) :: ans

WRITE(*,'(/A)')' Enter the sides of the cube : '
WRITE(*,'(/A)',ADVANCE = 'NO')' xSide = ' 
READ(*,*)xSide
WRITE(*,'(/A)',ADVANCE = 'NO')' ySide = ' 
READ(*,*)ySide
WRITE(*,'(/A)',ADVANCE = 'NO')' zSide = ' 
READ(*,*)zSide
WRITE(*,'(/A)',ADVANCE='NO')' Enter distance between adjacent lattice points : '
READ(*,*)width
WRITE(*,'(/A)',ADVANCE='NO')' Enter the scaling factor(E0) of energy function : '
READ(*,*)E0
WRITE(*,'(/A)',ADVANCE='NO')                                                    &
' Enter the scaling factor(A) of misorientation angle : '
READ(*,*)A
WRITE(*,'(/A)',ADVANCE='NO')' Enter the incremental time step : '
READ(*,*)delT
WRITE(*,'(/A)',ADVANCE='NO')' Enter the mobility(mu) : '
READ(*,*)mu
WRITE(*,'(/A)',ADVANCE='NO')' Enter the radius bounding the gaussian minima : '
READ(*,*)radius

WRITE(*,'(/A)',ADVANCE='NO')' Do you want to use symmetry operations(y/n)? : '
READ(*,*)ans
symm = 0
IF(ans=="y") symm = 1

WRITE(*,'(/A)',ADVANCE='NO')' Want to read pot minima from file(y/n)? : '
READ(*,*)ans
potOption = 0
IF(ans=="y") potOption = 1
WRITE(*,'(/A)',ADVANCE='NO')' Enter the total number of potential minima : ' 

IF(potOption == 1) totalMinEulerAngles = 8
IF(potOption == 0) totalMinEulerAngles = 2

CALL SetInverseGaussian

RETURN
!===============================================================================
                   END SUBROUTINE SetInputs
!===============================================================================


!===============================================================================
                   SUBROUTINE SetConstants
!===============================================================================
IMPLICIT NONE

pi = 4.0D0*DATAN(1.0D0)
factor = 180.0D0/pi      !factor for converting radian to degrees
WRITE(*,*) pi, factor

RETURN
!===============================================================================
                   END SUBROUTINE SetConstants
!===============================================================================


!===============================================================================
                         SUBROUTINE SetInverseGaussian
!This subroutine is meant for reading the depth & Std. Deviation of minimas
!===============================================================================
IMPLICIT NONE

!Local
INTEGER  ::  i

ALLOCATE(depth(totalMinEulerAngles), SD(totalMinEulerAngles))

DO i = 1,totalMinEulerAngles
   WRITE(*,'(A)',ADVANCE='NO')' Enter depth of inverse Gaussian : '
   READ(*,*)depth(i)
   WRITE(*,'(A)',ADVANCE='NO')' Enter SD of inverse Gaussian : '
   READ(*,*)SD(i)
ENDDO

sdNoise = DSQRT(2.D0/(delT *mu))     

RETURN
!===============================================================================
                   END SUBROUTINE SetInverseGaussian
!===============================================================================


!===============================================================================
                   SUBROUTINE FileError(fileName,error,modul)
!Subroutine for file errors
!===============================================================================
IMPLICIT NONE

!Shared
INTEGER::error
CHARACTER(LEN=*)::fileName
CHARACTER(LEN=*),OPTIONAL :: modul

IF(error/=0)THEN
  WRITE(*,1)
  WRITE(*,'(" ===>> FATAL ERROR :: file not found : ")')fileName
  IF(PRESENT(modul))WRITE(*,'(1x,"           Module :: ",A20,/)')modul
  STOP
ENDIF

1 FORMAT(1x,'*==================================================',&
            '==================*') 
!===============================================================================
                   END SUBROUTINE FileError
!===============================================================================


!===============================================================================
                        SUBROUTINE Fatal(routine,string,modul)
!Subroutine for giving the ouput error in the code
!===============================================================================
IMPLICIT NONE

!Shared
CHARACTER(LEN=*)  :: routine, string
CHARACTER(LEN=*),OPTIONAL :: modul


WRITE(*,'(1x, "TATA failed")')

WRITE(*,1)
WRITE(*,*)string
WRITE(*,'(1x,"===>> FATAL ERROR in the routine :",A15)')routine
IF(PRESENT(modul))WRITE(*,'(1x,">>>MODULE<<< ",A30,/)')modul
STOP

1 FORMAT(1x,' ==================================================',&
             '================== ') 
!===============================================================================
                           END SUBROUTINE Fatal
!===============================================================================

                       END MODULE moduleGlobal
!===============================================================================
