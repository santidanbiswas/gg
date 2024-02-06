!===============================================================================
                          MODULE moduleInput
!===============================================================================
USE moduleGlobal 

IMPLICIT NONE

! Euler angles:
TYPE, PUBLIC :: eAngle
    REAL(DP) :: a, b, c
END TYPE eAngle
TYPE(eAngle), ALLOCATABLE :: euler(:)

! Symmetry matrices:
TYPE, PUBLIC :: symmMatrix
    REAL(DP), POINTER :: S(:,:),B(:,:)
END TYPE symmMatrix
TYPE(symmMatrix), ALLOCATABLE :: symmMat(:)

! Potential minima euler angles:
TYPE, PUBLIC :: minEangle
    REAL(DP) :: a1, b1, c1
END TYPE minEangle
TYPE(minEangle), ALLOCATABLE :: minEuler(:)

!===============================================================================
!Public subroutines
PUBLIC  :: ReadEulerAngles  !moduleSites  Subroutine: GetQuaternions
PUBLIC  :: ReadSymmMatrices !moduleSites  Subroutine: GetQuaternions
PUBLIC  :: ReadContinueQuat !moduleSites  Subroutine: GetQuatContinued
PUBLIC  :: ReadMinEulerAngles !modulePotMinima
!===============================================================================
                                 CONTAINS

!===============================================================================
                         SUBROUTINE ReadEulerAngles

!To read the euler angles from the below mentioned file				!===============================================================================
IMPLICIT NONE

!Local variables:
REAL(DP) :: a, b, c
INTEGER :: i, counter
CHARACTER(LEN=132) :: line
!CAUTION: Please Write 'end' at the end of the below mentioned file

fileName='../data/200_100_3.dat'
WRITE(*,*)" Opening file ",fileName," to read euler angles"
OPEN(UNIT=eang, FILE=fileName,FORM='formatted',STATUS='old',IOSTAT=error)
CALL fileError(fileName,error,'moduleInput')

! Find the number of lines
counter = 0
DO
  READ(eang,'(A132)')line 
  IF(line(1:3) == 'end') EXIT
  counter = counter + 1
ENDDO

totalEulerAngles = counter
REWIND(eang)

ALLOCATE(euler(totalEulerAngles))

! Now read the file
DO i = 1, totalEulerAngles
  READ(eang,*)a,b,c 
  euler(i)%a = a
  euler(i)%b = b
  euler(i)%c = c
ENDDO
totalEulerAngles = counter
REWIND(eang)
CLOSE(eang)
WRITE(*,*)" Euler angles has been populated successfully"
RETURN
!===============================================================================
                   END SUBROUTINE ReadEulerAngles
!===============================================================================
!===============================================================================
                   SUBROUTINE ReadContinueQuat
!This subroutine is only used when for some reason the execution got killed
!and you want to continue from that time onwardby taking the final set of
!Quaternion as the initial one.
!CAUTION:Please remember to rename your old output files otherwise the new execution will overwrite them and hence the old output datas will be lost
!===============================================================================
IMPLICIT NONE

!Local variables:
INTEGER :: a,b,i,j,k, counter
CHARACTER(LEN=132) :: line
!please edit the below mentioned file by writing an 'end' at the EOF
fileName='../data/output.dat' 
WRITE(*,*)" Opening file ",fileName," to read euler angles"
OPEN(UNIT=outpt, FILE=fileName,FORM='formatted',STATUS='old',IOSTAT=error)
CALL fileError(fileName,error,'moduleInput')
counter = 0

DO
  READ(outpt,'(A132)')line 
  IF(line(1:3) == 'end') EXIT
  counter = counter + 1
ENDDO

totalEulerAngles = counter
REWIND(outpt)

ALLOCATE(euler(totalEulerAngles))

! Now read the file
DO a = 1, totalEulerAngles
  READ(outpt,*)i,j,k!,(site(i,j,k)%Q(b), b =1,4)
ENDDO
REWIND(outpt)
CLOSE(outpt)
WRITE(*,*)" Quaternions has been populated successfully"

   
RETURN
!===============================================================================
                   END SUBROUTINE ReadContinueQuat
!===============================================================================
!===============================================================================
                         SUBROUTINE ReadSymmMatrices
!To read symmetry elements and create 24 symmetry matrices
!===============================================================================
IMPLICIT NONE

!Local variables:
INTEGER :: i, j, k, m, n
REAL(DP), ALLOCATABLE, DIMENSION(:) :: element

fileName='../data/symm.dat'
WRITE(*,*)" Opening file ",fileName," to read the symmetry matrices"
OPEN(UNIT=sym, FILE=fileName,FORM='formatted',STATUS='old',IOSTAT=error)
CALL fileError(fileName,error,'moduleInput')

ALLOCATE(element(totalSymmElements))
ALLOCATE(symmMat(totSymm))

DO k = 1, totSymm
  READ(sym,*)(element(m), m = 1, totalSymmElements)
  ALLOCATE(symmMat(k)%S(dim,dim))
  ALLOCATE(symmMat(k)%B(dim,dim))
  n = 0
  DO i = 1, dim
    DO j = 1, dim 
      n = n + 1
      symmMat(k)%B(i,j) = element(n)
    ENDDO
  ENDDO
    symmMat(k)%S = TRANSPOSE(symmMat(k)%B)
  DEALLOCATE(symmMat(k)%B)
ENDDO

DEALLOCATE(element) 

WRITE(*,*)" Symmetry matrices has been populated successfully"

RETURN
!===============================================================================
                   END SUBROUTINE ReadSymmMatrices
!===============================================================================
  
!===============================================================================
                         SUBROUTINE ReadMinEulerAngles
!To read the minima euler angles
!===============================================================================
IMPLICIT NONE

!Local variables:
REAL(DP) :: a1, b1, c1
INTEGER :: i

CALL SetConstants

fileName='../data/minimaEulerAngles.dat'

WRITE(*,*)" Openning file ",fileName," to read the minimum euler angles"

OPEN(UNIT=mineang, FILE=fileName,FORM='formatted',STATUS='old',IOSTAT=error)
CALL fileError(fileName,error,'moduleInput')

ALLOCATE(minEuler(totalMinEulerAngles))

! Now read the file
DO i = 1, totalMinEulerAngles
  READ(mineang,*)a1,b1,c1
  minEuler(i)%a1 = a1/factor  !To convert from degrees to radian
  minEuler(i)%b1 = b1/factor
  minEuler(i)%c1 = c1/factor
ENDDO
WRITE(*,*)" Minimum euler angles has been populated successfully"

RETURN
!===============================================================================
                   END SUBROUTINE ReadMinEulerAngles
!===============================================================================
                      
                         END MODULE moduleInput
!===============================================================================
 
