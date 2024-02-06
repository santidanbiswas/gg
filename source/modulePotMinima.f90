!===============================================================================
                           MODULE modulePotMinima
!===============================================================================
USE moduleInput
!USE moduleSites

IMPLICIT NONE
!===============================================================================

! Free energy potential minima:
TYPE, PUBLIC :: potentialminima
 REAL(DP), POINTER :: Quat(:)
END TYPE potentialminima

TYPE(potentialminima), ALLOCATABLE :: potMinima(:)
!===============================================================================

!Public subroutines
PUBLIC  :: GetPotMinQuaternions
PUBLIC  :: GetPotMinQuat
PUBLIC  :: CheckPotMinQuaternions
!PUBLIC  :: GetAnglePotMinima
PRIVATE  :: MakePotQuatUnique
PRIVATE :: ResetPotQuatMagnitude
PRIVATE :: PotStereoProjection
!===============================================================================
                                 CONTAINS
!===============================================================================

!===============================================================================
                   SUBROUTINE GetPotMinQuaternions
!To read potential minima euler angles from the file
!===============================================================================
IMPLICIT NONE
          
!Local variables: 
INTEGER :: i
REAL(DP) :: add, diff, b

CALL ReadMinEulerAngles

ALLOCATE(potMinima(totalMinEulerAngles))
DO i = 1,totalMinEulerAngles 
       
       add = 0.5D0*(minEuler(i)%a1 + minEuler(i)%c1)
       diff= 0.5D0*(minEuler(i)%a1 - minEuler(i)%c1)
       b = 0.5D0*minEuler(i)%b1

       ALLOCATE(potMinima(i)%Quat(4))

       potMinima(i)%Quat(1)=DCOS(diff)*DSIN(b)
       potMinima(i)%Quat(2)=DSIN(diff)*DSIN(b)
       potMinima(i)%Quat(3)=DSIN(add)*DCOS(b)
       potMinima(i)%Quat(4)=DCOS(add)*DCOS(b)
       CALL  MakePotQuatUnique(i)
       CALL  CheckPotMinQuaternions(i)
       CALL  ResetPotQuatMagnitude(i)
!WRITE(*,*)"PotMinima", potMinima(i)%Quat
ENDDO 

!Deallocating the euler angle array
DEALLOCATE(minEuler)
CALL GetAnglePotMinima
CALL PotStereoProjection
RETURN

!===============================================================================
                   END SUBROUTINE GetPotMinQuaternions
!===============================================================================
!===============================================================================
                   SUBROUTINE GetPotMinQuat
!CAUTION: DON"T CALL THIS SUBROUTINE
!===============================================================================
IMPLICIT NONE
          
!Local variables: 
INTEGER :: i,j

ALLOCATE(potMinima(totalMinEulerAngles))

DO i = 1, totalMinEulerAngles
  ALLOCATE(potMinima(i)%Quat(4))
  DO j = 1, 4
    IF( i == 1) THEN
!      potMinima(i)%Quat(j) = site(1,1,1)%Q(j)
    ELSE  
!      potMinima(i)%Quat(j) = site(13,1,10)%Q(j)
    ENDIF
  ENDDO

  CALL  CheckPotMinQuaternions(i)

ENDDO

RETURN
!===============================================================================
                   END SUBROUTINE GetPotMinQuat
!===============================================================================
!===============================================================================
                   SUBROUTINE CheckPotMinQuaternions(i)
!To check Quaternion(of the minima) magnitude
!===============================================================================
IMPLICIT NONE

!Shared
INTEGER, INTENT(IN) :: i

!Locals
INTEGER :: m

DO m = 1, 4
  IF(potMinima(i)%Quat(m) > 1)&
  CALL Fatal('CheckPotMinQuaternions','Magnitude > 1','modulePotMinima')
ENDDO

RETURN
!===============================================================================
                   END SUBROUTINE CheckPotMinQuaternions
!===============================================================================
!===============================================================================
                   SUBROUTINE GetAnglePotMinima
!To get the angle between the minimas.
!===============================================================================
IMPLICIT NONE
!Local
Real(DP) :: theta,cosTheta,add
INTEGER  :: i,j,m
WRITE(100,*)'no.','   potMinima'
DO i = 1, totalMinEulerAngles
  WRITE(100,'(I1,1X,4F16.13)') i,potMinima(i)%Quat
ENDDO
WRITE(100,*)'Dot pdt b/w ','  Angle in deg','  costheta'
DO i = 1, (totalMinEulerAngles - 1)
  DO j = (i+1) ,totalMinEulerAngles
    add = 0.D0
    DO m = 1,4
      add = add + (potMinima(i)%Quat(m) - potMinima(j)%Quat(m))**2
    ENDDO
    cosTheta = 1.D0 - 2.D0*add + 0.5D0*add*add
    theta = (DACOS(cosTheta))*factor
    IF(add > 2.D0) theta = 360 - theta
    WRITE(100,'(2I1,5X,2(F16.13,1x))')i,j,theta,cosTheta
  ENDDO 
ENDDO





RETURN
!===============================================================================
                   END SUBROUTINE GetAnglePotMinima
!===============================================================================
!===============================================================================
                   SUBROUTINE ResetPotQuatMagnitude(x)
!To reset the Quaternion magnitude
!===============================================================================
IMPLICIT NONE

!Shared
INTEGER, INTENT(IN) :: x

!Local
INTEGER :: i
REAL(DP)::magQ
   
magQ = 0.0D0
DO i=1,4
  magQ=magQ + potMinima(x)%Quat(i)*potMinima(x)%Quat(i)
ENDDO
magQ=DSQRT(magQ)

IF( magQ >1.0D0 .OR. magQ < 1.0D0) THEN
  DO i=1,4
    potMinima(x)%Quat(i) = potMinima(x)%Quat(i)/magQ
  ENDDO
ENDIF

RETURN
!===============================================================================
                   END SUBROUTINE ResetPotQuatMagnitude
!===============================================================================
!===============================================================================
                   SUBROUTINE PotStereoProjection
! To take the stereographic projections
!===============================================================================
IMPLICIT NONE
!Local
INTEGER :: i
REAL(DP)::x,y,z,coeff
fileName = '../output/stPrjcnPot.dat'
OPEN(UNIT = potstereo,FILE = fileName)
DO i = 1, totalMinEulerAngles
  IF( potMinima(i)%Quat(4) == -1.D0) CYCLE
  coeff = 1.D0/(1.D0 + potMinima(i)%Quat(4)) 
  x =coeff* potMinima(i)%Quat(1)
  y =coeff* potMinima(i)%Quat(2)
  z =coeff* potMinima(i)%Quat(3)
  WRITE(potstereo,*)x,y,z
  WRITE(potstereo,*)''
  WRITE(potstereo,*)''
ENDDO
REWIND potstereo
CLOSE(potstereo)
RETURN
!===============================================================================
                   END SUBROUTINE PotStereoProjection
!===============================================================================
!===============================================================================
                   SUBROUTINE MakePotQuatUnique(i)
!To make the Potential minima Quaternions unique
!===============================================================================
IMPLICIT NONE
          
!Shared: 
INTEGER,INTENT(IN)  :: i
!Local:
IF( potMinima(i)%Quat(4) > 0.D0)RETURN
IF( potMinima(i)%Quat(4) < 0.D0)THEN
  potMinima(i)%Quat = -potMinima(i)%Quat
ELSEIF( potMinima(i)%Quat(4) == 0.D0)THEN
    IF( potMinima(i)%Quat(1)<0.D0)THEN
       potMinima(i)%Quat = -potMinima(i)%Quat
ELSEIF( potMinima(i)%Quat(1)==0.D0)THEN
    IF( potMinima(i)%Quat(2)<0.D0)THEN
       potMinima(i)%Quat = -potMinima(i)%Quat
ELSEIF( potMinima(i)%Quat(2)==0.D0)THEN
    IF( potMinima(i)%Quat(3)<0.D0)THEN
      potMinima(i)%Quat = -potMinima(i)%Quat
ELSEIF( potMinima(i)%Quat(3)==0.D0) THEN
      potMinima(i)%Quat(4)=1.D0
ENDIF
ENDIF
ENDIF
ENDIF
RETURN
!===============================================================================
                       END SUBROUTINE MakePotQuatUnique
!===============================================================================


                         END MODULE modulePotMinima
!===============================================================================
