!===============================================================================
                           MODULE moduleSites
!===============================================================================
!MODIFIED 21.3.07 LINE 652 ADDED SYMMETRY CHECK CONDITION

USE moduleInput
USE modulePotMinima
USE moduleRand
IMPLICIT NONE
!===============================================================================
! Lattice sites:
TYPE, PUBLIC :: sites
  REAL(DP), POINTER :: rotMat(:,:)
  REAL(DP), POINTER :: invRotMat(:,:)
  REAL(DP), POINTER :: Q(:),gradX(:),gradY(:),gradZ(:)
  REAL(DP), POINTER :: misOrienAng(:)
  REAL(DP)          :: delQX,delQY,delQZ
END TYPE sites

TYPE(sites), ALLOCATABLE :: site(:,:,:)


! Matrices product in different directions
TYPE, PUBLIC :: prdir
  REAL(DP), POINTER :: prod(:,:), sMat(:,:)
  REAL(DP)          :: maxTr, tr
  INTEGER           :: symNo, ss
END TYPE prdir
TYPE(prdir), ALLOCATABLE :: dir(:)

!===============================================================================

!Public subroutines
PUBLIC  :: AllocateSite        ! Referred in tata.f90
PUBLIC  :: GetQuaternions      ! tata.f90
PUBLIC  :: GetSlantIntfaceQuat ! tata.f90
PUBLIC  :: GetQuatContinued
PUBLIC  :: GetBiCrystalQuat
PUBLIC  :: GetBandsOfQuat
PUBLIC  :: TempSmallLatticeQuat
PUBLIC  :: GetColumnsOfQuat    !26/02/07
PUBLIC  :: GetGrad             ! tata.f90
PUBLIC  :: GetMaxTrace
PUBLIC  :: GetMisOrienAng      !moduleCalculation
PUBLIC  :: GetQuatDummy
PUBLIC  :: MakeQuatUnique
PUBLIC  :: ResetQuatMagnitude
PRIVATE :: CheckQuaternions    
PRIVATE :: ResetGradMagnitude
PRIVATE :: GetRotMat
PRIVATE :: GetInvRotMat
PRIVATE :: StereoProjection
!===============================================================================
                                 CONTAINS
!===============================================================================


!===============================================================================
                        SUBROUTINE AllocateSite
! To allocate all the variables used in the program
!===============================================================================
IMPLICIT NONE

!Locals
INTEGER :: i, j, k


! Allocating site and its components:
ALLOCATE(site(xSide, ySide, zSide))

DO i = 1, xSide
  DO j = 1, ySide
    DO k = 1, zSide
       ALLOCATE(site(i,j,k)%rotMat(dim,dim),   &
                site(i,j,k)%invRotMat(dim,dim),&
                site(i,j,k)%misOrienAng(dim),  &
                site(i,j,k)%Q(dim+1),          &
                site(i,j,k)%gradX(dim+1),      &
                site(i,j,k)%gradY(dim+1),      &
                site(i,j,k)%gradZ(dim+1)       )
    ENDDO
  ENDDO
ENDDO

RETURN
!===============================================================================
                        END SUBROUTINE AllocateSite
!===============================================================================


!===============================================================================
                   SUBROUTINE GetQuaternions
!To convert the euler angles to Quaternions from the data file of euler angles
! by calling the subroutine ReadEulerAngles which is in moduleInput.f90
!===============================================================================
IMPLICIT NONE
          
!Local variables: 
INTEGER :: i, j, k, e
REAL(DP) :: add, diff, b

CALL ReadEulerAngles
DO k = 1, zSide
  e=0
  DO j = 1, ySide
    DO i = 1, xSide
       e = e+1
       IF( e > totalEulerAngles )&
       CALL Fatal('GetQuaternions','total euler angles exceeded','moduleSites')

       add = 0.5*(euler(e)%a + euler(e)%c)
       diff = 0.5*(euler(e)%a - euler(e)%c)
       b = 0.5*euler(e)%b
!      Bunge euler angles to quaternions conversion
       site(i,j,k)%Q = 0.D0
       site(i,j,k)%Q(1)=DCOS(diff)*DSIN(b)
       site(i,j,k)%Q(2)=DSIN(diff)*DSIN(b)
       site(i,j,k)%Q(3)=DSIN(add)*DCOS(b)
       site(i,j,k)%Q(4)=DCOS(add)*DCOS(b)
       CALL  MakeQuatUnique(i,j,k)
       CALL  CheckQuaternions(i,j,k)
       CALL  ResetQuatMagnitude(i,j,k) 
!       CALL  StereoProjection(i,j,k) !tangent at(1,0,0,0)
    ENDDO
  ENDDO
ENDDO
!CALL TempSmallLatticeQuat! CAUTION: this reduces the lattice size
CALL ReadSymmMatrices

!Deallocating the euler angle array
DEALLOCATE(euler)

RETURN
!===============================================================================
                   END SUBROUTINE GetQuaternions
!===============================================================================

!===============================================================================
                   SUBROUTINE TempSmallLatticeQuat
! Temporarily created to just read the long lattice size from the above 
!mentioned subroutine but only take out a small part of the lattice and observe !the small lattice
!CAUTION: please do not call this subroutine when you wish to execute a Large 
! lattice size on the whole
!===============================================================================
IMPLICIT NONE
          
!Local variables: 
INTEGER :: i, j, k, x, y
x = 0
y = 0

DO k = 1,3
  DO j = 26,75
    y = MOD(j,25)
    IF(j>=50) y = MOD(j,25)+25
    IF(j == 75) y = 50
    DO i = 51,150
      x = MOD(i,50)
      IF( i>=100) x=MOD(i,50)+50
      IF(i == 150) y = 100
      site(x,y,k)%Q = site(i,j,k)%Q
    ENDDO
  ENDDO
ENDDO
xSide = 100; ySide = 50; zSide = 3
RETURN
!===============================================================================
                   END SUBROUTINE TempSmallLatticeQuat
!===============================================================================

!===============================================================================
                   SUBROUTINE GetQuatContinued
!To continue from the middle of an execution which got killed somehow
!===============================================================================
IMPLICIT NONE
          
!Local variables: 

INTEGER :: a,b,i,j,k, counter
CHARACTER(LEN=132) :: line
!please edit the below mentioned file by writing an 'end' at the EOF
fileName='../data/cntdataQuat.dat' 
WRITE(*,*)" Opening file ",fileName," to read euler angles"
OPEN(UNIT=cnt, FILE=fileName,FORM='formatted',STATUS='old',IOSTAT=error)
CALL fileError(fileName,error,'moduleInput')
counter = 0

DO
  READ(cnt,'(A132)')line 
  IF(line(1:3) == 'end') EXIT
  counter = counter + 1
ENDDO

totalEulerAngles = counter
REWIND(cnt)


! Now read the file
DO a = 1, totalEulerAngles
  READ(cnt,*)i,j,k,site(i,j,k)%Q(1),site(i,j,k)%Q(2),site(i,j,k)%Q(3),&
               site(i,j,k)%Q(4)
ENDDO
REWIND(cnt)
CLOSE(cnt)
WRITE(*,*)" Quaternions has been populated successfully"

DO k = 1, zSide
  DO j = 1, ySide
    DO i = 1, xSide
      CALL MakeQuatUnique(i,j,k)
      CALL CheckQuaternions(i,j,k)
      CALL ResetQuatMagnitude(i,j,k)
    ENDDO
  ENDDO
ENDDO 
CALL ReadSymmMatrices
RETURN
!===============================================================================
                   END SUBROUTINE GetQuatContinued
!===============================================================================
!===============================================================================
                    SUBROUTINE StereoProjection(i,j,k)
!To show the Quaternions stereographic projections
!===============================================================================
IMPLICIT NONE
!Shared
INTEGER, INTENT(IN) :: i,j,k

!Locals
REAL(DP) :: x,y,z,coeff
IF(site(i,j,k)%Q(4) == -1.D0) RETURN

coeff = 1.0D0/(1.D0 + site(i,j,k)%Q(4))

x = coeff*site(i,j,k)%Q(1)
y = coeff*site(i,j,k)%Q(2)
z = coeff*site(i,j,k)%Q(3)
WRITE(5000,*)x,y,z



RETURN
!===============================================================================
                   END SUBROUTINE StereoProjection
!===============================================================================

!===============================================================================
                    SUBROUTINE GetSlantIntfaceQuat
!To produce a lattice bicrystal having slant interface between them 
!===============================================================================
IMPLICIT NONE
          
!Local variables: 
INTEGER :: i, j, k
REAL(DP):: inp1, inp2, inp4

WRITE(*,'(/A)', ADVANCE='NO')                                                 &
' Give the 1st,2nd and 4th comp of q on the other half of the cube?'
READ(*,*)inp1,inp2,inp4

DO k=1,zSide
  DO j=1,ySide
    DO i=1,xSide
      site(i,j,k)%Q = 0.D0
      IF( (REAL(i)/REAL(xSide)+REAL(k)/REAL(zSide)) < 1.0D0) THEN
        site(i,j,k)%Q(1) = inp1!0.0d0
        site(i,j,k)%Q(2) = inp2!0.0d0
        site(i,j,k)%Q(3) = DSQRT(1.0d0-(inp4*inp4+inp1*inp1+inp2*inp2))!0.0d0
        site(i,j,k)%Q(4) = inp4
      ELSE
        site(i,j,k)%Q(1) = 0.0D0!inp1
        site(i,j,k)%Q(2) = 0.0D0 !inp2
        site(i,j,k)%Q(3) = 0.0D0 !dsqrt(1.0d0 - inp1*inp1 - inp2*inp2) !0.99d0
        site(i,j,k)%Q(4) = 1.0D0
      ENDIF
     
      CALL CheckQuaternions(i,j,k)

    ENDDO
  ENDDO
ENDDO

CALL ReadSymmMatrices

RETURN
!===============================================================================
                   END SUBROUTINE GetSlantIntfaceQuat
!===============================================================================
!===============================================================================
                    SUBROUTINE GetBiCrystalQuat
!To produce a bicrystal lattice
!===============================================================================
IMPLICIT NONE
          
!Local variables: 
INTEGER :: i, j, k
REAL(DP):: inp1, inp2, inp4

WRITE(*,'(/A)', ADVANCE='NO')                                                 &
' Give the 1st,2nd and 4th comp of q on the other half of the cube?'
READ(*,*)inp1,inp2,inp4

DO k=1,zSide
  DO j=1,ySide
    DO i=1,xSide/3
      site(i,j,k)%Q = 0.D0

      site(i,j,k)%Q(1) = inp1!0.0d0
      site(i,j,k)%Q(2) = inp2!0.0d0
      site(i,j,k)%Q(3) = DSQRT(1.0d0-(inp4*inp4+inp1*inp1+inp2*inp2))!0.0d0
      site(i,j,k)%Q(4) = inp4

      CALL CheckQuaternions(i,j,k)

    ENDDO
  ENDDO
ENDDO

DO k=1,zSide
  DO j=1,ySide
    DO i=xSide/3 + 1, xSide

      site(i,j,k)%Q(1) = 0.0D0!inp1
      site(i,j,k)%Q(2) = 0.0D0 !inp2
      site(i,j,k)%Q(3) = 0.0D0 !dsqrt(1.0d0 - inp1*inp1 - inp2*inp2) !0.99d0
      site(i,j,k)%Q(4) = 1.0D0

      CALL CheckQuaternions(i,j,k)

    ENDDO
  ENDDO
ENDDO

CALL ReadSymmMatrices

RETURN
!===============================================================================
                   END SUBROUTINE GetBiCrystalQuat
!===============================================================================
!===============================================================================
                    SUBROUTINE GetBandsOfQuat
!To produce Bands of different orientations in a lattice
!===============================================================================
IMPLICIT NONE
          
!Local variables: 
INTEGER :: i, j, k, n, m, p
REAL(DP):: x

  x = 0.0 ; i = 0; n = 0; k = 1
  DO m = 1, xSide/6
!    x = grnd()
!    IF(0.00 <= x .AND. x < 0.25) n = 1
!    IF(0.25 <= x .AND. x < 0.5)  n = 2
!    IF(0.50 <= x .AND. x < 0.75) n = 3
!    IF(0.75 <= x .AND. x < 1.00) n = 4  
   n = 2*n + 1
   IF(n > 3) n = 1 !3 because bicrystal 
!    DO k = 1, zSide 
      DO j = 1,ySide
        DO p = 1, 6
         i = p + (m-1)*6
         IF( grnd() > 0.5) THEN
          site(i,j,k)%Q =  potMinima(n)%Quat
            
         ELSEIF( n == 1) THEN
          site(i,j,k)%Q =  potMinima(2)%Quat !for bicrystal
!         site(i,j,k)%Q =  0.D0
!         site(i,j,k)%Q(3)= 1.D0
         ELSEIF( n == 3) THEN
          site(i,j,k)%Q =  potMinima(4)%Quat !for bicrystal
!         site(i,j,k)%Q(1) = 0.4523329753602607
!         site(i,j,k)%Q(2) =-0.15591509359944885
!         site(i,j,k)%Q(3) = 0.348380931307752
!         site(i,j,k)%Q(4) =0.818427478663149
!         ELSEIF( n == 3) THEN
!         site(i,j,k)%Q(1) = 0.320128993850940
!         site(i,j,k)%Q(2) = -.0922472436022094
!         site(i,j,k)%Q(3) = 0.409599404412543
!         site(i,j,k)%Q(4) = 0.889846868907011
!         ELSEIF( n == 4) THEN
!         site(i,j,k)%Q(1) =3.011600557016635E-002
!         site(i,j,k)%Q(2) =0.533096924003365 
!         site(i,j,k)%Q(3) =-1.413735856398096E-002
!         site(i,j,k)%Q(4) =  0.845399805369910
         ENDIF
        CALL  ResetQuatMagnitude(i,j,k)
      ENDDO
!    ENDDO
  ENDDO
ENDDO
DO k = 2, zSide 
  DO j = 1,ySide
      DO i = 1, xSide
        site(i,j,k)%Q = site(i,j,1)%Q 
      ENDDO
  ENDDO
ENDDO



!site(8,20,1)%Q =  potMinima(1)%Quat
!site(8,20,2)%Q =  potMinima(1)%Quat
!site(8,20,3)%Q =  potMinima(1)%Quat
CALL ReadSymmMatrices
RETURN
!===============================================================================
                   END SUBROUTINE GetBandsOfQuat
!===============================================================================


!===============================================================================
                    SUBROUTINE GetColumnsOfQuat         !modified 26.2.07
!To produce columns of different orientations with random thickness
!===============================================================================
IMPLICIT NONE
          
!Local variables: 
INTEGER :: i, j, k, n, m, p, xCol, xFilled, c
REAL:: x

x = 0.0 ; i = 0; n = 0; k = 1
xFilled = 0
DO 
  x = grnd()
  IF(0.00 <= x .AND. x < 0.25) xCol = 4
  IF(0.25 <= x .AND. x < 0.5)  xCol = 5
  IF(0.50 <= x .AND. x < 0.75) xCol = 6
  IF(0.75 <= x .AND. x < 1.00) xCol = 7
  c = 0  
  n = n + 1
  IF(n == 20) xCol = xSide - xFilled
  IF( MOD(n,2) /= 0) c = 2 
  DO j = 1,ySide
    DO p = 1, xCol
      i = p + xFilled
      x = grnd()
      IF(x >= 0.75) site(i,j,k)%Q =  potMinima(2*c+1)%Quat
      IF(x >= 0.5 .AND. x < 0.75) site(i,j,k)%Q =  potMinima(2*c+2)%Quat 
      IF(x >= 0.25 .AND. x < 0.5) site(i,j,k)%Q =  potMinima(2*c+3)%Quat
      IF(x < 0.25) site(i,j,k)%Q =  potMinima(2*c+4)%Quat
      CALL  ResetQuatMagnitude(i,j,k)
    ENDDO
  ENDDO
  xFilled = xFilled + xCol
  IF (xFilled == xSide) EXIT
ENDDO
DO k = 2, zSide 
  DO j = 1,ySide
      DO i = 1, xSide
        site(i,j,k)%Q = site(i,j,1)%Q 
      ENDDO
  ENDDO
ENDDO
CALL ReadSymmMatrices

RETURN
!===============================================================================
                   END SUBROUTINE GetColumnsOfQuat !modified 26.2.07
!===============================================================================


!===============================================================================
                   SUBROUTINE CheckQuaternions(i,j,k)
!To check the magnitude of each quaternion
!===============================================================================
IMPLICIT NONE

!Shared
INTEGER, INTENT(IN) :: i,j,k

!Locals
INTEGER :: m

DO m = 1, 4
  IF(site(i,j,k)%Q(m) > 1.)&
  CALL Fatal('CheckQuaternions','Magnitude > 1','moduleSites')
ENDDO

RETURN
!===============================================================================
                   END SUBROUTINE CheckQuaternions
!===============================================================================


!===============================================================================
                           SUBROUTINE GetGrad
!To get the gradients of the quaternions
!===============================================================================
IMPLICIT NONE
          
!Local: 
INTEGER  :: i, j, k, x, y, z, comp
!INTEGER, INTENT(IN)::symm

DO k = 1, zSide
! Periodic BC
  z = k
  IF( k == zSide ) z = 0
  DO j = 1, ySide
!   Periodic BC
    y =j
    IF( j == ySide ) y = 0
    DO i = 1, xSide
!     Periodic BC
      x = i
      IF( i == xSide ) x = 0

      CALL GetMaxTrace(i,j,k,x,y,z)
       
      site(i,j,k)%delQX=0.D0; site(i,j,k)%delQY=0.D0; site(i,j,k)%delQY=0.D0

      DO comp = 1, 4
        site(i,j,k)%gradX(comp) = (site(x+1,j,k)%Q(comp) - &
                                   site(i,j,k)%Q(comp))/width
        site(i,j,k)%gradY(comp) = (site(i,y+1,k)%Q(comp) - &
                                   site(i,j,k)%Q(comp))/width
        site(i,j,k)%gradZ(comp) = (site(i,j,z+1)%Q(comp) - &
                                    site(i,j,k)%Q(comp) )/width

        site(i,j,k)%delQX = site(i,j,k)%delQX + site(i,j,k)%gradX(comp)**2
        site(i,j,k)%delQY = site(i,j,k)%delQY + site(i,j,k)%gradY(comp)**2
        site(i,j,k)%delQZ = site(i,j,k)%delQZ + site(i,j,k)%gradZ(comp)**2
      ENDDO
      CALL ResetGradMagnitude(i,j,k)
    ENDDO
  ENDDO
ENDDO
RETURN
!===============================================================================
                       END SUBROUTINE GetGrad
!===============================================================================
!===============================================================================
                   SUBROUTINE MakeQuatUnique(i,j,k)
!To make the quaternion rotation unique because Q = -Q , so we restrict it to 
!plus Q
!===============================================================================
IMPLICIT NONE
          
!Shared: 
INTEGER,INTENT(IN)  :: i, j, k
!Local:
IF( site(i,j,k)%Q(4) > 0.D0) RETURN
IF( site(i,j,k)%Q(4) < 0.D0)THEN
  site(i,j,k)%Q = -site(i,j,k)%Q
ELSEIF( site(i,j,k)%Q(4) == 0.D0)THEN
    IF( site(i,j,k)%Q(1)<0.D0)THEN
       site(i,j,k)%Q = -site(i,j,k)%Q
ELSEIF( site(i,j,k)%Q(1)==0.D0)THEN
    IF( site(i,j,k)%Q(2)<0.D0)THEN
       site(i,j,k)%Q = -site(i,j,k)%Q
ELSEIF( site(i,j,k)%Q(2)==0.D0)THEN
    IF( site(i,j,k)%Q(3)<0.D0)THEN
       site(i,j,k)%Q = -site(i,j,k)%Q
ELSEIF( site(i,j,k)%Q(3)==0.D0) THEN
      site(i,j,k)%Q(4)=1.D0
ENDIF 
ENDIF
ENDIF      
ENDIF        
RETURN
!===============================================================================
                       END SUBROUTINE MakeQuatUnique
!===============================================================================

!===============================================================================
                   SUBROUTINE ResetGradMagnitude(i,j,k)
! To scale back the delq's 
!===============================================================================
IMPLICIT NONE

!Shared
INTEGER, INTENT(IN) :: i,j,k

IF(site(i,j,k)%delQX > 4.D0)site(i,j,k)%delQX = 4.D0
IF(site(i,j,k)%delQY > 4.D0)site(i,j,k)%delQY = 4.D0
IF(site(i,j,k)%delQZ > 4.D0)site(i,j,k)%delQZ = 4.D0

RETURN
!===============================================================================
                   END SUBROUTINE ResetGradMagnitude
!===============================================================================


!===============================================================================
                   SUBROUTINE GetRotMat(i,j,k)
!To get the rotation matrix related to a Quaternion 
!===============================================================================
IMPLICIT NONE
          
!Local: 
INTEGER, INTENT(IN)  :: i, j, k
REAL(DP) ::Q11,Q12,Q13,Q14,Q22,Q23,Q24,Q33,Q34,Q44

       Q11 =  site(i,j,k)%Q(1)*site(i,j,k)%Q(1)   
       Q22 =  site(i,j,k)%Q(2)*site(i,j,k)%Q(2)   
       Q33 =  site(i,j,k)%Q(3)*site(i,j,k)%Q(3)
       Q44 =  site(i,j,k)%Q(4)*site(i,j,k)%Q(4)   

       Q12 =  site(i,j,k)%Q(1)*site(i,j,k)%Q(2)
       Q13 =  site(i,j,k)%Q(1)*site(i,j,k)%Q(3)
       Q14 =  site(i,j,k)%Q(1)*site(i,j,k)%Q(4)
       Q23 =  site(i,j,k)%Q(2)*site(i,j,k)%Q(3)
       Q24 =  site(i,j,k)%Q(2)*site(i,j,k)%Q(4)
       Q34 =  site(i,j,k)%Q(3)*site(i,j,k)%Q(4)
     
       site(i,j,k)%rotMat(1,1) = (Q44 + Q11) -(Q22 + Q33)
       site(i,j,k)%rotMat(2,2) = (Q44 + Q22) -(Q11 + Q33)
       site(i,j,k)%rotMat(3,3) = (Q33 + Q44) -(Q11 + Q22)
       site(i,j,k)%rotMat(1,2) = 2*(Q12 + Q34)
       site(i,j,k)%rotMat(2,1) = 2*(Q12 - Q34)
       site(i,j,k)%rotMat(1,3) = 2*(Q13 - Q24)
       site(i,j,k)%rotMat(3,1) = 2*(Q13 + Q24)
       site(i,j,k)%rotMat(2,3) = 2*(Q23 + Q14)
       site(i,j,k)%rotMat(3,2) = 2*(Q23 - Q14)

RETURN
!===============================================================================
                    END SUBROUTINE GetRotMat
!===============================================================================


!===============================================================================
                   SUBROUTINE GetInvRotMat(i,j,k)
!To get back the inverse rotation matrix
!===============================================================================
IMPLICIT NONE
          
!Shared: 
INTEGER, INTENT(IN)  :: i, j, k
!Local
REAL(DP) ::Q11,Q12,Q13,Q14,Q22,Q23,Q24,Q33,Q34,Q44

       Q11 =  site(i,j,k)%Q(1)*site(i,j,k)%Q(1)   
       Q22 =  site(i,j,k)%Q(2)*site(i,j,k)%Q(2)   
       Q33 =  site(i,j,k)%Q(3)*site(i,j,k)%Q(3)
       Q44 =  site(i,j,k)%Q(4)*site(i,j,k)%Q(4)   

       Q12 =  site(i,j,k)%Q(1)*site(i,j,k)%Q(2)
       Q13 =  site(i,j,k)%Q(1)*site(i,j,k)%Q(3)
       Q14 =  site(i,j,k)%Q(1)*site(i,j,k)%Q(4)
       Q23 =  site(i,j,k)%Q(2)*site(i,j,k)%Q(3)
       Q24 =  site(i,j,k)%Q(2)*site(i,j,k)%Q(4)
       Q34 =  site(i,j,k)%Q(3)*site(i,j,k)%Q(4)

       site(i,j,k)%invRotMat(1,1) =  (Q44 + Q11) -(Q22 + Q33)
       site(i,j,k)%invRotMat(2,2) =  (Q44 + Q22) -(Q11 + Q33)
       site(i,j,k)%invRotMat(3,3) =  (Q33 + Q44) -(Q11 + Q22)
       site(i,j,k)%invRotMat(1,2) =  2*(Q12 - Q34)
       site(i,j,k)%invRotMat(2,1) =  2*(Q12 + Q34)
       site(i,j,k)%invRotMat(1,3) =  2*(Q13 + Q24)
       site(i,j,k)%invRotMat(3,1) =  2*(Q13 - Q24)
       site(i,j,k)%invRotMat(2,3) =  2*(Q23 - Q14)
       site(i,j,k)%invRotMat(3,2) =  2*(Q23 + Q14)

RETURN
!===============================================================================
                    END SUBROUTINE GetInvRotMat
!===============================================================================


!===============================================================================
                   SUBROUTINE GetMaxTrace(i,j,k,x,y,z)
!To get the Maximum trace for finding minimum misorientation
!===============================================================================
IMPLICIT NONE

!Shared:
INTEGER, INTENT(IN)  :: i, j, k, x, y, z

!Local: 
INTEGER ::  sss, d

IF(symm /= 0 .AND. symm /= 1)CALL Fatal("GetMaxTrace","Wrong symmetry option")


ALLOCATE( dir(dim) )

dir(1)%maxTr = -100.0D0
dir(2)%maxTr = -100.0D0
dir(3)%maxTr = -100.0D0

CALL GetRotMat(i,j,k)

direction:DO d = 1, dim

  ALLOCATE(dir(d)%prod(dim,dim))

  IF( d==1)THEN
    CALL GetInvRotMat(x+1,j,k)
    dir(d)%prod = MATMUL(site(x+1,j,k)%invRotMat,site(i,j,k)%rotMat)
  ENDIF
  IF( d==2)THEN
    CALL GetInvRotMat(i,y+1,k)
    dir(d)%prod = MATMUL(site(i,y+1,k)%invRotMat,site(i,j,k)%rotMat)
  ENDIF
  IF( d==3)THEN
    CALL GetInvRotMat(i,j,z+1)
    dir(d)%prod = MATMUL(site(i,j,z+1)%invRotMat,site(i,j,k)%rotMat)
  ENDIF
      
  SELECT CASE(symm)

!   No symmetry case
    CASE(0)
      IF( dir(d)%maxTr < 3.0d0)THEN
        dir(d)%tr = dir(d)%prod(1,1)+dir(d)%prod(2,2)+dir(d)%prod(3,3)
        IF( dir(d)%maxTr < dir(d)%tr) dir(d)%maxTr = dir(d)%tr
      ENDIF

!   Symmetry case
    CASE(1)
      dir(d)%symNo=0; dir(d)%ss = 0;

      AllSymm: DO  sss=1,totSymm

        ALLOCATE(dir(d)%sMat(dim,dim))

        dir(d)%ss=sss
!        IF( j.NE.1 .OR. k.NE.1) dir(1)%ss = 1
!        IF( k .NE. 1) dir(2)%ss = 1

        IF( dir(d)%maxTr < 3.0d0)THEN
          dir(d)%sMat = MATMUL(dir(d)%prod,symmMat(dir(d)%ss)%S)      
          dir(d)%tr = dir(d)%sMat(1,1)+dir(d)%sMat(2,2)+dir(d)%sMat(3,3)
          IF( dir(d)%maxTr < dir(d)%tr) THEN
            dir(d)%maxTr = dir(d)%tr
            dir(d)%symNo = dir(d)%ss
          ENDIF
        ENDIF

      DEALLOCATE(dir(d)%sMat)

      IF( dir(d)%maxTr == 3.0d0)CYCLE

     ENDDO AllSymm

  END SELECT
     
  DEALLOCATE (dir(d)%prod)

ENDDO direction 

CALL GetMisOrienAng(i,j,k,symm)
!IF(symm == 1) CALL GetQuatDummy(i,j,k)   

DEALLOCATE(dir)

RETURN
!===============================================================================
                    END SUBROUTINE GetMaxTrace
!===============================================================================


!===============================================================================
            SUBROUTINE GetMisOrienAng(i,j,k,symm)
!To get the misorientation angle
!===============================================================================
IMPLICIT NONE
!Shared
INTEGER, INTENT(IN) :: i, j, k, symm
!Local
REAL(DP),ALLOCATABLE,DIMENSION(:)::deg
INTEGER :: d

! Calculation of the minimum misorientation angles

ALLOCATE (deg(dim))

DO d = 1,dim
  IF( dir(d)%maxTr < 3.0d0) THEN
    IF( dir(d)%maxTr > -1.0d0) THEN
      site(i,j,k)%misOrienAng(d)=DACOS(0.5d0*(dir(d)%maxTr- 1.0d0))
    ELSE
     site(i,j,k)%misOrienAng(d)=pi 
    ENDIF
  ELSE
     site(i,j,k)%misOrienAng(d)=0.0d0
  ENDIF
!  Converting from radian to degrees  minimum angles
   IF(symm == 0) CYCLE
   deg(d) = site(i,j,k)%misOrienAng(d) * factor
   IF(deg(d) > 62.8)&
   CALL Fatal('GetMisOrienAng','DEGREE >62.8;SYMM FAILS','moduleSites')
  
ENDDO
        
DEALLOCATE(deg)

RETURN
!==============================================================================
                    END SUBROUTINE GetMisOrienAng
!==============================================================================


!==============================================================================
           SUBROUTINE GetQuatDummy(i,j,k)
!To get  the symmetric Quaternion i.e the quaternion update after symmetry
!is performed on the quaternion
! This is still in preparation 
!CAUTION: DON'T USE THIS SUBROUTINE AT PRESENT
!==============================================================================
IMPLICIT NONE

TYPE ::matRot 
  REAL(DP), POINTER :: R(:,:)
END TYPE matRot

TYPE(matRot), ALLOCATABLE :: axis(:)

INTEGER,INTENT(IN)::i,j,k
REAL(DP)::a,b,c,d,denom,tr1,s1,magQ
INTEGER::x,y,z,m


CALL Fatal('GetQuatDummy','dont use this subroutine','moduleSites')
ALLOCATE (axis(dim))        

DO m =1,dim
   IF( m == 1) THEN
     x = i + 1
     y = j
     z = k
     IF(i == xSide) x = 1
     IF( j /= 1 .OR. k /= 1) CYCLE
   ELSEIF( m == 2) THEN
     x = i 
     y = j + 1
     z = k
     IF(j == ySide) y = 1
     IF( k /= 1) CYCLE
   ELSEIF( m == 3) THEN
     x = i 
     y = j 
     z = k + 1
     IF(k == zSide) z = 1
   ENDIF

   ALLOCATE (axis(m)%R(dim,dim))
   CALL GetInvRotMat(x,y,z)
   axis(m)%R = MATMUL(symmMat(dir(m)%symNo)%S,site(x,y,z)%invRotMat)
   tr1 = axis(m)%R(1,1)+axis(m)%R(2,2)+axis(m)%R(3,3)+1.0D0

   IF( tr1.gt.0.0D0)THEN
     s1 = 0.5D0/DSQRT(tr1)
     site(x,y,z)%Q(1) = 0.25D0/s1
     site(x,y,z)%Q(2) = (axis(m)%R(3,2) -axis(m)%R(2,3)) * s1 
     site(x,y,z)%Q(3) = (axis(m)%R(1,3) - axis(m)%R(3,1)) * s1  
     site(x,y,z)%Q(4) = (axis(m)%R(2,1) - axis(m)%R(1,2)) * s1 
   ELSEIF( axis(m)%R(1,1)>axis(m)%R(2,2).AND.axis(m)%R(1,1)>axis(m)%R(3,3))THEN
     s1 = DSQRT(1.0D0 + axis(m)%R(1,1)-(axis(m)%R(2,2) + axis(m)%R(3,3)))*2.0D0
     b = axis(m)%R(1,2)+axis(m)%R(2,1)
     c = axis(m)%R(1,3)+axis(m)%R(3,1)
     d = axis(m)%R(3,2)-axis(m)%R(2,3)
     site(x,y,z)%Q(1) = d/s1
     site(x,y,z)%Q(2) = 0.25D0 * s1
     site(x,y,z)%Q(3) = b/s1
     site(x,y,z)%Q(4) = c/s1
   ELSEIF( axis(m)%R(2,2)>axis(m)%R(3,3))THEN
     s1 =DSQRT(1.0D0 + axis(m)%R(2,2) - (axis(m)%R(1,1) + axis(m)%R(3,3)))*2.0D0
     b = axis(m)%R(1,2)+axis(m)%R(2,1)
     c = axis(m)%R(2,3)+axis(m)%R(3,2)
     d = axis(m)%R(1,3)-axis(m)%R(3,1)
     site(x,y,z)%Q(1) = d/s1
     site(x,y,z)%Q(2) = b/s1
     site(x,y,z)%Q(3) = 0.25D0*s1
     site(x,y,z)%Q(4) = c/s1
   ELSE
     s1 =DSQRT(1.0D0 + axis(m)%R(3,3)-(axis(m)%R(2,2) + axis(m)%R(1,1)))*2.0D0 
     b = axis(m)%R(1,3)+axis(m)%R(3,1)
     c = axis(m)%R(2,3)+axis(m)%R(3,2)
     d = axis(m)%R(2,1)-axis(m)%R(1,2)
     site(x,y,z)%Q(1) = d/s1
     site(x,y,z)%Q(2) = b/s1
     site(x,y,z)%Q(3) = c/s1
     site(x,y,z)%Q(4) = 0.25D0 * s1
   ENDIF
   
   DEALLOCATE(axis(m)%R)
   CALL ResetQuatMagnitude(x,y,z)
ENDDO

DEALLOCATE(axis)

RETURN
!==============================================================================
                    END SUBROUTINE GetQuatDummy
!==============================================================================

!===============================================================================
                   SUBROUTINE ResetQuatMagnitude(x,y,z)
! To rescale the magnitude to 1
!===============================================================================
IMPLICIT NONE

!Shared
INTEGER, INTENT(IN) :: x,y,z

!Local
INTEGER :: i
REAL(DP)::magQ
   
magQ = 0.0D0
DO i=1,4
  magQ=magQ + site(x,y,z)%Q(i)*site(x,y,z)%Q(i)
ENDDO
magQ=DSQRT(magQ)

IF( magQ >1.0D0 .OR. magQ < 1.0D0) THEN
  DO i=1,4
    site(x,y,z)%Q(i)=site(x,y,z)%Q(i)/magQ
  ENDDO
ENDIF

RETURN
!===============================================================================
                   END SUBROUTINE ResetQuatMagnitude
!===============================================================================


                       END MODULE moduleSites
!==============================================================================


