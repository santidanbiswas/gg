!===============================================================================
                           MODULE moduleCalculation
!===============================================================================
USE moduleSites
USE modulePotMinima
IMPLICIT NONE
!===============================================================================
TYPE,PUBLIC :: temporaryQ
  REAL(DP), POINTER :: dummyQ(:)
END TYPE temporaryQ

TYPE(temporaryQ), ALLOCATABLE :: temp(:,:,:)

!Public subroutines
PUBLIC  :: GetQuatEvolution ! tata.f90
PUBLIC  :: GetEnergy
PUBLIC  :: GetPotContribution
PUBLIC  :: GetQuatUpdate
PUBLIC  :: GetOutput        !tata.f90
PUBLIC  :: SetRandom
PUBLIC  :: GetRandNoise
PUBLIC  :: NumberQuat
!Private subroutines
PRIVATE :: GetCosDelta
!Private functions
PRIVATE :: ranF
PRIVATE :: metric
!===============================================================================
                                 CONTAINS
!===============================================================================


!===============================================================================
                   SUBROUTINE GetQuatEvolution(divsor,iter)
!To calculate the quaternion evolution
!===============================================================================
IMPLICIT NONE

!Local:
INTEGER  :: i, j, k, x, y, z, m, indx
REAL(DP), INTENT(IN) ::divsor
INTEGER, INTENT(IN):: iter
REAL(DP) :: sumnoise, delQ, gradQ, minusDelQ, minusGradQ, P1
REAL(DP) :: T1, T2, sumF1, noiseterm 
REAL(DP), ALLOCATABLE, DIMENSION(:) :: F1, rate, noise, term1 
REAL(DP), ALLOCATABLE, DIMENSION(:) :: V
REAL(DP) :: sfn, potbarrier

ALLOCATE (F1(4) , rate(4), noise(4), term1(dim))
ALLOCATE (temp(xSide,ySide,zSide))
ALLOCATE (V(4))
IF (MOD(iter,500) == 0)WRITE(enrgy,*) 'iter' , iter

DO k = 1, zSide
  DO j = 1, ySide
    DO i = 1, xSide


      ALLOCATE (temp(i,j,k)%dummyQ(4))
      temp(i,j,k)%dummyQ = 0.0
      CALL GetRandNoise(divsor, noise) 
      sumnoise=0.0d0
      DO indx=1,4
        sumnoise=sumnoise+noise(indx)*site(i,j,k)%Q(indx)
      ENDDO
      DO indx=1,4
        F1(indx) = 0.D0
        DO m =1,dim
          IF( m == 1) THEN
             x = i - 1
             y = j
             z = k
            IF(i == 1) x = xSide
            delQ = site(i,j,k)%delQX
            minusDelQ = site(x,y,z)%delQX
            gradQ = site(i,j,k)%gradX(indx)
            minusGradQ =site(x,y,z)%gradX(indx) 
          ELSEIF( m == 2) THEN
            x = i 
            y = j - 1 
            z = k
            IF(j == 1) y = ySide
            delQ = site(i,j,k)%delQY
            minusDelQ = site(x,y,z)%delQY
            gradQ = site(i,j,k)%gradY(indx)
            minusGradQ =site(x,y,z)%gradY(indx) 
          ELSEIF( m == 3) THEN
            x = i 
            y = j 
            z = k - 1
            IF(k == 1) z = zSide
            delQ = site(i,j,k)%delQZ
            minusDelQ = site(x,y,z)%delQZ
            gradQ = site(i,j,k)%gradZ(indx)
            minusGradQ =site(x,y,z)%gradZ(indx) 
          ENDIF
          
          CALL GetCosDelta(delQ,P1)

          IF( DABS(P1) < 1.0d0)THEN
            IF(delQ > 2.D0) &
            site(i,j,k)%misorienAng(m)=2.D0*pi - site(i,j,k)%misorienAng(m)

            T1 = ( 1.D0-(DTANH(A*site(i,j,k)%misOrienAng(m)))**2.D0 )*        &
                 (4.0d0-2.0d0*delQ)/DSQRT(1.0d0 - P1*P1)*gradQ
          ELSE 
            T1 = 0.D0
          ENDIF

          CALL GetCosDelta(minusDelQ,P1)

          IF( DABS(P1) < 1.0d0)THEN
            IF(minusDelQ > 2.D0) &
            site(i,j,k)%misorienAng(m)=2.D0*pi - site(i,j,k)%misorienAng(m)

            T2 = ( 1.D0-(DTANH(A*site(x,y,z)%misOrienAng(m)))**2.D0 )*        &
                      (4.0d0-2.0d0*minusDelQ)/DSQRT(1.0d0 - P1*P1)            &
                       * minusGradQ
          ELSE
            T2 = 0.D0
          ENDIF

          term1(m) = T1 - T2

          F1(indx) = F1(indx) + term1(m)
        ENDDO
      ENDDO

      sumF1=0.0d0
      DO indx = 1, 4
        sumF1=sumF1 + site(i,j,k)%Q(indx)*F1(indx)
      ENDDO 
      CALL GetPotContribution(i,j,k,V,sfn,iter)
        
!      V = 0.0D0 ! caution remove after check
      DO indx = 1, 4
        potbarrier = site(i,j,k)%Q(indx)*sfn
!        potbarrier = 0.0D0 ! caution remove after check
        noiseterm = site(i,j,k)%Q(indx)*sumnoise - noise(indx)
!        noiseterm = 0.D0 ! caution remove after check
        
        rate(indx) = mu*E0*A*((F1(indx) + potbarrier) - &
                     (site(i,j,k)%Q(indx)*sumF1 + V(indx))) - noiseterm
        temp(i,j,k)%dummyQ(indx) = site(i,j,k)%Q(indx) + rate(indx)*delT
      ENDDO
    ENDDO
  ENDDO
ENDDO
IF (MOD(iter,500) == 0)WRITE(enrgy,*) ''

DEALLOCATE (F1, rate, noise, term1, V)

CALL GetQuatUpdate

RETURN
!===============================================================================
                   END SUBROUTINE GetQuatEvolution
!===============================================================================


!===============================================================================
                   SUBROUTINE GetEnergy(iter)
!To get the Energy from the Gaussian minimas
!===============================================================================
IMPLICIT NONE

!Local:
INTEGER  :: i, j, k
INTEGER, INTENT(IN):: iter
REAL(DP), ALLOCATABLE, DIMENSION(:) :: V
REAL(DP) :: sfn, potbarrier

ALLOCATE (V(4))
WRITE(enrgy,*) 'iter' , iter
DO k = 1, zSide
  DO j = 1, ySide
    DO i = 1, xSide

      CALL GetPotContribution(i,j,k,V,sfn,iter)!chk
    
    ENDDO
  ENDDO
ENDDO

!===============================================================================
                   END SUBROUTINE GetEnergy
!===============================================================================


!===============================================================================
                   SUBROUTINE GetCosDelta(tempP1,P1)
!===============================================================================
IMPLICIT NONE

!Shared
REAL(DP),INTENT(IN)  ::tempP1
REAL(DP),INTENT(OUT) ::P1

P1 =1.0D0 - 2.0D0 *tempP1 + 0.5D0*tempP1*tempP1

RETURN
!===============================================================================
                   END SUBROUTINE GetCosDelta
!===============================================================================

!===============================================================================
                   SUBROUTINE GetPotContribution(i,j,k,V,sfn,iter)
!To calculate the potential contribution
!===============================================================================
IMPLICIT NONE

!Shared
REAL(DP), INTENT(OUT) :: V(4), sfn
INTEGER, INTENT(IN)   :: i, j, k, iter

!Local
REAL(DP) ::sf, tt, tExp, sdTerm, tDummy, expContrbn,gaussianMinima
INTEGER  ::m, n

sfn = 0.D0
V = 0.D0
gaussianMinima = 0.D0
DO m = 1,totalMinEulerAngles
  sf  = 0.D0; tExp = 0.D0
  DO n = 1,4
    sf = sf + site(i,j,k)%Q(n)*potMinima(m)%Quat(n)
    tt = site(i,j,k)%Q(n) - potMinima(m)%Quat(n)
    tExp=tExp + tt*tt
  ENDDO
  sf=1.0d0 - sf    
  tExp = - tExp
  sdTerm = SD(m)*SD(m)
  expContrbn = DEXP(tExp/(2.0D0*sdTerm))
  tDummy = depth(m) * expContrbn/sdTerm
  gaussianMinima = gaussianMinima +  (depth(m) * expContrbn)
  sfn = sfn + tDummy*sf
  DO n = 1,4
    V(n) = V(n) + tDummy * (site(i,j,k)%Q(n) -potMinima(m)%Quat(n))
  ENDDO
ENDDO
IF (MOD(iter,500) == 0) WRITE(enrgy,*) i,j,k, -gaussianMinima
!WRITE(400,*) site(i,j,k)%Q(3), site(i,j,k)%Q(4), -gaussianMinima

RETURN
!===============================================================================
                   END SUBROUTINE GetPotContribution
!===============================================================================

!===============================================================================
                   SUBROUTINE SetRandom
!To generate random number for the random noise
!===============================================================================
IMPLICIT NONE

REAL(DP) :: C,CD,CM,U(97)
REAL(DP) :: S,T
INTEGER  ::I97,J97,I,J,K,M,L,II,JJ
COMMON /xt101/ U
COMMON /xt102/ C,CD,CM
COMMON /xt103/ I97,J97

I=12
J=34
K=56
L=78
DO II = 1,97
      S = 0.0d0
      T = 0.5d0
      DO JJ = 1,24
          M = MOD(MOD(I*J,179)*K,179)
          I = J
          J = K
          K = M
          L = MOD(53*L+1,169)
          IF(MOD(L*M,64)>=32) S=S+T
          T = 0.5d0*T
      ENDDO
      U(II) = S
ENDDO
C  =   362436.0d0/16777216.0d0
CD =  7654321.0d0/16777216.0d0
CM = 16777213.0d0/16777216.0d0
I97 = 97
J97 = 33
RETURN
!===============================================================================
                   END SUBROUTINE SetRandom
!===============================================================================


!===============================================================================
                   SUBROUTINE GetRandNoise(divsor, rf)
!To generate the random noise
!===============================================================================
IMPLICIT NONE
          
REAL(DP),INTENT(OUT):: rf(4)
REAL(DP),INTENT(IN)::divsor
REAL(DP)::t,deltaT,rk,rsq,v1,v2,fac,deviate,gset
INTEGER::iset,j
SAVE iset,gset
DATA iset/0/

DO j = 1,4
!  Returns a normally distributed deviate with zero mean and variance sigma
   IF( iset == 0) THEN
891   v1 = 2.0d0*ranF()-1.0d0
      v2 = 2.0d0*ranF()-1.0d0
      rsq = v1**2.0d0+v2**2.0d0
      IF (rsq >= 1.0d0 .OR. rsq == 0.0d0) GOTO 891
      fac=DSQRT(-2.0d0*(DLOG(rsq))/rsq)
      gset=v1*fac
      deviate=v2*fac
      iset=1
   ELSE
      deviate=gset
      iset=0
   ENDIF
   rf(j)=deviate*sdNoise/divsor
ENDDO

RETURN
!===============================================================================
                      END SUBROUTINE GetRandNoise
!===============================================================================       

!===============================================================================
                      FUNCTION ranF
!===============================================================================
IMPLICIT NONE

REAL(DP):: C,CD,CM,U(97),UNI,ranF
INTEGER :: I97,J97
COMMON /xt101/ U              
COMMON /xt102/ C,CD,CM        
COMMON /xt103/ I97,J97        

UNI = U(I97)-U(J97)
IF(UNI < 0.D0) UNI=UNI + 1.D0
U(I97) = UNI
I97 = I97-1
IF(I97 == 0) I97=97
J97 = J97-1
IF(J97 == 0) J97=97
C = C-CD
IF(C < 0.D0) C=C+CM
UNI = UNI-C
IF(UNI < 0.D0) UNI=UNI+1.D0
ranF = UNI

RETURN
!===============================================================================
                           END FUNCTION ranF
!===============================================================================


!===============================================================================
                      SUBROUTINE GetQuatUpdate
!To do the Quaternionic update after the time evolution
!===============================================================================
IMPLICIT NONE

INTEGER:: i, j, k
DO k = 1,zSide
  DO j = 1,ySide
    DO i = 1,xSide
!      site(i,j,k)%Q = 0.0D0
      site(i,j,k)%Q = temp(i,j,k)%dummyQ
      CALL MakeQuatUnique(i,j,k)  
      CALL ResetQuatMagnitude(i,j,k)
      DEALLOCATE(temp(i,j,k)%dummyQ)
    ENDDO
  ENDDO
ENDDO

DEALLOCATE(temp)

RETURN
!===============================================================================
                      END SUBROUTINE GetQuatUpdate
!===============================================================================


!===============================================================================
                      SUBROUTINE GetOutput(iter)
!To print the output file
! The output.dat file in the data directory contains the Quaternion value of 
!only the last iteration.
!===============================================================================
IMPLICIT NONE
!Local:
INTEGER  :: i, j, k
INTEGER,INTENT(IN) :: iter
!fileName='../data/output.dat'

!OPEN(UNIT = outpt, FILE = fileName )

WRITE(ggdata,*)'iter',iter
!REWIND outpt
WRITE(outpt,*)'iter',iter
DO k = 1, zSide
  DO j = 1, ySide
! DO j = 20,40
    DO i = 1, xSide
!   DO i = 120,140
      WRITE(outpt,'(3(I5,1x),4(E23.15,1X))')i,j,k,site(i,j,k)%Q 
     
     CALL NumberQuat(i,j,k,iter)   
    ENDDO
  ENDDO
ENDDO

!END FILE outpt
!REWIND outpt
!CLOSE(outpt)



RETURN                   
!===============================================================================
                      END SUBROUTINE GetOutput
!===============================================================================

!===============================================================================
                   SUBROUTINE NumberQuat(i,j,k,iter)
!IMPORTANT:
!This subroutine assigns to each Quaternion 4-vector a single scalar value
!by calculating the metric distance of the quaternions from the minimas
!===============================================================================
IMPLICIT NONE
!Local
REAL(DP) :: geodesic, x
INTEGER, INTENT(IN)  :: i,j,k,iter
INTEGER  :: n
REAL     :: nn


geodesic = 100.D0
nn = 0.0


DO n = 1,totalMinEulerAngles
   x = metric(i,j,k,n)
   IF(x < geodesic) THEN 
     geodesic = x
     nn =n !6.4*real(n-1) + geodesic
   ENDIF
ENDDO
IF(MOD(nn,2.)==0.) nn = nn - 1.0
IF(geodesic > 0.1772)nn = nn 
IF(geodesic > 0.1772.AND. geodesic <= 1.0472) nn = nn + 0.5
IF(geodesic > 1.0472)nn = nn + 1.0
WRITE(ggdata,*)i,j,k,nn

RETURN
!===============================================================================
                   END SUBROUTINE NumberQuat
!===============================================================================

!===============================================================================
                   FUNCTION  metric(i,j,k,n)
!To calculate metric distance
!===============================================================================
IMPLICIT NONE
!Local
Real(DP) :: theta,cosTheta,summ,metric
INTEGER, INTENT(IN)  :: i,j,k,n
INTEGER  :: m
summ =0.D0
DO m = 1,4
  summ = summ + (site(i,j,k)%Q(m) - potMinima(n)%Quat(m))**2
ENDDO
metric = DACOS(1.D0 - 2.D0*summ + 0.5D0*summ*summ)
IF(summ >= 2.D0) metric= 2.D0*pi - metric
RETURN
!===============================================================================
                   END FUNCTION metric
!===============================================================================


                    END MODULE moduleCalculation
!==============================================================================

