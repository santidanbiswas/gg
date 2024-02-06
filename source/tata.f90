!===============================================================================
                              PROGRAM TATA
!===============================================================================
USE moduleSites
USE modulePotMinima
USE moduleCalculation

IMPLICIT NONE

INTEGER :: i, j, k, l, iter
REAL(DP):: time, finalT, divsor
INTEGER(4):: rate, timeStart,timeEnd
INTEGER :: seconds, totalTime, latOption

!SYSTEM CLOCK 
CALL SYSTEM_CLOCK(count=timeStart,count_rate=rate)

CALL SetInputs
WRITE(*,'(/A)',ADVANCE='NO') 'Option: read lattice orientn from file &
or slant interface or bicrsytal or bands or columns or continued code&
 (1/2/3/4/5/6)? :' 
READ(*,*)latOption

WRITE(*,'(/A)',ADVANCE='NO')' How long do you want to evolve the quaternions?:'
READ(*,*)finalT

WRITE(*,'(/A)',ADVANCE='NO') ' Enter divisor for noise : '
READ(*,*)divsor

CALL AllocateSite
IF(potOption == 1)CALL GetPotMinQuaternions !this subroutine reads from datafile
IF(potOption == 0)CALL GetPotMinQuat     !for slant interface(doesn't work now)

SELECT CASE(latOption)

 CASE(1)
 CALL GetQuaternions

 CASE(2)
 CALL GetSlantIntfaceQuat

 CASE(3)
 CALL GetBiCrystalQuat

 CASE(4)
 CALL GetBandsOfQuat

 CASE(5)
 CALL GetColumnsOfQuat

 CASE(6)
 CALL GetQuatContinued

END SELECT

CALL SetRandom

CALL SYSTEM_CLOCK(count=timeEnd,count_rate=rate)
seconds = (timeEnd-timeStart)/rate 
WRITE(*,'(1x,"Elapsed CPU time so far = ",I6," seconds")')seconds 

!IMPORTANT PART OF THE PROGRAM 
time = 0.D0
iter = 0

fileName='../output/graingrowth.dat'
OPEN(UNIT = ggdata, FILE = fileName )

fileName='../output/energy.dat'
OPEN(UNIT = enrgy, FILE = fileName )

fileName='../data/output.dat'
OPEN(UNIT = outpt, FILE = fileName )
CALL GetEnergy(iter)
CALL GetOutput(iter)
totalTime = 0.0
WRITE(*,'(/1x,"Iteration",6x,"time",2x,"TimeElapsed")')
DO
  IF (finalT < time) EXIT

  IF(MOD(iter,20)==0)CALL SYSTEM_CLOCK(count=timeStart,count_rate=rate)

  time = time + delT
  iter = iter + 1   
  CALL GetGrad
  CALL GetQuatEvolution(divsor,iter)
  IF (MOD(iter,500) == 0) CALL GetOutput(iter)
!  IF( iter > 300 .AND. iter < 500)THEN
!    IF( MOD(iter,10) == 0) CALL GetOutput(iter)
!  ENDIF

  IF(MOD(iter,20)==0)THEN 
    CALL SYSTEM_CLOCK(count=timeEnd,count_rate=rate)
    seconds = (timeEnd-timeStart)/rate 
    totalTime = totalTime + seconds
    WRITE(*,'(1x,I9,1x, F9.4, 2x, I10)')iter,time,totalTime
  ENDIF
ENDDO
END FILE enrgy
REWIND enrgy
CLOSE  (enrgy)
END FILE ggdata
REWIND ggdata
CLOSE(ggdata)
END FILE outpt
REWIND outpt
CLOSE(outpt)


WRITE(*,'(1x,"Elapsed CPU time = ",I6," seconds")')totalTime

END
