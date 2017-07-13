      module clouddata
        implicit none
        integer :: i, nlines
        real*8, dimension(1 : 128) :: INTEMP, INPRES ! PHOENIX pressure/temperature profile [dyn/cm^2, K]
C        integer,parameter :: count = 21
C        integer,parameter :: count2 = 1400
C        real*8, dimension (1 : count) :: teplota,telak ! temperature, pressure
C        real*8, dimension (1 : count2) :: stredni,cislo ! mean radius, number density
C        data teplota/ count * 0.0/
C        data telak/ count * 0.0/
C        data stredni/ count2 * 0.0/
C        data cislo/ count2 * 0.0/
      end module clouddata
