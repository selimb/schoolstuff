! ************************************************
!
! 	MECH 539 Project 3
!       Generate Grid to Solve TSD Equations
!
!       Siva Nadarajah
!       3 March 2016
!
! ************************************************

program generategrid
implicit none
integer            :: i,j
integer, parameter :: imax=60,jmax=40 ! imax must be divisible by 2
double precision, parameter :: del_min=0.025D0 ! Minimum Grid Spacing
double precision x(1:imax), y(1:jmax)

!  Generate distribution of x-coordinates
x(imax/2)=0.0D0

!  Loop over from imax/2 to imax.
do i=(imax/2+1),imax
  if (x(i-1)<=0.50D0) then
    ! If x is on the airfoil surface then use constant delta
    x(i)=x(i-1)+del_min
  else
    ! If x is downstream then stretch the grid until the right farfield boundary
    x(i)=x(i-1)+(x(i-1)-x(i-2))*1.5
  end if

  if (i.lt.imax) x(imax-i)=-x(i) ! Mirror the grid about imax/2.
end do

x=x+20.5D0 ! Move the midpoint of the grid to 20.5

!  Generate distribution of y-coordinates
!- Y-Grid

! Initialize first and second y coordinates.
y(1)=0.0D0
y(2)=del_min
do i=3,jmax
  ! Stretch the grid until the upper farfield boundary
  y(i)=y(i-1)+(y(i-1)-y(i-2))*1.1**1.25D0
end do

open(20, file="grid")
write(20, *) imax, jmax
write(20, *) (x(j), j=1,imax)
write(20, *) (y(j), j=1,jmax)

end program generategrid
