!> \file
!! Table lookup subroutine `hunt` from Numerical Recipes.

!******************************************************************************
!> Interface module for subroutine `hunt`.
!!
MODULE mhunt
  INTERFACE
      SUBROUTINE hunt(xx,x,jlo)
        USE precision
        INTEGER, INTENT(INOUT) :: jlo
        REAL(DP), INTENT(IN) :: x
        REAL(DP), DIMENSION(:), INTENT(IN) :: xx
      END SUBROUTINE hunt
  END INTERFACE
END MODULE mhunt

!******************************************************************************
!> Table lookup subroutine `hunt` from Numerical Recipes.
!!
!! Given an array `xx(1:N)`, and given a value `x`, returns a value `jlo` such 
!! that `x` is between `xx(jlo)` and `xx(jlo+1)`. `xx` must be monotonic, 
!! either increasing or decreasing. `jlo=0` or `jlo=N` is returned to indicate 
!! that `x` is out of range. `jlo` on input is taken as the initial guess for 
!! `jlo` on output.
!!
!! @author W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery: 
!! Numerical Recipes in Fortran 90
!!
SUBROUTINE hunt(xx,x,jlo)
  USE precision
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: jlo
  REAL(DP), INTENT(IN) :: x
  REAL(DP), DIMENSION(:), INTENT(IN) :: xx
  INTEGER :: n,inc,jhi,jm
  LOGICAL :: ascnd
  n=size(xx)
  ascnd = (xx(n) >= xx(1))
  if (jlo <= 0 .or. jlo > n) then
      jlo=0
      jhi=n+1
  else
      inc=1
      if (x >= xx(jlo) .eqv. ascnd) then
          do
              jhi=jlo+inc
              if (jhi > n) then
                  jhi=n+1
                  exit
              else
                  if (x < xx(jhi) .eqv. ascnd) exit
                  jlo=jhi
                  inc=inc+inc
              end if
          end do
      else
          jhi=jlo
          do
              jlo=jhi-inc
              if (jlo < 1) then
                  jlo=0
                  exit
              else
                  if (x >= xx(jlo) .eqv. ascnd) exit
                  jhi=jlo
                  inc=inc+inc
              end if
          end do
      end if
  end if
  do
      if (jhi-jlo <= 1) then
          if (x == xx(n)) jlo=n-1
          if (x == xx(1)) jlo=1
          exit
      else
          jm=(jhi+jlo)/2
          if (x >= xx(jm) .eqv. ascnd) then
              jlo=jm
          else
              jhi=jm
          end if
      end if
  end do
END SUBROUTINE hunt
