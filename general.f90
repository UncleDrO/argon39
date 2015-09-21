!> \file
!! General subroutines and functions


!******************************************************************************
!> Interpolation of `type` between tabulated values of function `fa` at `xi`
!!
!! type=0: linear-linear
!! type=1: log-log
!! 
REAL(DP) FUNCTION interpolate_tab_func(xa,fa,xi,type)
  USE precision
  USE mhunt
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(IN) :: xa,fa
  REAL(DP), INTENT(IN) :: xi
  INTEGER, INTENT(IN) :: type
  INTEGER :: jlo
  REAL(DP) :: x1,x2,f1,f2,zz
  call hunt(xa(:),xi,jlo)
  jlo = lbound(xa,1) - 1 + jlo
  if (jlo == lbound(xa,1)-1 .or. jlo == ubound(xa,1)) then
      write (*,*) 'Out of bounds in hunt from interpolate_tab_func!'
      write (*,*) xa(lbound(xa,1)), ' to ', xa(ubound(xa,1))
      write (*,*) 'desired: ',xi
      STOP 'TERMINATION.'
  end if
  x1 = xa(jlo)
  x2 = xa(jlo+1)
  f1 = fa(jlo)
  f2 = fa(jlo+1)
  zz = xi
  select case (type)
  case (0)
      interpolate_tab_func = ((x2-zz)*f1 + (zz-x1)*f2) / (x2-x1)
  case (1)
      x1 = log(x1)
      x2 = log(x2)
      f1 = log(f1)
      f2 = log(f2)
      zz = log(zz)
      interpolate_tab_func = ((x2-zz)*f1 + (zz-x1)*f2) / (x2-x1)
      interpolate_tab_func = exp(interpolate_tab_func)
  end select
END FUNCTION interpolate_tab_func


