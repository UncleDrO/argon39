!> \file
!! Stopping power of alphas

!******************************************************************************
!> Read in stopping power from input file
!!
SUBROUTINE read_stopping_power
  USE all_data
  IMPLICIT NONE
  INTEGER :: unit, i
  CHARACTER :: filename*255

  filename = trim(input_directory)//'/'//trim(file_stopping_power)
  unit = 25
  call open_file(unit, trim(filename), 'old', 'read')
  write (*,*) 'Reading stopping power'
  write (outp,'("Reading stopping power from file: ",a)') trim(filename)

  !< get number of entries
  nsp = 0
  do 
      read (unit, *, end=99)
      nsp = nsp + 1
  end do
99 rewind (unit)

  !< read stopping power
  allocate (stoppow(nsp,2))
  do i=1,nsp
      read (unit,*) stoppow(i,1:2)
  end do
  close (unit)

  !< output stopping power
  write (outp,'(4x,a,4x,a)') 'E in MeV', 'SP in MeV cm2/g'
  do i=1,nsp
      write (outp,'(f12.4,f12.2)') stoppow(i,1:2)
      if (stoppow(i,1) > maxEalpha) exit
  end do
  write (outp,*)

  if (stoppow(nsp,1) < maxEalpha) &
      STOP 'Max stopping power energy lower than maxEalpha!!!'

END SUBROUTINE read_stopping_power


!******************************************************************************
!> Calculate and output alpha particle range
!!
SUBROUTINE calculate_alpha_range
  USE all_data
  IMPLICIT NONE
  INTERFACE
      REAL(DP) FUNCTION func_range(E)
        USE all_data
        USE mhunt
        REAL(DP), INTENT(IN) :: E
      END FUNCTION func_range
  END INTERFACE
  INTEGER :: unit
  REAL(DP) :: a,b,s,range
  CHARACTER :: filename*255

  write (*,*) 'Calculating alpha range'
  filename = 'alpha_range.out'
  unit = 37
  call open_file(unit, trim(filename), 'unknown', 'write')
  write (outp,'("<<OUTFILE>> Range of alpha''s written to: ",a)') trim(filename)
  write (unit,'("#  E[MeV]  Range[micrometers]")')

  range = zero
  b = zero
  do
      a = b
      b = b + dEan
      if (a < stoppow(1,1)) cycle
      if (b > maxEalpha) exit
      call qtrap(func_range,a,b,s)
      range = range + s
      write (unit,'(2f12.3)') b, range*1e4_dp 
  end do

  close (unit)

END SUBROUTINE calculate_alpha_range


!******************************************************************************
!> Function 1 / linear_stopping_power provided for integration of alpha range
!!
REAL(DP) FUNCTION func_range(E)
  USE all_data
  USE mhunt
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: E
  INTERFACE
      REAL(DP) FUNCTION interpolate_tab_func(xa,fa,xi,type)
        USE precision
        USE mhunt
        REAL(DP), DIMENSION(:), INTENT(IN) :: xa,fa
        REAL(DP), INTENT(IN) :: xi
        INTEGER, INTENT(IN) :: type
      END FUNCTION interpolate_tab_func
  END INTERFACE
  REAL(DP) :: lspow
  !< get linear stopping power: log-log interpolation
  lspow = interpolate_tab_func(stoppow(:,1),stoppow(:,2),E,1) * rock_density
  !< integrand for neutron production function (including target atom density)
  func_range = one / lspow
END FUNCTION func_range


