!> \file
!! Deals with input/output files


!******************************************************************************
!> Starts outputfile after calling get_inputs
!!
SUBROUTINE get_inputs_start_output
  USE all_data
  IMPLICIT NONE
  CHARACTER :: argv*127, command*255, date*8, time*10, line*127
  INTEGER :: len

  call get_inputs(argv,len,command)

  !< open and initiate output file
  output_file = 'argon39_'//trim(adjustl(model_name))//'.out'
  call open_file(outp, trim(output_file), 'replace', 'write')
  call date_and_time(date, time)
  write (outp,'(a)') 'Starting program ARGON39'
  write (outp,'(6a,/)') 'Command line "',trim(command),'" executed on ',date,' at ',time

  !< copy input file to output file
  write (outp,'(a,a)') argv(1:len),':'
  call open_file(11, argv(1:len), 'old', 'read')
  do
      read (11, '(a)', end=99) line
      write (outp,'(a)') trim(line)
  end do
99 close (11)
  write (outp,*)

END SUBROUTINE get_inputs_start_output


!******************************************************************************
!> Gets inputs from inputfile
!!
SUBROUTINE get_inputs(argv,len,command)
  USE all_data
  IMPLICIT NONE
  CHARACTER, INTENT(OUT) :: argv*127, command*255
  INTEGER, INTENT(OUT) :: len

  !< get input file name
  call get_inputfile_from_command_line(argv,len,command)

  !< read inputs
  call open_file(11, argv(1:len), 'old', 'read')
  read (11, nml=inputs)
  close (11)

END SUBROUTINE get_inputs


!******************************************************************************
!> Subroutine that extracts inputfile name from commnad line
!!
!! Assumes the program is executed with one command line argument
!! which is the inputfile name.
!!
!! Uses Fortran 2003 features and may not work with all F90 compilers.
!! Works with gfortran (gcc version 4.8.0).
!!
SUBROUTINE get_inputfile_from_command_line(argv,len,command)
  USE all_data
  IMPLICIT NONE
  CHARACTER(127), INTENT(OUT) :: argv    !< inputfile name (first argument)
  INTEGER, INTENT(OUT)        :: len     !< length of inputfile name
  CHARACTER(255), INTENT(OUT)  :: command !< program execution command
  INTEGER :: status, argc

  call get_command (command, len, status)
  if (status /= 0) then
      write (*,*) 'ERROR: get_command failed with status = ', status
      STOP
  end if
  write (*,*) 'command line = ', command (1:len)

  argc = command_argument_count ()
  if (argc /= 1) then
      write (*,*) 'USAGE:  ./program inputfile'
      STOP
  end if

  call get_command_argument (1, argv, len, status)
  if (status /= 0) then
      write (*,*) 'ERROR: get_command_argument failed with status = ', status
      STOP
  end if

END SUBROUTINE get_inputfile_from_command_line


!******************************************************************************
!> Opens external file including error handling
!!
SUBROUTINE open_file(u,filename,st,act)
  INTEGER, INTENT(IN) :: u
  CHARACTER(*), INTENT(IN) :: filename, st, act
  INTEGER :: ios
  open (unit=u, file=filename, status=st, action=act, iostat=ios)
  if (ios /= 0) then
      write (*,*) 'Problem opening file ',filename
      STOP
  end if
END SUBROUTINE open_file


!******************************************************************************
!> Closes output file and prints final message on stdout
!!
SUBROUTINE finalize_output
  USE all_data
  IMPLICIT NONE
  INTEGER :: it, ic, unit
  REAL(DP) :: rysumit, rysumtot
  REAL(DP), ALLOCATABLE :: rytotic(:)
  CHARACTER(127) :: filename

  allocate (rytotic(1:ndecc))
  do ic=1,ndecc
      rytotic(ic) = sum(results(:,ic,3))
  end do
  rysumtot = sum(rytotic(:))

  write (outp,'(a)') "SUMMARY OUTPUT TABLES"

  !< neutron yield [per decay]
  filename = 'RESULTS_nyield_per_decay.out'
  unit = 60
  call open_file(unit, trim(filename), 'unknown', 'write')
  write (outp,111) "neutron yield [per decay]"
  write (outp,'("<<OUTFILE>> Also written to file: ",a,/)') trim(filename)

  write (outp,222) "target", dchain(:)%tpnucid
  write (unit,222) "target", dchain(:)%tpnucid
  do it=1,ntgts
      write (outp,333) targets(it)%nucid, results(it,:,1)
      write (unit,333) targets(it)%nucid, results(it,:,1)
  end do
  write (outp,333) 'SF   ', results(ntgts+1,:,1)
  write (unit,333) 'SF   ', results(ntgts+1,:,1)
  close (unit)

  !< neutron prod. rate / (yr * kg-rock * wtf-parent-elem * wtf-target-elem)
  filename = 'RESULTS_nrate_yr-kg-wtf-wtf.out'
  unit = 61
  call open_file(unit, trim(filename), 'unknown', 'write')
  write (outp,111) "neutron prod. rate / (yr * kg-rock * wtf-parent-elem * wtf-target-elem)"
  write (outp,'("<<OUTFILE>> Also written to file: ",a,/)') trim(filename)

  write (outp,222) "target", dchain(:)%tpnucid
  write (unit,222) "target", dchain(:)%tpnucid
  do it=1,ntgts
      write (outp,333) targets(it)%nucid, results(it,:,4)
      write (unit,333) targets(it)%nucid, results(it,:,4)
  end do
  write (outp,333) 'SF   ', results(ntgts+1,:,4)
  write (unit,333) 'SF   ', results(ntgts+1,:,4)
  close (unit)

  !< neutron prod. rate / (yr * kg-rock)
  filename = 'RESULTS_nrate_yr-kg.out'
  unit = 62
  call open_file(unit, trim(filename), 'unknown', 'write')
  write (outp,111) "neutron prod. rate / (yr * kg-rock)"
  write (outp,'("<<OUTFILE>> Also written to file: ",a,/)') trim(filename)

  write (outp,222) "target", dchain(:)%tpnucid, "SUM  "
  write (unit,222) "target", dchain(:)%tpnucid, "SUM  "
  do it=1,ntgts
      rysumit = sum(results(it,:,3))
      write (outp,333) targets(it)%nucid, results(it,:,3), rysumit
      write (unit,333) targets(it)%nucid, results(it,:,3), rysumit
  end do
  rysumit = sum(results(ntgts+1,:,3))
  write (outp,333) 'SF   ', results(ntgts+1,:,3), rysumit
  write (unit,333) 'SF   ', results(ntgts+1,:,3), rysumit
  write (outp,444) "TOTAL", rytotic, rysumtot
  write (unit,444) "TOTAL", rytotic, rysumtot
  close (unit)
  deallocate (rytotic)

111 format (/,a)
222 format (a,1x,4a12)
333 format (a5,4x,4es12.3)
444 format (a5,4x,4es12.3,/)

  write (outp,111) 'PROGRAM FINISHED SUCCESFULLY.'
  close (outp)

  write (*,*) '************************************'
  write (*,*) 'See output in file:  ',trim(output_file)
  write (*,*) '************************************'
END SUBROUTINE finalize_output


