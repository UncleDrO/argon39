!******************************************************************************
!> This program prepares input files and batch execution file for running TALYS
!!
!! Used in conjucntion with ARGON39 and TALYS codes.
!!
!! @author  Ondřej Šrámek, University of Maryland, sramek@umd.edu
!! @date    5 December 2013
!!
PROGRAM TALYSprep
  USE all_data
  USE nucl_type
  IMPLICIT NONE
  CHARACTER :: argv*127, command*255, fb*127, filename*127, dir*255, elemsymb*2
  INTEGER :: len, ntargets, batch=2, it, A

  print *
  print *, 'TALYSprep'
  print *

  call get_inputs(argv,len,command)

  print *
  print *, 'This program may create new directories in your file system'
  print *, 'and overwrite some existing files.'
  print *, 'You can kill it now (ctrl-c). Otherwise press ENTER.'
  read *
  print *

  call open_file(outp, 'TALYSprep.out', 'replace', 'write')

  call read_elements_nuclides

  call get_requested_an_targets(ntargets)

  fb = 'TALYSbatch_'//trim(adjustl(model_name))//'.sh'
  write (*,*) 'Batch file name: ', trim(adjustl(fb))
  call open_file(batch, fb, 'unknown', 'write')
  write (batch,'(a,/)') "#! /bin/bash"
  write (batch,'(a,/)') 'homedir=`pwd`'
  write (*,*) 'TALYS directory: ', trim(adjustl(talys_directory))

  !< create talys_directory if needed
  write (command,'(5a)') 'if [ ! -d ', trim(adjustl(talys_directory)), &
      ' ] ; then mkdir ', trim(adjustl(talys_directory)),' ; fi'
  call execute_command_line(trim(command))

  do it=1,ntargets
      write (*,*) trim(adjustl(targets(it)%nucid))
      write (batch,'("### ",a)') trim(adjustl(targets(it)%nucid))
      dir = trim(adjustl(talys_directory))//'/'//trim(adjustl(targets(it)%nucid))
      write (outp,*)
      write (outp,*) 'nuclide directory:  ', trim(adjustl(dir))
      write (*,*)    'nuclide directory:  ', trim(adjustl(dir))

      !< create nuclide subdirectory if needed
      write (command,'(5a)') 'if [ ! -d ', trim(adjustl(dir)), &
          ' ] ; then mkdir ', trim(adjustl(dir)),' ; fi'
      call execute_command_line(trim(command))

      !< execute TALYS
      write (batch,'(2a)') 'cd ', trim(adjustl(dir))
      write (batch,'(a)') 'pwd'
      write (batch,'(a)') 'talys < input > output'
      write (batch,'(a,/)') 'cd $homedir'

      !< TALYS input file
      filename = trim(adjustl(dir))//'/input'
      write (outp,*) 'TALYS input file:   ', trim(adjustl(filename))
      write (*,*)    'TALYS input file:   ', trim(adjustl(filename))
      elemsymb = element( targets(it)%Z )%symb
      A = targets(it)%A
      call write_talys_input (elemsymb, A, filename)

      !< TALYS energy file
      filename = trim(adjustl(dir))//'/erange'
      write (outp,*) 'Energy discr. file: ', trim(adjustl(filename))
      write (*,*)    'Energy discr. file: ', trim(adjustl(filename))
      call write_talys_erange(dEtalys, filename)
  end do

  close (batch)

  !< change permissions to executable
  write (command,'(2a)') 'chmod 755 ', trim(adjustl(fb))
  call execute_command_line(trim(command))

  print *
  print *, 'PROGRAM FINISHED SUCCESFULLY.'
  print *
  print *, 'NOW EXECUTE SHELL SCRIPT:   ', trim(adjustl(fb))
  print *

  write (outp,*)
  write (outp,*) 'PROGRAM FINISHED SUCCESFULLY.'
  write (outp,*)
  write (outp,*) 'NOW EXECUTE SHELL SCRIPT:   ', trim(adjustl(fb))
  write (outp,*)

END PROGRAM TALYSprep


!******************************************************************************
!> Writes TALYS input file
!!
SUBROUTINE write_talys_input (elemsymb, A, filename)
  USE all_data
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN) :: elemsymb !< element symbol
  CHARACTER(*), INTENT(IN) :: filename !< filename
  INTEGER, INTENT(IN)      :: A        !< mass number
  INTEGER :: u=22

  call open_file(u, filename, 'unknown', 'write')
  write (u,'(a)') '#'
  write (u,'(a)') '# General'
  write (u,'(a)') '#'
  write (u,'(a)') 'projectile a'
  write (u,'(2a)') 'element ', trim(adjustl(elemsymb))
  write (u,'(a,i3)') 'mass ', A
  write (u,'(a)') 'energy erange'
  write (u,'(a)') '#'
  write (u,'(a)') '# Compound nucleus'
  write (u,'(a)') '#'
  write (u,'(a)') 'widthfluc y'
  write (u,'(a)') '#'
  write (u,'(a)') '# Output'
  write (u,'(a)') '#'
  write (u,'(a)') 'outspectra y'
  write (u,'(a)') 'filespectrum n'
  close (u)

END SUBROUTINE write_talys_input


!******************************************************************************
!> Writes TALYS energy file
!!
SUBROUTINE write_talys_erange (dE, filename)
  USE all_data
  IMPLICIT NONE
  REAL(DP), INTENT(IN)     :: dE       !< energy stepping
  CHARACTER(*), INTENT(IN) :: filename !< filename
  INTEGER :: u=22
  REAL(DP) :: Emax=8.8_dp, E

  call open_file(u, filename, 'unknown', 'write')
  E = zero
  do
      E = E + dE
      if (E > Emax+0.5*dE) exit
      write (u,'(f7.3)') E
  end do
  close (u)

END SUBROUTINE write_talys_erange
