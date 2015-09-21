!> \file
!! MCNP6prep program

!******************************************************************************
!> Module for input/output paramters pulled from ARGON39 code
!!
MODULE prepare
  USE inpout
  INTEGER :: nchains, ntargets, batch
  CHARACTER(5), ALLOCATABLE :: chains(:), targets(:)
  !< ******** START PARAMETERS ********
  REAL(DP) :: K39atomcm3 = zero
  REAL(DP) :: Kwtfrac = zero
  REAL(DP) :: Mg24atomcm3 = zero
  REAL(DP) :: Mgwtfrac = zero
  !< ******** END PARAMETERS ********
END MODULE prepare


!******************************************************************************
!> This program prepares input files and batch execution file for running MCNP6
!!
!! Used in conjucntion with ARGON39 and MCNP6 codes.
!!
!! @author  Ondřej Šrámek, University of Maryland, sramek@umd.edu
!! @date    1 October 2013
!!
PROGRAM MCNP6prep
  USE prepare
  IMPLICIT NONE
  INTEGER :: ic, it
  CHARACTER(255) :: fbatch, mcnp6file, command

  call get_inputs
  write (*,*)

  write (*,*) 'ENTER THE ATOMIC DENSITY OF 39K'
  write (*,*) '- run MCNP6 using the “I” option with PRINT card in the input file'
  write (*,*) '    i.e.:  $ mcnp6 i name=modelname'
  write (*,*) '- get atom fraction of 39K from table 40  =: AF'
  write (*,*) '- get the total cell atom density from table 50  =: AD'
  write (*,*) '- enter here AF*AD value:'
  read (*,*) K39atomcm3
  write (*,*) 'Using:', K39atomcm3
  write (*,*)

  write (*,*) 'ENTER THE WEIGHT FRACTION OF K'
  write (*,*) '- sum of mass fractions of all K nuclides from MCNP6 table 40:'
  read (*,*) Kwtfrac
  write (*,*) 'Using:', Kwtfrac
  write (*,*)

  write (*,*) 'ENTER THE ATOMIC DENSITY OF 24Mg'
  write (*,*) '- get atom fraction of 24Mg from table 40  =: AF'
  write (*,*) '- get the total cell atom density from table 50  =: AD'
  write (*,*) '- enter here AF*AD value:'
  read (*,*) Mg24atomcm3
  write (*,*) 'Using:', Mg24atomcm3
  write (*,*)

  write (*,*) 'ENTER THE WEIGHT FRACTION OF Mg'
  write (*,*) '- sum of mass fractions of all Mg nuclides from MCNP6 table 40:'
  read (*,*) Mgwtfrac
  write (*,*) 'Using:', Mgwtfrac
  write (*,*)

  call get_requested_decay_chains

  call get_requested_an_targets

  fbatch = 'MCNP6batch_'//trim(adjustl(model_name))//'.sh'
  write (*,*) 'Batch file name: ', trim(adjustl(fbatch))
  call open_file(batch, fbatch, 'unknown', 'write')
  write (batch,'(a,/)') "#! /bin/bash"

  !< (a,n) neutrons
  do it=1,ntargets
      write (batch,'("### ",a)') trim(adjustl(targets(it)))
      do ic=1,nchains
          mcnp6file = 'mcnp6_'//&
              trim(adjustl(model_name))//'_'//&
              trim(adjustl(chains(ic)))//'_'//&
              trim(adjustl(targets(it)))
          write (batch,'("mcnp6 name=",a)') trim(adjustl(mcnp6file))
          call write_mcnp6_input(it,ic,mcnp6file)
      end do
  end do

  !< spontaneous fission neutrons
  write (batch,'("### ",a)') 'SF'
  do ic=1,nchains
      mcnp6file = 'mcnp6_'//&
          trim(adjustl(model_name))//'_'//&
          trim(adjustl(chains(ic)))//'_SF'
      write (batch,'("mcnp6 name=",a)') trim(adjustl(mcnp6file))
      call write_mcnp6_input(ntargets+1,ic,mcnp6file)
  end do

  close (batch)

  !< change permissions to executable
  write (command,'(2a)') 'chmod 755 ', trim(adjustl(fbatch))
  call execute_command_line(trim(command))

  print *
  print *, 'NOW EXECUTE SHELL SCRIPT:   ', trim(adjustl(fbatch))
  print *


END PROGRAM MCNP6prep


!******************************************************************************
!> Writes MCNP6 input file
!!
SUBROUTINE write_mcnp6_input(it,ic,filename)
  USE prepare
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: it,ic
  CHARACTER(*), INTENT(IN) :: filename
  INTEGER :: u
  u = 22
  call open_file(u, filename, 'unknown', 'write')

  write(u,'(a,a)') 'Production of 39Ar by 39K(n,p) in a rock - ',trim(adjustl(model_name))
  write(u,'(a)') 'c'
  write(u,'(a)') 'c ******** BLOCK 1: CELL CARDS ************************************************'
  write(u,'(a,f5.2,a)') '  1 10 ', -rock_density, ' -1    $ rock of -rho[g/cm3] inside of sphere'
  write(u,'(a)') '  2  0       1     $ void outside of sphere'
  write(u,'(a)') ''
  write(u,'(a)') 'c ******** BLOCK 2: SURFACE CARDS *********************************************'
  write(u,'(a,f5.0,a)') '  1 SO ', mcnp6_radcm, '     $ sphere or radius[cm] centered at origin'
  write(u,'(a)') ''
  write(u,'(a)') 'c ******** BLOCK 3: DATA CARDS ************************************************ '
  write(u,'(a)') 'c ======== MATERIAL CARDS ========'
  write(u,'(a)') 'M0 nlib=70c    $ default data library'
  write(u,'(a)') 'c '
  write(u,'(a)') 'READ file=zaid_atomfrac_'//trim(adjustl(model_name))//'.mcnp6'
  write(u,'(a)') 'c '
  write(u,'(a)') 'M91    19039 1.0'
  write(u,'(a)') 'M92    12024 1.0'
  write(u,'(a)') 'c '
  write(u,'(a)') 'c ======== PHYSICS CARDS ========'
  write(u,'(a)') 'PHYS:N 9. 0 0 J J J 0 99.     $ neutron transport, max E 9 MeV, no models'
  write(u,'(a)') 'c '
  write(u,'(a)') 'c ======== SOURCE CARDS ========'
  write(u,'(a)') 'SDEF  pos=0 0 0  erg=d1       $ position, energy'
  if (it == ntargets+1) then
      write(u,'(a)') 'READ file=sdef_si_'//trim(adjustl(chains(ic)))//'_SF.mcnp6'
      write(u,'(a)') 'READ file=sdef_sp_'//trim(adjustl(chains(ic)))//'_SF.mcnp6'
  else
      write(u,'(a)') 'READ file=sdef_si_'//trim(adjustl(chains(ic)))//'_'//trim(adjustl(targets(it)))//'.mcnp6'
      write(u,'(a)') 'READ file=sdef_sp_'//trim(adjustl(chains(ic)))//'_'//trim(adjustl(targets(it)))//'.mcnp6'
  end if
  write(u,'(a)') 'c '
  write(u,'(a)') 'c ======== TALLIES ========'
  write(u,'(a)') 'F2:n 1'
  write(u,'(a)') 'Fc2 >>>>>>>>>>>>> TALLY 2: FLUX OF NEUTRONS LEAVING CELL 1 <<<<<<<<<<<<<'
  write(u,'(a)') 'Sd2 1'
  write(u,'(a)') 'c '
  write(u,'(a)') 'F14:n 1'
  write(u,'(a)') 'Fc14 >>>>>>>>>>>>> TALLY 14: # OF PROTONS FROM 39K(n,p) ONLY <<<<<<<<<<<<<'
  write(u,'(a,es13.7,a)') 'Fm14 ', K39atomcm3, ' 91 103     $ to account for atom density of 39K'
  write(u,'(a,f12.10,a)') 'Sd14 ', Kwtfrac, '             $ to get 39Ar per wt.frac. K'
  write(u,'(a)') 'c '
  write(u,'(a)') 'F24:n 1'
  write(u,'(a)') 'Fc24 >>>>>>>>>>>>> TALLY 24: # OF ALPHAS FROM 24Mg(n,a) ONLY <<<<<<<<<<<<<'
  write(u,'(a,es13.7,a)') 'Fm24 ', Mg24atomcm3, ' 92 107     $ to account for atom density of 24Mg'
  write(u,'(a,f12.10,a)') 'Sd24 ', Mgwtfrac, '             $ to get 21Ne per wt.frac. Mg'
  write(u,'(a)') 'c '
  write(u,'(a)') 'c ======== RELATED TO VARRED ========'
  write(u,'(a)') 'IMP:n 1 0                     $ neutron importance in cells'
  write(u,'(a)') 'c '
  write(u,'(a)') 'c ======== TERMINATION, OUTPUT, OTHER ========'
  write(u,'(a,f5.3,a)') 'STOP f14=', mcnp6_accu, ' nps=1000000    $ termination criteria'
  write(u,'(a)') 'PRINT                         $ print everything about the calculation'
  close (u)
END SUBROUTINE write_mcnp6_input


!******************************************************************************
!> Gets inputs from inputfile
!!
SUBROUTINE get_inputs
  USE prepare
  IMPLICIT NONE
  CHARACTER :: argv*127, command*255
  INTEGER :: len
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
  USE prepare
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
!> Gets decays chains from input file and allocates `chains(:)` array
!!
SUBROUTINE get_requested_decay_chains
  USE prepare
  IMPLICIT NONE
  INTEGER :: unit,i
  CHARACTER :: filename*255, nucid*5
  unit = 22
  filename = trim(input_directory)//'/'//trim(file_decay_chains)
  !< get number of entries (=lines)
  call open_file(unit, trim(filename), 'old', 'read')
  nchains = 0
  do
      read(unit, *, end=99)
      nchains = nchains + 1
  end do
99 continue
  if (nchains <= 0) STOP 'number of decay chains must be > 0 !'
  allocate (chains(nchains))
  !< get the entries (NUCIDs of top parent nuclide)
  rewind (unit)
  write (*,*) 'Number of decay chains =',nchains
  do i=1,nchains
      read(unit,'(a)') nucid
      call validate_NUCID_format(nucid) 
      chains(i) = nucid
      write (*,*) chains(i)
  end do
  write (*,*)
  close (unit)
END SUBROUTINE get_requested_decay_chains


!******************************************************************************
!> Gets (alpha,n) targets from input file and allocates `targets(:)` array
!!
SUBROUTINE get_requested_an_targets
  USE prepare
  IMPLICIT NONE
  INTEGER :: unit,i
  CHARACTER :: filename*255, nucid*5
  unit = 34
  filename = trim(input_directory)//'/'//trim(file_an_targets)
  !< get number of entries (=lines)
  call open_file(unit, trim(filename), 'old', 'read')
  ntargets = 0
  do
      read(unit, *, end=99)
      ntargets = ntargets + 1
  end do
99 continue
  if (ntargets < 0) STOP 'number of (alpha,n) targets must be >= 0 !'
  if (ntargets == 0) then
      write (*,*) 'WARNING: no (alpha,n) targets specified'
  end if
  allocate (targets(ntargets))
  !< get the entries (NUCIDs of top parent nuclide)
  rewind (unit)
  write (*,*) 'Number of (alpha,n) targets =',ntargets
  do i=1,ntargets
      read(unit,'(a5,1x,a255)') nucid
      call validate_NUCID_format(nucid)
      targets(i) = nucid
      write (*,*) targets(i)
  end do
  write (*,*)
  close (unit)
END SUBROUTINE get_requested_an_targets


!******************************************************************************
!> Checks NUCID format. Fixes positioning, or terminates if invalid.
!!
SUBROUTINE validate_NUCID_format(nucid)
  USE prepare
  IMPLICIT NONE
  CHARACTER(5), INTENT(INOUT) :: nucid !< NUCID :)
  CHARACTER(5) :: tmp
  INTEGER :: i,j,ver,plt,pbl
  tmp = adjustl(nucid)
  call upcase(tmp) !< just in case
  ver = verify(tmp,'ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890 ')
  if (ver /= 0) goto 99
  plt = verify(tmp,'1234567890')
  if (plt > 4 .or. plt < 2) goto 99
  pbl = plt - 1 + verify(tmp(plt:),'ABCDEFGHIJKLMNOPQRSTUVWXYZ')
  if (pbl == plt-1) pbl = 6
  if (pbl > plt+2) goto 99
  do i=pbl,5
      if (tmp(i:i) /= ' ') goto 99
  end do
  nucid = '     '
  nucid(4:5) = tmp(plt:plt+1)
  do i=plt-1,1,-1
      j = i-plt+4
      nucid(j:j) = tmp(i:i)
  end do
  RETURN
99 write (*,*) "ERROR: invalid NUCID '",tmp,"'"
STOP 'termination.'
END SUBROUTINE validate_NUCID_format


!******************************************************************************
!> Converts all alphanumeric characters in the string to upper case
!!
SUBROUTINE upcase(string)
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT) :: string
  INTEGER :: i,ic,len
  len = len_trim(string)
  do i=1,len
      ic = ichar(string(i:i))
      if (ic >= 97 .AND. ic <= 122) string(i:i) = char(ic-32)
  end do
END SUBROUTINE upcase


