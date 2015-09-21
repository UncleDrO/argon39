!> \file
!! MCNP6extract program

!******************************************************************************
!> Module for input/output paramters pulled from ARGON39 code
!!
MODULE extract
  USE inpout
  INTEGER :: nchains, ntargets
  CHARACTER(5), ALLOCATABLE :: chains(:), targets(:)
CONTAINS
  CHARACTER(len(trim(adjustl(string)))) FUNCTION trimal(string)
    CHARACTER(len=*), INTENT(IN) :: string
    trimal = trim(adjustl(string))
  END FUNCTION trimal
END MODULE extract


!******************************************************************************
!> This program extracts 39Ar values from MCNP6 output files
!!
!! Used in conjucntion with ARGON39 and MCNP6 codes.
!!
!! @author  Ondřej Šrámek, University of Maryland, sramek@umd.edu
!! @date    1 October 2013
!!
PROGRAM MCNP6prep
  USE extract
  IMPLICIT NONE
  INTEGER :: out0=3, out1=1, out2=2, ic, it
  CHARACTER :: filename0*127, filename1*127, filename2*127, stat*3, trgt*5
  REAL(DP) :: tal0, acc0, tal1, acc1, tal2, acc2
  REAL(DP), ALLOCATABLE :: tallies0(:,:,:), tallies1(:,:,:), tallies2(:,:,:)
  LOGICAL :: success

  call get_inputs

  call get_requested_decay_chains

  call get_requested_an_targets

  filename0 = 'MCNP6extract_'//trimal(model_name)//'_tally0.out'
  filename1 = 'MCNP6extract_'//trimal(model_name)//'_tally1.out'
  filename2 = 'MCNP6extract_'//trimal(model_name)//'_tally2.out'
  write (*,'(2a)') 'Output file name: ', trimal(filename0)
  write (*,'(2a)') 'Output file name: ', trimal(filename1)
  write (*,'(2a,/)') 'Output file name: ', trimal(filename2)

  write (*,*) '...Hit ENTER...'
  read (*,*)

  allocate (tallies0(ntargets+1,nchains,2), tallies1(ntargets+1,nchains,2), tallies2(ntargets+1,nchains,2))
  tallies0 = zero
  tallies1 = zero
  tallies2 = zero

  do it=1,ntargets+1
      do ic=1,nchains
          call read_from_mcnp6_output(it,ic,tal0,acc0,tal1,acc1,tal2,acc2,success)
          if (.not.success) then
              write (*,*) 'Problem getting the tallies...'
              read (*,*)
          end if
          tallies0(it,ic,1) = tal0
          tallies0(it,ic,2) = acc0
          tallies1(it,ic,1) = tal1
          tallies1(it,ic,2) = acc1
          tallies2(it,ic,1) = tal2
          tallies2(it,ic,2) = acc2
      end do
  end do

  call open_file(out0, filename0, 'unknown', 'write')
  call open_file(out1, filename1, 'unknown', 'write')
  call open_file(out2, filename2, 'unknown', 'write')
  write (out0,'(a,6a12)') 'target', chains, chains
  write (out1,'(a,6a12)') 'target', chains, chains
  write (out2,'(a,6a12)') 'target', chains, chains
  do it=1,ntargets+1
      stat = 'ok.'
      if (maxval(tallies1(it,:,2)) > mcnp6_accu) stat = '!!!'
      if (it==ntargets+1) then
          trgt = 'SF   '
      else
          trgt = targets(it)
      end if
      write (out0,222) trgt, tallies0(it,:,1), tallies0(it,:,2)
      write (out1,111) trgt, tallies1(it,:,1), tallies1(it,:,2), stat
      write (out2,222) trgt, tallies2(it,:,1), tallies2(it,:,2), maxval(tallies2(it,:,2))
  end do
  close (out0)
  close (out1)
  close (out2)
111 format (a5,4x,6es12.4,3x,a3)
222 format (a5,4x,6es12.4,f8.4)

  deallocate (tallies0, tallies1, tallies2)

  write (*,'(/,2a)') 'Output file name: ', trimal(filename0)
  write (*,'(2a)')   'Output file name: ', trimal(filename1)
  write (*,'(2a,/)') 'Output file name: ', trimal(filename2)

END PROGRAM MCNP6prep


!******************************************************************************
!> Writes MCNP6 input file
!!
SUBROUTINE read_from_mcnp6_output(it, ic, tally0, accu0, tally1, accu1, tally2, accu2, success)
  USE extract
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: it,ic
  REAL(DP), INTENT(OUT) :: tally0, accu0, tally1, accu1, tally2, accu2
  LOGICAL, INTENT(OUT) :: success
  CHARACTER :: mcnp6file*127, line*160
  INTEGER :: u, count, indx, i

  tally0 = -99.999_dp
  accu0  = -99.999_dp
  tally1 = -99.999_dp
  accu1  = -99.999_dp
  tally2 = -99.999_dp
  accu2  = -99.999_dp

  if (it == ntargets+1) then
      mcnp6file = 'mcnp6_'//trimal(model_name)//'_'//trimal(chains(ic))//'_SFo'
  else
      mcnp6file = 'mcnp6_'//trimal(model_name)//'_'//trimal(chains(ic))//'_'//&
          trimal(targets(it))//'o'
  end if
  write (*,'("Reading file ",a)') trimal(mcnp6file)
  u = 22
  call open_file(u, mcnp6file, 'old', 'read')

  !< WORK ON TALLY "0"
  count = 0
  do
      read (u,'(a)',end=99) line
      indx = index(line,'TALLY 2: FLUX OF NEUTRONS LEAVING CELL')
      if (indx > 0) count = count + 1
      if (count == 3) exit
  end do

  do i=1,16
      read (u,'(a)',end=99) line
      if (index(line,'surface  1') > 0) goto 11
  end do

99 close (u) !< SOME/ALL TALLY RESULTS NOT FOUND IN FILE :(
  success = .false.
  RETURN

11 read (u,*) tally0, accu0

  !< WORK ON TALLY "1"
  do
      read (u,'(a)',end=99) line
      if (index(line,'TALLY 14: # OF PROTONS FROM 39K(n,p) ONLY') > 0) exit
  end do

  do i=1,16
      read (u,'(a)',end=99) line
      if (index(line,'multiplier bin:') > 0) goto 22
  end do
  goto 99

22 read (u,*) tally1, accu1

  !< WORK ON TALLY "2"
  do
      read (u,'(a)',end=99) line
      if (index(line,'TALLY 24: # OF ALPHAS FROM 24Mg(n,a) ONLY') > 0) exit
  end do

  do i=1,16
      read (u,'(a)',end=99) line
      if (index(line,'multiplier bin:') > 0) goto 33
  end do
  goto 99

33 read (u,*) tally2, accu2
  close (u)

  success = .true.

END SUBROUTINE read_from_mcnp6_output


!******************************************************************************
!> Gets inputs from inputfile
!!
SUBROUTINE get_inputs
  USE extract
  IMPLICIT NONE
  CHARACTER :: argv*127, command*127
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
  USE extract
  IMPLICIT NONE
  CHARACTER(127), INTENT(OUT) :: argv    !< inputfile name (first argument)
  INTEGER, INTENT(OUT)        :: len     !< length of inputfile name
  CHARACTER(127), INTENT(OUT)  :: command !< program execution command
  INTEGER :: status, argc

  call get_command (command, len, status)
  if (status /= 0) then
      write (*,*) 'ERROR: get_command failed with status = ', status
      STOP
  end if
  write (*,'(/,2a,/)') 'command line = ', command (1:len)

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
      write (*,*) 'Problem opening file ',trim(adjustl(filename))
      write (*,*) 'Termination.'
      STOP
  end if
END SUBROUTINE open_file


!******************************************************************************
!> Gets decays chains from input file and allocates `chains(:)` array
!!
SUBROUTINE get_requested_decay_chains
  USE extract
  IMPLICIT NONE
  INTEGER :: unit,i
  CHARACTER :: filename*127, nucid*5
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
  write (*,'(a,i4)') 'Number of decay chains =',nchains
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
  USE extract
  IMPLICIT NONE
  INTEGER :: unit,i
  CHARACTER :: filename*127, nucid*5
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
  write (*,'(a,i4)') 'Number of (alpha,n) targets =',ntargets
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
  USE extract
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


