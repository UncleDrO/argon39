!> \file
!! ENSDF related


!******************************************************************************
!> Interface module for ENSDF related functions
!!
MODULE ensdf
  INTERFACE
      CHARACTER(5) FUNCTION AZ_to_NUCID(A,Z)
        USE all_data
        INTEGER, INTENT(IN) :: A,Z
      END FUNCTION AZ_to_NUCID
  END INTERFACE
END MODULE ensdf


!******************************************************************************
!> Returns NUCID from Z and A compliant to ENSDF standard, i.e., mass number 
!! in cols 1-3 right-justied, and element symbol in col 4-5 left-justied.
!!
CHARACTER(5) FUNCTION AZ_to_NUCID(A,Z)
  USE all_data
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: A !< mass (nucleon) number
  INTEGER, INTENT(IN) :: Z !< atomic (proton) number
  write(AZ_to_NUCID,'(i3,a2)') A,element(Z)%symb
END FUNCTION AZ_to_NUCID


!******************************************************************************
!> Gets A and Z from NUCID. Assumes ENSDF standard-compliant NUCID.
!!
SUBROUTINE NUCID_to_AZ(nucid,A,Z)
  USE all_data
  IMPLICIT NONE
  CHARACTER(5), INTENT(IN) :: nucid !< NUCID :)
  INTEGER, INTENT(OUT) :: A !< mass number
  INTEGER, INTENT(OUT) :: Z !< atomic number
  CHARACTER(2) :: symb
  A = -1
  Z = -1
  read (nucid(1:3),*) A
  read (nucid(4:5),*) symb
  do Z=1,nelem
      if (element(Z)%symb == symb) RETURN
  end do
  STOP 'ERROR: NUCID_to_AZ failed.'
END SUBROUTINE NUCID_to_AZ


!******************************************************************************
!> Checks NUCID format. Fixes positioning, or terminates if invalid.
!!
SUBROUTINE validate_NUCID_format(nucid)
  USE all_data
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
!> Reads ENSDF decay dataset.
!!
!! @attention Do I need to worry about 234PA IT DECAY?? (probably not...)
!!
SUBROUTINE read_ENSDF_decay_file(ensdf, idec, adec, bdec)
  USE inpout
  USE abra_type
  USE bbra_type
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN) :: ensdf
  INTEGER, INTENT(OUT) :: idec
  TYPE(abra), INTENT(OUT) :: adec
  TYPE(bbra), INTENT(OUT) :: bdec
  INTEGER :: unit,j,istat
  CHARACTER :: record*80, dnucid*5, dsid*30, pnucid*5, jpar*18
  REAL(DP) :: edl,t12,qval,br,eal,ial
  LOGICAL :: boolA, boolB

  boolA=.false.
  boolB=.false.
  allocate (adec%alevel(adec%nlev))
  unit = 24
  call open_file(unit, ensdf, 'old', 'read')
  !< look for Identification record
  read (unit, '(a80)', err=88) record
  !write (*,'(a)') record
  if (len_trim(record) == 0 .or. record(6:9) /= '    ') then
      write (*,*) 'ERROR(read_ENSDF_decay_file): ',&
          ensdf,': first record must be ID record!!!'
      STOP
  end if
  mrec: do
!      write (*,'(a)') record
      call ENSDF_ID_record(record,dnucid,dsid)
      if (index(dsid,'A DECAY') > 0 .or. index(dsid,'B- DECAY') > 0) then
          !< look for Parent or Normalization record
          do j=1,2
              do
                  read (unit, '(a80)', err=88) record
                  if (record(6:8) == '  P' .or. record(6:8) == '  N') exit
              end do
!              write (*,'(a)') record
              if (record(6:8) == '  P') then
                  call ENSDF_P_record(record,pnucid,edl,jpar,t12,qval)
              else if (record(6:8) == '  N') then
                  call ENSDF_N_record(record,br)
              end if
          end do
          !< look for A|B records
          if (index(dsid,'A DECAY') > 0) then !< look for Alpha records
              boolA = .true.
              adec%pnucid = pnucid
              adec%dnucid = dnucid
              adec%dsid = dsid
              adec%t12 = t12
              adec%Qval = qval * 1.e-3_dp !< keV to MeV
              adec%br = br
              AA: do
                  read (unit, '(a80)', err=88) record
                  if (record(6:8) == '  A') then
!                      write (*,'(a)') record
                      call ENSDF_A_record(record,eal,ial,istat)
                      if (istat==0) then
                          call append_level_to_abranch(eal,ial,adec)
                      else
                          write (outp,'(6x,a)') 'WARNING: Funny A record:'
                          write (outp,'(10x,a)') record
                      end if
                  else if (len_trim(record) == 0) then !< END record
                      exit AA
                  end if
              end do AA
          else if (index(dsid,'B- DECAY') > 0) then !< look for Beta records
              boolB = .true.
              bdec%pnucid = pnucid
              bdec%dnucid = dnucid
              bdec%dsid = dsid
              bdec%t12 = t12
              bdec%Qval = qval * 1.e-3_dp !< keV to MeV
              bdec%br = br
              BB: do
                  read (unit, '(a80)', err=88) record
                  if (record(7:8) == ' B') then
!                      write (*,'(a)') record
                      continue !< do nothing
                  else if (len_trim(record) == 0) then !< END record
                      exit BB
                  end if
              end do BB
          else
              write (*,*) 'Should not be here...'
              STOP
          end if
      else
          XX: do
              read (unit, '(a80)', err=88) record
              if (len_trim(record) == 0) exit XX !< END record
          end do XX
      end if
!      write (*,'(a)') record
      do
          read (unit, '(a80)', end=99) record
!          write (*,'(a)') record
          if (record(6:9) == '    ') cycle mrec !< another data set (ID record)
      end do
  end do mrec

99 close (unit)
  !< set idec
  if (boolA .and. boolB) then
      idec = 2
  else if (boolA) then
      idec = 0
  else if (boolB) then
      idec = 1
  else
      idec = -1
  end if
  RETURN

88 write (*,*) 'ERROR(read_ENSDF_decay_file): trouble reading records in',ensdf
  STOP

END SUBROUTINE read_ENSDF_decay_file


!******************************************************************************
!> Reads ENSDF Identification record. 
!!
!! Output names correspond to ENSDF standard.
!! @note Not all defined fields are being read and processed.
!!
SUBROUTINE ENSDF_ID_record(rec,nucid,dsid)
  IMPLICIT NONE
  CHARACTER(80), INTENT(IN) :: rec
  CHARACTER, INTENT(OUT) :: nucid*5, dsid*30 !< NUCID, DSID :)
  nucid(1:5) = rec(1:5)
  call upcase(nucid) !< just in case
  dsid(1:30) = rec(10:39)
END SUBROUTINE ENSDF_ID_record


!******************************************************************************
!> Reads ENSDF Parent record. 
!!
!! Output names correspond to ENSDF standard. Half-life output in seconds.
!! If stable isotope, t set to -1.
!! @note Not all defined fields are being read and processed.
!!
SUBROUTINE ENSDF_P_record(rec,nucid,e,j,t,qp)
  USE precision
  IMPLICIT NONE
  CHARACTER(80), INTENT(IN)  :: rec    !< ENSDF 80-char record
  CHARACTER(5), INTENT(OUT)  :: nucid  !< NUCID :)
  REAL(DP), INTENT(OUT)      :: e      !< energy of decaying level in keV
  CHARACTER(18), INTENT(OUT) :: j      !< parity
  REAL(DP), INTENT(OUT)      :: t      !< half-life in seconds
  REAL(DP), INTENT(OUT)      :: qp     !< gs-gs Qvalue in keV
  CHARACTER :: ce*10, cj*18, ct*10, cqp*10, unit*2
  e = zero
  t = zero
  qp = zero
  nucid(1:5) = rec(1:5)
  ce(1:10) = rec(10:19)
  cj(1:18) = rec(22:39)
  ct(1:10) = rec(40:49)
  cqp(1:10) = rec(65:74)
  call remove_pX_from_E_field(ce)
  read (ce,*) e
  j = adjustl(cj)
  read (cqp,*) qp
  if (trim(adjustl(ct)) == 'STABLE') then
      t = -1._dp
  else
      read (ct,*) t,unit
      select case (unit)
      case ('Y')
          t = t * secinyr
      case ('D')
          t = t * secinday
      case ('H')
          t = t * 3600._dp
      case ('M')
          t = t * 60._dp
      case ('S')
          continue
      case ('MS')
          t = t * 1.e-3_dp
      case ('US')
          t = t * 1.e-6_dp
      case ('NS')
          t = t * 1.e-9_dp
      case ('PS')
          t = t * 1.e-12_dp
      case ('FS')
          t = t * 1.e-15_dp
      case ('AS')
          t = t * 1.e-18_dp
      case default
          write (*,*) 'ENSDF_P_record: Something funky while reading half-life:'
          write (*,*) rec
          STOP 
      end select
  end if
END SUBROUTINE ENSDF_P_record


!******************************************************************************
!> Removes +X from E entries in ENSDF record
!!
!! @attention But what does +X mean in ENSDF energy values??
!!
SUBROUTINE remove_pX_from_E_field(ce)
  USE inpout
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT) :: ce
  INTEGER :: ipx
  ipx = index(ce,'+X')
  if (ipx > 0) then 
      write (outp,'(6x,a,x,a)') 'WARNING:',ce
      ce(ipx:ipx+1) = '  '
      !!write (*,*) '    ',ce
      !!read (*,*)
  end if
END SUBROUTINE remove_pX_from_E_field

!******************************************************************************
!> Reads ENSDF Normalization record. 
!!
!! Output names correspond to ENSDF standard.
!! @note Not all defined fields are being read and processed.
!!
SUBROUTINE ENSDF_N_record(rec,br)
  USE precision
  IMPLICIT NONE
  CHARACTER(80), INTENT(IN)  :: rec    !< ENSDF 80-char record
  REAL(DP), INTENT(OUT)      :: br     !< branching ratio
  CHARACTER :: cbr*8
  br = zero
  cbr(1:8) = rec(32:39)
  read (cbr,*) br
END SUBROUTINE ENSDF_N_record


!******************************************************************************
!> Reads ENSDF Alpha record. 
!!
!! Output names correspond to ENSDF standard.
!! @note Not all defined fields are being read and processed.
!!
SUBROUTINE ENSDF_A_record(rec,e,ia,istat)
  USE precision
  IMPLICIT NONE
  CHARACTER(80), INTENT(IN)  :: rec !< ENSDF 80-char record
  REAL(DP), INTENT(OUT)      :: e   !< alpha energy in keV
  REAL(DP), INTENT(OUT)      :: ia  !< intensity of level as % of total alpha
  !!REAL(DP), INTENT(OUT)      :: hf  !< hindrance factor
  INTEGER, INTENT(OUT)       :: istat !< 
  CHARACTER :: ce*10, cia*8 !!, chf*8
  e = zero
  ia = zero
  !!hf = zero
  ce(1:10) = rec(10:19)
  cia(1:8) = rec(22:29)
  !!chf(1:8) = rec(32:39)
  read (ce,*,end=99) e
  read (cia,*,end=99) ia
  !!read (chf,*,end=99) hf
  if (e > zero .and. ia > zero) then
      istat = 0
      RETURN
  end if
99 istat = -1
END SUBROUTINE ENSDF_A_record


!******************************************************************************
!> Appends alpha level to a variable of type abra (alpha decay branch)
!!
SUBROUTINE append_level_to_abranch(eal,ial,adec)
  USE abra_type
  IMPLICIT NONE
  REAL(DP), INTENT(IN)      :: eal  !< alpha energy in keV
  REAL(DP), INTENT(IN)      :: ial  !< intensity of level as % of total alpha
  TYPE(abra), INTENT(INOUT) :: adec !< alpha decay branch which is updated
  TYPE(alev), ALLOCATABLE :: tmp(:)
  INTEGER :: nl
  nl = adec%nlev
  allocate (tmp(nl))
  tmp(1:nl) = adec%alevel(1:nl) !< copy alevel array to tmp
  deallocate (adec%alevel)
  allocate (adec%alevel(nl+1))  !< expand alevel array by one element
  adec%alevel(1:nl) = tmp(1:nl) !< copy back from tmp
  deallocate (tmp)
  nl = nl + 1                   !< add new level data
  adec%alevel(nl)%E = eal * 0.001_dp !< keV to MeV
  adec%alevel(nl)%I = ial * 0.01_dp  !< % to frac
  adec%nlev = nl
END SUBROUTINE append_level_to_abranch


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


