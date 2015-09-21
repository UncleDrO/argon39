!> \file
!! Sets up decay chains 


!******************************************************************************
!> Top subroutine to set up decay chains
!!
SUBROUTINE setup_decay_chains
  USE all_data
  IMPLICIT NONE
  INTEGER :: ic

  call get_requested_decay_chains

  do ic=1,ndecc
      call setup_decay_chain(ic)
      call output_decay_chain(ic)
  end do
  minEalpha = minval(dchain(:)%minEa)
  maxEalpha = maxval(dchain(:)%maxEa)
  write (outp,'(/,a,f10.2,a)') 'Overall min E alpha: ', minEalpha*1.e3, ' keV'
  write (outp,'(a,f10.2,a,/)') 'Overall max E alpha: ', maxEalpha*1.e3, ' keV'

  write (outp,'(/,a)') "ALPHA YIELD"
  write (outp,111)
  do ic=1,ndecc
      write (outp,222) dchain(ic)%tpnucid, dchain(ic)%nalpha, &
          dchain(ic)%nalpha * dchain(ic)%num2rate * dchain(ic)%tpwtf
  end do
  write (outp,*)

  do ic=1,ndecc
      call output_alpha_E_spectra(ic)
  end do
  write (outp,*)

111 format (4x,"chain",7x,"per chain",4x,"per sec per kg")
222 format (4x,a,i12,8x,es12.3)
END SUBROUTINE setup_decay_chains


!******************************************************************************
!> Gets decays chains from input file and allocates `dchain(:)` array
!!
SUBROUTINE get_requested_decay_chains
  USE all_data
  IMPLICIT NONE
  INTEGER :: unit,i
  CHARACTER :: filename*255, nucid*5

  unit = 22
  filename = trim(input_directory)//'/'//trim(file_decay_chains)

  !< get number of entries (=lines)
  call open_file(unit, trim(filename), 'old', 'read')
  ndecc = 0
  do
      read(unit, *, end=99)
      ndecc = ndecc + 1
  end do
99 continue

  if (ndecc <= 0) STOP 'number of decay chains must be > 0 !'
  
  allocate (dchain(ndecc))

  !< get the entries (NUCIDs of top parent nuclide)
  rewind (unit)
  do i=1,ndecc
      read(unit,'(a)') nucid
      call validate_NUCID_format(nucid) 
      dchain(i)%tpnucid = nucid
      allocate (dchain(i)%Abranch(dchain(i)%nAbr)) !< initial allocation
      allocate (dchain(i)%Bbranch(dchain(i)%nBbr)) !< initial allocation
  end do
  close (unit)

  !< output info
  write (*,*) 'Number of decay chains =',ndecc
  write (outp,'(/,a,i2)') 'Number of decay chains =',ndecc
  do i=1,ndecc
      write (outp,'(4x,a)') dchain(i)%tpnucid
  end do
  write (outp,*)

END SUBROUTINE get_requested_decay_chains


!******************************************************************************
!> Sets up one decay chain.
!! Reads files `NUCID.chain` to get nuclides in the chain, 
!! then reads ENSDF decay data files for each nuclide to get decay data.
!! @note Could be made more automated: `get_requested_decay_chains` could read
!! top parent and stable daughter nuclides for each chain, then 
!! `setup_decay_chain` would go from the stable daughter up using ENSDF info
!! (rather than reading in the `NUCID.chain` file with all nuclides involved).
!!
SUBROUTINE setup_decay_chain(ic)
  USE all_data
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ic
  CHARACTER :: filename*255, nucid*5, ensdf*255
  INTEGER :: unit,idec,Z,A,inuc,i,iel
  TYPE(abra) :: adec
  TYPE(bbra) :: bdec
  REAL(DP) :: lambda,isotfrac,molarmass

  write (*,'(" Setting up chain ",a)') dchain(ic)%tpnucid
  write (outp,'(/"Setting up chain ",a)') dchain(ic)%tpnucid

  unit = 23
  filename = trim(input_directory)//'/'//trim(adjustl(dchain(ic)%tpnucid))//'.chain'
  call open_file(unit, trim(filename), 'old', 'read')
  read (unit, *) nucid !< skip top parent
  do
      read (unit, *, end=99) nucid
      ensdf = trim(input_directory)//'/'//trim(adjustl(nucid))//'.ensdf'
      !write (*   ,'(5x,"Reading dataset",1x,a)') trim(ensdf)
      write (outp,'(2x,"Reading dataset",1x,a)') trim(ensdf)
      call read_ENSDF_decay_file(trim(ensdf), idec, adec, bdec)
      select case (idec)
      case (0)
          call append_alpha_decay_to_chain(adec,ic)
      case (1)
          call append_beta_decay_to_chain(bdec,ic)
      case (2)
          call append_alpha_decay_to_chain(adec,ic)
          call append_beta_decay_to_chain(bdec,ic)
      case (-1)
          STOP 'no decay??? (in setup_decay_chain)'
          continue !< nothing to do
      end select
      call calculate_total_branching(ic,nucid)
  enddo
99 close (unit)
  write (outp,*)

  !< get number of alpha's emitted
  call NUCID_to_AZ(dchain(ic)%tpnucid,A,Z)
  dchain(ic)%nalpha = A
  call NUCID_to_AZ(dchain(ic)%Abranch(dchain(ic)%nAbr)%dnucid,A,Z)
  dchain(ic)%nalpha = (dchain(ic)%nalpha - A)/4

  !< calcualte "number-to-rate" conversion factor (# -> # per sec per kg)
  lambda = log(2._dp)/dchain(ic)%Abranch(1)%t12
  call NUCID_to_AZ(dchain(ic)%tpnucid,A,Z)
  call find_nuclide_index(A,Z,inuc)
  !< find abundance (wt.frac)
  iel = -1
  do i=1,nabund
      if (abund(i)%Z == Z) then
          iel = i
          exit
      end if
  end do
  if (iel == -1) STOP 'ERROR: Abundance of decay chain parent not found.'
  dchain(ic)%tpwtf = abund(iel)%wtf
  isotfrac = nuclide(inuc)%nab
  molarmass = nuclide(inuc)%wt 
  dchain(ic)%num2rate = lambda * isotfrac * Avogadro / molarmass * 1.e3_dp

END SUBROUTINE setup_decay_chain


!******************************************************************************
!> Appends alpha branch `adec` to `dchain(ic)%Abranch` array
!!
SUBROUTINE append_alpha_decay_to_chain(adec,ic)
  USE all_data
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ic
  TYPE(abra), INTENT(IN) :: adec
  TYPE(abra), ALLOCATABLE :: tmp(:)
  INTEGER :: na
  na = dchain(ic)%nAbr
  allocate (tmp(na))
  tmp(1:na) = dchain(ic)%Abranch(1:na) !< copy Abranch array to tmp
  deallocate (dchain(ic)%Abranch)
  allocate (dchain(ic)%Abranch(na+1))   !< expand Abranch array by one element
  dchain(ic)%Abranch(1:na) = tmp(1:na) !< copy back from tmp
  deallocate (tmp)
  na = na + 1
  dchain(ic)%Abranch(na) = adec        !< add new alpha branch data
  dchain(ic)%nAbr = na
END SUBROUTINE append_alpha_decay_to_chain


!******************************************************************************
!> Appends beta branch `bdec` to `dchain(ic)%Bbranch` array
!!
SUBROUTINE append_beta_decay_to_chain(bdec,ic)
  USE all_data
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ic
  TYPE(bbra), INTENT(IN) :: bdec
  TYPE(bbra), ALLOCATABLE :: tmp(:)
  INTEGER :: nb
  nb = dchain(ic)%nBbr
  allocate (tmp(nb))
  tmp(1:nb) = dchain(ic)%Bbranch(1:nb) !< copy Bbranch array to tmp
  deallocate (dchain(ic)%Bbranch)
  allocate (dchain(ic)%Bbranch(nb+1))   !< expand Abranch array by one element
  dchain(ic)%Bbranch(1:nb) = tmp(1:nb) !< copy back from tmp
  deallocate (tmp)
  nb = nb + 1
  dchain(ic)%Bbranch(nb) = bdec        !< add new alpha branch data
  dchain(ic)%nBbr = nb
END SUBROUTINE append_beta_decay_to_chain


!******************************************************************************
!> Calculates the total, including up-the-chain, branching of each decay branch
!!
!! @atention Total branching ratio of 234PA decay fixed to 1
!!
SUBROUTINE calculate_total_branching(ic,dnucid)
  USE all_data
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ic
  CHARACTER(5), INTENT(IN) :: dnucid
  CHARACTER(5) :: dnucid2
  LOGICAL :: found
  INTEGER :: ia,ib,ia2,ib2
  REAL(DP) :: lbr,ubr,tbr

  call find_Adecay_to_daughter(ic,dnucid,found,ia)
  if (found) then
      lbr = dchain(ic)%Abranch(ia)%br
      dnucid2 = dchain(ic)%Abranch(ia)%pnucid
      ubr = zero
      call find_Adecay_to_daughter(ic,dnucid2,found,ia2)
      if (found) ubr = ubr + dchain(ic)%Abranch(ia2)%tbr
      call find_Bdecay_to_daughter(ic,dnucid2,found,ib2)
      if (found) ubr = ubr + dchain(ic)%Bbranch(ib2)%tbr
      if (dnucid2 == dchain(ic)%tpnucid) ubr = one !< first alpha decay: tbr=1
      tbr = lbr*ubr
      if (tbr > one .or. tbr < zero) STOP 'total branching out of bounds???'
      dchain(ic)%Abranch(ia)%tbr = tbr
  end if

  call find_Bdecay_to_daughter(ic,dnucid,found,ib)
  if (found) then
      lbr = dchain(ic)%Bbranch(ib)%br
      dnucid2 = dchain(ic)%Bbranch(ib)%pnucid
      ubr = zero
      call find_Adecay_to_daughter(ic,dnucid2,found,ia2)
      if (found) ubr = ubr + dchain(ic)%Abranch(ia2)%tbr
      call find_Bdecay_to_daughter(ic,dnucid2,found,ib2)
      if (found) ubr = ubr + dchain(ic)%Bbranch(ib2)%tbr
      tbr = lbr*ubr
      if (dnucid == '234U ') tbr = one !< fix to 234PA decay!
      if (tbr > one .or. tbr < zero) STOP 'total branching out of bounds???'
      dchain(ic)%Bbranch(ib)%tbr = tbr
  end if

END SUBROUTINE calculate_total_branching


!******************************************************************************
!> Finds alpha decay in decay chain `ic` with daughter `dnucid`
!!
SUBROUTINE find_Adecay_to_daughter(ic,dnucid,found,ia)
  USE all_data
  IMPLICIT NONE
  INTEGER, INTENT(IN)      :: ic
  CHARACTER(5), INTENT(IN) :: dnucid
  LOGICAL, INTENT(OUT)     :: found
  INTEGER, INTENT(OUT)     :: ia
  INTEGER :: i
  found = .false.
  do i=1,dchain(ic)%nAbr
      if (dchain(ic)%Abranch(i)%dnucid == dnucid) then
          ia = i
          found = .true.
      end if
  end do
END SUBROUTINE find_Adecay_to_daughter


!******************************************************************************
!> Finds beta decay in decay chain `ic` with daughter `dnucid`
!!
SUBROUTINE find_Bdecay_to_daughter(ic,dnucid,found,ib)
  USE all_data
  IMPLICIT NONE
  INTEGER, INTENT(IN)      :: ic
  CHARACTER(5), INTENT(IN) :: dnucid
  LOGICAL, INTENT(OUT)     :: found
  INTEGER, INTENT(OUT)     :: ib
  INTEGER :: i
  found = .false.
  do i=1,dchain(ic)%nBbr
      if (dchain(ic)%Bbranch(i)%dnucid == dnucid) then
          ib = i
          found = .true.
      end if
  end do
END SUBROUTINE find_Bdecay_to_daughter


!******************************************************************************
!> Outputs parameters of decay chain `ic` to output file.
!! @attention 219AT B- decay not available from NNDC (branching 0.03)
!! @attention 215PO B- decay not available from NNDC (but branching only 2.3e-6)
SUBROUTINE output_decay_chain(ic)
  USE all_data
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ic
  INTEGER :: ia, nlev, ilev, ib
  REAL(DP) :: Elev

  write (outp,'(a,x,a,/)') dchain(ic)%tpnucid,'DECAY CHAIN'

  !< output alpha decay branches
  write (outp,'(2x,a,i3)') 'Number of ALPHA branches:',dchain(ic)%nAbr
  do ia=1,dchain(ic)%nAbr
      nlev = dchain(ic)%Abranch(ia)%nlev
      write (outp,111) &
          dchain(ic)%Abranch(ia)%pnucid, &
          dchain(ic)%Abranch(ia)%dnucid
      write (outp,222) &
          dchain(ic)%Abranch(ia)%t12, &
          dchain(ic)%Abranch(ia)%Qval * 1.d3, &
          dchain(ic)%Abranch(ia)%br, &
          dchain(ic)%Abranch(ia)%tbr
      write (outp,'(8x,"levels =",i3," :")') nlev
      do ilev=1,nlev
          Elev = dchain(ic)%Abranch(ia)%alevel(ilev)%E
          if (Elev < dchain(ic)%minEa) dchain(ic)%minEa = Elev
          if (Elev > dchain(ic)%maxEa) dchain(ic)%maxEa = Elev
          write (outp,333) &
              Elev * 1.d3, &
              dchain(ic)%Abranch(ia)%alevel(ilev)%I
      end do
  end do
  write (outp,'(//,2x,"Min alpha energy:",f10.2," keV")') dchain(ic)%minEa*1.e3
  write (outp,'(2x,"Max alpha energy:",f10.2," keV",//)') dchain(ic)%maxEa*1.e3

  !< output alpha decay branches
  write (outp,'(2x,a,i3)') 'Number of BETA  branches:',dchain(ic)%nBbr
  do ib=1,dchain(ic)%nBbr
      write (outp,111) &
          dchain(ic)%Bbranch(ib)%pnucid, &
          dchain(ic)%Bbranch(ib)%dnucid
      write (outp,222) &
          dchain(ic)%Bbranch(ib)%t12, &
          dchain(ic)%Bbranch(ib)%Qval * 1.d3, &
          dchain(ic)%Bbranch(ib)%br, &
          dchain(ic)%Bbranch(ib)%tbr
  end do
  write (outp,*)

111 format (4x,a," to ",a)
222 format (8x,"t12 =",es10.3," sec  Qval =",f8.2," keV  branch =",f11.8,"/",f10.8)
333 format (12x,f7.2,4x,f11.8)
END SUBROUTINE output_decay_chain


!******************************************************************************
!> Outputs energy spectra of alphas produced in chain `ic`
!!
SUBROUTINE output_alpha_E_spectra(ic)
  USE all_data
  USE mhunt
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ic
  REAL(DP) :: Elev,Aint,Asum,Esum
  INTEGER :: ia,nlev,ilev,unit,unit2
  CHARACTER :: filename*255

  filename = 'alpha_energy_spectrum_'//trim(adjustl(dchain(ic)%tpnucid))//'.out'
  write (outp,'("<<OUTFILE>> Alpha energy spectrum written to file: ",a)') &
      trim(filename)
  unit = 32
  call open_file(unit, trim(filename), 'unknown', 'write')

  filename = 'alpha_energy_table_'//trim(adjustl(dchain(ic)%tpnucid))//'.txt'
  write (outp,'("<<OUTFILE>> Alpha energy table written to file: ",a)') &
      trim(filename)
  unit2 = 33
  call open_file(unit2, trim(filename), 'unknown', 'write')

  Asum = zero
  Esum = zero
  do ia=1,dchain(ic)%nAbr
      nlev = dchain(ic)%Abranch(ia)%nlev
      write (unit2,666) &
          dchain(ic)%Abranch(ia)%pnucid, &
          dchain(ic)%Abranch(ia)%dnucid, &
          dchain(ic)%Abranch(ia)%Qval * 1.d3, &
          dchain(ic)%Abranch(ia)%tbr
      do ilev=1,nlev
          Elev = dchain(ic)%Abranch(ia)%alevel(ilev)%E
          Aint = dchain(ic)%Abranch(ia)%tbr * dchain(ic)%Abranch(ia)%alevel(ilev)%I
          if (Aint > zero) &
              write (unit,'(f8.4,es12.4)') Elev, Aint
          if (Aint > 1e-4) &
              write (unit2,777) Elev * 1.d3, &
              dchain(ic)%Abranch(ia)%alevel(ilev)%I

          Asum = Asum + Aint
          Esum = Esum + Elev*Aint
      enddo
  end do

  close (unit)
  close (unit2)

  write (outp,'(a,": Spectrum sums to",f6.3,", should be",i2)') &
      dchain(ic)%tpnucid, Asum, dchain(ic)%nalpha
  write (outp,'(a,": Mean alpha energy (MeV) =",f6.3,/)') &
      dchain(ic)%tpnucid, Esum/Asum

!666 format (a," to ",a,T17,"&",T20,"*",f7.2,T30,"&",T33,"*",f8.6)
!777 format (T17,"&",T21,f7.2,T30,"&",T34,f8.6)
666 format (a," to ",a,3x,"Qval = ",f7.2," MeV",3x,"tot.branch. = ",f8.6)
777 format (T20,"Ea = ",f7.2," MeV",3x,"int = ",f8.6)
END SUBROUTINE output_alpha_E_spectra
