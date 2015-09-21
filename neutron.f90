!> \file
!! Neutron production: (a,n)


!******************************************************************************
!> Outputs possible (alpha,n) targets, where "possible" means a particular 
!! target nuclide is present.
!!
SUBROUTINE output_possible_targets
  USE all_data
  IMPLICIT NONE
  INTEGER :: i, j, j4, Z, A, Z1, A1
  REAL(DP) :: m1, m2, m3, m4, nab2, atf2, Qval, Ethresh, Bcoul, third
  CHARACTER(5) :: nucid2, nucid4
  write (outp,'(/,a)') 'Possible targets (i.e., nuclide is present):'
  write (outp,'(x,a,6x,a,6x,3(a,3x),a,4x,a,12x,a,12x,a,6x,a,7x,a,8x,a,9x,a)') &
      '2', '4', 'Z2', 'A2', 'Z4', 'A4', 'm2', 'm4', 'Qval', 'Eth', 'Vc', &
      'nab2', 'atf2'
  call find_nuclide_index(4,2,j)
  m1 = nuclide(j)%wt * u_in_MeV
  m3 = mass_neutron * u_in_MeV
  do i=1,nabund
      Z = abund(i)%Z
      do j=1,nnucl
          nab2 = nuclide(j)%nab
          if (nuclide(j)%Z == Z .and. nab2 > zero) then
              nucid2 = nuclide(j)%nucid
              atf2 = abund(i)%atf * nab2
              m2 = nuclide(j)%wt * u_in_MeV
              A = nuclide(j)%A
              call find_nuclide_index(A+3,Z+2,j4)
              m4 = nuclide(j4)%wt * u_in_MeV
              nucid4 = nuclide(j4)%nucid
              Qval = m1 + m2 - m3 - m4
              if (Qval < zero) then
                  Ethresh = -(m1+m2)/m2*Qval
              else
                  Ethresh = zero
              end if
              Z1 = 2 !< Z of alpha
              A1 = 4 !< A of alpha
              third = 1._dp/3._dp
              Bcoul = 1.44 * Z * Z1 / (1.2 * (A**third + A1**third))
              m2 = m2/u_in_MeV
              m4 = m4/u_in_MeV
              write (outp,111) nucid2, nucid4, Z, A, Z+2, A+3, m2, m4, Qval, &
                  Ethresh, Bcoul, nab2, atf2
          end if
      end do
  end do
111 format (a,2x,a,4i5,2f14.8,3f10.4,f14.8,es12.3)
END SUBROUTINE output_possible_targets


!******************************************************************************
!> Gets (alpha,n) targets from input file and allocates `targets(:)` array
!!
SUBROUTINE get_requested_an_targets(ntargets)
  USE all_data
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: ntargets
  INTEGER :: unit,i,A,Z
  CHARACTER :: filename*255, nucid*5, postfix*255

  unit = 34
  filename = trim(input_directory)//'/'//trim(file_an_targets)

  !< get number of entries (=lines)
  call open_file(unit, trim(filename), 'old', 'read')
  ntgts = 0
  do
      read(unit, *, end=99)
      ntgts = ntgts + 1
  end do
99 continue

  if (ntgts < 0) STOP 'number of (alpha,n) targets must be >= 0 !'
  if (ntgts == 0) then
      write (*,*) 'WARNING: no (alpha,n) targets specified'
      write (outp,'(/,a,/)') 'WARNING: no (alpha,n) targets specified'
  end if
  ntargets = ntgts
  
  allocate (targets(ntgts), target_fnames(ntgts), results(ntgts+1,ndecc,4))
  results = zero

  !< get the entries (NUCIDs of top parent nuclide)
  rewind (unit)
  do i=1,ntgts
      read(unit,'(a5,1x,a255)') nucid, postfix
      call validate_NUCID_format(nucid)
      call NUCID_to_AZ(nucid,A,Z)
      targets(i)%nucid = nucid
      targets(i)%A = A
      targets(i)%Z = Z
      target_fnames(i) = trim(adjustl(nucid))//trim(adjustl(postfix))
  end do
  close (unit)

  !< output info
  write (*,*) 'Number of (alpha,n) targets =',ntgts
  write (outp,'(/,a,i2)') 'Number of (alpha,n) targets =',ntgts
  do i=1,ntgts
      write (outp,'(4x,a,2i6)') targets(i)%nucid, targets(i)%A, targets(i)%Z
  end do
  write (outp,*)

END SUBROUTINE get_requested_an_targets


!******************************************************************************
!> Sets up values and pointers to target nuclide properties
!!
SUBROUTINE setup_target_nuclide(it)
  USE all_data
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: it
  INTEGER :: Z,A,inuc,i,iel
  CHARACTER(5) :: nucid2,nucid4

  Z = targets(it)%Z
  A = targets(it)%A

  !< 1: alpha
  call find_nuclide_index(4,2,inuc)
  m1an = nuclide(inuc)%wt * u_in_MeV
  !< 2: target nuclide
  call find_nuclide_index(A,Z,inuc)
  target_nucl => nuclide(inuc)
  nucid2 = target_nucl%nucid
  m2an = target_nucl%wt * u_in_MeV
  !< 3: neutron
  m3an = mass_neutron * u_in_MeV
  !< 4: product nuclide
  A = target_nucl%A+3
  Z = target_nucl%Z+2
  call find_nuclide_index(A,Z,inuc)
  m4an = nuclide(inuc)%wt * u_in_MeV
  nucid4 = nuclide(inuc)%nucid

  QvalAN = m1an + m2an - m3an - m4an
  if (QvalAN < zero) then
      EthreshAN = -(m1an+m2an)/m2an*QvalAN
  else
      EthreshAN = zero
  end if

  !< target nuclide atomic density (in #/cm3)
  Z = target_nucl%Z
  iel = -1
  do i=1,nabund
      if (abund(i)%Z == Z) then
          iel = i
          exit
      end if
  end do
  if (iel == -1) then
      write (*,*) 'ERROR: Abundance of target element not found, Z =', Z
      STOP 'Termination.'
  end if
  target_nucl_atcm3 = abund(iel)%atcm3 * target_nucl%nab
  target_elem_wtf = abund(iel)%wtf

  write (*,*) trim(adjustl(nucid2)),'(a,n)',trim(adjustl(nucid4)),' reaction'
  write (outp,'(/,3a)') trim(adjustl(nucid2)),'(a,n)',trim(adjustl(nucid4))
  write (outp,'(4x,2a)') 'Target nuclide: ',nucid2
  write (outp,'(4x,"wt.frac =",f12.8,4x,"atoms/cm3 =",es11.3)') &
      abund(iel)%wtf * target_nucl%nab, target_nucl_atcm3
  write (outp,'(4x,"m1 =",f12.8,"u")') m1an/u_in_MeV
  write (outp,'(4x,"m2 =",f12.8,"u")') m2an/u_in_MeV
  write (outp,'(4x,"m3 =",f12.8,"u")') m3an/u_in_MeV
  write (outp,'(4x,"m4 =",f12.8,"u")') m4an/u_in_MeV
  write (outp,'(4x,"Q value =",f12.8," MeV")') QvalAN
  write (outp,'(4x,"Ethresh =",f12.8" MeV",/)') EthreshAN

END SUBROUTINE setup_target_nuclide


!******************************************************************************
!> Reads (a,n) cross section
!!
SUBROUTINE read_an_cross_section(it)
  USE all_data
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: it
  INTEGER :: unit, i, mini
  CHARACTER :: filename*255, line*80

  select case (which_an_calc)
  case (1)
      filename = trim(local_directory)//'/'//trim(adjustl(target_fnames(it)))
  case (2)
      filename = trim(talys_directory)//'/'//trim(adjustl(targets(it)%nucid))//'/nprod.tot'
  end select

  unit = 26
  call open_file(unit, trim(filename), 'old', 'read')
  write (*,*) '  Reading (a,n) cross section'
  write (outp,'(/,"Reading ",a,"(a,n) cross section from file: ",a)') &
      trim(adjustl(target_nucl%nucid)), trim(filename)

  !< get number of entries
  nEaCS = 0
  do 
      read (unit, '(a80)', end=99) line
      if (line(1:1) == '#') cycle
      nEaCS = nEaCS + 1
  end do
99 rewind (unit)

  !< read cross section
  allocate (CS(nEaCS), EaCS(nEaCS))
  i = 1
  do
      read (unit, '(a80)') line
      if (line(1:1) /= '#') then
          read (line,*) EaCS(i), CS(i)
          i = i+1
      end if
      if (i > nEaCS) exit
  end do
  close (unit)

  !< get Emin
  mini = 1
  do i=1,nEaCS
      if (CS(i) > zero) exit
      mini = i
  end do
  minE_CS = EaCS(mini)
  write (outp,'(4x,"Min E of non-zero cross section =",f12.8," MeV")') minE_CS

  !< output cross section
  write (outp,'(4x,a,4x,a)') 'E in MeV', 'XS in mbarn'
  do i=mini,nEaCS
      write (outp,'(2f12.3)') EaCS(i), CS(i)
      if (EaCS(i) > maxEalpha) exit
  end do
  write (outp,*)

  if (EaCS(nEaCS) < maxEalpha) then
      write (*,*) 'Max (a,n) xsec energy lower than maxEalpha!'
      write (*,*) EaCS(nEaCS), ' < ', maxEalpha
      STOP 'Termination.'
  end if

  if (stoppow(1,1) > minE_CS) &
      STOP 'Min stopping power energy higher than minE_CS!'

END SUBROUTINE read_an_cross_section


!==============================================================================
!--- CALCULATE YIELD ----------------------------------------------------------

!******************************************************************************
!> Top subroutine to calculate (a,n) neutron yield from each alpha decay branch of each chain
!!
SUBROUTINE calculate_an_neutron_yield(it)
  USE all_data
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: it
  INTERFACE
      REAL(DP) FUNCTION interpolate_tab_func(xa,fa,xi,type)
        USE precision
        USE mhunt
        REAL(DP), DIMENSION(:), INTENT(IN) :: xa,fa
        REAL(DP), INTENT(IN) :: xi
        INTEGER, INTENT(IN) :: type
      END FUNCTION interpolate_tab_func
  END INTERFACE
  REAL(DP) :: yield, Ekl, fkl, npf, tbr
  INTEGER :: ic,k,l

  write (*,*) '  Calculating neutron production function'
  call tabulate_neutron_prod_func

  call output_neutron_prod_func_table

  write (*,*) '  Calculating (a,n) neutron yields from decay chains'
  do ic=1,ndecc
      dchain(ic)%nyieldAN = zero
      do k=1,dchain(ic)%nAbr
          yield = zero
          do l=1,dchain(ic)%Abranch(k)%nlev
              Ekl = dchain(ic)%Abranch(k)%alevel(l)%E
              fkl = dchain(ic)%Abranch(k)%alevel(l)%I
              npf = interpolate_tab_func(Enpfc,NPFc,Ekl,0) !< lin-lin
              yield = yield + fkl*npf
          end do
          tbr = dchain(ic)%Abranch(k)%tbr
          dchain(ic)%Abranch(k)%nyieldAN = tbr * yield
          dchain(ic)%nyieldAN = dchain(ic)%nyieldAN + dchain(ic)%Abranch(k)%nyieldAN
      end do
  end do

  call output_neutron_yield(it)

END SUBROUTINE calculate_an_neutron_yield


!******************************************************************************
!> Tabulates neutron production function
!! 
SUBROUTINE tabulate_neutron_prod_func
  USE all_data
  IMPLICIT NONE
  INTERFACE
      REAL(DP) FUNCTION func_npf(E)
        USE all_data
        USE mhunt
        REAL(DP), INTENT(IN) :: E
      END FUNCTION func_npf
  END INTERFACE
  REAL(DP) :: minE, maxE, a, b, s
  INTEGER :: i

  minE = zero
  maxE = ceiling(maxEalpha/dEan)*dEan
  Nnpfc = nint((maxE-minE)/dEan)
  allocate (Enpfc(0:Nnpfc), NPFc(0:Nnpfc))
  Enpfc(0:Nnpfc) = (/ ( minE+i*dEan,i=0,Nnpfc) /)

  NPFc(0) = zero
  do i=1,Nnpfc
      a = Enpfc(i-1)
      b = Enpfc(i)
      if (a < minE_CS) then
          s = zero
      else
          call qtrap(func_npf,a,b,s)
      end if
      NPFc(i) = NPFc(i-1) + s
  end do

END SUBROUTINE tabulate_neutron_prod_func


!******************************************************************************
!> Function (a,n)_cross_section / linear_stopping_power provided for 
!! integration of neutron production function.
!!
REAL(DP) FUNCTION func_npf(E)
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
  REAL(DP) :: xsec, lspow
  !< get cross section: lin-lin interpolation
  xsec = interpolate_tab_func(EaCS(:),CS(:),E,0)
  !< get linear stopping power: log-log interpolation
  lspow = interpolate_tab_func(stoppow(:,1),stoppow(:,2),E,1) * rock_density
  !< integrand for neutron production function (including target atom density)
  func_npf = target_nucl_atcm3 * xsec*1.e-27_dp / lspow
END FUNCTION func_npf


!==============================================================================
!--- CALCULATE SPECTRA --------------------------------------------------------

!******************************************************************************
!> Top subroutine to calculate (a,n) neutron energy spectra from each decay chain
!!
SUBROUTINE calculate_an_neutron_spectrum
  USE all_data
  USE mhunt
  IMPLICIT NONE
  REAL(DP) :: tbr, Ekl, fkl
  REAL(DP), ALLOCATABLE :: nspec1(:)
  INTEGER :: ic, k, l

  write (*,*) '  Calculating (a,n) neutron enegy spectra'

  call tabulate_differential_neutron_prod_func

  !< calculate differential neutron yield
  allocate (nspec1(0:Nnspec))
  do ic=1,ndecc
      allocate (dchain(ic)%nspecAN(0:Nnspec))
      dchain(ic)%nspecAN(:) = zero
      do k=1,dchain(ic)%nAbr
          tbr = dchain(ic)%Abranch(k)%tbr
          do l=1,dchain(ic)%Abranch(k)%nlev
              Ekl = dchain(ic)%Abranch(k)%alevel(l)%E
              fkl = dchain(ic)%Abranch(k)%alevel(l)%I
              select case (which_an_calc)
              case (1)
                  call local_an_differential_npf_for_alpha_E(Ekl,nspec1)
              case (2)
                  call talys_an_differential_npf_for_alpha_E(Ekl,nspec1)
              end select
              nspec1(:) = nspec1(:) * tbr * fkl
              dchain(ic)%nspecAN(:) = dchain(ic)%nspecAN(:) + nspec1(:)
          end do
      end do
  end do

  call cleanup_check_neutron_spectra

  call output_neutron_spectra

  deallocate(nspec1, CS, EaCS, Enpfc, NPFc, Enspec)
  if (which_an_calc==2) deallocate (DCS, DNPFc)
  do ic=1,ndecc
      deallocate (dchain(ic)%nspecAN)
  end do

END SUBROUTINE calculate_an_neutron_spectrum


!******************************************************************************
!> Tabulates TALYS differential neutron production function
!!
!! @todo FIX integration of DNPF...
SUBROUTINE tabulate_differential_neutron_prod_func
  USE all_data
  IMPLICIT NONE
  INTERFACE
      REAL(DP) FUNCTION func_dnpf_nj(E)
        USE all_data
        USE mhunt
        REAL(DP), INTENT(IN) :: E
      END FUNCTION func_dnpf_nj
  END INTERFACE
  INTEGER :: i, j, k
  REAL(DP) :: minE, Efile, a, b, s
  REAL(DP), ALLOCATABLE :: DCStal(:)

  !< prepare neutron spectrum energy array
  minE = zero
  En_cutoff = maxEalpha + QvalAN + 1.0
  write (outp,'(a,a,f7.3,a,/)') trim(adjustl(target_nucl%nucid)), &
      '(a,n) energy cutoff for neutron spectra:', En_cutoff, ' MeV'
  flush (outp)
  En_cutoff = ceiling(En_cutoff/dEan)*dEan
  Nnspec = nint((En_cutoff-minE)/dEan)
  allocate (Enspec(0:Nnspec))
  Enspec(0:Nnspec) = (/ ( minE+i*dEan,i=0,Nnspec) /)

  if (which_an_calc==1) RETURN !< "local" calculation

  !< tabulate differential cross section array
  allocate (DCS(nEaCS,0:Nnspec), DCStal(0:Nnspec))
  DCS(:,:) = zero
  do i=1,nEaCS
      if (CS(i) == zero) cycle
      Efile = EaCS(i)
      call retrieve_interpolate_talys_spectrum(Efile, DCStal)
      DCS(i,:) = DCStal(:)
  end do
  deallocate (DCStal)

  !< tabulate differential neutron production function
  allocate (DNPFc(0:Nnpfc,0:Nnspec))
  DNPFc(0,:) = zero
  do j=1,Nnspec
      indx_dnpf_nj = j
      do i=1,Nnpfc
          a = Enpfc(i-1)
          b = Enpfc(i)
          if (a < minE_CS) then
              s = zero
          else
              !!call qtrap(func_dnpf_nj, a, b, s)
              do k=1,6 !! preset number of trapzd steps - FIX...
                  call trapzd(func_dnpf_nj, a, b, s, k)
              enddo
          end if
          DNPFc(i,j) = DNPFc(i-1,j) + s
      end do
  end do

END SUBROUTINE tabulate_differential_neutron_prod_func


!******************************************************************************
!> Function (a,n)_cross_section_spectrum / linear_stopping_power at a given
!! neutron energy provided for integration of diff. neutron production function.
!!
REAL(DP) FUNCTION func_dnpf_nj(E)
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
  REAL(DP) :: dxsec, lspow
  INTEGER :: j
  !< neutron energy index
  j = indx_dnpf_nj
  !< get differential cross section: lin-lin interpolation
  dxsec = interpolate_tab_func(EaCS(:), DCS(:,j), E, 0)
  !< get linear stopping power: log-log interpolation
  lspow = interpolate_tab_func(stoppow(:,1),stoppow(:,2), E, 1) * rock_density
  !< integrand for neutron production function (including target atom density)
  func_dnpf_nj = target_nucl_atcm3 * dxsec*1.e-27_dp / lspow
END FUNCTION func_dnpf_nj


!******************************************************************************
!> Retrieves TALYS differential neutron production cross section (n spectrum)
!! and interpolates to `Enspec` energy grid
!!
SUBROUTINE retrieve_interpolate_talys_spectrum(Efile, DCStal)
  USE all_data
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: Efile
  REAL(DP), INTENT(OUT) :: DCStal(0:Nnspec)
  INTERFACE
      REAL(DP) FUNCTION interpolate_tab_func(xa,fa,xi,type)
        USE precision
        USE mhunt
        REAL(DP), DIMENSION(:), INTENT(IN) :: xa,fa
        REAL(DP), INTENT(IN) :: xi
        INTEGER, INTENT(IN) :: type
      END FUNCTION interpolate_tab_func
  END INTERFACE
  CHARACTER :: filename*127, line*80
  INTEGER :: unit=55, ne, i
  REAL(DP) :: Ei, de
  REAL(DP), ALLOCATABLE :: etal(:), stal(:)

  write (filename, '(a,"/",a,"/nspec",i3.3,f4.3,".tot")') &
      trim(adjustl(talys_directory)), trim(adjustl(target_nucl%nucid)), &
      int(Efile), Efile-int(Efile)
  call open_file(unit, trim(filename), 'old', 'read')
  read (unit,*)
  read (unit,*)
  read (unit,*)
  read (unit,'(a80)') line
  read (line(15:80),*) ne
  read (unit,*)

  !< read in TALYS spectrum
  allocate (etal(0:ne+2), stal(0:ne+2))
  etal(:) = zero
  stal(:) = zero
  do i=1,ne
      read (unit,*) etal(i), stal(i)
  end do
  close (unit)
  !< to make sure that interpolation works (if TALYS spectrum short)
  de = etal(ne) - etal(ne-1)
  etal(ne+1) = etal(ne) + de !< next energy step spectrum = 0
  etal(ne+2) = max(etal(ne+1) + de, Enspec(Nnspec) + de ) !< also = 0 above maxE

  !< interpolate TALYS spectrum
  DCStal(:) = zero

  do i=1,Nnspec
      Ei = Enspec(i)
      DCStal(i) = interpolate_tab_func(etal, stal, Ei, 0)
  end do

  deallocate (etal, stal)

END SUBROUTINE retrieve_interpolate_talys_spectrum


!******************************************************************************
!> Interpolates differential neutron production function to alpha energy `Ea`
!! 
SUBROUTINE talys_an_differential_npf_for_alpha_E(Ea,nspec)
  USE all_data
  USE mhunt
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: Ea
  REAL(DP), DIMENSION(0:Nnspec), INTENT(OUT) :: nspec
  INTEGER :: jlo
  REAL(DP) :: w0, w1, E0, E1

  !< lin-lin interpolate tabulated DNPFc to alpha energy E0
  call hunt(Enpfc(:), Ea, jlo)
  jlo = lbound(Enpfc,1) - 1 + jlo
  E0 = Enpfc(jlo)
  E1 = Enpfc(jlo+1)
  w0 = (E1-Ea)/(E1-E0)  !< weight of y0
  w1 = (Ea-E0)/(E1-E0)  !< weight of y1
  nspec(:) = w0 * DNPFc(jlo,:) + w1 * DNPFc(jlo+1,:)

END SUBROUTINE talys_an_differential_npf_for_alpha_E


!******************************************************************************
!> Calculates neutron spectrum `nspec` from (a,n) with alpha of energy `E0`
!! uses simple reaction kinematics, ground state to ground state
!! 
SUBROUTINE local_an_differential_npf_for_alpha_E(E0,nspec)
  USE all_data
  USE mhunt
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: E0
  REAL(DP), DIMENSION(0:Nnspec), INTENT(OUT) :: nspec
  INTERFACE
      REAL(DP) FUNCTION interpolate_tab_func(xa,fa,xi,type)
        USE precision
        USE mhunt
        REAL(DP), DIMENSION(:), INTENT(IN) :: xa,fa
        REAL(DP), INTENT(IN) :: xi
        INTEGER, INTENT(IN) :: type
      END FUNCTION interpolate_tab_func
  END INTERFACE
  REAL(DP) :: npf0,dnpf,Ectr,E3min,E3max, dnspec(0:Nnspec)
  INTEGER :: jlo,jstp,jlomin,jhimax,jspan

  nspec(:) = zero
  call hunt(Enpfc(:),E0,jlo)
  jlo = lbound(Enpfc,1) - 1 + jlo
  if (jlo == lbound(Enpfc,1)-1 .or. jlo == ubound(Enpfc,1)) &
      STOP 'Out of bounds in hunt - an_neutron_spectrum_for_incoming_alpha!'

  !< step down energy discretization
  do jstp=jlo,0,-1
      if (jstp==jlo) then
          npf0 = interpolate_tab_func(Enpfc,NPFc,E0,0) !< lin-lin
          dnpf = npf0 - NPFc(jstp)
          Ectr = 0.5 * (E0 + Enpfc(jstp))
      else
          dnpf = NPFc(jstp+1) - NPFc(jstp)
          Ectr = 0.5 * (Enpfc(jstp+1) + Enpfc(jstp))
      end if

      !< calculate spectrum contribution at jstp/Ectr energy bin
      dnspec(:) = zero
      if (Ectr >= EthreshAN) then
          call nonrel_reaction_kinematics(m1an,m2an,m3an,m4an,Ectr,E3min,E3max)
          if (E3max > E3min) then
              call hunt(Enpfc(:),E3min,jlomin)
              jlomin = lbound(Enpfc,1) - 1 + jlomin
              call hunt(Enpfc(:),E3max,jhimax)
              jhimax = lbound(Enpfc,1) - 1 + jhimax
              jspan = jhimax - jlomin
              dnspec(jlomin:jhimax) = dnpf / (jspan*dEan)
          end if
      end if

      !< add contribution to spectrum
      nspec(:) = nspec(:) + dnspec(:)

  end do

END SUBROUTINE local_an_differential_npf_for_alpha_E


!******************************************************************************
!> Calculated nonrelativistic reaction kinematics according to 
!! Marmier & Sheldon, 1969, Appendix C.4, pp.612–614
!! 
SUBROUTINE nonrel_reaction_kinematics(m1,m2,m3,m4,E1,E3min,E3max)
  USE precision
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: m1,m2,m3,m4,E1
  REAL(DP), INTENT(OUT) :: E3min,E3max
  REAL(DP) :: Qval,E,Ethresh,m1234,fqe,A13,A14,A23,A24,th,thmax

  Qval = m1 + m2 - m3 - m4
  E = E1+Qval
  m1234 = (m1+m2)*(m3+m4)
  A13 = m1*m3/m1234 * E1/E
  A14 = m1*m4/m1234 * E1/E
  fqe = one + m1/m2*Qval/E
  A23 = m2*m3/m1234 * fqe
  A24 = m2*m4/m1234 * fqe

  E3min = zero
  E3max = zero
  if (Qval >= zero) then
      Ethresh = zero
  else
      Ethresh = -m1/m2*Qval !< kinematic threshold for E (not E1)!!!
  end if
  if (E < Ethresh) then  !< kinematic threshold not satisfied
      write (*,*) 'nonrel_reaction_kinematics: kinematic treshold not satisfied'
      STOP 'termination.'
  end if

  if (A13 > A24) then
      thmax = asin(sqrt(A24/A13))
      !< th = zero
      E3max = E * A13 * ( one - sqrt(A24/A13) )**2
      th = thmax
      E3min = E * A13 * ( cos(th) )**2
  else
      !< th = zero
      E3max = E * A13 * ( one + sqrt(A24/A13) )**2
      !< th = PI
      E3min = E * A13 * (-one + sqrt(A24/A13) )**2
  end if

END SUBROUTINE nonrel_reaction_kinematics


!******************************************************************************
!> Cleans up and checks neutron spectra :)
!!
SUBROUTINE cleanup_check_neutron_spectra
  USE all_data
  IMPLICIT NONE
  REAL(DP) :: tiny=1.e-8_dp, thresh, ints, rerr
  INTEGER :: ic

  !< cleanup
  do ic=1,ndecc
      thresh = maxval(dchain(ic)%nspecAN(:)) * tiny
      where (dchain(ic)%nspecAN(:) < thresh) &
          dchain(ic)%nspecAN(:) = zero
  end do

  !< check
  write (outp,'("Integrated ",a,"(a,n) energy spectra")') &
      trim(adjustl(target_nucl%nucid))
  do ic=1,ndecc
      ints = sum(dchain(ic)%nspecAN(:)) * dEan
      ints = ints - 0.5_dp * dEan * (dchain(ic)%nspecAN(0) + dchain(ic)%nspecAN(Nnspec))
      rerr = (ints / dchain(ic)%nyieldAN - one) * 100._dp
      write (outp,'(a,es12.3,"  (% diff w/ yield calc =",f8.3,")")') dchain(ic)%tpnucid,ints,rerr
      if (abs(rerr) > one) then
          write (*,'(6x,"WARNING: (a,n) spectrum diff wrt yield for ",a," = ",f6.1," %")') trim(adjustl(dchain(ic)%tpnucid)),rerr
          write (outp,'(50x,"WARNING: large difference!")') 
      end if
  end do
  write (outp,*)

END SUBROUTINE cleanup_check_neutron_spectra


!==============================================================================
!--- OUTPUT -------------------------------------------------------------------

!******************************************************************************
!> Outputs neutron production table :)
!!
SUBROUTINE output_neutron_prod_func_table
  USE all_data
  IMPLICIT NONE
  REAL(DP) :: Esave=0.1_dp
  INTEGER :: i,isave,unit
  CHARACTER :: filename*255

  filename = 'neutron_prod_func_'//trim(adjustl(target_nucl%nucid))//'.out'
  unit = 31
  call open_file(unit, trim(filename), 'unknown', 'write')

  isave = nint((Esave)/dEan)
  write (outp,'(/,a,"(a,n) neutron production function")') &
      trim(adjustl(target_nucl%nucid))
  write (outp,'(4x,"<<OUTFILE>> Also written to file: ",a)') trim(filename)
  write (outp,'(4x,a,4x,a)') 'E in MeV', 'n.p.f.'

  do i=lbound(Enpfc,1),ubound(Enpfc,1)
      write (unit,'(f8.4,es12.4)') Enpfc(i), NPFc(i)
      if (mod(i,isave)==0) &
          write (outp,'(f12.3,es12.3)') Enpfc(i), NPFc(i)
  end do
  close (unit)
  write (outp,*)
END SUBROUTINE output_neutron_prod_func_table


!******************************************************************************
!> Outputs neutron yield :)
!!
SUBROUTINE output_neutron_yield(it)
  USE all_data
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: it
  INTEGER :: ic,k
  REAL(DP) :: nyield, yd, ydw, ry, ryww, totrate
  write (outp,'(/,a,a)') trim(adjustl(target_nucl%nucid)),'(a,n) NEUTRON YIELD/PRODUCTION RATE'
  write (outp,111)
  totrate = zero
  do ic=1,ndecc
      write (outp,'(a," decay chain")') dchain(ic)%tpnucid
      do k=1,dchain(ic)%nAbr
          nyield = dchain(ic)%Abranch(k)%nyieldAN
          call calc_neutron_output(ic, nyield, yd, ydw, ry, ryww)
          write (outp,222) dchain(ic)%Abranch(k)%pnucid, &
              dchain(ic)%Abranch(k)%dnucid, yd, ydw, ry, ryww
      end do
      call calc_neutron_output(ic, dchain(ic)%nyieldAN, yd, ydw, ry, ryww)
      write (outp,333) yd, ydw, ry, ryww
      write (outp,*)
      totrate = totrate + ry
      !< save results
      results(it,ic,1:4) = (/ yd, ydw, ry, ryww /)
  end do
  write (outp,444) ndecc, totrate

111 format (3x,"parent -> daught.  per decay   per dec-wtf   per yr-kg   per yr-kg-wtf-wtf",/)
222 format (4x,a," -> ",a,x,2es12.3,2x,2es12.3)
333 format (4x,"TOTAL",10x,2es12.3,2x,2es12.3)
444 format (4x,"TOTAL (",i1," chains)",25x,es12.3,/)
END SUBROUTINE output_neutron_yield


!******************************************************************************
!> Calculates various outputs
!!
SUBROUTINE calc_neutron_output(ic, nyield, yd, ydw, ry, ryww)
  USE all_data
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ic
  REAL(DP), INTENT(IN) :: nyield
  REAL(DP), INTENT(OUT) :: yd, ydw, ry, ryww
  !< neutron yield per decay
  yd = nyield
  !< neutron yield per decay per wtf target element
  ydw = yd / target_elem_wtf
  !< neutron production rate per yr per kg of rock
  ry = yd * dchain(ic)%num2rate * dchain(ic)%tpwtf * secinyr
  !< neutron prod. rate per yr per wtf parent per wtf target elem per kg of rock
  ryww = ry / (dchain(ic)%tpwtf * target_elem_wtf)
END SUBROUTINE calc_neutron_output


!******************************************************************************
!> Outputs neutron spectra, user-friendly xy file + for MCNP6 input
!!
SUBROUTINE output_neutron_spectra
  USE all_data
  IMPLICIT NONE
  INTEGER :: unit, unitSI, unitSP, ic, j
  CHARACTER :: filename*255

  do ic=1,ndecc

      filename = 'neutron_spectrum_'//trim(adjustl(dchain(ic)%tpnucid))//'_'//trim(adjustl(target_nucl%nucid))//'.out'
      write (outp,'("<<OUTFILE>> Neutron (a,n) spectra written to file: ",a)') trim(filename)
      unit = 33
      call open_file(unit, trim(filename), 'unknown', 'write')

      do j=0,Nnspec
          write (unit,'(f8.4,es12.4)') Enspec(j), dchain(ic)%nspecAN(j)
      end do

      close (unit)

      !< input for MCNP6
      filename = 'sdef_si_'//trim(adjustl(dchain(ic)%tpnucid))//'_'//trim(adjustl(target_nucl%nucid))//'.mcnp6'
      write (outp,'("<<OUTFILE>> MCNP6 input written to file: ",a)') trim(filename)
      unitSI = 34
      call open_file(unitSI, trim(filename), 'unknown', 'write')

      filename = 'sdef_sp_'//trim(adjustl(dchain(ic)%tpnucid))//'_'//trim(adjustl(target_nucl%nucid))//'.mcnp6'
      write (outp,'("<<OUTFILE>> MCNP6 input written to file: ",a,/)') trim(filename)
      unitSP = 35
      call open_file(unitSP, trim(filename), 'unknown', 'write')

      write (unitSI,'(a)') 'SI1 a'
      write (unitSI,'(6x,7es10.3)') Enspec(0:Nnspec)
      write (unitSP,'(a)') 'SP1'
      write (unitSP,'(6x,7es10.3)') dchain(ic)%nspecAN(0:Nnspec)

      close (unitSI)
      close (unitSP)

  end do
  write (outp,*)

END SUBROUTINE output_neutron_spectra


