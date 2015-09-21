!> \file
!! Elements & nuclides setup


!******************************************************************************
!> Sets up element and nuclide data structures from input files
!!
SUBROUTINE read_elements_nuclides
  USE all_data
  IMPLICIT NONE
  CHARACTER :: filename*255, symb*2
  INTEGER :: unit,Z,A
  REAL(DP) :: nwt,nab,awt
  TYPE(elem) :: el
  TYPE(nucl) :: nu

  nnucl = 0
  allocate (nuclide(nnucl))
  filename = trim(input_directory)//'/'//trim(file_atomic_mass)
  write (outp,'("Reading elements/nuclides from file: ",a)') trim(filename)
  unit = 21
  call open_file(unit, trim(filename), 'old', 'read')
  do
      read (unit, *, end=99) Z,symb,A,nwt,nab,awt
      if (Z <= nelem) then
          element(Z)%Z = Z
          call upcase(symb) !< convert symbol to all upper case
          if (Z==1) symb = 'H ' !< special fix for hydrogen to prevent D or T
          element(Z)%symb = symb
          element(Z)%wt = awt
          call append_nuclide(Z,A,nwt,nab)
      end if
  end do
99 close (unit)

  call fill_element_names

  !< some outputs
  if (all(element(:)%Z > 0)) then
      write (*,*) 'Element inputs read =',nelem
      write (outp,'(a,i5)') 'Element inputs read =',nelem
      do Z=1,nelem
          el = element(Z)
          write (outp,'(4x,i3,a4,2x,a,f12.6)') el%Z, el%symb, el%name, el%wt
      end do
      write (outp,*)
  else
      write (outp,*) 'ERROR: Problem reading elements...'
      write (*   ,*) 'ERROR: Problem reading elements...'
  end if
  write (*,*) 'Number of nuclide inputs read =',nnucl
  write (outp,'(a,i5)') 'Number of nuclide inputs read =',nnucl
  write (outp,'(a)') 'Writing only naturally occuring nuclides'
  do Z=1,nnucl
      nu = nuclide(Z)
      if (nu%nab > zero) &
          write (outp,'(4x,a5,i4,i5,f12.6,f13.9)') &
          nu%nucid, nu%Z, nu%A, nu%wt, nu%nab
  enddo
  write (outp,*)
END SUBROUTINE read_elements_nuclides


!******************************************************************************
!> Adds nuclide to the nuclide array by expanding its size
!!
SUBROUTINE append_nuclide(Z,A,wt,nab)
  USE all_data
  USE ensdf
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: Z,A
  REAL(DP), INTENT(IN) :: wt,nab
  TYPE(nucl), ALLOCATABLE :: tmp(:)
  allocate (tmp(nnucl))
  tmp(1:nnucl) = nuclide(1:nnucl) !< copy nuclide array to tmp
  deallocate (nuclide)
  allocate (nuclide(nnucl+1))     !< expand nuclide array by one element
  nuclide(1:nnucl) = tmp(1:nnucl) !< copy back from tmp
  deallocate (tmp)
  nnucl = nnucl + 1               !< add new nuclide data
  nuclide(nnucl)%Z = Z
  nuclide(nnucl)%A = A
  nuclide(nnucl)%nucid = AZ_to_NUCID(A,Z)
  nuclide(nnucl)%wt = wt
  nuclide(nnucl)%nab = nab
END SUBROUTINE append_nuclide


!******************************************************************************
!> Simply fills names of elements in element array
!!
SUBROUTINE fill_element_names
  USE all_data
  IMPLICIT NONE
  element(1:94)%name = (/ &
      'Hydrogen    ', 'Helium      ', 'Lithium     ', 'Beryllium   ', &
      'Boron       ', 'Carbon      ', 'Nitrogen    ', 'Oxygen      ', &
      'Fluorine    ', 'Neon        ', 'Sodium      ', 'Magnesium   ', &
      'Aluminium   ', 'Silicon     ', 'Phosphorus  ', 'Sulfur      ', &
      'Chlorine    ', 'Argon       ', 'Potassium   ', 'Calcium     ', &
      'Scandium    ', 'Titanium    ', 'Vanadium    ', 'Chromium    ', &
      'Manganese   ', 'Iron        ', 'Cobalt      ', 'Nickel      ', &
      'Copper      ', 'Zinc        ', 'Gallium     ', 'Germanium   ', &
      'Arsenic     ', 'Selenium    ', 'Bromine     ', 'Krypton     ', &
      'Rubidium    ', 'Strontium   ', 'Yttrium     ', 'Zirconium   ', &
      'Niobium     ', 'Molybdenum  ', 'Technetium  ', 'Ruthenium   ', &
      'Rhodium     ', 'Palladium   ', 'Silver      ', 'Cadmium     ', &
      'Indium      ', 'Tin         ', 'Antimony    ', 'Tellurium   ', &
      'Iodine      ', 'Xenon       ', 'Cesium      ', 'Barium      ', &
      'Lanthanum   ', 'Cerium      ', 'Praseodymium', 'Neodymium   ', &
      'Promethium  ', 'Samarium    ', 'Europium    ', 'Gadolinium  ', &
      'Terbium     ', 'Dysprosium  ', 'Holmium     ', 'Erbium      ', &
      'Thulium     ', 'Ytterbium   ', 'Lutetium    ', 'Hafnium     ', &
      'Tantalum    ', 'Tungsten    ', 'Rhenium     ', 'Osmium      ', &
      'Iridium     ', 'Platinum    ', 'Gold        ', 'Mercury     ', &
      'Thallium    ', 'Lead        ', 'Bismuth     ', 'Polonium    ', &
      'Astatine    ', 'Radon       ', 'Francium    ', 'Radium      ', &
      'Actinium    ', 'Thorium     ', 'Protactinium', 'Uranium     ', &
      'Neptunium   ', 'Plutonium   ' /)
END SUBROUTINE fill_element_names


!******************************************************************************
!> Reads rock composition from input file.
!! \note
!! Modification to exclude 18O output into ZAID-atf file because of 
!! incomplete data  libraries on MCNP6
SUBROUTINE read_abundances
  USE all_data
  IMPLICIT NONE
  CHARACTER :: filename*255
  INTEGER :: unit, i, j, Z, ZAID, ZAID18O=8018, ZAID12C=6012, ZAID13C=6013
  REAL(DP) :: wtfrac, sumwt, summol, sum3, sum4, atf1, atf18O, atf12C, atf13C

  filename = trim(input_directory)//'/'//trim(file_abundances)
  !filename = trim(file_abundances)
  unit = 21
  call open_file(unit, trim(filename), 'old', 'read')
  write (*,*) 'Reading element abundances'
  write (outp,'("Reading element abundances from file: ",a)') trim(filename)

  !< get number of entries
  nabund = 0
  do 
      read (unit, *, end=99)
      nabund = nabund + 1
  end do
99 rewind (unit)

  !< read composition
  sumwt = zero
  summol = zero
  allocate (abund(nabund))
  do i=1,nabund
      read (unit, *) Z, wtfrac
      abund(i)%Z = Z
      abund(i)%wtf = wtfrac
      abund(i)%atf = wtfrac/element(Z)%wt
      abund(i)%atgram = wtfrac/element(Z)%wt*Avogadro
      abund(i)%atcm3 = abund(i)%atgram*rock_density
      sumwt = sumwt + wtfrac
      summol = summol + wtfrac/element(Z)%wt
  end do
  close (unit)
  abund(:)%atf = abund(:)%atf/summol
  summol = sum(abund(:)%atf)
  sum3 = sum(abund(:)%atgram)
  sum4 = sum(abund(:)%atcm3)

  !< output elemental composition
  write (outp,'(26x,a)') 'wt.frac     at.frac     atoms/g     atoms/cm3'
  do i=1,nabund
      write (outp,'(4x,a,i4,2x,4es12.2)') element(abund(i)%Z)%name, abund(i)
  end do
  write (outp,'(4x,"TOTAL",13x,2f12.8,2es12.2,/)') sumwt, summol, sum3, sum4

  write (outp,'("Rock density =",f6.3," g/cm3",/)') rock_density
  
  !< output isotopic composition for MCNP6 input ("ZAID atom.frac." pairs) 
  filename = 'zaid_atomfrac_'//trim(adjustl(model_name))//'.mcnp6'
  write (outp,'("<<OUTFILE>> ZAID-atomfrac composition written to file: ",a)') trim(filename)
  unit = 19
  call open_file(unit, trim(filename), 'unknown', 'write')
  write (unit,'(a)') 'M10'
  do i=1,nabund
      Z = abund(i)%Z
      do j=1,nnucl
          if (nuclide(j)%Z == Z) then
              ZAID = Z*1000 + nuclide(j)%A
              atf1 = abund(i)%atf * nuclide(j)%nab
              if (atf1 > zero) then
                  if (ZAID == ZAID18O) then !< MODIF FOR MCNP6 (no data for 18O)
                      atf18O = atf1
                  else if (ZAID == ZAID12C) then !< MODIF FOR MCNP6 (no data)
                      atf12C = atf1
                  else if (ZAID == ZAID13C) then !< MODIF FOR MCNP6 (no data)
                      atf13C = atf1
                  else
                      write (unit,'(4x,i6,2x,es10.4)') ZAID, atf1
                  end if
              end if
          end if
      end do
  end do
  write (unit,'("c ",4x,i6,2x,es10.4)') ZAID12C, atf12C
  write (unit,'("c ",4x,i6,2x,es10.4)') ZAID13C, atf13C
  write (unit,'("c ",4x,i6,2x,es10.4)') ZAID18O, atf18O
  
  close (unit)

END SUBROUTINE read_abundances


!******************************************************************************
!> Returns index `inuc` of nuclide `Z,A`
!!
SUBROUTINE find_nuclide_index(A,Z,inuc)
  USE all_data
  USE ensdf
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: A    !< mass number
  INTEGER, INTENT(IN)  :: Z    !< atomic number
  INTEGER, INTENT(OUT) :: inuc !< nuclide index
  INTEGER :: i
  CHARACTER(5) :: nucid
  do i=1,nnucl
      if (nuclide(i)%Z /= Z .or. nuclide(i)%A /= A) cycle
      inuc = i
      RETURN
  end do
  nucid = AZ_to_NUCID(A,Z)
  write (*,*) nucid, " - could not find nuclide in find_nuclide_index :("
  STOP 'ERROR termination.'
END SUBROUTINE find_nuclide_index
