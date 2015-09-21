!> \file
!! Neutron production by spontaneous fission


!******************************************************************************
!> 
!!
SUBROUTINE calculate_SF_neutron_production
  USE all_data
  IMPLICIT NONE
  INTEGER :: ic, unit, i, j, unitSI, unitSP
  CHARACTER :: tpnucid*5, filename*128
  REAL(DP) :: Ej, Sj, a, b, minE, maxE, ints

  write (*,*) '  Calculating spontaneous fission neutron yields and spectra'
  do ic=1,ndecc
      dchain(ic)%nyieldSF = zero
      tpnucid = dchain(ic)%tpnucid

      !< READ SF INPUTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      filename = trim(input_directory)//'/'//trim(adjustl(tpnucid))//'.sf'
      unit = 38
      write (outp,'("Reading ",a," spont. fission parameters from file: ",a)') &
          trim(adjustl(tpnucid)), trim(filename)
      call open_file(unit, trim(filename), 'old', 'read')
      read (unit,*) dchain(ic)%SFbran
      read (unit,*) dchain(ic)%SFneut
      read (unit,*) dchain(ic)%WattA
      read (unit,*) dchain(ic)%WattB
      close (unit)
      write (outp,'("SF branching              = ",es8.2)') dchain(ic)%SFbran
      write (outp,'("SF neutrons per fission   = ",f4.2)') dchain(ic)%SFneut
      write (outp,'("Watt spectrum parameter a = ",f6.4," MeV")') dchain(ic)%WattA
      write (outp,'("Watt spectrum parameter b = ",f6.4," MeV")') dchain(ic)%WattB

      !< CALCULATE SF NEUTRON YIELD AND PRODUCTION RATE >>>>>>>>>>>>>>>>>>>>>>>
      dchain(ic)%nyieldSF = dchain(ic)%SFbran * dchain(ic)%SFneut

      !< SF neutron yield per decay
      results(ntgts+1,ic,1) = dchain(ic)%nyieldSF
      !< SF neutron production rate per yr per kg of rock
      results(ntgts+1,ic,3) = dchain(ic)%nyieldSF * dchain(ic)%num2rate &
          * dchain(ic)%tpwtf * secinyr
      !< SF neutron prod. rate per yr per wtf parent per per kg of rock
      results(ntgts+1,ic,4) = dchain(ic)%nyieldSF * dchain(ic)%num2rate &
          * secinyr

      !< CALCULATE SF NEUTRON SPECTRA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      minE = zero
      maxE = 10._dp !< MeV
      Nnspec = nint((maxE-minE)/dEan)
      allocate (Enspec(0:Nnspec), dchain(ic)%nspecSF(0:Nnspec))
      Enspec(0:Nnspec) = (/ ( minE+i*dEan,i=0,Nnspec) /)

      !< Watt fission spectrum
      a = dchain(ic)%WattA
      b = dchain(ic)%WattB
      do j=0,Nnspec
          Ej = Enspec(j)
          Sj = exp(-Ej/a) * sinh(sqrt(b*Ej))
          dchain(ic)%nspecSF(j) = Sj
      end do
      !< normalize spectrum
      ints = sum(dchain(ic)%nspecSF(:)) * dEan &
          - 0.5_dp * dEan * (dchain(ic)%nspecSF(1) + dchain(ic)%nspecSF(Nnspec))
      dchain(ic)%nspecSF(:) = dchain(ic)%nspecSF(:) / ints * dchain(ic)%nyieldSF

      !< output spectrum
      filename = 'neutron_spectrum_'//trim(adjustl(dchain(ic)%tpnucid))//'_SF.out'
      write (outp,'(/,"<<OUTFILE>> Neutron SF spectra written to file: ",a)') trim(filename)
      unit = 39
      call open_file(unit, trim(filename), 'unknown', 'write')
      do j=0,Nnspec
          write (unit,'(f8.4,es12.4)') Enspec(j), dchain(ic)%nspecSF(j)
      end do
      close (unit)

      !< input for MCNP6
      filename = 'sdef_si_'//trim(adjustl(dchain(ic)%tpnucid))//'_SF.mcnp6'
      write (outp,'("<<OUTFILE>> MCNP6 input written to file: ",a)') trim(filename)
      unitSI = 34
      call open_file(unitSI, trim(filename), 'unknown', 'write')

      filename = 'sdef_sp_'//trim(adjustl(dchain(ic)%tpnucid))//'_SF.mcnp6'
      write (outp,'("<<OUTFILE>> MCNP6 input written to file: ",a,//)') trim(filename)
      unitSP = 35
      call open_file(unitSP, trim(filename), 'unknown', 'write')

      write (unitSI,'(a)') 'SI1 a'
      write (unitSI,'(6x,7es10.3)') Enspec(0:Nnspec)
      write (unitSP,'(a)') 'SP1'
      write (unitSP,'(6x,7es10.3)') dchain(ic)%nspecSF(0:Nnspec)

      close (unitSI)
      close (unitSP)

      deallocate(Enspec)
  end do

END SUBROUTINE calculate_SF_neutron_production
