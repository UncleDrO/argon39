!> \file
!! Main program file

!******************************************************************************
!> Main program
!!
PROGRAM ARGON39
  IMPLICIT NONE
  INTEGER :: ntargets,it

  write (*,*) '************************************'
  write (*,*) 'Starting program ARGON39.'

  call get_inputs_start_output

  call read_elements_nuclides

  call read_abundances

  call setup_decay_chains

  call read_stopping_power

  call calculate_alpha_range


  call output_possible_targets

  call get_requested_an_targets(ntargets)

  do it=1,ntargets

      call setup_target_nuclide(it)

      call read_an_cross_section(it)

      call calculate_an_neutron_yield(it)

      call calculate_an_neutron_spectrum

  end do

  call calculate_SF_neutron_production


  call finalize_output

END PROGRAM ARGON39
