!> \file
!! Contains modules

MODULE precision
  INTEGER, PARAMETER :: SP = kind(1.e0) !< Single precision.
  INTEGER, PARAMETER :: DP = kind(1.d0) !< Double precision.
  REAL(DP), PARAMETER :: zero = 0.000000000000000000000000000000000000000_dp
  REAL(DP), PARAMETER :: one  = 1.000000000000000000000000000000000000000_dp
  REAL(DP), PARAMETER :: half = 0.500000000000000000000000000000000000000_dp
  REAL(DP), PARAMETER :: PI   = 3.141592653589793238462643383279502884197_dp
  REAL(DP), PARAMETER :: infty  = huge(one)
  !
  REAL(DP), PARAMETER :: secinday = 24._dp*3600._dp
  REAL(DP), PARAMETER :: secinyr = 365.25_dp*secinday
  REAL(DP), PARAMETER :: Avogadro = 6.02214129e23_dp
  REAL(DP), PARAMETER :: mass_neutron = 1.008664916_dp !< in u (=g/mol)
  REAL(DP), PARAMETER :: u_in_MeV = 931.494061_dp
END MODULE precision

!******************************************************************************
!> Module for input/output paramters
!!
MODULE inpout
  USE precision
  !< INPUT ********
  CHARACTER(127) :: model_name !< Model label to be used in filenamins
  CHARACTER(127) :: file_abundances !< name of file with rock composition
  !! Assumes the input file contains 2 columns:
  !!   1. Z of element
  !!   2. mass fraction of element
  !! Include major elements + K,Th,U
  CHARACTER(127) :: input_directory !< Directory where input files are located.
  INTEGER        :: which_an_calc=0 !< selector for (a,n) calculation options
  CHARACTER(127) :: local_directory !< Directory where local files are located.
  CHARACTER(127) :: talys_directory !< Directory where TALYS files are located.
  CHARACTER(127) :: file_atomic_mass !< Name of input file with atomic masses.
  !!
  !! Assumes the input file contains 6 columns:
  !!   1. Atomic number
  !!   2. Atomic symbol
  !!   3. Mass number
  !!   4. Relative nuclide mass in u (=g/mol)
  !!   5. Natural isotopic composition (isotopic fraction)
  !!   6. Standard atomic weight in u (=g/mol)
  !!
  !! Nuclide mass, atomic weight are in unified atomic mass units (u or Da),
  !! i.e, normalized to carbon-12 and equivalent to g/mol. If isotopic
  !! composition is unknown, put zero in the input file. If standard weight 
  !! is unknown, put mass number.
  CHARACTER(127) :: file_decay_chains !< Name of input file with NUCIDs of
  !! top parent in each decay chain to be included (e.g., 232TH). One NUCID
  !! (chain) per line in the input file.
  CHARACTER(127) :: file_stopping_power !< Name of input file with alpha 
  !! stopping power
  !! Assumes the input file contains 2 columns:
  !!   1. Energy in MeV
  !!   2. Mass stopping power in Mev cm2/g
  CHARACTER(127) :: file_an_targets !< Name of input file with NUCIDs of
  !! (a,n) target nuclides to be included (e.g., "18O") and "postfixes" of the
  !! cross section file name (e.g., "_tendl-2012.xy"). The cross section 
  !! filename will be NUCID//postfix (e.g., "18O_tendl-2012.xy").
  !! Assumes the input files contains 2 columns:
  !!   1. NUCID
  !!   2. postfix
  !! Assumes the cross section library file contains 2 columns:
  !!   1. Energy in MeV
  !!   2. Cross section in mbarn
  REAL(DP) :: rock_density !< rock density in g/cm3
  REAL(DP) :: dEtalys      !< energy stepping for TALYS cross sections in MeV
  REAL(DP) :: dEan         !< energy stepping for neutron production in MeV
  REAL(DP) :: mcnp6_accu   !< desired accuracy of MCNP6 tally
  REAL(DP) :: mcnp6_radcm  !< radius of spherical volume for MCNP6 calculation
  !
  !< NOT IN NAMELIST
  CHARACTER(127) :: output_file            !< Name of main output file.
  INTEGER, PARAMETER :: outp=1             !< output file unit number
  !
  NAMELIST /inputs/          &
      model_name,            &
      input_directory,       &
      which_an_calc,         &
      local_directory,       &
      talys_directory,       &
      file_atomic_mass,      &
      file_abundances,       &
      file_decay_chains,     &
      file_stopping_power,   &
      file_an_targets,       &
      rock_density,          &
      dEtalys,               &
      dEan,                  &
      mcnp6_accu,            &
      mcnp6_radcm
END MODULE inpout

!******************************************************************************
!> Derived type definition of type for chemical element.
!!
MODULE elem_type
  USE precision
  TYPE elem !< type for chemical element
      INTEGER       :: Z    !< atomic number
      CHARACTER(2)  :: symb !< symbol for element
      CHARACTER(12) :: name !< name of element
      REAL(DP)      :: wt   !< standard atomic weight in u (=g/mol)
  END type elem
END MODULE elem_type

!******************************************************************************
!> Derived type definition of type for nuclide.
!!
MODULE nucl_type
  USE precision
  TYPE nucl !< type for nuclide
      INTEGER      :: Z     !< atomic number
      INTEGER      :: A     !< mass number
      CHARACTER(5) :: nucid !< NUCID (e.g., "232TH", " 14C ", ...)
      REAL(DP)     :: wt    !< nuclide mass in u (=g/mol)
      REAL(DP)     :: nab   !< natural abundance (isotopic fraction)
  END type nucl
END MODULE nucl_type

!******************************************************************************
!> Derived type definition of type for alpha level.
!!
MODULE alev_type
  USE precision
  TYPE alev !< type for alpha level
      REAL(DP) :: E  !< level energy in MeV
      REAL(DP) :: I  !< intensity (as a fraction, not percent)
  END type alev
END MODULE alev_type

!******************************************************************************
!> Derived type definition of type for alpha decay branch.
!!
MODULE abra_type
  USE alev_type
  TYPE abra !< type for alpha decay branch
      CHARACTER(5)            :: pnucid    !< parent nuclide NUCID
      CHARACTER(5)            :: dnucid    !< daughter nuclide NUCID
      CHARACTER(30)           :: dsid      !< DSID of ENSDF
      REAL(DP)                :: t12       !< half-life in seconds
      REAL(DP)                :: Qval      !< ground state Q value in MeV
      REAL(DP)                :: br        !< branching ratio
      REAL(DP)                :: tbr       !< total branching ratio
      INTEGER                 :: nlev=0    !< number of levels
      TYPE(alev), ALLOCATABLE :: alevel(:) !< levels
      REAL(DP)                :: nyieldAN  !< (a,n) neutron yield
  END type abra
END MODULE abra_type

!******************************************************************************
!> Derived type definition of type for beta- decay branch.
!!
MODULE bbra_type
  USE precision
  TYPE bbra !< type for beta- decay branch
      CHARACTER(5)            :: pnucid    !< parent nuclide NUCID
      CHARACTER(5)            :: dnucid    !< daughter nuclide NUCID
      CHARACTER(30)           :: dsid      !< DSID of ENSDF
      REAL(DP)                :: t12       !< half-life in seconds
      REAL(DP)                :: Qval      !< ground state Q value in MeV
      REAL(DP)                :: br        !< branching ratio
      REAL(DP)                :: tbr       !< total branching ratio
  END type bbra
END MODULE bbra_type

!******************************************************************************
!> Derived type definition of type for decay chain.
!!
MODULE decch_type
  USE abra_type
  USE bbra_type
  TYPE decch !< type for decay chain
      CHARACTER(5)            :: tpnucid    !< top parent NUCID
      REAL(DP)                :: tpwtf      !< top parent wt.frac.
      INTEGER                 :: nAbr=0     !< number of alpha branches
      TYPE(abra), ALLOCATABLE :: Abranch(:) !< alpha branches
      INTEGER                 :: nBbr=0     !< number of beta branches
      TYPE(bbra), ALLOCATABLE :: Bbranch(:) !< beta branches
      INTEGER                 :: nalpha     !< number of alpha's emitted
      REAL(DP)                :: maxEa=zero !< max alpha energy
      REAL(DP)                :: minEa=infty!< min alpha energy
      REAL(DP)                :: nyieldAN   !< (a,n) neutron yield
      REAL(DP)                :: nyieldSF   !< spontaneous fission neutron yield
      REAL(DP)                :: num2rate   !< dec-to-rate conversion factor "C"
      REAL(DP), ALLOCATABLE   :: nspecAN(:) !< (a,n) neutron prod. energy spec.
      REAL(DP), ALLOCATABLE   :: nspecSF(:) !< spont.fiss. neutron energy spec.
      REAL(DP)                :: SFbran=zero!< branching ratio for spont.fission
      REAL(DP)                :: SFneut=zero!< neutrons per spont.fission
      REAL(DP)                :: WattA=zero !< Watt spectrum parameter "a"
      REAL(DP)                :: WattB=zero !< Watt spectrum parameter "b"
  END type decch
END MODULE decch_type

!******************************************************************************
!> Derived type definition of type for chemical abundance.
!!
MODULE abun_type
  USE precision
  TYPE abun !< type for abundance
      INTEGER  :: Z      !< atomic number
      REAL(DP) :: wtf    !< weight fraction
      REAL(DP) :: atf    !< atomic fraction
      REAL(DP) :: atgram !< atoms per gram of rock
      REAL(DP) :: atcm3  !< atoms per cm3 of rock
  END type abun
END MODULE abun_type

!******************************************************************************
!> Module for all data
!!
MODULE all_data
  USE inpout
  USE elem_type
  USE nucl_type
  USE abun_type
  USE decch_type
  INTEGER, PARAMETER :: nelem=94        !< number of chemical elements
  TYPE(elem) :: element(nelem)          !< chemical element array
  INTEGER :: nnucl=0                    !< number of nuclides
  TYPE(nucl), TARGET, ALLOCATABLE :: nuclide(:) !< nuclide array
  INTEGER :: nabund=0                   !< number of elem. components
  TYPE(abun), ALLOCATABLE :: abund(:)   !< chem. abundance array
  INTEGER :: ndecc=0                    !< number of decay chains
  TYPE(decch), ALLOCATABLE :: dchain(:) !< decay chains
  REAL(DP) :: maxEalpha=zero            !< maximum alpha particle energy
  REAL(DP) :: minEalpha=infty           !< maximum alpha particle energy
  REAL(DP) :: En_cutoff                 !< cuttoff energy of neutron spectra
  !
  INTEGER :: nsp=0                      !< number of stopping power data points
  REAL(DP), ALLOCATABLE :: stoppow(:,:) !< stopping power array
  !
  INTEGER :: ntgts                      !< number of (a,n) targets
  TYPE(nucl), ALLOCATABLE :: targets(:) !< (a,n) targets
  CHARACTER(255), TARGET, ALLOCATABLE :: target_fnames(:) !< xsec filenames
  REAL(DP), ALLOCATABLE :: results(:,:,:) !< summary results array
  REAL(DP) :: m1an,m2an,m3an,m4an       !< masses of particles in (a,n) in MeV
  REAL(DP) :: QvalAN                    !< Q value of (a,n) in MeV
  REAL(DP) :: EthreshAN                 !< kinematic energy threshold
  INTEGER :: nEaCS=0                    !< number of cross section data points in Ealpha
  INTEGER :: nEnDCS=0                   !< number of diff. cross section data points in Eneutron
  REAL(DP), ALLOCATABLE ::   CS(:)      !< (a,n) cross section
  REAL(DP), ALLOCATABLE ::   EaCS(:)    !< alpha energy grid for (a,n) CS
!!$  REAL(DP), ALLOCATABLE ::   EnDCS(:)   !< neutron energy grid for (a,n) diff.CS
  REAL(DP), ALLOCATABLE ::   DCS(:,:)   !< (a,n) differential cross section
  REAL(DP) :: minE_CS                   !< lower E bound of (a,n) xsec
  TYPE(nucl), POINTER :: target_nucl    !< (a,n) target nuclide
  REAL(DP) :: target_nucl_atcm3         !< atomic density of target nuclide
  REAL(DP) :: target_elem_wtf           !< wt.frac. of target element
  INTEGER :: Nnpfc=0                    !< number of points in npftab
  INTEGER :: Nnspec=0                   !< number of points in n-spectrum
  REAL(DP), ALLOCATABLE :: Enpfc(:)     !< neutron production energy array
  REAL(DP), ALLOCATABLE :: NPFc(:)      !< neutron production P_i array
  REAL(DP), ALLOCATABLE :: Enspec(:)    !< neutron spectra energy array
  REAL(DP), ALLOCATABLE :: DNPFc(:,:)   !< diff. neutron production DP_i array
  INTEGER :: indx_dnpf_nj               !< neutron energy index for func_dnpf_nj
END MODULE all_data
