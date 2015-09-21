Detailed Description    {#detail}
====================

# Publication #    {#publication}

O. Šrámek, L. Stevens, W. F. McDonough, S. Mukhopadhyay, and J. Peterson: “Subterranean production of neutrons, <sup>39</sup>Ar and <sup>21</sup>Ne: Rates and uncertainties.” In preparation.

# Overview #    {#overview}

Software `argon39` handles the calculation of nucleogenic production rates of neutrons, <sup>21</sup>Ne and <sup>39</sup>Ar in underground geological settings. Prepares input files for calculations of cross sections and spectra using [TALYS](http://www.talys.eu) code. Prepares input files for neutron history simulations using [MCNP6](http://mcnp.lanl.gov) code and reads results from its output files.

* Alpha particles are produced in decay chains of long-lived natural radionuclides (mainly <sup>232</sup>Th and <sup>238</sup>U).
*  Neutrons are produced by (&alpha;,n) reactions on light nuclides and by spontaneous fission of <sup>238</sup>U.
* <sup>21</sup>Ne is produced by (&alpha;,n) on <sup>18</sup>O.
* <sup>39</sup>Ar is produced by (n,p) reaction on <sup>39</sup>K.

# Four programs #    {#fourprog}

Contains 4 separate programs; the corresponding executables are 

* `argon39.x`
* `talysprep.x`
* `mcnp6extract.x`
* `mcnp6prep.x`

## Basic calculation sequence ##    {#sequence}

__0a.__ (optional) run `talysprep.x` to prepare input files for TALYS.

__0b.__ (optional) Run TALYS to calculate (&alpha;,n) cross sections and spectra.

__1.__ Run `argon39.x` to calculate &alpha; production rates and spectra, &alpha; range, neutron production functions, neutron yields and production rates and spectra

__2.__ Run `mcnp6prep.x` to prepare input files for MCNP6.

__3.__ Run MCNP6 simulations.

__4.__ Run `mcnp6extract.x` to extract results from MCNP6 output files.

As TALYS-calculated (&alpha;,n) cross sections and spectra are included in the distribution, steps 0a and 0b are optional.


## argon39.x ##   {#argon39}

### Description ###    {#desc1}

`argon39.x` calculates natural &alpha; particle emission and neutron production. Uses (&alpha;,n) cross sections and spectra calculated by [TALYS](http://www.talys.eu). List of the main subroutines with brief description follows:

* `get_inputs_start_output`
* `read_elements_nuclides` ... Reads in atomic & nuclide mass data and nuclide natural abundances (data from [NIST](http://www.nist.gov/pml/data/comp.cfm))
* `read_abundances` ... Reads in rock composition (see publication for references)
* `setup_decay_chains` ... Sets up requested decay chains by reading in parameters from [ENSDF](http://www.nndc.bnl.gov) decay datasets
* `read_stopping_power` ... Reads in stopping power calculated externally (e.g., by [SRIM](http://srim.org))
* `calculate_alpha_range`
* `output_possible_targets`
* `get_requested_an_targets`
* `setup_target_nuclide`
* `read_an_cross_section` ... Reads in (&alpha;,n) cross section, calculated externally (e.g., by [TALYS](http://www.talys.eu)), for requested target nuclides
* `calculate_an_neutron_yield` ... Calculates (&alpha;,n) neutron yield 
* `calculate_an_neutron_spectrum` ... Calculates neutron energy spectra
* `calculate_SF_neutron_production`
* `finalize_output`

### Execution ###    {#exec1}

In subdirectory `example/`, run

    ../argon39.x argon39_example.inp

### Input file ###    {#inputfile}

The name of the input file (e.g., `argon39_example.inp`) is the one required command line parameter. The input file contains the following inputs:

    model_name = "example",

Main output file `argon39_MODELNAME.out` is created upon execution.

    input_directory = "../inputdata",

Directory where input data are located. All following input files are assumed to reside in this directory.

    file_abundances = "elem_wtfrac_RG03UC+C",

Name of 2-column input file with elemental rock composition:
 1. Atomic number (Z)
 2. Weight fraction of element


    file_atomic_mass = "atom_wt_nist_full.94",

Name of a 6-column input file with atomic masses:
 1. Atomic number (Z)
 2. Atomic symbol
 3. Mass number (A)
 4. Relative nuclide mass in u (= g/mol)
 5. Natural isotopic composition (isotopic fraction)
 6. Standard atomic weight in u (= g/mol)

If isotopic composition is unknown, put zero in the input file. If standard weight is unknown, put mass number.

    file_decay_chains = "decay.chains",

Name of input file with NUCIDs of top parent in each decay chain to be included (e.g., 232TH), one NUCID per line.

    file_stopping_power = "srim_tot_RG03UC_14ThU",

Name of a 2-column input file with mass stopping power of &alpha; particle in the rock:
 1. Energy of &alpha; particle in MeV
 2. Mass stopping power in MeV.cm2/g


    file_an_targets = "targets_tendl-2012.an",

Name of input file with NUCIDs of (&alpha;,n) target nuclides (e.g., 27Al), one per line.

    which_an_calc = 2,    !! 1=local(kinematic), 2=talys
    local_directory = "../localxs",
    talys_directory = "../talysxs",

Leave `which_an_calc = 2`. `talys_directory` is the name of the directory which contains the structure of subdirectories with (&alpha;,n) cross sections and spectra calculated by [TALYS](http://www.talys.eu) for the desired target nuclides. These cross sections/spectra are included in this distribution, but may also be calculated anew (see @ref talysprep).

    rock_density = 2.70,    !CRUST2.0 average density in layers 3-5

Material density of the rock in g/cm<sup>3</sup>.

    dEan = 0.01,

Desired energy stepping (in MeV) of neutron production calculations.

    dEtalys = 0.1,

Energy stepping (in MeV) of TALYS cross sections, relevant for @ref talysprep. 

    mcnp6_accu = 0.005,
    mcnp6_radcm = 300.,

Accuracy criterion for MCNP6 simulations and the radius (in cm) of spherical domain of the MCNP6 simulations. Relevant for @ref mcnp6prep and @ref mcnp6extract.

### Output ###    {#outp1}

* `argon39_MODELNAME.out` ... main output file, description of what has been done.
* `alpha_energy_spectrum_PPP.out` ... spectrum of &alpha; particles from PPP (e.g., `232TH`) decay chains
* `alpha_energy_table_PPP.txt` ... table of &alpha; particle emission from PPP decay chains
* `alpha_range.out` ... range of &alpha; particles
* `neutron_prod_func_TTT.out` ... thick target neutron production function for target nuclide TTT (e.g.. `27AL`)
* `neutron_spectrum_232TH_13C.out` ... (&alpha;,n) neutron spectrum with &alpha;s from PPP chain and TTT target
* `RESULTS_nrate_yr-kg-wtf-wtf.out` ... table of (&alpha;,n) and spontaneous fission neutron production rate per year per kg of rock per wt.frac. of PPP element (e.g., thorium) per wt.frac. of TTT element (e.g., aluminum)
* `RESULTS_nrate_yr-kg.out` ... table of (&alpha;,n) and spontaneous fission neutron production rate per year per kg of rock
* `RESULTS_nyield_per_decay.out` ... table of (&alpha;,n) and spontaneous fission neutron yield
* `sdef_si_PPP_TTT.mcnp6` and `sdef_sp_PPP_TTT.mcnp6` ... neutron spectrum input files for MCNP6 simulation
* `zaid_atomfrac_example.mcnp6` ... input file for MCNP6 simulation, nuclide atomic fractions in the rock


## talysprep.x ##    {#talysprep}

### Description ###    {#desc2}

`talysprep.x` can be used to prepare the calculation of (&alpha;,n) cross sections and spectra by [TALYS](http://www.talys.eu), if desired and if TALYS is installed on your machine.

Uses same @ref inputfile as @ref argon39.

### Execution ###    {#exec2}

In subdirectory `example/`, run

    ../talysprep.x argon39_example.inp

### Output ###    {#outp2}

* `TALYSprep.out` ... main output file with basic information.
* `TALYSbatch_example.sh` ... bash script to run TALYS calculations.

## mcnp6prep.x ##    {#mcnp6prep}

### Description ###    {#desc3}

`mcnp6prep.x` prepares input files for MCNP6 simulations.

Uses same @ref inputfile as @ref argon39.

### Execution ###    {#exec3}

In subdirectory `example/`, run

    ../mcnp6prep.x argon39_example.inp

The program interactively asks for some values. You can just enter 0s the first time. Then run MCNP6 in the "initiate" mode (`mcnp6 i name=inputfile`), then inspect the output file (input file name with appended `o`) to calculate input quantities asked for during the second run of `mcnp6prep.x`. _[Of course, this could be automated too...]_

### Output ###    {#outp3}

* `MCNP6batch_example.sh` ... ... bash script to run TALYS calculations.
* `mcnp6_MODELNAME_PPP_TTT` ... series of input files for MCNP6, where `PPP` is NUCID of decay chain parent and `TTT` is NUCID of (&alpha;,n) target, e.g., `mcnp6_example_232TH_27AL`.

## mcnp6extract.x ##    {#mcnp6extract}

### Description ###    {#desc4}

Extracts results from MCNP6 output files, after MCNP simulations have been run.

Uses same @ref inputfile as @ref argon39.

### Execution ###    {#exec4}

In subdirectory `example/`, run

    ../mcnp6extract.x argon39_example.inp

### Output ###    {#outp4}

* `MCNP6extract_MODELNAME_tally0.out` ... tally of neutron flux out of the calculation domain; should be all 0s.
* `MCNP6extract_MODELNAME_tally1.out` ... tally of (n,p) <sup>39</sup>Ar yield per neutron per wt.frac. K (first 3 columns) and corresponding statistical uncertainty (next 3 columns)
* `MCNP6extract_MODELNAME_tally2.out` ... tally of (n,&alpha;) <sup>21</sup>Ne yield per neutron per wt.frac. Mg (first 3 columns) and corresponding statistical uncertainty (next 3 columns)

# Input Data #    {#inputdata}

_to be completed..._

