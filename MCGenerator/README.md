# How to use the MC Generator

## Generate cardfiles
### Create energy to momenta conversion files
Use the binary file `ene2mon` to generate energy to momenta conversions in keV
steps between two integer values. E.g.:

``` sh
./ene2mon 100 200
```

generates 100 files containing the momentum of an electron in 1 keV steps
between 100 keV and 200 keV. The output files are created in the `ene`
directory.

### Create cardfiles
Using the shell script `make_card.csh` the cardfiles, named `gedeo#####.card`,
that takes the momenta lines from the `ene` directory and creates a cardfile
that can be used to simulate electrons of that energy. One cardfile will be
created for each file in the `ene` directory.

## Input calibration output
Modify the `uhinit.f` file using the output from the calibration:

``` fortran
parameter (offset=#######)
parameter (slope=#######)
```

## Compile and execute gedeo
Compile `gedeo.f`, this is the main program for `uhinit.f`:

``` sh
make gedeo
```

execute `make_gedeo.csh`, this will generate jobscripts in the `script`
directory.

TODO rewrite make_gedeo for the new job submission system

## Edit and execute nqs.gedeo.sh

Make sure that the output directories are set properly:

``` sh
set f_err
set f_out
```

and that the energy range is appropriate for the range you are simulating:

``` sh
set i = #####
while ($i <- #####)
```

Executing this will submit all of the jobs in the `scripts` directory, and
output simulation in the form of `.hbk` files in the `hbk` directory. Use the
utility `h2root` to convert these to ROOT files. There is a script that does
this called `h2root.sh`.

## Fit the beam data
This step compares the beam data with the MC sample and finds the best fit for
the initial electron energy, before traversing through the Ti window at the end
of the beam pipe and the various components of the germanium detector.

### Modify fit_linac_data.C
The ROOT macro `fit_linac_data.C` uses the calibration constants from the first
part -- note that this ROOT macro uses the `po` and `p1` values from the fit,
not the slope and intercept that is used in a previous step -- in the array
`calib_const`. These need updating to be representative of the Ni calibration
file -- I think.

### Create the list of runs for a measurement of the energy scale.
A `.txt` file is used to contain the following information in columns:
- SK run number
- Approximate LINAC energy
- 0
- Approximate X-position
- Approximate Z-position
- Filename of the Ge data
- Lower bound of MC input electron energy
- Upper bound of MC input electron energy
- Lower bound of Ge energy deposit for chi2 calculation
- Upper bound of Ge energy deposit for chi2 calculation
- Attenuator flag - 1: without attenuator, 2: with attenuator
- Ignored and can be used for commenting

Any lines that begin with a `#` are also ignored by the macro

**Important note**
The first three columns need to be consistent with the information stored in
`/home/sklowe/linac/const/linac_sk#_runsum.dat`

### Fit the Ge spectrum and determine the best fit
