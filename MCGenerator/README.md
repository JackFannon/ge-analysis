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
