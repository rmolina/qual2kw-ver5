# qual2kw-ver5
## QUAL2KW version 5

This is the Fortran source code for version 5 of QUAL2Kw. 

It was slightly modified to get it working on Linux.

The original source code is available at https://github.com/gjpelletier/qual2kw-ver5/blob/master/qual2kw51b52_f.zip

The Windows executable version is available at
https://ecology.wa.gov/Research-Data/Data-resources/Models-spreadsheets/Modeling-the-environment/Models-tools-for-TMDLs


### Building

This uses the Fortran Package Manager (FPM) for building.
See: https://fpm.fortran-lang.org/

```
(fortran) ruben@pop-os:~/github/qual2kw-ver5$ fpm build
nrtype.f90                             done.
phsolve.f90                            done.
classwaterquality.f90                  done.
classreach13.f90                       done.
classlightheat.f90                     done.
classdate.f90                          done.
ClassIntegralData.f90                  done.
classSolarPosition.f90                 done.
classDownstream.f90                    done.
systemparams.f90                       done.
classrivertopo13.f90                   done.
classhydraulics.f90                    done.
classheadwater.f90                     done.
classmeteo.f90                         done.
classstoch.f90                         done.
classSolarCalc.f90                     done.
ClassSrcIn.f90                         done.
classreadfile.f90                      done.
classoutput.f90                        done.
classintegral.f90                      done.
libqual2kw-ver5.a                      done.
q2kmain.f90                            done.
q2kmain                                done.
[100%] Project compiled successfully.
```

### Running

Our `message.dat` list the names of the input and output files:

```
(fortran) ruben@pop-os:~/github/qual2kw-ver5$ cat message.dat 
BC_1987-08-21.q2k
BC_1987-08-21.out
(fortran) ruben@pop-os:~/github/qual2kw-ver5$ ls -lh BC_1987-08-21.*
-rw-rw-r-- 1 ruben ruben 30K feb 26  2015 BC_1987-08-21.q2k
```

Here we use `BC_1987-08-21.q2k` as input file.
(This file is distributed with the Windows executable, but a copy is available in this repository.)

Run the model:

```
(fortran) ruben@pop-os:~/github/qual2kw-ver5$ fpm run
Project is up to date

 QUAL2Kw version 5.1
 Department of Ecology and Tufts University

 G.J. Pelletier, S.C. Chapra, and Hua Tao

 Program is running, please wait...

 Integrating: Euler method.
 day           1
 day           2
 day           3
 day           4
 day           5
 elapsed time:   0.314000010      seconds
```

Check for an output file:

```
(fortran) ruben@pop-os:~/github/qual2kw-ver5$ ls BC_1987-08-21.* -lh
-rw-rw-r-- 1 ruben ruben 3,3M ene 25 13:35 BC_1987-08-21.out
-rw-rw-r-- 1 ruben ruben  30K feb 26  2015 BC_1987-08-21.q2k
(fortran) ruben@pop-os:~/github/qual2kw-ver5$ head BC_1987-08-21.out 
 ** Hydraulics Summary **
Downstrea  Hydraulics          E'           H           B          Ac           U   trav time       slope  Reaeration      Reaeration formulas
 distance      Q,m3/s        m3/s           m           m          m2         mps           d                ka,20,/d               water/wind
  13.6000     0.71348     0.35674     0.20934    12.50000     2.61669     0.27267     0.00000     0.00400     0.00000                                   
  13.1750     1.47910     0.73955     0.32654    12.50000     4.08175     0.36237     0.01357     0.00400    11.83131     Specified/No wind             
  12.7500     1.49473     0.49824     0.32865    12.50000     4.10809     0.36385     0.02709     0.00400    11.76145     Specified/No wind             
  11.9000     1.52598     0.76299     0.33284    12.50000     4.16047     0.36678     0.05392     0.00400    11.62504     Specified/No wind             
  11.0500     1.55723     0.77861     0.33700    12.50000     4.21245     0.36967     0.08053     0.00400    11.49285     Specified/No wind             
  10.2000     1.58848     0.79424     0.34112    12.50000     4.26403     0.37253     0.10694     0.00400    11.36465     Specified/No wind             
   9.3500     2.20973     1.10486     0.43530    12.50000     5.44123     0.40611     0.13116     0.00350     8.50250     Specified/No wind             
(fortran) ruben@pop-os:~/github/qual2kw-ver5$ 
```
