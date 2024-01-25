# qual2kw-ver5
## QUAL2KW version 5

This is the Fortran source code for version 5 of QUAL2Kw. 

It was slightly modified to get it working on Linux.

The original source code is available at https://github.com/gjpelletier/qual2kw-ver5/blob/master/qual2kw51b52_f.zip

The Windows executable version is available at
https://ecology.wa.gov/Research-Data/Data-resources/Models-spreadsheets/Modeling-the-environment/Models-tools-for-TMDLs


### Building

Use the Fortran Package Manager (FPM) for building:

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