# TItle

to run the pps, make an shell script to do:

   - find the pps.trees files
     - make a folder
     - move trees file to folder
     - make xml files for 100 trees from the trees file
     - make a slurm script that would run beast sequentially on each of these xml files. Remember to request sufficient time.
     - go up and restart loop

   - Then submit all slurm scripts

I also need to run the mismatched modes analyses. For this we need the trees used in each analyses. 

  - bd cc --analyses running
  - bd ce
  - cc ce
  - cc bd
  - ce cc
  - ce bd

