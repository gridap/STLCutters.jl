#!/bin/bash
#PBS -q {{q}} 
#PBS -l walltime={{walltime}}
#PBS -l ncpus={{ncpus}}
#PBS -l mem={{mem}}
#PBS -N {{{name}}}
#PBS -l wd
#PBS -o {{{o}}}
#PBS -e {{{e}}} 

module load julia/1.5.3

julia --project={{{projectdir}}} \
      -O3 --check-bounds=no -e\
      'using STLCutters;
       using STLCutters.Tests;
       {{{includes}}}
       {{func}}({{{args}}})'
