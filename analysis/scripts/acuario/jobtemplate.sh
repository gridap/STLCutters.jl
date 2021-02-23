#!/bin/bash
#SBATCH --job-name={{{name}}}
#SBATCH --output={{{o}}}
#SBATCH --error={{{e}}}
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks=1
##SBATCH --mem-per-cpu={{mem}} #MB

##SBATCH --time={{walltime}}

{{{julia}}} --project={{{projectdir}}} \
      -O3 --check-bounds=no -e\
      'using STLCutters;
       using STLCutters.Tests;
       {{{includes}}}
       {{func}}({{{args}}})'
