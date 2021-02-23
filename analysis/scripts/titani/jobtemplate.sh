#!/bin/bash

#SBATCH --partition={{q}}
#SBATCH --job-name={{{name}}}
#SBATCH --output={{{o}}}
#SBATCH --error={{{e}}}
#SBATCH --mem={{mem}} #MB
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

{{{julia}}} --project={{{projectdir}}} \
      -O3 --check-bounds=no -e\
      'using STLCutters;
       using STLCutters.Tests;
       {{{includes}}}
       {{func}}({{{args}}})'
