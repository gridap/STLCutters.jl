#!/bin/bash

#SBATCH --partition={{q}}
#SBATCH --job-name={{{name}}}
#SBATCH --output={{{o}}}
#SBATCH --error={{{e}}}
#SBATCH --mem={{mem}} #MB
#SBATCH --cpus-per-task={{ncpus}}
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

{{{julia}}} --project={{{projectdir}}} \
      -O3 --check-bounds=no --threads {{ncpus}} -e\
      'using STLCuttersAnalysis;
       {{{includes}}}
       {{func}}({{{args}}})'
