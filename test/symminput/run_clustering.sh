module load python3/pyrmsd/4.1.gita558b8a
module load imp-fast/last_ok_build

export top_dir=/home/rakesh/WORK/Projects/3-SPOTSComplex/IMPModeling/AnalysisModules/imp-sampcon/test/symminput
export analys_dir=$top_dir/analys/
export name=SPOTS

python /home/rakesh/WORK/Projects/3-SPOTSComplex/IMPModeling/AnalysisModules/imp-sampcon/pyext/src/Master_Sampling_Exhaustiveness_Analysis.py --sysname $name --path $analys_dir --ambiguity symm_groups.txt --mode cuda --cores 4 --align --density density.txt --gridsize 5.0 --gnuplot --scoreA A_models_clust-1.txt --scoreB B_models_clust-1.txt --rmfA A_models_clust-1.rmf3 --rmfB B_models_clust-1.rmf3
