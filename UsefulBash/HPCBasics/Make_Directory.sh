#! /bin/bash -l

#SBATCH --partition=panda   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=make_dir
#SBATCH --time=00:60:00   # HH/MM/SS
#SBATCH --mem=32G   # memory requested, units available: K,M,G,T
#SBATCH --output make_dir-%j.out
#SBATCH --error make_dir-%j.err

source ~/.bashrc

echo "Starting at:" `date` >> make_dir_output.txt
echo "This is job #:" $SLURM_JOB_ID >> make_dir_output.txt
echo "Running on node:" `hostname` >> make_dir_output.txt
echo "Running on cluster:" $SLURM_CLUSTER_NAME >> make_dir_output.txt
echo "This job was assigned the temporary (local) directory:" $TMPDIR >> make_dir_output.txt

maindir=/projects/TF_$SLURM_JOB_ID
mkdir -p "$maindir" && cd "$maindir" || exit -1
touch test.txt

exit
