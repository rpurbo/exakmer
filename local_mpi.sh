
export UCX_WARN_UNUSED_ENV_VARS=n
ulimit -s 10240
mpiexec -np 6 julia bcast_load_genomes.jl 2>&1 | tee log

