
export UCX_WARN_UNUSED_ENV_VARS=n
ulimit -s 10240
mpiexec -np 7 julia exakmer.jl 2>&1 | tee log
#mpiexec -np 7 julia exakmer.jl 2>&1 > log
