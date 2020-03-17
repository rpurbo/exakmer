# Comm Pattern:
# 1) broadcast (genomes to load?)
# 2) scatter (put data on active rank)
# 3) gather (only active rank put data, other send null)
# 4) broadcast data to every rank
# 5) move to next rank => 2

import MPI

MPI.Init()

comm = MPI.COMM_WORLD
root = 0
host = gethostname()
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)
banner = "[rank" * string(rank) * ":" * string(host) * "] "
recv = 0
msg = 42
tag = 10

# 
# set active workers except root, i=1 
# 
for active in 2:size

	# INIT BUFFERS

	if(rank == root)
		send_mesg = Array{Int8}(undef, size)
		recv_mesg = Array{Int8}(undef, size)
		bcast_mesg = Array{Int8}(undef, 1)

		fill!(send_mesg,0)
		send_mesg[active] = 1
	else
		buff = Array{Int8}(undef, 1)
	end

	# === SCATTER DATA TO WORKERS ===
	isActive = 0
	if(rank == root)
		MPI.Scatter_in_place!(send_mesg,1,root,comm)
		for i in 2:(size)
			out = string(i-1) * " " * string(send_mesg[i])
			# println(out)
		end 
	else
		buff = Array{Int8}(undef, 1)
		MPI.Scatter_in_place!(buff,1,root,comm)
		isActive = buff[1]
	end

	MPI.Barrier(comm)


	# === GATHER DATA FROM ACTIVE WORKER ===
	if(rank == root)
		MPI.Gather_in_place!(recv_mesg,1,root,comm)
		bcast_mesg = recv_mesg[active]
		println(recv_mesg)
	else
		buff = Array{Int8}(undef, 1)
		if(isActive == 1)
			println(banner,"is active!")
			buff[1] = rank
		else
			buff[1] = 0
		end

		MPI.Gather_in_place!(buff,1,root,comm)
	end

	# === BROADCAST DATA TO OTHER WORKERS
	if(rank == root)
		MPI.bcast(bcast_mesg, root, comm)
	else
		recv_mesg = Array{Int8}(undef, 1)
		recv_mesg = MPI.bcast(C_NULL, root, comm)
		# println(recv_mesg)
	end
	

end











