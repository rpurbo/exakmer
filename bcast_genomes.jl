import MPI
using BioSequences
using FASTX
using GZip
using CodecZlib
using Dictionaries
using DataStructures


mutable struct genome_list
	kmer_size::Int
	size::Int
	fasta::Dict{UInt8, String}
	names::Dict{UInt8, String}
	partition::Dict{UInt8, UInt8}
end


#########################################
# Main Function
#########################################

function main()
	# Init parameters
	file = "list.txt"
	ksize = 27

	# Init MPI
	MPI.Init()
	comm = MPI.COMM_WORLD
	root = 0
	mpi_host = gethostname()
	mpi_rank = MPI.Comm_rank(comm)
	mpi_size = MPI.Comm_size(comm)
	mpi_last_worker = mpi_size-1
	banner = "[rank" * string(mpi_rank) * ":" * string(mpi_host) * "] "

	# Read genome fasta paths from list file
	if(mpi_rank == root)	
		glist = get_genome_list(file,ksize)

		# Assign genomes to workers. non-root worker id starts from 2
		# comment: currently the genomes are spread out fifo. more advanced way
		# is to distribute based on the size of total genomes per worker
		worker_id = 1	
		for (id,worker) in sort(glist.partition)
			glist.partition[id] = worker_id
			worker_id += 1
			if(worker_id > mpi_last_worker)
				worker_id = 1
			end
		end
	end

	# Init buffers
	if(mpi_rank == root)
		msg_buf = nothing
	else
		msg_buf = nothing
	end

	# Broadcast genome list to workers
	# comment: no option but to send everything as scatter doesnt work for Strings yet	
	if(mpi_rank == root)
		msg_buf = MPI.bcast(glist, root, comm)
	else
		msg_buf = MPI.bcast(C_NULL, root, comm)
		
		for (id,wid) in msg_buf.partition
			if wid == mpi_rank
				name = msg_buf.names[id]
				path = msg_buf.fasta[id]
				# println("",banner," ", name, " ",path)
			end
		end

		# Read genomes from files
		# load to kmer tables

	end

end



#########################################
# function: get_genome_list
# Read input file (name path) and package them as genome_list struct 
########################################

function get_genome_list(input::String,kmer_size::Int)
	fasta = Dict{UInt8, String}()
	names = Dict{UInt8, String}()
	parts = Dict{UInt8, UInt8}()	
	
	idx = 1
	open(input, "r") do file
		for line in eachline(file)
			arr = split(line," ")
			fasta[idx] = arr[2]
			names[idx] = arr[1]
			parts[idx] = 0
			idx = idx + 1
		end
	end
	genome_tot = idx - 1

	glist = genome_list(kmer_size, genome_tot, fasta, names, parts)
	return glist
end


main()

