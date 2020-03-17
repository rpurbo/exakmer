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
	# comment: need to check if all files are readable
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
		kmer_size = nothing
		msg_buf = nothing
		gpath = Dict{String, String}()
		gname2idx = Dict{String, UInt8}() 
		gidx2name = Dict{UInt8, String}()
		kmer_dicts = nothing
	end

	# Broadcast genome list to workers
	# comment: no option but to send everything as scatter doesnt work for Strings yet	
	if(mpi_rank == root)
		msg_buf = MPI.bcast(glist, root, comm)
	else
		msg_buf = MPI.bcast(C_NULL, root, comm)
		kmer_size = msg_buf.kmer_size		

		idx = 1
		for (id,wid) in msg_buf.partition
			if(wid == mpi_rank)
				name = msg_buf.names[id]
				path = msg_buf.fasta[id]
				gpath[name] = path
				idx += 1
			end
		end

		msg_buf = nothing
	end

	# On workers, read the genome fasta files and load them to tables
	if(mpi_rank != root) 
		genomes_num = length(keys(gpath))
		kmer_dicts = Array{ Dict{DNAMer{kmer_size}, Int8}}(undef,genomes_num)
		
		idx = 1
		for (name,path) in gpath
			# println("",banner," ", name, " ",path)
			gtable = Dict{DNAMer{kmer_size}, Int8}()
			load_kmer!(path, gtable, kmer_size)
			kmer_dicts[idx] = gtable

			gname2idx[name] = idx
			gidx2name[idx] = name
			idx += 1	
		end	
	end

	# Barrier to ensure all workers loaded their assigned genomes
	MPI.Barrier(comm)
	sleep(10)

	# Next, 1) encode kmers to UInt64, 2) scatter & gather kmers from worker, 3) broadcast the kmers. 
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

#########################################
# function: load_kmer!
# Read fasta from path and load all the kmers to a dictionary
# comment: need to check if it's zipped or is in fasta format
#########################################
function load_kmer!(path::String, db::Dict{DNAMer{27}, Int8}, kmer_size::Int)
	ON::Int8 = 1
	reader = FASTA.Reader(GzipDecompressorStream(open(path)))

	for record in reader
		for iter in each(DNAMer{27}, FASTX.FASTA.sequence(record))
			db[iter.fw] = ON
			db[iter.bw] = ON
		end
	end
	close(reader)
end

#########################################
# function: tag_kmer!
# Tag kmer (set to 0) from a kmer dictionary. Need to explore if can be done
# with automatic memory recycling
#########################################
function tag_kmer!(refdb::Dict{DNAMer{27}, Int8}, subjdb::Dict{DNAMer{27}, Int8})
	for (mer, count) in refdb
		subjdb[mer] = 0
	end

end


main()

