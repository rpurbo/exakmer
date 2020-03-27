import MPI
using BioSequences
using FASTX
using GZip
using CodecZlib
using Dictionaries
using DataStructures
using Printf

mutable struct composite_kmer_table
	kmer_size::Int8
	occ_idx::Int8
	tag_idx::Int8
	session_cnt::Int8
	genome_cnt::Int8
	curr_idx::Int8
	genome2idx::Dict{String, UInt64}
	idx2genome::Dict{UInt64, String}
	table::Dict{UInt64, BitArray}

	function composite_kmer_table(member_size::Int, kmer_size::Int)
		occ_idx = 1
		tag_idx = 2
		session_cnt = 0
		genome_cnt = member_size
		curr_idx = 1
		genome2idx = Dict{String, UInt64}()
		idx2genome = Dict{UInt64, String}()
		table = Dict{UInt64, BitArray}() # save DNAMer as UInt64
		new(kmer_size, occ_idx, tag_idx, session_cnt, genome_cnt, curr_idx, genome2idx, idx2genome, table)
	end
end

mutable struct genome_list
	kmer_size::Int
	size::Int
	fasta::Dict{UInt8, String}
	names::Dict{UInt8, String}
	partition::Dict{UInt8, UInt8}
end

#########################################
# Main Function
# Comment:
# NEXT: 
# 	1) create benchmark/baseline output w/ perl
# 	2) test on PBS cluster
# 	3) test on more genomes
#	4) benchmarking
#	5) aggresive mode = instead of ring comm, use all-to-all comm
#########################################

function main()
	# Init parameters
	file = "list.txt"
	ksize = 27
	outdir = "/gpfs1/scratch/admin/projects/julia/repos/output"

	##########################################################################
	## STAGE 1
	## Set up MPI comms and load genomes as kmer table to memory
        ##########################################################################

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
	ckt=composite_kmer_table(0, ksize)
	if(mpi_rank != root) 
		genomes_num = length(keys(gpath))
		ckt = composite_kmer_table(genomes_num, ksize)
		
		for (name,path) in gpath
			# println("",banner," ", name, " ",path)
			gtable = Dict{NucleotideSeq, Int8}()
			load_kmer!(path, gtable, kmer_size)
			CKT_add_kmer_table!(ckt, name, gtable)
			gtable = nothing
		end	
	end

	# Barrier to ensure all workers loaded their assigned genomes
	MPI.Barrier(comm)

	##########################################################################
	## STAGE 2
	## Run the scatter & gather communication
	## 1) root assign active worker (consecutively based on rank) - scatter
	## 2) active worker send its kmer table to root
	## 3) root broadcast kmer table - broadcast
	## 4) workers tag the received kmer table - broadcast
	##########################################################################

	# For all non-root workers. For mpi_size cluster, root(rank=0)'s data is in array[1].
	# workers(rank=i)'s data is in array[i+1].
	for active in 2:mpi_size

		# Stage 2.1 - Init communication buffers
		if(mpi_rank == root)
			send_mesg = Array{Int8}(undef, mpi_size)
			recv_mesg = Array{Int8}(undef, mpi_size)
			bcast_mesg = Array{Int8}(undef, 1)

			fill!(send_mesg,0)
			send_mesg[active] = 1
		else
			buff = Array{Int8}(undef, 1)
		end

		# Stage 2.2 - Root send active assignment to workers consecutively
		isActive = 0
		if(mpi_rank == root)
			MPI.Scatter_in_place!(send_mesg,1,root,comm)
		else
			buff = Array{Int8}(undef, 1)
			MPI.Scatter_in_place!(buff,1,root,comm)
			isActive = buff[1]
		end

		# Stage 2.3 - Active worker send the kmer table to root directly
		ckt_arr = UInt64[]
		active_rank = active-1

		if(mpi_rank == root)
			# Get incoming kmer table size
			buff = Array{Int64}(undef, 1)
			MPI.Recv!(buff, active_rank, 1, comm)
			arr_size = buff[1]

			# Receive kmer table
			ckt_arr = Array{UInt64}(undef, arr_size)
			MPI.Recv!(ckt_arr, active_rank, 2, comm)

			
			#println(string(active_rank),":",string(length(ckt_arr)))
			#CHECKPOINT - root received correct numbers of kmer from workers
		else
			if(isActive == 1)
				# Send the size of kmer table for memory allocation
				CKT_get_all_encoded_genome_kmer_arr!(ckt,ckt_arr)
				arr_size = length(ckt_arr)
				buff = Array{Int64}(undef, 1)
				buff[1] = arr_size		
				MPI.Send(buff, root, 1, comm)

				# Send the kmer table
				MPI.Send(ckt_arr, root, 2, comm)
			end
		end

		# Stage 2.4 - Broadcast kmer table to all workers (except active)
		recv_arr = nothing
		if(mpi_rank == root)
			arr_size = length(ckt_arr)
			MPI.bcast(arr_size, root, comm)
			MPI.bcast(ckt_arr, root, comm)
			#println(string(active_rank),":",string(length(ckt_arr)))
			#CHECKPOINT - root broadcasted correct numbers of kmer for workers

		else	
			arr_size = MPI.bcast(C_NULL, root, comm)
			#println("from ",string(active_rank)," to ",string(mpi_rank),":",string(arr_size))

			recv_arr = Array{UInt64}(undef, arr_size)
			recv_arr = MPI.bcast(C_NULL, root, comm)
			#println("from ",string(active_rank)," to ",string(mpi_rank),":",string(length(recv_arr)))
			#CHECKPOINT - workers received correct numbers of kmer for workers


		end

		# Stage 2.5 - Tag non-active workers's CKT kmer table using the received kmers
		if(mpi_rank != root)
			if(isActive == 0)
				#println("from ",string(active_rank)," to ",string(mpi_rank),":",string(length(recv_arr)))
				CKT_tag_kmer_arr!(ckt, recv_arr)
				#println("from ",string(active_rank)," to ",string(mpi_rank),":",string(length(collect(keys(ckt.table)))))
				#CHECKPOINT - both recv_arr and ckt are intact
			end
		end

	end

	
	#println("from ",string(mpi_rank),":",string(length(collect(keys(ckt.table)))))
	#CHECKPOINT - ckt is intact

	##########################################################################
	## STAGE 3 
	## Write un-tagged (unique) kmers to file
	##########################################################################
	if(mpi_rank != root)
		genomes = CKT_get_genomes(ckt)

		#println("from ",string(mpi_rank),":",string(length(collect(keys(ckt.table)))))

		for genome in genomes
			#println("from ",string(mpi_rank),":",genome)
			outfile = outdir * "/" * genome * ".uniq"
			io = open(outfile, "w")

			gtable = Dict{NucleotideSeq, Int8}()
			CKT_get_genome_kmer_table_untagged!(ckt, gtable, genome)

			#println("from ",string(mpi_rank),":",string(length(collect(keys(gtable)))))

			for (kmer,count) in gtable
				println(io,string(kmer))
			end

			close(io)
		end
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

#########################################
# function: load_kmer!
# Read fasta from path and load all the kmers to a dictionary
# comment: need to check if it's zipped or is in fasta format
#########################################
function load_kmer!(path::String, db::Dict{NucleotideSeq, Int8}, kmer_size::Int)
	ON::Int8 = 1
	reader = FASTA.Reader(GzipDecompressorStream(open(path)))

	for record in reader
		for iter in each(DNAMer{kmer_size}, FASTX.FASTA.sequence(record))
			db[iter.fw] = ON
			db[iter.bw] = ON
		end
	end
	close(reader)
end

#########################################
# function: CKT_tag_kmer!
# Set the column tag_idx of the composite kmer table to
# true, indicating that the kmer occurs in other genomes
#########################################
function CKT_tag_kmer!(table::composite_kmer_table, kmer::NucleotideSeq)
	kmer_int = BioSequences.encoded_data(kmer)
	table.table[kmer_int][:, table.tag_idx] .= true
end

#########################################
## function: CKT_tag_kmer_arr!
## Set the column tag_idx of the composite kmer table to
## true from an encoded kmer table, indicating that the kmer occurs in other genomes
##########################################
function CKT_tag_kmer_arr!(table::composite_kmer_table, arr::Array{UInt64})
	for kmer_int in arr
		if(haskey(table.table, kmer_int))
			table.table[kmer_int][:, table.tag_idx] .= true
		end
	end	
end

#########################################
## function: CKT_get_all_encoded_genome_kmer_arr!
## Get all kmer table from CKT data structure
## in encoded (UInt64) data type. The resulting
## array can be sent via MPI comm.
##########################################
function CKT_get_all_encoded_genome_kmer_arr!(table::composite_kmer_table, arr::Array{UInt64})
	for kmer_int in collect(keys(table.table))
		push!(arr, kmer_int)
	end
end


#########################################
# function: CKT_get_genome_kmer_table!
# Get a genome's kmer table from CKT data structure
# comment: no error-handling on non-existing genome
#########################################
function CKT_get_genome_kmer_table_untagged!(table::composite_kmer_table, kmer_table::Dict{NucleotideSeq, Int8}, name::String)
	idx = table.genome2idx[name]

	for kmer_int in collect(keys(table.table))
		hasKmer = table.table[kmer_int][idx, table.occ_idx]
		isTagged = table.table[kmer_int][idx, table.tag_idx]

		if(hasKmer == true && isTagged == false)
			key = DNAMer{table.kmer_size}(kmer_int)
			kmer_table[key] = 1
			nottag += 1	
		end

	end

end

#########################################
# function: CKT_get_genome_kmer_table!
# Get a genome's kmer table from CKT data structure
# comment: no error-handling on non-existing genome 
#########################################
function CKT_get_genome_kmer_table!(table::composite_kmer_table, kmer_table::Dict{NucleotideSeq, Int8}, name::String)
	idx = table.genome2idx[name]
	for kmer_int in collect(keys(table.table))
		hasKmer = table.table[kmer_int][idx, table.occ_idx]

		if hasKmer == true
			key = DNAMer{table.kmer_size}(kmer_int)
			kmer_table[key] = 1
		end
	end
end

#########################################
# function: CKT_get_encoded_genome_kmer_arr!
# Get a genome's kmer table from CKT data structure
# in encoded (UInt64) data type. The resulting
# array can be sent via MPI comm.
#########################################
function CKT_get_encoded_genome_kmer_arr!(table::composite_kmer_table, arr::Array{UInt64}, name::String)
	idx = table.genome2idx[name]
	for kmer_int in collect(keys(table.table))
		hasKmer = table.table[kmer_int][idx, table.occ_idx]

		if hasKmer == true
			push!(arr, kmer_int)
		end
	end
end

#########################################
# function: CKT_get_genomes
# Get all genome names contained in the CKT data structure
#
#########################################
function CKT_get_genomes(table::composite_kmer_table)
	return collect(keys(table.genome2idx))
end

#########################################
# function: CKT_print_table
# Print CKT data structure for troubleshooting
#
#########################################
function CKT_print_table(table::composite_kmer_table)
	names = collect(keys(table.genome2idx))
	@printf("#kmer")
	for name in names
		@printf(" %s", name)
	end
	println("")

	for kmer_int in collect(keys(table.table))
		@printf("%s", string(DNAMer{table.kmer_size}(kmer_int)))
		for name in names
			idx = table.genome2idx[name]
			occ = table.table[kmer_int][idx, table.occ_idx]
			tag = table.table[kmer_int][idx, table.tag_idx]
			@printf(" %s/%s",string(occ),string(tag))
		end
		println("")
	end
end

#########################################
# function: CKT_add_kmer_table!
# Add a genome's kmer table to the CKT data structure.
#
#########################################
function CKT_add_kmer_table!(table::composite_kmer_table, name::String, kmer_table::Dict{NucleotideSeq, Int8})
	idx = table.curr_idx
	table.genome2idx[name] = idx
	table.idx2genome[idx] = name

	for (kmer,cnt) in kmer_table
		key = BioSequences.encoded_data(kmer)

		if(haskey(table.table, key))
			bits = table.table[key]
			bits[idx, table.occ_idx] = true
			bits[idx, table.tag_idx] = false
			table.table[key] = bits
		else
			bits = falses(table.genome_cnt,2)
			bits[idx, table.occ_idx] = true
			bits[idx, table.tag_idx] = false
			table.table[key] = bits
		end
	end
	table.curr_idx += 1
end


main()

