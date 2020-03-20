import MPI
using BioSequences
using FASTX
using GZip
using CodecZlib
using Dictionaries
using DataStructures


mutable struct composite_kmer_table
	occ_idx::Int8
	tag_idx::Int8
	session_cnt::Int8
	genome_cnt::Int8
	curr_idx::Int8
	genome2idx::Dict{String, UInt64}
	idx2genome::Dict{UInt64, String}
	table::Dict{UInt64, BitArray}

	function composite_kmer_table(member_size::Int)
		occ_idx = 1
		tag_idx = 2
		session_cnt = 0
		genome_cnt = member_size
		curr_idx = 1
		genome2idx = Dict{String, UInt64}()
		idx2genome = Dict{UInt64, String}()
		table = Dict{UInt64, BitArray}() # save DNAMer as UInt64
		new(occ_idx, tag_idx, session_cnt, genome_cnt, curr_idx, genome2idx, idx2genome, table)
	end
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
	# steps
	# 1) root: scatter active worker
	# 2) active worker: encoded kmer, send back via gather. passive worker, send null
	# 3) root: broadcast encoded kmers
	# 4) workers: tag kmers (dont self-tag)
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

#########################################
# function: encode_kmer_table!
# Encode kmer table dictionary to array of UInt64 (for MPI comm)
#########################################  
function encode_kmer_table!(db::Dict{DNAMer{27}, Int8}, arr::Array{UInt64})
	idx=1
	for (mer, count) in db
		enc_mer = BioSequences.encoded_data(mer)
		arr[idx] = enc_mer
		idx += 1
	end
end


function CKT_tag_kmer!(table::composite_kmer_table, kmer::DNAMer{27})
	kmer_int = BioSequences.encoded_data(kmer)
	table.table[kmer_int][:, table.tag_idx] .= true
end


function CKT_get_genome_kmer_table!(table::composite_kmer_table, kmer_table::Dict{DNAMer{27}, Int8}, name::String)
	idx = table.genome2idx[name]
	for kmer_int in collect(keys(table.table))
		hasKmer = table.table[kmer_int][idx, table.occ_idx]

		if hasKmer == true
			key = DNAMer{27}(kmer_int)
			kmer_table[key] = 1
		end
	end
end

function CKT_get_encoded_genome_kmer_arr!(table::composite_kmer_table, arr::Array{UInt64}, name::String)
	idx = table.genome2idx[name]
	for kmer_int in collect(keys(table.table))
		hasKmer = table.table[kmer_int][idx, table.occ_idx]

		if hasKmer == true
			push!(arr, kmer_int)
		end
	end
end


function CKT_get_genomes(table::composite_kmer_table)
	return collect(keys(table.genome2idx))
end

function CKT_print_table(table::composite_kmer_table)
	names = collect(keys(table.genome2idx))
	@printf("#kmer")
	for name in names
		@printf(" %s", name)
	end
	println("")

	for kmer_int in collect(keys(table.table))
		@printf("%s", string(DNAMer{27}(kmer_int)))
		for name in names
			idx = table.genome2idx[name]
			occ = table.table[kmer_int][idx, table.occ_idx]
			tag = table.table[kmer_int][idx, table.tag_idx]
			@printf(" %s/%s",string(occ),string(tag))
		end
		println("")
	end
end

function CKT_add_kmer_table!(table::composite_kmer_table, name::String, kmer_table::Dict{DNAMer{27}, Int8})
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

