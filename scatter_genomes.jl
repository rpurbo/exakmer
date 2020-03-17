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
	fasta::Dict{String, String}
	gidx::Dict{String, UInt8}
end


#########################################
# Main Function
#########################################

function main()
	file = "list.txt"
	ksize = 27

	#########################################
	# Init MPI
	#########################################

	MPI.Init()
	comm = MPI.COMM_WORLD
	root = 0
	host = gethostname()
	rank = MPI.Comm_rank(comm)
	size = MPI.Comm_size(comm)
	banner = "[rank" * string(rank) * ":" * string(host) * "] "

	if(rank == root)	
		glist = get_genome_list(file,ksize)
	end


	
	# For all workers
	for active in 2:size
		# Init buffers
		if(rank == root)
			send_mesg = Array{Char}(undef, size)
			recv_mesg = Array{Int8}(undef, size)
			fill!(send_mesg,0)
			send_mesg[active] = 'A'	
		else
			buff = Array{String}(undef, 1)
		end

		# Scatter data to workers	
		isActive = ""
		if(rank == root)
			MPI.Scatter_in_place!(send_mesg,1,root,comm)
		else
			buff = Array{Char}(undef, 1)
			MPI.Scatter_in_place!(buff,1,root,comm)
			isActive = buff[1]
		end
	
		if(isActive == 'A')
			println("",banner,"is active")
		end

	end


end



#########################################
# Get Genome List
########################################

function get_genome_list(input::String,kmer_size::Int)
	genomes = Dict{String, String}()
	genome2arr = Dict{String, UInt8}()

	idx = 1
	open(input, "r") do file
		for line in eachline(file)
			arr = split(line," ")
			genomes[arr[2]] = arr[3]
			genome2arr[arr[2]] = parse(UInt8,arr[1])
			idx = idx + 1
		end
	end
	genome_tot = idx - 1

	glist = genome_list(kmer_size, genome_tot, genomes, genome2arr)
	return glist
end


main()

