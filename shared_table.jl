using BioSequences
using FASTX
using GZip
using CodecZlib
using Dictionaries
using DataStructures

mutable struct composite_kmer_table
	session_cnt::Int8
	genome_cnt::Int8
	curr_idx::Int8
	genome2idx::Dict{String, UInt64}
	idx2genome::Dict{UInt64, String}
	table::Dict{UInt64, BitArray}

	function composite_kmer_table(member_size::Int)
		session_cnt = 0
		genome_cnt = member_size
		curr_idx = 1
		genome2idx = Dict{String, UInt64}()
		idx2genome = Dict{UInt64, String}()
		table = Dict{UInt64, BitArray}() # save DNAMer as UInt64
		new(session_cnt, genome_cnt, curr_idx, genome2idx, idx2genome, table)
	end
end


function main()
	KMER_SIZE =27
	input="list2.txt"
	genomes = Dict{String, String}()
	genome2arr = Dict{String, UInt8}()

	idx = 1
	open(input, "r") do file
		for line in eachline(file)
			arr = split(line," ")
			genomes[arr[1]] = arr[2]
			genome2arr[arr[1]] = idx
			idx = idx + 1
		end
	end

	genome_tot = idx - 1

	# load fasta file and add them to shared table
	db = composite_kmer_table(genome_tot)
	for (name, path) in genomes
		idx =  genome2arr[name]
		merdb = Dict{DNAMer{KMER_SIZE}, Int8}()
		load_kmer!(path, merdb,KMER_SIZE)
		CKT_add_kmer_table!(db, name, merdb)
		merdb = nothing
	end

	gens = CKT_get_genomes(db)
	test = gens[1]

	merdb = Dict{DNAMer{KMER_SIZE}, Int8}()
	CKT_get_genome_kmer_table!(db, merdb, test)

	#for (mer,count) in merdb
	#	println(mer)
	#end
	
	arr = UInt64[]
	CKT_get_encoded_genome_kmer_arr!(db, arr, test)
	#for kmer_int in arr
	#	println(kmer_int,"\t", DNAMer{27}(kmer_int))
	#end	

end

function CKT_get_genome_kmer_table!(table::composite_kmer_table, kmer_table::Dict{DNAMer{27}, Int8}, name::String)
	idx = table.genome2idx[name]
	for kmer_int in collect(keys(table.table))
		hasKmer = table.table[kmer_int][idx]

		if hasKmer == true
			key = DNAMer{27}(kmer_int)
			kmer_table[key] = 1
		end
	end
	
end

function CKT_get_encoded_genome_kmer_arr!(table::composite_kmer_table, arr::Array{UInt64}, name::String)
	idx = table.genome2idx[name]
	for kmer_int in collect(keys(table.table))
		hasKmer = table.table[kmer_int][idx]

		if hasKmer == true
			push!(arr, kmer_int)
		end
	end
end


function CKT_get_genomes(table::composite_kmer_table)
	return collect(keys(table.genome2idx))
end

function CKT_add_kmer_table!(table::composite_kmer_table, name::String, kmer_table::Dict{DNAMer{27}, Int8})
	idx = table.curr_idx

	table.genome2idx[name] = idx
	table.idx2genome[idx] = name

	for (kmer,cnt) in kmer_table
		key = BioSequences.encoded_data(kmer)
		
		if(haskey(table.table, key))
			bits = table.table[key]
			bits[idx] = true
			table.table[key] = bits
		else
			bits = falses(table.genome_cnt)
			bits[idx] = true
			table.table[key] = bits
		end

	end

	table.curr_idx += 1
end


function encode_kmer_table!(db::Dict{DNAMer{27}, Int8}, arr::Array{UInt64})
	idx=1
	for (mer, count) in db
		enc_mer = BioSequences.encoded_data(mer)
		arr[idx] = enc_mer
		idx += 1		
	end	
end

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

function tag_kmer!(refdb::Dict{DNAMer{27}, Int8}, subjdb::Dict{DNAMer{27}, Int8})

	for (mer, count) in refdb
		subjdb[mer] = 0
	end

end


main()
