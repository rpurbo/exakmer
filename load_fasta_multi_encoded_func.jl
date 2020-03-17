using BioSequences
using FASTX
using GZip
using CodecZlib
using Dictionaries
using DataStructures

function main()
	KMER_SIZE =27
	input="list2.txt"
	genomes = Dict{String, String}()
	genome2arr = Dict{String, UInt8}()

	idx = 1
	open(input, "r") do file
		for line in eachline(file)
			arr = split(line,"\t")
			genomes[arr[1]] = arr[2]
			genome2arr[arr[1]] = idx
			idx = idx + 1
		end
	end

	genome_tot = idx - 1
	db = Array{ Dict{DNAMer{KMER_SIZE}, Int8}}(undef,genome_tot)

	for (name, path) in genomes
		# println("loading ",name, " ", genome2arr[name])
		idx =  genome2arr[name]
		merdb = Dict{DNAMer{KMER_SIZE}, Int8}()
		load_kmer!(path, merdb,KMER_SIZE)
		db[idx] = merdb
	end

	#for (name, path) in genomes
	#	idx =  genome2arr[name]
	#	act_db = db[idx]
	#	for (mer, count) in act_db
	#		# println(name,"\t",mer,"\t",count)
	#		enc_mer = BioSequences.encoded_data(mer)
	#		dec_mer = DNAMer{27}(enc_mer)
	#		# println(name,"\t",mer,"\t",string(enc_mer),"\t",string(dec_mer))
	#		if(string(mer) != string(dec_mer))
	#			println(name,"\t",mer," doesnt matches")
	#		end
	#	end	
	#end

	act_db = db[1]
	arr = Array{UInt64}(undef,length(keys(act_db)))
	encode_kmer_table!(act_db, arr)	

	for i in 1:length(arr)
		str = string(DNAMer{27}(arr[i]))
		println(arr[i],"\t",str)
	end

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
