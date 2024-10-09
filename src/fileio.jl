"""
Preprocesses the barcode file by modifying sequences based on specific criteria.
"""
function preprocess_bc_file(bc_file::String, complement::Bool)
	mime_type = readchomp(`file --mime-type -b $bc_file`)
	delim = "\t"
	if occursin(r"\.csv$", bc_file) || mime_type == "text/csv"
		delim = ","
	end
	bc_df = CSV.read(bc_file, DataFrame, delim = delim)
	for i in 1:nrow(bc_df)
		prefix_region = 1:findfirst('B', bc_df.Full_annotation[i])-1
		suffix_region = findlast('B', bc_df.Full_annotation[i])+1:length(bc_df.Full_annotation[i])
		prefix = SubString(bc_df.Full_seq[i], prefix_region)
		suffix = SubString(bc_df.Full_seq[i], suffix_region)
		bc_df.Full_seq[i] = SubString(bc_df.Full_seq[i], length(prefix)+1)
		bc_df.Full_seq[i] = SubString(bc_df.Full_seq[i], 1, length(bc_df.Full_seq[i])-length(suffix))
	end
	bc_df.Full_seq = uppercase.(bc_df.Full_seq)
	bc_df.Full_seq = replace.(bc_df.Full_seq, "U" => "T")
	if complement == true
		bc_df.Full_seq = reverse.(bc_df.Full_seq)
		bc_df.Full_seq = replace.(bc_df.Full_seq, "A" => "T", "T" => "A", "G" => "C", "C" => "G")
	end
	return bc_df
end

function preprocess_fastq(file_path::String)
    # Decompress if gzipped
    mine_type = readchomp(`file --mime-type -b $file_path`)
	is_gzipped = mine_type == "application/x-gzip"
    if is_gzipped
        decompressed_path = replace(file_path, r"\.gz$" => "")
        run(`gunzip -c $file_path \> $decompressed_path`)
    else
        decompressed_path = file_path
    end

    # Extract the prefix (part of the filename before .fastq)
    prefix = replace(basename(decompressed_path), r"\.fastq$" => "")
    return decompressed_path, prefix
end

"""
Divides a pair of FASTQ files into smaller parts for parallel processing.
It calculates the number of reads per worker and uses the split command to divide the files.
"""
function divide_fastq(FASTQ_file1::String, FASTQ_file2::String, output_dir::String, workers::Int)
	divided_dir = mktempdir(output_dir)
	total_lines = countlines(FASTQ_file1)
	total_reads = total_lines ÷ 4
	reads_per_worker = cld(total_reads, workers)
	lines_per_worker = reads_per_worker * 4
	run(`split -l $lines_per_worker -a 5 -d $FASTQ_file1 $divided_dir/file1_ --additional-suffix=.fastq`)
	run(`split -l $lines_per_worker -a 5 -d $FASTQ_file2 $divided_dir/file2_ --additional-suffix=.fastq`)
	return divided_dir
end

"""
Divides a single FASTQ file for parallel processing.
"""
function divide_fastq(FASTQ_file::String, output_dir::String,workers::Int)
	divided_dir = mktempdir(output_dir)
	total_lines = countlines(FASTQ_file)
	total_reads = total_lines ÷ 4
	reads_per_worker = cld(total_reads, workers)
	lines_per_worker = reads_per_worker * 4
	run(`split -l $lines_per_worker -a 5 -d $FASTQ_file $divided_dir/ --additional-suffix=.fastq`)
	return divided_dir
end

function multi_demultiplex(thread_num::Int, bc_df::DataFrame, divided_dir::String, output_prefix1::String, output_prefix2::String, max_error_rate::Float64, min_delta::Float64, mismatch::Int, indel::Int, classify_both::Bool)
	FASTQ_file1 = divided_dir * "/file1_" * lpad((thread_num - 1), 5, "0") * ".fastq"
	FASTQ_file2 = divided_dir * "/file2_" * lpad((thread_num - 1), 5, "0") * ".fastq"
	mkdir(divided_dir * "/thread" * string(thread_num))
	output_dir = divided_dir * "/thread" * string(thread_num)
	classify_sequences(FASTQ_file1, FASTQ_file2, bc_df, output_dir, output_prefix1, output_prefix2, max_error_rate, min_delta, mismatch, indel, classify_both)
end

function multi_demultiplex(thread_num::Int, bc_df::DataFrame, divided_dir::String, output_prefix::String, max_error_rate::Float64, min_delta::Float64, mismatch::Int, indel::Int)
	FASTQ_file = divided_dir * "/" * lpad((thread_num - 1), 5, "0") * ".fastq"
	mkdir(divided_dir * "/thread" * string(thread_num))
	output_dir = divided_dir * "/thread" * string(thread_num)
	classify_sequences(FASTQ_file, bc_df, output_dir, output_prefix, max_error_rate, min_delta, mismatch, indel)
end


function merge_fastq_files(paths::Vector, bc_df::DataFrame, output_dir::String, prefix::String)
	paths_unknown = filter(x -> occursin(r".*unknown.fastq", x), paths)
	if paths_unknown != []
		run(pipeline(`cat $paths_unknown`, stdout = "$output_dir/$prefix.unknown.fastq"))
	end

	paths_ambiguous_classification = filter(x -> occursin(r".*/ambiguous_classification.fastq", x), paths)
	if paths_ambiguous_classification != []
		run(pipeline(`cat $paths_ambiguous_classification`, stdout = "$output_dir/$prefix.ambiguous_classification.fastq"))
	end

	for i in 1:nrow(bc_df)
		regex = string(bc_df.ID[i]) * ".fastq"
		paths_matched = filter(x -> occursin(regex, x), paths)
		if paths_matched != []
			run(pipeline(`cat $paths_matched`, stdout = "$output_dir/$prefix.$(bc_df.ID[i]).fastq"))
		end
	end
end