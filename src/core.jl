"""
Orchestrates the entire demultiplexing process for FASTQ files.
Handles the preprocessing, dividing, demultiplexing, and merging of files.
"""
function execute_demultiplexing(file_R1::String, file_R2::String, bc_file::String, output_dir::String; max_error_rate::Float64 = 0.22, min_delta::Float64 = 0.1, mismatch::Int = 1, indel::Int = 1, classify_both::Bool = false, bc_rev::Bool = true)
	if isdir(output_dir)
		error("Output directory already exists")
	end
	mkdir(output_dir)
	workers = nworkers()
	bc_df = preprocess_bc_file(bc_file, bc_rev)
	if workers == 1
		classify_sequences(file_R1, file_R2, bc_df, output_dir, max_error_rate, min_delta, mismatch, indel, classify_both)
	else
		divide_fastq(file_R1, file_R2, output_dir, workers)
		pmap(x -> mlt_demltplex(x, bc_df, output_dir, max_error_rate, min_delta, mismatch, indel, classify_both), 1:workers)
		paths = []
		for (root, dirs, files) in walkdir(output_dir)
			for file in files
				push!(paths, joinpath(root, file))
			end
		end
		if classify_both
			paths_R1 = filter(x -> occursin(r"R1", x), paths)
			paths_R2 = filter(x -> occursin(r"R2", x), paths)
			mkdir(output_dir * "/R1")
			mkdir(output_dir * "/R2")
			merge_fastq_files(paths_R1, bc_df, output_dir * "/R1")
			merge_fastq_files(paths_R2, bc_df, output_dir * "/R2")
		else
			merge_fastq_files(paths, bc_df, output_dir)
		end

		rm(joinpath(output_dir, "divided_fastq"), recursive = true)
		for i in 1:workers
			rm(joinpath(output_dir, "thread" * string(i)), recursive = true)
		end
	end
end

function execute_demultiplexing(file_R1::String, bc_file::String, output_dir::String; max_error_rate::Float64 = 0.2, min_delta::Float64 = 0.1, mismatch::Int = 1, indel::Int = 1, bc_rev::Bool = true)
	if isdir(output_dir)
		error("Output directory already exists")
	end
	mkdir(output_dir)

	workers = nworkers()
	bc_df = preprocess_bc_file(bc_file, bc_rev)
	if workers == 1
		classify_sequences(file_R1, bc_df, output_dir, max_error_rate, min_delta, mismatch, indel)
	else
		divide_fastq(file_R1, output_dir, workers)
		pmap(x -> mlt_demltplex(x, bc_df, output_dir, max_error_rate, min_delta, mismatch, indel), 1:workers)
		paths = []
		for (root, dirs, files) in walkdir(output_dir)
			for file in files
				push!(paths, joinpath(root, file))
			end
		end
		merge_fastq_files(paths, bc_df, output_dir)

		rm(joinpath(output_dir, "divided_fastq"), recursive = true)
		for i in 1:workers
			rm(joinpath(output_dir, "thread" * string(i)), recursive = true)
		end
	end
end