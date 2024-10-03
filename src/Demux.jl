__precompile__()

module Demux

export
	semiglobal_alignment,
	find_best_matching_bc,
	determine_filemname,
	write_fastq_entry,
	classify_sequences,

	preprocess_bc_file,
	preprocess_fastq,
	divide_fastq,
	mlt_demltplex,
	merge_fastq_files,
	execute_demultiplexing


using DataFrames, CSV, Distributed
include("classification.jl")
include("fileio.jl")
include("core.jl")

end#module
