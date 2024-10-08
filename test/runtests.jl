using Demux
using Test

const demo_FASTQ_file1 = "demo_file1.fastq"
const demo_FASTQ_file2 = "demo_file2.fastq"
const demo_FASTQ_file1_mutated = "demo_file1_mutated.fastq"
const demo_FASTQ_file2_mutated = "demo_file2_mutated.fastq"
const demo_bc_file = "demo_bc_file.csv"

function check_output_files(output_dir::String, ideal_dir::String)
    ideal_files = readdir(ideal_dir)

    for ideal_file in ideal_files
        output_file_path = joinpath(output_dir, ideal_file)
        ideal_file_path = joinpath(ideal_dir, ideal_file)

        @test isfile(output_file_path)
        @test isfile(ideal_file_path)   

        diff_result = readchomp(`diff $output_file_path $ideal_file_path`)
        @test diff_result == ""
    end
end

@testset "Demux.jl" begin
    @testset "Single thread with 1 FASTQ file" begin
        mktempdir() do output_dir
            ideal_dir = "path/to/ideal_outputs"

            execute_demultiplexing(demo_FASTQ_file1, demo_bc_file, output_dir)
            check_output_files_with_diff(output_dir, ideal_dir)
        end
    end

    @testset "Single thread with 2 FASTQ files" begin
        mktempdir() do output_dir
            ideal_dir = "path/to/ideal_outputs"

            execute_demultiplexing(demo_FASTQ_file1, demo_FASTQ_file2, demo_bc_file, output_dir)
            check_output_files_with_diff(output_dir, ideal_dir)
        end
    end

    @testset "Multi-thread with 1 FASTQ file" begin
        mktempdir() do output_dir
            ideal_dir = "path/to/ideal_outputs"

            addprocs(4)
            @everywhere using Demux
            execute_demultiplexing(demo_FASTQ_file1, demo_bc_file, output_dir)
            check_output_files_with_diff(output_dir, ideal_dir)
        end
        rmprocs(nworkers())
    end

    @testset "Multi-thread with 2 FASTQ files" begin
        mktempdir() do output_dir
            ideal_dir = "path/to/ideal_outputs"

            addprocs(4)
            @everywhere using Demux
            execute_demultiplexing(demo_FASTQ_file1, demo_FASTQ_file2, demo_bc_file, output_dir)
            check_output_files_with_diff(output_dir, ideal_dir)
        end
        rmprocs(nworkers())
    end

    @testset "Multi-thread with 2 FASTQ files with options" begin
        mktempdir() do output_dir
            ideal_dir = "path/to/ideal_outputs"

            addprocs(16)
            @everywhere using Demux
            execute_demultiplexing(demo_FASTQ_file1_mutated, demo_FASTQ_file2_mutated, demo_bc_file, output_dir,
                max_error_rate=0.1, min_delta=0.05, mismatch=2, indel=2, classify_both=true, bc_complement=true,
                output_prefix1="test_prefix1", output_prefix2="test_prefix2")
            
            check_output_files_with_diff(output_dir, ideal_dir)
        end
        rmprocs(nworkers())
    end
end