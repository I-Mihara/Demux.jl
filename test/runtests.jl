using Demux
using Test

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
            ideal_dir = "results/demo1_R1"

            files = readdir("FASTQ_files/demo1_R1")
            for file in files
                execute_demultiplexing("FASTQ_files/demo1_R1/$file", "reference_files/demo1.tsv", "output_dir")
            end
            check_output_files(output_dir, ideal_dir)
        end
    end

    @testset "Single thread with 2 FASTQ files" begin
        mktempdir() do output_dir
            ideal_dir = "results/demo1_R2"

            files_R1 = readdir("FASTQ_files/demo1_R1")
            files_R2 = readdir("FASTQ_files/demo1_R2")
            for (file_R1, file_R2) in zip(files_R1, files_R2)
                execute_demultiplexing("FASTQ_files/demo1_R1/$file", "FASTQ_files/demo1_R2/$file", "reference_files/demo1.tsv", "output_dir")
            end
            check_output_files(output_dir, ideal_dir)
        end
    end

    @testset "Multi-thread with 1 FASTQ file" begin
        mktempdir() do output_dir
            ideal_dir = "results/demo1_R1"

            addprocs(4)
            @everywhere using Demux
            files = readdir("FASTQ_files/demo1_R1")
            for file in files
                execute_demultiplexing("FASTQ_files/demo1_R1/$file", "reference_files/demo1.tsv", "output_dir")
            end
            check_output_files(output_dir, ideal_dir)
        end
        rmprocs(nworkers())
    end

    @testset "Multi-thread with 2 FASTQ files" begin
        mktempdir() do output_dir
            ideal_dir = "results/demo1_R2"

            addprocs(4)
            @everywhere using Demux
            files_R1 = readdir("FASTQ_files/demo1_R1")
            files_R2 = readdir("FASTQ_files/demo1_R2")
            for (file_R1, file_R2) in zip(files_R1, files_R2)
                execute_demultiplexing("FASTQ_files/demo1_R1/$file", "FASTQ_files/demo1_R2/$file", "reference_files/demo1.tsv", "output_dir")
            end
            check_output_files(output_dir, ideal_dir)
        end
        rmprocs(nworkers())
    end

    @testset "Multi-thread with 2 FASTQ files with options" begin
        mktempdir() do output_dir
            ideal_dir = "results/demo2"

            addprocs(4)
            @everywhere using Demux
            files_R1 = readdir("FASTQ_files/demo2_R1")
            files_R2 = readdir("FASTQ_files/demo2_R2")
            for (file_R1, file_R2) in zip(files_R1, files_R2)
                execute_demultiplexing("FASTQ_files/demo2_R1/$file_R1", "FASTQ_files/demo2_R2/$file_R2", "reference_files/demo2.csv", "output_dir",
                    max_error_rate=0.25, min_delta=0.15, mismatch=1, indel=2, classify_both=true, bc_complement=true,
                    output_prefix="test_prefix1", output_prefix2="test_prefix2")
            end
            check_output_files(output_dir, ideal_dir)
        end
        rmprocs(nworkers())
    end
end