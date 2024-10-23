using BammannChiesaJohnsons
using CSV
using DataFrames
using Test

@testset verbose=true "BammannChiesaJohnsons.jl" begin
    df_file     = CSV.read("Data_Tension_e0002_T295.csv", DataFrame;
    # df_file     = CSV.read("Data_Tension_e570_T295.csv", DataFrame;
        header=true, delim=',', types=[Float64, Float64, Float64, Float64, String])
    params      = begin
        df          = CSV.read("Props_BCJ_4340_fit.csv", DataFrame; header=true, delim=',', types=[String, Float64])
        rowsofconstants = findall(occursin.(r"C\d{2}", df[!, "Comment"]))
        C_0         = Vector{Float64}(undef, 20)
        C_0[rowsofconstants] .= df[!, "For Calibration with vumat"][rowsofconstants]
        bulk_mod    = df[!, "For Calibration with vumat"][findfirst(occursin("Bulk Mod"), df[!, "Comment"])]
        shear_mod   = df[!, "For Calibration with vumat"][findfirst(occursin("Shear Mod"), df[!, "Comment"])]
        Dict( # collect as dictionary
            "C01"       => C_0[ 1],
            "C02"       => C_0[ 2],
            "C03"       => C_0[ 3],
            "C04"       => C_0[ 4],
            "C05"       => C_0[ 5],
            "C06"       => C_0[ 6],
            "C07"       => C_0[ 7],
            "C08"       => C_0[ 8],
            "C09"       => C_0[ 9],
            "C10"       => C_0[10],
            "C11"       => C_0[11],
            "C12"       => C_0[12],
            "C13"       => C_0[13],
            "C14"       => C_0[14],
            "C15"       => C_0[15],
            "C16"       => C_0[16],
            "C17"       => C_0[17],
            "C18"       => C_0[18],
            "C19"       => C_0[19],
            "C20"       => C_0[20],
            "bulk_mod"  => bulk_mod,
            "shear_mod" => shear_mod)
    end
    bcj_ref     = BCJ_metal(295., 2e-3, 1., 200, 1, params)
    # bcj_ref     = BCJ_metal(295., 570., 0.15, 200, 1, params)
    bcj_current = BCJ_metal_currentconfiguration_init(bcj_ref, DK)
    solve!(bcj_current)
    σ__ = bcj_current.σ__ # [1, :]
    # σvM     = σ__
    σvM     = sum(map.(x->x^2., [σ__[1, :] - σ__[2, :], σ__[2, :] - σ__[3, :], σ__[3, :] - σ__[1, :]])) + (
        6sum(map.(x->x^2., [σ__[4, :], σ__[5, :], σ__[6, :]])))
    σvM     = sqrt.(σvM .* 0.5)
    idx = []
    for t in df_file[!, "Strain"]
        j = findlast(bcj_current.ϵ__[1, :] .<= t)
        push!(idx, !isnothing(j) ? j : findfirst(bcj_current.ϵ__[1, :] .>= t))
    end
    err = sum(((df_file[!, "Stress"] .* 1e6) - σvM[idx]) .^ 2.)
    x = sum((df_file[!, "Stress"] .* 1e6) .^ 2)
    # println(x)
    err /= if x > 1e9^2
        1e9^2
    elseif x > 1e6^2
        1e6^2
    elseif x > 1e3^2
        1e3^2
    else
        1.
    end
    # println(err)
    # p = scatter(df_file[!, "Strain"], df_file[!, "Stress"] .* 1e6, label="Data", ylims=(0., 2e9))
    # plot!(p, bcj_current.ϵ__[1, :], σvM, label="Model")
    # display(p)
    @test isapprox(err, 0.020599315626415465; atol=1e-6)
end
