"""
Script to sample the validation and landcover rasters and save them again as geoparquet (.pq) files.
"""

using SpaceLiDAR
using ProgressMeter
using GeoDataFrames
using DataFrames
using GeoArrays
using Statistics
using ProgressMeter
using StaticArrays
using GeoParquet

areas = ("sw", "nz", "nl")

demdir = "/mnt/d3770dc8-0b9d-4a21-84c7-d0a3dc5946c9/dtm/dtm/"
idatadir = "validation/data/"
odatadir = "/mnt/ec66e171-5639-4c62-9d2c-08e81c462669/validation/"

function sample!(df, ga, symbol, buffer)
    return df[!, symbol] = sample_buffer.(Ref(ga), df.longitude, df.latitude, buffer)
end

function sample_diff!(df, ga, symbol, buffer, refsymbol, reducer=median)
    return df[!, symbol] =
        sample_buffer.(Ref(ga), df.longitude, df.latitude, buffer, reducer) .-
        df[!, refsymbol]
end

unzip(a) = map(x -> getfield.(a, x), fieldnames(eltype(a)))

function sample_diff_shift!(df, ga, symbol, dist, refsymbol, angle=0)
    x, y = unzip(SpaceLiDAR.shift.(df.longitude, df.latitude, df.angle .+ angle, dist))
    return df[!, symbol] =
        sample.(Ref(ga), x, y) .- df[!, refsymbol]
end

function smooth(A, Δ=30)
    B = copy(A)
    R = LinearIndices(A)
    I_first, I_last = first(R), last(R)
    @inbounds @simd for I in R
        patch = max(I_first, I - Δ):min(I_last, I + Δ)
        B[I] = median(view(A, patch))
    end
    return B
end


const strategy = GeoArrays.Center()

function sample_buffer(ga::GeoArray, x::Real, y::Real, buffer=0, reducer=median)::Float32
    I = SVector{2}(x, y)
    i, j = indices(ga, I, strategy)
    0 < i - buffer <= Base.size(ga.A)[1] || return NaN32
    0 < j - buffer <= Base.size(ga.A)[2] || return NaN32
    0 < i + buffer <= Base.size(ga.A)[1] || return NaN32
    0 < j + buffer <= Base.size(ga.A)[2] || return NaN32
    data = collect(skipmissing(ga.A[i-buffer:i+buffer, j-buffer:j+buffer, 1]))
    mask = isfinite.(data) .& (data .!= -9999) .& (data .!= 3.4028235f38)
    if sum(mask) == 0
        return Inf32
    else
        return Float32(reducer(data[mask]))
    end
end

function sample(ga::GeoArray, x::Real, y::Real)::Float32
    I = SVector{2}(x, y)
    i, j = indices(ga, I, strategy)
    0 < i <= Base.size(ga.A)[1] || return NaN32
    0 < j <= Base.size(ga.A)[2] || return NaN32
    0 < i <= Base.size(ga.A)[1] || return NaN32
    0 < j <= Base.size(ga.A)[2] || return NaN32
    data = ga[i, j, 1]

    if !ismissing(data) && isfinite(data) && (data != -9999) && (data != 3.4028235f38)
        return Float32(data)
    else
        return Inf32
    end
end

for area in areas
    demfn = joinpath(demdir, area, "dtm.tif")
    granules = readdir(joinpath(idatadir, area), join=true)
    filter!(x -> endswith(x, ".pq"), granules)
    @info length(granules), isfile(demfn)
    ga = GeoArrays.read(demfn, masked=false)

    coverfn = "/mnt/ec66e171-5639-4c62-9d2c-08e81c462669/esa-worldcover/$area/esa-worldcover.tif"
    coverga = GeoArrays.read(coverfn)

    @showprogress "Sampling $area" 1 for granule in granules
        fn = joinpath(odatadir, area, basename(granule))
        # isfile(fn) && continue

        if occursin("gedi", granule)
            buffer = round(Int, 12.5 / (ga.f.linear[1] * 110_000))
        else
            buffer = round(Int, 9.0 / (ga.f.linear[1] * 110_000))
        end
        df = GeoParquet.read(granule)
        df.track = String.(df.track)
        sample_diff!(df, ga, :ref_center, 0, :height)
        sample_diff!(df, ga, :ref_median, buffer, :height, median)
        sample_diff!(df, ga, :ref_mean, buffer, :height, mean)
        sample!(df, coverga, :landcover, 0)
        filter!(:ref_center => isfinite, df)
        filter!(:landcover => isfinite, df)
        df.landcover = Int16.(df.landcover)

        if area == "sw"
            df.angle = smooth(df.angle, 30)

            gas = GeoArrays.read(
                "/mnt/d3770dc8-0b9d-4a21-84c7-d0a3dc5946c9/dtm/dtm/sw/dtm_wgs84_slope.tif",
                masked=false,
            )
            sample_diff_shift!(df, ga, :ref_f25, 2.5, :height)
            sample_diff_shift!(df, ga, :ref_f5, 5, :height)
            sample_diff_shift!(df, ga, :ref_f75, 7.5, :height)
            sample_diff_shift!(df, ga, :ref_f10, 10, :height)
            sample_diff_shift!(df, ga, :ref_b25, -2.5, :height)
            sample_diff_shift!(df, ga, :ref_b5, -5, :height)
            sample_diff_shift!(df, ga, :ref_b75, -7.5, :height)
            sample_diff_shift!(df, ga, :ref_b10, -10, :height)
            sample_diff_shift!(df, ga, :ref_l25, 2.5, :height, -90)
            sample_diff_shift!(df, ga, :ref_l5, 5, :height, -90)
            sample_diff_shift!(df, ga, :ref_l75, 7.5, :height, -90)
            sample_diff_shift!(df, ga, :ref_l10, 10, :height, -90)
            sample_diff_shift!(df, ga, :ref_r25, 2.5, :height, 90)
            sample_diff_shift!(df, ga, :ref_r5, 5, :height, 90)
            sample_diff_shift!(df, ga, :ref_r75, 7.5, :height, 90)
            sample_diff_shift!(df, ga, :ref_r10, 10, :height, 90)
            sample!(df, gas, :slope, 0)
        end
        size(df)[1] == 0 && continue
        GeoParquet.write(fn, df)
    end
end
