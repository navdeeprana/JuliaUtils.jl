using CSV
using DataFrames
using HDF5

read_data(fname, dsets; slice = (:, :)) = [h5read(fname, dset, slice) for dset in dsets]

function write_data(fname; mode = "cw", overwrite = true, kwargs...)
    h5open(fname, mode) do fid
        for (k, v) in pairs(kwargs)
            ds = string(k)
            if haskey(fid, ds) & overwrite
                delete_object(fid, ds)
            end
            fid[ds] = v
        end
    end
end

function read_xydata(fname; delim = " ", header = false, fcol = :Column1, ffun = (x -> true), kwargs...)
    df = CSV.read(
        fname,
        DataFrame;
        delim = delim,
        ignorerepeated = true,
        header = header,
        ntasks = 1,
        kwargs...
    )
    return eachcol(filter(fcol => ffun, df))
end

function write_xydata(fname, data; kwargs...)
    return CSV.write(
        fname,
        DataFrame(data, :auto),
        delim = " ",
        writeheader = false,
        kwargs...
    )
end
