using CSV
using DataFrames
using HDF5
using DelimitedFiles

"""
!!! note "Efficiency"
   It is noted that CSV is faster compared to Readdlm package for large files.
   The loadtxt function defined below uses readdlm.
   Check the timings using @time macro for your needs.
   see ["stackoverflow"](https://stackoverflow.com/questions/26810171/is-there-an-equivalent-or-close-to-numpy-loadtxt-for-julia)

!!! warning "Not tested"
   Read-write using CSV and HDF5 are not tested fully. Use [`loadtxt`](@ref) function for the time being.
"""

"""
Read from HDF5 file.
$(TYPEDSIGNATURES)
Args:\\
[Req] fname; (filename)
[Opt] slice = (:, :) (Read HDF5 documentation)\\
[Opt] dsets = ["phi1", "phi2"] (Read HDF5 documentation)
"""
function read_data(fname; slice = (:, :), dsets = ["phi1", "phi2"])
   return [h5read(fname, dset, slice) for dset in dsets]
end
"""
Write to HDF5 file.
$(TYPEDSIGNATURES)
Args:\\
[Req] fname; (filename)\\
[Opt] mode = "cw", (Read HDF5 documentation)\\
[Opt] overwrite = true, (Read HDF5 documentation)\\
[Opt] kwargs...   (Read HDF5 documentation)\\
"""
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
"""
Read and XY format data file usig CSV package.
$(TYPEDSIGNATURES)
Args:\\
[Req] fname; (filename)\\
[Opt] delim = " ",   (read documentation for CSV.read)
[Opt] header = false,   (read documentation for CSV.read)
[Opt] fcol = :Column1,  (read documentation for CSV.read)
[Opt] ffun = (x -> true),  (read documentation for CSV.read)
[Opt] kwargs...   (read documentation for CSV.read)
"""
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
"""
Wrtie XY format data to a file using CSV package.
$(TYPEDSIGNATURES)
Args:\\
[Req] fname (filename)\\
[Opt] data (data to be written)
[Opt] kwargs...   (read documentation for CSV.read)
"""
function write_xydata(fname, data; kwargs...)
   return CSV.write(
      fname,
      DataFrame(data, :auto),
      delim = " ",
      writeheader = false,
      kwargs...
   )
end
"""
Read XY format data file using Readdlm package.
$(TYPEDSIGNATURES)
Args:\\
[Req] fname (filename)\\
[Opt] skiprow (# of header to be skiped)\\
[Opt] delim::AbstractChar = '\t',   (seperator)\\
[Opt] eol::AbstractChar='\n', (end of line character)\\
[Opt] kwargs...   (read documentation forDelimitedFiles.readdlm)
"""
function loadtxt(fname; skiprows=0, delim::AbstractChar = '\t',eol::AbstractChar='\n', kwargs...)
   return readdlm(
      fname,
      delim,
      eol,
      skipstart=skiprows,
      skipblanks = true,
      kwargs...)
end
"""
Write XY format data file using Readdlm package.
$(TYPEDSIGNATURES)
Args:\\
[Req] fname (filename)\\
[Req] data (data to be written, Vector/Matrix or an iterable collection of iterable rows)\\
[Opt] mode = "w" (mode of opening file, "w" = opwn for write)
[Opt] delim::AbstractChar = '\t',   (seperator)\\
[Opt] kwargs...   (read documentation for DelimitedFiles.writedlm)\\
"""
function savetxt(fname, data; mode="w", delim::AbstractChar = '\t', kwargs...)
   open(fname, mode) do fid
      writedlm(fid, data,delim,kwargs...)
  end
end
