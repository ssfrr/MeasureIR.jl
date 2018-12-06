using MeasureIR
using Plots, PortAudio
using Unitful: s, Hz, kHz
using DSP
using FileIO
using LibSndFile: load, save

"""
    chunkminmax(data, nchunks)

Break `data` into `nchunks` chunks, and return a 3-element `Tuple` of
where the 1st item is a `Range` with the indices of each chunk, and the
second and third are the lower and upper bounds of each chunk, respectively.

## Example

```julia
sig = randn(1_000_000) .* cos.(linspace(0,2π, 1_000_000))

i, l, u = chunkminmax(sig, 1000)
plot(i, l, fillrange=u)
```
"""
function chunkminmax(data, nchunks)
    N = size(data, 1)
    C = size(data, 2)
    nchunks < N || return 1:N, data, data
    lower = zeros(nchunks, C)
    upper = zeros(nchunks, C)

    chunksize = floor(Int, N/nchunks)
    for i in 1:nchunks
        offset = (i-1) * chunksize
        n = min(chunksize, N-offset)
        lower[i, :] = minimum(view(data, (1:n) + offset, :), 1)
        upper[i, :] = maximum(view(data, (1:n) + offset, :), 1)
    end

    (0:(nchunks-1))*chunksize+1, lower, upper
end

function sigplot(sig, sr=1; kwargs...)
    i,l,u = chunkminmax(sig, 1000)
    plot(i/sr, l, fillrange=u, alpha=0.8; kwargs...)
end

function sigplot!(sig, sr=1; kwargs...)
    i,l,u = chunkminmax(sig, 1000)
    plot!(i/sr, l, fillrange=u, alpha=0.8; kwargs...)
end

function measure(measurement)
    # wait so we don't get the sound of the keypress
    sleep(1)
    str = PortAudioStream(2, 2; synced=true)
    sig = stimulus(measurement) * 0.5
    @async write(str, sig)
    resp = read(str, length(sig))
    close(str)
    resp
end

resps = []
stims = []
irs = []
for f in [expsweep, golay, impulse, mls, schroeder]
    meas = f(5*48000)
    resp = measure(meas)
    push!(resps, resp)
    push!(irs, analyze(meas, resp))
end
for (f, resp) in zip([expsweep, golay, impulse, mls, schroeder], resps)
    meas = f(5*48000)
    push!(irs, analyze(meas, resp))
end

sigplot(irs[5])
sigplot(stimulus(schroeder(5*48000)))
sigplot(resp)
sigplot(analyze(schroeder(5*48000), resp[:, 2]))
sigplot(xcorr(Float64.(resp[:, 2].data), stimulus(schroeder(5*48000))))

for (resp, stim) in zip(resps, stims)
    LibSndFile.save("/home/sfr/.julia/v0.6/MeasureIR/$stim.wav", resp)
end

resp = LibSndFile.load("/home/sfr/.julia/v0.6/MeasureIR/schroeder.wav")


stims = ["expsweep", "golay", "impulse", "mls", "schroeder"]

testir = randn(1000) .* exp(linspace(0,-10,1000))
m = schroeder(5*48000)
sigplot(analyze(m, conv(stimulus(m), testir))[1:1000])
sqrt(sum(stimulus(m).^2))
maximum(testir) / maximum(analyze(m, conv(stimulus(m), testir))[1:1000]/(5*48000))
plot!(testir)

sigplot(stimulus(m))
sigplot(analyze(meas, resp))
