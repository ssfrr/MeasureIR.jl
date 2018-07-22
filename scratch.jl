using MeasureIR
using Plots, PortAudio
using Unitful: s, Hz, kHz
using DSP

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
    @show nchunks
    for i in 1:nchunks
        offset = (i-1) * chunksize
        n = min(chunksize, N-offset)
        @show i
        @show offset
        lower[i, :] = minimum(view(data, (1:n) + offset, :), 1)
        upper[i, :] = maximum(view(data, (1:n) + offset, :), 1)
    end

    (0:(nchunks-1))*chunksize+1, lower, upper
end

function sigplot(sig, sr=1)
    i,l,u = chunkminmax(sig, 1000)
    plot(i/sr, l, fillrange=u, alpha=0.8)
end

function measure(measurement)
    # wait so we don't get the sound of the keypress
    sleep(1)
    str = PortAudioStream(2, 2; synced=true)
    sig = stimulus(measurement)
    @async write(str, sig)
    resp = read(str, length(sig))
    close(str)
    resp
end

swp = expsweep(5*44100; nonlinear=true)
swp = golay(1*44100)
swp = mls(44100)
resp = convert(Array{Float64}, measure(swp))
ir = analyze(swp, resp)
sigplot(resp, 44100)
sigplot(ir, 44100)
size(ir)
plot(linspace(0,22050, size(ir,1)÷2+1), abs.(rfft(ir[:, 2])), legend=false)



testir = [zeros(100); (rand(100)*2-1) .* e.^(linspace(0,-10,100))]
plot(testir)

g = mls(44100)
stim = stimulus(g)
ir = analyze(g, conv(stim, testir))
plot(ir[1:500]*2.2)
plot!(testir)

stim = rand(8)
# resp = stim
resp = [stim; zeros(10)]
# plot(stim)
# plot(resp)

length(timexcorr(resp,stim))
plot(xcorr(resp,stim))
plot(causalxcorr(resp,stim))
vline!([length(stim)])

causalxcorr(resp, stim) = xcorr(resp, stim)[max(length(resp), length(stim)):end]

function timexcorr(u,v)
    out = zeros(length(u)+length(v)-1)
    i = 1
    for τ in -length(v)+1:length(u)-1
        low = max(0, τ)
        high = min(length(u)-1, length(v)+τ-1)
        out[i] = sum((u[t+1]*conj(v[t-τ+1]) for t in low:high))
        i += 1
    end

    out
end
