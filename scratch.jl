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
    for i in 1:nchunks
        offset = (i-1) * chunksize
        n = min(chunksize, N-offset)
        lower[i, :] = minimum(view(data, (1:n) + offset, :), 1)
        upper[i, :] = maximum(view(data, (1:n) + offset, :), 1)
    end

    (0:(nchunks-1))*chunksize+1, lower, upper
end

function sigplot(sig, sr=1)
    i,l,u = chunkminmax(sig, 1000)
    plot(i/sr, l, fillrange=u, alpha=0.8)
end

function sigplot!(sig, sr=1)
    i,l,u = chunkminmax(sig, 1000)
    plot!(i/sr, l, fillrange=u, alpha=0.8)
end

function measure(measurement)
    # wait so we don't get the sound of the keypress
    sleep(1)
    str = PortAudioStream(1, 2; synced=true)
    sig = stimulus(measurement)
    @async write(str, [sig zeros(length(sig))])
    resp = read(str, length(sig))
    close(str)
    resp
end

fs = [expsweep, golay, impulse, mls, rpms]
fnames = [split(string(f), ".")[2] for f in fs]
ms = [f(5*44100) for f in fs]
resps = [measure(m) for m in ms]
irs = [analyze(m, r) for (m, r) in zip(ms, resps)]
plot(irs[5][1:25000], ylims=(-0.3, 0.3))
plot!(ir[1:25000], ylims=(-0.3, 0.3))

ir = analyze(rpms(5*44100), resps[5])
plot()

swp = expsweep(5*44100; nonlinear=true)
swp = golay(1*44100)
swp = mls(44100)
resp = convert(Array{Float64}, measure(swp))
ir = analyze(swp, resp)
sigplot(resp, 44100)
sigplot(ir, 44100)
size(ir)
plot(linspace(0,22050, size(ir,1)÷2+1), abs.(rfft(ir[:, 2])), legend=false)

m1 = rpms(5*44100)
m2 = rpms(5*44100)
m1.sig == m2.sig

testir = (rand(100)*2-1) .* e.^(linspace(0,-10,100))
plot(stimulus(m))
plot!(testir)
stim = stimulus(m)
resp = measure(m)
sigplot(resp)
sigplot(analyze(m, resp)[1:50000])


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

# quick-and-dirty true-peak estimation
truepeak(x::AbstractVector) = maximum(abs.(resample(x, 8//1)))
truepeak(x::AbstractMatrix) = mapslices(truepeak, x, 1)

N = 2^25
Ns = 2.^(5:25)
2^24/48000
truepeak(expsweep(N).sig)
truepeak([golay(N).A; golay(N).B])
truepeak(mls(N).sig)
rpmspeaks = [truepeak(rpms(N).sig * sqrt(N/2)) for N in Ns]
golaypeaks = [truepeak([golay(N).A;golay(N).B]) for N in Ns]
mlspeaks = [truepeak(mls(N).sig) for N in Ns]
plot(rpmspeaks)
plot(golaypeaks)
plot(mlspeaks)
g = golay(2^18)
histogram(resample([g.A;g.B], 8//1))
histogram(resample(mls(2^18).sig, 8//1))
histogram(resample(rpms(2^18).sig*sqrt(2^18), 8//1))
var(rpms(2^19).sig*sqrt(2^19))
N = 524288
truepeak(rpms(N).sig * sqrt(N/2))
sigplot(rpms(N).sig * sqrt(N/2))

import LibSndFile
filt = sinc.((-800:800)/8)/8
filt .*= gaussian(length(filt), 0.15)
plot(20log10.(abs.(rfft(filt))))
x = sinc.(((-100:100)+0.5))
LibSndFile.save("c:\\Users\\sfr\\Desktop\\sinc.wav", Float32.(x))
sticks(x)
x .*= gaussian(length(x), 0.15)
sticks(x)
xz = MeasureIR.zerostuff(x, 8) * 8
x8 = conv(xz,filt)
maximum(x8)
plot(x)
x8 = resample(x, 32//1)
20log10(maximum(x8))
db2amp(-1)
plot(x)
plot!((0:length(x8)-1)/32+1, x8)
xlims!(48,53)
maximum(x8)
sinc(0)
f = resample_filter(2.0, 32, 1.0, 80)
plot(filt(FIRFilter(f, 2.0), x))
plot(f)

function percentilethresh(x, p)
    xs = sort(abs.(x))
    xs[ceil(Int, p/100*length(x))]
end
