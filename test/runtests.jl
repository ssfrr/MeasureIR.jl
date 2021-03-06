using MeasureIR
using Unitful: s, kHz
using SampledSignals: SampleBuf, samplerate
using Suppressor
using DSP
using Test

"""
Simulate running the given stimulus through a real system, including a linear
IR, nonlinear distortion, noise, and clock skew.

If `dist` is given, it should be a multiplier that scales the stimulus before
passing into `tanh`. The distorted signal is scaled to have the same energy as
the original. A multiplier of 1 gives pretty mild distortion, and it goes up
from there.

Returns the simulated response, as well as the IR used.
"""
function simIR(stim; ir=nothing, dist=nothing, noisepow=-Inf, skew=1.0)
    if dist !== nothing
        stimenergy = sum(x->x^2, stim)
        stim .= tanh.(stim)
        distenergy = sum(x->x^2, stim)

        stim .*= sqrt(sigenergy/distenergy)
    end
    if ir === nothing
        ir = [zeros(100); 1; randn.() .* exp.(range(-1.5, -8, length=7500))]
    end
    resp = conv(stim, ir)
    if noisepow > -Inf
        resp .+= db2amp(noisepow) .* randn.()
    end
    if skew != 1.0
        resp = resample(resp, skew)
    end

    resp, ir
end

"""
Generates a broadband random signal of length `N` that's periodic with period
`P`. `P` does not need to be an integer.
"""
function randperiodic(N, P)
    f = 1/P # fundamental freq
    # we add eps() here so we don't get hit by floating-point error when f
    # divides evenly into 1
    H = Int(fld(0.5+eps(), f)) # number of harmonics we can fit below nyquist
    sum(cos.(2π.*(h.*f.*(0:N-1).+rand()))./H for h in 1:H)
end

function testmeasure(measfunc)
    meas = measfunc(8192)
    stim = stimulus(meas)
    funcname = split(string(measfunc), ".")[2]
    @testset "$funcname - pure impulse IR" begin
        ir = analyze(meas, stim)
        maxval, maxidx = findmax(ir)
        @test maxidx == 1
        @test maximum(ir[10:end].^2) < 0.01*maxval^2
    end
    @testset "$funcname - pure delay IR" begin
        delay = 150
        ir = analyze(meas, [zeros(delay); stim])
        maxval, maxidx = findmax(ir)
        @test maxidx == delay+1
        @test maximum(ir[1:delay-10].^2) < 0.01*maxval^2
        @test maximum(ir[delay+10:end].^2) < 0.01*maxval^2
    end
    @testset "$funcname - convolved IR" begin
        delay = 250
        L = 250
        testir = [zeros(delay); randn(L) .* exp.(range(0,-8, length=L))]
        # filter the impulse response so it's within the passband of the sweep.
        # We should be able to pretty much perfectly recover it
        testir = filt(digitalfilter(Bandpass(0.05/π, 0.5),
                                    FIRWindow(hanning(128))),
                      testir)
        resp = conv(stim, testir)
        ir = analyze(meas, resp)
        @test sum((ir[1:delay+L]-testir).^2)/sum(testir.^2) < 0.01
    end

    # FIXME: is this still supported at all?
    # @testset "$funcname - noisefloor" begin
    #     @test isapprox(noisefloor(meas, stim)[], 0.0, atol=0.01)
    #     nf = noisefloor(meas, stim .+ randn.() .* 0.2)[]
    #     @test isapprox(nf,  0.2^2, rtol=0.05)
    # end

    @testset "$funcname - multichannel response" begin
        delays = [150,200]
        resp1 = [zeros(delays[1]); stim[1:end-delays[1]]]
        resp2 = [zeros(delays[2]); stim[1:end-delays[2]]]
        resp = [resp1 resp2]

        ir = analyze(meas, resp)
        for ch in 1:2
            maxval, maxidx = findmax(ir[:,ch])
            @test maxidx == delays[ch]+1
            @test maximum(ir[1:delays[ch]-10, ch].^2) < 0.01*maxval^2
            @test maximum(ir[delays[ch]+10:end, ch].^2) < 0.01*maxval^2
        end
    end

    @testset "$funcname - vector/matrix type handling" begin
        meas = measfunc(16)
        stim = stimulus(meas)
        @test stim isa Vector
        @test analyze(meas, stim) isa Vector
        @test analyze(meas, [stim stim]) isa Matrix
    end

    # FIXME: this seems to be having a problem in `similar`
    # @testset "$funcname - analysis works with SampleBuf" begin
    #     meas = measfunc(16)
    #     stim = stimulus(meas)
    #     stimbuf = SampleBuf(stim, 44100)
    #     irbuf = analyze(meas, stimbuf)
    #     @test irbuf isa SampleBuf
    #     @test samplerate(irbuf) == 44100
    # end
end

@testset "MeasureIR" begin
# don't indent contents here, so we can run tests individually in Juno more
# easily

@testset "Standard Measurements" begin
    for f in [expsweep, golay, mls, rpms, impulse]
        testmeasure(f)
    end
end

@testset "Noncausal exponential sweep" begin
    L = 4096
    meas = expsweep(L)
    stim = stimulus(meas)
    @test indmax(analyze(meas, stim, warn=false)) == 1
    @test indmax(analyze(meas, stim, noncausal=true)) == L+1
end

# only the golay measurement truncates and warns - for the others theres
# no reason to truncate - though this may be reworked later
@testset "golay - truncation warning" begin
    L = 4096
    meas = golay(L)
    ir = randn(L÷2) .* exp.(linspace(0,-8, L÷2))
    stim = stimulus(meas)
    out = @capture_err analyze(meas, conv(stim, ir))
    @test out == ""
    ir = randn(3L÷2) .* exp.(linspace(0,-8, 3L÷2))
    out = @capture_err analyze(meas, conv(stim, ir))
    @test out == string("Warning: nonsilent samples past end of analysis ",
                        "window. Check your test signal is long enough for ",
                        "the response\n")
end

# not working right now. See source for `analyze(::ExpSweep, ...) for details
# @testset "Nonlinearity detection" begin
#     L = 4096
#     meas = expsweep(L)
#     stim = stimulus(meas)
#     out = @capture_err analyze(meas, stim)
#     @test out == ""
#     out = @capture_err analyze(meas, tanh.(stim*2))
#     @test out == "Warning: Energy in noncausal IR is above noise floor. Check for nonlinearity\n"
#
#     # now try with some noise and a more realistic IR
#     testir = randn(512) .* exp.(linspace(0,-8, 512)) / 2
#     resp = conv(stim, testir)[1:length(stim)]
#     out = @capture_err analyze(meas, stim.+randn.().*0.1)
#     @test out == ""
#     out = @capture_err analyze(meas, tanh.(stim*2).+randn.().*0.1)
#     @test out == "Warning: Energy in noncausal IR is above noise floor. Check for nonlinearity\n"
# end

@testset "Golay Sequences" begin
    for L in (8, 16)
        seq = golay(L)
        @test length(seq.A) == L
        @test length(seq.B) == L
        @test stimulus(seq) == [zeros(L);
                                seq.gain .* seq.A;
                                zeros(L);
                                seq.gain .* seq.B;
                                zeros(L)]
    end

    # should upgrade non-powers-of-two
    @test golay(13) == golay(16)

    # test equality
    @test golay(16) == golay(16)
    @test golay(16) != golay(8)

    @test maximum(stimulus(golay(16, gain=1.0))) ≈ 1.0

    # can set the bandwidth
    bwlim = golay(1024, upsample=2)
    spec = rfft(stimulus(bwlim))
    @test maximum(abs.(spec[round(Int, 0.6*length(spec)):end])) < 0.02
    resp = (analyze(bwlim, stimulus(bwlim)))
    # response energy should be compact into a bandlimited impulse
    @test sum(resp[20:60].^2) > 1.95
    @test maximum(resp[80:end]) < 1e16
end

@testset "Unitful Golay initialization" begin
    # unitful sampling rate is optional
    @test golay(2s, 1kHz) == golay(2000)
    @test golay(2s, 1000) == golay(2000)
end

@testset "Unitful Impulse" begin
    @test impulse(0.5s, 48_000) == impulse(24_000)
    @test impulse(0.5s, 48kHz) == impulse(24_000)
end


@testset "snr" begin
    sig = randn(100000)
    noise = randn(100000)

    @test isapprox(snr(sig*2+noise, sig), 4; rtol=0.02)
end

# sr = 48kHz
# f = 1kHz
# dur = 1s
# sig = [zeros(50000); pilot(dur, f, samplerate=sr); zeros(48000*10)]
# sig .+= sqrt(2) .* randn.()
# plt = plot(welch_pgram(sig[(1:10000) .+ 50000], 256).power)
@testset "findpilot" begin
    sr = 48kHz
    dur = 1s
    @testset "f=$f" for f in 0.99kHz:0.005kHz:1.01kHz
        sig = [zeros(50000); pilot(dur, f, samplerate=sr); zeros(48000*10)]
        sig .+= 8 .* randn.() # SNR of -21.1dB - not bad!
        onset = findpilot(sig, 1kHz, 1s; samplerate=sr)
        @test isapprox(onset, 50001, rtol=0.1)
    end
end

@testset "pilotsync" begin
    dur = 5s
    sr = 48kHz
    f = 1kHz
    L = 5*48000
    @testset "skew=$skew" for skew in 0.99:0.005:1.01
        p = pilot(dur, f; samplerate=sr)
        sig = [zeros(50000);
               isapprox(skew, 1) ? p : resample(p, skew);
               zeros(48000*15)]
        sig .+= randn.() ./ sqrt(2) # SNR of 0dB
        onset = findpilot(sig, f, dur; samplerate=sr)
        @test isapprox(pilotsync(sig[(0:L-1) .+ onset], f; samplerate=sr),
                       1/skew, rtol=5e-7)
    end
end

@testset "getperiod" begin
    N = 65535
    meas = mls(N, 20)

    # test at a variety of skews in ppm
    skews = 1 .+ (-400:80:400) ./ 1e6

    @testset "skew=$skew" for skew in skews
        resp, ir = simIR(stimulus(meas), skew=skew, noisepow=-8.2+20) # -20dB SNR
        period = getperiod(resp, N)
        @test isapprox(period, skew*N, rtol=1e-7)
    end
end

end # @testset MeasureIR
