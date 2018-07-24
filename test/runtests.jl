using MeasureIR
using TestSetExtensions
using Unitful: s, kHz
using SampledSignals: SampleBuf, samplerate
using Suppressor
using DSP
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
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
        testir = [zeros(delay); randn(L) .* exp.(linspace(0,-8, L))]
        # filter the impulse response so it's within the passband of the sweep.
        # We should be able to pretty much perfectly recover it
        testir = filt(digitalfilter(Bandpass(0.05/π, 0.5),
                                    FIRWindow(hanning(128))),
                      testir)
        resp = conv(stim, testir)
        ir = analyze(meas, resp)
        @test sum((ir[1:delay+L]-testir).^2)/sum(testir.^2) < 0.01
    end

    @testset "$funcname - noisefloor" begin
        @test isapprox(noisefloor(meas, stim)[], 0.0, atol=0.01)
        nf = noisefloor(meas, stim .+ randn.() .* 0.2)[]
        @test isapprox(nf,  0.2^2, rtol=0.05)
    end

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

    @testset "$funcname - analysis works with SampleBuf" begin
        meas = measfunc(16)
        stim = stimulus(meas)
        stimbuf = SampleBuf(stim, 44100)
        irbuf = analyze(meas, stimbuf)
        @test irbuf isa SampleBuf
        @test samplerate(irbuf) == 44100
    end
end

@testset "MeasureIR" ExtendedTestSet begin
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
    @test indmax(analyze(meas, stim)) == 1
    @test indmax(analyze(meas, stim, noncausal=true)) == L+1
end

# only the golay measurement truncates and warns - for the others theres
# no reason to truncate - though this may be reworked later
@testset "golay - truncation warning" begin
    L = 4096
    meas = golay(L)
    ir = randn(L÷2) .* exp.(linspace(0,-8, L÷2))
    stim = stimulus(meas)
    out = @color_output false @capture_err analyze(meas, conv(stim, ir))
    @test out == ""
    ir = randn(3L÷2) .* exp.(linspace(0,-8, 3L÷2))
    out = @color_output false @capture_err analyze(meas, conv(stim, ir))
    @test out == string("Warning: nonsilent samples past end of analysis ",
                        "window. Check your test signal is long enough for ",
                        "the response\n")
end

# not working right now. See source for `analyze(::ExpSweep, ...) for details
# @testset "Nonlinearity detection" begin
#     L = 4096
#     meas = expsweep(L)
#     stim = stimulus(meas)
#     out = @color_output false @capture_err analyze(meas, stim)
#     @test out == ""
#     out = @color_output false @capture_err analyze(meas, tanh.(stim*2))
#     @test out == "Warning: Energy in noncausal IR is above noise floor. Check for nonlinearity\n"
#
#     # now try with some noise and a more realistic IR
#     testir = randn(512) .* exp.(linspace(0,-8, 512)) / 2
#     resp = conv(stim, testir)[1:length(stim)]
#     out = @color_output false @capture_err analyze(meas, stim.+randn.().*0.1)
#     @test out == ""
#     out = @color_output false @capture_err analyze(meas, tanh.(stim*2).+randn.().*0.1)
#     @test out == "Warning: Energy in noncausal IR is above noise floor. Check for nonlinearity\n"
# end

@testset "Golay Sequences" begin
    for L in (8, 16)
        seq = golay(L)
        @test length(seq.A) == L
        @test length(seq.B) == L
        @test stimulus(seq) == [zeros(L);
                                seq.amp .* seq.A;
                                zeros(L);
                                seq.amp .* seq.B;
                                zeros(L)]
    end

    # should upgrade non-powers-of-two
    @test golay(13) == golay(16)

    # test equality
    @test golay(16) == golay(16)
    @test golay(16) != golay(8)

    # default amplitude is 1/2.2
    @test maximum(stimulus(golay(16))) ≈ 1/2.2
    @test maximum(stimulus(golay(16, amp=1.0))) ≈ 1.0

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

end # @testset MeasureIR
