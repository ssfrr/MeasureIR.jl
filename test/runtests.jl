using MeasureIR
using TestSetExtensions
using Unitful: s, kHz
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

@testset "MeasureIR" ExtendedTestSet begin
    @testset "Golay Sequences" begin
        for L in (8, 16)
            seq = golay(L)
            @test length(seq.A) == L
            @test length(seq.B) == L
            @test stimulus(seq) == [seq.amp .* seq.A; zeros(L); seq.amp .* seq.B; zeros(L)]
            # should get pure impulse
            @test isapprox(analyze(seq, stimulus(seq)),
                           [1, zeros(L-1)...],
                           atol=sqrt(eps()))
        end

        local seq
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

        # test reconstruction with random impulse response
        kernel = rand(16)
        seq = golay(16)
        measured = conv(stimulus(seq), kernel)
        @test isapprox(analyze(seq, measured), kernel,
                       atol=sqrt(eps()))

        # just make sure we're not crazy slow with realistic measurement lengths
        @test @elapsed(golay(nextpow2(48000*6))) < 0.5
    end

    @testset "Unitful Golay initialization" begin
        # unitful sampling rate is optional
        @test golay(2s, 1kHz) == golay(2000)
        @test golay(2s, 1000) == golay(2000)
    end

    @testset "Impulse" begin
        m = impulse(4)
        @test length(stimulus(m)) == 4
        @test stimulus(m) == [1.0, 0.0, 0.0, 0.0]
        @test analyze(m, stimulus(m)) == [1.0, 0.0, 0.0, 0.0]
        @test analyze(m, [2.0, 1.2, 6.5, 1.4]) == [2.0, 1.2, 6.5, 1.4]
    end

    @testset "Unitful Impulse" begin
        @test impulse(0.5s, 48000) == impulse(24000)
        @test impulse(0.5s, 48kHz) == impulse(24000)
    end

    @testset "Exponential Sweep" begin
        m = expsweep(1024)
        stim = stimulus(m)
        @test length(stim) == 1024
        @test analyze(m, stim) == [1.0; ones(1023)]
    end
end
