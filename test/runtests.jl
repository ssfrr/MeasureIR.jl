using MeasureIR
using TestSetExtensions
using Unitful: s, Hz
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

function padequal(x::AbstractArray{T}, y::AbstractArray{T}) where {T}
    if length(x) > length(y)
        y = [y; zeros(T, length(x)-length(y))]
    elseif length(y) > length(x)
        x = [x; zeros(T, length(y)-length(x))]
    end
    x, y
end
@testset "MeasureIR" begin
    @testset "Golay Sequences" ExtendedTestSet begin
        for L in (8, 16)
            seq = golay(L)
            @test length(seq.A) == L
            @test length(seq.B) == L
            @test stimulus(seq) == [seq.A; zeros(L); seq.B; zeros(L)]
            # should get pure impulse
            @test isapprox(analyze(seq, stimulus(seq)),
                           [1, zeros(L-1)...],
                           atol=sqrt(eps()))
        end

        local seq
        # should upgrade non-powers-of-two
        seq = golay(13)
        @test length(seq.A) == 16
        @test length(seq.B) == 16

        # test equality
        @test golay(16) == golay(16)
        @test golay(16) != golay(8)

        kernel = rand(16)
        seq = golay(16)
        measured = conv(stimulus(seq), kernel)
        @test isapprox(analyze(seq, measured), kernel,
                       atol=sqrt(eps()))

        # just make sure we're not crazy slow with realistic measurement lengths
        @test @elapsed(golay(nextpow2(48000*6))) < 0.5
    end

    @testset "Unitful initialization" begin
        # unitful sampling rate is optional
        @test golay(2s, 48000Hz) == golay(96000)
        @test golay(2s, 48000) == golay(96000)
    end
end
