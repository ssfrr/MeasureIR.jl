struct MLS{AT} <: IRMeasurement
    sig::AT
end

function mls(L)
    N = Int(log2(nextpow2(L)))
    1 <= N <= length(mlspolys) || throw(ArgumentError("N must be positive and <= $(length(mlspolys))"))
    poly = mlspolys[N]
    # where we'll keep the state used to compute the next value. The
    # rightmost bit is the most recent value. Nothing super special about this
    # initialization value, as long as there's at least 1 nonzero bit in the
    # filter span
    register = 0x0000000000000001
    # bit mask for the filter taps
    taps = reduce(|, (0x0000000000000001 << (N-1-shift) for shift in poly))
    out = zeros(2^N-1)
    # initialize our output array to match the register above - with only the
    # most recent sample as 1.
    out[1:N-1] = -1
    out[N] = 1
    for i in N+1:length(out)
        nextval = UInt64(count_ones(register & taps)) & 0x0000000000000001
        register = (register << 1) | nextval
        out[i] = float(nextval) * 2 - 1
    end

    MLS(out)
end

stimulus(m::MLS) = m.sig

function analyze(m::MLS, response::AbstractArray)
    L = length(m.sig)
    mapslices(response, 1) do v
        xcorr(v, m.sig)[L:end] ./ L
    end
end

# these polynomial coefficients represent the tap indices of the filter, where
# zero is the oldest value (delay of N)
# they are from:
# Stahnke, Wayne 1973. Primitive Binary Polynomials.
# Mathematics of Computation 27(124): 977–980.
# the paper has higher orders, but this seems like more than enough for
# reasonable audio purposes. Feel free to submit a PR with more if they are
# useful to you.

const mlspolys = [
    [0],            # 01
    [1,  0],        # 02
    [1,  0],        # 03
    [1,  0],        # 04
    [2,  0],        # 05
    [1,  0],        # 06
    [1,  0],        # 07
    [6,  5,  1, 0], # 08
    [4,  0],        # 09
    [3,  0],        # 10
    [2,  0],        # 11
    [7,  4,  3, 0], # 12
    [4,  3,  1, 0], # 13
    [12, 11, 1, 0], # 14
    [1,  0],        # 15
    [5,  3,  2, 0], # 16
    [3,  0],        # 17
    [7,  0],        # 18
    [6,  5,  1, 0], # 19
    [3,  0],        # 20
    [2,  0],        # 21
    [1,  0],        # 22
    [5,  0],        # 23
    [4,  3,  1, 0], # 24
    [3,  0],        # 25
    [8,  7,  1, 0], # 26
    [8,  7,  1, 0], # 27
    [3,  0],        # 28
    [2,  0],        # 29
    [16, 15, 1, 0]  # 30
]
