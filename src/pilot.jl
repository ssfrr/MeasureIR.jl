function pilot(length, freq; samplerate=1)
    L = round(Int, uconvert(NoUnits, length*samplerate))
    dθ = uconvert(NoUnits, 2π*freq/samplerate)

    cos.((0:L-1) * dθ)
end
