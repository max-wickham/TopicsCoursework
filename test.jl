using Images, FileIO, Wavelets, ImageView
using Plots
img_path = "MRI_of_Human_Brain.jpg"
X0 = load(img_path)
X0 = X0[1:256,1:256]
X0 = Float64.(X0)

using FFTW
using StatsBase
F0 = fftshift((fft(Float64.(X0))))

using Random
Random.seed!(1);

sampling_indexes = sort(sample(collect(1:size(X0)[1]*size(X0)[1]), round(Int64, 0.3 * size(X0)[1] * size(X0)[2])))

wavlet_filt = wavelet(WT.db2, WT.Filter)

using LinearAlgebra

############ Q9

function P(x)
    result = zeros(ComplexF64,256,256)
    for index in sampling_indexes
        row = (index % 255) + 1
        column = ceil(Int64, index / 256)
        result[row, column] = x[row,column]
    end
    return result
end

y = P(F0);

function ∇f(x)
    Fx = fftshift(fft(x))
    vecF0 = P(Fx)
    normal = vecF0 - y
    normal = ifft(ifftshift(normal))
    return normal
end

function f(x)
    f = 0.5*norm(y - P(fftshift(fft(x))))
    return f
end

function fhat(x,z, γ)
    result = f(z)
    result += tr((x-z)'*real.(∇f(z)))
    result += (1/(2*γ))*(norm(x-y, 2)^2)
    return result
end


################################

using Statistics
MINX = zeros(256,256)/100
function soft(a, c)
    if abs(a) > c
        return sign(a) * (abs(a) - c)
    end
        return 0.0
end
for i in 1:10
    global MINX
    η = 0.9
    # γ = 0.00675
    γ = 0.01
    λ = 0.2
    c = γ * λ
    z = []
    while true
        z = ∇f(MINX)
        z = MINX - γ*real.(z)
        z = dwt(z, wavlet_filt, 2)
        z = soft.(z, c)
        z = idwt(z, wavlet_filt, 2)
        if f(z) <= real.(fhat(z,MINX,γ))
            break
        end
        γ = η * γ
        c = γ * λ
    end
    MINX = z
end

# plot(Gray.(100*MINX))
print(maximum(10*MINX))
