export gen_prop_cut,

@inline function CUDA.exp(x::T) where T <: Complex#complex exponential for GPU
    scale = CUDA.exp( x.re )
    return ComplexF64( scale * CUDA.cos(x.im), scale * CUDA.sin(x.im) )
end

@inline function c_sqrt(x::T) where T <: Complex#complex sqrt
    abs = CUDA.sqrt(x.re^2 + x.im^2)
    cos_x = x.re / abs
    out_re = CUDA.sqrt(abs*(cos_x + 1)/2)
    out_im = CUDA.sqrt(abs*(1 - cos_x)/2)

    out_im *= -sign(x.im)
    return ComplexF64(out_re, out_im)
end

function gen_prop_cut(fields, k0_prop, surface_prop, surface_prop_out,
            k0, dist, R, tilt_x, tilt_y, surfaces, X, Y, Kx, Ky)

    id_x = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    id_y = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    id_z = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    stride_x = blockDim().x * gridDim().x
    stride_y = blockDim().y * gridDim().y
    stride_z = blockDim().z * gridDim().z

    for i in id_x:stride_x:length(X)
        for j in id_y:stride_y:length(Y)

            #disk cutoff along X and Y
            disk_diff::Bool = (X[i]^2 + Y[j]^2 < R^2)

            for z in id_z:stride_z:length(k0)

                #perform initial cut
                fields[z,1,i,j] *= disk_diff
                fields[z,2,i,j] *= disk_diff

                x = CUDA.exp(-1im*k0[z]*tilt_x[z]*X[i])
                y = CUDA.exp(-1im*k0[z]*tilt_y[z]*Y[j])
                s = CUDA.exp(-1im*k0[z]*surfaces[z,i,j])

                #combine surface misalignments and disk cutoff
                surface_prop[z,i,j] = x * y * s * disk_diff

                #New indices to perfom fftshift on k0_prop
                m = (i - 1 + length(Kx)÷2)%length(Kx) + 1
                l = (j - 1 + length(Ky)÷2)%length(Ky) + 1

                χ = conj(c_sqrt(k0[z]^2 - Kx[m]^2 - Ky[l]^2))
                k0_prop[z,i,j] = CUDA.exp(-1im*dist[z]*χ)
            end

            #Calculate phase shifts for the last region without cutoff
            z = length(k0)

            x = CUDA.exp(-1im*k0[z]*tilt_x[z]*X[i])
            y = CUDA.exp(-1im*k0[z]*tilt_y[z]*Y[j])
            s = CUDA.exp(-1im*k0[z]*surfaces[z,i,j])

            surface_prop_out[i,j] = x * y * s
        end
    end
    return nothing
end
