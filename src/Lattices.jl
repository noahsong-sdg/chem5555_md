module Lattices
# A module for generating the three crystal lattices with a cubic unit cell.

export bcc_lattice, fcc_lattice, sc_lattice

function bcc_lattice(N::Int64, L::Float64)
    # BCC is a cubic lattice with 2 atoms per unit cell.
    nside = round(Int64, (N/2)^(1/3) )
    @show 2*nside^3 N
    if !isapprox(2*nside^3, N)
        xerr = 2*nside^3
        error("N/2 must be a perfect cube. Try $xerr")
    end
    positions = zeros(Float64, 3, N)
    
    a = L/nside  # lattice constant
    count = 1
    
    for i in 0:nside-1
        for j in 0:nside-1
            for k in 0:nside-1
                if count > N break end
                # Corner atom
                positions[:, count] = [i*a, j*a, k*a]
                count += 1
                
                if count > N break end
                # Center atom
                positions[:, count] = [(i + 0.5)*a, (j + 0.5)*a, (k + 0.5)*a]
                count += 1
            end
        end
    end
    return positions
end

function sc_lattice(N::Int64, L::Float64)
    nside = round(Int64, N^(1/3) )
    if !isapprox(nside^3, N)
        xerr = nside^3
        error("N must be a perfect cube. Try $xerr")
    end
    positions = zeros(Float64, 3, N)
    
    a = L/nside  # lattice constant
    count = 1
    
    for i in 0:nside-1
        for j in 0:nside-1
            for k in 0:nside-1
                if count > N break end
                positions[:, count] = [i*a, j*a, k*a]
                count += 1
            end
        end
    end
    return positions
end

function fcc_lattice(N::Int64, L::Float64)
    nside = round(Int64, (N/4)^(1/3) )
    if !isapprox(4*nside^3, N)
        xerr = 4*nside^3
        error("N/4 must be a perfect cube. Try $xerr")
    end
    positions = zeros(Float64, 3, N)
    
    a = L/nside  # lattice constant
    count = 1
    
    for i in 0:nside-1
        for j in 0:nside-1
            for k in 0:nside-1
                if count > N break end
                # Corner atom
                positions[:, count] = [i*a, j*a, k*a]
                count += 1
                
                if count > N break end
                # Face-centered atoms
                positions[:, count] = [(i + 0.5)*a, (j + 0.5)*a, k*a]
                count += 1
                if count > N break end
                positions[:, count] = [(i + 0.5)*a, j*a, (k + 0.5)*a]
                count += 1
                if count > N break end
                positions[:, count] = [i*a, (j + 0.5)*a, (k + 0.5)*a]
                count += 1
            end
        end
    end
    return positions
end

end