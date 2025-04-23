module MSD

import ..Accessors: number_of_particles, positions
import ..Calculations: Calculation, initialize!, compute!, finalize!

mutable struct msd <: Calculation
    r₀::Matrix{Float64}
    r²::Vector{Float64}
    t::Int64
    T::Int64
    averages::Int64
    function msd(; T::Int64=100)  # Constructor only sets T
        new(Matrix{Float64}(undef,0,0), 
            Vector{Float64}(undef,0), 
            1, T, 1)
    end
end

export msd, msd_initialize!, msd_compute!, msd_finalize!

function msd_initialize!(calc::msd, system_::systemtype) where systemtype
    r = positions(system_)
    calc.r₀ = deepcopy(r)  # Store initial positions
    calc.r² = zeros(calc.T)  # Use T from constructor
end

@inline function running(; μ::typename, xᵢ::typename, nₐ::Int64) where typename
    return ((nₐ-1)*μ + xᵢ) / nₐ
end

function msd_compute!(calc::msd, system::systemtype) where systemtype
    r = positions(system)
    N = number_of_particles(system)
    
    if calc.t == calc.T + 1
        calc.r₀ = deepcopy(r)
        calc.t = 1
        calc.averages += 1
    end

    r² = 0.0
    Δr = r - calc.r₀
    
    for Δrᵢ in eachcol(Δr)
        r² += sum(Δrᵢ.*Δrᵢ)
    end
    r² /= N
    
    calc.r²[calc.t] = running(μ = calc.r²[calc.t], xᵢ = r², nₐ=calc.averages)
    
    calc.t += 1
end

function msd_finalize!(calc::msd) 
    open("msd.dat", "w") do file
        for value in calc.r²
            println(file, value)
        end
    end
end

end
