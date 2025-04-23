module RDF

import ..Accessors: number_of_particles, positions, box_length
import ..Calculations: Calculation, initialize!, compute!, finalize!
import ..MolecularDynamics: VelocityVerlet, vv_integrate!
using UnPack

mutable struct rdf <: Calculation
    hist::Vector{Float64}    # Histogram of pair distances
    bins::Int64             # Number of histogram bins
    dr::Float64             # Bin width
    r_max::Float64          # Maximum radius to consider
    averages::Int64         # Number of samples taken
    N::Int64                 # Number of particles
    L::Float64                 # Box length
    function rdf(; bins::Int64=100, r_max::Float64=0.0)
        new(Vector{Float64}(undef, bins), 
            bins,
            0.0,  # dr will be set in initialize!
            r_max,
            0,
            0,
            0.0)
    end
end

export rdf, rdf_initialize!, rdf_compute!, rdf_finalize!

function rdf_initialize!(calc::rdf, system_::systemtype) where systemtype
    # Set maximum radius to box length if not specified
    calc.L = box_length(system_)
    calc.r_max = box_length(system_)
    calc.N = number_of_particles(system_)
    calc.dr = calc.r_max / calc.bins
    calc.hist .= 0.0
end

function rdf_compute!(calc::rdf, system::systemtype) where systemtype
    r = positions(system)
    @unpack N,L, bins = calc
    
    # Create temporary histogram for this sample
    temp_hist = zeros(Float64, bins)
    R = calc.r_max
    # Loop over unique pairs of particles
    for i in 1:N-1
        rᵢ = r[:,i]
        for j in i+1:N
            Δr = rᵢ - r[:,j]
            Δr = Δr .- L .* round.(Δr ./ L)
            dᵢⱼ = sqrt(sum(Δr.^2))
            
            if dᵢⱼ < R
                bin = floor(Int, dᵢⱼ/calc.dr) + 1
                temp_hist[bin] += 2.0
            end
        end
    end
    
    # Add normalized histogram to running average
    @. calc.hist = (calc.averages * calc.hist + temp_hist) / (calc.averages + 1)
    calc.averages += 1
end

function rdf_finalize!(calc::rdf)
    # Volume of the entire simulation box
    V = calc.L^3
    ρ = calc.N / V  # Number density
    # Normalize histogram to get g(r)
    for i in 1:calc.bins
        r = (i - 0.5) * calc.dr        # Center of bin
        dV = 4π * r^2 * calc.dr        # Volume of spherical shell
        ideal = ρ * dV * calc.N        # Expected number of particles in shell for ideal gas
        
        if ideal > 0
            calc.hist[i] /= ideal
        else
            calc.hist[i] = 0.0
        end
    end
    
    # Optional: Write results to file
    open("rdf.dat", "w") do file
        for i in 1:calc.bins
            r = (i - 0.5) * calc.dr
            println(file, r, " ", calc.hist[i])
        end
    end
end















end
