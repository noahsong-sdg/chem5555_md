module MonteCarlo

# Perform a Monte Carlo simulation of a fluid. Compute functions from the configuraitons 
# and keep averages.

using ..HardSpheres
using ..Accessors
import ..Calculations: initialize!, compute!, finalize!

export sweep!, mc_simulation!, acceptance_criterion

@inline function dot(x,y)
    sum(x.*y)
end

function acceptance_criterion(dᵢⱼ²::Float64)
    # Return 1.0 if the move is accepted and 0.0 otherwise.
    if dᵢⱼ² < 1.0
        return false
    else
        return true
    end
end

function pbc(i::Int64, j::Int64, r::Array{Float64,2}, L::Float64)
    # Calculate the distance between particles i and j with periodic boundary conditions.
    rᵢⱼ = r[:, i] - r[:, j]
    # Apply periodic boundary conditions
    rᵢⱼ = rᵢⱼ .- L * round.(Int64, rᵢⱼ ./L)
    # Compute the squared distance
    dᵢⱼ² = dot(rᵢⱼ, rᵢⱼ)
    return rᵢⱼ, dᵢⱼ²
end

function move!(i::Int64, r::Array{Float64,2}, L::Float64, N::Int64, dr::Float64)
    # Make a trial move on particle i. Accept it or reject it according to Metropolis MC.
    # Compute the distances between particle i and all of its neighbors. If any *Squared* distance is less than dr^2, reject the move.
    # Otherwise, accept it.
    rᵢ = copy(r[:, i])
    Δ = (rand(3) .- 0.5)*dr
    @. r[:, i] += Δ
    accept = 1.0
    # Loop over particles j ≠ i, compute rᵢⱼ and dᵢⱼ². If dᵢⱼ² < 1.0, reject the move, otherwise accept it.
    @inbounds for j in 1:N
        if i != j
            _, dᵢⱼ² = pbc(i, j, r, L)
            if acceptance_criterion(dᵢⱼ²) == false
                # Reject the move
                accept = 0.0
                r[:, i] = rᵢ
                return accept
            end
        end
    end
    # Set accept = 1 if the move is accepted and 0 otherwise.
    return accept
end

function sweep!(sys::hard_spheres)
    # Perform a Monte Carlo sweep, choosing N particles to move at random and attempting moves for each.
    r = positions(sys)
    L = box_length(sys)
    N = number_of_particles(sys)
    acceptance_rate = 0.0
    N = number_of_particles(sys)
    dr = sys.dr
    pick = rand(1:N,N)
    for i ∈ pick
        acceptance_rate += move!(i,r,L,N,dr)
    end
    acceptance_rate /= N
    return acceptance_rate
end

function adjust_dr!(sys::hard_spheres, δ::Float64, mean_rate::Float64)
    sys.dr = sys.dr * (1.0 + δ * (mean_rate - sys.acceptance))
end

function mc_simulation!(; sys::hard_spheres, 
                         sweeps::Int64=10000,
                         calculations::Vector{T}=T[], # Default to empty vector
                         calculate_freq::Int64=100, 
                         adjust_interval::Int64=50, 
                         equil_time::Int64=1000,
                         show_rate::Bool=false,
                         δ::Float64=0.1) where T  # Changed default to 0.1
    # Rest of the function remains the same
    for calc in calculations
        initialize!(calc, sys)
    end
    
    mean_rate = 0.0

    for t in 1:sweeps
        mean_rate += sweep!(sys)

        if t % calculate_freq == 0 && t >= equil_time
            for c in calculations
                compute!(c, sys)
            end
        end

        if t % adjust_interval == 0 && t <= equil_time
            mean_rate = mean_rate/adjust_interval
            if show_rate 
                @show mean_rate
            end
            adjust_dr!(sys, δ, mean_rate)
            mean_rate = 0.0
        end
    end

    for calc in calculations
        finalize!(calc)
    end
end

end