module ForceFields

using ..CellLists: CellList, update_cells!, get_cell_neighbors
using ..Accessors: number_of_particles, positions, forces, box_length
export ForceField, LennardJones, compute_pair_interaction, compute_forces!, update_energy!

abstract type ForceField end

mutable struct LennardJones <: ForceField
    base::Any
    parameters::Dict{Symbol, Any}
    cells::Union{Nothing, CellList}
    potential_energy::Float64
    virial::Float64
end

function LennardJones(base, parameters; cells=nothing)
    LennardJones(base, parameters, cells, 0.0, 0.0)
end

function compute_forces!(ff::LennardJones)
    # Use Accessors functions to get positions, forces, etc.
    r = positions(ff.base)
    f = forces(ff.base)
    N = number_of_particles(ff.base)
    L = box_length(ff.base)
    
    V = W = 0.0
    f .= 0.0
    rij = zeros(3)
    rcut² = isnothing(ff.parameters[:rcut]) ? Inf : ff.parameters[:rcut]^2
    if isa(ff.cells, CellList)
        update_cells!(ff.cells, r, L)
        for i in 1:N
            ri = r[:,i]
            for j in get_cell_neighbors(ff.cells, i, r, L, rcut²)
                rij .= ri - r[:,j]
                rij .= rij .- L .* round.(Int64, rij ./ L)
                r2 = sum(rij .* rij)
                if r2 < rcut²
                    fij, vij = compute_pair_interaction(ff, r2)
                    f[:,i] .+= fij.*rij
                    f[:,j] .-= fij.*rij
                    V += vij
                    W += 0.5*sum(fij.*rij)
                    #W₀ = sum(fij.*rij)
                    #W += i < j ? W₀ : -W₀
                end
            end
        end
    else
        for i in 1:(N-1)
            ri = r[:,i]
            for j in (i+1):N
                rij .= ri - r[:,j]
                rij .= rij .- L .* round.(Int64, rij ./ L)
                r2 = sum(rij .* rij)
                if r2 < rcut²
                    fij, vij = compute_pair_interaction(ff, r2)
                    f[:,i] .+= fij.*rij
                    f[:,j] .-= fij.*rij
                    V += vij
                    W += 0.5*sum(fij.*rij)
                end
            end
        end
    end
    ff.potential_energy = V
    ff.virial = -W
    return
end

function compute_pair_interaction(ff::LennardJones, r2::Float64)
    # Lennard-Jones potential: U(r) = 4ε[(σ/r)¹² - (σ/r)⁶]
    # In reduced units where energy is measured in units of ε and distance in units of σ:
    ε = ff.parameters[:ε]  # Should be 1.0 in reduced units
    σ = ff.parameters[:σ]  # Should be 1.0 in reduced units
    
    r2_inv = σ^2 / r2
    r6_inv = r2_inv^3
    r12_inv = r6_inv^2
    
    # Force magnitude: f = 48ε[(σ/r)¹² - 0.5(σ/r)⁶]/r²
    fij = 48 * ε * (r12_inv - 0.5 * r6_inv) / r2
    
    # Potential energy: U = 4ε[(σ/r)¹² - (σ/r)⁶]
    vij = 4 * ε * (r12_inv - r6_inv)
    
    return fij, vij
end

function update_energy!(ff::LennardJones)
    # Use Accessors functions to get positions
    r = positions(ff.base)
    N = number_of_particles(ff.base)
    L = box_length(ff.base)
    
    V = 0.0
    rij = zeros(3)
    rcut² = isnothing(ff.parameters[:rcut]) ? Inf : ff.parameters[:rcut]^2
    if isa(ff.cells, CellList)
        update_cells!(ff.cells, r, L)
        for i in 1:N
            ri = r[:,i]
            for j in get_cell_neighbors(ff.cells, i, r, L, rcut²)
                rij .= ri - r[:,j]
                rij .= rij .- L .* round.(Int64, rij ./ L)
                r2 = sum(rij .* rij)
                if r2 < rcut²
                    _, vij = compute_pair_interaction(ff, r2)
                    V += vij
                end
            end
        end
    else
        for i in 1:(N-1)
            ri = r[:,i]
            for j in (i+1):N
                rij .= ri - r[:,j]
                rij .= rij .- L .* round.(Int64, rij ./ L)
                r2 = sum(rij .* rij)
                if r2 < rcut²
                    _, vij = compute_pair_interaction(ff, r2)
                    V += vij
                end
            end
        end
    end
    ff.potential_energy = V
    return
end

end
