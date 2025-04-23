module CellLists

export CellList, update_cells!, get_cell_neighbors, cell_index  # Added cell_index to exports

mutable struct CellList
    head::Vector{Int64}    # First particle in each cell
    next::Vector{Int64}    # Next particle in linked list
    nc::Vector{Int64}      # Number of cells in each dimension
    cell_size::Float64     # Size of each cell
    n_cells::Int64        # Total number of cells
    neighbors::Vector{Int64}  # Pre-allocated neighbor list
    cell_neighbors::Vector{Vector{Int64}}  # Cached neighboring cells for each cell
end

function CellList(L::Float64, rcut::Float64, N::Int64)
    cell_size = rcut
    nc = floor.(Int, [L/cell_size, L/cell_size, L/cell_size])
    total_cells = prod(nc)
    
    if total_cells <= 27
        @warn "Only $total_cells cells created (≤ 27). Cell lists would be inefficient - reverting to direct calculation."
        return nothing
    end
    
    cell_size = L / nc[1]
    head = zeros(Int64, total_cells)
    next = zeros(Int64, N)
    neighbors = Vector{Int64}()
    
    # Pre-compute neighboring cells (only J > I)
    cell_neighbors = Vector{Vector{Int64}}(undef, total_cells)
    for cell in 1:total_cells
        # Convert linear index to 3D coordinates
        ix = (cell-1) % nc[1] + 1
        iy = ((cell-1) ÷ nc[1]) % nc[2] + 1
        iz = ((cell-1) ÷ (nc[1]*nc[2])) + 1
        
        cell_neighbors[cell] = Int64[]
        
        # Check neighboring cells, only cache those with higher indices
        for dx in -1:1, dy in -1:1, dz in -1:1
            jx = ix + dx
            if jx <= 0 
                jx = nc[1] + jx
            elseif jx > nc[1]
                jx = jx - nc[1]
            end
            
            jy = iy + dy
            if jy <= 0
                jy = nc[2] + jy
            elseif jy > nc[2]
                jy = jy - nc[2]
            end
            
            jz = iz + dz
            if jz <= 0
                jz = nc[3] + jz
            elseif jz > nc[3]
                jz = jz - nc[3]
            end
            
            neighbor_cell = jx + nc[1]*((jy-1) + nc[2]*(jz-1))
            if neighbor_cell > cell  # Only store higher-indexed cells
                push!(cell_neighbors[cell], neighbor_cell)
            end
        end
    end
    
    return CellList(head, next, nc, cell_size, total_cells, neighbors, cell_neighbors)
end

@inline function cell_index(cl::CellList, r::AbstractVector{Float64}, L::Float64)
    # Wrap coordinates LAMMPS-style 
    r_wrapped = zeros(3)
    for k in 1:3
        r_wrapped[k] = r[k]
        while r_wrapped[k] >= L
            r_wrapped[k] -= L
        end
        while r_wrapped[k] < 0
            r_wrapped[k] += L
        end
    end

    # Simple binning to get 1-based indices
    ix = floor(Int, r_wrapped[1]/cl.cell_size)
    iy = floor(Int, r_wrapped[2]/cl.cell_size)
    iz = floor(Int, r_wrapped[3]/cl.cell_size)

    indx = 1 + ix + iy*cl.nc[1] + iz*cl.nc[1]*cl.nc[2]
    return indx
end

function update_cells!(cells::CellList, r::Matrix{Float64}, L::Float64)
    fill!(cells.head, 0)
    for i in 1:size(r,2)
        cell = cell_index(cells, view(r,:,i), L)
        cells.next[i] = cells.head[cell]
        cells.head[cell] = i
    end
end

function get_cell_neighbors(cells::CellList, i::Int64, r::Matrix{Float64}, L::Float64, cutoff²::Float64)
    # This version simply collects candidate neighbor indices.
    empty!(cells.neighbors)
    ri = r[:, i]
    cell_i = cell_index(cells, ri, L)
    
    # For particles in the same cell, add only if j > i.
    j = cells.head[cell_i]
    while j != 0
        if j > i
            push!(cells.neighbors, j)
        end
        j = cells.next[j]
    end
    
    # For neighbor cells, add all particles.
    for cell_j in cells.cell_neighbors[cell_i]
        j = cells.head[cell_j]
        while j != 0
            push!(cells.neighbors, j)
            j = cells.next[j]
        end
    end
    
    return cells.neighbors
end

end
