module MD_Simulation

# Include files from the correct location
try
    include(joinpath(@__DIR__, "src", "MonoAtomic.jl"))
    include(joinpath(@__DIR__, "src", "HardSpheres.jl"))  # This is directly in Module_4
    include(joinpath(@__DIR__, "src", "Accessors.jl"))
    include(joinpath(@__DIR__, "src", "Calculations.jl"))
    include(joinpath(@__DIR__, "src", "MSD.jl"))
    include(joinpath(@__DIR__, "src", "MonteCarlo.jl"))
    include(joinpath(@__DIR__, "src", "CellLists.jl"))
    include(joinpath(@__DIR__, "src", "ForceFields.jl"))
    include(joinpath(@__DIR__, "src", "MolecularDynamics.jl"))
    include(joinpath(@__DIR__, "src", "RDF.jl"))
    include(joinpath(@__DIR__, "src", "utils.jl"))
catch e
    println("Error including files: $e")
    error("Failed to include required module files")
end

using .MonoAtomic
using .HardSpheres
using .Accessors
using .Calculations
using .MSD
using .RDF
using .MonteCarlo
using .CellLists
using .ForceFields
using .MolecularDynamics
using .MD_utils
# Export the types and functions needed
export msd, rdf

# from accessors, celllists, forcefields, [lattices unneeded, its called by other mods],
# then moleculardynamics, 
end
