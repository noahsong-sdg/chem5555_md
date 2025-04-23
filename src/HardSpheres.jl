module HardSpheres
#A Monte Carlo simulation of the hard sphere fluid

using ..MonoAtomic
using UnPack, TOML

export hard_spheres
# Usage here would be the following:
# sys = MonoAtomic.system(N=1000, ϕ=0.5)
# hs = hard_spheres(sys)

#Or:
# hs = hard_spheres("parameters.toml")

#Or:
# hs = hard_spheres(N=1000, ϕ=0.4, lattice="sc")


function packing_fraction_to_density(packing_fraction::Float64)
    # phi to rho
    return 6*packing_fraction/π
end

function packing_fraction(density::Float64)
    # rho to phi
    return π*density/6
end
mutable struct hard_spheres
    system::MonoAtomic.system
    dr::Float64
    acceptance::Float64
    function hard_spheres(sys::MonoAtomic.system)
        sys_ = deepcopy(sys)
        dr_ = 0.05
        acceptance_ = 0.5
        new(sys_, dr_, acceptance_)
    end
    function hard_spheres(; N::Int64=1000, ϕ::Float64=0.4,lattice::String="sc",dr=0.05, acceptance::Float64=0.5)
        ρ_ = packing_fraction_to_density(ϕ)
        sys_ = MonoAtomic.system(N=N, ρ =ρ_, T = 0.0,lattice=lattice)
        new(sys_, dr, acceptance)
    end

    function hard_spheres(filename::String)
        tmp = MonoAtomic.system(filename)
        ϕ_ = tmp.ρ#Interpreting the number density as a packing fraction for hard spheres.
        sys_ = hard_spheres(N=tmp.N, ϕ = ϕ_, lattice=tmp.lattice)
        parameters = TOML.parsefile(filename)
        @unpack dr, acceptance = parameters
        new(sys_.system, dr, acceptance)
    end

end

end
