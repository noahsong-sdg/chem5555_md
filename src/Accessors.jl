module Accessors

export positions, velocities, forces, box_length, number_of_particles, temperature, kinetic_energy, potential_energy, virial

# Accessing positions for hard_spheres type
function positions(sys::Any)
    # Check if the system is a hard_spheres object
    if hasfield(typeof(sys), :system) && hasfield(typeof(sys.system), :atoms)
        return sys.system.atoms.r
    elseif hasfield(typeof(sys), :atoms)
        return sys.atoms.r
    else
        error("Cannot access positions from system of type $(typeof(sys))")
    end
end

# Similar functions for other accessors
function velocities(sys::Any)
    if hasfield(typeof(sys), :system) && hasfield(typeof(sys.system), :atoms)
        return sys.system.atoms.v
    elseif hasfield(typeof(sys), :atoms)
        return sys.atoms.v
    else
        error("Cannot access velocities from system of type $(typeof(sys))")
    end
end

function forces(sys::Any)
    if hasfield(typeof(sys), :system) && hasfield(typeof(sys.system), :atoms)
        return sys.system.atoms.f
    elseif hasfield(typeof(sys), :atoms)
        return sys.atoms.f
    else
        error("Cannot access forces from system of type $(typeof(sys))")
    end
end

function box_length(sys::Any)
    if hasfield(typeof(sys), :system) && hasfield(typeof(sys.system), :L)
        return sys.system.L
    elseif hasfield(typeof(sys), :L)
        return sys.L
    else
        error("Cannot access box length from system of type $(typeof(sys))")
    end
end

function number_of_particles(sys::Any)
    if hasfield(typeof(sys), :system) && hasfield(typeof(sys.system), :N)
        return sys.system.N
    elseif hasfield(typeof(sys), :N)
        return sys.N
    else
        error("Cannot access number of particles from system of type $(typeof(sys))")
    end
end

function temperature(sys::Any)
    # Get kinetic energy and calculate temperature
    ke = kinetic_energy(sys)
    N = number_of_particles(sys)
    return 2.0 * ke / (3.0 * N)  # T = 2K/3N in reduced units
end

function kinetic_energy(sys::Any)
    v = velocities(sys)
    ke = 0.0
    for i in 1:size(v, 2)
        ke += 0.5 * sum(v[:,i].^2)
    end
    return ke
end

function potential_energy(ff::Any)
    return ff.potential_energy
end

function virial(ff::Any)
    return ff.virial
end

end
