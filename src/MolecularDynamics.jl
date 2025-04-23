module MolecularDynamics

using ..Accessors: positions, velocities, forces, box_length, number_of_particles, temperature
using ..ForceFields: compute_forces!, update_energy!

export VelocityVerlet, vv_initialize!, vv_integrate!, vv_finalize!, vv_rescale_velocities!

mutable struct VelocityVerlet
    dt::Float64
    system::Any
    forcefield::Any
    target_temperature::Float64
    thermostat_frequency::Int64
    step::Int64
    ke::Float64
    pe::Float64
    virial::Float64
end

function VelocityVerlet(dt, system, forcefield, target_temperature, thermostat_frequency)
    VelocityVerlet(dt, system, forcefield, target_temperature, thermostat_frequency, 0, 0.0, 0.0, 0.0)
end

function vv_initialize!(integrator::VelocityVerlet, system, forcefield, dt, target_temperature, thermostat_frequency)
    # Initial force calculation
    compute_forces!(forcefield)
    
    # Initialize velocities if they're all zero
    v = velocities(system)
    if all(v .== 0.0)
        # Initialize velocities from Maxwell-Boltzmann distribution
        # In reduced units where T* = k_BT/ε, velocities should be sampled
        # from a distribution with variance T* (since mass=1)
        N = number_of_particles(system)
        for i in 1:N
            for j in 1:3
                v[j,i] = sqrt(target_temperature) * randn()
            end
        end
        
        # Remove center of mass motion
        vcm = sum(v, dims=2) / N
        for i in 1:N
            v[:,i] .-= vcm
        end
        
        # Scale to the target temperature
        vv_rescale_velocities!(system, target_temperature)
    end
    
    integrator.dt = dt
    integrator.target_temperature = target_temperature
    integrator.thermostat_frequency = thermostat_frequency
    integrator.step = 0
end

function vv_rescale_velocities!(system, target_temperature)
    v = velocities(system)
    current_temp = temperature(system)
    
    # Avoid division by zero
    if current_temp > 0
        factor = sqrt(target_temperature / current_temp)
        v .*= factor
    end
end

function vv_integrate!(integrator::VelocityVerlet, steps::Int64)
    dt = integrator.dt
    system = integrator.system
    forcefield = integrator.forcefield
    
    r = positions(system)
    v = velocities(system)
    f = forces(system)
    
    for i in 1:steps
        integrator.step += 1
        
        # First half of velocity verlet
        v .+= 0.5 * dt * f
        r .+= dt * v
        
        # Apply PBC
        L = box_length(system)
        r .= r .- L .* floor.(r ./ L)
        
        # Recompute forces at new positions
        compute_forces!(forcefield)
        
        # Second half of velocity verlet
        v .+= 0.5 * dt * f
        
        # Apply thermostat only at specified intervals, not every step
        if integrator.thermostat_frequency > 0 && 
           mod(integrator.step, integrator.thermostat_frequency) == 0
            vv_rescale_velocities!(system, integrator.target_temperature)
        end
        
        
        # Record energy and virial
        integrator.ke = 0.5 * sum(v.^2)
        integrator.pe = forcefield.potential_energy
        integrator.virial = forcefield.virial
        
        # Print progress
        if mod(integrator.step, 100) == 0
            current_temp = temperature(system)
            println("Step: $(integrator.step), Temp: $(current_temp), KE: $(integrator.ke), PE: $(integrator.pe), Virial: $(integrator.virial)")
        end
    end
end

function vv_finalize!(integrator::VelocityVerlet)
    # Final calculations, if any
    println("Simulation completed after $(integrator.step) steps")
end

end
