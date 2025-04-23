module MD_utils

import ..Accessors: number_of_particles, positions, box_length
import ..Calculations: Calculation, initialize!, compute!, finalize!
import ..RDF: rdf, rdf_initialize!, rdf_compute!, rdf_finalize!
import ..MolecularDynamics: VelocityVerlet, vv_integrate!
using UnPack, PlotlyBase, Plots, DelimitedFiles

export plot_particles_with_box
################ Visualizations ############################
function plot_particles_with_box(positions, L)
    x = positions[1, :]
    y = positions[2, :]
    z = positions[3, :]
    
    p = Plots.scatter3d(x, y, z, 
        markersize=2,
        label="Particles",
        xlabel="X", ylabel="Y", zlabel="Z",
        xlim=(0,L), ylim=(0,L), zlim=(0,L))
    
    # to make it more interactive
    plot!(p, camera=(45, 45))  # rotate view
end

function plot_lammps_rdf(rdf_file="rdf.txt")
    # Read the RDF data
    data = readdlm(rdf_file, comments=true)
    
    # LAMMPS RDF compute output format:
    # First line of each chunk has: Timestep Number_of_bins
    # Then Number_of_bins lines with: Bin_index R g(R)
    
    # Find the last complete RDF dataset
    timestep_rows = findall(x -> length(x) == 2, eachrow(data))
    if isempty(timestep_rows)
        error("No RDF data found in file")
    end
    
    # Get the last timestep row
    last_timestep_row = timestep_rows[end]
    num_bins = Int(data[last_timestep_row, 2])
    
    # Extract the RDF data for the last timestep
    rdf_start = last_timestep_row + 1
    rdf_end = rdf_start + num_bins - 1
    
    if rdf_end > size(data, 1)
        error("Incomplete RDF data at the end of file")
    end
    
    # Extract R values and g(R) values
    r_values = data[rdf_start:rdf_end, 2]
    gr_values = data[rdf_start:rdf_end, 3]
    
    # Plot the RDF
    plt = plot(
        r_values, gr_values,
        xlabel="Distance r (Å)",
        ylabel="g(r)",
        title="Radial Distribution Function",
        linewidth=2,
        legend=false,
        grid=true
    )
    
    # Add vertical lines at expected peaks for LJ fluid
    vline!([1.0, 2.0], linestyle=:dash, color=:gray, alpha=0.5, label="")
    
    display(plt)
    savefig(plt, "rdf_plot.png")
    
    return plt, r_values, gr_values
end


function run_md_with_rdf!(integrator::VelocityVerlet, rdf_calc::rdf, n_steps::Int, rdf_steps::Int)
    for i in 1:n_steps
        vv_integrate!(integrator, 1)
        if i % rdf_steps == 0
            rdf_compute!(rdf_calc, integrator.system)
        end
    end
end
end
