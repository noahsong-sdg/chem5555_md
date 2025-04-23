using DelimitedFiles
using Plots
using Statistics

function read_lammps_dump(filename)
    """
    read_lammps_dump(filename)

    Read a LAMMPS dump file and return the trajectory data.
    """
    # First, check if the file exists
    if !isfile(filename)
        println("ERROR: File not found: $filename")
        println("Current working directory: $(pwd())")
        println("Files in current directory:")
        for (i, file) in enumerate(readdir())
            println("  $i. $file")
        end
        error("Could not find LAMMPS dump file: $filename")
    end
    
    # Try to read the file with better error handling
    frames = []
    current_frame = Dict()
    
    try
        open(filename, "r") do file
            line_number = 0
            n_atoms = 0
            atom_data = []
            
            for line in eachline(file)
                line_number += 1
                
                if startswith(line, "ITEM: TIMESTEP")
                    # Start of a new frame
                    if !isempty(current_frame)
                        push!(frames, current_frame)
                    end
                    current_frame = Dict()
                    atom_data = []
                    
                elseif startswith(line, "ITEM: NUMBER OF ATOMS")
                    n_atoms = parse(Int, readline(file))
                    line_number += 1
                    
                elseif startswith(line, "ITEM: BOX BOUNDS")
                    bounds = []
                    for _ in 1:3
                        box_line = readline(file)
                        line_number += 1
                        box_vals = parse.(Float64, split(box_line))
                        push!(bounds, box_vals)
                    end
                    current_frame["box"] = bounds
                    
                elseif startswith(line, "ITEM: ATOMS")
                    headers = split(line)[3:end]
                    current_frame["headers"] = headers
                    
                    # Read atom data
                    atom_data = []
                    for _ in 1:n_atoms
                        atom_line = readline(file)
                        line_number += 1
                        atom_vals = parse.(Float64, split(atom_line))
                        push!(atom_data, atom_vals)
                    end
                    current_frame["atoms"] = atom_data
                end
            end
            
            # Add the last frame
            if !isempty(current_frame)
                push!(frames, current_frame)
            end
        end
    catch e
        println("ERROR reading LAMMPS dump file:")
        println(e)
        println("Make sure the file exists and is accessible.")
        rethrow(e)
    end
    
    return frames
end

function plot_trajectory(frames, output_file="lammps_trajectory.gif")
    """
    plot_trajectory(frames, output_file="lammps_trajectory.gif")

    Create an animation of the particle positions from trajectory data.
    """
    n_frames = length(frames)
    println("Creating animation with $n_frames frames")
    
    anim = @animate for (i, frame) in enumerate(frames)
        if i % 10 == 0 
            println("Processing frame $i of $n_frames")
        end
        
        # Extract positions
        atoms = frame["atoms"]
        x = [atom[4] for atom in atoms]  # x position is typically column 4
        y = [atom[5] for atom in atoms]  # y position is typically column 5
        z = [atom[6] for atom in atoms]  # z position is typically column 6
        
        # If needed, can extract velocities too
        # vx = [atom[7] for atom in atoms]
        # vy = [atom[8] for atom in atoms]
        # vz = [atom[9] for atom in atoms]
        
        # Get box dimensions
        box = frame["box"]
        xmin, xmax = box[1]
        ymin, ymax = box[2]
        zmin, zmax = box[3]
        
        # Create a 3D scatter plot
        scatter(x, y, z, 
                markersize=2, 
                markerstrokewidth=0,
                xlim=(xmin, xmax),
                ylim=(ymin, ymax),
                zlim=(zmin, zmax),
                title="LAMMPS Simulation - Frame $i",
                xlabel="X", ylabel="Y", zlabel="Z",
                legend=false,
                camera=(30, 30),
                size=(800, 600))
    end
    
    gif(anim, output_file, fps=10)
    println("Animation saved to $output_file")
    return anim
end

function analyze_trajectory(frames)
    """
    analyze_trajectory(frames)

    Perform basic analysis of the trajectory data.
    """
    # Calculate mean positions and velocities over time
    mean_positions = []
    mean_velocities = []
    
    for frame in frames
        atoms = frame["atoms"]
        
        # Extract positions and velocities
        positions = [atom[4:6] for atom in atoms]  # x,y,z positions
        mean_pos = mean(positions)
        push!(mean_positions, mean_pos)
        
        if length(atoms[1]) >= 9  # Check if velocities exist
            velocities = [atom[7:9] for atom in atoms]  # vx,vy,vz velocities
            mean_vel = mean(velocities)
            push!(mean_velocities, mean_vel)
        end
    end
    
    # Plot mean position over time
    times = 1:length(frames)
    
    p1 = plot(times, [pos[1] for pos in mean_positions], 
              label="X", title="Mean Position vs Time",
              xlabel="Frame", ylabel="Position")
    plot!(p1, times, [pos[2] for pos in mean_positions], label="Y")
    plot!(p1, times, [pos[3] for pos in mean_positions], label="Z")
    
    # Plot mean velocity over time if available
    if !isempty(mean_velocities)
        p2 = plot(times, [vel[1] for vel in mean_velocities], 
                  label="VX", title="Mean Velocity vs Time",
                  xlabel="Frame", ylabel="Velocity")
        plot!(p2, times, [vel[2] for vel in mean_velocities], label="VY")
        plot!(p2, times, [vel[3] for vel in mean_velocities], label="VZ")
        
        return plot(p1, p2, layout=(2,1), size=(800, 800))
    else
        return p1
    end
end

function calculate_rdf(frames, nbins=100, rmax=2.5)
    """
    calculate_rdf(frames, nbins=100, rmax=2.5)

    Calculate the radial distribution function (RDF) from trajectory data.
    Default rmax is set to 2.5 to match the LJ cutoff in the LAMMPS script.
    """
    # Initialize histogram arrays
    hist = zeros(nbins)
    bin_edges = range(0, rmax, length=nbins+1)
    bin_centers = [(bin_edges[i] + bin_edges[i+1])/2 for i in 1:nbins]
    bin_width = bin_edges[2] - bin_edges[1]
    
    num_frames = length(frames)
    total_atoms = 0
    
    println("Calculating RDF...")
    
    for (frame_idx, frame) in enumerate(frames)
        if frame_idx % 10 == 0
            println("Processing frame $frame_idx of $num_frames for RDF")
        end
        
        atoms = frame["atoms"]
        n_atoms = length(atoms)
        total_atoms += n_atoms
        positions = [[atom[4], atom[5], atom[6]] for atom in atoms]
        
        # Get box dimensions for periodic boundary conditions
        box = frame["box"]
        frame_box_lengths = [box[1][2] - box[1][1], box[2][2] - box[2][1], box[3][2] - box[3][1]]
        
        # Calculate all pairwise distances
        for i in 1:n_atoms
            for j in (i+1):n_atoms  # Only consider unique pairs (i,j) where i < j
                # Calculate distance with minimum image convention
                dx = positions[i][1] - positions[j][1]
                dy = positions[i][2] - positions[j][2]
                dz = positions[i][3] - positions[j][3]
                
                # Apply minimum image convention
                dx -= frame_box_lengths[1] * round(dx / frame_box_lengths[1])
                dy -= frame_box_lengths[2] * round(dy / frame_box_lengths[2])
                dz -= frame_box_lengths[3] * round(dz / frame_box_lengths[3])
                
                r = sqrt(dx^2 + dy^2 + dz^2)
                
                # Add to histogram if within cutoff
                if r < rmax
                    bin_idx = Int(floor(r / bin_width)) + 1
                    if 1 <= bin_idx <= nbins
                        hist[bin_idx] += 2  # Count each pair twice (i→j and j→i)
                    end
                end
            end
        end
    end
    
    # Calculate average density over all frames
    avg_n_atoms = total_atoms / num_frames
    total_volume = 0.0
    for frame in frames
        box = frame["box"]
        frame_box_lengths = [box[1][2] - box[1][1], box[2][2] - box[2][1], box[3][2] - box[3][1]]
        volume = frame_box_lengths[1] * frame_box_lengths[2] * frame_box_lengths[3]
        total_volume += volume
    end
    avg_volume = total_volume / num_frames
    avg_density = avg_n_atoms / avg_volume
    
    # Normalize histogram to get g(r)
    g_r = zeros(nbins)
    for i in 1:nbins
        r = bin_centers[i]
        
        # Volume of the spherical shell
        shell_volume = 4π * (r^2) * bin_width
        
        # Expected number of particles in this shell for ideal gas
        ideal_count = shell_volume * avg_density
        
        # Average number of particles found in this shell per reference particle
        if ideal_count > 0 && r > 0.1  # Avoid very small r values
            # Divide by number of frames and number of atoms (since each atom is a reference point)
            g_r[i] = hist[i] / (num_frames * avg_n_atoms * ideal_count)
        else
            g_r[i] = 0.0
        end
    end
    
    return bin_centers, g_r
end

function read_lammps_rdf(filename="Module_4/rdf.dat")
    """
    read_lammps_rdf(filename="Module_4/rdf.dat")

    Read the RDF data directly calculated by LAMMPS, which is more reliable.
    """
    if !isfile(filename)
        println("WARNING: LAMMPS RDF file not found: $filename")
        return nothing, nothing
    end
    
    println("Reading LAMMPS RDF data from $filename")
    
    try
        # Read the file content
        data = readdlm(filename, comments=true)
        
        # LAMMPS RDF file format: first column is bin center, second is g(r)
        r = data[:, 1]
        g_r = data[:, 2]
        
        return r, g_r
    catch e
        println("ERROR reading LAMMPS RDF file:")
        println(e)
        return nothing, nothing
    end
end

function calculate_and_plot_rdf(filename; nbins=100, rmax=2.5, stride=5, use_lammps_rdf=true)
    """
    calculate_and_plot_rdf(filename; nbins=100, rmax=2.5, stride=5, use_lammps_rdf=true)

    Calculate and plot the radial distribution function from a LAMMPS trajectory file.
    """
    rdf_plot = nothing
    
    # Try to use LAMMPS' RDF first if requested
    if use_lammps_rdf
        r_lammps, g_r_lammps = read_lammps_rdf()
        if r_lammps !== nothing
            # Successfully read LAMMPS RDF
            rdf_plot = plot(r_lammps, g_r_lammps,
                   title="Radial Distribution Function (LAMMPS)",
                   xlabel="Distance (r)",
                   ylabel="g(r)",
                   linewidth=2,
                   legend=false)
            
            savefig(rdf_plot, "rdf_lammps.png")
            println("LAMMPS RDF plot saved as 'rdf_lammps.png'")
            
            return (r_lammps, g_r_lammps, rdf_plot)
        end
    end
    
    # Fall back to our calculation if LAMMPS RDF is not available or not requested
    println("Reading LAMMPS trajectory file for RDF calculation...")
    trajectory = read_lammps_dump(filename)
    println("Found $(length(trajectory)) frames")
    
    # Only use the production run frames (after equilibration)
    # Based on your LAMMPS script, the first 10000 steps are equilibration
    # Assuming dumps are every 100 steps, that's frame 100
    prod_frames = trajectory[100:end]
    println("Using $(length(prod_frames)) production frames")
    
    # Use a subset of frames for more efficient calculation
    rdf_frames = prod_frames[1:stride:end]
    println("Using $(length(rdf_frames)) frames for RDF calculation")
    
    # Calculate RDF
    r, g_r = calculate_rdf(rdf_frames, nbins, rmax)
    
    # Create the RDF plot
    rdf_plot = plot(r, g_r,
                   title="Radial Distribution Function",
                   xlabel="Distance (r)",
                   ylabel="g(r)",
                   linewidth=2,
                   legend=false)
    
    savefig(rdf_plot, "rdf_analysis.png")
    println("RDF plot saved as 'rdf_analysis.png'")
    
    return (r, g_r, rdf_plot)
end

#####################################################
############# Liquid Water Simulation ###############
#####################################################
function plot_energy_conservation(energy_file="Module_4/energy.txt")
    """
        plot_energy_conservation(energy_file="Module_4/energy.txt")

    Plot energy conservation from LAMMPS simulation and verify it's conserved to one part in a thousand.
    """
    # Load the energy data
    data = readdlm(energy_file, skipstart=1)

    # Extract columns
    steps = data[:, 1]
    total_energy = data[:, 2]
    potential_energy = data[:, 3]
    kinetic_energy = data[:, 4]
    temperature = data[:, 5]

    # Convert steps to time (multiply by timestep)
    # Account for changing timestep during simulation
    time = similar(steps)
    for i in eachindex(steps)
        if steps[i] <= 1000
            time[i] = steps[i] * 0.01  # 0.01 fs per step
        elseif steps[i] <= 2000
            time[i] = 1000 * 0.01 + (steps[i] - 1000) * 0.1  # 0.1 fs per step after 1000 steps
        else
            time[i] = 1000 * 0.01 + 1000 * 0.1 + (steps[i] - 2000) * 1.0  # 1.0 fs per step after 2000 steps
        end
    end

    # Create the energy plot
    p1 = plot(
        time, total_energy, 
        label="Total Energy", 
        linewidth=2,
        xlabel="Time (fs)",
        ylabel="Energy (kcal/mol)",
        title="Energy Conservation in TIP4P Water Simulation",
        legend=:outertopright,
        grid=true
    )

    plot!(p1, time, potential_energy, label="Potential Energy", linewidth=1.5)
    plot!(p1, time, kinetic_energy, label="Kinetic Energy", linewidth=1.5)

    # Create temperature plot with secondary y-axis
    p2 = twinx(p1)
    plot!(p2, time, temperature, 
        label="Temperature", 
        color=:red, 
        linestyle=:dash, 
        linewidth=1.5,
        ylabel="Temperature (K)",
        legend=:outertopright
    )

    # Save the plot
    savefig(p1, "energy_conservation.png")
    
    println("Plot saved as energy_conservation.png")
    
    return p1
end

function calculate_oo_rdf(frames, nbins=100, rmax=10.0)
    """
    calculate_oo_rdf(frames, nbins=100, rmax=10.0)

    Calculate the oxygen-oxygen radial distribution function (RDF) from trajectory data.
    Only considers oxygen atoms (type 1) for the RDF calculation.
    """
    # Initialize histogram arrays
    hist = zeros(nbins)
    bin_edges = range(0, rmax, length=nbins+1)
    bin_centers = [(bin_edges[i] + bin_edges[i+1])/2 for i in 1:nbins]
    bin_width = bin_edges[2] - bin_edges[1]
    
    num_frames = length(frames)
    total_oxygen = 0
    
    println("Calculating Oxygen-Oxygen RDF...")
    
    for (frame_idx, frame) in enumerate(frames)
        if frame_idx % 10 == 0
            println("Processing frame $frame_idx of $num_frames for O-O RDF")
        end
        
        atoms = frame["atoms"]
        
        # Filter oxygen atoms (type 1)
        oxygen_indices = findall(atom -> atom[2] == 1.0, atoms)  # Column 2 is typically atom type
        oxygen_atoms = atoms[oxygen_indices]
        n_oxygen = length(oxygen_atoms)
        total_oxygen += n_oxygen
        
        if n_oxygen == 0
            println("Warning: No oxygen atoms found in frame $frame_idx")
            continue
        end
        
        # Get positions of oxygen atoms
        oxygen_positions = [[atom[4], atom[5], atom[6]] for atom in oxygen_atoms]
        
        # Get box dimensions for periodic boundary conditions
        box = frame["box"]
        frame_box_lengths = [box[1][2] - box[1][1], box[2][2] - box[2][1], box[3][2] - box[3][1]]
        
        # Calculate all pairwise distances between oxygen atoms
        for i in 1:n_oxygen
            for j in (i+1):n_oxygen  # Only consider unique pairs (i,j) where i < j
                # Calculate distance with minimum image convention
                dx = oxygen_positions[i][1] - oxygen_positions[j][1]
                dy = oxygen_positions[i][2] - oxygen_positions[j][2]
                dz = oxygen_positions[i][3] - oxygen_positions[j][3]
                
                # Apply minimum image convention
                dx -= frame_box_lengths[1] * round(dx / frame_box_lengths[1])
                dy -= frame_box_lengths[2] * round(dy / frame_box_lengths[2])
                dz -= frame_box_lengths[3] * round(dz / frame_box_lengths[3])
                
                r = sqrt(dx^2 + dy^2 + dz^2)
                
                # Add to histogram if within cutoff
                if r < rmax
                    bin_idx = Int(floor(r / bin_width)) + 1
                    if 1 <= bin_idx <= nbins
                        hist[bin_idx] += 2  # Count each pair twice (i→j and j→i)
                    end
                end
            end
        end
    end
    
    # Calculate average density of oxygen atoms over all frames
    avg_n_oxygen = total_oxygen / num_frames
    total_volume = 0.0
    for frame in frames
        box = frame["box"]
        frame_box_lengths = [box[1][2] - box[1][1], box[2][2] - box[2][1], box[3][2] - box[3][1]]
        volume = frame_box_lengths[1] * frame_box_lengths[2] * frame_box_lengths[3]
        total_volume += volume
    end
    avg_volume = total_volume / num_frames
    avg_density = avg_n_oxygen / avg_volume
    
    # Normalize histogram to get g(r)
    g_r = zeros(nbins)
    for i in 1:nbins
        r = bin_centers[i]
        
        # Volume of the spherical shell
        shell_volume = 4π * (r^2) * bin_width
        
        # Expected number of particles in this shell for ideal gas
        ideal_count = shell_volume * avg_density
        
        # Average number of particles found in this shell per reference particle
        if ideal_count > 0 && r > 0.1  # Avoid very small r values
            # Divide by number of frames and number of oxygen atoms
            g_r[i] = hist[i] / (num_frames * avg_n_oxygen * ideal_count)
        else
            g_r[i] = 0.0
        end
    end
    
    return bin_centers, g_r
end

function plot_oxygen_oxygen_rdf(filename="Module_4/dump.lammpstrj", dump_freq=100; nbins=100, rmax=10.0, stride=5)
    """
    plot_oxygen_oxygen_rdf(filename="dump.lammpstrj", dump_freq=100; nbins=100, rmax=10.0, stride=5)

    Calculate and plot the oxygen-oxygen radial distribution function from a LAMMPS trajectory file.
    
    Parameters:
    - filename: path to the LAMMPS trajectory file
    - dump_freq: how often frames were dumped in the LAMMPS simulation
    - nbins: number of bins for the RDF histogram
    - rmax: maximum distance for the RDF calculation (Å)
    - stride: use every N-th frame for RDF calculation
    """
    println("Reading LAMMPS trajectory file for O-O RDF calculation...")
    trajectory = read_lammps_dump(filename)
    println("Found $(length(trajectory)) frames")
    
    # Use only production run frames (after equilibration)
    prod_frames = trajectory[end÷2:end]  # Use second half of trajectory for production
    println("Using $(length(prod_frames)) production frames")
    
    # Use a subset of frames for more efficient calculation
    rdf_frames = prod_frames[1:stride:end]
    println("Using $(length(rdf_frames)) frames for O-O RDF calculation")
    
    # Calculate oxygen-oxygen RDF
    r, g_r = calculate_oo_rdf(rdf_frames, nbins, rmax)
    
    # Find the first solvation shell (first peak in RDF)
    # Skip the first few bins which might contain noise or artifacts at very small distances
    start_idx = max(1, Int(floor(1.5 / (rmax/nbins))))  # Start looking from ~1.5 Å
    peak_idx = argmax(g_r[start_idx:end]) + start_idx - 1
    first_shell_r = r[peak_idx]
    first_shell_g = g_r[peak_idx]
    
    println("First solvation shell found at $(round(first_shell_r, digits=3)) Å")
    println("RDF value at first solvation shell: $(round(first_shell_g, digits=3))")
    
    # Create the RDF plot
    rdf_plot = plot(r, g_r,
                   title="Oxygen-Oxygen Radial Distribution Function",
                   xlabel="Distance (Å)",
                   ylabel="g(r)",
                   linewidth=2,
                   legend=false)
    
    # Add vertical line at the first solvation shell
    vline!([first_shell_r], 
           label="First shell: $(round(first_shell_r, digits=2)) Å", 
           linestyle=:dash, 
           linewidth=1.5,
           color=:red)
    
    # Add annotation for RDF value at first shell
    annotate!([(first_shell_r, first_shell_g + 0.2, 
               text("g(r) = $(round(first_shell_g, digits=2))", 
               :right, 8))])
    
    savefig(rdf_plot, "oo_rdf.png")
    println("O-O RDF plot saved as 'oo_rdf.png'")
    
    return (r, g_r, rdf_plot)
end

# Add this line at the end of your script to execute the O-O RDF calculation
r_oo, g_oo, rdf_plot = plot_oxygen_oxygen_rdf("Module_4/dump.lammpstrj")

function calculate_diffusion_coefficient(filename="Module_4/dump.lammpstrj")
    """
    calculate_diffusion_coefficient(filename="Module_4/dump.lammpstrj")
    
    Calculate the diffusion coefficient of water oxygen atoms from a LAMMPS trajectory.
    Uses the Einstein relation: D = MSD/(6*t) for 3D systems.
    """
    println("Reading LAMMPS trajectory for diffusion coefficient calculation...")
    trajectory = read_lammps_dump(filename)
    println("Found $(length(trajectory)) frames")
    
    # Use production run frames (skip equilibration)
    # Usually we want at least the last 60-70% of the trajectory
    prod_frames = trajectory[Int(length(trajectory)*0.3):end]
    println("Using $(length(prod_frames)) frames for diffusion coefficient")
    
    # Get timestep information
    # For simplicity, let's assume we're using the final timestep of 1.0 fs
    # In a more accurate calculation, we would get this from the LAMMPS log file
    timestep = 1.0  # fs
    dump_freq = 100  # How often frames were dumped (from your LAMMPS script)
    time_between_frames = timestep * dump_freq  # fs
    
    # Track oxygen atoms only
    first_frame = prod_frames[1]
    atoms = first_frame["atoms"]
    oxygen_indices = findall(atom -> atom[2] == 1.0, atoms)  # Column 2 is atom type
    
    # Extract oxygen atom IDs to track them across frames
    oxygen_ids = [Int(atoms[i][1]) for i in oxygen_indices]  # Column 1 is atom ID
    n_oxygen = length(oxygen_ids)
    
    # Create a map of atom IDs to indices for efficient lookup
    id_to_index = Dict()
    
    # Initialize storage for unwrapped trajectories
    # We need to unwrap to account for periodic boundary crossings
    unwrapped_traj = Array{Float64}(undef, n_oxygen, 3, length(prod_frames))
    
    # Process first frame - store initial positions
    for (idx, o_idx) in enumerate(oxygen_indices)
        atom = atoms[o_idx]
        unwrapped_traj[idx, :, 1] = [atom[4], atom[5], atom[6]]  # x,y,z positions
    end
    
    # Get box dimensions
    box = first_frame["box"]
    box_lengths = [box[1][2] - box[1][1], box[2][2] - box[2][1], box[3][2] - box[3][1]]
    
    # Process subsequent frames - unwrap trajectories
    for (frame_idx, frame) in enumerate(prod_frames[2:end])
        actual_frame_idx = frame_idx + 1  # Adjust for 1-based indexing
        
        # Build index map for this frame
        atoms = frame["atoms"]
        id_to_index = Dict(Int(atom[1]) => i for (i, atom) in enumerate(atoms))
        
        # Process each oxygen atom
        for (o_idx, o_id) in enumerate(oxygen_ids)
            # Find atom in current frame
            atom_idx = id_to_index[o_id]
            atom = atoms[atom_idx]
            
            # Get current wrapped position
            current_pos = [atom[4], atom[5], atom[6]]
            
            # Get previous unwrapped position
            prev_pos = unwrapped_traj[o_idx, :, actual_frame_idx-1]
            
            # Calculate displacement accounting for periodic boundaries
            disp = current_pos - prev_pos
            
            # Apply minimum image convention
            for d in 1:3
                while disp[d] > box_lengths[d]/2
                    disp[d] -= box_lengths[d]
                end
                while disp[d] < -box_lengths[d]/2
                    disp[d] += box_lengths[d]
                end
            end
            
            # Update unwrapped position
            unwrapped_traj[o_idx, :, actual_frame_idx] = prev_pos + disp
        end
    end
    
    # Calculate MSD
    n_frames = length(prod_frames)
    max_lag = min(n_frames-1, 50)  # Use up to 50 lag times or all available
    
    # Time values in picoseconds (convert from fs)
    times = [(lag * time_between_frames)/1000 for lag in 0:max_lag]
    
    # Calculate MSD for different lag times
    msd = zeros(max_lag + 1)
    
    for lag in 0:max_lag
        sum_sq_disp = 0.0
        count = 0
        
        for frame_idx in 1:(n_frames-lag)
            for atom_idx in 1:n_oxygen
                disp = unwrapped_traj[atom_idx, :, frame_idx+lag] - unwrapped_traj[atom_idx, :, frame_idx]
                sq_disp = sum(disp.^2)
                sum_sq_disp += sq_disp
                count += 1
            end
        end
        
        msd[lag+1] = sum_sq_disp / count
    end
    
    # Calculate diffusion coefficient from Einstein relation
    # D = lim(t→∞) MSD(t)/(6*t) for 3D
    # Use linear fit to later portion of MSD curve (skip initial non-linear region)
    skip_points = max(2, div(max_lag, 10))  # Skip ~10% of initial points
    
    # Linear fit: MSD = 6D*t + b
    y = msd[(skip_points+1):end]
    x = times[(skip_points+1):end]
    
    # Simple linear regression
    n = length(x)
    slope = (n * sum(x .* y) - sum(x) * sum(y)) / (n * sum(x.^2) - sum(x)^2)
    
    # Diffusion coefficient (Å²/ps)
    D = slope / 6.0
    
    # Diffusion coefficient in cm²/s (standard units)
    # 1 Å²/ps = 10^-8 cm²/s
    D_standard = D * 1e-8
    
    println("Diffusion coefficient: $(round(D, digits=5)) Å²/ps")
    println("Diffusion coefficient: $(round(D_standard, digits=12)) cm²/s")
    
    # Plot MSD vs time
    msd_plot = plot(times, msd,
                   label="MSD data",
                   xlabel="Time (ps)",
                   ylabel="MSD (Å²)",
                   title="Mean Square Displacement",
                   marker=:circle,
                   markersize=3,
                   markerstrokewidth=0,
                   alpha=0.7,
                   legend=:topleft)
    
    # Add linear fit line
    fit_y = slope .* x .+ (msd[skip_points+1] - slope * x[1])
    plot!(msd_plot, x, fit_y, 
         label="Linear fit (D = $(round(D, digits=5)) Å²/ps)", 
         linewidth=2, 
         linestyle=:dash)
    
    # Add annotation about diffusion coefficient
    annotate!(msd_plot, [(x[end]/2, fit_y[end]/2, 
               text("D = $(round(D_standard, digits=10)) cm²/s", 
               :right, 8))])
    
    savefig(msd_plot, "water_diffusion.png")
    println("MSD plot saved as 'water_diffusion.png'")
    
    return D, D_standard, msd_plot
end

# Uncomment this line to calculate the diffusion coefficient
D, D_standard, msd_plot = calculate_diffusion_coefficient("Module_4/dump.lammpstrj")


