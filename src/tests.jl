# Test suite for MD equilibration validation

using Test

function run_md_validation_tests(data, system, target_temperature)
    @testset "MD Equilibration Tests" begin
        # 1. Temperature stability
        temp_mean = mean(data["temperature"][div(end,2):end])
        temp_std = std(data["temperature"][div(end,2):end])
        @test isapprox(temp_mean, target_temperature, rtol=0.05) # Within 5% of target
        @test temp_std / temp_mean < 0.1  # Fluctuations less than 10%
        
        # 2. Energy conservation
        energy_drift = (data["total_energy"][end] - data["total_energy"][div(end,2)])/data["total_energy"][div(end,2)]
        @test abs(energy_drift) < 0.01  # Energy drift less than 1%
        
        # 3. Energy fluctuations
        energy_mean = mean(data["total_energy"][div(end,2):end])
        energy_std = std(data["total_energy"][div(end,2):end])
        @test energy_std / abs(energy_mean) < 0.05  # Energy fluctuations less than 5%
        
        # 4. Equipartition theorem
        ke_per_dof = mean(data["kinetic_energy"][div(end,2):end]) / (3*system.N/2)
        @test isapprox(ke_per_dof, target_temperature, rtol=0.1)
        
        # 5. Virial theorem (approximate for LJ)
        if hasfield(typeof(system), :virial)
            virial_mean = mean(system.virial[div(end,2):end])
            pe_mean = mean(data["potential_energy"][div(end,2):end])
            # LJ doesn't perfectly follow virial theorem, but should be approximately -3*PE
            @test isapprox(virial_mean, -3*pe_mean, rtol=0.2)
        end
    end
end

# Run the validation tests
run_md_validation_tests(validation_data, system, target_temperature)

#################################### RDF Module Tests ##############################3#######
include("MD_Simulation.jl")
using .MD_Simulation

@testset "RDF Module Tests" begin
    # Test RDF initialization
    @testset "RDF Initialization" begin
        # Create test system
        ϕ = π * 0.8 / 6
        test_system = hard_spheres(N=32, ϕ=ϕ, lattice="fcc")
        forcefield = LennardJones(test_system, Dict(:ε => 1.0, :σ => 1.0, :rcut => 2.5))
        
        # Initialize RDF calculator
        rdf_calc = rdf(bins=50, r_max=0.5 * box_length(test_system))
        rdf_initialize!(rdf_calc, test_system)
        
        # Check initialization state
        @test rdf_calc.bins == 50
        @test rdf_calc.r_max ≈ box_length(test_system)
        @test rdf_calc.dr ≈ box_length(test_system) / 50
        @test rdf_calc.averages == 0
        @test rdf_calc.N == number_of_particles(test_system)
        @test all(rdf_calc.hist .== 0.0)
    end
    
    # Test RDF computation for a simple case
    @testset "RDF Computation" begin
        # Create a simple system with known structure
        ϕ = π * 0.3 / 6  # Dilute system
        test_system = hard_spheres(N=32, ϕ=ϕ, lattice="fcc")
        forcefield = LennardJones(test_system, Dict(:ε => 1.0, :σ => 1.0, :rcut => 2.5))
        
        # Initialize integrator and equilibrate briefly
        integrator = VelocityVerlet(0.001, test_system, forcefield, 1.0, 100)
        vv_initialize!(integrator, test_system, forcefield, 0.001, 1.0, 100)
        vv_integrate!(integrator, 100)  # Short equilibration
        
        # Create and initialize RDF calculator
        rdf_calc = rdf(bins=100, r_max=0.5 * box_length(test_system))
        rdf_initialize!(rdf_calc, test_system)
        
        # Take some RDF samples
        for i in 1:5
            rdf_compute!(rdf_calc, test_system)
            vv_integrate!(integrator, 10)
        end
        
        # Finalize RDF calculation
        rdf_finalize!(rdf_calc)
        
        # Basic tests on the RDF result
        @test length(rdf_calc.hist) == 100
        @test all(rdf_calc.hist .>= 0.0)  # g(r) should be non-negative
    end
    
    # Test RDF for a system with known peaks
    @testset "RDF Peak Structure" begin
        # Create FCC lattice which has known peak structure
        ϕ = π * 0.45 / 6
        test_system = hard_spheres(N=108, ϕ=ϕ, lattice="fcc")  # 108 particles for 3x3x3 FCC cells
        forcefield = LennardJones(test_system, Dict(:ε => 0.1, :σ => 1.0, :rcut => 3.0))
        
        # Create RDF calculator
        rdf_calc = rdf(bins=200, r_max=0.5 * box_length(test_system))
        rdf_initialize!(rdf_calc, test_system)
        
        # Compute RDF for the crystalline structure (without dynamics)
        rdf_compute!(rdf_calc, test_system)
        rdf_finalize!(rdf_calc)
        
        # In an FCC lattice, we expect peaks at approximately r/σ = √2, √3, 2, etc.
        # Find the positions of the peaks
        peak_indices = findall(i -> i > 1 && i < length(rdf_calc.hist) && 
                              rdf_calc.hist[i] > rdf_calc.hist[i-1] && 
                              rdf_calc.hist[i] > rdf_calc.hist[i+1], 
                              1:length(rdf_calc.hist))
        
        # If we have sufficient peaks, check the pattern
        if length(peak_indices) >= 3
            peak_positions = [(i - 0.5) * rdf_calc.dr for i in peak_indices[1:3]]
            # We expect roughly these peak positions, but with some tolerance
            # due to finite bin size and lattice deformation
            @test peak_positions[1] > 0.9 && peak_positions[1] < 1.5  # Around r ≈ 1.0-1.4
            @test peak_positions[2] > 1.3 && peak_positions[2] < 2.0  # Around r ≈ 1.4-2.0
        end
    end
    
    # Test RDF averaging over multiple samples
    @testset "RDF Averaging" begin
        ϕ = π * 0.6 / 6
        test_system = hard_spheres(N=64, ϕ=ϕ, lattice="fcc")
        forcefield = LennardJones(test_system, Dict(:ε => 1.0, :σ => 1.0, :rcut => 2.5))
        integrator = VelocityVerlet(0.001, test_system, forcefield, 1.5, 50)
        vv_initialize!(integrator, test_system, forcefield, 0.001, 1.5, 50)
        
        # Equilibrate
        vv_integrate!(integrator, 200)
        
        # Create RDF calculator and test accumulation
        rdf_calc = rdf(bins=50, r_max=0.4 * box_length(test_system))
        rdf_initialize!(rdf_calc, test_system)
        
        # First sample
        rdf_compute!(rdf_calc, test_system)
        first_sample = copy(rdf_calc.hist)
        @test rdf_calc.averages == 1
        
        # Move system and take second sample
        vv_integrate!(integrator, 50)
        rdf_compute!(rdf_calc, test_system)
        second_sample = copy(rdf_calc.hist)
        @test rdf_calc.averages == 2
        
        # Check that the histogram has been averaged
        @test !all(first_sample .== second_sample)  # Samples should differ due to dynamics
        @test all(rdf_calc.hist .>= 0.0)  # Averages should still be non-negative
    end
    
    # Test integration with MD run
    @testset "RDF with MD Integration" begin
        # Helper function to run MD with RDF calculation
        function run_md_with_rdf!(integrator, rdf_calc, steps, rdf_interval)
            for i in 1:steps
                vv_integrate!(integrator, 1)
                if mod(i, rdf_interval) == 0
                    rdf_compute!(rdf_calc, integrator.system)
                end
            end
        end
        
        ϕ = π * 0.7 / 6
        test_system = hard_spheres(N=128, ϕ=ϕ, lattice="fcc")
        forcefield = LennardJones(test_system, Dict(:ε => 1.0, :σ => 1.0, :rcut => 2.5))
        integrator = VelocityVerlet(0.001, test_system, forcefield, 2.0, 100)
        vv_initialize!(integrator, test_system, forcefield, 0.001, 2.0, 100)
        
        # Equilibrate
        vv_integrate!(integrator, 300)
        
        # Create RDF calculator
        rdf_calc = rdf(bins=100, r_max=0.5 * box_length(test_system))
        rdf_initialize!(rdf_calc, test_system)
        
        # Run MD with RDF sampling
        run_md_with_rdf!(integrator, rdf_calc, 1000, 50)
        rdf_finalize!(rdf_calc)
        
        # Basic checks on final RDF
        @test rdf_calc.averages > 0
        @test all(rdf_calc.hist .>= 0.0)
        
        # For liquid state, we expect g(r) to approach 1.0 at large r
        # Take average of g(r) for the last 20% of bins
        long_range_avg = mean(rdf_calc.hist[Int(round(0.8*length(rdf_calc.hist))):end])
        @test isapprox(long_range_avg, 1.0, atol=0.2)  # Should approach 1 with some tolerance
    end
end

#################################### MSD Module Tests #####################################3
@testset "MSD Module Tests" begin
    # Test MSD initialization
    @testset "MSD Initialization" begin
        # Create test system
        ϕ = π * 0.8 / 6
        test_system = hard_spheres(N=32, ϕ=ϕ, lattice="fcc")
        
        # Initialize MSD calculator
        msd_calc = msd(T=100)
        msd_initialize!(msd_calc, test_system)
        
        # Check initialization state
        @test size(msd_calc.r₀) == size(positions(test_system))
        @test length(msd_calc.r²) == 100
        @test all(msd_calc.r² .== 0.0)
        @test msd_calc.t == 1
        @test msd_calc.averages == 1
    end
    
    # Test MSD computation for a simple case
    @testset "MSD Computation" begin
        # Create a simple system
        ϕ = π * 0.3 / 6  # Dilute system
        test_system = hard_spheres(N=32, ϕ=ϕ, lattice="fcc")
        forcefield = LennardJones(test_system, Dict(:ε => 1.0, :σ => 1.0, :rcut => 2.5))
        
        # Initialize integrator
        integrator = VelocityVerlet(0.001, test_system, forcefield, 1.0, 100)
        vv_initialize!(integrator, test_system, forcefield, 0.001, 1.0, 100)
        
        # Initialize MSD calculator with small T for easier testing
        msd_calc = msd(T=10)
        msd_initialize!(msd_calc, test_system)
        
        # Store initial positions
        initial_positions = deepcopy(positions(test_system))
        
        # Compute MSD at t=0 (should be 0)
        msd_compute!(msd_calc, test_system)
        @test msd_calc.r²[1] ≈ 0.0
        @test msd_calc.t == 2
        
        # Move particles and compute MSD again
        vv_integrate!(integrator, 10)
        msd_compute!(msd_calc, test_system)
        
        # MSD should be positive after movement
        @test msd_calc.r²[2] > 0.0
        @test msd_calc.t == 3
    end
    
    # Test MSD reset and averaging
    @testset "MSD Reset and Averaging" begin
        ϕ = π * 0.5 / 6
        test_system = hard_spheres(N=64, ϕ=ϕ, lattice="fcc")
        forcefield = LennardJones(test_system, Dict(:ε => 1.0, :σ => 1.0, :rcut => 2.5))
        integrator = VelocityVerlet(0.001, test_system, forcefield, 1.0, 100)
        vv_initialize!(integrator, test_system, forcefield, 0.001, 1.0, 100)
        
        # Very small T for testing reset
        msd_calc = msd(T=5)
        msd_initialize!(msd_calc, test_system)
        
        # Fill the MSD array
        for i in 1:5
            msd_compute!(msd_calc, test_system)
            vv_integrate!(integrator, 5)
        end
        
        # Now t should be 6, triggering a reset on next compute
        @test msd_calc.t == 6
        
        # The next compute should reset
        first_r₀ = deepcopy(msd_calc.r₀)
        msd_compute!(msd_calc, test_system)
        
        # Check that r₀ was updated and t was reset
        @test msd_calc.t == 2
        @test msd_calc.averages == 2
        @test !all(msd_calc.r₀ .== first_r₀)
    end
    
    # Test MSD diffusive behavior
    @testset "MSD Diffusive Behavior" begin
        ϕ = π * 0.4 / 6  # Moderately dilute liquid
        test_system = hard_spheres(N=100, ϕ=ϕ, lattice="fcc")
        forcefield = LennardJones(test_system, Dict(:ε => 1.0, :σ => 1.0, :rcut => 2.5))
        
        # High temperature to get more movement
        integrator = VelocityVerlet(0.001, test_system, forcefield, 2.0, 100)
        vv_initialize!(integrator, test_system, forcefield, 0.001, 2.0, 100)
        
        # Equilibrate
        vv_integrate!(integrator, 200)
        
        # Initialize MSD calculator
        msd_calc = msd(T=50)
        msd_initialize!(msd_calc, test_system)
        
        # Run for more steps and compute MSD
        for i in 1:60
            msd_compute!(msd_calc, test_system)
            vv_integrate!(integrator, 5)
        end
        
        # For diffusive behavior, MSD should generally increase with time
        # Test that MSD at later times is larger than at earlier times
        @test msd_calc.r²[30] > msd_calc.r²[10]
        @test msd_calc.r²[45] > msd_calc.r²[25]
    end
    
    # Test MSD file output
    @testset "MSD File Output" begin
        ϕ = π * 0.5 / 6
        test_system = hard_spheres(N=32, ϕ=ϕ, lattice="fcc")
        forcefield = LennardJones(test_system, Dict(:ε => 1.0, :σ => 1.0, :rcut => 2.5))
        integrator = VelocityVerlet(0.001, test_system, forcefield, 1.0, 100)
        vv_initialize!(integrator, test_system, forcefield, 0.001, 1.0, 100)
        
        # Initialize MSD calculator
        msd_calc = msd(T=10)
        msd_initialize!(msd_calc, test_system)
        
        # Generate some data
        for i in 1:10
            msd_compute!(msd_calc, test_system)
            vv_integrate!(integrator, 5)
        end
        
        # Finalize to write file
        msd_finalize!(msd_calc)
        
        # Check if file was created
        @test isfile("msd.dat")
        
        # Read file and check contents
        data = readlines("msd.dat")
        @test length(data) == 10
        @test all(tryparse.(Float64, data) .!== nothing)
        
        # Clean up
        rm("msd.dat")
    end
end
