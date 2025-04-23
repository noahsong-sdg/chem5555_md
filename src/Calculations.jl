module Calculations

export Calculation, initialize!, compute!, finalize!

abstract type Calculation end

# Define interface functions
function initialize!(calc::Calculation, system)
    error("initialize! not implemented for $(typeof(calc))")
end

""" This is a generic interface for calculations that can be performed on a system.
    It defines the basic structure for initializing, computing, and finalizing calculations.
    
    The `Calculation` type is an abstract type that serves as a base for all specific calculation types.
    The `compute!` function is used to perform the calculation on a given system.
    The `initialize!` function is used to set up the calculation before it starts.
    The `finalize!` function is used to clean up or finalize the calculation after it is done.
    These functions are designed to be overridden by specific calculation types to provide the actual implementation.

    Example usage:
    ```julia
    mutable struct MyCalculation <: Calculation
        # Define specific fields for MyCalculation
    end
    function initialize!(calc::MyCalculation, system)
        # Initialization code for MyCalculation
    end
    function compute!(calc::MyCalculation, system)
        # Computation code for MyCalculation
    end
    function finalize!(calc::MyCalculation)
        # Finalization code for MyCalculation
    end
    ```
    This interface allows for flexibility and extensibility in adding new types of calculations to the system.
    """

function compute!(calc::Calculation, system)
    error("compute! not implemented for $(typeof(calc))")
end

function finalize!(calc::Calculation)
    nothing
end

end