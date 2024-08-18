module OrdinaryDiffEqRosenbrock

import OrdinaryDiffEqCore: alg_order, alg_adaptive_order, isWmethod, isfsal, _unwrap_val,
                           DEFAULT_PRECS, OrdinaryDiffEqRosenbrockAlgorithm, @cache,
                           alg_cache, initialize!, @unpack,
                           calculate_residuals!, OrdinaryDiffEqMutableCache,
                           OrdinaryDiffEqConstantCache, _ode_interpolant, _ode_interpolant!,
                           _vec, _reshape, perform_step!, trivial_limiter!,
                           OrdinaryDiffEqRosenbrockAdaptiveAlgorithm,
                           OrdinaryDiffEqRosenbrockAlgorithm, generic_solver_docstring,
                           namify, initialize!, perform_step!,
                           constvalue, only_diagonal_mass_matrix,
                           calculate_residuals, has_stiff_interpolation, ODEIntegrator,
                           resize_non_user_cache!, _ode_addsteps!, full_cache,
                           DerivativeOrderNotPossibleError
using MuladdMacro, FastBroadcast, RecursiveArrayTools
import MacroTools
using MacroTools: @capture
using DiffEqBase: @def
import LinearSolve
import LinearSolve: UniformScaling
import ForwardDiff
using FiniteDiff
using LinearAlgebra: mul!, diag, diagm, I, Diagonal, norm
import ADTypes: AutoForwardDiff
import OrdinaryDiffEqCore

using OrdinaryDiffEqDifferentiation: TimeDerivativeWrapper, TimeGradientWrapper,
                                     UDerivativeWrapper, UJacobianWrapper,
                                     wrapprecs, calc_tderivative, build_grad_config,
                                     build_jac_config, issuccess_W, jacobian2W!,
                                     resize_jac_config!, resize_grad_config!,
                                     calc_W, calc_rosenbrock_differentiation!, build_J_W,
                                     UJacobianWrapper, dolinsolve

import OrdinaryDiffEqNonlinearSolve # Required for DAE initialization

using Reexport
@reexport using DiffEqBase

import OrdinaryDiffEqCore: alg_autodiff
import OrdinaryDiffEqCore

function rosenbrock_wanner_docstring(description::String,
        name::String;
        references::String = "",
        extra_keyword_description = "",
        extra_keyword_default = "",
        with_step_limiter = false)
    keyword_default = """
    autodiff = Val{true}(),
    concrete_jac = nothing,
    linsolve = nothing,
    precs = DEFAULT_PRECS,
    """ * extra_keyword_default

    keyword_default_description = """
    - `autodiff`: boolean to control if the Jacobian should be computed via AD or not
    - `concrete_jac`: function of the form `jac!(J, u, p, t)`
    - `linsolve`: custom solver for the inner linear systems
    - `precs`: custom preconditioner for the inner linear solver
    """ * extra_keyword_description

    if with_step_limiter
        keyword_default *= "step_limiter! = OrdinaryDiffEq.trivial_limiter!,\n"
        keyword_default_description *= "- `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`\n"
    end

    generic_solver_docstring(
        description, name, "Rosenbrock-Wanner Method. ", references,
        keyword_default_description, keyword_default
    )
end

function rosenbrock_docstring(description::String,
        name::String;
        references::String = "",
        extra_keyword_description = "",
        extra_keyword_default = "",
        with_step_limiter = false)
    keyword_default = """
    autodiff = Val{true}(),
    concrete_jac = nothing,
    linsolve = nothing,
    precs = DEFAULT_PRECS,
    """ * extra_keyword_default

    keyword_default_description = """
    - `autodiff`: boolean to control if the Jacobian should be computed via AD or not
    - `concrete_jac`: function of the form `jac!(J, u, p, t)`
    - `linsolve`: custom solver for the inner linear systems
    - `precs`: custom preconditioner for the inner linear solver
    """ * extra_keyword_description

    if with_step_limiter
        keyword_default *= "step_limiter! = OrdinaryDiffEq.trivial_limiter!,\n"
        keyword_default_description *= "- `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`\n"
    end

    generic_solver_docstring(
        description, name, "Rosenbrock Method. ", references,
        keyword_default_description, keyword_default
    )
end

include("algorithms.jl")
include("alg_utils.jl")
include("generic_rosenbrock.jl")
include("rosenbrock_caches.jl")
include("rosenbrock_tableaus.jl")
include("interp_func.jl")
include("rosenbrock_interpolants.jl")
include("stiff_addsteps.jl")
include("rosenbrock_perform_step.jl")
include("integrator_interface.jl")

import PrecompileTools
import Preferences
PrecompileTools.@compile_workload begin
    lorenz = OrdinaryDiffEqCore.lorenz
    lorenz_oop = OrdinaryDiffEqCore.lorenz_oop
    solver_list = [Rosenbrock23(), Rodas5P()]
    prob_list = []

    if Preferences.@load_preference("PrecompileDefaultSpecialize", true)
        push!(prob_list, ODEProblem(lorenz, [1.0; 0.0; 0.0], (0.0, 1.0)))
        push!(prob_list, ODEProblem(lorenz, [1.0; 0.0; 0.0], (0.0, 1.0), Float64[]))
    end

    if Preferences.@load_preference("PrecompileAutoSpecialize", false)
        push!(prob_list,
            ODEProblem{true, SciMLBase.AutoSpecialize}(lorenz, [1.0; 0.0; 0.0],
                (0.0, 1.0)))
        push!(prob_list,
            ODEProblem{true, SciMLBase.AutoSpecialize}(lorenz, [1.0; 0.0; 0.0],
                (0.0, 1.0), Float64[]))
    end

    if Preferences.@load_preference("PrecompileFunctionWrapperSpecialize", false)
        push!(prob_list,
            ODEProblem{true, SciMLBase.FunctionWrapperSpecialize}(lorenz, [1.0; 0.0; 0.0],
                (0.0, 1.0)))
        push!(prob_list,
            ODEProblem{true, SciMLBase.FunctionWrapperSpecialize}(lorenz, [1.0; 0.0; 0.0],
                (0.0, 1.0), Float64[]))
    end

    if Preferences.@load_preference("PrecompileNoSpecialize", false)
        push!(prob_list,
            ODEProblem{true, SciMLBase.NoSpecialize}(lorenz, [1.0; 0.0; 0.0], (0.0, 1.0)))
        push!(prob_list,
            ODEProblem{true, SciMLBase.NoSpecialize}(lorenz, [1.0; 0.0; 0.0], (0.0, 1.0),
                Float64[]))
    end

    for prob in prob_list, solver in solver_list
        solve(prob, solver)(5.0)
    end

    prob_list = nothing
    solver_list = nothing
end

export Rosenbrock23, Rosenbrock32, RosShamp4, Veldd4, Velds4, GRK4T, GRK4A,
       Ros4LStab, ROS3P, Rodas3, Rodas23W, Rodas3P, Rodas4, Rodas42, Rodas4P, Rodas4P2,
       Rodas5, Rodas5P, Rodas5Pe, Rodas5Pr,
       RosenbrockW6S4OS, ROS34PW1a, ROS34PW1b, ROS34PW2, ROS34PW3, ROS34PRw, ROS3PRL,
       ROS3PRL2, ROK4a,
       ROS2, ROS2PR, ROS2S, ROS3, ROS3PR, Scholz4_7

end
