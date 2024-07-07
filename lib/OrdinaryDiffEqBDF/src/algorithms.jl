"""
QBDF1: Multistep Method

An alias of `QNDF1` with κ=0.
"""
QBDF1(; kwargs...) = QNDF1(; kappa = 0, kwargs...)

"""
QNDF2: Multistep Method
An adaptive order 2 quasi-constant timestep L-stable numerical differentiation function (NDF) method.

See also `QNDF`.
"""

"""
QBDF2: Multistep Method

An alias of `QNDF2` with κ=0.
"""
QBDF2(; kwargs...) = QNDF2(; kappa = 0, kwargs...)

"""
QNDF: Multistep Method
An adaptive order quasi-constant timestep NDF method.
Utilizes Shampine's accuracy-optimal kappa values as defaults (has a keyword argument for a tuple of kappa coefficients).

@article{shampine1997matlab,
title={The matlab ode suite},
author={Shampine, Lawrence F and Reichelt, Mark W},
journal={SIAM journal on scientific computing},
volume={18},
number={1},
pages={1--22},
year={1997},
publisher={SIAM}
}
"""

"""
QBDF: Multistep Method

An alias of `QNDF` with κ=0.
"""
QBDF(; kwargs...) = QNDF(; kappa = tuple(0 // 1, 0 // 1, 0 // 1, 0 // 1, 0 // 1), kwargs...)

"""
FBDF: Fixed leading coefficient BDF

An adaptive order quasi-constant timestep NDF method.
Utilizes Shampine's accuracy-optimal kappa values as defaults (has a keyword argument for a tuple of kappa coefficients).

@article{shampine2002solving,
title={Solving 0= F (t, y (t), y′(t)) in Matlab},
author={Shampine, Lawrence F},
year={2002},
publisher={Walter de Gruyter GmbH \\& Co. KG}
}
"""

struct FBDF{MO, CS, AD, F, F2, P, FDT, ST, CJ, K, T, StepLimiter} <:
    OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
 max_order::Val{MO}
 linsolve::F
 nlsolve::F2
 precs::P
 κ::K
 tol::T
 extrapolant::Symbol
 controller::Symbol
 step_limiter!::StepLimiter
end

function FBDF(; max_order::Val{MO} = Val{5}(), chunk_size = Val{0}(),
    autodiff = Val{true}(), standardtag = Val{true}(), concrete_jac = nothing,
    diff_type = Val{:forward},
    linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(), κ = nothing,
    tol = nothing,
    extrapolant = :linear, controller = :Standard, step_limiter! = trivial_limiter!) where {MO}
FBDF{MO, _unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
    typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
    _unwrap_val(concrete_jac),
    typeof(κ), typeof(tol), typeof(step_limiter!)}(
    max_order, linsolve, nlsolve, precs, κ, tol, extrapolant,
    controller, step_limiter!)
end

TruncatedStacktraces.@truncate_stacktrace FBDF

"""
Uri M. Ascher, Steven J. Ruuth, Brian T. R. Wetton. Implicit-Explicit Methods for Time-
Dependent Partial Differential Equations. 1995 Society for Industrial and Applied Mathematics
Journal on Numerical Analysis, 32(3), pp 797-823, 1995. doi: https://doi.org/10.1137/0732037
"""
struct SBDF{CS, AD, F, F2, P, FDT, ST, CJ, K, T} <:
       OrdinaryDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    κ::K
    tol::T
    extrapolant::Symbol
    order::Int
    ark::Bool
end

function SBDF(order; chunk_size = Val{0}(), autodiff = Val{true}(),
    standardtag = Val{true}(), concrete_jac = nothing, diff_type = Val{:forward},
    linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(), κ = nothing,
    tol = nothing,
    extrapolant = :linear, ark = false)
SBDF{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve), typeof(nlsolve),
    typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac),
    typeof(κ), typeof(tol)}(linsolve,
    nlsolve,
    precs,
    κ,
    tol,
    extrapolant,
    order,
    ark)
end

# All keyword form needed for remake
function SBDF(; chunk_size = Val{0}(), autodiff = Val{true}(), standardtag = Val{true}(),
    concrete_jac = nothing, diff_type = Val{:forward},
    linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(), κ = nothing,
    tol = nothing,
    extrapolant = :linear,
    order, ark = false)
SBDF{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve), typeof(nlsolve),
    typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac),
    typeof(κ), typeof(tol)}(linsolve,
    nlsolve,
    precs,
    κ,
    tol,
    extrapolant,
    order,
    ark)
end

"""
    SBDF2(;kwargs...)

The two-step version of the IMEX multistep methods of

  - Uri M. Ascher, Steven J. Ruuth, Brian T. R. Wetton.
    Implicit-Explicit Methods for Time-Dependent Partial Differential Equations.
    Society for Industrial and Applied Mathematics.
    Journal on Numerical Analysis, 32(3), pp 797-823, 1995.
    doi: [https://doi.org/10.1137/0732037](https://doi.org/10.1137/0732037)

See also `SBDF`.
"""
SBDF2(; kwargs...) = SBDF(2; kwargs...)

"""
    SBDF3(;kwargs...)

The three-step version of the IMEX multistep methods of

  - Uri M. Ascher, Steven J. Ruuth, Brian T. R. Wetton.
    Implicit-Explicit Methods for Time-Dependent Partial Differential Equations.
    Society for Industrial and Applied Mathematics.
    Journal on Numerical Analysis, 32(3), pp 797-823, 1995.
    doi: [https://doi.org/10.1137/0732037](https://doi.org/10.1137/0732037)

See also `SBDF`.
"""
SBDF3(; kwargs...) = SBDF(3; kwargs...)

"""
    SBDF4(;kwargs...)

The four-step version of the IMEX multistep methods of

  - Uri M. Ascher, Steven J. Ruuth, Brian T. R. Wetton.
    Implicit-Explicit Methods for Time-Dependent Partial Differential Equations.
    Society for Industrial and Applied Mathematics.
    Journal on Numerical Analysis, 32(3), pp 797-823, 1995.
    doi: [https://doi.org/10.1137/0732037](https://doi.org/10.1137/0732037)

See also `SBDF`.
"""
SBDF4(; kwargs...) = SBDF(4; kwargs...)