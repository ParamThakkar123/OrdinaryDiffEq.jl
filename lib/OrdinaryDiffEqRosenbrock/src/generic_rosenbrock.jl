@doc rosenbrock_wanner_docstring(
"""
An Order 3/2 A-Stable Rosenbrock-W method which is good for mildly stiff equations without oscillations at low tolerances. Note that this method is prone to instability in the presence of oscillations, so use with caution. 2nd order stiff-aware interpolation.
""",
"Rosenbrock32",
references = """
- Shampine L.F. and Reichelt M., (1997) The MATLAB ODE Suite, SIAM Journal of
Scientific Computing, 18 (1), pp. 1-22.
""",
with_step_limiter = true) Rosenbrock32

@doc rosenbrock_wanner_docstring(
"""
An Order 2/3 L-Stable Rosenbrock-W method which is good for very stiff equations with oscillations at low tolerances. 2nd order stiff-aware interpolation.
""",
"Rosenbrock23",
references = """
- Shampine L.F. and Reichelt M., (1997) The MATLAB ODE Suite, SIAM Journal of
Scientific Computing, 18 (1), pp. 1-22.
""",
with_step_limiter = true) Rosenbrock23

"""
    RosShamp4Tableau()

L. F. Shampine, Implementation of Rosenbrock Methods,
ACM Transactions on Mathematical Software (TOMS), 8: 2, 93-113.
doi:10.1145/355993.355994
"""
function RosShamp4Tableau()
    a=[0      0     0;
       2      0     0;
       48//25 6//25 0]
    C=[ 0         0        0    0;
       -8         0        0    0;
        372//25   12//5    0    0;
       -112//125 -54//125 -2//5 0]
    b=[19//9,1//2,25//108,125//108]
    btilde=[17//54,7//36,0,125//108]
    gamma=1//2
    c=[0,1,3//5]
    d=[1//2,-3//2,242//100,116//1000]#2.42,0.116
    RosenbrockAdaptiveTableau(a,C,b,btilde,gamma,d,c)
end

@doc rosenbrock_docstring(
"""
An A-stable 4th order Rosenbrock method.
""",
"RosShamp4",
references = """
- L. F. Shampine, Implementation of Rosenbrock Methods, ACM Transactions on
  Mathematical Software (TOMS), 8: 2, 93-113. doi:10.1145/355993.355994
""") RosShamp4

"""
    Veldd4Tableau()

van Veldhuizen, D-stability and Kaps-Rentrop-methods,
M. Computing (1984) 32: 229. doi:10.1007/BF02243574
"""
function Veldd4Tableau()
    a=[0                 0                 0;
       2                 0                 0;
       4.812234362695436 4.578146956747842 0]
    C=[ 0                  0                  0                 0;
       -5.333333333333331  0                  0                 0;
        6.100529678848254  1.804736797378427  0                 0;
       -2.540515456634749 -9.443746328915205 -1.988471753215993 0]
    b=[4.289339254654537,5.036098482851414,0.6085736420673917,1.355958941201148]
    btilde=[2.175672787531755,2.950911222575741,-0.7859744544887430,-1.355958941201148]
    gamma=0.2257081148225682
    c=[0,0.4514162296451364,0.8755928946018455]
    d=[0.2257081148225682,-0.04599403502680582,0.5177590504944076,-0.03805623938054428]
    RosenbrockAdaptiveTableau(a,C,b,btilde,gamma,d,c)
end

@doc rosenbrock_docstring(
"""
A 4th order D-stable Rosenbrock method.
""",
"Veldd4",
references = """
- van Veldhuizen, D-stability and Kaps-Rentrop-methods, M. Computing (1984) 32: 229.
  doi:10.1007/BF02243574
""",
with_step_limiter=true) Veldd4

"""
    Velds4Tableau()

van Veldhuizen, D-stability and Kaps-Rentrop-methods,
M. Computing (1984) 32: 229. doi:10.1007/BF02243574
"""
function Velds4Tableau()
    a=[0    0    0;
       2    0    0;
       7//4 1//4 0]
    C=[ 0     0    0 0;
       -8     0    0 0;
       -8    -1    0 0;
        1//2 -1//2 2 0]
    b=[4//3,2//3,-4//3,4//3]
    btilde=[-1//3,-1//3,0,-4//3]
    gamma=1//2
    c=[0,1,1//2]
    d=[1//2,-3//2,-3//4,1//4]
    RosenbrockAdaptiveTableau(a,C,b,btilde,gamma,d,c)
end

@doc rosenbrock_docstring(
"""
A 4th order A-stable Rosenbrock method.
""",
"Velds4",
references = """
- van Veldhuizen, D-stability and Kaps-Rentrop-methods, M. Computing (1984) 32: 229.
  doi:10.1007/BF02243574
""",
with_step_limiter=true) Velds4

"""
    GRK4TTableau()

Kaps, P. & Rentrop, Generalized Runge-Kutta methods of order four
with stepsize control for stiff ordinary differential equations.
P. Numer. Math. (1979) 33: 55. doi:10.1007/BF01396495
"""
function GRK4TTableau()
    a=[0                 0                 0;
       2                 0                 0;
       4.524708207373116 4.163528788597648 0]
    C=[ 0                  0                   0                 0;
       -5.071675338776316  0                   0                 0;
        6.020152728650786  0.1597506846727117  0                 0;
       -1.856343618686113 -8.505380858179826  -2.084075136023187 0]
    b=[3.957503746640777,4.624892388363313,0.6174772638750108,1.282612945269037]
    btilde=[2.302155402932996,3.073634485392623,-0.8732808018045032,-1.282612945269037]
    gamma=0.231
    c=[0,0.462,0.8802083333333334]
    d=[0.2310000000000000,-0.03962966775244303,0.5507789395789127,-0.05535098457052764]
    RosenbrockAdaptiveTableau(a,C,b,btilde,gamma,d,c)
end

@doc rosenbrock_docstring(
"""
An efficient 4th order Rosenbrock method.
""",
"GRK4T",
references = """
- Kaps, P. & Rentrop, Generalized Runge-Kutta methods of order four with stepsize control
  for stiff ordinary differential equations. P. Numer. Math. (1979) 33: 55. doi:10.1007/BF01396495
""",
with_step_limiter=true) GRK4T

"""
    GRK4ATableau()

Kaps, P. & Rentrop, Generalized Runge-Kutta methods of order four
with stepsize control for stiff ordinary differential equations.
P. Numer. Math. (1979) 33: 55. doi:10.1007/BF01396495
"""
function GRK4ATableau()
    a=[0                 0                  0;
       1.108860759493671 0                  0;
       2.377085261983360 0.1850114988899692 0]
    C=[ 0                 0                  0                 0;
       -4.920188402397641 0                  0                 0;
        1.055588686048583 3.351817267668938  0                 0;
        3.846869007049313 3.427109241268180 -2.162408848753263 0]
    b=[1.845683240405840,0.1369796894360503,0.7129097783291559,0.6329113924050632]
    btilde=[0.04831870177201765,-0.6471108651049505,0.2186876660500240,-0.6329113924050632]
    gamma=0.395
    c=[0,0.438,0.87]
    d=[0.395,-0.3726723954840920,0.06629196544571492,0.4340946962568634]
    RosenbrockAdaptiveTableau(a,C,b,btilde,gamma,d,c)
end

@doc rosenbrock_docstring(
"""
An A-stable 4th order Rosenbrock method. Essentially "anti-L-stable" but efficient.
""",
"",
references = """
- Kaps, P. & Rentrop, Generalized Runge-Kutta methods of order four with stepsize control
  for stiff ordinary differential equations. P. Numer. Math. (1979) 33: 55. doi:10.1007/BF01396495
""",
with_step_limiter=true) GRK4A

"""
    Ros4LSTableau()

E. Hairer, G. Wanner, Solving ordinary differential equations II,
stiff and differential-algebraic problems. Computational mathematics (2nd revised ed.), Springer (1996)
"""
function Ros4LSTableau()
    a=[0                 0                  0;
       2                 0                  0;
       1.867943637803922 0.2344449711399156 0]
    C=[ 0                  0                   0                  0;
       -7.137615036412310  0                   0                  0;
        2.580708087951457  0.6515950076447975  0                  0;
       -2.137148994382534 -0.3214669691237626 -0.6949742501781779 0]
    b=[2.255570073418735,0.2870493262186792,0.4353179431840180,1.093502252409163]
    btilde=[-0.2815431932141155,-0.07276199124938920,-0.1082196201495311,-1.093502252409163]
    gamma=0.5728200000000000
    c=[0,1.145640000000000,0.6552168638155900]
    d=[0.5728200000000000,-1.769193891319233,0.7592633437920482,-0.1049021087100450]
    RosenbrockAdaptiveTableau(a,C,b,btilde,gamma,d,c)
end

@doc rosenbrock_docstring(
"""
A 4th order A-stable stiffly stable Rosenbrock method with a stiff-aware 3rd order interpolant
""",
"Ros4LStab",
references = """
- E. Hairer, G. Wanner, Solving ordinary differential equations II, stiff and
  differential-algebraic problems. Computational mathematics (2nd revised ed.), Springer (1996)
""",
with_step_limiter=true) Ros4LStab

@doc rosenbrock_docstring(
"""
3rd order A-stable and stiffly stable Rosenbrock method. Keeps high accuracy on discretizations of nonlinear parabolic PDEs.
""",
"ROS3P",
references = """
- Lang, J. & Verwer, ROS3P—An Accurate Third-Order Rosenbrock Solver Designed for
  Parabolic Problems J. BIT Numerical Mathematics (2001) 41: 731. doi:10.1023/A:1021900219772
""",
with_step_limiter = true) ROS3P

"""
    ROS3PRTableau()

3nd order stiffly accurate Rosenbrock-Wanner method with 3 internal stages with B_PR consistent of order 3, which is strongly A-stable with Rinf~=-0.73.

Rang, Joachim (2014): The Prothero and Robinson example: 
Convergence studies for Runge-Kutta and Rosenbrock-Wanner methods. https://doi.org/10.24355/dbbs.084-201408121139-0
"""
function ROS3PRTableau() # 3rd order
    gamma=7.88675134594813e-01
    Alpha=[0                       0                      0;
           2.36602540378444e+00    0                      0;
           0.00000000000000e+00    1.0000000000000e+00    0]
    Gamma=[gamma                  0                        0;
           -2.36602540378444e+00   gamma                   0;
           -2.84686425165674e-01   -1.08133897861876e+00   gamma]
    B=[2.92663844023951e-01,-8.13389786187641e-02, 7.88675134594813e-01]
    Bhat=[1.11324865405187e-01, 1.00000000000000e-01, 7.88675134594813e-01]
    a,C,b,btilde,d,c=_transformtab(Alpha,Gamma,B,Bhat)
    RosenbrockAdaptiveTableau(a,C,b,btilde,gamma,d,c)
end

@doc rosenbrock_docstring(
"""
3nd order stiffly accurate Rosenbrock method with 3 internal stages with B_PR consistent of order 3, which is strongly A-stable with Rinf~=-0.73.
""",
"ROS3PR",
references = """
- Rang, Joachim (2014): The Prothero and Robinson example:
  Convergence studies for Runge-Kutta and Rosenbrock-Wanner methods.
  https://doi.org/10.24355/dbbs.084-201408121139-0
""") ROS3PR

"""
    ROS3PRLTableau()

3rd order stiffly accurate Rosenbrock-Wanner method with 4 internal stages,
B_PR consistent of order 2 with Rinf=0.
The order of convergence decreases if medium stiff problems are considered, but it has good results for very stiff cases.

Rang, Joachim (2014): The Prothero and Robinson example: 
Convergence studies for Runge-Kutta and Rosenbrock-Wanner methods. https://doi.org/10.24355/dbbs.084-201408121139-0
"""
function ROS3PRLTableau() # 3rd order
    gamma=4.3586652150845900e-01
    Alpha=[0                       0                       0                       0;
           5.00000000000000e-01    0                       0                       0;
           5.00000000000000e-01    5.00000000000000e-01    0                       0;
           5.00000000000000e-01    5.00000000000000e-01    0                       0]
    Gamma=[ gamma                  0                      0                        0;
           -5.00000000000000e-01   gamma                  0                        0;
           -7.91564804204642e-01   3.52442167927514e-01   gamma                    0;
           -4.97889699145187e-01   3.86075154415805e-01   -3.24051976779077e-01    gamma]
    B=[2.11030085481324e-03, 8.86075154415805e-01, -3.24051976779077e-01, 4.35866521508459e-01]
    Bhat=[0.5, 3.87524229532982e-01, -2.09492263150452e-01, 3.21968033617470e-01]
    a,C,b,btilde,d,c=_transformtab(Alpha,Gamma,B,Bhat)
    RosenbrockAdaptiveTableau(a,C,b,btilde,gamma,d,c)
end

@doc rosenbrock_docstring(
"""
3rd order stiffly accurate Rosenbrock method with 4 internal stages,
B_PR consistent of order 2 with Rinf=0.
The order of convergence decreases if medium stiff problems are considered, but it has good results for very stiff cases.
""",
"ROS3PRL",
references = """
- Rang, Joachim (2014): The Prothero and Robinson example:
  Convergence studies for Runge-Kutta and Rosenbrock-Wanner methods.
  https://doi.org/10.24355/dbbs.084-201408121139-0
""") ROS3PRL

"""
    ROS3PRL2Tableau()

3rd order stiffly accurate Rosenbrock method with 4 internal stages,
B_PR consistent of order 3.
The order of convergence does NOT decreases if medium stiff problems are considered as it does for [`ROS3PRL`](@ref).

Rang, Joachim (2014): The Prothero and Robinson example: 
Convergence studies for Runge-Kutta and Rosenbrock-Wanner methods. https://doi.org/10.24355/dbbs.084-201408121139-0
"""
function ROS3PRL2Tableau() # 3rd order
    gamma=4.35866521508459e-01
    Alpha=[0                       0                       0                       0;
           1.30759956452538e+00    0                       0                       0;
           5.00000000000000e-01    5.00000000000000e-01    0                       0;
           5.00000000000000e-01    5.00000000000000e-01    0                       0]
    Gamma=[gamma                  0                       0                       0;
           -1.30759956452538e+00   gamma                  0                        0;
           -7.09885758609722e-01   -5.59967359602778e-01  gamma                    0;
           -1.55508568075521e-01   -9.53885165751122e-01  6.73527212318184e-01    gamma]
    B=[3.44491431924479e-01,-4.53885165751122e-01, 6.73527212318184e-01, 4.35866521508459e-01]
    Bhat=[5.00000000000000e-01, -2.57388120865221e-01, 4.35420087247750e-01, 3.21968033617470e-01]
    a,C,b,btilde,d,c=_transformtab(Alpha,Gamma,B,Bhat)
    RosenbrockAdaptiveTableau(a,C,b,btilde,gamma,d,c)
end

@doc rosenbrock_docstring(
"""
3rd order stiffly accurate Rosenbrock method with 4 internal stages,
B_PR consistent of order 3.
The order of convergence does NOT decreases if medium stiff problems are considered as it does for [`ROS3PRL`](@ref).
""",
"ROS3PRL2",
references = """
- Rang, Joachim (2014): The Prothero and Robinson example:
  Convergence studies for Runge-Kutta and Rosenbrock-Wanner methods.
  https://doi.org/10.24355/dbbs.084-201408121139-0
""") ROS3PRL2

@doc rosenbrock_docstring(
"""
3rd order A-stable and stiffly stable Rosenbrock method.
""",
"Rodas3",
references = """
- Steinebach G. Construction of Rosenbrock–Wanner method Rodas5P and numerical benchmarks 
  within the Julia Differential Equations package.
  In: BIT Numerical Mathematics, 63(2), 2023
""",
with_step_limiter=true) Rodas3

@doc rosenbrock_docstring(
"""
3rd order A-stable and stiffly stable Rosenbrock method with a stiff-aware 3rd order interpolant
and additional error test for interpolation. Keeps accuracy on discretizations of linear parabolic PDEs.
""",
"Rodas3P",
references = """
- Steinebach G., Rodas23W / Rodas32P - a Rosenbrock-type method for DAEs with additional error estimate
  for dense output and Julia implementation,
  In progress.
""",
with_step_limiter=true) Rodas3P

@doc rosenbrock_wanner_docstring(
"""
An Order 2/3 L-Stable Rosenbrock-W method for stiff ODEs and DAEs in mass matrix form. 2nd order stiff-aware interpolation and additional error test for interpolation.
""",
"Rodas23W",
references = """
- Steinebach G., Rodas23W / Rodas32P - a Rosenbrock-type method for DAEs with additional error estimate for dense output and Julia implementation,
  In progress.
""",
with_step_limiter = true) Rodas23W

@doc rosenbrock_wanner_docstring(
"""
4th order A-stable stiffly stable Rosenbrock method with a stiff-aware 3rd order interpolant. 4th order
on linear parabolic problems and 3rd order accurate on nonlinear parabolic problems (as opposed to
lower if not corrected).
""",
"Rodas4P",
references = """
- Steinebach G., Rodas23W / Rodas32P - a Rosenbrock-type method for DAEs with additional error estimate 
  for dense output and Julia implementation,
  In progress.
""",
with_step_limiter=true) Rodas4P

@doc rosenbrock_wanner_docstring(
"""
A 4th order L-stable stiffly stable Rosenbrock method with a stiff-aware 3rd order interpolant. 4th order
on linear parabolic problems and 3rd order accurate on nonlinear parabolic problems. It is an improvement
of Roadas4P and in case of inexact Jacobians a second order W method.
""",
"Rodas4P2",
references = """
- Steinebach G., Rodas23W / Rodas32P - a Rosenbrock-type method for DAEs with additional error estimate 
  for dense output and Julia implementation,
  In progress.
""",
with_step_limiter=true) Rodas4P2

@doc rosenbrock_docstring(
"""
A 4th order L-stable Rosenbrock method.
""",
"Rodas4",
references = """
- E. Hairer, G. Wanner, Solving ordinary differential equations II, stiff and
  differential-algebraic problems. Computational mathematics (2nd revised ed.), Springer (1996)
""",
with_step_limiter=true) Rodas4

@doc rosenbrock_docstring(
"""
A 4th order A-stable stiffly stable Rosenbrock method with a stiff-aware 3rd order interpolant
""",
"Rodas42",
references = """
- E. Hairer, G. Wanner, Solving ordinary differential equations II, stiff and
  differential-algebraic problems. Computational mathematics (2nd revised ed.), Springer (1996)
""",
with_step_limiter=true) Rodas42

@doc rosenbrock_docstring(
"""
A 5th order A-stable stiffly stable Rosenbrock method with a stiff-aware 4th order interpolant.
""",
"Rodas5",
references = """
- Di Marzo G. RODAS5(4) – Méthodes de Rosenbrock d’ordre 5(4) adaptées aux problemes
  différentiels-algébriques. MSc mathematics thesis, Faculty of Science,
  University of Geneva, Switzerland.
""",
with_step_limiter=true) Rodas5

@doc rosenbrock_docstring(
"""
A 5th order A-stable stiffly stable Rosenbrock method with a stiff-aware 4th order interpolant.
Has improved stability in the adaptive time stepping embedding.
""",
"Rodas5P",
references = """
- Steinebach G. Construction of Rosenbrock–Wanner method Rodas5P and numerical benchmarks 
  within the Julia Differential Equations package.
  In: BIT Numerical Mathematics, 63(2), 2023
""",
with_step_limiter=true) Rodas5P

@doc rosenbrock_docstring(
"""
A 5th order A-stable stiffly stable Rosenbrock method with a stiff-aware 4th order interpolant.
Has improved stability in the adaptive time stepping embedding.
""",
"Rodas5Pr",
references = """
- Steinebach G. Rosenbrock methods within OrdinaryDiffEq.jl - Overview, recent developments and applications -
  Preprint 2024
  https://github.com/hbrs-cse/RosenbrockMethods/blob/main/paper/JuliaPaper.pdf
""",
with_step_limiter=true) Rodas5Pr

@doc rosenbrock_docstring(
"""
A 5th order A-stable stiffly stable Rosenbrock method with a stiff-aware 4th order interpolant.
Has improved stability in the adaptive time stepping embedding.
""",
"Rodas5Pe",
references = """
- Steinebach G. Rosenbrock methods within OrdinaryDiffEq.jl - Overview, recent developments and applications -
  Preprint 2024
  https://github.com/hbrs-cse/RosenbrockMethods/blob/main/paper/JuliaPaper.pdf
""",
with_step_limiter=true) Rodas5Pe

#! format: off

abstract type RosenbrockTableau{T,T2} end
struct RosenbrockFixedTableau{T,T2}<:RosenbrockTableau{T,T2}
    a::Array{T,2}
    C::Array{T,2}
    b::Array{T,1}
    gamma::T2
    d::Array{T,1}
    c::Array{T2,1}
end

struct RosenbrockAdaptiveTableau{T,T2}<:RosenbrockTableau{T,T2}
    a::Array{T,2}
    C::Array{T,2}
    b::Array{T,1}
    btilde::Array{T,1}
    gamma::T2
    d::Array{T,1}
    c::Array{T2,1}
end

"""
    @_bitarray2boolarray RosenbrockTableau(tab.a.!=0,...)

Transform BitArray (in the form of `xs.!=0` ) into 1D-Array of Bools by
`[i for i in xs.!=0]` to satisfy the type constraint of RosenbrockTableau
"""
macro _bitarray2boolarray(expr)
    args=[:([i for i in $arg]) for arg in expr.args[2:end]]
    args[end-2]=:(tab.gamma!=0)
    esc(:($(expr.args[1])($(args...))))
end

"""
    _masktab(tab)

Convert normal tableau into a dummy tableau consisting of Bools. We use dummy tableaus
where we only care about whether values in the tableau are zeros.
"""
_masktab(tab::RosenbrockFixedTableau)=@_bitarray2boolarray RosenbrockFixedTableau(tab.a.!=0,tab.C.!=0,tab.b.!=0,tab.gamma!=0,tab.d.!=0,tab.c.!=0)
_masktab(tab::RosenbrockAdaptiveTableau)=@_bitarray2boolarray RosenbrockAdaptiveTableau(tab.a.!=0,tab.C.!=0,tab.b.!=0,tab.btilde.!=0,tab.gamma!=0,tab.d.!=0,tab.c.!=0)


"""
    _common_nonzero_vals(tab::RosenbrockTableau)

Return the common nonzero symbols in the tableau. Typical return value:
`[[:a21,:a31,:a32],[:C21,:C31,:C32],[:b1,:b2,:b3],:gamma,[:d1,:d2,:d3],[:c1,:c2,:c3]]`
"""
function _common_nonzero_vals(tab::RosenbrockTableau)
    nzvals=[]
    push!(nzvals,[Symbol(:a,ind[1],ind[2]) for ind in findall(!iszero,tab.a)])
    push!(nzvals,[Symbol(:C,ind[1],ind[2]) for ind in findall(!iszero,tab.C)])
    push!(nzvals,[Symbol(:b,ind) for ind in findall(!iszero,tab.b)])
    push!(nzvals,:gamma)
    push!(nzvals,[Symbol(:d,ind) for ind in findall(!iszero,tab.d)])
    push!(nzvals,[Symbol(:c,ind) for ind in findall(!iszero,tab.c)])
    nzvals
end

"""
    _nonzero_vals(tab::RosenbrockFixedTableau)

Return all the nonzero symbols in the tableau. Typical return value:
`[:a21,:a31,:a32,:C21,:C31,:C32,:b1,:b2,:b3,:gamma,:d1,:d2,:d3,:c1,:c2,:c3]`
"""
function _nonzero_vals(tab::RosenbrockFixedTableau)
    nzvals=_common_nonzero_vals(tab)
    vcat(nzvals...)
end

"""
    _nonzero_vals(tab::RosenbrockAdaptiveTableau)

Typical return value:
`[:a21,:a31,:a32,:C21,:C31,:C32,:b1,:b2,:b3,:btilde1,:btilde2,:btilde3,:gamma,:d1,:d2,:d3,:c1,:c2,:c3]`
"""
function _nonzero_vals(tab::RosenbrockAdaptiveTableau)
    nzvals=_common_nonzero_vals(tab)
    push!(nzvals,[Symbol(:btilde,ind) for ind in findall(!iszero,tab.btilde)])
    vcat(nzvals...)
end

"""
    _push_assigns!(valexprs,inds,name,type)

Insert a series of field statements like `[:(c2::T2),:(c3::T2)]` into the array `valexprs`.

# Arguments
- `valexpr::Array{Expr,1}`: the array to be inserted
- `inds`: an iterator that gives indices
- `name::Symbol`: the prefix name of the values
- `type::Symbol`: type in the statements
"""
function _push_assigns!(valexprs,inds,name,type::Symbol)
    for ind in inds
        push!(valexprs,:($(Symbol(name,"$(Tuple(ind)...)"))::$type))
    end
end

"""
    gen_tableau_struct(tab::RosenbrockTableau,tabstructname::Symbol)

Generate the tableau struct expression from a given tableau emulating those in
`tableaus/rosenbrock_tableaus.jl`. The statements of `aij`,`Cij` and `ci` are generated
according to the nonzero values of `tab.a`,`tab.C` and `tab.c` while others are generated
from their indices. One may choose to pass in a dummy tableau with type `<:RosenbrockTalbeau{Bool,Bool}`
to fully control the tableau struct.
"""
function gen_tableau_struct(tab::RosenbrockTableau,tabstructname::Symbol)
    valexprs=Array{Expr,1}()
    _push_assigns!(valexprs,findall(!iszero,tab.a),:a,:T)
    _push_assigns!(valexprs,findall(!iszero,tab.C),:C,:T)
    _push_assigns!(valexprs,eachindex(tab.b),:b,:T)
    if typeof(tab)<:RosenbrockAdaptiveTableau
        _push_assigns!(valexprs,eachindex(tab.btilde),:btilde,:T)
    end
    push!(valexprs,:(gamma::T2))
    _push_assigns!(valexprs,eachindex(tab.d),:d,:T)
    _push_assigns!(valexprs,findall(!iszero,tab.c),:c,:T2)
    quote struct $tabstructname{T,T2}
        $(valexprs...)
        end
    end
end

"""
    gen_tableau(tab::RosenbrockTableau,tabstructexpr::Expr,tabname::Symbol)

Generate the tableau function expression emulating those in `tableaus/rosenbrock_tableaus.jl`.
It takes in the tableau struct expression (generated by gen_tableau_struct(...) or written by hand)
to make sure the actual values of the tableau are organized in the right order.
"""
function gen_tableau(tab::RosenbrockTableau,tabstructexpr::Expr,tabname::Symbol)
    @capture(tabstructexpr, struct T_ fields__ end) || error("incorrect tableau expression")
    tabstructname = namify(T)
    valsym2tabdict=Dict("a"=>tab.a,"C"=>tab.C,"gamma"=>tab.gamma,"c"=>tab.c,"d"=>tab.d,"b"=>tab.b)
    if typeof(tab)<:RosenbrockAdaptiveTableau
        valsym2tabdict["btilde"]=tab.btilde
    end
    pattern=r"^([a-zA-Z]+)([1-9]{0,2})$"
    assignexprs = Expr[]
    valsyms = Symbol[]
    for field in fields
        if @capture(field, valsym_Symbol::valtype_)
            push!(valsyms, valsym)
            m = match(pattern, String(valsym))
            val = valsym2tabdict[m[1]][(parse(Int, i) for i in m[2])...]
            push!(assignexprs, :($valsym = convert($valtype, $val)))
        end
    end
    quote function $tabname(T, T2)
            $(assignexprs...)
            $tabstructname($(valsyms...))
        end
    end
end

"""
    gen_cache_struct(tab::RosenbrockTableau,cachename::Symbol,constcachename::Symbol)

Generate cache struct expression emulating those in `caches/rosenbrock_caches.jl`.
The length of k1,k2,... in the mutable cache struct is determined by the length of `tab.b`
because in the end of Rosenbrock's method, we have: `y_{n+1}=y_n+ki*bi`.
"""
function gen_cache_struct(tab::RosenbrockTableau,cachename::Symbol,constcachename::Symbol)
    kstype=[:($(Symbol(:k,i))::rateType) for i in 1:length(tab.b)]
    constcacheexpr=quote struct $constcachename{TF,UF,Tab,JType,WType,F} <: OrdinaryDiffEqConstantCache
            tf::TF
            uf::UF
            tab::Tab
            J::JType
            W::WType
            linsolve::F
        end
    end
    cacheexpr=quote
        @cache mutable struct $cachename{uType,rateType,uNoUnitsType,JType,WType,TabType,TFType,UFType,F,JCType,GCType} <: RosenbrockMutableCache
            u::uType
            uprev::uType
            du::rateType
            du1::rateType
            du2::rateType
            $(kstype...)
            fsalfirst::rateType
            fsallast::rateType
            dT::rateType
            J::JType
            W::WType
            tmp::rateType
            atmp::uNoUnitsType
            weight::uNoUnitsType
            tab::TabType
            tf::TFType
            uf::UFType
            linsolve_tmp::rateType
            linsolve::F
            jac_config::JCType
            grad_config::GCType
        end
    end
    constcacheexpr,cacheexpr
end

"""
    gen_algcache(cacheexpr::Expr,constcachename::Symbol,algname::Symbol,tabname::Symbol)

Generate expressions for `alg_cache(...)` emulating those in `caches/rosenbrock_caches.jl`.
"""
function gen_algcache(cacheexpr::Expr,constcachename::Symbol,algname::Symbol,tabname::Symbol)
    @capture(cacheexpr, @cache mutable struct T_ fields__ end) || error("incorrect cache expression")
    cachename = namify(T)
    ksinit = Expr[]
    valsyms = Symbol[]
    for field in fields
        if @capture(field, valsym_Symbol::valtype_)
            push!(valsyms, valsym)

            if match(r"^k[1-9]+$", String(valsym)) !== nothing
                push!(ksinit, :($valsym = zero(rate_prototype)))
            end
        end
    end

    quote
        function alg_cache(alg::$algname,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
            tf = TimeDerivativeWrapper(f,u,p)
            uf = UDerivativeWrapper(f,t,p)
            J,W = build_J_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits,Val(false))
            $constcachename(tf,uf,$tabname(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits)),J,W,nothing)
        end
        function alg_cache(alg::$algname,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
            du = zero(rate_prototype)
            du1 = zero(rate_prototype)
            du2 = zero(rate_prototype)
            $(ksinit...)
            fsalfirst = zero(rate_prototype)
            fsallast = zero(rate_prototype)
            dT = zero(rate_prototype)
            J,W = build_J_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits,Val(true))
            tmp = zero(rate_prototype)
            atmp = similar(u, uEltypeNoUnits)
            weight = similar(u, uEltypeNoUnits)
            tab = $tabname(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))

            tf = TimeGradientWrapper(f,uprev,p)
            uf = UJacobianWrapper(f,t,p)
            linsolve_tmp = zero(rate_prototype)
            linprob = LinearProblem(W,_vec(linsolve_tmp); u0=_vec(tmp))
            linsolve = init(linprob,alg.linsolve,alias_A=true,alias_b=true,
                            Pl = LinearSolve.InvPreconditioner(Diagonal(_vec(weight))),
                            Pr = Diagonal(_vec(weight)))
            grad_config = build_grad_config(alg,f,tf,du1,t)
            jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,du2)
            $cachename($(valsyms...))
        end
    end
end

"""
    gen_initialize(cachename::Symbol,constcachename::Symbol)

Generate expressions for `initialize!(...)` in `perform_step/rosenbrock_perform_step.jl`.
It only generates a default version of `initialize!` which support 3rd-order Hermite interpolation.
"""
function gen_initialize(cachename::Symbol,constcachename::Symbol)
    quote
        function initialize!(integrator, cache::$constcachename)
            integrator.kshortsize = 2
            integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
            integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
            integrator.stats.nf += 1

            # Avoid undefined entries if k is an array of arrays
            integrator.fsallast = zero(integrator.fsalfirst)
            integrator.k[1] = integrator.fsalfirst
            integrator.k[2] = integrator.fsallast
          end

          function initialize!(integrator, cache::$cachename)
            integrator.kshortsize = 2
            @unpack fsalfirst,fsallast = cache
            integrator.fsalfirst = fsalfirst
            integrator.fsallast = fsallast
            resize!(integrator.k, integrator.kshortsize)
            integrator.k .= [fsalfirst,fsallast]
            integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
            integrator.stats.nf += 1
          end
    end
end

"""
    gen_constant_perform_step(tabmask::RosenbrockTableau{Bool,Bool},cachename::Symbol,n_normalstep::Int,specialstepexpr=:nothing)

Generate non-inplace version of `perform_step!` expression emulating those in `perform_step/rosenbrock_perform_step.jl`.
The `perform_step!` function calculates `k1,k2,k3,...` defined by `(-W)ki=f(u+aij*kj,t+ci*dt)+di*dt*dT+Cij*kj*dt, i=1,2,...,n_normalstep`
and then gives the result by `y_{n+1}=y_n+ki*bi`. Terms with 0s (according to tabmask) are skipped in the expressions.
Special steps can be added before calculating `y_{n+1}`. The non-inplace `perform_step!` assumes the mass_matrix==I.
"""
function gen_constant_perform_step(tabmask::RosenbrockTableau{Bool,Bool},cachename::Symbol,n_normalstep::Int,specialstepexpr=:nothing)
    unpacktabexpr=:(@unpack ()=cache.tab)
    unpacktabexpr.args[3].args[1].args=_nonzero_vals(tabmask)
    dtCijexprs=[:($(Symbol(:dtC,Cind[1],Cind[2]))=$(Symbol(:C,Cind[1],Cind[2]))/dt) for Cind in findall(!iszero,tabmask.C)]
    dtdiexprs=[:($(Symbol(:dtd,dind))=dt*$(Symbol(:d,dind))) for dind in eachindex(tabmask.d)]
    iterexprs=[]
    for i in 1:n_normalstep
        aijkj=[:($(Symbol(:a,i+1,j))*$(Symbol(:k,j))) for j in findall(!iszero,tabmask.a[i+1,:])]
        Cijkj=[:($(Symbol(:dtC,i+1,j))*$(Symbol(:k,j))) for j in findall(!iszero,tabmask.C[i+1,:])]
        push!(iterexprs,
        quote
            $(Symbol(:k,i)) = _reshape(W\-_vec(linsolve_tmp), axes(uprev))
            integrator.stats.nsolve += 1
            u=+(uprev,$(aijkj...))
            du = f(u, p, t+$(Symbol(:c,i+1))*dt)
            integrator.stats.nf += 1
            if mass_matrix === I
                linsolve_tmp=+(du,$(Symbol(:dtd,i+1))*dT,$(Cijkj...))
            else
                linsolve_tmp=du+$(Symbol(:dtd,i+1))*dT+mass_matrix*(+($(Cijkj...)))
            end
        end)
    end
    push!(iterexprs,specialstepexpr)
    n=length(tabmask.b)
    biki=[:($(Symbol(:b,i))*$(Symbol(:k,i))) for i in 1:n]
    push!(iterexprs,
    quote
        $(Symbol(:k,n))=_reshape(W\-_vec(linsolve_tmp), axes(uprev))
        integrator.stats.nsolve += 1
        u=+(uprev,$(biki...))
    end)

    adaptiveexpr=[]
    if typeof(tabmask)<:RosenbrockAdaptiveTableau
        btildeiki=[:($(Symbol(:btilde,i))*$(Symbol(:k,i))) for i in findall(!iszero,tabmask.btilde)]
        push!(adaptiveexpr,quote
            if integrator.opts.adaptive
                utilde =  +($(btildeiki...))
                atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
                                        integrator.opts.reltol,integrator.opts.internalnorm,t)
                integrator.EEst = integrator.opts.internalnorm(atmp,t)
            end
        end)
    end
    quote
        @muladd function perform_step!(integrator, cache::$cachename, repeat_step=false)
            @unpack t,dt,uprev,u,f,p = integrator
            @unpack tf,uf = cache
            $unpacktabexpr

            $(dtCijexprs...)
            $(dtdiexprs...)
            dtgamma = dt*gamma

            mass_matrix = integrator.f.mass_matrix

            # Time derivative
            tf.u = uprev
            dT = ForwardDiff.derivative(tf, t)

            W = calc_W(integrator, cache, dtgamma, repeat_step, true)
            linsolve_tmp = integrator.fsalfirst + dtd1*dT #calc_rosenbrock_differentiation!

            $(iterexprs...)

            integrator.fsallast = f(u, p, t + dt)
            integrator.stats.nf += 1

            integrator.k[1] = integrator.fsalfirst
            integrator.k[2] = integrator.fsallast
            integrator.u = u

            $(adaptiveexpr...)
        end
    end

end

"""
    gen_perform_step(tabmask::RosenbrockTableau{Bool,Bool},cachename::Symbol,n_normalstep::Int,specialstepexpr=:nothing)

Generate inplace version of `perform_step!` expression emulating those in `perform_step/rosenbrock_perform_step.jl`.
The inplace `perform_step!` produces the same result as the non-inplace version except that it treats the mass_matrix appropriately.
"""
function gen_perform_step(tabmask::RosenbrockTableau{Bool,Bool},cachename::Symbol,n_normalstep::Int,specialstepexpr=:nothing)
    unpacktabexpr=:(@unpack ()=cache.tab)
    unpacktabexpr.args[3].args[1].args=_nonzero_vals(tabmask)
    dtCij=[:($(Symbol(:dtC,"$(Cind[1])$(Cind[2])"))=$(Symbol(:C,"$(Cind[1])$(Cind[2])"))/dt) for Cind in findall(!iszero,tabmask.C)]
    dtdi=[:($(Symbol(:dtd,dind[1]))=dt*$(Symbol(:d,dind[1]))) for dind in eachindex(tabmask.d)]
    iterexprs=[]
    for i in 1:n_normalstep
        ki=Symbol(:k,i)
        dtdj=Symbol(:dtd,i+1)
        aijkj=[:($(Symbol(:a,i+1,j))*$(Symbol(:k,j))) for j in findall(!iszero,tabmask.a[i+1,:])]
        dtCijkj=[:($(Symbol(:dtC,i+1,j))*$(Symbol(:k,j))) for j in findall(!iszero,tabmask.C[i+1,:])]
        repeatstepexpr=[]
        if i==1
            repeatstepexpr=[:(!repeat_step)]
        end
        push!(iterexprs,quote

            if $(i==1)
                # Must be a part of the first linsolve for preconditioner step
                linres = dolinsolve(integrator, linsolve; A = !repeat_step ? W : nothing, b = _vec(linsolve_tmp))
            else
                linres = dolinsolve(integrator, linsolve; b = _vec(linsolve_tmp))
            end

            linsolve = linres.cache
            vecu = _vec(linres.u)
            vecki = _vec($ki)

            @.. broadcast=false vecki = -vecu

            integrator.stats.nsolve += 1
            @.. broadcast=false u = +(uprev,$(aijkj...))
            f( du,  u, p, t+$(Symbol(:c,i+1))*dt)
            integrator.stats.nf += 1
            if mass_matrix === I
                @.. broadcast=false linsolve_tmp = +(du,$dtdj*dT,$(dtCijkj...))
            else
                @.. broadcast=false du1 = +($(dtCijkj...))
                mul!(_vec(du2),mass_matrix,_vec(du1))
                @.. broadcast=false linsolve_tmp = du + $dtdj*dT + du2
            end
        end)
    end
    push!(iterexprs,specialstepexpr)
    n=length(tabmask.b)
    ks=[Symbol(:k,i) for i in 1:n]
    klast=Symbol(:k,n)
    biki=[:($(Symbol(:b,i))*$(Symbol(:k,i))) for i in 1:n]
    push!(iterexprs,quote

        linres = dolinsolve(integrator, linsolve; b = _vec(linsolve_tmp))
        linsolve = linres.cache
        vecu = _vec(linres.u)
        vecklast = _vec($klast)
        @.. broadcast=false vecklast = -vecu

        integrator.stats.nsolve += 1
        @.. broadcast=false u = +(uprev,$(biki...))
    end)

    adaptiveexpr=[]
    if typeof(tabmask)<:RosenbrockAdaptiveTableau
        btildeiki=[:($(Symbol(:btilde,i))*$(Symbol(:k,i))) for i in findall(!iszero,tabmask.btilde)]
        push!(adaptiveexpr,quote
            utilde=du
            if integrator.opts.adaptive
                @.. broadcast=false utilde = +($(btildeiki...))
                calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol,
                                    integrator.opts.reltol,integrator.opts.internalnorm,t)
                integrator.EEst = integrator.opts.internalnorm(atmp,t)
            end
        end)
    end
    quote
        @muladd function perform_step!(integrator, cache::$cachename, repeat_step=false)
            @unpack t,dt,uprev,u,f,p = integrator
            @unpack du,du1,du2,fsallast,dT,J,W,uf,tf,$(ks...),linsolve_tmp,jac_config,atmp,weight = cache
            $unpacktabexpr

            # Assignments
            sizeu  = size(u)
            uidx = eachindex(integrator.uprev)
            mass_matrix = integrator.f.mass_matrix

            # Precalculations
            $(dtCij...)
            $(dtdi...)
            dtgamma = dt*gamma

            calculate_residuals!(weight, fill!(weight, one(eltype(u))), uprev, uprev,
                                 integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)

            calc_rosenbrock_differentiation!(integrator, cache, dtd1, dtgamma, repeat_step, true)

            linsolve = cache.linsolve

            $(iterexprs...)

            f( fsallast,  u, p, t + dt)
            integrator.stats.nf += 1

            $(adaptiveexpr...)
        end
    end
end

"""
    RosenbrockW6S4OSTableau()

Rahunanthan, A., & Stanescu, D. (2010). High-order W-methods.
Journal of computational and applied mathematics, 233(8), 1798-1811.
"""
function RosenbrockW6S4OSTableau()
    a=[0                  0                  0                  0                  0;
       0.5812383407115008 0                  0                  0                  0;
       0.9039624413714670 1.8615191555345010 0                  0                  0;
       2.0765797196750000 0.1884255381414796 1.8701589674910320 0                  0;
       4.4355506384843120 5.4571817986101890 4.6163507880689300 3.1181119524023610 0;
       10.791701698483260 -10.05691522584131 14.995644854284190 5.2743399543909430 1.4297308712611900]
    C=[0                  0                  0                  0                  0;
       -2.661294105131369 0                  0                  0                  0;
       -3.128450202373838 0.0000000000000000 0                  0                  0;
       -6.920335474535658 -1.202675288266817 -9.733561811413620 0                  0;
       -28.09530629102695 20.371262954793770 -41.04375275302869 -19.66373175620895 0;
       9.7998186780974000 11.935792886603180 3.6738749290132010 14.807828541095500 0.8318583998690680]
    b=[6.4562170746532350,-4.853141317768053,9.7653183340692600,2.0810841772787230,0.6603936866352417,0.6000000000000000]
    gamma=0.2500000000000000
    d=[0.2500000000000000,0.0836691184292894,0.0544718623516351,-0.3402289722355864,0.0337651588339529,-0.0903074267618540]
    c=[0                 ,0.1453095851778752,0.3817422770256738,0.6367813704374599,0.7560744496323561,0.9271047239875670]
    RosenbrockFixedTableau(a,C,b,gamma,d,c)
end

"""
    @RosenbrockW6S4OS(part)

Generate code for the RosenbrockW6S4OS method.
`part` should be one of `:tableau`, `:cache`, `:init`, `:performstep`.
`@RosenbrockW6S4OS(:tableau)` should be placed in `tableaus/rosenbrock_tableaus.jl`.
`@RosenbrockW6S4OS(:cache)` should be placed in `caches/rosenbrock_caches.jl`.
`@RosenbrockW6S4OS(:init)` and `@RosenbrockW6S4OS(:performstep)` should be
placed in `perform_step/rosenbrock_perform_step.jl`.
"""
macro RosenbrockW6S4OS(part)
    tab=RosenbrockW6S4OSTableau()
    tabmask=_masktab(tab)
    algname=:RosenbrockW6S4OS
    tabname=:RosenbrockW6S4OSTableau
    tabstructname=:RosenbrockW6STableau
    cachename=:RosenbrockW6SCache
    constcachename=:RosenbrockW6SConstantCache
    n_normalstep=length(tab.b)-1
    if part.value==:tableau
        #println("Generating const cache")
        tabstructexpr=gen_tableau_struct(tabmask,tabstructname)
        tabexpr=gen_tableau(tab,tabstructexpr,tabname)
        return esc(quote $([tabstructexpr,tabexpr]...) end)
    elseif part.value==:cache
        #println("Generating cache")
        constcacheexpr,cacheexpr=gen_cache_struct(tabmask,cachename,constcachename)
        algcacheexpr=gen_algcache(cacheexpr,constcachename,algname,tabname)
        return esc(quote $([constcacheexpr,cacheexpr,algcacheexpr]...) end)
    elseif part.value==:init
        #println("Generating initialize")
        return esc(gen_initialize(cachename,constcachename))
    elseif part.value==:performstep
        #println("Generating perform_step")
        constperformstepexpr=gen_constant_perform_step(tabmask,constcachename,n_normalstep)
        performstepexpr=gen_perform_step(tabmask,cachename,n_normalstep)
        return esc(quote $([constperformstepexpr,performstepexpr]...) end)
    else
        throw(ArgumentError("Unknown parameter!"))
        nothing
    end
end

"""
    Ros4dummyTableau()

Generate a dummy tableau for ROS4 methods. It can be considered as performing elementwise OR to the masks
of those specific tableaus: `Ros4dummyTableau()==_masktab(RosShamp4Tableau()) OR _masktab(Veldd4Tableau()) OR ...`
ROS4 methods have the property of a4j==a3j so a is a 3*3 matrix instead of a 4*4 matrix and c is a 1*3 vector instead of a 1*4 vector.
"""
function Ros4dummyTableau()#create a tabmask for all ROS4 methods where false->0,true->non-0
    a=[false false false;
       true  false false;
       true  true  false]
    C=[false false false false;
       true  false false false;
       true  true  false false;
       true  true  true  false]
    b=[true,true,true,true]
    btilde=[true,true,true,true]
    gamma=true
    c=[false,true,true]
    d=[true,true,true,true]
    RosenbrockAdaptiveTableau(a,C,b,btilde,gamma,d,c)
end

"""
    @Rosenbrock4(part)

Generate code for the Rosenbrock4 methods: RosShamp4, Veldd4, Velds4, GRK4A, GRK4T, Ros4LStab.
`part` should be one of `:tableau`, `:cache`, `:performstep`.
`@Rosenbrock4(:tableau)` should be placed in `tableaus/rosenbrock_tableaus.jl`.
`@Rosenbrock4(:cache)` should be placed in `caches/rosenbrock_caches.jl`.
`@Rosenbrock4(:performstep)` should be placed in `perform_step/rosenbrock_perform_step.jl`.
The `initialize!` function for Rosenbrock4 methods is already included in `rosenbrock_perform_step.jl`.
The special property of ROS4 methods that a4j==a3j requires a special step in `perform_step!` that
calculates `linsolve_tmp` from the previous `du` which reduces a function call.
"""
macro Rosenbrock4(part)
    tabmask=Ros4dummyTableau()#_masktab(tab)
    cachename=:Rosenbrock4Cache
    constcachename=:Rosenbrock4ConstantCache
    RosShamp4tabname=:RosShamp4Tableau
    Veldd4tabname=:Veldd4Tableau
    Velds4tabname=:Velds4Tableau
    GRK4Ttabname=:GRK4TTableau
    GRK4Atabname=:GRK4ATableau
    Ros4LStabname=:Ros4LStabTableau
    n_normalstep=2 #for the third step a4j=a3j which reduced one function call
    if part.value==:tableau
        #println("Generating tableau for Rosenbrock4")
        tabstructexpr=gen_tableau_struct(tabmask,:Ros4Tableau)
        tabexprs=Array{Expr,1}()
        push!(tabexprs,tabstructexpr)
        push!(tabexprs,gen_tableau(RosShamp4Tableau(),tabstructexpr,RosShamp4tabname))
        push!(tabexprs,gen_tableau(Veldd4Tableau(),tabstructexpr,Veldd4tabname))
        push!(tabexprs,gen_tableau(Velds4Tableau(),tabstructexpr,Velds4tabname))
        push!(tabexprs,gen_tableau(GRK4TTableau(),tabstructexpr,GRK4Ttabname))
        push!(tabexprs,gen_tableau(GRK4ATableau(),tabstructexpr,GRK4Atabname))
        push!(tabexprs,gen_tableau(Ros4LSTableau(),tabstructexpr,Ros4LStabname))
        return esc(quote $(tabexprs...) end)
    elseif part.value==:cache
        #println("Generating cache for Rosenbrock4")
        constcacheexpr,cacheexpr=gen_cache_struct(tabmask,cachename,constcachename)
        cacheexprs=Array{Expr,1}([constcacheexpr,cacheexpr])
        push!(cacheexprs,gen_algcache(cacheexpr,constcachename,:RosShamp4,RosShamp4tabname))
        push!(cacheexprs,gen_algcache(cacheexpr,constcachename,:Veldd4,Veldd4tabname))
        push!(cacheexprs,gen_algcache(cacheexpr,constcachename,:Velds4,Velds4tabname))
        push!(cacheexprs,gen_algcache(cacheexpr,constcachename,:GRK4T,GRK4Ttabname))
        push!(cacheexprs,gen_algcache(cacheexpr,constcachename,:GRK4A,GRK4Atabname))
        push!(cacheexprs,gen_algcache(cacheexpr,constcachename,:Ros4LStab,Ros4LStabname))
        return esc(quote $(cacheexprs...) end)
    elseif part.value==:performstep
        #println("Generating perform_step for Rosenbrock4")
        specialstepconst=quote
            k3 = _reshape(W\-_vec(linsolve_tmp), axes(uprev))
            integrator.stats.nsolve += 1
            #u = uprev  + a31*k1 + a32*k2 #a4j=a3j
            #du = f(u, p, t+c3*dt) #reduced function call
            if mass_matrix === I
                linsolve_tmp =  du + dtd4*dT + dtC41*k1 + dtC42*k2 + dtC43*k3
            else
                linsolve_tmp = du + dtd4*dT + mass_matrix * (dtC41*k1 + dtC42*k2 + dtC43*k3)
            end
        end
        specialstep=quote

            linres = dolinsolve(integrator, linsolve; b = _vec(linsolve_tmp))
            linsolve = linres.cache
            cache.linsolve = linsolve
            vecu = _vec(linres.u)
            veck3 = _vec(k3)
            @.. broadcast=false veck3 = -vecu

            integrator.stats.nsolve += 1
            #@.. broadcast=false u = uprev + a31*k1 + a32*k2 #a4j=a3j
            #f( du,  u, p, t+c3*dt) #reduced function call
            if mass_matrix === I
                @.. broadcast=false linsolve_tmp = du + dtd4*dT + dtC41*k1 + dtC42*k2 + dtC43*k3
            else
                @.. broadcast=false du1 = dtC41*k1 + dtC42*k2 + dtC43*k3
                mul!(du2,mass_matrix,du1)
                @.. broadcast=false linsolve_tmp = du + dtd4*dT + du2
            end
        end
        constperformstepexpr=gen_constant_perform_step(tabmask,constcachename,n_normalstep,specialstepconst)
        performstepexpr=gen_perform_step(tabmask,cachename,n_normalstep,specialstep)
        return esc(quote $([constperformstepexpr,performstepexpr]...) end)
    else
        throw(ArgumentError("Unknown parameter!"))
        nothing
    end
end

#ROS2, ROS23 and ROS34PW methods (Rang and Angermann, 2005)

"""
    Ros34dummyTableau()

Generate a dummy tableau for ROS34W methods proposed by Rang and Angermann. This type of methods has 4 steps.
"""
function Ros34dummyTableau()
    a=[false false false false;
       true  false false false;
       true  true  false false;
       true  true  true  false]
    C=[false false false false;
       true  false false false;
       true  true  false false;
       true  true  true  false]
    b=[true,true,true,true]
    btilde=[true,true,true,true]
    gamma=true
    c=[false,true,true,true]
    d=[true,true,true,true]
    RosenbrockAdaptiveTableau(a,C,b,btilde,gamma,d,c)
end

"""
    Ros23dummyTableau()

Generate a dummy tableau for ROS23 methods proposed by Rang. This type of methods has 3 steps.
"""
function Ros23dummyTableau()
    a=[false false false;
       true  false false;
       true  true  false]
    C=[false false false;
       true  false false;
       true  true  false;]
    b=[true,true,true]
    btilde=[true,true,true]
    gamma=true
    c=[false,true,true]
    d=[true,true,true]
    RosenbrockAdaptiveTableau(a,C,b,btilde,gamma,d,c)
end

"""
    Ros2dummyTableau()

Generate a dummy tableau for ROS2 methods. This type of methods has 2 steps.
"""
function Ros2dummyTableau()
    a=[false false;
       true  false]
    C=[false false;
       true  false]
    b=[true,true]
    btilde=[true,true]
    gamma=true
    c=[false,true]
    d=[true,true]
    RosenbrockAdaptiveTableau(a,C,b,btilde,gamma,d,c)
end


"""
    _transformtab(Alpha,Gamma,B,Bhat)

Transform the tableau from values in the paper into values used in OrdinaryDiffEq according to p112 in Hairer and Wanner.

E. Hairer, G. Wanner, Solving ordinary differential equations II, stiff and
differential-algebraic problems. Computational mathematics (2nd revised ed.), Springer (1996)
"""
function _transformtab(Alpha,Gamma,B,Bhat)
    invGamma=inv(Gamma)
    a=Alpha*invGamma
    C=diagm(0=>diag(invGamma))-invGamma
    b=[(transpose(B)*invGamma)...]# [2Darray...]=>1Darray
    btilde=[(transpose(B-Bhat)*invGamma)...]
    gamma=Gamma[1,1]#Gamma11==Gamma22==...==Gammass
    d=[sum(Gamma,dims=2)...]#di=sum_j Gamma_ij
    c=[sum(Alpha,dims=2)...]#ci=sum_j Alpha_ij
    (a,C,b,btilde,d,c)
end




# 2 step ROS Methods
"""
    ROS2Tableau()
2nd order stiffly accurate Rosenbrock method with 2 internal stages with (Rinf=0).
The embedded method is taken from Kinetic PreProcessor (KPP).
J. G. Verwer et al. (1999): A second-order Rosenbrock method applied to photochemical dispersion problems
https://doi.org/10.1137/S1064827597326651
"""
function ROS2Tableau() # 2nd order
    gamma=1.7071067811865475 # 1+1/sqrt(2)
    Alpha=[0     0;
           1.    0]
    Gamma=[gamma                   0;
           -3.414213562373095   gamma]
    B=[0.5, 0.5]
    Bhat=[1, 0]
    a,C,b,btilde,d,c=_transformtab(Alpha,Gamma,B,Bhat)
    RosenbrockAdaptiveTableau(a,C,b,btilde,gamma,d,c)
end

@doc rosenbrock_wanner_docstring(
"""
A 4th order L-stable Rosenbrock-W method.
""",
"ROS34PW1a") ROS34PW1a

@doc rosenbrock_wanner_docstring(
"""
A 4th order L-stable Rosenbrock-W method.
""",
"ROS34PW1b") ROS34PW1b

@doc rosenbrock_wanner_docstring(
"""
A 4th order stiffy accurate Rosenbrock-W method for PDAEs.
""",
"ROS34PW2") ROS34PW2

@doc rosenbrock_wanner_docstring(
"""
A 4th order strongly A-stable (Rinf~0.63) Rosenbrock-W method.
""",
"ROS34PW3") ROS34PW3

"""
    @ROS2(part)

Generate code for the 2 step ROS methods: ROS2
`part` should be one of `:tableau`, `:cache`, `:init`, `:performstep`.
`@ROS2(:tableau)` should be placed in `tableaus/rosenbrock_tableaus.jl`.
`@ROS2(:cache)` should be placed in `caches/rosenbrock_caches.jl`.
`@ROS2(:init)` and `@ROS2(:performstep)` should be placed in
`perform_step/rosenbrock_perform_step.jl`.
"""
macro ROS2(part)
    tabmask=Ros2dummyTableau()
    cachename=:ROS2Cache
    constcachename=:ROS2ConstantCache
    ROS2tabname=:ROS2Tableau
    n_normalstep=length(tabmask.b)-1
    if part.value==:tableau
        tabstructexpr=gen_tableau_struct(tabmask,:Ros2Tableau)
        tabexprs=Array{Expr,1}([tabstructexpr])
        push!(tabexprs,gen_tableau(ROS2Tableau(),tabstructexpr,ROS2tabname))
        return esc(quote $(tabexprs...) end)
    elseif part.value==:cache
        constcacheexpr,cacheexpr=gen_cache_struct(tabmask,cachename,constcachename)
        cacheexprs=Array{Expr,1}([constcacheexpr,cacheexpr])
        push!(cacheexprs,gen_algcache(cacheexpr,constcachename,:ROS2,ROS2tabname))
        return esc(quote $(cacheexprs...) end)
    elseif part.value==:init
        return esc(gen_initialize(cachename,constcachename))
    elseif part.value==:performstep
        performstepexprs=Array{Expr,1}()
        push!(performstepexprs,gen_constant_perform_step(tabmask,constcachename,n_normalstep))
        push!(performstepexprs,gen_perform_step(tabmask,cachename,n_normalstep))
        return esc(quote $(performstepexprs...) end)
    else
        throw(ArgumentError("Unknown parameter!"))
        nothing
    end
end

@doc rosenbrock_docstring(
"""
A 2nd order L-stable Rosenbrock method with 2 internal stages.
""",
"ROS2",
references = """
- J. G. Verwer et al. (1999): A second-order Rosenbrock method applied to photochemical dispersion problems
  https://doi.org/10.1137/S1064827597326651
""") ROS2

# 3 step ROS Methods
"""
    ROS2PRTableau()

2nd order stiffly accurate Rosenbrock method with 3 internal stages with (Rinf=0).
For problems with medium stiffness the convergence behaviour is very poor and it is recommended to use 
[`ROS2S`](@ref) instead.

Rang, Joachim (2014): The Prothero and Robinson example: 
Convergence studies for Runge-Kutta and Rosenbrock-Wanner methods. https://doi.org/10.24355/dbbs.084-201408121139-0
"""
function ROS2PRTableau() # 2nd order
    gamma=2.28155493653962e-01
    Alpha=[0                       0                      0;
           1.00000000000000e+00    0                      0;
           0.00000000000000e+00    1.0000000000000e+00    0]
    Gamma=[gamma                  0                        0;
           -2.28155493653962e-01   gamma                   0;
            6.47798871261042e-01   -8.75954364915004e-01   gamma]
    B=[6.47798871261042e-01,1.24045635084996e-01, 2.28155493653962e-01]
    Bhat=[7.71844506346038e-01, 2.28155493653962e-01, 0.00000000000000e+00]
    a,C,b,btilde,d,c=_transformtab(Alpha,Gamma,B,Bhat)
    RosenbrockAdaptiveTableau(a,C,b,btilde,gamma,d,c)
end

@doc rosenbrock_docstring(
"""
2nd order stiffly accurate Rosenbrock method with 3 internal stages with (Rinf=0).
For problems with medium stiffness the convergence behaviour is very poor and it is recommended to use 
[`ROS2S`](@ref) instead.
""",
"ROS2PR",
references = """
- Rang, Joachim (2014): The Prothero and Robinson example:
  Convergence studies for Runge-Kutta and Rosenbrock-Wanner methods.
  https://doi.org/10.24355/dbbs.084-201408121139-0
""")
ROS2PR



"""
    ROS2STableau()

2nd order stiffly accurate Rosenbrock-Wanner W-method with 3 internal stages with B_PR consistent of order 2 with (Rinf=0).
More Information at https://doi.org/10.24355/dbbs.084-201408121139-0

Rang, Joachim (2014): The Prothero and Robinson example: 
Convergence studies for Runge-Kutta and Rosenbrock-Wanner methods. https://doi.org/10.24355/dbbs.084-201408121139-0
"""
function ROS2STableau() # 2nd order
    gamma=2.92893218813452e-01
    Alpha=[0                       0                      0;
           5.85786437626905e-01    0                      0;
           0.00000000000000e+00    1.0000000000000e+00    0]
    Gamma=[gamma                  0                        0;
           -5.85786437626905e-01   gamma                   0;
            3.53553390593274e-01   -6.46446609406726e-01   gamma]
    B=[3.53553390593274e-01,3.53553390593274e-01, 2.92893218813452e-01]
    Bhat=[3.33333333333333e-01, 3.33333333333333e-01, 3.33333333333333e-01]
    a,C,b,btilde,d,c=_transformtab(Alpha,Gamma,B,Bhat)
    RosenbrockAdaptiveTableau(a,C,b,btilde,gamma,d,c)
end

@doc rosenbrock_wanner_docstring(
"""
2nd order stiffly accurate Rosenbrock-Wanner W-method with 3 internal stages with B_PR consistent of order 2 with (Rinf=0).
""",
"ROS2S",
references = """
- Rang, Joachim (2014): The Prothero and Robinson example:
  Convergence studies for Runge-Kutta and Rosenbrock-Wanner methods.
  https://doi.org/10.24355/dbbs.084-201408121139-0
""")
ROS2S


"""
    ROS3Tableau()
E. Hairer, G. Wanner, Solving ordinary differential equations II,
stiff and differential-algebraic problems. Computational mathematics (2nd revised ed.), Springer (1996)
With coefficients from https://doi.org/10.1016/S1352-2310(97)83212-8

"""
function ROS3Tableau() # 3rd order
    gamma=0.435866521508459
    Alpha=[0                    0                      0;
           0.435866521508459    0                      0;
           0.435866521508459    0                      0]
    Gamma=[gamma                  0                        0;
           -0.19294655696029095   gamma                   0;
           0                      1.7492714812579468   gamma]
    B=[-0.7545741238540432,1.9410040706196443, -0.18642994676560104]
    Bhat=[-1.5335874578414959, 2.817451311486258, -0.28386385364476185]
    a,C,b,btilde,d,c=_transformtab(Alpha,Gamma,B,Bhat)
    RosenbrockAdaptiveTableau(a,C,b,btilde,gamma,d,c)
end

@doc rosenbrock_docstring(
"""
3rd order L-stable Rosenbrock method with 3 internal stages with an embedded strongly
A-stable 2nd order method.
""",
"ROS3",
references = """
- E. Hairer, G. Wanner, Solving ordinary differential equations II, stiff and
  differential-algebraic problems. Computational mathematics (2nd revised ed.), Springer (1996)
""") ROS3


"""
    Scholz4_7Tableau()

3nd order stiffly accurate Rosenbrock method with 3 internal stages with B_PR consistent of order 3, which is strongly A-stable with Rinf~=-0.73
Convergence with order 4 for the stiff case, but has a poor accuracy.

Rang, Joachim (2014): The Prothero and Robinson example: 
Convergence studies for Runge-Kutta and Rosenbrock-Wanner methods. https://doi.org/10.24355/dbbs.084-201408121139-0
"""
function Scholz4_7Tableau() # 3rd order
    gamma=7.88675134594813e-01
    Alpha=[0                       0                      0;
           2.36602540378444e+00    0                      0;
           2.50000000000000e-01    1.0000000000000e+00    0]
    Gamma=[gamma                  0                        0;
           -2.36602540378444e+00   gamma                   0;
           -6.13414364537605e-01   -1.10383267558217e+00   gamma]
    B=[4.95076910424059e-01,-1.12898126628685e-01, 6.17821216204626e-01]
    Bhat=[3.33333333333333e-01, 3.33333333333333e-01, 3.33333333333333e-01]
    a,C,b,btilde,d,c=_transformtab(Alpha,Gamma,B,Bhat)
    RosenbrockAdaptiveTableau(a,C,b,btilde,gamma,d,c)
end

@doc rosenbrock_docstring(
"""
3nd order stiffly accurate Rosenbrock method with 3 internal stages with B_PR consistent of order 3, which is strongly A-stable with Rinf~=-0.73.
Convergence with order 4 for the stiff case, but has a poor accuracy.
""",
"Scholz4_7",
references = """
- Rang, Joachim (2014): The Prothero and Robinson example:
  Convergence studies for Runge-Kutta and Rosenbrock-Wanner methods.
  https://doi.org/10.24355/dbbs.084-201408121139-0
""") Scholz4_7


"""
    @ROS23(part)

Generate code for the 3 step ROS methods: ROS2PR, ROS2S, ROS3, ROS3PR, Scholz4_7
`part` should be one of `:tableau`, `:cache`, `:init`, `:performstep`.
`@ROS23(:tableau)` should be placed in `tableaus/rosenbrock_tableaus.jl`.
`@ROS23(:cache)` should be placed in `caches/rosenbrock_caches.jl`.
`@ROS23(:init)` and `@ROS23(:performstep)` should be placed in
`perform_step/rosenbrock_perform_step.jl`.
"""
macro ROS23(part)
    tabmask=Ros23dummyTableau()
    cachename=:ROS23Cache
    constcachename=:ROS23ConstantCache
    ROS2PRtabname=:ROS2PRTableau
    ROS2Stabname=:ROS2STableau
    ROS3tabname=:ROS3Tableau
    ROS3PRtabname=:ROS3PRTableau
    Scholz4_7tabname=:Scholz4_7Tableau
    n_normalstep=length(tabmask.b)-1
    if part.value==:tableau
        tabstructexpr=gen_tableau_struct(tabmask,:Ros23Tableau)
        tabexprs=Array{Expr,1}([tabstructexpr])
        push!(tabexprs,gen_tableau(ROS2PRTableau(),tabstructexpr,ROS2PRtabname))
        push!(tabexprs,gen_tableau(ROS2STableau(),tabstructexpr,ROS2Stabname))
        push!(tabexprs,gen_tableau(ROS3Tableau(),tabstructexpr,ROS3tabname))
        push!(tabexprs,gen_tableau(ROS3PRTableau(),tabstructexpr,ROS3PRtabname))
        push!(tabexprs,gen_tableau(Scholz4_7Tableau(),tabstructexpr,Scholz4_7tabname))
        return esc(quote $(tabexprs...) end)
    elseif part.value==:cache
        constcacheexpr,cacheexpr=gen_cache_struct(tabmask,cachename,constcachename)
        cacheexprs=Array{Expr,1}([constcacheexpr,cacheexpr])
        push!(cacheexprs,gen_algcache(cacheexpr,constcachename,:ROS2PR,ROS2PRtabname))
        push!(cacheexprs,gen_algcache(cacheexpr,constcachename,:ROS2S,ROS2Stabname))
        push!(cacheexprs,gen_algcache(cacheexpr,constcachename,:ROS3,ROS3tabname))
        push!(cacheexprs,gen_algcache(cacheexpr,constcachename,:ROS3PR,ROS3PRtabname))
        push!(cacheexprs,gen_algcache(cacheexpr,constcachename,:Scholz4_7,Scholz4_7tabname))
        return esc(quote $(cacheexprs...) end)
    elseif part.value==:init
        return esc(gen_initialize(cachename,constcachename))
    elseif part.value==:performstep
        performstepexprs=Array{Expr,1}()
        push!(performstepexprs,gen_constant_perform_step(tabmask,constcachename,n_normalstep))
        push!(performstepexprs,gen_perform_step(tabmask,cachename,n_normalstep))
        return esc(quote $(performstepexprs...) end)
    else
        throw(ArgumentError("Unknown parameter!"))
        nothing
    end
end




# 4 step ROS Methods
"""
    ROS34PW1aTableau()

L-Stable, Rosenbrock-W method with order 3 and 4 inner steps

Rang, J., & Angermann, L. (2005). New Rosenbrock W-methods of order 3 for partial
differential algebraic equations of index 1. BIT Numerical Mathematics, 45(4), 761-787.
"""
function ROS34PW1aTableau()
    gamma=4.358665215084590e-1
    Alpha=[0                 0                    0   0;
           2.218787467653286 0                    0   0;
           0                 0                    0   0; # can reduce one function call with specialized perform_step
           1.208587690772214 7.511610241919324e-2 0.5 0]
    Gamma=[ gamma                 0                    0     0;
           -2.218787467653286     gamma                0     0;
           -9.461966143940745e-2 -7.913526735718213e-3 gamma 0;
           -1.870323744195384    -9.624340112825115e-2 2.726301276675511e-1 gamma]
    B=[3.285609536316354e-1,-5.785609536316354e-1,0.25,1]
    Bhat=[-0.25,0,0.25,1]#B-Bhat[3:4]==[0,0]
    a,C,b,btilde,d,c=_transformtab(Alpha,Gamma,B,Bhat)
    RosenbrockAdaptiveTableau(a,C,b,btilde,gamma,d,c)
end

"""
    ROS34PW1bTableau()

L-Stable, Rosenbrock-W method with order 3 and 4 inner steps

Rang, J., & Angermann, L. (2005). New Rosenbrock W-methods of order 3 for partial
differential algebraic equations of index 1. BIT Numerical Mathematics, 45(4), 761-787.
"""
function ROS34PW1bTableau()
    gamma=4.358665215084590e-1
    Alpha=[0                 0 0   0;
           2.218787467653286 0 0   0;
           2.218787467653286 0 0   0; # can reduce one function call with specialized perform_step
           1.453923375357884 0 0.1 0]
    Gamma=[ gamma              0                    0     0;
           -2.218787467653286  gamma                0     0;
           -2.848610224639349 -5.267530183845237e-2 gamma 0;
           -1.128167857898393 -1.677546870499461e-1 5.452602553351021e-2 gamma]
    B=[5.495647928937977e-1,-5.507258170857301e-1,0.25,7.511610241919324e-1]
    Bhat=[-1.161024191932427e-3,0,0.25,7.511610241919324e-1]#B-Bhat[3:4]==[0,0]
    a,C,b,btilde,d,c=_transformtab(Alpha,Gamma,B,Bhat)
    RosenbrockAdaptiveTableau(a,C,b,btilde,gamma,d,c)
end

"""
    ROS34PW2Tableau()

A stiffy accurate Rosenbrock-W method with order 3 and 4 inner steps whose
embedded method is strongly A-stable with Rinf~=0.48

Rang, J., & Angermann, L. (2005). New Rosenbrock W-methods of order 3 for partial
differential algebraic equations of index 1. BIT Numerical Mathematics, 45(4), 761-787.
"""
function ROS34PW2Tableau()
    gamma=4.3586652150845900e-1
    Alpha=[0                      0                     0 0;
           8.7173304301691801e-1  0                     0 0;
           8.4457060015369423e-1 -1.1299064236484185e-1 0 0;
           0                      0                     1 0]
    Gamma=[ gamma                  0                     0     0;
           -8.7173304301691801e-1  gamma                 0     0;
           -9.0338057013044082e-1  5.4180672388095326e-2 gamma 0;
            2.4212380706095346e-1 -1.2232505839045147    5.4526025533510214e-1 gamma]
    B=[2.4212380706095346e-1,-1.2232505839045147,1.5452602553351020,4.3586652150845900e-1]
    Bhat=[3.7810903145819369e-1,-9.6042292212423178e-2,0.5,2.1793326075422950e-1]
    a,C,b,btilde,d,c=_transformtab(Alpha,Gamma,B,Bhat)
    RosenbrockAdaptiveTableau(a,C,b,btilde,gamma,d,c)
end

"""
    ROS34PW3Tableau()

an A-stable (Rinf~=0.63), Rosenbrock-W method with order 4 and 4 inner steps.

Rang, J., & Angermann, L. (2005). New Rosenbrock W-methods of order 3 for partial
differential algebraic equations of index 1. BIT Numerical Mathematics, 45(4), 761-787.
"""
function ROS34PW3Tableau()#4th order
    gamma=1.0685790213016289
    Alpha=[0                      0                     0 0;
           2.5155456020628817     0                     0 0;
           5.0777280103144085e-1  0.75                  0 0;
           1.3959081404277204e-1 -3.3111001065419338e-1 8.2040559712714178e-1 0]
    Gamma=[ gamma                  0                      0     0;
           -2.5155456020628817     gamma                  0     0;
           -8.7991339217106512e-1 -9.6014187766190695e-1  gamma 0;
           -4.1731389379448741e-1  4.1091047035857703e-1 -1.3558873204765276 gamma]
    B=[2.2047681286931747e-1,2.7828278331185935e-3,7.1844787635140066e-3,7.6955588053404989e-1]
    Bhat=[3.1300297285209688e-1,-2.8946895245112692e-1,9.7646597959903003e-1,0]
    a,C,b,btilde,d,c=_transformtab(Alpha,Gamma,B,Bhat)
    RosenbrockAdaptiveTableau(a,C,b,btilde,gamma,d,c)
end

"""
    ROS34PRwTableau()

3rd order stiffly accurate Rosenbrock-Wanner W-method with 4 internal stages,
B_PR consistent of order 2.
The order of convergence decreases if medium stiff problems are considered.

Joachim Rang, Improved traditional Rosenbrock-Wanner methods for stiff ODEs and DAEs,
Journal of Computational and Applied Mathematics: https://doi.org/10.1016/j.cam.2015.03.010
"""
function ROS34PRwTableau() # 3rd order
    gamma=4.3586652150845900e-01
    Alpha=[0                         0                         0                       0;
           8.7173304301691801e-01    0                         0                       0;
           1.4722022879435914e+00    -3.1840250568090289e-01   0                       0;
           8.1505192016694938e-01    0.5                       -3.1505192016694938e-01 0]
    Gamma=[ gamma                    0                        0                        0;
           -8.7173304301691801e-01   gamma                    0                        0;
           -1.2855347382089872e+00   5.0507005541550687e-01   gamma                    0;
           -4.8201449182864348e-01   2.1793326075422950e-01   -1.7178529043404503e-01 gamma]
    B=[3.3303742833830591e-01, 7.1793326075422947e-01, -4.8683721060099439e-01, 4.3586652150845900e-01]
    Bhat=[0.25, 7.4276119608319180e-01, -3.1472922970066219e-01, 3.2196803361747034e-01]
    a,C,b,btilde,d,c=_transformtab(Alpha,Gamma,B,Bhat)
    RosenbrockAdaptiveTableau(a,C,b,btilde,gamma,d,c)
end

@doc rosenbrock_wanner_docstring(
"""
3rd order stiffly accurate Rosenbrock-Wanner W-method with 4 internal stages,
B_PR consistent of order 2.
The order of convergence decreases if medium stiff problems are considered.
""",
"ROS34PRw",
references = """
- Joachim Rang, Improved traditional Rosenbrock–Wanner methods for stiff ODEs and DAEs,
  Journal of Computational and Applied Mathematics,
  https://doi.org/10.1016/j.cam.2015.03.010
""") ROS34PRw


"""
    ROK4aTableau()

4rd order L-stable Rosenbrock-Krylov method with 4 internal stages,
with a 3rd order embedded method which is strongly A-stable with Rinf~=0.55. (when using exact Jacobians)
Tranquilli, Paul and Sandu, Adrian (2014): Rosenbrock--Krylov Methods for Large Systems of Differential Equations 
https://doi.org/10.1137/130923336
"""
function ROK4aTableau() # 4rd order
    gamma=0.572816062482135
    Alpha=[0                        0                       0                       0;
           1                        0                       0                       0;
           0.10845300169319391758   0.39154699830680608241  0                       0;
           0.43453047756004477624   0.14484349252001492541 -0.07937397008005970166  0]
    Gamma=[gamma                    0                       0                       0;
           -1.91153192976055097824  gamma                   0                       0;
            0.32881824061153522156  0.0                     gamma                   0;
            0.03303644239795811290 -0.24375152376108235312 -0.17062602991994029834  gamma]
    B=   [0.16666666666666666667, 0.16666666666666666667, 0.0, 0.66666666666666666667]
    Bhat=[0.50269322573684235345, 0.27867551969005856226, 0.21863125457309908428, 0.0]
    a,C,b,btilde,d,c=_transformtab(Alpha,Gamma,B,Bhat)
    RosenbrockAdaptiveTableau(a,C,b,btilde,gamma,d,c)
end

@doc rosenbrock_wanner_docstring(
"""
4rd order L-stable Rosenbrock-Krylov method with 4 internal stages,
with a 3rd order embedded method which is strongly A-stable with Rinf~=0.55. (when using exact Jacobians)
""",
"ROK4a",
references = """
- Tranquilli, Paul and Sandu, Adrian (2014): 
  Rosenbrock--Krylov Methods for Large Systems of Differential Equations
  https://doi.org/10.1137/130923336
""") ROK4a

"""
    @ROS34PW(part)

Generate code for the 4 steps ROS34PW methods: ROS34PW1a, ROS34PW1b, ROS34PW2, ROS34PW3, ROS34PRw, ROS3PRL, ROS3PRL2, ROK4a.
`part` should be one of `:tableau`, `:cache`, `:init`, `:performstep`.
`@ROS34PW(:tableau)` should be placed in `tableaus/rosenbrock_tableaus.jl`.
`@ROS34PW(:cache)` should be placed in `caches/rosenbrock_caches.jl`.
`@ROS34PW(:init)` and `@ROS34PW(:performstep)` should be placed in
`perform_step/rosenbrock_perform_step.jl`.
"""
macro ROS34PW(part)
    tabmask=Ros34dummyTableau()
    cachename=:ROS34PWCache
    constcachename=:ROS34PWConstantCache
    ROS34PW1atabname=:ROS34PW1aTableau
    ROS34PW1btabname=:ROS34PW1bTableau
    ROS34PW2tabname=:ROS34PW2Tableau
    ROS34PW3tabname=:ROS34PW3Tableau
    ROS34PRwtabname=:ROS34PRwTableau
    ROS3PRLtabname=:ROS3PRLTableau
    ROS3PRL2tabname=:ROS3PRL2Tableau
    ROK4atabname=:ROK4aTableau
    n_normalstep=length(tabmask.b)-1
    if part.value==:tableau
        tabstructexpr=gen_tableau_struct(tabmask,:Ros34Tableau)
        tabexprs=Array{Expr,1}([tabstructexpr])
        push!(tabexprs,gen_tableau(ROS34PW1aTableau(),tabstructexpr,ROS34PW1atabname))
        push!(tabexprs,gen_tableau(ROS34PW1bTableau(),tabstructexpr,ROS34PW1btabname))
        push!(tabexprs,gen_tableau(ROS34PW2Tableau(),tabstructexpr,ROS34PW2tabname))
        push!(tabexprs,gen_tableau(ROS34PW3Tableau(),tabstructexpr,ROS34PW3tabname))
        push!(tabexprs,gen_tableau(ROS34PRwTableau(),tabstructexpr,ROS34PRwtabname))
        push!(tabexprs,gen_tableau(ROS3PRLTableau(),tabstructexpr,ROS3PRLtabname))
        push!(tabexprs,gen_tableau(ROS3PRL2Tableau(),tabstructexpr,ROS3PRL2tabname))
        push!(tabexprs,gen_tableau(ROK4aTableau(),tabstructexpr,ROK4atabname))
        return esc(quote $(tabexprs...) end)
    elseif part.value==:cache
        constcacheexpr,cacheexpr=gen_cache_struct(tabmask,cachename,constcachename)
        cacheexprs=Array{Expr,1}([constcacheexpr,cacheexpr])
        push!(cacheexprs,gen_algcache(cacheexpr,constcachename,:ROS34PW1a,ROS34PW1atabname))
        push!(cacheexprs,gen_algcache(cacheexpr,constcachename,:ROS34PW1b,ROS34PW1btabname))
        push!(cacheexprs,gen_algcache(cacheexpr,constcachename,:ROS34PW2,ROS34PW2tabname))
        push!(cacheexprs,gen_algcache(cacheexpr,constcachename,:ROS34PW3,ROS34PW3tabname))
        push!(cacheexprs,gen_algcache(cacheexpr,constcachename,:ROS34PRw,ROS34PRwtabname))
        push!(cacheexprs,gen_algcache(cacheexpr,constcachename,:ROS3PRL,ROS3PRLtabname))
        push!(cacheexprs,gen_algcache(cacheexpr,constcachename,:ROS3PRL2,ROS3PRL2tabname))
        push!(cacheexprs,gen_algcache(cacheexpr,constcachename,:ROK4a,ROK4atabname))
        return esc(quote $(cacheexprs...) end)
    elseif part.value==:init
        return esc(gen_initialize(cachename,constcachename))
    elseif part.value==:performstep
        performstepexprs=Array{Expr,1}()
        push!(performstepexprs,gen_constant_perform_step(tabmask,constcachename,n_normalstep))
        push!(performstepexprs,gen_perform_step(tabmask,cachename,n_normalstep))
        return esc(quote $(performstepexprs...) end)
    else
        throw(ArgumentError("Unknown parameter!"))
        nothing
    end
end

#=========================================================================================
# How to add a new method
1. `OrdinaryDiffEq.jl`: export <Algorithm_name>
2. `alg_utils.jl`: alg_order(alg::<Algorithm_name>)=<Algorithm_order>
    if the method is a W-method, add isWmethod(alg::<Algorithm_name>) = true as well
3. `algorithms.jl`: add algorithm struct
4. `generic_rosenbrock.jl`:
    a. write dummy tableau function (or generate from actual tableau using _masktab) for generating
       table struct, cache struct and perform_step
    b. write tableau function. When only `Alpha, Gamma, B, Bhat` are given, use _transformtab
    c. write macro with :tableau, :cache, :init and :performstep
    d. put the macros in the right places.
5. test\algconvergence\ode_rosenbrock_tests.jl: add a test for your method
# How to refactor methods into generic ones
RUN CONVERGENCE TESTS BETWEEN ANY OF THE TWO STEPS!
1. write tableau function and macro definition in this file
2. replace the tableau function (usually named with `XXXConstCache()`) using `gen_tableau()`
   and the original tableau struct expression in `tableaus/rosenbrock_tableaus.jl`
3. replace the `perform_step!` methods in `perform_step/rosenbrock_perform_step.jl` using
   `gen_perform_step()` and `gen_constant_perform_step()`
4. replace cache struct and `alg_cache` in `caches/rosenbrock_caches.jl` using `gen_cache_struct()`
   and `gen_algcache()`
5. If the method only have 3rd-order Hermite interpolation, you can replace `initialize!()`
   in `perform_step/rosenbrock_perform_step.jl` with `gen_initialize()`
DONE

# How to debug
Use macroexpand like `macroexpand(OrdinaryDiffEq,:(@ROS34PW(:performstep)))` and check the
generated codes.

`Revise.jl` is not compatible with macros. One may want to manually re-eval files that use
the macro like `@eval OrdinaryDiffEq include(...)`

# If you want to refactor Rosenbrock methods ...
You need to change respective places in this file.
1. `perform_step/rosenbrock_perform_step.jl` -> `gen_perform_step()`, `gen_constant_perform_step()`,
    `gen_initialize()` and special step expressions in macro definitions
2. `caches/rosenbrock_caches.jl` ->  `gen_algcache()`, `gen_cache_struct()`
3. `tableaus/rosenbrock_tableaus.jl` -> `gen_tableau_struct()` and `gen_tableau()`
=========================================================================================#