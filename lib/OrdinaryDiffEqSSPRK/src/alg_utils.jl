isfsal(alg::SSPRK53_2N1) = false
isfsal(alg::SSPRK53_2N2) = false
isfsal(alg::SSPRK22) = false
isfsal(alg::SSPRK33) = false
isfsal(alg::SSPRK53) = false
isfsal(alg::SSPRK53_H) = false
isfsal(alg::SSPRK63) = false
isfsal(alg::SSPRK73) = false
isfsal(alg::SSPRK83) = false
isfsal(alg::SSPRK43) = false
isfsal(alg::SSPRK432) = false
isfsal(alg::SSPRK932) = false
isfsal(alg::SSPRK54) = false
isfsal(alg::SSPRK104) = false

alg_order(alg::KYKSSPRK42) = 2
alg_order(alg::SSPRKMSVS32) = 2
alg_order(alg::SSPRK33) = 3
alg_order(alg::SSPRK53_2N1) = 3
alg_order(alg::SSPRK53_2N2) = 3
alg_order(alg::SSPRK22) = 2
alg_order(alg::SSPRK53) = 3
alg_order(alg::SSPRK53_H) = 3
alg_order(alg::SSPRK63) = 3
alg_order(alg::SSPRK73) = 3
alg_order(alg::SSPRK83) = 3
alg_order(alg::SSPRK43) = 3
alg_order(alg::SSPRK432) = 3
alg_order(alg::SSPRKMSVS43) = 3
alg_order(alg::SSPRK932) = 3
alg_order(alg::SSPRK54) = 4
alg_order(alg::SSPRK104) = 4

"""
    ssp_coefficient(alg)

Return the SSP coefficient of the ODE algorithm `alg`. If one time step of size
`dt` with `alg` can be written as a convex combination of explicit Euler steps
with step sizes `cᵢ * dt`, the SSP coefficient is the minimal value of `1/cᵢ`.

# Examples

```julia-repl
julia> ssp_coefficient(SSPRK104())
6
```
"""

ssp_coefficient(alg::SSPRK53_2N1) = 2.18
ssp_coefficient(alg::SSPRK53_2N2) = 2.148
ssp_coefficient(alg::SSPRK53) = 2.65
ssp_coefficient(alg::SSPRK53_H) = 2.65
ssp_coefficient(alg::SSPRK63) = 3.518
ssp_coefficient(alg::SSPRK73) = 4.2879
ssp_coefficient(alg::SSPRK83) = 5.107
ssp_coefficient(alg::SSPRK43) = 2
ssp_coefficient(alg::SSPRK432) = 2
ssp_coefficient(alg::SSPRK932) = 6
ssp_coefficient(alg::SSPRK54) = 1.508
ssp_coefficient(alg::SSPRK104) = 6
ssp_coefficient(alg::SSPRK33) = 1
ssp_coefficient(alg::SSPRK22) = 1
ssp_coefficient(alg::SSPRKMSVS32) = 0.5
ssp_coefficient(alg::SSPRKMSVS43) = 0.33
ssp_coefficient(alg::KYKSSPRK42) = 2.459
ssp_coefficient(alg::KYK2014DGSSPRK_3S2) = 0.8417