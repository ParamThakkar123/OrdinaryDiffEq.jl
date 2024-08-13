function initialize!(integrator, cache::Union{RosenbrockCache})
    integrator.kshortsize = 2
    @unpack k₁, k₂, fsalfirst, fsallast = cache
    integrator.fsalfirst = fsalfirst
    integrator.fsallast = fsallast
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = k₁
    integrator.k[2] = k₂
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    integrator.stats.nf += 1
end

function initialize!(integrator,
        cache::Union{RosenbrockConstantCache})
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    integrator.stats.nf += 1

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = zero(integrator.fsalfirst)
    integrator.k[2] = zero(integrator.fsalfirst)
end

@muladd function perform_step!(integrator, cache::RosenbrockCache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p, opts = integrator
    @unpack k1, k2, k3, du, du1, du2, f₁, fsalfirst, fsallast, dT, J, W, tmp, uf, tf, linsolve_tmp, jac_config, atmp, weight, stage_limiter!, step_limiter! = cache
    @unpack a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, c₃₂, d, a32, C21, C31, C32, C41, C42, C43, C51, C52, C53, C54, C61, C62, C63, C64, C65, C71, C72, C73, C74, C75, C76, C81, C82, C83, C84, C85, C86, C87, b1, b2, b3, btilde1, btilde2, btilde3, gamma, c2, c3, c4, c5, d1, d2, d3, d4, d5 = cache.tab

    # Assignments
    sizeu = size(u)
    mass_matrix = integrator.f.mass_matrix
    utilde = du
    uidx = eachindex(integrator.uprev)
    mass_matrix = integrator.f.mass_matrix

    dtC21 = C21 / dt
    dtC31 = C31 / dt
    dtC32 = C32 / dt
    dtC41 = C41 / dt
    dtC42 = C42 / dt
    dtC43 = C43 / dt
    dtC51 = C51 / dt
    dtC52 = C52 / dt
    dtC53 = C53 / dt
    dtC54 = C54 / dt
    dtC61 = C61 / dt
    dtC62 = C62 / dt
    dtC63 = C63 / dt
    dtC64 = C64 / dt
    dtC65 = C65 / dt
    dtC71 = C71 / dt
    dtC72 = C72 / dt
    dtC73 = C73 / dt
    dtC74 = C74 / dt
    dtC75 = C75 / dt
    dtC76 = C76 / dt
    dtC81 = C81 / dt
    dtC82 = C82 / dt
    dtC83 = C83 / dt
    dtC84 = C84 / dt
    dtC85 = C85 / dt
    dtC86 = C86 / dt
    dtC87 = C87 / dt

    dtd1 = dt * d1
    dtd2 = dt * d2
    dtd3 = dt * d3
    dtd4 = dt * d4
    dtd5 = dt * d5
    dtgamma = dt * gamma

    f(cache.fsalfirst, uprev, p, t) # used in calc_rosenbrock_differentiation!
    integrator.stats.nf += 1

    calc_rosenbrock_differentiation!(integrator, cache, dtd1, dtgamma, repeat_step, true)
    calc_rosenbrock_differentiation!(integrator, cache, dtd4, dtgamma, repeat_step, true)

    # Precalculations
    γ = dt * d
    dto2 = dt / 2
    dto6 = dt / 6

    if repeat_step
        f(integrator.fsalfirst, uprev, p, t)
        integrator.stats.nf += 1
    end

    calc_rosenbrock_differentiation!(integrator, cache, γ, γ, repeat_step, false)

    calculate_residuals!(weight, fill!(weight, one(eltype(u))), uprev, uprev,
        integrator.opts.abstol, integrator.opts.reltol,
        integrator.opts.internalnorm, t)

    if repeat_step
        linres = dolinsolve(
            integrator, cache.linsolve; A = nothing, b = _vec(linsolve_tmp),
            du = integrator.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = γ))
    else
        linres = dolinsolve(integrator, cache.linsolve; A = W, b = _vec(linsolve_tmp),
            du = integrator.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = γ))
    end

    vecu = _vec(linres.u)
    veck₁ = _vec(k₁)

    @.. broadcast=false $(_vec(k1))=-linres.u
    integrator.stats.nsolve += 1

    @.. broadcast=false veck₁=-vecu
    integrator.stats.nsolve += 1

    @.. broadcast=false u=uprev + dto2 * k₁
    stage_limiter!(u, integrator, p, t + dto2)
    f(f₁, u, p, t + dto2)
    integrator.stats.nf += 1

    @.. broadcast=false u=uprev + a21 * k1
    stage_limiter!(u, integrator, p, t + c2 * dt)
    f(du, u, p, t + c2 * dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        copyto!(tmp, k₁)
        tmp .= k₁
    else
        mul!(_vec(tmp), mass_matrix, _vec(k₁))
    end

    @.. broadcast=false linsolve_tmp=f₁ - tmp

    if mass_matrix === I
        @.. broadcast=false linsolve_tmp=du + dtd2 * dT + dtC21 * k1
    else
        @.. broadcast=false du1=dtC21 * k1
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast=false linsolve_tmp=du + dtd2 * dT + du2
    end

    @.. broadcast=false u=uprev + a31 * k1 + a32 * k2
    stage_limiter!(u, integrator, p, t + c3 * dt)
    f(du, u, p, t + c3 * dt)
    integrator.stats.nf += 1

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    vecu = _vec(linres.u)
    veck2 = _vec(k2)
    @.. broadcast=false veck2=-vecu
    integrator.stats.nsolve += 1

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    @.. broadcast=false $(_vec(k2))=-linres.u
    integrator.stats.nsolve += 1

    @.. broadcast=false u=uprev + a31 * k1 + a32 * k2
    stage_limiter!(u, integrator, p, t + c3 * dt)
    f(du, u, p, t + c3 * dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        @.. broadcast=false linsolve_tmp=cache.fsalfirst + dtd3 * dT +
                                         (dtC31 * k1 + dtC32 * k2)
    else
        @.. broadcast=false du1=dtC31 * k1 + dtC32 * k2
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast=false linsolve_tmp=cache.fsalfirst + dtd3 * dT + du2
    end

    @.. broadcast=false veck2=-vecu
    integrator.stats.nsolve += 1

    @.. broadcast=false k₂+=k₁
    @.. broadcast=false u=uprev + dt * k₂
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)

    if mass_matrix === I
        @.. broadcast=false linsolve_tmp=du + dtd3 * dT + dtC31 * k1 + dtC32 * k2
    else
        @.. broadcast=false du1=dtC31 * k1 + dtC32 * k2
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast=false linsolve_tmp=du + dtd3 * dT + du2
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    vecu = _vec(linres.u)
    veck3 = _vec(k3)
    @.. broadcast=false veck3=-vecu
    integrator.stats.nsolve += 1

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    @.. broadcast=false $(_vec(k3))=-linres.u
    integrator.stats.nsolve += 1

    @.. broadcast=false u=uprev + a41 * k1 + a42 * k2 + a43 * k3
    stage_limiter!(u, integrator, p, t + c2 * dt)
    f(du, u, p, t + dt)
    integrator.stats.nf += 1

    @.. broadcast=false veck3=-vecu
    integrator.stats.nsolve += 1

    if mass_matrix === I
        @.. broadcast=false linsolve_tmp=du + dtd4 * dT +
                                         (dtC41 * k1 + dtC42 * k2 + dtC43 * k3)
    else
        @.. broadcast=false du1=dtC41 * k1 + dtC42 * k2 + dtC43 * k3
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast=false linsolve_tmp=du + dtd4 * dT + du2
    end

    if mass_matrix === I
        @.. broadcast=false linsolve_tmp=du +
                                         (dtC41 * k1 + dtC42 * k2 + dtC43 * k3)
    else
        @.. broadcast=false du1=dtC41 * k1 + dtC42 * k2 + dtC43 * k3
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast=false linsolve_tmp=du + du2
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    veck4 = _vec(k4)
    @.. broadcast=false veck4=-vecu
    integrator.stats.nsolve += 1

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    @.. broadcast=false $(_vec(k4))=-linres.u
    integrator.stats.nsolve += 1

    @.. broadcast=false u=uprev + a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4
    stage_limiter!(u, integrator, p, t + dt)
    f(du, u, p, t + dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        @.. broadcast=false linsolve_tmp=du + dtd5 * dT +
                                         (dtC52 * k2 + dtC54 * k4 + dtC51 * k1 + dtC53 * k3)
    else
        @.. broadcast=false du1=dtC52 * k2 + dtC54 * k4 + dtC51 * k1 + dtC53 * k3
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast=false linsolve_tmp=du + dtd5 * dT + du2
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    veck5 = _vec(k5)
    @.. broadcast=false veck5=-vecu
    integrator.stats.nsolve += 1

    if mass_matrix === I
        @.. broadcast=false linsolve_tmp=du +
                                         (dtC52 * k2 + dtC54 * k4 + dtC51 * k1 + dtC53 * k3)
    else
        @.. broadcast=false du1=dtC52 * k2 + dtC54 * k4 + dtC51 * k1 + dtC53 * k3
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast=false linsolve_tmp=du + du2
    end

    @.. broadcast=false u=uprev + a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5
    stage_limiter!(u, integrator, p, t + dt)
    f(du, u, p, t + dt)
    integrator.stats.nf += 1

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    @.. broadcast=false $(_vec(k5))=-linres.u
    integrator.stats.nsolve += 1

    du = u + k4 #-- p=2 solution
    u .+= k5

    @.. broadcast=false u=uprev + b1 * k1 + b2 * k2 + b3 * k3
    step_limiter!(u, integrator, p, t + dt)

    f(fsallast, u, p, t + dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        @.. broadcast=false linsolve_tmp=du + (dtC61 * k1 + dtC62 * k2 + dtC65 * k5 +
                                          dtC64 * k4 + dtC63 * k3)
    else
        @.. broadcast=false du1=dtC61 * k1 + dtC62 * k2 + dtC65 * k5 + dtC64 * k4 +
                                dtC63 * k3
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast=false linsolve_tmp=du + du2
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    veck6 = _vec(k6)
    @.. broadcast=false veck6=-vecu
    integrator.stats.nsolve += 1

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    @.. broadcast=false $(_vec(k6))=-linres.u
    integrator.stats.nsolve += 1

    u .+= k6
    step_limiter!(u, integrator, p, t + dt)
    f(du, u, p, t + dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        @.. broadcast=false linsolve_tmp=du + (dtC71 * k1 + dtC72 * k2 + dtC73 * k3 +
                                          dtC74 * k4 + dtC75 * k5 + dtC76 * k6)
    else
        @.. broadcast=false du1=dtC71 * k1 + dtC72 * k2 + dtC73 * k3 + dtC74 * k4 +
                                dtC75 * k5 + dtC76 * k6
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast=false linsolve_tmp=du + du2
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    veck7 = _vec(k7)
    @.. broadcast=false veck7=-vecu
    integrator.stats.nsolve += 1

    u .+= k7
    f(du, u, p, t + dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        @.. broadcast=false linsolve_tmp=du + (dtC81 * k1 + dtC82 * k2 + dtC83 * k3 +
                                          dtC84 * k4 + dtC85 * k5 + dtC86 * k6 + dtC87 * k7)
    else
        @.. broadcast=false du1=dtC81 * k1 + dtC82 * k2 + dtC83 * k3 + dtC84 * k4 +
                                dtC85 * k5 + dtC86 * k6 + dtC87 * k7
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast=false linsolve_tmp=du + du2
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    veck8 = _vec(k8)
    @.. broadcast=false veck8=-vecu
    integrator.stats.nsolve += 1

    du .= k8
    u .+= k8

    step_limiter!(u, integrator, p, t + dt)

    if integrator.opts.adaptive
        if (integrator.alg isa Rodas5Pe)
            @. du = 0.2606326497975715 * k1 - 0.005158627295444251 * k2 +
                    1.3038988631109731 * k3 + 1.235000722062074 * k4 +
                    -0.7931985603795049 * k5 - 1.005448461135913 * k6 -
                    0.18044626132120234 * k7 + 0.17051519239113755 * k8
        end
        calculate_residuals!(atmp, du, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    EEst = 0.0
    if integrator.opts.calck
        @unpack h21, h22, h23, h24, h25, h31, h32, h33, h34, h35, h2_21, h2_22, h2_23, h2_24, h2_25 = cache.tab
        @.. broadcast=false integrator.k[1]=h21 * k1 + h22 * k2 + h23 * k3 + h24 * k4 +
                                            h25 * k5
        @.. broadcast=false integrator.k[2]=h31 * k1 + h32 * k2 + h33 * k3 + h34 * k4 +
                                            h35 * k5
        @.. broadcast=false integrator.k[3]=h2_21 * k1 + h2_22 * k2 + h2_23 * k3 +
                                            h2_24 * k4 + h2_25 * k5
        if integrator.opts.adaptive
            calculate_interpoldiff!(
                du1, du2, uprev, du, u, integrator.k[1], integrator.k[2], integrator.k[3])
            calculate_residuals!(atmp, du2, uprev, du1, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm, t)
            EEst = max(EEst, integrator.opts.internalnorm(atmp, t))  #-- role of t unclear
        end

        if (integrator.alg isa Rodas5Pr) && integrator.opts.adaptive &&
            (integrator.EEst < 1.0)
             k2 = 0.5 * (uprev + u +
                   0.5 * (integrator.k[1] + 0.5 * (integrator.k[2] + 0.5 * integrator.k[3])))
             du1 = (0.25 * (integrator.k[2] + integrator.k[3]) - uprev + u) / dt
             f(du, k2, p, t + dt / 2)
             integrator.stats.nf += 1
             if mass_matrix === I
                 du2 = du1 - du
             else
                 mul!(_vec(du2), mass_matrix, _vec(du1))
                 du2 = du2 - du
             end
             EEst = norm(du2) / norm(integrator.opts.abstol .+ integrator.opts.reltol .* k2)
             integrator.EEst = max(EEst, integrator.EEst)
         end
    end

    if (integrator.alg isa Rodas23W)
        du1[:] = u[:]
        u[:] = du[:]
        du[:] = du1[:]
        if integrator.opts.calck
            integrator.k[1][:] = integrator.k[3][:]
            integrator.k[2][:] .= 0.0
        end
    end

    if integrator.opts.adaptive
        f(fsallast, u, p, t + dt)
        integrator.stats.nf += 1

        if mass_matrix === I
            @.. broadcast=false linsolve_tmp=fsallast - c₃₂ * (k₂ - f₁) -
                                             2(k₁ - fsalfirst) + dt * dT
        else
            @.. broadcast=false du2=c₃₂ * k₂ + 2k₁
            mul!(_vec(du1), mass_matrix, _vec(du2))
            @.. broadcast=false linsolve_tmp=fsallast - du1 + c₃₂ * f₁ + 2fsalfirst +
                                             dt * dT
        end

        linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
        vecu = _vec(linres.u)
        veck3 = _vec(k₃)
        @.. broadcast=false veck3=-vecu

        integrator.stats.nsolve += 1

        @.. broadcast=false u=uprev + dto6 * (k₁ + 4k₂ + k₃)
        step_limiter!(u, integrator, p, t + dt)

        if integrator.opts.adaptive
            calculate_residuals!(atmp, u - du, uprev, u, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm, t)
            integrator.EEst = max(EEst, integrator.opts.internalnorm(atmp, t))
        end

        if mass_matrix === I
            @.. broadcast=false tmp=dto6 * (k₁ - 2 * k₂ + k₃)
        else
            veck₁ = _vec(k₁)
            veck₂ = _vec(k₂)
            veck₃ = _vec(k₃)
            vectmp = _vec(tmp)
            @.. broadcast=false vectmp=ifelse(cache.algebraic_vars,
                false, dto6 * (veck₁ - 2 * veck₂ + veck₃))
        end
        calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)

        if mass_matrix !== I
            algvar = reshape(cache.algebraic_vars, size(u))
            @.. broadcast=false atmp=ifelse(algvar, fsallast, false) /
                                     integrator.opts.abstol
            integrator.EEst += integrator.opts.internalnorm(atmp, t)
        end
    end

    if integrator.opts.adaptive
        calculate_residuals!(atmp, k6, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    if integrator.opts.calck
        @unpack h21, h22, h23, h24, h25, h31, h32, h33, h34, h35 = cache.tab
        @.. broadcast=false integrator.k[1]=h21 * k1 + h22 * k2 + h23 * k3 + h24 * k4 +
                                            h25 * k5
        @.. broadcast=false integrator.k[2]=h31 * k1 + h32 * k2 + h33 * k3 + h34 * k4 +
                                            h35 * k5
    end
    
    cache.linsolve = linres.cache
end

@muladd function perform_step!(integrator, cache::RosenbrockConstantCache,
        repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack c₃₂, d, tf, uf = cache
    @unpack a21, a31, a32, a41, a42, a43, C21, C31, C32, C41, C42, C43, C51, C52, C53, C54, b1, b2, b3, btilde1, btilde2, btilde3, gamma, c2, c3, d1, d2, d3 = cache.tab

    # Precalculations
    γ = dt * d
    dto2 = dt / 2
    dto6 = dt / 6

    dtC21 = C21 / dt
    dtC31 = C31 / dt
    dtC32 = C32 / dt
    dtC41 = C41 / dt
    dtC42 = C42 / dt
    dtC43 = C43 / dt
    dtC51 = C51 / dt
    dtC52 = C52 / dt
    dtC53 = C53 / dt
    dtC54 = C54 / dt

    dtd1 = dt * d1
    dtd2 = dt * d2
    dtd3 = dt * d3
    dtgamma = dt * gamma

    if repeat_step
        integrator.fsalfirst = f(uprev, p, t)
        integrator.stats.nf += 1
    end

    mass_matrix = integrator.f.mass_matrix

    # Time derivative
    tf.u = uprev
    dT = calc_tderivative(integrator, cache)

    W = calc_W(integrator, cache, dtgamma, repeat_step, true)
    if !issuccess_W(W)
        integrator.EEst = 2
        return nothing
    end

    W = calc_W(integrator, cache, γ, repeat_step)
    if !issuccess_W(W)
        integrator.EEst = 2
        return nothing
    end

    du = f(uprev, p, t)
    integrator.stats.nf += 1
    k3 = copy(du)  #-- save for stage 3

    linsolve_tmp = du + dtd1 * dT

    k1 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    u = uprev + a21 * k1
    du = f(u, p, t + c2 * dt)
    integrator.stats.nf += 1

    linsolve_tmp = integrator.fsalfirst + dtd1 * dT

    k1 = _reshape(W \ -_vec((integrator.fsalfirst + γ * dT)), axes(uprev))
    integrator.stats.nsolve += 1
    u = uprev + a21 * k1
    du = f(u, p, t + c2 * dt)
    f₁ = f(uprev + dto2 * k1, p, t + dto2)
    integrator.stats.nf += 1

    if mass_matrix === I
        linsolve_tmp = du + dtd2 * dT + dtC21 * k1
    else
        linsolve_tmp = du + dtd2 * dT + mass_matrix * (dtC21 * k1)
    end

    k2 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    u = uprev + a31 * k1 + a32 * k2
    du = f(u, p, t + c3 * dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        k₂ = _reshape(W \ -_vec(f₁ - k₁), axes(uprev)) + k₁
    else
        k₂ = _reshape(W \ -_vec(f₁ - mass_matrix * k₁), axes(uprev)) + k₁
    end
    integrator.stats.nsolve += 1
    u = uprev + dt * k₂

    if mass_matrix === I
        linsolve_tmp = du + dtd3 * dT + dtC31 * k1 + dtC32 * k2
    else
        linsolve_tmp = du + dtd3 * dT + mass_matrix * (dtC31 * k1 + dtC32 * k2)
    end

    k3 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    u = uprev + b1 * k1 + b2 * k2 + b3 * k3
    u = uprev + a41 * k1 + a42 * k2 + a43 * k3
    integrator.fsallast = f(u, p, t + dt)
    du = f(u, p, t + dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        linsolve_tmp = du + (dtC41 * k1 + dtC42 * k2 + dtC43 * k3)
    else
        linsolve_tmp = du + mass_matrix * (dtC41 * k1 + dtC42 * k2 + dtC43 * k3)
    end

    k4 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1

    if mass_matrix === I
        linsolve_tmp = du + (dtC52 * k2 + dtC54 * k4 + dtC51 * k1 + dtC53 * k3)
    else
        linsolve_tmp = du +
                       mass_matrix * (dtC52 * k2 + dtC54 * k4 + dtC51 * k1 + dtC53 * k3)
    end

    k5 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    du = u + k4 #-- solution p=2
    u = u + k5 #-- solution p=3

    EEst = 0.0
    if integrator.opts.calck
        @unpack h21, h22, h23, h24, h25, h31, h32, h33, h34, h35, h2_21, h2_22, h2_23, h2_24, h2_25 = cache.tab
        integrator.k[1] = h21 * k1 + h22 * k2 + h23 * k3 + h24 * k4 + h25 * k5
        integrator.k[2] = h31 * k1 + h32 * k2 + h33 * k3 + h34 * k4 + h35 * k5
        integrator.k[3] = h2_21 * k1 + h2_22 * k2 + h2_23 * k3 + h2_24 * k4 + h2_25 * k5
        if integrator.opts.adaptive
            if isa(linsolve_tmp, AbstractFloat)
                u_int, u_diff = calculate_interpoldiff(
                    uprev, du, u, integrator.k[1], integrator.k[2], integrator.k[3])
            else
                u_int = linsolve_tmp
                u_diff = linsolve_tmp .+ 0
                calculate_interpoldiff!(u_int, u_diff, uprev, du, u, integrator.k[1],
                    integrator.k[2], integrator.k[3])
            end
            atmp = calculate_residuals(u_diff, uprev, u_int, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm, t)
            EEst = max(EEst, integrator.opts.internalnorm(atmp, t))  #-- role of t unclear
        end
    end

    if (integrator.alg isa Rodas23W)
        k1 = u .+ 0
        u = du .+ 0
        du = k1 .+ 0
        if integrator.opts.calck
            integrator.k[1] = integrator.k[3] .+ 0
            integrator.k[2] = 0 * integrator.k[2]
        end
    end

    if integrator.opts.adaptive
        atmp = calculate_residuals(u - du, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = max(EEst, integrator.opts.internalnorm(atmp, t))
    end

    if integrator.opts.adaptive
        utilde = btilde1 * k1 + btilde2 * k2 + btilde3 * k3
        atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    if integrator.opts.adaptive
        integrator.fsallast = f(u, p, t + dt)
        integrator.stats.nf += 1

        if mass_matrix === I
            k₃ = _reshape(
                W \
                -_vec((integrator.fsallast - c₃₂ * (k₂ - f₁) -
                       2 * (k₁ - integrator.fsalfirst) + dt * dT)),
                axes(uprev))
        else
            linsolve_tmp = integrator.fsallast - mass_matrix * (c₃₂ * k₂ + 2 * k₁) +
                           c₃₂ * f₁ + 2 * integrator.fsalfirst + dt * dT
            k₃ = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
        end
        integrator.stats.nsolve += 1
        u = uprev + dto6 * (k₁ + 4k₂ + k₃)

        if integrator.opts.adaptive
            utilde = dto6 * (k₁ - 2k₂ + k₃)
            atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm, t)
            integrator.EEst = integrator.opts.internalnorm(atmp, t)
    
            if mass_matrix !== I
                atmp = @. ifelse(!integrator.differential_vars, integrator.fsallast, false) ./
                          integrator.opts.abstol
                integrator.EEst += integrator.opts.internalnorm(atmp, t)
            end
        end

        if u isa Number
            utilde = dto6 * f.mass_matrix[1, 1] * (k₁ - 2 * k₂ + k₃)
        else
            utilde = dto6 * f.mass_matrix * (k₁ - 2 * k₂ + k₃)
        end
        atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)

        if mass_matrix !== I
            atmp = @. ifelse(!integrator.differential_vars, integrator.fsallast, false) ./
                      integrator.opts.abstol
            integrator.EEst += integrator.opts.internalnorm(atmp, t)
        end
    end
    integrator.k[1] = k₁
    integrator.k[2] = k₂
    integrator.u = u
    return nothing
end

@muladd function perform_step!(integrator, cache::RosenbrockConstantCache,
        repeat_step = false)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
    return nothing
end

################################################################################


@muladd function perform_step!(integrator, cache::RosenbrockCache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack du, du1, du2, fsalfirst, fsallast, k1, k2, k3, k4, dT, J, W, uf, tf, linsolve_tmp, jac_config, atmp, weight, stage_limiter!, step_limiter! = cache
    @unpack A, C, b, btilde, gamma, c, d = cache.tab  # Coefficients from the tableau

    # Assignments
    sizeu = size(u)
    mass_matrix = integrator.f.mass_matrix
    utilde = du

    # Precalculations
    dtC = C ./ dt
    dtd = dt .* d
    dtgamma = dt * gamma

    # Differentiation
    calc_rosenbrock_differentiation!(integrator, cache, dtd[1], dtgamma, repeat_step, true)
    
    calculate_residuals!(weight, fill!(weight, one(eltype(u))), uprev, uprev,
        integrator.opts.abstol, integrator.opts.reltol,
        integrator.opts.internalnorm, t)

    # Linear solve setup
    if repeat_step
        linres = dolinsolve(integrator, cache.linsolve; A = nothing, b = _vec(linsolve_tmp),
            du = integrator.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = dtgamma))
    else
        linres = dolinsolve(integrator, cache.linsolve; A = W, b = _vec(linsolve_tmp),
            du = integrator.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = dtgamma))
    end

    vecu = _vec(linres.u)
    k = Vector{typeof(k1)}(undef, length(A))
    
    # Stage calculation loop
    for i in 1:size(A)
        @.. broadcast=false k[i] = -vecu
        integrator.stats.nsolve += 1

        if i < length(A)  # Skip the last iteration for k computation
            u = uprev + sum(A[i,j] * k[j] for j in 1:(i-1))
            du = f(u, p, t + c[i] * dt)
            integrator.stats.nf += 1

            if mass_matrix === I
                linsolve_tmp = du + dtd[i+1] * dT + dtC[i] * k[i]
            else
                du1 = dtC[i] * k[i]
                mul!(_vec(du2), mass_matrix, _vec(du1))
                linsolve_tmp = du + dtd[i+1] * dT + du2
            end

            linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
            vecu = _vec(linres.u)
        end
    end

    # Combine stages to get the final solution
    u = uprev + sum(b[i] * k[i] for i in 1:size(b))
    step_limiter!(u, integrator, p, t + dt)

    integrator.fsallast = f(u, p, t + dt)
    integrator.stats.nf += 1

    if integrator.opts.adaptive
        utilde = sum(btilde[i] * k[i] for i in 1:size(btilde))
        calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    cache.linsolve = linres.cache
    return nothing
end

################################################################################

#### ROS2 type method

@ROS2(:init)
@ROS2(:performstep)

################################################################################

#### ROS23 type method

@ROS23(:init)
@ROS23(:performstep)

################################################################################

#### ROS34PW type method

@ROS34PW(:init)
@ROS34PW(:performstep)

################################################################################

#### ROS4 type method

@Rosenbrock4(:performstep)

################################################################################

#### Rodas3P type method

function initialize!(integrator, cache::Union{RosenbrockConstantCache})
    integrator.kshortsize = 3
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    # Avoid undefined entries if k is an array of arrays
    integrator.k[1] = zero(integrator.u)
    integrator.k[2] = zero(integrator.u)
    integrator.k[3] = zero(integrator.u)
end

function initialize!(integrator, cache::Union{RosenbrockCache})
    integrator.kshortsize = 3
    @unpack dense1, dense2, dense3 = cache
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = dense1
    integrator.k[2] = dense2
    integrator.k[3] = dense3
end

function calculate_interpoldiff(uprev, up2, up3, c_koeff, d_koeff, c2_koeff)
    u_int = 0.0
    u_diff = 0.0
    a1 = up3 + c_koeff - up2 - c2_koeff
    a2 = d_koeff - c_koeff + c2_koeff
    a3 = -d_koeff
    dis = a2^2 - 3 * a1 * a3
    u_int = up3
    u_diff = 0.0
    if dis > 0.0 #-- Min/Max occurs
        tau1 = (-a2 - sqrt(dis)) / (3 * a3)
        tau2 = (-a2 + sqrt(dis)) / (3 * a3)
        if tau1 > tau2
            tau1, tau2 = tau2, tau1
        end
        for tau in (tau1, tau2)
            if (tau > 0.0) && (tau < 1.0)
                y_tau = (1 - tau) * uprev +
                        tau * (up3 + (1 - tau) * (c_koeff + tau * d_koeff))
                dy_tau = ((a3 * tau + a2) * tau + a1) * tau
                if abs(dy_tau) > abs(u_diff)
                    u_diff = dy_tau
                    u_int = y_tau
                end
            end
        end
    end
    return u_int, u_diff
end

function calculate_interpoldiff!(u_int, u_diff, uprev, up2, up3, c_koeff, d_koeff, c2_koeff)
    for i in eachindex(up2)
        a1 = up3[i] + c_koeff[i] - up2[i] - c2_koeff[i]
        a2 = d_koeff[i] - c_koeff[i] + c2_koeff[i]
        a3 = -d_koeff[i]
        dis = a2^2 - 3 * a1 * a3
        u_int[i] = up3[i]
        u_diff[i] = 0.0
        if dis > 0.0 #-- Min/Max occurs
            tau1 = (-a2 - sqrt(dis)) / (3 * a3)
            tau2 = (-a2 + sqrt(dis)) / (3 * a3)
            if tau1 > tau2
                tau1, tau2 = tau2, tau1
            end
            for tau in (tau1, tau2)
                if (tau > 0.0) && (tau < 1.0)
                    y_tau = (1 - tau) * uprev[i] +
                            tau * (up3[i] + (1 - tau) * (c_koeff[i] + tau * d_koeff[i]))
                    dy_tau = ((a3 * tau + a2) * tau + a1) * tau
                    if abs(dy_tau) > abs(u_diff[i])
                        u_diff[i] = dy_tau
                        u_int[i] = y_tau
                    end
                end
            end
        end
    end
end

#### Rodas4 type method

function initialize!(integrator, cache::RosenbrockConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    # Avoid undefined entries if k is an array of arrays
    integrator.k[1] = zero(integrator.u)
    integrator.k[2] = zero(integrator.u)
end

@muladd function perform_step!(integrator, cache::RosenbrockConstantCache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack tf, uf = cache
    @unpack gamma, c, d, C, a = cache.tab

    # Precalculations
    dtC = C / dt
    dtd = dt * d
    dtgamma = dt * gamma

    mass_matrix = integrator.f.mass_matrix

    # Time derivative
    tf.u = uprev
    dT = calc_tderivative(integrator, cache)

    W = calc_W(integrator, cache, dtgamma, repeat_step, true)
    if !issuccess_W(W)
        integrator.EEst = 2
        return nothing
    end

    du = f(uprev, p, t)
    integrator.stats.nf += 1

    k = Array{typeof(uprev)}(undef, length(d))

    linsolve_tmp = du + dtd1 * dT

    k[1] = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    u = uprev + a[2, 1] * k[1]
    du = f(u, p, t + c[2] * dt)
    integrator.stats.nf += 1

    for i in 2:length(d)
        if mass_matrix === I
            linsolve_tmp = du + dtd[i] * dT + sum(dtC[i, 1:i-1] .* k[1:i-1])
        else
            linsolve_tmp = du + dtd[i] * dT + mass_matrix * sum(dtC[i, 1:i-1] .* k[1:i-1])
        end

        k[i] = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
        integrator.stats.nsolve += 1
        u += sum(a[i+1, 1:i] .* k[1:i])
        du = f(u, p, t + c[i+1] * dt)
        integrator.stats.nf += 1
    end

    u += k[end]
    du = f(u, p, t + dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        linsolve_tmp = du + sum(dtC[end, :] .* k)
    else
        linsolve_tmp = du + mass_matrix * sum(dtC[end, :] .* k)
    end

    k[end] = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    u += k[end]

    if integrator.opts.adaptive
        atmp = calculate_residuals(k[end], uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    if integrator.opts.calck
        @unpack h = cache.tab
        integrator.k[1] = sum(h[1, :] .* k)
        integrator.k[2] = sum(h[2, :] .* k)
    end
    integrator.u = u
    return nothing
end

function initialize!(integrator, cache::RosenbrockCache)
    integrator.kshortsize = 2
    @unpack dense1, dense2 = cache
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = dense1
    integrator.k[2] = dense2
end

###############################################################################

### Rodas5 Method

function initialize!(integrator, cache::RosenbrockConstantCache)
    integrator.kshortsize = 3
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    # Avoid undefined entries if k is an array of arrays
    integrator.k[1] = zero(integrator.u)
    integrator.k[2] = zero(integrator.u)
    integrator.k[3] = zero(integrator.u)
end

@muladd function perform_step!(integrator, cache::RosenbrockConstantCache,
        repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack tf, uf = cache
    @unpack A, C, b, d, c, gamma = cache.tab  # Coefficients from the tableau

    # Precalculations
    dtC = C ./ dt
    dtd = dt .* d
    dtgamma = dt * gamma

    mass_matrix = integrator.f.mass_matrix

    # Time derivative
    dT = calc_tderivative(integrator, cache)

    W = calc_W(integrator, cache, dtgamma, repeat_step, true)
    if !issuccess_W(W)
        integrator.EEst = 2
        return nothing
    end

    du1 = f(uprev, p, t)
    integrator.stats.nf += 1

    linsolve_tmp = du1 + dtd[1] * dT
    k = Vector{typeof(du1)}(undef, length(A))

    k[1] = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1

    for i in 2:length(A)
        u = uprev + sum(A[i-1, j] * k[j] for j in 1:(i-1))
        du = f(u, p, t + c[i] * dt)
        integrator.stats.nf += 1

        if mass_matrix === I
            linsolve_tmp = du + dtd[i] * dT + dtC[i] * k[i]
        else
            linsolve_tmp = du + dtd[i] * dT + mass_matrix * (dtC[i] * k[i])
        end

        k[i] = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
        integrator.stats.nsolve += 1
    end

    u = uprev + sum(A[end, i] * k[i] for i in 1:length(k))
    du = f(u, p, t + dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        linsolve_tmp = du + sum(dtC[i] * k[i] for i in 1:length(k))
    else
        linsolve_tmp = du + mass_matrix * sum(dtC[i] * k[i] for i in 1:length(k))
    end

    k[end+1] = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1

    u = u + k[end+1]
    du = f(u, p, t + dt)
    integrator.stats.nf += 1

    if integrator.opts.adaptive
        if (integrator.alg isa Rodas5Pe)
            linsolve_tmp = sum(btilde[i] * k[i] for i in 1:length(k))
        end
        atmp = calculate_residuals(linsolve_tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    if integrator.opts.calck
        @unpack h = cache.tab
        integrator.k[1] = sum(h[i] * k[i] for i in 1:length(k))
        integrator.k[2] = sum(h[length(k)+i] * k[i] for i in 1:length(k))
        integrator.k[3] = sum(h[2*length(k)+i] * k[i] for i in 1:length(k))
        if (integrator.alg isa Rodas5Pr) && integrator.opts.adaptive &&
           (integrator.EEst < 1.0)
            k2 = 0.5 * (uprev + u +
                  0.5 * (integrator.k[1] + 0.5 * (integrator.k[2] + 0.5 * integrator.k[3])))
            du1 = (0.25 * (integrator.k[2] + integrator.k[3]) - uprev + u) / dt
            du = f(k2, p, t + dt / 2)
            integrator.stats.nf += 1
            if mass_matrix === I
                du2 = du1 - du
            else
                du2 = mass_matrix * du1 - du
            end
            EEst = norm(du2) / norm(integrator.opts.abstol .+ integrator.opts.reltol .* k2)
            integrator.EEst = max(EEst, integrator.EEst)
        end
    end

    integrator.u = u
    return nothing
end

function initialize!(integrator, cache::RosenbrockCache)
    integrator.kshortsize = 3
    @unpack dense1, dense2, dense3 = cache
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = dense1
    integrator.k[2] = dense2
    integrator.k[3] = dense3
end

@RosenbrockW6S4OS(:init)
@RosenbrockW6S4OS(:performstep)