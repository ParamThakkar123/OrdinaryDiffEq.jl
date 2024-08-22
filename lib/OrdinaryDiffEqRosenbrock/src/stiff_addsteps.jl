function _ode_addsteps!(k, t, uprev, u, dt, f, p,
        cache::Union{Rosenbrock23ConstantCache,
            Rosenbrock32ConstantCache},
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false)
    if length(k) < 2 || always_calc_begin
        @unpack tf, uf, d = cache
        γ = dt * d
        tf.u = uprev
        if cache.autodiff isa AutoForwardDiff
            dT = ForwardDiff.derivative(tf, t)
        else
            dT = FiniteDiff.finite_difference_derivative(tf, t, dir = sign(dt))
        end

        mass_matrix = f.mass_matrix
        if uprev isa AbstractArray
            J = ForwardDiff.jacobian(uf, uprev)
            W = mass_matrix - γ * J
        else
            J = ForwardDiff.derivative(uf, uprev)
            W = 1 - γ * J
        end
        k₁ = W \ (f(uprev, p, t) + dt * d * dT)
        f₁ = f(uprev + dt * k₁ / 2, p, t + dt / 2)
        k₂ = W \ (f₁ - k₁) + k₁
        copyat_or_push!(k, 1, k₁)
        copyat_or_push!(k, 2, k₂)
    end
    nothing
end

function _ode_addsteps!(k, t, uprev, u, dt, f, p,
        cache::Union{Rosenbrock23Cache, Rosenbrock32Cache},
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false)
    if length(k) < 2 || always_calc_begin
        @unpack k₁, k₂, k₃, du1, du2, f₁, fsalfirst, fsallast, dT, J, W, tmp, uf, tf, linsolve_tmp, weight = cache
        @unpack c₃₂, d = cache.tab
        uidx = eachindex(uprev)

        # Assignments
        sizeu = size(u)
        mass_matrix = f.mass_matrix
        γ = dt * d
        dto2 = dt / 2
        dto6 = dt / 6

        @.. broadcast=false linsolve_tmp=@muladd fsalfirst + γ * dT

        ### Jacobian does not need to be re-evaluated after an event
        ### Since it's unchanged
        jacobian2W!(W, mass_matrix, γ, J, false)

        linsolve = cache.linsolve

        linres = dolinsolve(cache, linsolve; A = W, b = _vec(linsolve_tmp),
            reltol = cache.reltol)

        vecu = _vec(linres.u)
        veck₁ = _vec(k₁)

        @.. broadcast=false veck₁=-vecu

        @.. broadcast=false tmp=uprev + dto2 * k₁
        f(f₁, tmp, p, t + dto2)

        if mass_matrix === I
            tmp .= k₁
        else
            mul!(_vec(tmp), mass_matrix, _vec(k₁))
        end

        @.. broadcast=false linsolve_tmp=f₁ - tmp

        linres = dolinsolve(cache, linres.cache; b = _vec(linsolve_tmp),
            reltol = cache.reltol)
        vecu = _vec(linres.u)
        veck2 = _vec(k₂)

        @.. broadcast=false veck2=-vecu

        @.. broadcast=false k₂+=k₁

        copyat_or_push!(k, 1, k₁)
        copyat_or_push!(k, 2, k₂)
        cache.linsolve = linres.cache
    end
    nothing
end

function _ode_addsteps!(
        k, t, uprev, u, dt, f, p, cache::Union{Rodas23WConstantCache, Rodas3PConstantCache},
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false)
    if length(k) < 2 || always_calc_begin
        @unpack tf, uf = cache
        @unpack a21, a41, a42, a43, C21, C31, C32, C41, C42, C43, C51, C52, C53, C54, gamma, c2, c3, d1, d2, d3 = cache.tab

        # Precalculations
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
        mass_matrix = f.mass_matrix

        # Time derivative
        tf.u = uprev
        if cache.autodiff isa AutoForwardDiff
            dT = ForwardDiff.derivative(tf, t)
        else
            dT = FiniteDiff.finite_difference_derivative(tf, t, dir = sign(dt))
        end

        # Jacobian
        uf.t = t
        if uprev isa AbstractArray
            J = ForwardDiff.jacobian(uf, uprev)
            W = mass_matrix / dtgamma - J
        else
            J = ForwardDiff.derivative(uf, uprev)
            W = 1 / dtgamma - J
        end

        du = f(uprev, p, t)
        k3 = copy(du)

        linsolve_tmp = du + dtd1 * dT

        k1 = W \ linsolve_tmp
        u = uprev + a21 * k1
        du = f(u, p, t + c2 * dt)

        linsolve_tmp = du + dtd2 * dT + dtC21 * k1

        k2 = W \ linsolve_tmp

        linsolve_tmp = k3 + dtd3 * dT + (dtC31 * k1 + dtC32 * k2)

        k3 = W \ linsolve_tmp
        u = uprev + a41 * k1 + a42 * k2 + a43 * k3
        du = f(u, p, t + dt)

        linsolve_tmp = du + (dtC41 * k1 + dtC42 * k2 + dtC43 * k3)

        k4 = W \ linsolve_tmp

        linsolve_tmp = du + (dtC52 * k2 + dtC54 * k4 + dtC51 * k1 + dtC53 * k3)

        k5 = W \ linsolve_tmp

        @unpack h21, h22, h23, h24, h25, h31, h32, h33, h34, h35, h2_21, h2_22, h2_23, h2_24, h2_25 = cache.tab
        k₁ = h21 * k1 + h22 * k2 + h23 * k3 + h24 * k4 + h25 * k5
        k₂ = h31 * k1 + h32 * k2 + h33 * k3 + h34 * k4 + h35 * k5
        #k₃ = h2_21 * k1 + h2_22 * k2 + h2_23 * k3 + h2_24 * k4 + h2_25 * k5
        copyat_or_push!(k, 1, k₁)
        copyat_or_push!(k, 2, k₂)
        #copyat_or_push!(k, 3, k₃)
    end
    nothing
end

function _ode_addsteps!(
        k, t, uprev, u, dt, f, p, cache::Union{Rodas23WCache, Rodas3PCache},
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false)
    if length(k) < 2 || always_calc_begin
        @unpack du, du1, du2, tmp, k1, k2, k3, k4, k5, dT, J, W, uf, tf, linsolve_tmp, jac_config, fsalfirst, weight = cache
        @unpack a21, a41, a42, a43, C21, C31, C32, C41, C42, C43, C51, C52, C53, C54, gamma, c2, c3, d1, d2, d3 = cache.tab

        # Assignments
        sizeu = size(u)
        uidx = eachindex(uprev)
        mass_matrix = f.mass_matrix

        # Precalculations
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

        @.. broadcast=false linsolve_tmp=@muladd fsalfirst + dtgamma * dT

        ### Jacobian does not need to be re-evaluated after an event
        ### Since it's unchanged
        jacobian2W!(W, mass_matrix, dtgamma, J, true)

        linsolve = cache.linsolve

        linres = dolinsolve(cache, linsolve; A = W, b = _vec(linsolve_tmp),
            reltol = cache.reltol)
        vecu = _vec(linres.u)
        veck1 = _vec(k1)

        @.. broadcast=false veck1=-vecu
        @.. broadcast=false tmp=uprev + a21 * k1
        f(du, tmp, p, t + c2 * dt)

        if mass_matrix === I
            @.. broadcast=false linsolve_tmp=du + dtd2 * dT + dtC21 * k1
        else
            @.. broadcast=false du1=dtC21 * k1
            mul!(du2, mass_matrix, du1)
            @.. broadcast=false linsolve_tmp=du + dtd2 * dT + du2
        end

        linres = dolinsolve(cache, linres.cache; b = _vec(linsolve_tmp),
            reltol = cache.reltol)
        vecu = _vec(linres.u)
        veck2 = _vec(k2)
        @.. broadcast=false veck2=-vecu

        if mass_matrix === I
            @.. broadcast=false linsolve_tmp=fsalfirst + dtd3 * dT +
                                             (dtC31 * k1 + dtC32 * k2)
        else
            @.. broadcast=false du1=dtC31 * k1 + dtC32 * k2
            mul!(du2, mass_matrix, du1)
            @.. broadcast=false linsolve_tmp=fsalfirst + dtd3 * dT + du2
        end

        linres = dolinsolve(cache, linres.cache; b = _vec(linsolve_tmp),
            reltol = cache.reltol)
        vecu = _vec(linres.u)
        veck3 = _vec(k3)
        @.. broadcast=false veck3=-vecu
        @.. broadcast=false tmp=uprev + a41 * k1 + a42 * k2 + a43 * k3
        f(du, tmp, p, t + dt)

        if mass_matrix === I
            @.. broadcast=false linsolve_tmp=du + (dtC41 * k1 + dtC42 * k2 + dtC43 * k3)
        else
            @.. broadcast=false du1=dtC41 * k1 + dtC42 * k2 + dtC43 * k3
            mul!(du2, mass_matrix, du1)
            @.. broadcast=false linsolve_tmp=du + du2
        end

        linres = dolinsolve(cache, linres.cache; b = _vec(linsolve_tmp),
            reltol = cache.reltol)
        vecu = _vec(linres.u)
        veck4 = _vec(k4)
        @.. broadcast=false veck4=-vecu

        if mass_matrix === I
            @.. broadcast=false linsolve_tmp=du + (dtC52 * k2 + dtC54 * k4 + dtC51 * k1 +
                                              dtC53 * k3)
        else
            @.. broadcast=false du1=dtC52 * k2 + dtC54 * k4 + dtC51 * k1 + dtC53 * k3
            mul!(du2, mass_matrix, du1)
            @.. broadcast=false linsolve_tmp=du + du2
        end

        linres = dolinsolve(cache, linres.cache; b = _vec(linsolve_tmp),
            reltol = cache.reltol)
        vecu = _vec(linres.u)
        veck5 = _vec(k5)
        @.. broadcast=false veck5=-vecu
        @unpack h21, h22, h23, h24, h25, h31, h32, h33, h34, h35, h2_21, h2_22, h2_23, h2_24, h2_25 = cache.tab
        @.. broadcast=false du=h21 * k1 + h22 * k2 + h23 * k3 + h24 * k4 + h25 * k5
        copyat_or_push!(k, 1, copy(du))

        @.. broadcast=false du=h31 * k1 + h32 * k2 + h33 * k3 + h34 * k4 + h35 * k5
        copyat_or_push!(k, 2, copy(du))
    end
    nothing
end

function _ode_addsteps!(k, t, uprev, u, dt, f, p, cache::RosenbrockConstantCache,
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false)
    
    if length(k) < 2 || always_calc_begin
        @unpack tf, uf = cache
        @unpack a, C, gamma, c, d = cache.tab

        # Precalculations
        dtC = C ./ dt
        dtd = dt .* d
        dtgamma = dt * gamma
        mass_matrix = f.mass_matrix

        # Time derivative
        tf.u = uprev
        dT = (cache.autodiff isa AutoForwardDiff) ? 
             ForwardDiff.derivative(tf, t) : 
             FiniteDiff.finite_difference_derivative(tf, t, dir = sign(dt))

        # Jacobian
        uf.t = t
        J = (uprev isa AbstractArray) ? 
            ForwardDiff.jacobian(uf, uprev) : 
            ForwardDiff.derivative(uf, uprev)
        
        W = (uprev isa AbstractArray) ? 
            mass_matrix / dtgamma - J : 
            1 / dtgamma - J

        # Initialize k values
        du = f(uprev, p, t)
        linsolve_tmp = du + dtd[1] * dT
        k_values = Vector{eltype(uprev)}(undef, 5)  # Store k1 to k5

        for i in 1:5
            if i > 1
                u = uprev + sum(a[i][j] * k_values[j] for j in 1:(i-1))
                du = f(u, p, t + c[i] * dt)
                linsolve_tmp = du + dtd[i] * dT + sum(dtC[i][j] * k_values[j] for j in 1:(i-1))
            end
            
            k[i] = W \ linsolve_tmp
        end

        @unpack h21, h22, h23, h24, h25, h31, h32, h33, h34, h35 = cache.tab
        k₁ = h21 * k[1] + h22 * k[2] + h23 * k[3] + h24 * k[4] + h25 * k[5]
        k₂ = h31 * k[1] + h32 * k[2] + h33 * k[3] + h34 * k[4] + h35 * k[5]
        
        copyat_or_push!(k, 1, k₁)
        copyat_or_push!(k, 2, k₂)
    end
    nothing
end

function _ode_addsteps!(k, t, uprev, u, dt, f, p, cache::RosenbrockCache,
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false)
    
    if length(k) < 2 || always_calc_begin
        @unpack dus, tmp, ks, dT, J, W, uf, tf, linsolve_tmp, jac_config, fsalfirst, weight = cache
        @unpack a, C, gamma, c, d = cache.tab

        # Assignments
        sizeu = size(u)
        uidx = eachindex(uprev)
        mass_matrix = f.mass_matrix

        # Precalculations
        dtC = C ./ dt
        dtd = dt .* d
        dtgamma = dt * gamma

        @.. broadcast=false linsolve_tmp = @muladd fsalfirst + dtgamma * dT

        # Jacobian does not need to be re-evaluated after an event
        jacobian2W!(W, mass_matrix, dtgamma, J, true)

        linsolve = cache.linsolve

        # Solve for k1 to k5 using a loop
        for i in 1:5
            linsolve_tmp = (i == 1) ? linsolve_tmp : 
                dus[1] + dtd[i] * dT + sum(dtC[i][j] * ks[j] for j in 1:(i-1))

            # Update dus based on the mass matrix
            if mass_matrix === I
                dus[2] = (i == 1) ? dus[1] : dus[2]
            else
                dus[2] = dtC[i][1] * ks[1] + sum(dtC[i][j] * ks[j] for j in 2:(i-1))
                mul!(dus[3], mass_matrix, dus[2])
            end

            # Solve the linear system
            linres = dolinsolve(cache, linsolve; b = _vec(linsolve_tmp),
                reltol = cache.reltol)
            vecu = _vec(linres.u)
            @.. broadcast=false ks[i] = -vecu

            # Update the temporary variable
            tmp = uprev + sum(a[i][j] * ks[j] for j in 1:i)
            f(dus[1], tmp, p, t + c[i] * dt)
        end

        # Final k6 calculations
        @unpack h21, h22, h23, h24, h25, h31, h32, h33, h34, h35 = cache.tab
        ks[6] = h21 * ks[1] + h22 * ks[2] + h23 * ks[3] + h24 * ks[4] + h25 * ks[5]
        copyat_or_push!(k, 1, copy(ks[6]))

        ks[6] = h31 * ks[1] + h32 * ks[2] + h33 * ks[3] + h34 * ks[4] + h35 * ks[5]
        copyat_or_push!(k, 2, copy(ks[6]))
    end
    nothing
end

function _ode_addsteps!(k, t, uprev, u, dt, f, p, cache::Rosenbrock5ConstantCache,
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false)
    if length(k) < 3 || always_calc_begin
        @unpack tf, uf = cache
        @unpack a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, C21, C31, C32, C41, C42, C43, C51, C52, C53, C54, C61, C62, C63, C64, C65, C71, C72, C73, C74, C75, C76, C81, C82, C83, C84, C85, C86, C87, gamma, d1, d2, d3, d4, d5, c2, c3, c4, c5 = cache.tab

        # Precalculations
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
        mass_matrix = f.mass_matrix

        # Time derivative
        tf.u = uprev
        #    if cache.autodiff isa AutoForwardDiff
        #      dT = ForwardDiff.derivative(tf, t)
        #    else
        dT = FiniteDiff.finite_difference_derivative(tf, t, dir = sign(dt))
        #    end

        # Jacobian
        uf.t = t
        if uprev isa AbstractArray
            J = ForwardDiff.jacobian(uf, uprev)
            W = mass_matrix / dtgamma - J
        else
            J = ForwardDiff.derivative(uf, uprev)
            W = 1 / dtgamma - J
        end

        du = f(uprev, p, t)

        linsolve_tmp = du + dtd1 * dT

        k1 = W \ linsolve_tmp
        u = uprev + a21 * k1
        du = f(u, p, t + c2 * dt)

        linsolve_tmp = du + dtd2 * dT + dtC21 * k1

        k2 = W \ linsolve_tmp
        u = uprev + a31 * k1 + a32 * k2
        du = f(u, p, t + c3 * dt)

        linsolve_tmp = du + dtd3 * dT + (dtC31 * k1 + dtC32 * k2)

        k3 = W \ linsolve_tmp
        u = uprev + a41 * k1 + a42 * k2 + a43 * k3
        du = f(u, p, t + c4 * dt)

        linsolve_tmp = du + dtd4 * dT + (dtC41 * k1 + dtC42 * k2 + dtC43 * k3)

        k4 = W \ linsolve_tmp
        u = uprev + a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4
        du = f(u, p, t + c5 * dt)

        linsolve_tmp = du + dtd5 * dT + (dtC52 * k2 + dtC54 * k4 + dtC51 * k1 + dtC53 * k3)

        k5 = W \ linsolve_tmp
        u = uprev + a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5
        du = f(u, p, t + dt)

        linsolve_tmp = du + (dtC61 * k1 + dtC62 * k2 + dtC63 * k3 + dtC64 * k4 + dtC65 * k5)

        k6 = W \ linsolve_tmp
        u = u + k6
        du = f(u, p, t + dt)

        linsolve_tmp = du +
                       (dtC71 * k1 + dtC72 * k2 + dtC73 * k3 + dtC74 * k4 + dtC75 * k5 +
                        dtC76 * k6)

        k7 = W \ linsolve_tmp

        u = u + k7
        du = f(u, p, t + dt)

        linsolve_tmp = du +
                       (dtC81 * k1 + dtC82 * k2 + dtC83 * k3 + dtC84 * k4 + dtC85 * k5 +
                        dtC86 * k6 + dtC87 * k7)

        k8 = W \ linsolve_tmp

        @unpack h21, h22, h23, h24, h25, h26, h27, h28, h31, h32, h33, h34, h35, h36, h37, h38, h41, h42, h43, h44, h45, h46, h47, h48 = cache.tab
        k₁ = h21 * k1 + h22 * k2 + h23 * k3 + h24 * k4 + h25 * k5 + h26 * k6 + h27 * k7 +
             h28 * k8
        k₂ = h31 * k1 + h32 * k2 + h33 * k3 + h34 * k4 + h35 * k5 + h36 * k6 + h37 * k7 +
             h38 * k8
        k₃ = h41 * k1 + h42 * k2 + h43 * k3 + h44 * k4 + h45 * k5 + h46 * k6 + h47 * k7 +
             h48 * k8
        copyat_or_push!(k, 1, k₁)
        copyat_or_push!(k, 2, k₂)
        copyat_or_push!(k, 3, k₃)
    end
    nothing
end

function _ode_addsteps!(k, t, uprev, u, dt, f, p, cache::RosenbrockCache,
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false)
    
    if length(k) < 3 || always_calc_begin
        @unpack dus, tmp, ks, dT, J, W, uf, tf, linsolve_tmp, jac_config, fsalfirst, weight = cache
        @unpack a, C, gamma, d, c = cache.tab

        # Assignments
        sizeu = size(u)
        uidx = eachindex(uprev)
        mass_matrix = f.mass_matrix
        tmp = k8 # integrator.tmp === linsolve_tmp, aliasing fails due to linsolve mutation

        # Precalculations
        dtC = C ./ dt
        dtd = dt .* d
        dtgamma = dt * gamma

        @.. broadcast=false linsolve_tmp = @muladd fsalfirst + dtgamma * dT

        # Jacobian does not need to be re-evaluated after an event
        jacobian2W!(W, mass_matrix, dtgamma, J, true)

        linsolve = cache.linsolve

        # Solve for k1 to k8 using a loop
        for i in 1:8
            if i == 1
                linres = dolinsolve(cache, linsolve; A = W, b = _vec(linsolve_tmp),
                    reltol = cache.reltol)
                vecu = _vec(linres.u)
                @.. broadcast=false ks[i] = -vecu
            else
                tmp = uprev + sum(a[i][j] * ks[j] for j in 1:(i-1))
                f(dus[1], tmp, p, t + c[i] * dt)

                if mass_matrix === I
                    linsolve_tmp = dus[1] + dtd[i] * dT + sum(dtC[i][j] * ks[j] for j in 1:(i-1))
                else
                    dus[1] = sum(dtC[i][j] * ks[j] for j in 1:(i-1))
                    mul!(dus[3], mass_matrix, dus[1])
                    linsolve_tmp = dus[1] + dtd[i] * dT + dus[3]
                end

                linres = dolinsolve(cache, linres.cache; b = _vec(linsolve_tmp),
                    reltol = cache.reltol)
                vecu = _vec(linres.u)
                @.. broadcast=false ks[i] = -vecu
            end
        end

        # Final calculations for k
        @unpack h21, h22, h23, h24, h25, h26, h27, h28, h31, h32, h33, h34, h35, h36, h37, h38, h41, h42, h43, h44, h45, h46, h47, h48 = cache.tab
        
        @.. broadcast=false k6=h21 * k1 + h22 * k2 + h23 * k3 + h24 * k4 + h25 * k5
        copyat_or_push!(k, 1, copy(k6))

        @.. broadcast=false k6=h31 * k1 + h32 * k2 + h33 * k3 + h34 * k4 + h35 * k5
        copyat_or_push!(k, 2, copy(k6))
    end
    nothing
end