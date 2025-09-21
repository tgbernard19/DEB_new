"""
    function calc_rate(su, rel_reserve::Tuple, turnover::Tuple, 
                       j_E_mai, y_E_Ea, y_E_Eb, y_V_E, κsoma, tstep)

Calculate growth rate using a numeric root-finder, also determining
wether the organ is alive or dead.

Returns a `Tuple` holding the rate and a Bool for alive/dead status
"""
function calc_rate end

function _secant_root(f, x0, x1; atol, maxiters)
    fx0 = f(x0)
    fx1 = f(x1)
    xprev, fprev = x0, fx0
    xcurr, fcurr = x1, fx1

    for _ in 1:maxiters
        denom = fcurr - fprev
        if iszero(denom)
            return xcurr, DEAD
        end

        xnext = xcurr - fcurr * (xcurr - xprev) / denom
        if abs(xnext - xcurr) <= atol
            return xnext, ALIVE
        end

        xprev, fprev = xcurr, fcurr
        xcurr = xnext
        fcurr = f(xcurr)
    end

    return xcurr, DEAD
end

function _bisection_root(f, lo, hi; atol, maxiters)
    flo = f(lo)
    fhi = f(hi)

    if iszero(flo)
        return lo, ALIVE
    elseif iszero(fhi)
        return hi, ALIVE
    end

    if signbit(flo) == signbit(fhi)
        return hi, DEAD
    end

    local mid, fmid
    for _ in 1:maxiters
        mid = (lo + hi) / 2
        if abs(hi - lo) <= atol
            return mid, ALIVE
        end

        fmid = f(mid)
        if iszero(fmid)
            return mid, ALIVE
        end

        if signbit(fmid) == signbit(flo)
            lo, flo = mid, fmid
        else
            hi, fhi = mid, fmid
        end
    end

    return mid, DEAD
end

function calc_rate(su, rel_reserve::Tuple, turnover::Tuple,
                   j_E_mai, y_E_Ea, y_E_Eb, y_V_E, κsoma, tstep)
    # Find the type of `rate` so that diff and units work with atol etc
    one_r = oneunit(eltype(turnover))
    # Get rate with a zero finder
    solver = let turnover=turnover, j_E_mai=j_E_mai, y_E_Ea=y_E_Ea, y_E_Eb=y_E_Eb, y_V_E=y_V_E, κsoma=κsoma
        atol = one_r * 1e-10
        maxiters = 200
        f(r) = rate_formula(r, su, rel_reserve, turnover, j_E_mai, y_E_Ea, y_E_Eb, y_V_E, κsoma)

        attempts = (
            (-2one_r, 1one_r),
            (zero(one_r), max(maximum(turnover), one_r))
        )

        for (x0, x1) in attempts
            r, state = _secant_root(f, x0, x1; atol=atol, maxiters=maxiters)
            if state == ALIVE && r >= zero(r)
                return r, ALIVE
            end
        end

        lo = zero(one_r)
        hi = maximum(turnover)
        if hi > lo
            r, state = _bisection_root(f, lo, hi; atol=atol, maxiters=maxiters)
            if state == ALIVE && r >= zero(r)
                return r, ALIVE
            end
        end

        zero(one_r), DEAD
    end

    r, state = solver

    if state == DEAD
        @warn "Root for rate not found at t = $tstep"
        return r, DEAD
    elseif r < zero(r)
        @warn "Rate is less than zero at t = $tstep"
        return zero(r), DEAD
    else
        return r, ALIVE
    end
end

"""
    rate_formula(r, ureserve::Tuple, turnover::Tuple, j_E_mai, y_V_E, κsoma)

Rate formulas for E, CN or CNE reserves
"""
function rate_formula end
rate_formula(r, su, rel_reserve::NTuple{1}, turnover::NTuple{1},
             j_E_mai, y_E_Ea, y_E_Eb, y_V_E, κsoma) = begin
    (j_E,) = rel_reserve .* (turnover .- r)
    y_V_E * (κsoma * j_E - j_E_mai) - r
end
rate_formula(r, su, rel_reserve::NTuple{2}, turnover::NTuple{2},
             j_E_mai, y_E_Ea, y_E_Eb, y_V_E, κsoma) = begin
    j_Ea, j_Eb = rel_reserve .* (turnover .- r)
    j_E = synthesizing_unit(su, j_Ea * y_E_Ea, j_Eb * y_E_Eb)
    y_V_E * (κsoma * j_E - j_E_mai) - r
end
rate_formula(r, su, rel_reserve::NTuple{3}, turnover::NTuple{3},
             j_E_mai, y_E_Ea, y_E_Eb, y_V_E, κsoma) = begin
    j_Ea, j_Eb, j_E = rel_reserve .* (turnover .- r)
    j_E += synthesizing_unit(su, j_Ea * y_E_Ea, j_Eb * y_E_Eb)
    y_V_E * (κsoma * j_E - j_E_mai) - r
end


