"""
    function calc_rate(su, rel_reserve::Tuple, turnover::Tuple, 
                       j_E_mai, y_E_Ea, y_E_Eb, y_V_E, κsoma, tstep)

Calculate growth rate using a numeric root-finder, also determining
wether the organ is alive or dead.

Returns a `Tuple` holding the rate and a Bool for alive/dead status
"""
function calc_rate end
function calc_rate(su, rel_reserve::Tuple, turnover::Tuple,
                   j_E_mai, y_E_Ea, y_E_Eb, y_V_E, κsoma, tstep)
    one_r = oneunit(eltype(turnover))
    f(r) = rate_formula(r, su, rel_reserve, turnover, j_E_mai, y_E_Ea, y_E_Eb, y_V_E, κsoma)

    r, success = _secant_rate(f, one_r)
    if success == ALIVE && r >= zero(r)
        return r, ALIVE
    end

    r, success = _bracketed_rate(f, turnover, one_r)
    if success == ALIVE && r >= zero(r)
        return r, ALIVE
    end

    @warn "Root for rate not found at t = $tstep"
    return zero(one_r), DEAD
end

_secant_rate(f, one_r) = begin
    atol = one_r * 1e-10
    maxiters = 200
    try
        findzero(Secant(-2one_r, one_r), atol, maxiters) do r
            f(r)
        end
    catch
        return zero(one_r), DEAD
    end
end

function _bracketed_rate(f, turnover, one_r)
    max_turn = _max_turnover(turnover, one_r)
    lo = zero(max_turn)
    flo = f(lo)
    _is_zero(flo) && return lo, ALIVE

    hi = max_turn
    step = max_turn
    for _ in 1:12
        fhi = f(hi)
        if _is_zero(fhi)
            return hi, ALIVE
        elseif _has_sign_change(flo, fhi)
            return _bisection(f, lo, hi, one_r * 1e-10, 200)
        end
        hi += step
        step *= 2
    end

    return zero(one_r), DEAD
end

function _max_turnover(turnover, one_r)
    max_turn = zero(one_r)
    for val in turnover
        aval = abs(val)
        max_turn = max(max_turn, aval)
    end
    max_turn == zero(max_turn) && return one_r
    return max_turn
end

_is_zero(x) = iszero(Unitful.ustrip(x))

_has_sign_change(a, b) = begin
    sa = Base.signbit(Unitful.ustrip(a))
    sb = Base.signbit(Unitful.ustrip(b))
    sa != sb
end

function _bisection(f, lo, hi, tol, maxiters)
    flo = f(lo)
    fhi = f(hi)
    mid = (lo + hi) / 2
    for _ in 1:maxiters
        mid = (lo + hi) / 2
        fmid = f(mid)
        if _is_zero(fmid) || abs(hi - lo) <= tol
            return mid, ALIVE
        end
        if _has_sign_change(flo, fmid)
            hi, fhi = mid, fmid
        else
            lo, flo = mid, fmid
        end
    end
    return mid, DEAD
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


