"""
    principal_angle(x)

principal value of angle `x` in [-π, π]"""
function principal_angle(x)
    return rem2pi(x, RoundNearest)
end

"""
    angle_in(x, x1, x2)

`true` if angle `x` is in the angle interval `[x1, x2]`

Angles `x`, `x1` and `x2` are supposed to be in [-π, π]
(use [`principal_angle`](@ref) if necessary).

Interval [x1, x2] is assumed to represent a _positive_ angular sector between 0 and 2π.
If `x1 <= x2`, then this function tests if `x ∈ [x1, x2]`.
However, since angles are defined module 2π, it may be that `x1 > x2`,
then this function tests if `x ∈ [-π, x2]` or `x ∈ [x1, +π]`,
expect that the -π and +π bounds are omitted from the test.
"""
function angle_in(x, x1, x2)
    if x1 <= x2 # easy case
        return x >= x1 && x <= x2
    else
        return x >= x1 || x <= x2
    end
end

"""
    solve_cossin(a, b, c)

Compute (x1, x2) the two principal solutions of `a`.cos(`x`) + `b`.sin(`x`) = `c`,
so that a.cos(x) + b.sin(x) ≥ c on [x1, x2]

Notice that since angles are defined module 2π, it is not guaranted that x1 ≤ x2

If the equation has no solution, there are two cases:
- if a.cos(x) + b.sin(x) ≥ c ∀ x, returns an interval of length 2π
- if a.cos(x) + b.sin(x) < c ∀ x, returns an interval of length 0
"""
function solve_cossin(a, b, c)
    # Special trivial case
    if a == 0 && b == 0
        if c > 0 # a.cos(x) + b.sin(x) < c ∀ x
            return (0., 0.)
        else
            return (-1π, 1π)
        end
    end

    # Normalization
    r = sqrt(a^2 + b^2)
    an = a/r
    bn = b/r
    cn = c/r
    
    # Transform to cos(x-x0) = cn
    # with x0 s.t. cos(x0) = an and sin(x0) = bn
    x0 = atan(bn, an)

    # Special cases:
    if cn > 1 # a.cos(x) + b.sin(x) < c ∀ x
        x1 = principal_angle(x0)
        x2 = x1
    elseif cn < -1 # a.cos(x) + b.sin(x) > c ∀ x
        x1 = -1π
        x2 = 1π
    else # If |cn| < 1, transform to cos(x-x0) = cos(Δx)
        Δx = acos(cn)
        x1 = principal_angle(x0 - Δx)
        x2 = principal_angle(x0 + Δx)
    end
    return (x1, x2)
end

solve_cossin(2., 0, 2*cos(0.1))
solve_cossin(2., 0, -10)