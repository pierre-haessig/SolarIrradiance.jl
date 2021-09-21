using SolarIrradiance
using Test

@testset "angle helpers" begin
    principal_angle = SolarIrradiance.principal_angle
    @test principal_angle(0.0) == 0.0
    @test principal_angle(1π) == 1π
    @test principal_angle(-1π) == -1π

    @test principal_angle(1π+2eps(Float64)) ≈ -1π atol=2eps(Float64)
    @test principal_angle(-1π-2eps(Float64)) ≈ 1π atol=2eps(Float64)

    @test principal_angle(2pi) ≈ 0.0 atol=2eps(Float64)

    angle_in = SolarIrradiance.angle_in
    # Normal angle interval [x1, x2] with x1 <= x2
    @test angle_in(0., -1., 1.)
    @test angle_in(-1., -1., 1.) # x=x1
    @test angle_in(+1., -1., 1.) # x=x2
    @test !angle_in(-1.1, -1., 1.)
    @test !angle_in(+1.1, -1., 1.)

    # Degenerate angle interval [x1, x2] with x1 == x2
    @test angle_in(0.2, 0.2, 0.2)
    @test !angle_in(1.1, 0.2, 0.2)
    @test !angle_in(-1.1, 0.2, 0.2)

    # Flipped values of angle interval: x1 > x2
    @test !angle_in(0., 1., -1.)
    @test angle_in(+1., 1., -1.) # x=x1
    @test angle_in(-1., 1., -1.) # x=x2
    @test angle_in(1.1, 1., -1.)
    @test angle_in(-1.1, 1., -1.)
end


@testset "a.cos(x) + b.sin(x) = c" begin
    solve_cossin = SolarIrradiance.solve_cossin
    atol = 10eps(Float64)
    # Let y>0:
    y = 0.1
    # 2.cos(x) = 2.cos(y) : output is (-y, +y)
    # (for small y only, due to output expressed in [-π,π])
    @test solve_cossin(2., 0, 2*cos(y))[1] ≈ -y atol=atol
    @test solve_cossin(2., 0, 2*cos(y))[2] ≈ +y atol=atol

    # 2.sin(x) = 2.cos(y) : output is (-y + π/2, +y + π/2)
    # (for small y only, due to output expressed in [-π,π])
    @test solve_cossin(0., 2., 2*cos(y))[1] ≈ -y + π/2 atol=atol
    @test solve_cossin(0., 2., 2*cos(y))[2] ≈ +y + π/2 atol=atol

    # 2.cos(x) = -2.cos(y) : output is (+y - π, -y + π)
    # (output seems flipped because it is expressed in [-π,π],
    #  but it is in fact a large angle interval)
    @test solve_cossin(2., 0, -2*cos(y))[1] ≈ +y - π atol=atol
    @test solve_cossin(2., 0, -2*cos(y))[2] ≈ -y + π atol=atol
end