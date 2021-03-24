using PyPlot
include("spherical-solver.jl")
include("utilities.jl")
SP = SphericalParticle

function PEPI_potential(
    inner_radius,
    outer_radius,
    Ks,
    ms,
    Kc,
    mc,
    theta0,
    V0s,
    V0c,
)

    solver =
        SP.SphericalSolver(inner_radius, outer_radius, Ks, ms, Kc, mc, theta0)
    return SP.PEPI19_potential(solver, V0s, V0c)
end

function my_potential_difference(
    inner_radius,
    outer_radius,
    Ks,
    ms,
    Kc,
    mc,
    theta0,
    V0s,
    V0c,
)

    solver =
        SP.SphericalSolver(inner_radius, outer_radius, Ks, ms, Kc, mc, theta0)

    cse = SP.core_strain_energy(solver, V0c)
    ccw = SP.core_compression_work(solver, V0c)

    sse = SP.shell_strain_energy(solver, inner_radius, V0s)
    scw = SP.shell_compression_work(solver, inner_radius, V0s)

    return (sse - scw) - (cse - ccw)
end

Ks, Kc = 247e9, 192e9
ms, mc = 126e9, 87e9
theta0 = -0.067
rhos, rhoc = 3.93e3, 3.68e3
V0s = 1.0 / rhos
V0c = 1.0 / rhoc

outer_radius = 1.0
dx = outer_radius / 1e3
inner_radius = dx:dx:outer_radius

pepipd =
    PEPI_potential.(
        inner_radius,
        outer_radius,
        Ks,
        ms,
        Kc,
        mc,
        theta0,
        V0s,
        V0c,
    )
mypd =
    my_potential_difference.(
        inner_radius,
        outer_radius,
        Ks,
        ms,
        Kc,
        mc,
        theta0,
        V0s,
        V0c,
    )

fig, ax = PyPlot.subplots()
ax.plot(inner_radius, pd)
ax.grid()
fig

# solver = SP.SphericalSolver(inner_radius,outer_radius,Ks,ms,Kc,mc,theta0)
#
# p1 = pressure(SP.core_stress(solver))
# shelldevstress = deviatoric_stress(SP.shell_stress(solver,inner_radius))
# srr = shelldevstress[1]
#
# t1 = p1*(V0s-V0c)
# t2 = -0.5*p1^2*(V0s/Ks - V0c/Kc)
# t3 = 0.5*V0s*(1/Ks+3/4ms)*srr^2
