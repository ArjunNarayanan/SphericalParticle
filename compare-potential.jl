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

Ks, Kc = 247., 192.
ms, mc = 126., 87.
theta0 = -0.067
rhos, rhoc = 3.93e3, 3.68e3
V0s = 1.0 / rhos
V0c = 1.0 / rhoc

ΔG0Jmol = -14351.0
molarmass = 0.147
ΔG0 = ΔG0Jmol/molarmass

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

err = maximum(abs.(mypd - pepipd))

pd = (mypd*1e9 .+ ΔG0)/abs(ΔG0)

fig, ax = PyPlot.subplots()
ax.plot(inner_radius, pd)
ax.grid()
fig
