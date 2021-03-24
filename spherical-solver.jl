module SphericalParticle

include("moduli-conversion.jl")
include("utilities.jl")

function analytical_coefficient_matrix(
    inner_radius,
    outer_radius,
    Ks,
    ms,
    Kc,
    mc,
)
    row1 = [inner_radius, -inner_radius, -1 / inner_radius^2]
    row2 = [3Kc, -3Ks, 4ms / inner_radius^3]
    row3 = [0, 3Ks, -4ms / outer_radius^3]

    return vcat(row1', row2', row3')
end

function analytical_coefficient_rhs(Ks, theta0)
    return [0, -Ks * theta0, Ks * theta0]
end

struct SphericalSolver
    inradius::Any
    outradius::Any
    A1c::Any
    A1s::Any
    A2s::Any
    Ks::Any
    ms::Any
    Kc::Any
    mc::Any
    theta0::Any
    function SphericalSolver(inradius, outradius, Ks, ms, Kc, mc, theta0)
        op = analytical_coefficient_matrix(inradius, outradius, Ks, ms, Kc, mc)
        R = analytical_coefficient_rhs(Ks, theta0)
        A1c, A1s, A2s = op \ R
        new(inradius, outradius, A1c, A1s, A2s, Ks, ms, Kc, mc, theta0)
    end
end

function core_strain(solver::SphericalSolver)
    return [solver.A1c, solver.A1c, solver.A1c]
end

function core_stress(solver::SphericalSolver)
    Kc = solver.Kc
    return 3Kc * solver.A1c * [1.0, 1.0, 1.0]
end

function core_strain_energy(solver, V0c)
    strain = core_strain(solver)
    stress = core_stress(solver)
    return 0.5 * V0c * sum(stress .* strain)
end

function core_compression_work(solver,V0c)
    strain = core_strain(solver)
    V = V0c*(1+sum(strain))
    srr = core_stress(solver)[1]
    return V*srr
end

function shell_strain(solver::SphericalSolver, r)
    t0 = solver.theta0 / 3
    err = solver.A1s - 2 * solver.A2s / r^3 - t0
    ett = solver.A1s + solver.A2s / r^3 - t0
    return [err, ett, ett]
end

function shell_stress(solver::SphericalSolver, r)
    K = solver.Ks
    m = solver.ms
    A1s = solver.A1s
    A2s = solver.A2s
    t0 = solver.theta0

    srr = 3K * A1s - 4m * A2s / r^3 - K * t0
    stt = 3K * A1s + 2m * A2s / r^3 - K * t0
    return [srr, stt, stt]
end

function shell_strain_energy(solver, r, V0s)
    strain = shell_strain(solver, r)
    stress = shell_stress(solver, r)
    return 0.5 * V0s * sum(strain .* stress)
end

function shell_compression_work(solver,r,V0s)
    strain = shell_strain(solver,r)
    srr = shell_stress(solver,r)[1]
    V = V0s*(1+sum(strain))
    return srr*V
end

function PEPI19_potential(solver::SphericalSolver, V0s, V0c)
    corestress = core_stress(solver)
    shellstress = shell_stress(solver, solver.inradius)

    p1 = pressure(corestress)
    shelldevstress = deviatoric_stress(shellstress)
    srr = shelldevstress[1]

    Ks, Kc = solver.Ks, solver.Kc
    ms = solver.ms

    return p1 * (V0s - V0c) - 0.5 * p1^2 * (V0s / Ks - V0c / Kc) +
           0.5 * V0s * (1 / Ks + 3 / 4ms) * srr^2

end

end
