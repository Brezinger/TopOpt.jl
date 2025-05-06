using TopOpt
using TopOpt.TopOptProblems: InpStiffness  # import the high-level INP reader
using TopOpt.TopOptProblems, TopOpt.TopOptProblems.InputOutput.INP.Parser
using Makie
using CairoMakie
# using GLMakie

using TimerOutputs

println("Start running.")
# https://github.com/KristofferC/TimerOutputs.jl
to = TimerOutput()
reset_timer!(to)

# Define the problem
E = 1.0 # Young’s modulus
v = 0.3 # Poisson’s ratio
f = 1.0; # downward force

# Parameter settings
V = 0.3 # volume fraction
# xmin = 0.001 # minimum density
xmin = 1e-6 # minimum density
rmin = 10.0 # density filter radius

filepath = joinpath(@__DIR__, "hive_blocked_RevB_Tet.inp")
content = extract_inp(filepath)
problem = InpStiffness(content; keep_load_cells=false)

# get black cell ids
black_elids = content.cellsets["Black"]
#cell_elids = problem.metadata.cells[:,1]
#black_cells = findall(elid -> elid in black_elids, cell_elids)

# build a mask of length n_cells
ncell = length(problem.black)
blackmask = falses(ncell)
blackmask[black_elids] .= true

# Define a finite element solver
@timeit to "penalty def" penalty = TopOpt.PowerPenalty(3.0)
@timeit to "solver def" begin
    solver = FEASolver(Direct, problem; xmin=xmin, penalty=penalty)
    filter = DensityFilter(solver, rmin=rmin)
end

# Define compliance objective
@timeit to "objective def" begin
    # Define compliance objective
    comp = Compliance(solver)
    obj = x -> comp(filter(PseudoDensities(x)))
end

# Define volume constraint
@timeit to "constraint def" begin
    volfrac = TopOpt.Volume(solver)
    constr = x -> volfrac(filter(PseudoDensities(x))) - V
end

@timeit to "define problem" begin
    x0 = fill(V, length(solver.vars))
    nvar = length(solver.vars)
    lb = zeros(nvar); ub = ones(nvar)
    # set the bounds of the black cells to 1.0 and set the initial guess to 1.0
    for i in black_elids
        lb[i] = 1. - 1.e-9; ub[i] = 1.
        x0[i] = 1.
    end    
    model = Model(obj)
    addvar!(model, lb, ub)

    add_ineq_constraint!(model, constr) 
    alg = MMA87()
    tol = Tolerance(x=1e-3, f=1e-6, fabs=1e-3, frel=0.0, kkt=1e-3, infeas=1e-3)
    convcriteria = GenericCriteria()
    options = MMAOptions(;
        maxiter=1000, tol=tol, convcriteria=convcriteria
    )
end

@timeit to "simp run" r = optimize(model, alg, x0; options)

# Print the timings in the default way
println()
show(to)

@show obj(r.minimizer)
@show constr(r.minimizer)

# Visualize the result using Makie.jl
fig = visualize(
    problem;
    topology=r.minimizer,
    default_exagg_scale=0.07,
    scale_range=10.0,
    display_supports=false,
    vector_linewidth=3,
    vector_arrowsize=0.005,
    default_support_scale=0.01,
    default_load_scale=0.01,
)
Makie.display(fig)

Makie.save("hive.png", fig)

TopOpt.save_mesh("hive.vtu", problem, r.minimizer)