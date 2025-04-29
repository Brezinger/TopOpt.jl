using Makie, CairoMakie, TopOpt, LinearAlgebra, StatsFuns
# using GLMakie


# Parameter setup
E = 1.0 # Young’s modulus
v = 0.3 # Poisson’s ratio
f = 1.0 # downward force
rmin = 4.0 # filter radius
xmin = 0.0001 # minimum density
problem_size = (30, 10, 10)
V = 0.5 # maximum volume fraction
p = 4.0 # penalty

x0 = fill(V, prod(problem_size)); # initial design


# Define new problem
using Ferrite
using TopOpt
using TopOpt.TopOptProblems: RectilinearGrid, Metadata
using TopOpt.TopOptProblems: left, right, bottom, middley, middlez,
    nnodespercell, nfacespercell, find_black_and_white, find_varind
using TopOpt.Utilities: @params

@params struct NewPointLoadCantilever{dim, T, N, M} <: StiffnessTopOptProblem{dim, T}
    rect_grid::RectilinearGrid{dim, T, N, M}
    E::T
    ν::T
    ch::ConstraintHandler{<:DofHandler{dim, <:Cell{dim,N,M}, T}, T}
    load_dict::Dict{Int, Vector{T}}
    black::AbstractVector
    white::AbstractVector
    varind::AbstractVector{Int}
    metadata::Metadata
end

function NewPointLoadCantilever(::Type{Val{CellType}}, nels::NTuple{dim,Int}, sizes::NTuple{dim}, 
    E = 1.0, ν = 0.3, force = 1.0) where {dim, CellType}
    iseven(nels[2]) && (length(nels) < 3 || iseven(nels[3])) || throw("Grid does not have an even number of elements along the y and/or z axes.")

    _T = promote_type(eltype(sizes), typeof(E), typeof(ν), typeof(force))
    if _T <: Integer
        T = Float64
    else
        T = _T
    end
    if CellType === :Linear || dim === 3
        rect_grid = RectilinearGrid(Val{:Linear}, nels, T.(sizes))
    else
        rect_grid = RectilinearGrid(Val{:Quadratic}, nels, T.(sizes))
    end

    if haskey(rect_grid.grid.facesets, "fixed_all") 
        pop!(rect_grid.grid.facesets, "fixed_all")
    end
    addnodeset!(rect_grid.grid, "fixed_all", x -> left(rect_grid, x));
    
    if haskey(rect_grid.grid.nodesets, "down_force") 
        pop!(rect_grid.grid.nodesets, "down_force")
    end
    if dim == 3
        addnodeset!(rect_grid.grid, "down_force", x -> right(rect_grid, x) && 
            bottom(rect_grid, x));
            #  && middlez(rect_grid, x));
    else
        addnodeset!(rect_grid.grid, "down_force", x -> right(rect_grid, x) && 
            right(rect_grid, x) && middley(rect_grid, x));
    end

    # Create displacement field u
    dh = DofHandler(rect_grid.grid)
    if CellType === :Linear || dim === 3
        push!(dh, :u, dim) # Add a displacement field
    else
        ip = Lagrange{2, RefCube, 2}()
        push!(dh, :u, dim, ip) # Add a displacement field        
    end
    close!(dh)
    
    ch = ConstraintHandler(dh)

    dbc = Dirichlet(:u, getnodeset(rect_grid.grid, "fixed_all"), (x,t) -> zeros(T, dim), collect(1:dim))
    add!(ch, dbc)
    close!(ch)
    t = T(0)
    Ferrite.update!(ch, t)

    metadata = Metadata(dh)
    load_dict = Dict{Int, Vector{T}}()
    for fnode in getnodeset(rect_grid.grid, "down_force")
    	load_dict[fnode] = [0, -force, 0]
    end

    N = nnodespercell(rect_grid)
    M = nfacespercell(rect_grid)

    black, white = find_black_and_white(dh)
    varind = find_varind(black, white)
    
    return NewPointLoadCantilever(rect_grid, E, ν, ch, load_dict, black, white, varind, metadata)
end

# used in FEA to determine default quad order
# we don't assume the problem struct has `rect_grid` to define its grid
TopOptProblems.nnodespercell(p::NewPointLoadCantilever) = nnodespercell(p.rect_grid)

# ! important, used for specification!
function TopOptProblems.getcloaddict(p::NewPointLoadCantilever{dim, T}) where {dim, T}
    # f = T[0, -p.force, 0]
    # fnode = Tuple(getnodeset(p.rect_grid.grid, "down_force"))[1]
    # return Dict{Int, Vector{T}}(fnode => f)
    return p.load_dict
end


problem = NewPointLoadCantilever(Val{:Linear}, problem_size, (1.0, 1.0, 1.0), E, v, f);

solver = FEASolver(Direct, problem, xmin = xmin);


solver()
u0 = solver.u


# visualize the deformation!
fig = visualize(problem, u=u0, default_exagg_scale=0.001)
Makie.display(fig)

# cg_assembly_solver = FEASolver(CG, Assembly, problem, xmin = xmin)
# cg_matrix_free_solver = FEASolver(CG, MatrixFree, problem, xmin = xmin)

cheqfilter = DensityFilter(solver, rmin = rmin)
stress = TopOpt.von_mises_stress_function(solver)
comp = TopOpt.Compliance(solver)

# minimize compliance
function obj(x)
    return comp(cheqfilter(PseudoDensities(x)))
end

# volume bound
function constr(x)
    return sum(cheqfilter(PseudoDensities(x))) / length(x) - V
end;


# Run optimization
m = Model(obj)
addvar!(m, zeros(length(x0)), ones(length(x0)))
Nonconvex.add_ineq_constraint!(m, constr)

options = MMAOptions(
    maxiter=300, tol = Tolerance(kkt = 1e-3, x=1e-3, f = 1e-3),
)
TopOpt.setpenalty!(solver, p)

@time r = Nonconvex.optimize(
    m, MMA87(dualoptimizer = ConjugateGradient()),
    x0, options = options,
);


# show the optimized results!
@show obj(r.minimizer)
@show constr(r.minimizer)
@show maximum(stress(cheqfilter(r.minimizer)))

topology = cheqfilter(PseudoDensities(r.minimizer)).x;
fig = visualize(problem; topology = topology)
Makie.display(fig)

save("result_3D.png", fig)
