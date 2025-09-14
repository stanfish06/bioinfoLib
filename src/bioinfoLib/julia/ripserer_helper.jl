using Ripserer
using Base.Threads
using Graphs
using SimpleWeightedGraphs
using Random

function boundary_mat(filtration::Ripserer.AbstractFiltration, thresh::Float64)
    n_threads = Threads.nthreads()
    thread_buffer = Vector{Vector{Tuple}}()
    for i in 1:n_threads
        push!(thread_buffer, Tuple[])
    end
    num_vertices = Ripserer.nv(filtration)
    Threads.@threads for i = 1:num_vertices
        Threads.@threads for j = (i+1):num_vertices
            Threads.@threads for k = (j+1):num_vertices
                # make vertices absolute just in case the simplex is oriented
                sim = abs(Ripserer.simplex(filtration, Val(2), (i, j, k)))
                if (sim != nothing && sim.birth <= thresh)
                    vs = Ripserer.vertices(sim)
                    push!(thread_buffer[Threads.threadid()], Tuple([vs[[1, 2]], vs[[1, 3]], vs[[2, 3]]]))
                end
            end
        end
    end
    return vcat(thread_buffer...)
end

function boundary_mat_d2(filtration_thresh::Ripserer.AbstractFiltration)
    edges_pool = Vector{Tuple}()
    trigs_pool = Vector{Tuple}()
    birth_t = Vector()
    edges_full = Ripserer.edges(filtration_thresh)
    triangles = Ripserer.columns_to_reduce(filtration_thresh, edges_full)
    for trig in triangles
        sim = abs(trig)
        vs = Ripserer.vertices(sim)
        push!(edges_pool, vs[[1, 2]]) 
        push!(edges_pool, vs[[1, 3]])
        push!(edges_pool, vs[[2, 3]])
        push!(trigs_pool, vs)
        push!(birth_t, sim.birth)
    end
    return (edges_pool, trigs_pool, birth_t)
end

# If you can get the column reduced matrix, then you can use fewer number of columns
# Get H1 and extracts its representatives
# birth time of each loop representative is the birth time of the birth edge
# so all other edges of the representative should be present
function reduced_boundary_mat_d2(filtration_thresh::Ripserer.AbstractFiltration)
    verts = Vector{Tuple}()
    e_rep_idx = Vector()
    birth_t = Vector()
    ripser = Ripserer.ripserer(filtration_thresh, reps=1, alg=:involuted)
    loops = ripser[2]
    i = 0
    for loop in loops
        es = Ripserer.representative(loop)
        for e in es
            vs = vertices(e)
            push!(verts, vs)
            push!(birth_t, birth(loop))
            push!(e_rep_idx, i)
        end
        i += 1
    end
    return (verts, birth_t, e_rep_idx)
end

function get_trigs(filtration_thresh::Ripserer.AbstractFiltration)
    verts = Vector{Tuple}()
    birth_t = Vector()
    edges = Ripserer.edges(filtration_thresh)
    triangles = Ripserer.columns_to_reduce(filtration_thresh, edges)
    for trig in triangles
        sim = abs(trig)
        push!(birth_t, sim.birth)
        vs = Ripserer.vertices(sim)
        push!(verts, vs)
    end
    return (verts, birth_t)
end

function inter_arrival_filtration(filtration_thresh::Ripserer.AbstractFiltration)
    birth = Vector()
    edges = Ripserer.edges(filtration_thresh)
    triangles = Ripserer.columns_to_reduce(filtration_thresh, edges)
    for trig in triangles
        sim = abs(trig)
        push!(birth, sim.birth)
    end
    return birth
end

function boundary_mat_fill_hole(filtration::Ripserer.AbstractFiltration, rep_cycle, thresh_trig::Float64, thresh_hole_birth::Float64, thresh_hole_life::Float64)
    n_threads = Threads.nthreads()
    thread_buffer = Vector{Vector{Tuple}}()
    for i in 1:n_threads
        push!(thread_buffer, Tuple[])
    end
    num_vertices = Ripserer.nv(filtration)
    Threads.@threads for i = 1:num_vertices
        Threads.@threads for j = (i+1):num_vertices
            Threads.@threads for k = (j+1):num_vertices
                sim = abs(Ripserer.simplex(filtration, Val(2), (i, j, k)))
                if (sim != nothing && sim.birth <= thresh_trig)
                    vs = Ripserer.vertices(sim)
                    push!(thread_buffer[Threads.threadid()], Tuple([vs[[1, 2]], vs[[1, 3]], vs[[2, 3]]]))
                end
            end
        end
    end
    Threads.@threads for rep in rep_cycle
        if (rep.birth < thresh_hole_birth && (rep.death - rep.birth) < thresh_hole_life)
            push!(thread_buffer[Threads.threadid()], Tuple([Ripserer.vertices(v) for v in rep.representative]))
        end
    end
    return vcat(thread_buffer...)
end

# Use Yen's k shortest path to get n best cycles
function reconstruct_n_loop_representatives(
    cocycles,
    filt,
    rep_idx,
    n,
    life_pct=0.1,
    n_force_deviate=4,
    n_reps_per_loop=8,
    loop_lower_pct=5,
    loop_upper_pct=95,
    n_max_cocycles=10
)
    rep = cocycles[2][length(cocycles[2])-rep_idx]
    filt_t = birth(rep) + (death(rep) - birth(rep)) * life_pct
    # get all cocycle representatives
    cocycles_filt = filter!(simplex.(representative(rep))) do sx
        birth(sx) <= filt_t
    end
    # get all existing edges
    edges_filt = filter!(edges(filt)) do sx
        birth(sx) <= filt_t
    end
    # create weighted graph and disconnect the cocycles by setting them to infinite weights
    sources = vcat(
        getindex.(vertices.(edges_filt), 1),
        getindex.(vertices.(cocycles_filt), 1)
    )
    destinations = vcat(
        getindex.(vertices.(edges_filt), 2),
        getindex.(vertices.(cocycles_filt), 2)
    )
    weights = vcat(
        birth.(edges_filt),
        birth.(cocycles_filt) * Inf
    )
    cycles_pool = Vector{Vector{Int64}}()
    cycles_dist = Vector{Float64}()
    for _ in 1:n_force_deviate
        cycles_pool_tmp = Vector{Vector{Int64}}()
        g = SimpleWeightedGraph(sources, destinations, weights; combine=max)
        n_cocycles_used = 0
        for (i, j) in Iterators.map(vertices, cocycles_filt)
            if n_cocycles_used == n_max_cocycles
                break
            end
            res = yen_k_shortest_paths(g, i, j, g.weights, n_reps_per_loop)
            for path in res.paths
                push!(cycles_pool, path)
                push!(cycles_pool_tmp, path)
            end
            append!(cycles_dist, res.dists)
            n_cocycles_used = n_cocycles_used + 1
        end
        # force deviation from the previous cycles by increasing the edge weights to infinity
        for path in cycles_pool_tmp
            append!(sources, path[1:end-1])
            append!(destinations, path[2:end])
            append!(weights, fill(Inf, length(path) - 1))
        end
    end
    cycles_pool = sort(collect(zip(cycles_dist, cycles_pool)), by=first)
    n_cycles_total = length(cycles_pool)
    n_cycles_return = min(n_cycles_total, n)
    step = (loop_upper_pct - loop_lower_pct) / (n_cycles_return - 1)
    cycles_idx_pick = Vector{Int64}()
    for i in 1:n_cycles_return
        dist_pct = (loop_lower_pct + step * (i - 1)) / 100
        push!(cycles_idx_pick, min(floor(n_cycles_total * dist_pct) + 1, n_cycles_total))
    end
    top_cycles = last.(cycles_pool[cycles_idx_pick])
    top_cycles_dist = first.(cycles_pool[cycles_idx_pick])
    mask = top_cycles_dist .< Inf
    return (top_cycles[mask], top_cycles_dist[mask])
end

function noisy_circle(n, r=1, noise=0.1)
    points = NTuple{2,Float64}[]
    for _ in 1:n
        θ = 2π * rand()
        push!(points, (r * sin(θ) + noise * rand(), r * cos(θ) + noise * rand()))
    end
    return points
end
