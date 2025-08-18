using Ripserer
using Base.Threads

function boundary_mat(filtration::Ripserer.AbstractFiltration, thresh::Float64)
    n_threads = Threads.nthreads()
    thread_buffer = Vector{Vector{Tuple}}()
    for i in 1:n_threads
        push!(thread_buffer, Tuple[])
    end
    num_vertices = Ripserer.nv(filtration)
    Threads.@threads for i = 1:num_vertices
        Threads.@threads for j = (i + 1):num_vertices
            Threads.@threads for k = (j +1):num_vertices
                # make vertices absolute just in case the simplex is oriented
                sim = abs(Ripserer.simplex(filtration, Val(2), (i, j, k)))
                if (sim != nothing && sim.birth <= thresh)
                    vs = Ripserer.vertices(sim)
                    push!(thread_buffer[Threads.threadid()],Tuple([vs[[1, 2]], vs[[1, 3]], vs[[2, 3]]]))
                end
            end
        end
    end
    return vcat(thread_buffer...)
end

function boundary_mat_d2(filtration_thresh::Ripserer.AbstractFiltration)
    faces = Vector{Tuple}()
    birth_t = Vector()
    edges = Ripserer.edges(filtration_thresh)
    triangles = Ripserer.columns_to_reduce(filtration_thresh, edges)
    for trig in triangles
        sim = abs(trig)
        push!(birth_t, sim.birth)
        vs = Ripserer.vertices(sim)
        push!(faces, vs[[1, 2]])
        push!(faces, vs[[1, 3]])
        push!(faces, vs[[2, 3]])
    end
    return (faces, birth_t)
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
        Threads.@threads for j = (i + 1):num_vertices
            Threads.@threads for k = (j + 1):num_vertices
                sim = abs(Ripserer.simplex(filtration, Val(2), (i, j, k)))
                if (sim != nothing && sim.birth <= thresh_trig)
                    vs = Ripserer.vertices(sim)
                    push!(thread_buffer[Threads.threadid()],Tuple([vs[[1, 2]], vs[[1, 3]], vs[[2, 3]]]))
                end
            end
        end
    end
    Threads.@threads for rep in rep_cycle
        if (rep.birth < thresh_hole_birth && (rep.death - rep.birth) < thresh_hole_life)
            push!(thread_buffer[Threads.threadid()],Tuple([Ripserer.vertices(v) for v in rep.representative]))
        end
    end
    return vcat(thread_buffer...)
end

function recontstruct_multiple_cycles(
    filtration::Ripserer.AbstractFiltration,
    interval
)
    
end
