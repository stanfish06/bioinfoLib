using Ripserer
using Base.Threads

function boundary_mat(filtration::Ripserer.AbstractFiltration, thresh::Float64)
    n_threads = Threads.nthreads()
    thread_buffer = Vector{Vector{Tuple{Int, Int, Int}}}()
    for i in 1:n_threads
        push!(thread_buffer, Tuple{Int, Int, Int}[])
    end
    num_vertices = Ripserer.nv(filtration)
    Threads.@threads for i = 1:num_vertices
        Threads.@threads for j = (i + 1):num_vertices
            Threads.@threads for k = (j +1):num_vertices
                sim = Ripserer.simplex(filtration, Val(2), (i, j, k))
                if (sim != nothing && sim.birth <= thresh)
                    push!(thread_buffer[Threads.threadid()],Ripserer.vertices(sim))
                end
            end
        end
    end
    return vcat(thread_buffer...)
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
                sim = Ripserer.simplex(filtration, Val(2), (i, j, k))
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
