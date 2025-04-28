from juliacall import Main as julia


def ripserer_helper():
    julia.seval("using Ripserer")
    julia.seval("using Base.Threads")
    julia.seval(
        " function boundary_mat(filtration::Ripserer.AbstractFiltration, dim::Integer, thresh::Float64)"
        "   n_threads = Threads.nthreads();"
        "   thread_buffer = Vector{Vector{Tuple{Int, Int, Int}}}();"
        "   for i in 1:n_threads"
        "       push!(thread_buffer, Tuple{Int, Int, Int}[])"
        "   end;"
        "   num_vertices = Ripserer.nv(filtration);"
        "   Threads.@threads for i = 1:num_vertices"
        "       Threads.@threads for j = (i + 1):num_vertices"
        "           Threads.@threads for k = (j +1):num_vertices"
        "               sim = Ripserer.simplex(filtration, Val(2), (i, j, k));"
        "               if (sim != nothing && sim.birth <= thresh)"
        "                   push!(thread_buffer[Threads.threadid()],Ripserer.vertices(sim));"
        "               end;"
        "           end;"
        "       end;"
        "   end;"
        "   return vcat(thread_buffer...)"
        " end"
    )
