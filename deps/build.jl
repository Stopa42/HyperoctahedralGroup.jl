using Permutations, StaticArrays, FileIO

println("Building of HyperoctahedralGroup.jl started.")

const xyp = SA[0 -1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1]
const xym = SA[0 1 0 0; -1 0 0 0; 0 0 1 0; 0 0 0 1]

const xzp = SA[0 0 -1 0; 0 1 0 0; 1 0 0 0; 0 0 0 1]
const xzm = SA[0 0 1 0; 0 1 0 0; -1 0 0 0; 0 0 0 1]

const xwp = SA[0 0 0 -1; 0 1 0 0; 0 0 1 0; 1 0 0 0]
const xwm = SA[0 0 0 1; 0 1 0 0; 0 0 1 0; -1 0 0 0]

const yzp = SA[1 0 0 0; 0 0 -1 0; 0 1 0 0; 0 0 0 1]
const yzm = SA[1 0 0 0; 0 0 1 0; 0 -1 0 0; 0 0 0 1]

const ywp = SA[1 0 0 0; 0 0 0 -1; 0 0 1 0; 0 1 0 0]
const ywm = SA[1 0 0 0; 0 0 0 1; 0 0 1 0; 0 -1 0 0]

const zwp = SA[1 0 0 0; 0 1 0 0; 0 0 0 -1; 0 0 1 0]
const zwm = SA[1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 -1 0]

function create_elements()
    elements = [SA[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]]
    counter = 0
    while length(elements) < 192 || counter == 16
        if true
            new_elements = [xyp*element for element in elements]
            for (i,e) in enumerate(new_elements)
                if e ∉ elements
                    push!(elements,e)
                end
            end
        end

        if true
            new_elements = [xzp*element for element in elements]
            for (i,e) in enumerate(new_elements)
                if e ∉ elements
                    push!(elements,e)
                end
            end
        end

        if true
            new_elements = [xwp*element for element in elements]
            for (i,e) in enumerate(new_elements)
                if e ∉ elements
                    push!(elements,e)
                end
            end
        end

        if false
            new_elements = [yzp*element for element in elements]
            for (i,e) in enumerate(new_elements)
                if e ∉ elements
                    push!(elements,e)
                end
            end
        end

        if false
            new_elements = [ywp*element for element in elements]
            for (i,e) in enumerate(new_elements)
                if e ∉ elements
                    push!(elements,e)
                end
            end
        end

        if false
            new_elements = [zwp*element for element in elements]
            for (i,e) in enumerate(new_elements)
                if e ∉ elements
                    push!(elements,e)
                end
            end
        end
        counter += 1
    end
    return elements
end

function verifyelements(elements)
    if length(unique(elements)) != 192
        return false
    end
    ones_1x4 = SA[1 1 1 1]
    ones_4x1 = reshape(ones_1x4,(4,1))
    for el in elements
        if size(el) != (4,4)
            return false
        end
        if !issubset(unique(el),[0,1,-1])
            return false
        end
        if map(abs,sum(el,dims=1)) != ones_1x4
            return false
        end
        if map(abs,sum(el,dims=2)) != ones_4x1
            return false
        end
    end
    return true
end

let elements_file = joinpath(@__DIR__,"elements.jld2")
    timestarted = time()
    verified = false
    howstring = "loaded"
    if isfile(elements_file)
        println("Loading group elements.")
        elements = load(elements_file,"elements")
        verified = verifyelements(elements)
    end
    if !verified
        howstring = "created"
        println("Elements could not be loaded. Creating group elements.")
        global elements = create_elements()
        save(elements_file,"elements",elements)
    end
    timeelapsed = time()-timestarted
    println("All elements $howstring after $timeelapsed s.")
end

function write_elements()
    f = open(joinpath(@__DIR__,"elements.jl"),"w")
    write(f,"elements = [")
    for el in elements[1:end-1]
        write(f,"SA")
        print(f,el)
        write(f,",\n")
    end
    write(f,"SA")
    print(f,last(elements))
    write(f,"]\n")
    close(f)
end

write_elements()

function order(e)
    r = e*e
    n = 1
    while r != e
        r = r*e
        n+=1
    end
    return n
end

function group_act_slow(a::Integer,b::Integer)
    if 192 >= a > 0 && 192 >= b > 0
        for i = 1:192
            if elements[i] == elements[a]*elements[b]
                return typeof(a)(i)
            end
        end

        println(elements[a])
        println(elements[b])
        println(elements[a]*elements[b])
        error("Not found? How is that even possible?")
    else
        error("Not defined.")
    end
    return zero(typeof(a))
end

const g = [group_act_slow(a,b) for a in 1:192, b in 1:192]

function ⊗(a::Integer,b::Integer)
    if 192 >= a > 0 && 192 >= b > 0
        return g[a,b]
    else
        error("$a o $b Not defined.")
    end
end

function insertsorted!(set,element;startindex = 1)
    if !issorted(set)
        error("Array is not sorted.")
    end
    i = startindex
    while  i <= length(set) && element > set[i]
        i += 1
    end
    i > length(set) ? push!(set, element) : insert!(set, i, element)
    return i
end

function complete!(set)
    added = -1
    queue = vec([(a,b) for a in set, b in set])
    while added != 0
        added = 0
        newElements = empty(set)
        for i in length(queue):-1:1
            a,b = queue[i]
            if ⊗(a,b) ∉ set
                push!(newElements,⊗(a,b))
                added = added+1
            end
            pop!(queue)
        end
        newElements = sort(unique(newElements))
        for a in newElements
            for b in set
                push!(queue,(a,b))
                push!(queue,(b,a))
            end
            for b in newElements
                push!(queue,(a,b))
            end
        end

        i = 1
        for a in newElements
            i = insertsorted!(set,a, startindex = i)
        end
    end
    return set
end

function completewith!(set,toadd)
    toadd = setdiff(toadd,set)
    isempty(toadd) && return set
    added = -1
    queue = vcat(vec([(a,b) for a in set, b in toadd]),vec([(b,a) for a in set, b in toadd]),vec([(a,b) for a in toadd, b in toadd]))
    for a in toadd
        insertsorted!(set,a)
    end
    while added != 0
        added = 0
        newElements = empty(set)
        for i in length(queue):-1:1
            a,b = queue[i]
            if ⊗(a,b) ∉ set
                push!(newElements,⊗(a,b))
                added = added+1
            end
            pop!(queue)
        end
        newElements = sort(unique(newElements))
        for a in newElements
            for b in set
                push!(queue,(a,b))
                push!(queue,(b,a))
            end
            for b in newElements
                push!(queue,(a,b))
            end
        end

        i = 1
        for a in newElements
            i = insertsorted!(set,a, startindex = i)
        end
    end
    return set
end

function getsubgroups(;maxgenerators = 4)
    subgroupdict = Dict()
    subgroupdict[(generators = 0,)] = [[1]]
    for n = 1:maxgenerators
        timestarted = time()
        println("Generating with $n generators.")
        subgroupbase = subgroupdict[(generators = n-1,)]
        subgroupstoadd_duplicates = vcat([[completewith!(copy(s),[a]) for s in subgroupbase] for a in 1:192]...)
        subgroupstoadd = setdiff(unique(subgroupstoadd_duplicates),vcat([subgroupdict[(generators = i,)] for i = 0:n-1]...))
        subgroupdict[(generators = n,)] = subgroupstoadd
        timeelapsed = time()-timestarted
        println("  Generated after $timeelapsed s.")
    end
    return subgroupdict
end

function addgenerator(subgroups,generator)
    newsubgroups = empty(subgroups)
    for subgroup in subgroups
        newsubgroup = copy(subgroup)
        generator ∈ subgroup && continue
        completewith!(newsubgroup,[generator])
        push!(newsubgroups,newsubgroup)
    end
    return newsubgroups
end

function isgroup(set)
    for a in set
        for b in set
            if ⊗(a,b) ∉ set
                return false
            end
        end
    end
    return true
end

function verifysubgroupdict(subgroupdict)
    #TODO: implement subgroupdict verification
    return true
end

let subgroupdict_file = joinpath(@__DIR__,"subgroupdict.jld2")
    timestarted = time()
    verified = false
    howstring = "loaded"
    if isfile(subgroupdict_file)
        println("Loading all subgroups.")
        subgroupdict = load(subgroupdict_file,"subgroupdict")
        verified = verifysubgroupdict(subgroupdict)
    end
    if !verified
        howstring = "generated"
        println("Generating all subgroups.")
        global subgroupdict = getsubgroups()
        for igens = 1:4
            permute!(subgroupdict[(generators = igens,)],sortperm(length.(subgroupdict[(generators = igens,)])))
        end
        subgroupdict[:allsubgroups] = vcat([subgroupdict[(generators = i,)] for i = 0:4]...)
        save(subgroupdict_file,"subgroupdict",subgroupdict)
    end
    timeelapsed = time()-timestarted
    println("All subgroups $howstring after $timeelapsed s.")
end

function write_subgroups()
    f = open(joinpath(@__DIR__,"subgroups.jl"),"w")
    write(f,"subgroupdict = Dict(")
    for key in keys(subgroupdict)
        print(f,key)
        write(f," => ")
        print(f,subgroupdict[key])
        write(f,",\n")
    end
    skip(f,-2)
    write(f,")\n")
    close(f)
end

write_subgroups()

const vertices = vec([SA[a,b,c,d] for a in [-1,1], b in [-1,1], c in [-1,1], d in [-1,1]])

function element_as_vertex_perm(element)
    permuted_vertices = [element*vertex for vertex in vertices]
    vertex_perm = [findfirst([v_ori==v_new for v_ori in vertices]) for v_new in permuted_vertices]
    return Permutation(vertex_perm)
end

elements_as_vertex_perms = [element_as_vertex_perm(e) for e in elements]

function allequal(array)
    previous = first(array)
    result = true
    for val in array
        result = result && isequal(previous,val)
        previous = val
        !result && continue
    end
    return result
end

function isequalcoloring(coloring1,coloring2)
    if length(unique(coloring1)) != length(unique(coloring2))
        return false
    end
    for c in unique(coloring1)
        v1 = coloring1 .== c
        if !allequal(coloring2[v1])
            return false
        end
    end
    return true
end

function get_subgroup_properties()
    subgroup_properties = [Dict() for s in subgroupdict[:allsubgroups]]
    for (i,sdict) in enumerate(subgroup_properties)
        sdict[:subgroup] = subgroupdict[:allsubgroups][i]
        sdict[:order] = length(sdict[:subgroup])
        j = i
        ngenerators = -1
        while j > 0
            ngenerators += 1
            with_ngenerators = length(subgroupdict[(generators = ngenerators,)])
            with_ngenerators == 0 && break
            j -= with_ngenerators
        end
        sdict[:ngenerators] = ngenerators
    end

    for (is,sdict) in enumerate(subgroup_properties)
        sdict[:coloring] = collect(1:16)
        for (iv,v) in enumerate(vertices)
            group_action_on_v = [elements[e]*v for e in sdict[:subgroup]]
            for connected_vertex in group_action_on_v
                i_con_vert = 1
                while vertices[i_con_vert] != connected_vertex && i_con_vert < 17
                    i_con_vert = i_con_vert+1
                end
                #if connected_vertex has lower color than v, leave it be (should not happen)
                sdict[:coloring][i_con_vert] = min(sdict[:coloring][i_con_vert],sdict[:coloring][iv])
            end
        end
    end

    for (is,sdict) in enumerate(subgroup_properties)
        subgroupcoloring = sdict[:coloring]
        newcolorvalue = 1
        newsubgroupcoloring = copy(subgroupcoloring)
        for colorvalue in sort(unique(subgroupcoloring))
            newsubgroupcoloring[subgroupcoloring.==colorvalue] .= newcolorvalue
            newcolorvalue = newcolorvalue+1
        end
        sdict[:coloring_sorted] = newsubgroupcoloring
    end

    for (is,sdict) in enumerate(subgroup_properties)
        sdict[:coloring_equiv] = zeros(Int,16)
    end

    for (is,sdict) in enumerate(subgroup_properties)
        sgcolor = sdict[:coloring_sorted]
        found = false
        for p in elements_as_vertex_perms
            sgcolorperm = sgcolor[[getindex(p,k) for k in axes(sgcolor,1)]]
            for other_sdict in subgroup_properties
                if isequalcoloring(sgcolorperm,other_sdict[:coloring_equiv])
                    found = true
                    sdict[:coloring_equiv] = other_sdict[:coloring_equiv]
                    continue
                end
            end
            found && continue
        end
        found && continue
        sdict[:coloring_equiv] = sgcolor
    end

    for (is,sdict) in enumerate(subgroup_properties)
        sdict[:coloring_equiv_index] = 0
    end

    coloring_equiv_index = 1
    for (is,sdict) in enumerate(subgroup_properties)
        found = false
        for other_is in 1:is-1
            if sdict[:coloring_equiv] == subgroup_properties[other_is][:coloring_equiv]
                found = true
                sdict[:coloring_equiv_index] = subgroup_properties[other_is][:coloring_equiv_index]
                continue
            end
        end
        found && continue
        sdict[:coloring_equiv_index] = coloring_equiv_index
        coloring_equiv_index += 1
    end

    for col_eq_i = 1:coloring_equiv_index-1
        subgs = [sdict[:coloring_equiv_index] == col_eq_i for sdict in subgroup_properties]
        if length(unique([sdict[:ngenerators] for sdict in subgroup_properties[subgs]])) > 1
            n = length(unique(first(subgroup_properties[subgs])[:coloring_equiv]))
            println("Varying number of generators for coloring_equiv_index = $col_eq_i with $n colors")
        end
    end

    mirrorx = [-1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
    elements_as_vertex_perms_mirror = vcat([element_as_vertex_perm(e) for e in elements],[element_as_vertex_perm(e*mirrorx) for e in elements])

    for (is,sdict) in enumerate(subgroup_properties)
        sdict[:coloring_equiv_mirror] = zeros(Int,16)
    end

    for (is,sdict) in enumerate(subgroup_properties)
        println("Processing $is-th subgroup.")
        sgcolor = sdict[:coloring_sorted]
        found = false
        for p in elements_as_vertex_perms_mirror
            sgcolorperm = sgcolor[[getindex(p,k) for k in axes(sgcolor,1)]]
            for other_sdict in subgroup_properties
                if isequalcoloring(sgcolorperm,other_sdict[:coloring_equiv_mirror])
                    found = true
                    sdict[:coloring_equiv_mirror] = other_sdict[:coloring_equiv_mirror]
                    continue
                end
            end
            found && continue
        end
        found && continue
        sdict[:coloring_equiv_mirror] = sgcolor
    end

    for (is,sdict) in enumerate(subgroup_properties)
        sdict[:coloring_equiv_mirror_index] = 0
    end

    coloring_equiv_mirror_index = 1
    for (is,sdict) in enumerate(subgroup_properties)
        found = false
        for other_is in 1:is-1
            if sdict[:coloring_equiv_mirror] == subgroup_properties[other_is][:coloring_equiv_mirror]
                found = true
                sdict[:coloring_equiv_mirror_index] = subgroup_properties[other_is][:coloring_equiv_mirror_index]
                continue
            end
        end
        found && continue
        sdict[:coloring_equiv_mirror_index] = coloring_equiv_mirror_index
        coloring_equiv_mirror_index += 1
    end

    for col_eq_i = 1:coloring_equiv_mirror_index-1
        subgs = [sdict[:coloring_equiv_mirror_index] == col_eq_i for sdict in subgroup_properties]
        if length(unique([sdict[:ngenerators] for sdict in subgroup_properties[subgs]])) > 1
            n = length(unique(first(subgroup_properties[subgs])[:coloring_equiv_mirror]))
            println("Varying number of generators for coloring_equiv_index = $col_eq_i with $n colors")
        end
    end
    return subgroup_properties
end

let subgroup_properties_file = joinpath(@__DIR__,"subgroup_properties.jld2")
    println("Calculating subgroup properties.")
    global subgroup_properties = get_subgroup_properties()
    save(subgroup_properties_file,"subgroup_properties", subgroup_properties)
end

function write_groupmatrix()
    f = open(joinpath(@__DIR__,"groupmatrix.jl"),"w")
    groupmatrixname = :g
    write(f,"const $(string(groupmatrixname)) = [")
    for row in axes(eval(groupmatrixname),1)
        for col in axes(eval(groupmatrixname),2)
            write(f,"$(g[row,col]) ")
        end
        write(f,";\n")
    end
    skip(f,-2)
    write(f,"]\n")
    close(f)
end

write_groupmatrix()