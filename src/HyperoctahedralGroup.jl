module HyperoctahedralGroup

using FileIO, Permutations, StaticArrays

export get_elements, get_subgroupdict, get_subgroup_properties
export ⊗

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

function get_elements()
    let elements_file = joinpath(@__DIR__,"../deps/elements.jld2")
        if isfile(elements_file)
            elements = load(elements_file,"elements")
            return elements
        else
            error("Build the package first!")
        end
    end
end

function order(e)
    r = e*e
    n = 1
    while r != e
        r = r*e
        n+=1
    end
    return n
end

include(joinpath(@__DIR__, "../deps/groupmatrix.jl"))

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

function get_subgroupdict()
    let subgroupdict_file = joinpath(@__DIR__,"../deps/subgroupdict.jld2")
        if isfile(subgroupdict_file)
            subgroupdict = load(subgroupdict_file,"subgroupdict")
            return subgroupdict
        else
            error("Build the package first!")
        end
    end
end

function get_subgroup_properties()
    let subgroup_properties_file = joinpath(@__DIR__,"../deps/subgroup_properties.jld2")
        if isfile(subgroup_properties_file)
            subgroup_properties = load(subgroup_properties_file,"subgroup_properties")
            return subgroup_properties
        else
            error("Build the package first!")
        end
    end
end

end #of module
