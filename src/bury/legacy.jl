# Original algorithm bodies retained as an independent compatibility backend.
# Public wrappers validate inputs before entering these unchanged bodies.

function _bury_legacy(g::Graphs.Graph, r0; verbose=false)
    r = r0

    # Node weight initialization
    weights = [0.0 for _ in 1:nv(g)]
    colored = [false for _ in 1:nv(g)]
    for (index,v) in enumerate(vertices(g))
        weights[index] = Graphs.degree(g,v) + 1
    end

    if verbose
        println("Initial weights and colors:")
        println("\nWeights table")
        println(weights)
        println("Colored table")
        println(colored)
    end

    # Global greedy selection
    while sum(colored)<r0
        cost, choice = findmin(weights)

        if verbose
            println("\n\n###########################\nStarting round. r=", r)
            println("Chose vertex ", choice, " for the price of ", cost)
        end

        if r>= cost
            nbd = vcat(choice,neighbors(g,choice))
            r -= cost
        else
            nbd = choice
            r -= 1
        end

        for v in nbd
            if !colored[v]
                weights[v] -= 1
                weights[neighbors(g,v)] .-= 1
                colored[v] = true
            end
        end
        weights[choice] = Inf
        if verbose
            println("\nUpdated weights table")
            println(weights)
            println("Updated Colored table")
            println(colored)
        end
    end

    new_reg = Dict()
    for i in 1:nv(g)
        new_reg[i] = colored[i]+1
    end
    return new_reg, colored
end

function _bury_seeded_legacy(g::Graphs.Graph, r0, seed; verbose=false)
    r = r0

    # Node weight initialization
    weights = [0.0 for _ in 1:nv(g)]
    colored = [false for _ in 1:nv(g)]

    region = [false for _ in 1:nv(g)]
    region[neighbors(g,seed)] .= true
    region[seed] = true

    for (index,v) in enumerate(vertices(g))
        weights[index] = Graphs.degree(g,v) + 1
    end

    if verbose
        println("Initial weights and colors:")
        println("\nWeights table")
        println(weights)
        println("Colored table")
        println(colored)
    end

    # Global greedy selections
    while sum(colored)<r0
        cost, choice = findmin(.!region*Inf + weights)
        if cost == Inf
            cost, choice = findmin(weights)
        end

        if verbose
            println("\n\n###########################\nStarting round. r=", r)
            println("Chose vertex ", choice, " for the price of ", cost)
        end

        if r>= cost
            nbd = vcat(choice,neighbors(g,choice))
            r -= cost
        else
            nbd = choice
            r -= 1
        end

        for v in nbd
            if !colored[v]
                weights[v] -= 1
                weights[neighbors(g,v)] .-= 1
                colored[v] = true
            end
            region[v] = true
        end
        weights[choice] = Inf

        if verbose
            println("\nUpdated weights table")
            println(weights)
            println("Updated Colored table")
            println(colored)
        end
    end

    new_reg = Dict()
    for i in 1:nv(g)
        new_reg[i] = colored[i]+1
    end
    return new_reg, colored
end
