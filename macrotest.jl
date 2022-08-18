macro x(expr); return trait_init(:SystemX, expr); end
macro u(expr); return trait_init(:SystemU, expr); end
macro y(expr); return trait_init(:SystemY, expr); end
macro s(expr); return trait_init(:SystemS, expr); end

function trait_init(trait::Symbol, expr::Expr)

    @assert expr isa Expr && expr.head == :(=) "Invalid syntax"

    target = expr.args[1]
    assignment = expr.args[2]

    rhs = :($(esc(assignment)))
    lhs = :($(esc(Systems.init))(::$(esc(trait))))

    if target isa Symbol
        push!(lhs.args, :(::$(esc(target))))
    elseif target isa Expr && target.head === :call
        args = target.args
        if length(args) == 1
            push!(lhs.args, :(::$(esc(args[1]))))
        elseif length(args) == 2
            push!(lhs.args, :($(esc(args[2]))::$(esc(args[1]))))
        else
            extra_args = args[3:end]
            for arg in extra_args
                push!(lhs.args, :($(esc(arg))))
            end
        end
    else
        error("Invalid syntax")
    end

    return :($lhs = $rhs)

end