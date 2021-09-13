module Utils

export pwf, no_extend_error, no_extend_warning

function pwf(s)#print with fieldnames
    for f in fieldnames(typeof(s))
        println("$f: $(getfield(s,f))")
    end
end

no_extend_error(f::Function, ::Type{S}) where {S} = error(
    "Function $f not implemented for type $S or incorrect call signature")

no_extend_warning(f::Function, ::Type{S}) where {S} = println(
    "Warning: Function $f not implemented for type $S or incorrect call signature")


end