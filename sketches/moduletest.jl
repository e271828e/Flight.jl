module Parent

#ChildA must be defined before the using .ChildA statement, that's why
#includes should come first
    module ChildA

    using ..Parent

    export child_function_a

    child_function_a(x) = 2*parent_function(x)

    end #ChildA

using Reexport
@reexport using .ChildA #exports child_function_a

export parent_function

const parent_constant = 1

parent_function(x) = parent_constant + x

end #Parent