using BlockArrays

a = [1:4;]
b = [5:7;]

block_a = BlockArray(a)
block_b = BlockArray(b)

block_ab = mortar([block_a, block_b])

println(typeof(block_ab))
println(block_ab)

block_ab[Block(2)][2] = 4985

println(block_b)

block_a[:] = 37ones(4)

println(block_ab)

#just need to initialize each system by getting the x BlockVectors of each
#of its children and assembling them into a BlockVector. this is done
#recursively until a Leaf system is reached. do we need to create an
#AbstractLeaf < AbstractSystem? or is the implementation of this initializer
#part of the concrete subsystem itself, depending on whether it is or not a
#Leaf. could also be a Trait.

#with this, the only piece missing is to define a namedtuple mapping Symbols
#representing children labels to Block numbers for quick access, and override
#the System's get property / getfield

MyNT = @NamedTuple{ldg_group::Int, pwp_group::Int}
children = (ldg_group  = Block(1), pwp_group = Block(2))
println(isa(children, MyNT))
child = :ldg_group
block_ab[children[child]]
