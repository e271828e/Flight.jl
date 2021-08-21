1. update Model for the new ContinuousSystem. test it both with EThruster and a group of EThrusters
1. finish TestAircraft and test it in Simulation
2. create a new branch new_arch_with_struct_y to change back to nested struct arrays for Y. Update
   Model and remove all Y methods. Remove field y from ContinuousSystem, and modify f_cont!() to
   return y. Modify NewContinuousModel to account for all this. benchmark and compare with the
   committed implementation.
3. if successful, return Wrench() to a immutable struct of SVector{3}, see the time improvements

if instead of having a ComponentVector Y as the only output, i define a discrete Yd and a continuous
Y and pack them in a struct to be saved by DiffEq on each time step, the resulting log would no
longer be a Vector of ComponentVectors, but a Vector of structs. this means i would have to extract
both y and yd for each t and assemble them into two separate VectorOfArrays. that is ugly to do by
hand. StructArrays to the rescue!!!

```
Base.@kwdef struct Output{Y, YD} #this is the type we pass to SavedValues
y::Y = ComponentVector(a = fill(1.0, 2), b = fill(2.0, 3))
yd::YD = ComponentVector(m = fill(4,3), n = fill(-2, 2))
end
log = collect(Output() for i in 1:5) #create some copies of it
sa = StructArray(log) #now all the y fields of log lie in a contiguous array
y = sa.y
y_voa = VectorOfArray(y) #now we have a vector of Y's that indexes like a matrix
y_mat = convert(Array, y_voa) #and a matrix of y's whose rows still preserve axis metadata
```

 with this, Y no longer needs to be an array! we could use nested immutable structs or NamedTuples
thereof. because StructArrays not only accepts Vectors of structs as inputs but also Vectors of
NamedTuples. and, if we keep everything immutable, it no longer has to be inefficient to allocate
it. the only problem is doing a deepcopy on each step.

but NO! this does not have to be. we can simply define y to be the output of f_cont!(sys). being an
immutable struct freshly generated from scratch on every iteration, it can go directly to saveval
without deepcopying, because it is NOT the same object as the one in the next step. #in fact, with
this approach, we no longer need to define y as a System field!!

ok, but what happens if we store things like RQuat or Wrench objects? Vectors, and therefore cannot
be handled by the VectorOfArray workflow. on the other hand, it doesn't have to. because we can
apply StructArray recursively. for example, if we get an array y_kin of YKin structs, each of them
containing a RQuat q_lb and another RQuat q_el, first we do sa = StructArray(y_kin), then v_q_lb =
sa.q_lb. now we have a Vector of RQuats. yep, but the pain does not stop here, because then RQuat
holds a UnitQuat, which in turn holds a Quat. this is not practical. so here is a self-imposed
limitation to facilitate simulation output post-processing: only use fields of SVectors of Floats,
Ints, Bools or non-deeply nested structs. A Wrench is acceptable because it is handled directly by
StructArrays to extract F and M. a WGS84 pos maybe too.

advantages of the nested struct approach:
1) we can hold different data types in a single struct, so we no longer need y and yd
2) for postprocessing we can use StructArrays to extract a Vector of outputs corresponding to a
   subsystem for plot dispatch, so we're still good there
3) we don't need to make Wrench a ComponentVector, which recovers some performance
4) because the overall y no longer needs to be allocated we don't need Y methods, the strct is
   directly instantiated in the fcont call we only need to define Y structs for Leaf components. in
   most other cases they will be compose using NTs
4) we can use ComponentGroups and other forms of composition without worrying about having to deal a
   particular component returning nothing for y. a NT of nothing is a NT!

disadvantages:
1) although they are stack allocated, these nested structs have still some allocation time. but this
   is probably comparable or even less than the cost of in-place assignment to the heap-preallocated
   yblocks
2) we're still restricted to storing SVectors, Floats and such things to make post-processing
   easier. however, we can still store them in y and pass them around for other uses


<!-- nah, it IS better to have y as a preallocated ComponentVector with views stored
throughout the hierarchy. especially if we can define another y_d for discrete
outputs if required. -->