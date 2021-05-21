module Attitude

using StaticArrays: MVector #this form of using does not allow method extension

# https://discourse.julialang.org/t/writing-functions-for-types-defined-in-another-module/31895
# https://discourse.julialang.org/t/what-is-the-preferred-way-to-use-multiple-files/30687/6

#do not do this:
# include("quaternions.jl")
# import .Quaternions
#This pastes the complete quaternions.jl, including the Quaternions module, in
#this file. But Flight.jl includes both quaternions.jl and attitude.jl.
#Therefore, in Flight.jl there will be two conflicting definitions of module
#Quaternions, one from the direct inclusion of quaternions.jl, and another
#through the inclusion of quaternions.jl in attitude.jl. this forces Julia to
#qualify all the types exported by the Quaternions module with their complete
#paths to keep them distinct and avoid ambiguity

#instead...
using ..Quaternions: UnitQuat, UnitQuat64

export Rotation

#implements a passive (alias) rotation, that is, a rotation describing the
#relative attitude between two reference frames.

#outer constructors must take their input arguments, generate the field values
#required by the type, and then pass them to the inner constructor
#rmat, rvec, axisangle


#Attitude(UnitQuat64()) #calls the inner constructor directly

# Attitude(UnitQuat32()) #also. here, the inner constructor automatically tries
#to convert() the input into a UnitQuat64. since this conversion is provided by
#the Quaternion module, it works silently. obviously, something like
#Rotation("Hello") fails, because there is no convert() defined from string to
#UnitQuat64.

#enables keyword argument "quat" and doubles as a zero argument constructor
# Rotation(; quat::UnitQuat = UnitQuat64()) = Rotation(quat)

#Julia DOES NOT dispatch on keyword argument names. that is, creating outer
#constructor with rmat, rvec, euler, etc as keyword arguments and expecting
#Julia to choose the right one will fail. it follows that for this purpose, a
#single keyword constructor is required which handles all possible keyword
#combinations. this might not be efficient, because the keywords must first be
#captured in a Dict and then checked to decide which actual constructor to call

#https://discourse.julialang.org/t/keyword-argument-constructor-breaks-incomplete-constructor/34198/3

#########################
abstract type Rotation end

# declare RMat, RVec, RQuat, AxAng, Euler types and work with them, providing
# convert() methods

#Rotation is an abstract type, of which these are concrete subtypes. not all
#operations are defined for every one of them.

#how do i make RVec(input) work for input<:Rotation? Simply make the constructor
#convert input to RQuat, then from RQuat to RVec
#convert(RQuat, input)
#convert(RVec, input)
#in fact, since we know that we will support conversions to Quat for all
#Rotation subtypes, we can simply do:
#RVec(input::Rotation) = RVec(RQ)
#Euler(1, 2, 3)

#a RMat can be constructed from a 3x3 Matrix. in that case, orthonormality is
#always enforced by applying QR factorization. if it comes from a conversion, it
#is not necessary

#i can admit * of rotation matrices without converting to RQuat, and also
#transforming vectors! only use promote when absolutely required!

#Euler can be a struct with fields upsi, utheta, uphi

#replace argument types with StaticArrays to fix lengths. also, allow different
#Real subtypes


#DO NOT TRY TO DISPATCH ON KEYWORDS, DISPATCH ON TYPES. JESUS CHRIST!
# Rotation(rmat::SVector33)
# Rotation(axis::SVector3, angle::Real)
# Rotation(rvec::SVector3, angle::Real)

#NOW create a keyword constructor for clarity, which then dispatchs to the
#previous ones

#go to Flight folder
#enter package manager
#activate .
#do not do activate Flight. in that case, the test and using commands do not
#work. they expect to be run from Flight's parent folder. why??

#Attitude(quat = qtest)
#Attitude(axis = [0, 0, 0], angle = 0)
    #conver to quaternion and then call the quaternion constructor

#must fail when passed a 1D array of size other than 3
#how do we do this? by defining StaticArrays types


end