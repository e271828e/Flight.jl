module System

using Flight.LBV

export AbstractContinuousSystem, get_x_type, get_u_type

abstract type AbstractSystem end

abstract type AbstractContinuousSystem{X, U <: Union{Nothing, AbstractLBV}} <: AbstractSystem end

# get_x_type(::Type{<:AbstractContinuousSystem{X, U}}) where {X, U} = X
# get_u_type(::Type{<:AbstractContinuousSystem{X, U}}) where {X, U} = U
# get_x_type(::AbstractContinuousSystem{X, U}) where {X, U} = get_x_type(AbstractContinuousSystem{X,U})
# get_u_type(::AbstractContinuousSystem{X, U}) where {X, U} = get_u_type(AbstractContinuousSystem{X,U})

end









# ejemplo: un Subsystem muReg contiene un x::Leaf{2}. cuando vaya a ensamblar el
# descriptor XLdg del parent system Ldg, defino n campo reg=typeof(sys.reg.x).

# Seria interesante que un System pudiera tener: x_disc x_cont

# x_disc puede representar estados discretos del sistema (por ejemplo, una
# maquina de estados) y tipicamente respondera a una ecuacion en diferencias.
# Podria definir por ejemplo un LBV AircraftXStateMachine con fields act, pwp,
# etc.

# Ahora, en step hace dos cosas: x_disc = f_disc(x_disc, x_cont, u) x_cont =
# integrate(f_cont(x_cont, x_disc, u))

# Todo sistema manejado por el scheduler deberá ser discreto. Un sistema
# discreto tendrá al menos un method init, otro step, otro shutdown. Si dentro
# contiene un Continuous System, lógicamente estará discretizado en un numerical
# method para exponer esos methods discretos.

# Y si tengo un hybrid system? Pues entonces necesitaré un subtype de Discrete
# Sys que además de integrar numéricamente el Continuous que envuelve, le puede
# hacer cosas en cada llamada a step

#para un Aircraft <: System, que deberia ser Aircraft? abstract type no,
#claramente. entre diferentes tipos de Aircraft no hay solo afinidad en cuanto a
#methods, sino tambien en cuanto a datos e implementacion. la duda es si definir
#solo un tipo Aircraft, y que la unica particularidad de cada subtipo este en
#los fields, o (y esta parece la opcion mas recomendable) un parametric type,
#que permite hacer dispatch para definir pwp_group, ldg_group, etc, con los
#dispatches correspondientes. p. ej, Aircraft{:Cessna172}