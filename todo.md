Aircraft type design concept:

We start from a parametric AbstractAircraft{Mass, Aerodynamics, Surfaces, Powerplant, LandingGear, Fuel, Electrics} type.

As soon as we have some aerodynamics, surfaces will have to be dealt with.
Their actuation dynamics and current position will be handled by the Surfaces type.
Note that as long as we keep the Aircraft struct immutable, we don't really need to store EVERYTHING in the type
parameters. The type parameters may contain only types, which are then referenced by the Aircraft struct fields.
For each of these parameters we can provide NoAerodynamics, NoSurfaces, NoPowerplant, NoLandingGear, etc. that provide
trivial implementations of.
we can exploit multiple dispatch for this. for example, we can define:
get_mass_data(mass::ConstantMassModel, args...) = mass.data
but:
get_mass_data(mass::C172MassModel, y::NoFuel) = #returns OEW mass data#
get_mass_data(mass::C172MassModel, y::C172FuelY) = #some more complex computation#

or:
get_wr_Ob_b(aero_model::C172Aero, y_srf::C172SrfY, y_ldg, trn::NoTerrain) = #consider only surfaces#
get_wr_Ob_b(aero_model::C172Aero, y_srf::C172SrfY, y_ldg, trn::Abstract) = #consider only surfaces#
#both of these work with the same call: get_mass_data(aircraft.mass, y_aircraft.fuel)


### Controls

For each Aircraft concrete type we must define:
1. A set of continuous and discrete inputs
2. A mapping of these inputs to the different Aircraft components' inputs

AbstractControlSet{Aircraft}
For 1. we need to define a mutable struct C172Controls <: AbstractControls
For 2. we define a method with interface@

map_controls(u::C172Controls, u_srf::AbstractSurfaceControls, pwp::AbstractPowerplant, ldg::AbstractLandingGear,
            fuel::AbstractFuelSystem, etc)

            #this is called before updating the aircraft components


U{ControlMap, Surfaces, Powerplant, LandingGear, Fuel, Electrics}
Each Controls subtype should provide a mutable struct (that may be defined by composing the U structs on each of the

aircraft's systems. then method that dispatches on the Controls type
itself, but also on the rest of the Aircraft systems: Srf, Ldg

a thought: the control inputs for each aircraft subsystem should reside statically in memory! they should not be
allocated from scratch on each call to f_cont!, no matter how cheap they'd be if they were immutable. it is extremely
inconvenient to do this all the time. even if it means a tiny performance cost when mutating heap-allocated data, they
should be mutable and preallocated somehwere. they can be ComponentVectors or mutable structs (more on this later), but
they should be mutable.

but where they should reside? we have established that aircraft inputs and aircraft subsystem inputs are different
things, and the mapping between them is NOT one to one. the aircraft inputs can be exposed to the outside world so they
can be modified (this is what we were doing with U). but how about system inputs in general? well, a possible design is
to force every system to include a :u field that

how do we define f_cont! in this case>


### GenericAircraft


There are some parts of an aircraft f_cont! that should be common:
- Updating x_kin! and getting y_kin at the beginning,
- Mapping its control inputs to the
updating x_dyn! and getting y_acc at the end, performing an update of its subsystems (with a call to f_cont_cmp!), and obtaining wr_Ob_b and h_Gc_b for all components and summing them up.

Now, the specifics are in f_cont_cmp! and f_disc_cmp! no generic implementation should be provided for these, because
they are totally dependent on


there is no performance harm in going to mutable structs for Y, as long as we assemble and allocate the whole Y
hierarchy once at the aircraft level and then mutate its elements as required. the problem with the nested struct
approach is that we no longer have the convenience of VectorOfArrays to generate a labelled matrix from which to pass
blocks to each system.

the only real reason to have nested structs for Y rather than plan ComponentVectors is the need for storing integers,
that is, discrete variables, rather than having everything flattened to Floats. to address this, we could probably
define a discrete Y and a discrete U for each system and System of systems. however, how can we save them separately
from the actual continuous y?




renormalization goes to Kinematics with tolerance as an optional argument and is called by the generic Aircraft f_disc!





Cambiar U de ComponentVectors a mutable structs. ComponentGroup las empaqueta en un NamedTuple. No importa que sea
immutable, porque para un ComponentGroup lo que queremos es poder fijar la U de cada componente por separado, y esas es
de esperar que, o bien sean ComponentVectors, o bien sean mutable structs cuyos campos podamos modificar. Cada subtype

Si un cierto Aircraft define su U a priori, en principio no deberia necesitar



Desde el punto de vista de pintar, sigue siendo interesante que Y sea un ComponentVector

Alternativa: es posible que cada subsystem defina localmente como struct field una x, xdot u e y? y que en vez de
llamar a sus f_cont! con argumentos lo hacemos sin ellos, asignando previamente sus blocks? cual es la diferencia? al
final el x grande del que vengan los blocks va a tener que ser allocatado en algun sitio fuera.

la ventaja de que cada sistema tenga internamente campos x, u, dx e y que apunten como views a blocks de vectores
externos es que puedo pasarle cada sistema, con su estado encapsulado, a diferentes methods, como un todo en uno. Ya no
necesito pasar al mismo tiempo los parametros del sistema, su u, y su u

que hago? dos alternativas:
1. Immutable structs con campo params separado, si queremos que sus dx, x, etc sean views a blocks de vectores externos,
   hay que instanciarlas desde el principio pasandoles esos blocks

A VER: Si quiero que cada subsistema contenga sus propios x, y, u, ya no hay separacion posible entre params y esos
vectores. Yo tendre un Aircraft con subsystems.
```
struct System
subsystems::S
dx::X
x::X
y::Y
u::U
end
```

This means that the types are not determined until the System is instantiated. this may lead to instability. but worst
case scenario, we can always define it as const once it is created. immutability should be enough but...

```
struct MyAircraft <: AbstractComponent
pwp::ConcretePwp
ldg::ConcreteLdg
srf::ConcreteSrf
end

function System(ac_params::MyAircraft, dx_ac = X(MyAircraft), x_ac = X(MyAicraft), y_ac =...)
pwp = System(ac_params.pwp, dx_ac.pwp, x_ac.pwp, y_ac.pwp, u_ac.pwp)
ldg = System(ac_params.ldg, dx_ac.ldg, x_ac.ldg, y_ac.ldg, u_ac.ldg)
srf = System(ac_params.srf, dx_ac.srf, x_ac.srf, y_ac.srf, u_ac.srf)
#although of course some of those blocks may not exist, for example, ldg might not have u
subsystems = (pwp = pwp, ldg = ldg, srf = srf)
System(subsystems, dx, x, y, u) #this dispatches to the generic inner constructor
end
```

But now, how does each of these inner constructor calls look?
```
struct ConcretePwp <: AbstractComponent
left::ConcreteThruster
right::ConcreteThruster
end

function System(pwp_params::ConcretePwp, dx_pwp = X(ConcretePwp), y_pwp = ...)
left = System(pwp_params.left, dx_pwp.left...)
right = System(pwp_params.right, dx_pwp.right...)
subsystems = (left = left, right = right)
System(subsystems, dx_pwp, x_)
end
```

note that we could define AbstractComponentGroup <: AbstractComponent and then define a @generated System constructor
with the interface:
@generated function System(cmp::C, dx = X(C), x = X(C)) where {C<:AbstractComponentGroup} that automates all the process
above if it is guaranteed that all members have states and outputs (which we have agreed that they must have)

if a specific Component for example has no inputs, it should extend the System constructor omitting those arguments.
For example, if we have an inputless LandingGear configuration, we could do
function System(ldg::InputlessLdg, dx::InputlessLdgX = X(InputlessLdg), x::Inputl..., y::Inputles...) #OMITTING U
end

the only constraint of this approach is that every parent system must be aware of the interface of each of the
subsystems it is using: if it has states or not, inputs or not, etc. THIS IS NO BIG ASK. This is absolutely sensible.


2. Mutable structs que permitan reasignar X Y U una vez creadas a las views que necesite

Esto me va a perjudicar en algun momento con DifferentialEquations? En principio no deberia. Porque aunque a priori las
structs que forman parte de la jerarquia de Aircraft sean mutables, y sus campos x, y, u tengan una cierta ambiguedad
porque no se sabe a priori si van a ser un ComponentArray como view o como vector, la realidad es que a
DifferentialEquations esto no deberia afectarle. El integrator trabaja con un vector X que es el root del sistema. En
ningun momento el tipo de X puede cambiar. siempre va a ser el que le especifique al principio en la condicion inicial.

Ahora, la historia es distinta dentro de las f_cont! de los diferentes subsistemas. Porque internamente esa f_cont! si
que va a echar mano de los x, y, u internos. Y si los tipos de esos no estan definidos desde el principio, entonces si
que voy a tener un problema. como podria solucionarlo? Puedo hacer:
data = SVector{length(XTemplate), Float64}(x); x = ComponentVector(data, XAxes). Y ya trabajar con esto como immutable.
es cierto que el compilador no tiene ninguna garantia
despues construir un ComponentVector

Probar este concepto con Powerplant a ver que pasa


---

Ahora, hay alguna forma de definir Systems concretos en los que X, U e Y esten restringidos por ejemplo a cierto tipo de
ComponentVectors?
Si. Por ejemplo, yo tengo:
System(params::P, dx = X(P), x = X(P), u = U(P), y = ...)
Pero me puedo definir:
System(params::EThrusterP, dx::EThrusterX, x::EThrusterX, u::EThrusterU, etc)

De manera que si alguien me intenta pasar un x que se ha definido en vez de apoyarse en el default, pero ese x no encaja
con mi EThrusterX, no se lo permito.

Yo tengo

struct ParametricAircraft{Pwp, Ldg, Srf}
pwp::Pwp
ldg::Ldg
srf::Srf
end


so there are two paths:
1. Give up on parametric aircraft and assume that any aircraft will have its own ways of defining its X and
   instantiating its corresponding System

struct MyAircraft
pwp::Powerplant #these do not need to be type parameters now
ldg::Ldg
srf::Srf
end

const MyAircraftXTemplate = assemble manually from individual subsystems by calling their
MyAircraftX = Component...
X(::Type{MyAircraftParams}) = copy(MyAircraftXTemplate)

f_kin!(aircraft.y.kin, aircraft.dx.pos, aircraft.x.kin)# this mutates aircraft

2. Use a parametric aircraft but give up on checking the passed state and input vectors for the right type

X(ac::Aircraft{Pwp,Ldg,Srf}) = X(Aircraft{Pwp,Ldg,Srf}) #fallback to type
X(::Type{Aircraft{Pwp,Ldg,Srf}}) where {...} = ComponentVector(pwp = X(Pwp), ldg = X(Ldg), srf = X(Srf))
#this is the default implementation, but if any of those returns nothing, we have to do it manually in an ad-hoc method
to which we dispatch based on the specific type parameters. alternatively, we could simply check with if statements. or
use traits HasStates() or HasNoStates(). or more simply, build a dictionary, push the return values of the X calls, drop
those that are nothing, and call ComponentVector. in any case, this method must return an appropriate ComponentVector



end










aaaaaaaaaaaaa