using Flight

#in this

#instantiate an aircraft System
aircraft = System(Cessna172R());

#create a default environment
env = System(SimpleEnvironment());


#trim and run a simulation. then program a simulation control callback to do
#something with the controls.

#then linearize
