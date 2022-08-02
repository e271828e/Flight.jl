using Flight


#This tutorial introduces the modeling approach used by the package.

#the core is a System:
# Represents a system with hybrid (continuous and discrete) dynamics. Continuous dynamics are defined by an ODE function f_cont!, and discrete dynamics by a difference equation function f_disc!

# System Descriptor has dual purpose. Dispatch for init methods that return the for System{desc} , hold any data that parameterized the behavior of the System

#we'll start by defining

import