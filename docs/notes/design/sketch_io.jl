#User-level listing for the §12 runtime-periphery machinery and the §13.3
#staging torture test. Same illustrative, non-committed syntax as sketch.jl /
#sketch_decoder.jl.
#
#Scenario: a first-order low-pass filter whose input is driven simultaneously
#by a fictitious single-axis joystick (100 Hz, complete writer) and a GUI
#panel (60 Hz, sparse writer), with 50 Hz boundaries at pace 1. The
#interference on u_cmd is deliberate — it is what the exercise probes.
#
#This listing is identical under all three candidate staging shapes (§12.3);
#that shape-invariance is itself a finding: the shapes differ only in
#framework internals and in behavior under concrete interleavings (§13.3).

######################## Model #################################################

struct LowPassFilter <: AbstractComponent end

init_x(::LowPassFilter) = (x = 0.0,)
#no z, no modes, no events — the minimal continuous primitive
#no g_s1: identity publication → the state x appears as a stage-1 port
#no g_s2: nothing beyond the law itself, and no published intermediates → the
#law may live directly in f (no duplicate computation site to drift, §9.4)

f(::LowPassFilter, y, u, t) = (x = (u.u_cmd - y.x) / u.τ,)

######################## Assembly ##############################################

asm = Assembly()
add!(asm, :filt, LowPassFilter())

#root input slots: published by no component, written only by the frame-top
#drain; sources to the scheduler; the periphery's entire write surface (§12.3).
#The model does NOT republish its input — u_cmd is a table slot, hence in
#every snapshot for free.
add_input!(asm, :u_cmd, 0.0)
add_input!(asm, :τ,     1.0)    #runtime-adjustable → a slot, not a parameter

connect!(asm, (:filt, :u, :u_cmd), (:u_cmd,))
connect!(asm, (:filt, :u, :τ),     (:τ,))

######################## Joystick device #######################################

struct ToyStick <: IODevice          #one device kind (§12.4); ~today's IODevices
    period::Float64                  #polling cadence, wall-clock seconds
end

init!(js::ToyStick)     = nothing    #open SDL, claim hardware, ...
shutdown!(js::ToyStick) = nothing

#runs on the device task each cycle; may block or pace itself — its business
poll(js::ToyStick) = (axis = read_hw_axis(js),)

struct StickToFilter <: IOMapping end

#pure, device data → slot writes (levels, §12.3); runs on the device task,
#never inside the loop's frame. A complete writer: whole write-set every poll.
#(Today's assign_input!(mdl, mapping, data), retargeted from mutation to value.)
map_input(data, ::StickToFilter) = (u_cmd = exp_axis_curve(data.axis; strength = 0.5),)

#framework-owned device task loop, for reference (the user never writes it):
#   init! → [ data = poll(dev); stage!(cell, map_input(data, mapping)) ]* → shutdown!

######################## GUI panel #############################################

#per-component extension (§12.5), reusable wherever a LowPassFilter appears.
#Widgets name the component's OWN PORTS; the build resolves each port to a
#root input slot (→ live widget) or a driving component (→ read-only
#rendering with provenance). Here both ports are root-driven, so both
#sliders are live; the SAME method renders u_cmd read-only in a
#configuration where another component drives it.
function GUI.draw!(ctx, ::LowPassFilter)
    s = view(ctx)                            #snapshot view scoped to this instance

    GUI.text("t = $(s.t)")
    GUI.display_bar("state x",       s.x)
    GUI.display_bar("applied u_cmd", s.u_cmd)   #what the sim actually integrated

    GUI.input_slider!(ctx, :u_cmd, -1.0, 1.0; label = "commanded input")
    GUI.input_slider!(ctx, :τ,      0.1, 10.0; label = "time constant")

    GUI.pause_button!(ctx)                   #control plane, not staging (§12.6)
end

#u_cmd appears on screen twice by design: the display bar shows the APPLIED
#value (from the snapshot); the slider knob shows own-pending-else-snapshot
#(§12.5 peek rule). Under joystick interference the two legitimately differ
#within a frame.

######################## Setup and run #########################################

sim = Simulation(asm; h = 0.02)             #50 Hz boundaries = drain/publish cadence

attach!(sim, ToyStick(0.01), StickToFilter())   #attachment order ⇒ conflict
attach!(sim, GUI.Renderer())                    #precedence: GUI attached later
                                                #⇒ wins same-frame conflicts (§12.3)
run!(sim; pace = 1.0)

trace = input_trace(sim)                    #frame → device-tagged slot batches
run!(sim; inputs = trace)                   #bit-identical replay, no devices needed
