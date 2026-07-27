#User-level listing for the §12 runtime periphery on the settled machinery
#(framework_design.md v0.16): root slots as exported input faces (§13.6), slot
#exclusivity and declarative bindings (§12.3, §17.4), derived GUI liveness and
#stage-on-interaction (§12.5), slot totality at init (§16.6), always-on trace
#and same-build replay (§12.3). Illustrative, non-committed syntax; not
#runnable.
#
#Lineage note: this file's original scenario — joystick and GUI writing u_cmd
#simultaneously, adjudicated by attachment order — was the §17.3 staging
#torture test. Slot exclusivity (§12.3, adjudicated by the §17.4 cast) made
#that configuration build-rejected by design: one writer per slot, claimed at
#attach. The listing now shows the settled behavior instead: the joystick
#claims u_cmd, the GUI's u_cmd widget renders read-only by derived liveness,
#and its τ widget stays live.

######################## Model #################################################

struct LowPassFilter <: AbstractComponent end

init_x(::LowPassFilter) = (x = 0.0,)
#no z, no modes, no events — the minimal continuous primitive

inputs(::LowPassFilter) = (u_cmd = Float64, τ = Float64)

#outputs is mandatory even here (§13.2): x names a state field no stage
#produces → auto-published at stage-1 position
outputs(::LowPassFilter, ::Type{T}) where {T<:Real} = (x = T,)

#no g_s1/g_s2: no published intermediates → the law lives directly in f
#(no duplicate computation site to drift, §9.4)
f(::LowPassFilter, x, m, y, u, t) = (x = (u.u_cmd - x.x) / u.τ,)

######################## Root assembly #########################################

#The root is an ordinary assembly (§13.5). At the root, an exported input
#face IS the write surface: "u_cmd" and "τ" become root slots with no extra
#vocabulary (§13.6) — the periphery's entire write surface, drained at frame
#top, in every snapshot (§12.3).
struct FilterRig <: AbstractComponent
    filt::LowPassFilter
end

connections(::FilterRig) = ()   #mandatory even when empty: the kind marker (§13.5)

exports(::FilterRig) = (
    "u_cmd" => "filt/u_cmd",    #input faces → root slots
    "τ"     => "filt/τ",        #runtime-adjustable → a slot, not a parameter
    "x"     => "filt/x",        #output face: the observation contract
)

######################## Joystick device #######################################

struct ToyStick <: IODevice          #one device kind (§12.4); ~today's IODevices
    period::Float64                  #polling cadence, wall-clock seconds
end

init!(js::ToyStick)     = nothing    #open SDL, claim hardware, ...
shutdown!(js::ToyStick) = nothing

#runs on the device task each cycle; may block or pace itself — its business
poll(js::ToyStick) = (axis = read_hw_axis(js),)

#The binding is declarative data, not shaping code (§17.4): device data key →
#face name + conditioning params. Conditioning (expo, deadzone) is device
#truth and lives here; command semantics stay in-model — the face carries
#post-conditioning, writer-independent meaning (§12.3). Another stick model
#is the same table with different keys, zero shaping code.
stick_binding = (axis = (face = "u_cmd", expo = 0.5, deadzone = 0.05),)

#framework-owned device task loop, for reference (the user never writes it):
#   init! → [ data = poll(dev); stage!(cell, condition(data, binding)) ]* → shutdown!
#Writes are levels, never deltas; the batch overwrites the device's own cell
#(inter-frame polls coalesce, ZOH-correct); the frame-top drain applies and
#traces it (§12.3).

######################## GUI panel #############################################

#Per-component extension (§12.5), reusable wherever a LowPassFilter appears.
#Widgets name the component's OWN PORTS, never root slots; the build resolves
#each port's feed chain through wires and exports. Liveness is derived, not
#configured: live iff the chain ends in a root slot that is currently
#unclaimed. Here, with the joystick attached, u_cmd is claimed → its slider
#renders read-only with provenance; τ is unclaimed → live. Unplug the stick
#→ the claim releases → the slider goes live at the last-held slot value.
function GUI.draw!(ctx, ::LowPassFilter)
    s = view(ctx)                            #snapshot view scoped to this instance

    GUI.text("t = $(s.t)")
    GUI.display_bar("state x",       s.x)
    GUI.display_bar("applied u_cmd", s.u_cmd)   #what the sim integrated; root
                                                #slots are in the snapshot (§12.2)

    GUI.input_slider!(ctx, :u_cmd, -1.0, 1.0; label = "commanded input")
    GUI.input_slider!(ctx, :τ,      0.1, 10.0; label = "time constant")

    GUI.pause_button!(ctx)                   #control plane, not staging (§12.6)
end

#Staging contract (§12.5): widgets stage on interaction events only — the τ
#slider stages the new level on edit, idempotent under the levels doctrine;
#no stage-every-pass. The slider knob shows own-pending-else-snapshot (peek
#rule); while paused, staged edits display indefinitely and apply at the
#un-pause drain.

######################## Setup and run #########################################

world = FilterRig(LowPassFilter())

sim = Simulation(world; algorithm = Heun(), h = 0.02, t_end = 60.0)
#the entire build pipeline runs here (§14): kind resolution, wiring and
#contract checks, face derivation, schedule, layouts, slot table

#Slot totality (§16.6): every root slot gets a condition value before
#boundary zero, or init! errors with the uncovered faces — no synthesized
#defaults, no third branch. t₀ is an init-service argument, not condition
#data (§16.5).
init!(sim, fragment(slots = (u_cmd = 0.0, τ = 1.0)); t0 = 0.0)

attach!(sim, ToyStick(0.01), stick_binding)
#at attach: faces resolved against the root contract (typo → did-you-mean),
#claim on "u_cmd" registered — a second writer on the same face errors HERE,
#not at runtime (§12.3 exclusivity; supersedes attachment-order precedence)

run!(sim; gui = true, pace = 1.0)            #50 Hz boundaries = drain/publish cadence

######################## Replay ################################################

trace = input_trace(sim)        #always-on (§12.7): header = initial (x, m, z)
                                #+ slot values; body = per-frame device-tagged
                                #slot batches — primary data, the log is
                                #recomputable from it

sim2 = Simulation(world; algorithm = Heun(), h = 0.02, t_end = 60.0)
run!(sim2; inputs = trace)      #bit-identical replay on the same build, no
                                #devices attached; also the state-trajectory
                                #inspector (snapshots carry no state stores)
