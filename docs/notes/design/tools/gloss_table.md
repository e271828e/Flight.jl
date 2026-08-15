# Gloss table — Phase 3 sweep 3 (first-use glosses)

One gloss per Appendix D entry, compressed from that entry's own text.
Column *class*: **A** = applied in the body at section first use (opaque
coinage); **B** = glossary link alone suffices (near-ordinary or
domain-standard vocabulary, or core Part-I terms the spec teaches head-on and
repeats on nearly every page); **—** = entry carries no body link (whitelisted
in `check_glossary.jl`), gloss recorded for completeness only.
*n* = body link count.

## D.1 Component model and declaration layer

| term | anchor | n | class | gloss |
|---|---|---|---|---|
| abstract entry | g-abstract-entry | 3 | A | an `input_types` entry admitting any concrete producer face |
| assembly | g-assembly | 24 | B | a component that only composes children, with no dynamics |
| auto-published port | g-auto-published-port | 3 | A | published by the framework from the state or mode store |
| class | g-class | 5 | A | a component's primitive-vs-assembly status, read off its declarations |
| component | g-component | 67 | B | a leaf primitive or an assembly of components |
| container children | g-container-children | 2 | A | a tuple field contributing its component elements as children |
| continuous component | g-continuous-component | 6 | B | the hybrid primitive: state, modes, flow, stages, events |
| contract | g-contract | 48 | B | a component's declared `input_types` and `output_types` |
| declaration inventory | g-declaration-inventory | 0 | — | the closed set of well-known functions a component defines |
| derived contract | g-derived-contract | 0 | — | the checkable surface an assembly or the `Build` derives instead of declaring |
| function family | g-function-family | 1 | A | which bundle fields a given function may legally receive |
| generic holding | g-generic-holding | 2 | A | a parent holding a child through a non-concrete field type |
| hybrid causal system | g-hybrid-causal-system | 1 | B | continuous flow, discrete dynamics, events and injected inputs |
| the letters | g-the-letters | 0 | — | *(resists compression: the entry is a naming table)* |
| periodic discrete component | g-periodic-discrete-component | 0 | — | a leaf updating at a declared rate, holding between ticks |
| rate scope | g-rate-scope | 3 | A | an assembly's `sample_times` declaration against the enclosing scope |
| schema authority | g-schema-authority | 2 | A | declarations define structure; evaluation only checks conformance |
| stage function / two-stage outputs | g-stage-function | 4 | A | `h_x` or `h_xu`, the two output stages every component provides |
| workspace | g-workspace | 6 | A | component-declared mutable scratch arriving as the `ws` bundle field |

## D.2 Signals and data homes

| term | anchor | n | class | gloss |
|---|---|---|---|---|
| buffer | g-buffer | 20 | B | the framework-owned flat vector backing all continuous state |
| bundle | g-bundle | 18 | A | the NamedTuple of zero-copy views a component function receives |
| cell | g-cell | 45 | B | one typed entry of the signal table, one per port |
| constant source | g-constant-source | 1 | B | a library component publishing a value its instance holds |
| entry | g-entry | 0 | — | *(resists compression: a disambiguation entry, never used bare)* |
| face | g-face | 42 | B | the name a port wears on its component's boundary |
| feedthrough | g-feedthrough | 14 | B | an instantaneous input→output dependence |
| field handle / function-valued signal | g-field-handle | 4 | A | an immutable query object consumers evaluate at their own arguments |
| immutable value semantics | g-immutable-value-semantics | 0 | — | immutability plus frozen references, so concurrent reads are safe |
| one home per datum | g-one-home-per-datum | 2 | A | each datum lives in exactly one store, mirrored nowhere |
| port | g-port | 42 | B | the addressable unit of the model: one name, one cell |
| signal table | g-signal-table | 12 | B | the framework-owned collection of cells holding every produced signal |
| slot | g-slot | 46 | B | a root input slot, the only thing the periphery writes |
| staging cell | g-staging-cell | 8 | A | where a device's pending write batch waits between drains |
| store | g-store | 9 | B | the typed home of `m` and a discrete leaf's `x` |
| summing junction | g-summing-junction | 2 | B | a library component doing N-to-1 aggregation through explicit wires |
| value-level constructor | g-value-level-constructor | 4 | A | the plain exported function building a field handle from the component and input values |
| view | g-view | 4 | B | a zero-copy reconstruction of a store handed through the bundle |
| `w` (private intermediates) | g-w | 0 | — | the optional private second slot of a stage's return |

## D.3 Evaluation and scheduling

| term | anchor | n | class | gloss |
|---|---|---|---|---|
| algebraic loop | g-algebraic-loop | 3 | B | a genuine cycle in the instantaneous dependency graph |
| flow / RHS | g-flow | 12 | B | `f`, the continuous derivative function |
| frame | g-frame | 6 | A | one iteration of the loop: drain, integrate, boundary, publication |
| projection | g-projection | 11 | B | the optional per-component hook `x ← project(x)` |
| schedule | g-schedule | 19 | B | the static evaluation order computed once at build time |
| sweep | g-sweep | 44 | B | one execution of the schedule against the current state |

## D.4 Time and events

| term | anchor | n | class | gloss |
|---|---|---|---|---|
| anchor | g-anchor | 2 | A | the exact `(T, τ)` pair an `Absolute` entry establishes |
| bound schedule | g-bound-schedule | 3 | A | the printable per-component `(D, Φ, Δt)` artifact deployment binding produces |
| boundary | g-boundary | 65 | B | a published consistency point where a snapshot goes out |
| boundary-detected | g-boundary-detected | 8 | A | checked for edges at step boundaries only, no root-finding |
| chattering / localization budget | g-chattering | 4 | A | the bounded per-frame localization allowance, exhaustion degrading rather than throwing |
| `Δt_base` | g-dt_base | 1 | A | the base tick period, an integer multiple `n·h` |
| due | g-due | 6 | A | admitted at this boundary by its compiled `(D, Φ)` pair |
| edge semantics / holding | g-edge-semantics | 7 | A | firing on not-holding → holding transitions, never bare sign changes |
| firing budget | g-firing-budget | 2 | A | the per-boundary cap on how often each event fires |
| guard | g-guard | 30 | B | the declared function defining an event's predicate |
| harmonic grid | g-harmonic-grid | 4 | A | every discrete period an integer multiple of `Δt_base` |
| input epoch | g-input-epoch | 3 | A | a maximal span of constant `u`, delimited by frame-top drains |
| interpolant | g-interpolant | 3 | A | the cubic Hermite continuous extension over the last completed step |
| localized | g-localized | 8 | A | the crossing instant bracketed by root-finding over trial sweeps |
| pacing / pacer debt | g-pacing | 8 | A | waits inserted between completed frames, never altering the boundary sequence |
| phase (`Φ`) | g-phase | 3 | A | a schedule's offset against its grid, compiled to base ticks |
| predicate | g-predicate | 14 | B | what a guard defines: a `Bool` form or a sign |
| prior | g-prior | 5 | A | the event's stored predicate sample from the previous boundary |
| quiescence | g-quiescence | 5 | A | the fixed point where a round of handlers fires nothing |
| remainder step | g-remainder-step | 1 | A | the integration from `t*` to the original grid target |
| `t*` | g-t | 0 | — | the localized event time, structurally strictly later than `tₙ` |
| tick | g-tick | 27 | B | an instant at which a discrete component's stages run |
| tier | g-tier | 29 | B | the continuous or discrete side of the hybrid formalism |

## D.5 Build pipeline

| term | anchor | n | class | gloss |
|---|---|---|---|---|
| activation | g-activation | 30 | A | a re-run of Stratum C at a given scalar type |
| always-on conformance check | g-always-on-conformance-check | 0 | — | one type test of a stage return at the table-write point |
| `Build` | g-build | 1 | B | the artifact `build(world)` produces: wires, faces, schedule, slots |
| chunking | g-chunking | 2 | A | splitting a large phase body into statically typed chunks |
| executable set | g-executable-set | 1 | A | the function set an activation can actually run, hence probes |
| executor | g-executor | 10 | A | the compiled execution form of the schedule |
| leaf walk | g-leaf-walk | 1 | A | the derivation of per-activation types from a declared nominal type |
| lens (`Getter`) | g-lens | 2 | A | the compiled navigation step of a condition entry |
| measurement seam / phase bodies | g-measurement-seam | 2 | A | `phase_bodies(sim)`, the compiled bodies bound over the simulation's buffers |
| nominal | g-nominal | 1 | A | the `Float64` activation, and a declaration's `Float64` face |
| probe | g-probe | 23 | B | the build's single evaluation of a user function with real values |
| probe value / input synthesis | g-probe-value | 2 | A | fabricated build-time values, synthesized at producerless root slots and flowing the probe chain |
| `ProbeDual` | g-probedual | 1 | B | the exported canonical concrete probe scalar |
| schema vs. layout | g-schema-vs-layout | 0 | — | *(resists compression: the entry contrasts a pair)* |
| stratum | g-stratum | 17 | A | one of the build's three phases: structure, schedule, activation |
| walked / pinned / exempt | g-walked | 21 | A | the eltype-genericity classes: follow the activation scalar, stay `Float64`, exempt |

## D.6 Runtime periphery

| term | anchor | n | class | gloss |
|---|---|---|---|---|
| bad datum | g-bad-datum | 1 | A | a datum unmappable for environmental reasons, tolerated not thrown |
| batch | g-batch | 6 | B | a device's staged set of face ⇒ value writes |
| binding | g-binding | 4 | B | the value passed at `attach!` that makes a device framework-legible |
| boundary counter | g-boundary-counter | 1 | B | the monotonic count of published boundaries |
| calling task | g-calling-task | 8 | A | the task that invoked `run!` |
| claim | g-claim | 14 | B | the set of faces a device may write |
| coalescing | g-coalescing | 5 | A | the CAS merge keeping one pending batch per device |
| control plane | g-control-plane | 6 | A | the separate atomic surface carrying pause, pace and stop |
| derived liveness | g-derived-liveness | 1 | A | widget liveness derived from the feed chain, never marked per port |
| device | g-device | 39 | B | any attached participant in the periphery |
| diagnostic cell | g-diagnostic-cell | 5 | A | the single-writer ring each writer owns for diagnostics and heartbeat |
| drain | g-drain | 20 | A | the frame-top swap that publishes staged device inputs into the root slots |
| framework status | g-framework-status | 5 | A | the frozen diagnostics value each snapshot carries beside the table |
| greedy claim | g-greedy-claim | 3 | A | the unclaimed complement, computed by the framework instead of returned |
| harness cell | g-harness-cell | 3 | A | the always-present staging cell of the harness register |
| harness register | g-harness-register | 2 | A | the framework-owned `stage!` write path of the calling task |
| `latest` | g-latest | 1 | B | the atomic reference a published snapshot is stored into |
| next-snapshot wait | g-next-snapshot-wait | 0 | — | the boundary counter plus one condition variable |
| operator interrupt | g-operator-interrupt | 3 | A | Ctrl-C read as a control-plane stop rather than a failure |
| orphaned claims | g-orphaned-claims | 1 | A | the claims of a device whose task died mid-run |
| peek | g-peek | 6 | A | showing a widget's own pending write, else the snapshot value |
| periphery | g-periphery | 12 | B | everything outside the loop that exchanges data with it |
| roster | g-roster | 17 | B | the list of attached device entries, read once at `run!` |
| scenario component | g-scenario-component | 4 | A | an ordinary periodic discrete component holding a sim-time script |
| selector (read-selector family) | g-selector | 6 | A | the closed family of deferred reads resolving against a source |
| `should_abort` | g-should_abort | 1 | B | whether a device's departure also requests a stop |
| snapshot | g-snapshot | 30 | B | the immutable per-boundary publication of the signal table |
| stage-on-interaction | g-stage-on-interaction | 2 | A | widgets stage on edit or activation, never per render pass |
| unattended run | g-unattended-run | 6 | A | a run with empty staging and no snapshot readers |
| write surface | g-write-surface | 2 | A | the set of faces a writer's batch entries may reach |

## D.7 Recording and replay

| term | anchor | n | class | gloss |
|---|---|---|---|---|
| decimation | g-decimation | 2 | A | the log's keep-every-kth retention policy |
| frame ordinal | g-frame-ordinal | 1 | A | the trace's key, the frame index a batch replays at |
| log | g-log | 1 | B | the retained sequence of published snapshots |
| recorders | g-recorders | 2 | B | the trace and the log jointly |
| replay | g-replay | 23 | B | the ordinary loop re-driven from the trace |
| run metadata | g-run-metadata | 1 | B | the trace header's deployment block |
| trace | g-trace | 29 | B | the primary record: drained, device-tagged batches per frame |
| trace header | g-trace-header | 13 | B | the trace's preamble: initial stores, slot values, schemas, deployment |
| trace record | g-trace-record | 0 | — | the retained form of a drained batch |
| what-if register | g-what-if-register | 2 | A | replaying a trace against the same structure with changed parameters |

## D.8 Stopped-sim services and the condition algebra

| term | anchor | n | class | gloss |
|---|---|---|---|---|
| `at` / `Scoped` | g-at | 1 | A | the scoping combinator storing a path prefix beside a node |
| baseline | g-baseline | 2 | B | an aircraft-shipped, full-coverage condition function |
| boundary zero | g-boundary-zero | 18 | A | the initialization boundary: the ordinary macro-sequence with an empty integrate |
| capture | g-capture | 1 | A | reading the current stores and slots back as a condition |
| component test rig | g-component-test-rig | 2 | B | a one-child assembly exporting the child's whole input face set |
| condition | g-condition | 12 | A | the path-addressed sparse overlay that sets a build's state |
| `design_world` | g-design_world | 1 | B | the shipped thin world that mounts an aircraft |
| fragment | g-fragment | 4 | B | the leaf node of the condition algebra |
| fragment tree | g-fragment-tree | 0 | — | the inert, lazy composition of condition nodes |
| merge | g-merge | 1 | B | the symmetric, collision-intolerant combinator over condition nodes |
| mounting | g-mounting | 2 | A | relocating a whole problem or tap set with `at(prefix, …)` |
| override | g-override | 1 | B | the ordered, asymmetric layering combinator: the patch wins |
| service lifecycle | g-service-lifecycle | 0 | — | the `Simulation` states and each service's legality against them |
| slot totality | g-slot-totality | 4 | A | the requirement that an application establishing a complete world cover every root slot |
| taps | g-taps | 3 | A | the three selector lists declaring what linearization seeds and reports |
| `TrimProblem` | g-trimproblem | 1 | B | the closed seven-field value specifying a trim |

## D.9 Error discipline and diagnostics

| term | anchor | n | class | gloss |
|---|---|---|---|---|
| carrier exception | g-carrier-exception | 0 | — | the single exception a set of diagnostics travels in |
| collect the checks, fail the evaluations fast | g-collect-the-checks-fail-the-evaluations-fast | 0 | — | declarative passes collect; the first user-code exception aborts |
| did-you-mean | g-did-you-mean | 16 | A | the offending name plus the list-in-hand it should have matched |
| error locality | g-error-locality | 1 | A | a mistake fails at the site of the mistake |
| execution cursor | g-execution-cursor | 1 | A | the mutable field recording where in the schedule execution is |
| feedthrough tracer | g-feedthrough-tracer | 1 | A | the set-propagation instrument classifying a rejected cycle |
| kind | g-kind | 3 | B | a diagnostic's identity in the closed Appendix C set |
| payload | g-payload | 8 | B | the structured data a diagnostic carries beside its kind |
| `stop_on` / termination is a state | g-stop_on | 1 | B | the deployment policy naming the faces the loop reads |
| warning streams | g-warning-streams | 1 | B | two streams, build and runtime, scoped separately |

## D.10 Meta-vocabulary

| term | anchor | n | class | gloss |
|---|---|---|---|---|
| blessed | g-blessed | 4 | A | the spec's marker for a practice it explicitly sanctions |
| the freeze | g-the-freeze | 4 | A | the roster freeze: `attach!`/`detach!` are stopped-sim operations |
| guarded addition | g-guarded-addition | 8 | A | a capability the design admits but does not build |
| normative / index, not a second home | g-normative | 1 | B | the spec norms; its appendices are indices |
| recorded, not built | g-recorded-not-built | 3 | A | a worked-out extension deliberately left unimplemented, its seams named |
| register | g-register | 10 | B | the spec's word for a mode or idiom, always compounded |
| row | g-row | 0 | — | a numbered entry of `framework_decisions.md` |
| seam | g-seam | 19 | B | a narrow, named interface kept deliberately thin |
| torture test | g-torture-test | 1 | B | an awkward existing artifact transliterated to validate a mechanism |
| worked (example) | g-worked | 4 | B | a full spelling of a mechanism against a real artifact |

## Application record (regenerated at the gloss ratification sitting, 2026-08-15)

179 inline gloss applications standing in the body (chapters 1-16): Part I 27,
Part II 60, Part III 34, Part IV 52, Part V 6. 64 of the 81 class-A terms are
applied at least once.

An *application* is a gloss parenthetical or appositive (comma, dash or colon)
attached to the term, carrying the gloss above or a close derivative of it. A
bare glossary link is not one, and neither is a copular definition sentence
("an X is ..."), which is the section's own definitional prose. Rows follow the
Appendix D grouping above.

| gloss | sites |
|---|---|
| abstract entry | §13.7 |
| auto-published port | §7.1, §14.4 |
| class | §11.2, §11.5, §12.1 |
| function family | §5.2 |
| generic holding | §11.5 |
| rate scope | §12.2 |
| schema authority | §11.3, §13.7 |
| stage function | §4.1, §4.3, §4.4, §11.1 |
| workspace | §7.5, §11.2, §12.1, §14.1 |
| bundle | §5.3, §7.3, §8.5, §8.6, §10.5, §12.3, §13.2, §15.5 |
| field handle | §7.5, §9.2 |
| staging cell | §10.7, §14.8 |
| value-level constructor | §14.1 |
| frame | §9.1 |
| anchor | §12.1 |
| bound schedule | §8.5, §12.1, §12.2 |
| boundary-detected | §2.2, §8.4, §10.5, §11.2, §16 |
| `Δt_base` | §5.5 |
| due | §12.7 |
| edge semantics | §8.6, §14.5 |
| firing budget | §5.3 |
| harmonic grid | §8.4, §14.5, §14.8 |
| interpolant | §8.4 |
| localized | §8.2, §8.4, §10.7, §11.2, §13.5 |
| pacing | §8.1, §8.5, §10.7 |
| phase (`Φ`) | §8.5 |
| prior | §8.4, §8.6, §14.5 |
| quiescence | §5.3, §8.4, §8.5, §8.6 |
| activation | §5.6, §6.1, §6.2, §7.3, §11.2, §11.3, §12.7, §13.7, §14.3, §14.8, §14.10 |
| chunking | §11.1 |
| executor | §5.2, §7.1, §11.1, §11.5, §12.5, §13.4 |
| leaf walk | §11.2 |
| lens | §14.3 |
| nominal | §6.1 |
| stratum | §5.6, §6.1, §11.2, §11.5, §13.1, §13.2, §13.3, §14.3 |
| calling task | §9.1, §9.2, §9.3, §9.6, §10.6 |
| coalescing | §9.4 |
| control plane | §10.4 |
| diagnostic cell | §10.4, §13.2 |
| drain | §8.4, §9.4, §9.5, §9.7, §14, §15.3 |
| framework status | §8.7, §10.2 |
| greedy claim | §15.4 |
| harness cell | §9.3, §9.4, §10.6 |
| peek | §9.2, §15.4 |
| scenario component | §10.5 |
| selector | §9.2, §14.7, §14.10 |
| stage-on-interaction | §14.5 |
| unattended run | §7.5, §9.8, §10.4, §13.4 |
| write surface | §9.3, §11.6 |
| decimation | §9.8 |
| frame ordinal | §10.7 |
| what-if register | §9.5, §10.7 |
| boundary zero | §8.4, §8.5, §10.4, §10.6, §10.7, §13.4, §14.4, §14.6, §14.8, §14.9, §14.10, §16 |
| capture | §14 |
| condition | §11.2, §13.3, §13.5, §14.4, §14.5 |
| mounting | §11.5, §14.4 |
| slot totality | §14 |
| taps | §14.4, §14.10 |
| did-you-mean | §5.2, §9.2, §10.7, §11.3, §11.4, §11.5, §12.1, §13.1, §13.2, §13.3, §14.3, §14.10 |
| error locality | §11.1 |
| execution cursor | §13.4 |
| feedthrough tracer | §7.2 |
| guarded addition | §4.3, §7.5, §9.3, §10.2 |
| recorded, not built | §14.7, §14.9, §14.10 |

The seventeen class-A terms with no application anywhere in the body — each
either defined copularly at its owning section and used bare thereafter, or
carrying only a differently-worded expansion: `at`/`Scoped`, bad datum,
blessed, chattering, container children, derived liveness, executable set,
harness register, input epoch, measurement seam, one home per datum, operator
interrupt, orphaned claims, probe value, remainder step, the freeze,
walked/pinned/exempt.
