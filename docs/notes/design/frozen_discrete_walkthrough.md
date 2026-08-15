# The frozen discrete tier: zero partials as the exact answer

*A companion explainer, not normative text. The ground truth is
`framework_spec.md` [§7.2][s7-2] (genericity scoping), [§8.5][s8-5] (gated tick scheduling),
[§11.2][s11-2] (the frozen-exact typing rule), [§12.5][s12-5] (the embedding guarantee) and
[§14.10][s14-10] (linearization and the sampled-data door), with decision D-070 and D-072
scoping the AD consumers. If this document and the spec ever disagree, the
spec wins. It was written during the 2026-08 declaration-signature
re-adjudication, to preserve the answer to one recurring question.*

Everything here answers one question. Take a continuous component A, a
discrete component B, and a continuous component C, wired in a chain: B's
`h_xu` reads one of A's outputs and computes an output wired to one of C's
inputs. Under a `Dual` activation, B's stages never run and B's cells are
frozen `Float64` constants with zero partials. **How can the AD chain work
without calling `h_xu` between A's and C's output stages?**

The answer: the instantaneous chain A → B → C that the question pictures does
not exist in the real system either. AD follows actual dataflow, and the
actual dataflow through B is temporal, not instantaneous. The frozen cell is
not an AD limitation on the signal path — it is the true zero of an
instantaneous dependence that the hybrid semantics never had.

## 1. What B does between ticks

B is a sampled-data component — the model of a digital device. Its execution
contract ([§8.5][s8-5]'s gated scheduling) is: at B's tick times $t_k$, its stages
run — `h_x`/`h_xu` read whatever their inputs *are at that instant*, compute
B's outputs, write them to B's cells, and `g` updates B's `x`. Between ticks, B's
stages do not run at all; its cells **hold** — the zero-order hold (ZOH) that
a real DAC performs. This is true in every nominal simulation sweep, before
AD ever enters the picture.

## 2. A timeline

Let B tick at 50 Hz, so $t_k = 0.02\,\mathrm{s}$, and consider an integrator
minor stage evaluating the model at $t = 0.023\,\mathrm{s}$:

- **A's output** is computed fresh from $x(0.023)$ — a live function of the
  current state.
- **B's cell** contains a *number*: the value computed at $t = 0.02$ from A's
  output at $t = 0.02$. The inputs that produced it were consumed then; they
  are gone.
- **C** reads that number.

C never sees a "current" B output between ticks. The signal chain A → B → C
exists only *at tick instants*; between them, C's input is a constant
determined by past samples.

## 3. The derivative of the execution story

Today's AD consumers — linearization (D-072) and trim (D-070) —
differentiate the instantaneous vector field $\dot{x} = f(x, \cdot)$ at a
frozen instant $t$, with the discrete state (and B's held cells, which are the same kind of
fact: sampled state) held. Seed $x$. What is the true sensitivity of C's
input to the seed? **Exactly zero** — C's input is a held constant whose
value was determined by past samples, and no perturbation of the present $x$
reaches back in time to change it.

So the `Dual` sweep that runs only the continuous chain, with B's cell
embedding as a zero-partial constant ([§12.5][s12-5]), is not skipping part of the
chain — it is computing the exact derivative of a dependence structure in
which B contributes *no dependence at this instant*. Zero partials is not an
approximation of a derivative too hard to get; it is the derivative.

## 4. The fictitious system

Suppose instead the sweep ran B's `h_xu` on A's current `Dual` output. The
partials would flow A → B → C, and the result would be the Jacobian of a
**fictitious system** — one in which B recomputes continuously, i.e. the
sample-and-hold has been deleted and the digital controller replaced by a
memoryless analog one. That Jacobian is *wrong* for the hybrid system being
simulated: the hold is not an implementation detail, it is dynamics — it is
why discretization matters in control design at all.

There is a structural echo of the same fact. A feedback loop closed through B
is schedulable precisely because B's output is held state rather than
feedthrough ([§5.3][s5-3]); treat `h_xu` as in-sweep feedthrough and that loop
becomes algebraic, which [§5.5][s5-5] would have to reject. The freeze and the
schedule are the same statement made twice.

## 5. The dependence that does exist

The intuition behind the question is not wrong — it is aimed at a different
object. $x(t_k)$ genuinely influences $x(t_{k+1})$ through B: A's output at
$t_k$ enters B's tick computation, and B's held output shapes $\dot{x}$ over
the whole interval $[t_k, t_{k+1})$. But that dependence lives in the
composition *across* the interval — the sampled-data step map

$$\Phi : (x_k, \mathrm{slots}) \mapsto x_{k+1}
\qquad \text{(over the whole state, continuous and discrete leaves alike:
integrate one period, then run the due ticks)}$$

and differentiating $\Phi$ is exactly where "call `h_xu` between A and C"
becomes correct: the `Dual`s must flow through B's stages *at tick position
in the step*, with the real leaves of B's `x` walked. That is verbatim [§14.10][s14-10]'s
recorded door — the sampled-data `Dual` activation, executable set
"continuous chain + `f` + discrete `h_x`/`h_xu` + `g`", with its honest
boundary ($\Phi$ is differentiable only where the event pattern is locally
constant; exactness across a firing needs saltation corrections).

## 6. The control-practice mapping

The frozen sweep answers "what is the plant Jacobian with the digital side
held" — it is what you use to obtain the continuous plant model that you then
discretize (`c2d`) and design against. The $\Phi$ Jacobian answers "what is
the closed-loop sensitivity through the controller" — the discrete-time
closed-loop analysis. They are different mathematical questions with
different machinery, not a shorter and a longer version of the same sweep.
Conflating them — differentiating "through" a ZOH as if it were a wire — is
precisely the error the type-level freeze makes unwritable.

<!-- citation link definitions — generated by tools/linkify.jl; do not edit -->
[s11-2]: framework_spec.md#112-the-declaration-inventory
[s12-5]: framework_spec.md#125-the-always-on-conformance-check
[s14-10]: framework_spec.md#1410-linearization-tap-selectors-one-seeded-pass-a-pure-query
[s5-3]: framework_spec.md#53-structural-feedthrough-stage-roles-schedule-and-step-boundaries
[s5-5]: framework_spec.md#55-algebraic-loop-policy-reject-at-build-time
[s7-2]: framework_spec.md#72-numeric-genericity-eltype
[s8-5]: framework_spec.md#85-multi-rate-tick-scheduling
