# QuantumControl.jl

```@eval
using Markdown
using Pkg

VERSION = Pkg.dependencies()[Base.UUID("8a270532-f23f-47a8-83a9-b33d10cad486")].version

github_badge = "[![Github](https://img.shields.io/badge/JuliaQuantumControl-QuantumControl.jl-blue.svg?logo=github)](https://github.com/JuliaQuantumControl/QuantumControl.jl)"

version_badge = "![v$VERSION](https://img.shields.io/badge/version-v$VERSION-green.svg)"

Markdown.parse("$github_badge $version_badge")
```

[QuantumControl.jl](https://github.com/JuliaQuantumControl/QuantumControl.jl) is a [Julia framework for quantum dynamics and control](https://github.com/JuliaQuantumControl).

Quantum optimal control [BrumerShapiro2003,BrifNJP2010,Shapiro2012,KochJPCM2016,SolaAAMOP2018,MorzhinRMS2019,Wilhelm2003.10132,KochEPJQT2022](@cite) attempts to steer a quantum system in some desired way by finding optimal control parameters or control fields inside the system Hamiltonian or Liouvillian. Typical control tasks are the preparation of a specific quantum state or the realization of a logical gate in a quantum computer (["pulse level control"](https://arxiv.org/abs/2004.06755)). Thus, quantum control theory is a critical part of realizing quantum technologies at the lowest level. Numerical methods of *open-loop* quantum control (methods that do not involve measurement feedback from a physical quantum device) such as [Krotov's method](https://github.com/JuliaQuantumControl/Krotov.jl) [Krotov1996,SomloiCP1993,BartanaJCP1997,PalaoPRA2003,ReichJCP2012,GoerzSPP2019](@cite) and [GRAPE](https://github.com/JuliaQuantumControl/GRAPE.jl) [KhanejaJMR2005,FouquieresJMR2011](@cite) address the control problem by [simulating the dynamics of the system](https://github.com/JuliaQuantumControl/QuantumPropagators.jl) and then iteratively improving the value of a functional that encodes the desired outcome.


##### Functional

Mathematically, the control problem is the minimization of a functional of the form

```math
J(\{ϵ_l(t)\})
    = J_T(\{|Ψ_k(T)⟩\})
    + λ_a \, \underbrace{∑_l \int_{0}^{T} g_a(ϵ_l(t)) \, dt}_{=J_a(\{ϵ_l(t)\})}
    + λ_b \, \underbrace{∑_k \int_{0}^{T} g_b(|Ψ_k(t)⟩) \, dt}_{=J_b(\{|Ψ_k(t)⟩\})}\,,
```

where ``\{ϵ_l(t)\}`` is a set of [control functions](@ref "Control Function") defined between the initial time ``t=0`` and the final time ``t=T``, and ``\{|Ψ_k(t)⟩\}`` is a set of ["trajectories"](@ref QuantumControl.Trajectory) evolving from a set of initial states ``\{|\Psi_k(t=0)⟩\}`` under the controls ``\{ϵ_l(t)\}``. The full functional consists of the final-time functional ``J_T``, a control-dependent running cost ``J_a`` weighted by ``λ_a``, and sometimes a state-dependent running cost ``J_b`` weighted by ``λ_b``. The states can be Hilbert space vectors or density matrices, and the equation of motion is implicit, typically the Schrödinger or Liouville equation.

The `QuantumControl.jl` package provides a single coherent [API](@ref QuantumControlAPI) for solving the quantum control problem with the [packages](https://github.com/JuliaQuantumControl#packages) in the [JuliaQuantumControl](https://github.com/JuliaQuantumControl) organization. Different [optimization methods](@ref "Control Methods") target specific variants of the above functional.


## Getting Started

* See the [installation instructions](https://github.com/JuliaQuantumControl/QuantumControl.jl#installation) on Github.

* Look at a [simple example for a state-to-state transition with GRAPE](@extref Examples :doc:`examples/simple_state_to_state/index`) to get a feeling for how the `QuantumControl` package is intended to be used, or look at the larger list of [Examples](@extref Examples :doc:`index`).

* Read the [Glossary](@ref) and [Overview](@ref) to understand the philosophy of the framework.

## Contents

```@contents
Pages = [
    "glossary.md",
    "overview.md",
    "methods.md",
    "howto.md",
]
Depth = 2
```

### Examples

```@contents
Pages = [
    "examples/index.md",
]
```

### API

```@contents
Pages = [
    "api/quantum_control.md",
]
Depth = 1
```

#### Sub-Packages

```@contents
Pages = [
    "api/quantum_propagators.md",
]
Depth = 1
```

### History

See the [Releases](https://github.com/JuliaQuantumControl/QuantumControl.jl/releases) on Github.
