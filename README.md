# TRIDENT – 2D Truss FEA Solver

MATLAB finite element solver for planar trusses. It performs a linear static analysis of a 3-bar concurrent truss under an inclined load: global stiffness assembly, boundary conditions by elimination, nodal displacements, support reactions and deformed-shape plot.

---

## Problem

Three bars meet at a single free node **A**. Nodes **B**, **C** and **D** are pinned. A load **P** is applied at A with an angle **α** to the x-axis.

```
        D (0, L)
          o
          |  \
          |    \   bar 3
          |      \
 C (0, 0) o--------o A (L, 0)   ← P, angle α
          |  bar 2 /
          |      /
          |    /   bar 1
          |  /
          o
        B (0, -L)
```

| Element | Nodes | Length |
|---------|-------|--------|
| Bar 1   | A – B | L√2    |
| Bar 2   | A – C | L      |
| Bar 3   | A – D | L√2    |

---

## Method

- **Element**: 2-node bar (truss) element, axial stiffness only, small displacements, linear elastic material.
- **Element stiffness** in global axes, with `c = cos θ` and `s = sin θ`:

```
            EA  [  c²   cs  -c²  -cs ]
   Ke  =   ---- [  cs   s²  -cs  -s² ]
             L  [ -c²  -cs   c²   cs ]
                [ -cs  -s²   cs   s² ]
```

- **Assembly** of the global stiffness matrix `K` (2 DOF per node).
- **Boundary conditions**: fixed DOFs are removed (row/column elimination), and imposed displacements are accounted for through a correction term on the load vector.
- **Solution**: `K_red · q_red = F_red`
- **Post-processing**: reactions are computed as `R = K · q`.

---

## Usage

1. Open `trident_truss_fea.m` in MATLAB (R2016b or later; no toolbox required).
2. Edit the input data at the top of the script if needed.
3. Run the script (`F5`).

### Input parameters

| Variable  | Description                        | Default  |
|-----------|------------------------------------|----------|
| `L`       | Reference length (m)               | 1        |
| `P`       | Applied load magnitude (N)         | 1000     |
| `alpha`   | Load angle (rad)                   | π/4      |
| `Module`  | Young's modulus E per bar (Pa)     | 200e9    |
| `Section` | Cross-section A per bar (m²)       | 1e-4     |
| `Coord`   | Nodal coordinates `[x, y]`         | see above|
| `Connec`  | Element connectivity `[n1, n2]`    | see above|

Geometry, materials, loads and boundary conditions are all defined as data. The solver therefore extends to any 2D truss: change `Nn`, `Nr`, `Coord`, `Connec` and the BC/load vectors.

---

## Outputs

- Element stiffness matrices `Ke` with length and direction cosines
- Global stiffness matrix `K` and reduced matrix `K_red`
- Nodal displacements `(Ux, Uy)`
- Reaction forces at supports
- Plot of the initial and deformed configurations (amplification ×1000)

### Expected results (default parameters)

Stiffness at node A: `Kxx = EA(1 + 1/√2)/L` and `Kyy = EA/(L√2)`, with the coupling terms cancelling by symmetry.

| DOF   | Value              |
|-------|--------------------|
| Ux(A) | ≈ 2.071 × 10⁻⁵ m   |
| Uy(A) | ≈ 5.000 × 10⁻⁵ m   |

These values provide an analytical check of the solver.

---

## Assumptions and limitations

- Linear elasticity, small displacements (geometric linearity)
- Bars carry axial load only (pin-jointed)
- Static loading; no buckling check
- 2D only

---

## Author

**Nathan Fagnou**, M1 Modélisation et Mécanique Numérique, Aerospace & Defense track
