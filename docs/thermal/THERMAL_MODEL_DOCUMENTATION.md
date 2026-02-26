# Thermal Model Documentation

## Table of Contents
1. [Overview](#overview)
2. [Physical Theory](#physical-theory)
3. [Architecture](#architecture)
4. [Algorithms](#algorithms)
5. [Material Properties](#material-properties)
6. [References](#references)

---

## Overview

The MKF thermal analysis module calculates steady-state temperature distribution in magnetic components (inductors, transformers) by modeling:

- **Heat Generation**: Core losses (hysteresis, eddy current), winding losses (DC resistance, proximity, skin effect)
- **Heat Transfer**: Conduction (through solids), convection (to air), radiation (at high temperatures)
- **Geometric Modeling**: Node-based thermal network with surface quadrant subdivision

### Key Features

- Per-turn thermal modeling for accurate hotspot detection
- Quadrant-based surface representation for directional heat transfer
- Support for toroidal and concentric (E, U, ETD, etc.) core types
- Natural and forced convection modeling
- Insulation layer modeling (inter-turn, turn-to-core)
- Automatic thermal schematic generation

---

## Physical Theory

### Thermal-Electrical Analogy

The thermal system is modeled using the well-established thermal-electrical analogy:

| Thermal Quantity | Symbol | Unit | Electrical Equivalent |
|-----------------|--------|------|----------------------|
| Temperature | T | °C | Voltage (V) |
| Heat Flow | Q | W | Current (I) |
| Thermal Resistance | R_th | K/W | Resistance (R) |
| Thermal Capacitance | C_th | J/K | Capacitance (C) |

**Ohm's Law for Heat Transfer:**
```
ΔT = Q × R_th
```

### Conduction Resistance

For heat conduction through a solid:

```
R_cond = L / (k × A)  [K/W]
```

Where:
- **L** = Length of heat path [m]
- **k** = Thermal conductivity [W/(m·K)]
- **A** = Cross-sectional area perpendicular to heat flow [m²]

**For composite materials** (e.g., wire with insulation):
```
R_total = R_wire + R_insulation1 + R_gap + R_insulation2
        = t_wire/k_wire/A + t_ins/k_ins/A + gap/k_air/A + t_ins/k_ins/A
```

### Convection Resistance

For heat transfer to surrounding fluid (air):

```
R_conv = 1 / (h × A)  [K/W]
```

Where:
- **h** = Convection coefficient [W/(m²·K)]
- **A** = Exposed surface area [m²]

#### Natural Convection (Churchill-Chu Correlation)

For vertical plates and cylinders:

```
Nu_L = {0.825 + 0.387 Ra_L^(1/6) / [1 + (0.492/Pr)^(9/16)]^(8/9)}²

h = Nu_L × k_air / L
```

Where:
- **Nu_L** = Nusselt number (dimensionless)
- **Ra_L** = Rayleigh number = g β (T_s - T_∞) L³ / (ν α)
- **Pr** = Prandtl number ≈ 0.71 for air
- **k_air** = Thermal conductivity of air [W/(m·K)]
- **L** = Characteristic length (surface area / perimeter) [m]

**Grashof Number:**
```
Gr = g β (T_s - T_∞) L³ / ν²
```

**Rayleigh Number:**
```
Ra = Gr × Pr
```

### Radiation Resistance

At high temperatures, radiation becomes significant:

```
R_rad = 1 / (h_rad × A)  [K/W]

h_rad = ε × σ × (T_s² + T_a²) × (T_s + T_a)
```

Where:
- **ε** = Surface emissivity (0-1, typically 0.9 for oxidized copper)
- **σ** = Stefan-Boltzmann constant = 5.67×10⁻⁸ W/(m²·K⁴)
- **T_s** = Surface temperature [K]
- **T_a** = Ambient temperature [K]

---

## Architecture

### Node-Based Thermal Network

The magnetic component is discretized into thermal nodes:

```
┌─────────────────────────────────────────────────────────────┐
│                    Thermal Network                          │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│   Core Nodes          Turn Nodes              Ambient      │
│   ┌─────┐            ┌───────┐                 ┌───┐       │
│   │Core0│◄──────────►│Turn0_i│◄───────────────►│   │       │
│   └──┬──┘            └───┬───┘                 │ A │       │
│      │                   │                     │ m │       │
│   ┌──▼──┐            ┌───▼───┐                 │ b │       │
│   │Core1│◄──────────►│Turn0_o│◄───────────────►│ i │       │
│   └─────┘            └───┬───┘                 │ e │       │
│                          │                     │ n │       │
│   Bobbin Nodes      ┌────▼────┐                │ t │       │
│   ┌───────┐         │Turn1_i  │◄──────────────►│   │       │
│   │Bobbin0│◄───────►└─────────┘                └───┘       │
│   └───────┘                                                 │
│                                                             │
└─────────────────────────────────────────────────────────────┘
```

### ThermalNode Structure

Each node contains:

```cpp
struct ThermalNetworkNode {
    // Identification
    std::string name;              // e.g., "Turn_0_Inner", "Core_Column_0"
    ThermalNodePartType part;      // CORE, TURN, BOBBIN, AMBIENT
    
    // Geometry
    std::vector<double> physicalCoordinates;  // [x, y, z] in meters
    Dimensions dimensions;         // width, height, depth
    
    // Thermal State
    double temperature;            // Current temperature [°C]
    double powerDissipation;       // Heat generation [W]
    
    // Surface Quadrants (4 quadrants for turns)
    std::array<ThermalNodeQuadrant, 4> quadrants;
};
```

### Quadrant System

Turn surfaces are divided into 4 quadrants for directional heat transfer modeling:

**Toroidal Cores (Polar Coordinates):**

```
                    TL (Tangential Left, +θ direction)
                              ↑
                              │
    RI (Radial Inner) ←───────●───────→ RO (Radial Outer)
    Toward center             │        Away from center
    (smaller radius)          │        (larger radius)
                              ↓
                    TR (Tangential Right, -θ direction)
```

**Concentric Cores (Cartesian Coordinates):**

```
                         T (Top, +Y)
                          ↑
                          │
    L (Left, -X) ←────────●────────→ R (Right, +X)
                          │
                          ↓
                         B (Bottom, -Y)
```

### ThermalResistanceElement Structure

Connections between node quadrants:

```cpp
struct ThermalResistanceElement {
    size_t nodeFromId;          // Source node index
    ThermalNodeFace quadrantFrom;  // Source quadrant face
    size_t nodeToId;            // Destination node index
    ThermalNodeFace quadrantTo;    // Destination quadrant face
    
    HeatTransferType type;      // CONDUCTION, CONVECTION, RADIATION
    double resistance;          // Thermal resistance [K/W]
    
    // For conduction calculations
    double length;              // Heat path length [m]
    double area;                // Contact area [m²]
    double thermalConductivity; // Material k [W/(m·K)]
};
```

---

## Algorithms

### Turn-to-Turn Connection Algorithm

This algorithm determines which turn quadrants should be connected via conduction.

```
ALGORITHM: CreateTurnToTurnConnections

INPUT:  turnNodes - list of turn node indices
        minConductionDist - threshold distance for connection
OUTPUT: thermalResistances - list of resistance elements

FOR each pair of turns (i, j) where i < j:
    
    // 1. Calculate center-to-center distance
    dx = turn[i].x - turn[j].x
    dy = turn[i].y - turn[j].y
    centerDist = sqrt(dx² + dy²)
    
    // 2. Calculate surface distance
    extent_i = min(turn[i].width, turn[i].height) / 2
    extent_j = min(turn[j].width, turn[j].height) / 2
    surfaceDist = centerDist - extent_i - extent_j
    
    // 3. Check connection criterion
    threshold = min(minDim_i, minDim_j) / 4
    IF surfaceDist > threshold:
        CONTINUE  // Too far apart
    
    // 4. Find facing quadrants (toroidal only)
    IF toroidalCore:
        angle_i = atan2(turn[i].y, turn[i].x)
        angle_j = atan2(turn[j].y, turn[j].x)
        
        // Direction from turn i to turn j
        dirX = -dx / centerDist
        dirY = -dy / centerDist
        
        // Find quadrant of turn i that best faces turn j
        bestDot = -1
        FOR each quadrant Q in turn[i]:
            qx, qy = GetQuadrantDirection(angle_i, Q.face)
            dot = qx*dirX + qy*dirY
            IF dot > bestDot:
                bestDot = dot
                face_i = Q.face
        
        // Find quadrant of turn j that best faces turn i
        bestDot = -1
        FOR each quadrant Q in turn[j]:
            qx, qy = GetQuadrantDirection(angle_j, Q.face)
            dot = qx*(-dirX) + qy*(-dirY)
            IF dot > bestDot:
                bestDot = dot
                face_j = Q.face
        
        // Create resistance between facing quadrants
        CREATE_RESISTANCE(turn[i], face_i, turn[j], face_j)
    
    ELSE:
        // Concentric: use standard quadrant based on relative position
        CREATE_RESISTANCE(turn[i], face_i, turn[j], face_j)
```

### Convection Blocking Algorithm

Determines which quadrants are exposed to air (convection) vs blocked by other turns.

```
ALGORITHM: DetermineConvectionExposure

INPUT:  node - turn node to check
        allNodes - list of all nodes
OUTPUT: exposedQuadrants - list of quadrants with convection to air

maxBlockingDist = max(wireWidth, wireHeight)
minRadialDiff = min(wireWidth, wireHeight) / 4

FOR each quadrant Q in node.quadrants:
    isExposed = TRUE
    
    IF Q.face == RADIAL_INNER:
        FOR each otherNode in allNodes:
            IF otherNode == node: CONTINUE
            
            otherR = sqrt(otherNode.x² + otherNode.y²)
            nodeR = sqrt(node.x² + node.y²)
            radialDiff = nodeR - otherR
            
            // Check if significantly closer AND within blocking range
            IF radialDiff > minRadialDiff AND radialDiff < maxBlockingDist:
                angleDiff = |atan2(node.y, node.x) - atan2(otherNode.y, otherNode.x)|
                IF angleDiff < 0.3 radians:  // ~17 degrees
                    isExposed = FALSE
                    BREAK
    
    ELSE IF Q.face == RADIAL_OUTER:
        // Similar logic for outward direction
        ...
    
    IF isExposed:
        CREATE_CONVECTION_RESISTANCE(node, Q.face, ambient)
```

### Thermal Circuit Solver

The thermal network is solved as a linear system:

```
ALGORITHM: SolveThermalCircuit

INPUT:  nodes - list of thermal nodes
        resistances - list of thermal resistances
OUTPUT: nodeTemperatures - temperature at each node

// Build conductance matrix G (inverse of resistance)
// Matrix size: N × N where N = number of nodes

FOR each node i:
    G[i][i] = sum of conductances connected to node i
    FOR each resistance R connected between i and j:
        conductance = 1 / R.resistance
        G[i][j] -= conductance
        G[i][i] += conductance

// Build power vector P
FOR each node i:
    P[i] = nodes[i].powerDissipation

// Solve system: G × T = P
// Using Gaussian elimination or sparse matrix solver
T = SOLVE_LINEAR_SYSTEM(G, P)

RETURN T
```

---

## Material Properties

### Ferrite Core Materials

| Material | k [W/(m·K)] | Source |
|----------|-------------|--------|
| 3C90 | 4.0-4.5 | Ferroxcube |
| 3C95 | 4.0-4.5 | Ferroxcube |
| 3C97 | 4.5-5.0 | Ferroxcube |
| 3F3 | 4.0 | Ferroxcube |
| 3F4 | 4.0 | Ferroxcube |
| N87 | 4.0-4.5 | TDK |
| N97 | 4.0-4.5 | TDK |
| N49 | 4.0 | TDK |
| PC40 | 4.0 | TDK |
| PC95 | 4.0 | TDK |

**Note:** Ferrites have low thermal conductivity compared to metals, leading to significant temperature gradients in the core.

### Copper Wire

| Property | Value | Notes |
|----------|-------|-------|
| Thermal Conductivity | 385-400 W/(m·K) | @ 20°C, decreases with temperature |
| Electrical Resistivity | 1.68×10⁻⁸ Ω·m | @ 20°C (annealed) |
| Temp. Coefficient (α) | 0.00393 /°C | For resistance calculation |
| Specific Heat | 385 J/(kg·K) | |
| Density | 8960 kg/m³ | |

### Wire Insulation (Enamel)

| Type | Thickness [μm] | k [W/(m·K)] | Max Temp [°C] |
|------|---------------|-------------|---------------|
| Polyurethane (UEW) | 20-30 | 0.25 | 130-155 |
| Polyester (PEW) | 25-35 | 0.24 | 155 |
| Polyimide (AIW) | 30-50 | 0.35 | 200-220 |

Insulation adds significant thermal resistance in the transverse direction.

### Air Properties

| Property | @ 25°C | @ 50°C | @ 100°C |
|----------|--------|--------|---------|
| Density [kg/m³] | 1.184 | 1.093 | 0.946 |
| k [W/(m·K)] | 0.0261 | 0.0278 | 0.0314 |
| Prandtl Number | 0.71 | 0.71 | 0.71 |
| Kinematic Viscosity [m²/s] | 15.7×10⁻⁶ | 17.9×10⁻⁶ | 23.1×10⁻⁶ |

### Bobbin Materials

| Material | k [W/(m·K)] | Max Temp [°C] |
|----------|-------------|---------------|
| Phenolic | 0.2-0.3 | 150 |
| PBT (Polybutylene Terephthalate) | 0.2-0.25 | 140 |
| PA66 (Nylon) | 0.25-0.3 | 130 |
| PET | 0.24 | 120 |

---

## References

### Heat Transfer Fundamentals

1. **Incropera, F. P., & DeWitt, D. P. (2002)**. Fundamentals of Heat and Mass Transfer (6th ed.). Wiley.
   - Chapter 1: Introduction
   - Chapter 3: One-Dimensional, Steady-State Conduction
   - Chapter 7: External Convection
   - Chapter 9: Natural Convection

2. **Churchill, S. W., & Chu, H. H. S. (1975)**. Correlating equations for laminar and turbulent free convection from a vertical plate. *International Journal of Heat and Mass Transfer*, 18(11), 1323-1329.
   - Source of the widely-used Churchill-Chu correlation for natural convection

### Magnetic Component Thermal Modeling

3. **Kyaw, P. T., et al. (2018)**. Thermal modeling of litz wire windings using homogenization techniques. *IEEE Transactions on Magnetics*, 54(11), 1-5.
   - Anisotropic thermal conductivity in windings
   - k_longitudinal ≈ 150-200 W/(m·K), k_transverse ≈ 0.5-2 W/(m·K)

4. **López, I., et al. (2019)**. Anisotropic thermal conductivity in magnetic component windings. *IEEE Transactions on Power Electronics*, 34(4), 3754-3763.
   - Comprehensive study of effective thermal conductivity
   - Gap effects and contact resistance

5. **Mukosiej, D., et al. (2022)**. Thermal modeling of high-frequency magnetic components using homogenization techniques. *IEEE Transactions on Power Electronics*, 37(2), 2156-2169.
   - Multi-scale modeling approaches
   - Litz wire thermal modeling

6. **Sullivan, C. R. (2022)**. Power Magnetics for Power Electronics: A Guide to Design and Applications. Springer.
   - Chapter 6: Thermal Design
   - Practical thermal management guidelines

### Industry Standards

7. **IEC 60317** Specifications for particular types of winding wires
   - Insulation thickness specifications
   - Thermal class ratings

8. **JEDEC JESD51** Methodology for the Thermal Measurement of Component Packages
   - Thermal measurement standards
   - Test board specifications

9. **MIL-STD-883** Test Methods and Procedures for Microelectronics
   - Method 1012: Thermal Characteristics

### Online Resources

10. **Ferroxcube Datasheets** - https://www.ferroxcube.com/
    - Thermal conductivity data for ferrite materials

11. **TDK Datasheets** - https://www.tdk.com/
    - Ferrite material properties

12. **NIST Chemistry WebBook** - https://webbook.nist.gov/
    - Air thermophysical properties

---

## Implementation Notes

### Numerical Stability

- Temperatures are solved in Celsius, but radiation calculations use Kelvin
- Convergence tolerance typically 0.1°C
- Maximum iterations: 100 (usually converges in 10-20)

### Solver Performance

- System is typically sparse (each node connects to few neighbors)
- For N turns with 2 nodes each: ~2N equations
- Solve time: O(N^1.5) for sparse direct solver

### Accuracy Considerations

1. **Conduction through contacts**: Contact resistance not modeled; assumes perfect contact
2. **Convection coefficients**: Empirical correlations have ±20% uncertainty
3. **Material properties**: Temperature-dependent k is approximated as constant
4. **Geometry**: Round wires approximated as rectangular for quadrant calculation

---

## File Organization

```
src/physical_models/
├── Temperature.h          # Main API and configuration
├── Temperature.cpp        # Implementation
├── ThermalNode.h          # Node and resistance structures
├── ThermalNode.cpp        # Node methods and utilities
└── ThermalResistance.h    # Heat transfer calculations

docs/thermal/
├── ARCHITECTURE.md              # High-level architecture
├── THERMAL_MODEL_DOCUMENTATION.md # This file
├── THERMAL_MODEL_ANALYSIS.md     # Literature review
└── thermal_resistance_models.md  # Resistance formulas
```

---

*Document Version: 1.0*
*Last Updated: February 2025*
*Maintainer: OpenMagnetics Team*
