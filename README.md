# Rate-and-State Friction Simulator GUI (Tkinter Version)

## Overview

The **Rate-and-State Friction Simulator GUI** is a Python application built with Tkinter that allows users to:
1. Configure key parameters for rate-and-state friction simulations.
2. Visualize results interactively.
3. Export simulation parameters and time-series data for further analysis.

---

## Features

### Parameter Input
- **a** (dimensionless): Direct effect parameter.
- **b** (dimensionless): Evolution effect parameter.
- **Dc** (µm): Critical slip distance.
- **k** (N/µm): Stiffness of the fault system.
- **Critical Stiffness**: Dynamically calculated as `(b - a) / Dc`.
- Friction laws:
  - Dieterich
  - Ruina
  - Dieterich with RD
  - Ruina with RD
- Velocity steps and times: Configurable via comma-separated input fields.

### Visualization
- Interactive plots for:
  - Friction vs. Displacement
  - Velocity vs. Displacement
  - Theta (State Variable) vs. Displacement
- Shared x-axis for consistent zooming and panning.
- Matplotlib toolbar for zoom, pan, and reset.

### Data Export
- Export parameters and time-series data to a `.txt` file, including:
  - Displacement
  - Friction
  - Velocity
  - Theta
  - Slip distance
  - Time

---

## How to Use

### Running a Simulation
1. Set parameter values (`a`, `b`, `Dc`, `k`) and choose a friction law.
2. Define velocity steps and times (comma-separated values).
3. Click **Run Simulation** to execute the simulation and update the plots.

### Exporting Results
1. Run a simulation.
2. Click **Export Data**.
3. Choose a file location and name.
4. The exported `.txt` file will include simulation parameters and all time-series data.

---

## Exported Data Format

### Example `.txt` File

```plaintext
Simulation Parameters:
a: 0.01
b: 0.012
Dc: 1.000000
k: 0.010000
Friction Law: Dieterich
Velocity Steps: 1.0, 10.0, 1.0
Velocity Times: 300.0, 300.0, 300.0

Simulation Data:
Time (s)        Displacement (µm)   Friction Coefficient   Velocity (µm/s)   Theta (s)   Slip Distance (µm)
0.000000        0.000000            0.600000               1.000000          1.000000    0.000000
0.100000        0.010000            0.601000               1.100000          0.990000    0.010000
...
