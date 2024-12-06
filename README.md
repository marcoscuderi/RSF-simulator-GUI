Documentation for Rate-and-State Friction Simulator
Overview
The Rate-and-State Friction Simulator is a Tkinter-based graphical user interface (GUI) designed to help users visualize the behavior of rate-and-state friction laws. The application allows users to:
1.	Input key parameters (a, b, Dc, k, velocity steps, and times).
2.	Choose a friction evolution law.
3.	Simulate and visualize:
·	Friction vs. Displacement
·	Velocity vs. Displacement
·	Theta (State Variable) vs. Displacement
1.	Interact with synchronized plots (shared x-axis) and zoom/pan functionality using a Matplotlib toolbar.
Application Features
Input Parameters
1.	Parameter Fields:
·	a (dimensionless): A small parameter representing the direct effect in the friction law.
·	b (dimensionless): A small parameter representing the evolution effect in the friction law.
·	Dc (µm): Critical slip distance for state variable evolution.
·	k (N/µm): Stiffness of the fault system.
1.	Velocity Steps and Times:
·	Velocity Steps: Comma-separated values defining different velocity steps during the simulation.
·	Velocity Times: Comma-separated values specifying the duration for each velocity step.
1.	Friction Evolution Law:
·	Choose from:
·	Dieterich
·	Ruina
·	Dieterich with RD
·	Ruina with RD
1.	Critical Stiffness Calculation:
·	Automatically calculates and displays (b - a) / Dc dynamically as inputs are updated.
Outputs and Visualization
1.	Interactive Plots:
·	Friction vs. Displacement: Shows how the friction coefficient evolves with slip displacement.
·	Velocity vs. Displacement: Visualizes velocity changes during slip.
·	Theta vs. Displacement: Tracks the evolution of the state variable over displacement.
1.	Shared X-Axis:
·	Zooming or panning one plot synchronizes the view for all three plots.
1.	Matplotlib Toolbar:
·	Provides tools for zooming, panning, and resetting the plots.
How to Use the Application
Step-by-Step Guide
1.	Launch the GUI:
·	Run the Python script containing the GUI.
1.	Input Parameters:
·	Enter values for a, b, Dc, and k in the respective fields.
·	Specify velocity steps and times as comma-separated values in the respective fields.
1.	Choose Friction Law:
·	Select one of the provided friction evolution laws using the radio buttons.
1.	Run the Simulation:
·	Click the Run Simulation button.
·	The application will compute and update the plots.
1.	Interact with Plots:
·	Use the Matplotlib toolbar (above the first plot) to:
·	Zoom: Click the magnifying glass icon and drag over the area of interest.
·	Pan: Click the hand icon and drag the plot to adjust the view.
·	Reset: Click the home icon to return to the default view.
Application Layout
Inputs Section (Left Pane)
·	Input fields for a, b, Dc, k.
·	Text boxes for velocity steps and times.
·	Dropdown for selecting the friction law.
·	Critical stiffness display.
Visualization Section (Right Pane)
·	Three vertically stacked plots:
1.	Friction vs. Displacement
2.	Velocity vs. Displacement
3.	Theta vs. Displacement
·	Shared Matplotlib toolbar at the top.
Error Handling
Input Validation
·	Ensure valid numerical input in the fields for a, b, Dc, k, velocity steps, and times.
·	Errors in input will raise an exception in the console output.
Debugging
·	If the simulation does not run as expected, check the console output for errors.
Dependencies
Required Libraries
·	Python 3.x
·	Matplotlib: For plotting.
·	Tkinter: For the GUI.
·	NumPy: For numerical calculations.
Installation Commands
bash
Copy code
pip install matplotlib numpy

Code Structure
Key Functions
1.	run_simulation:
·	Core function performing the simulation of rate-and-state friction laws.
·	Outputs displacement, friction coefficient, velocity, and theta arrays.
1.	update_plots:
·	Reads input parameters.
·	Runs the simulation using run_simulation.
·	Updates the plots dynamically with the simulation results.
1.	Dynamic Widgets:
·	Input fields dynamically trigger updates (e.g., recalculating critical stiffness).
Example Use Case
Simulating a Velocity-Step Experiment
1.	Set parameters:
·	a = 0.01
·	b = 0.012
·	Dc = 1 µm
·	k = 0.01 N/µm
1.	Define velocity steps:
·	Steps: 1, 10, 1
·	Times: 300, 300, 300
1.	Choose "Dieterich" as the friction law.
2.	Click Run Simulation.
3.	Use the toolbar to zoom in on areas of interest in the plots.


