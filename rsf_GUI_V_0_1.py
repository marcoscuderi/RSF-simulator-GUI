import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import matplotlib.pyplot as plt
import numpy as np
from tkinter import ttk, filedialog

# Mock function to simulate the model output
def run_simulation(a, b, Dc, k, law, v_steps, v_times):
    import sys
    import numpy as np
    import matplotlib.pylab as plt
    import math
    import datetime

    # adaptive step size control parameters
    EPSILON = 1e-13
    PGROW = -0.200
    PSHRINK = -0.250
    FCOR = 1.0 / 15.0
    SAFETY = 0.9
    ERRCON = 6.0e-4
    VPEAK_WARN_THRES = 0.01
    MUPEAK_WARN_THRES = 1e-4

    # These are used globally so define them here
    myConstants = [0, 0, 0.5, 0.5, 1]  # these are used in do_rk
    OneOverSix = 1 / 6  # multiplication is better than division millions of times in do_rk
    rk_J = np.zeros(5)  # arrays for Runge-Kutta. These are where we store slope estimates
    rk_K = np.zeros(5)  # J is theta_1, K is theta_2, M is friction
    rk_M = np.zeros(5)  # the first term in each array is always zero; it never changes.
    rk_V = np.zeros(5)  # velocity term for radiation damping

    #with radiation damping
    def do_rk_d_rd(mu, theta, t_Step, vel):		
        for ii in range(1, 5):
            temp_theta = theta + rk_J[ii - 1] * myConstants[ii]
            v_temp = vel	 + rk_V[ii - 1] * myConstants[ii]
            rk_J[ii] = (1 - (v_temp * temp_theta) * OneOverDc) * t_Step
            new_dtheta1_dt = rk_J[ii]/t_Step
            rk_V[ii] = t_Step * ( ( stiff*(v_lp - v_temp) - (fric_b*new_dtheta1_dt/temp_theta) ) / (fric_a/v_temp + Gprime) )
        # end of for loop

        theta += (rk_J[1] + 2 * rk_J[2] + 2 * rk_J[3] + rk_J[4]) * OneOverSix
        vel += (rk_V[1] + 2 * rk_V[2] + 2 * rk_V[3] + rk_V[4]) * OneOverSix
        mu = mu_o + fric_a*np.log(vel/v_o) + fric_b*np.log(v_oOnDc * theta) 

        return mu, theta, vel

    def do_rk_d(mu, theta, t_Step, vel):
        for ii in range(1, 5):
            temp_mu = mu + rk_M[ii - 1] * myConstants[ii]
            temp_theta = theta + rk_J[ii - 1] * myConstants[ii]

            v_temp = v_o * np.exp(
                (temp_mu - mu_o - (fric_b * np.log(v_oOnDc * temp_theta))) * OneOver_fric_a
            )  # do a partial update on v

            rk_M[ii] = stiff * (v_lp - v_temp) * t_Step
            rk_J[ii] = (1 - (v_temp * temp_theta) * OneOverDc) * t_Step

        # end of for loop

        theta += (rk_J[1] + 2 * rk_J[2] + 2 * rk_J[3] + rk_J[4]) * OneOverSix
        # theta += (rk_K[1] + 2*rk_K[2] + 2*rk_K[3] + rk_K[4])*OneOverSix   #2nd state variable
        mu += (rk_M[1] + 2 * rk_M[2] + 2 * rk_M[3] + rk_M[4]) * OneOverSix

        vel = v_o * np.exp( (mu - mu_o - (fric_b * np.log(v_oOnDc * theta))) / fric_a)

        return mu, theta, vel


    def do_rk_r(mu, theta, t_Step, vel):
        for ii in range(1, 5):
            temp_mu = mu + rk_M[ii - 1] * myConstants[ii]
            temp_theta = theta + rk_J[ii - 1] * myConstants[ii]

            v_temp = v_o * np.exp(
                (temp_mu - mu_o - (fric_b * np.log(v_oOnDc * temp_theta))) * OneOver_fric_a
            )  # do a partial update on v

            rk_M[ii] = stiff * (v_lp - v_temp) * t_Step
            rk_J[ii] = (
                -(v_temp * temp_theta * OneOverDc)
                * np.log(v_temp * temp_theta * OneOverDc)
                * t_Step
            )

        # end of for loop

        theta += (rk_J[1] + 2 * rk_J[2] + 2 * rk_J[3] + rk_J[4]) * OneOverSix
        # theta += (rk_K[1] + 2*rk_K[2] + 2*rk_K[3] + rk_K[4])*OneOverSix   #2nd state variable
        mu += (rk_M[1] + 2 * rk_M[2] + 2 * rk_M[3] + rk_M[4]) * OneOverSix

        vel = v_o * np.exp( (mu - mu_o - (fric_b * np.log(v_oOnDc * theta))) / fric_a)

        return mu, theta, vel


    # set RSF parameters
    fric_a = a
    fric_b = b
    Dc = Dc  # micrion
    stiff = k  # friction/micron
    #law = "d"  # Dieterich law
    # law = 'r'	 #Ruina law
    #law = "d_rd"  # Dieterich law w Radiation damping
    # law = 'r_rd'	 #Ruina law w Radiation damping
    law = law

    Gprime = 1e-6  #G/(sigma_n * beta) Radiation damping term in s/micron; this is for G=30 MPa, sigma=10MPa, beta=3km/s
    mu_o = 0.6  # ref friction
    v_o = 1.0  # ref velocity
    theta_o = Dc / v_o  # ref state

    OneOverDc = 1 / Dc
    OneOver_fric_a = 1 / fric_a
    v_oOnDc = v_o / Dc  # no need to do this calculation millions of times in do_rk

    vel = v_o
    mu_err_scale = mu_o  # set error scale for rk adaptive step size 
    theta_err_scale = Dc / vel
    vel_err_scale = vel*1e3

    # main calculation

    # indexes for vel sequence and arrays for data to plot 

    # for v_step this is the the velocity after the step; initial velocity is v_o
    # loadingVelocity = 1   #final velocity for v_step or reload velocity for SHS
    # hold_time = 1

    holdtimes = []			#for healing plots
    healingvalues = []			#for healing plots
    delta_mu_c = []			#for healing plots
    v_step_count = 0

    # hold_times = [1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000]

    v_steps = v_steps
    v_step_times = v_times
    #v_steps = [3, 10, 30, 100]
    #v_step_times = [100, 30, 10, 3]

    dataFor_thisV_step = np.zeros( len(v_steps))  # set up a counter so we know when to change velocity
    indexFor_V_step = np.zeros( len(v_steps))  # set up an index of v steps 

    if len(v_steps) != len(v_step_times):
        print("The number of v_steps and v_step times have to be equal\n")
        sys.exit()

    # for velocities
    v_step_count = 0
    arraysize = 2  # one point for initial state and one point at the first v_step
    step_size = 0.1  # this is for holds, assume that 0.1s is ok...

    for v_lp in v_steps:

        disp = v_lp * v_step_times[v_step_count]  # data for this v-step

        # keep 10 points per Dc.
        write_inc = Dc / 10  # how often should we save data to arrays?

        if v_lp == 0:
            dataFor_thisV_step[v_step_count] = int(v_step_times[v_step_count] / step_size)
        else:
            dataFor_thisV_step[v_step_count] = int( disp / write_inc)  # number of data points for this velocity

        arraysize += dataFor_thisV_step[v_step_count]
        indexFor_V_step[v_step_count] = int(arraysize+2) # index of where v steps occur, add 2 for first two points
        thisIndex = int(indexFor_V_step[v_step_count])
        print(
            "{} data points will be written for vel of {} of duration {}, index={}\n"
                .format( dataFor_thisV_step[v_step_count], v_lp, v_step_times[v_step_count],thisIndex)
        )
        v_step_count += 1
    # end of for loop: for v_lp in v_steps

    # these are arrays for output and plotting

    arraysize = int(arraysize)  # the += above makes this a float
    print(arraysize)
    mu_array = np.zeros(arraysize)
    x_lp_array = np.zeros(arraysize)
    theta_array = np.zeros(arraysize)
    vel_array = np.zeros(arraysize)
    slipDisp_array = np.zeros(arraysize)
    time_array = np.zeros(arraysize)

    # this is the approach for arrays defined by np.zeros

    x_lp_array[0] = slipDisp_array[0] = time_array[0] = -10
    x_lp_array[1] = slipDisp_array[1] = time_array[1] = 0
    mu_array[0] = mu_array[1] = 0.6
    theta_array[0] = theta_array[1] = Dc / v_o
    vel_array[0] = vel_array[1] = v_o
    slipDisp = slipDisp_array[1]

    # variables
    jj = 0
    dmu = 0
    dtheta = 0
    lpDisp = 0
    calcTime = 0
    # initial values so that we can track peaks and changes
    muPeak = 0
    mu = mu_array[1]
    theta = theta_array[1]
    vel = vel_array[1]

    i = 2  # I used values 0 and 1 above
    nextStep = timeStep = ( write_inc / v_o)  # initial time step. This will be changed using adaptive step size control
    v_step_count = 0  # set this back to zero for start

    current_time = datetime.datetime.now()
    print("Current time:", current_time)

    v_lp = v_steps[v_step_count]
    print("first velocity is = {:g} for duration {:g}\n".format(v_lp, v_step_times[0]))
    #print( "step size is {:g} and write_inc/v_lp = {:g}\n".format(step_size, write_inc / v_lp))
    #print("nextStep = {:g}\n".format(nextStep))
    thisIndex = int( indexFor_V_step[v_step_count])
    print("arraysize = {}, i={}, thisIndex={}, v_step_count={}\n".format(arraysize,i,thisIndex,v_step_count))
    print(indexFor_V_step[v_step_count])

    while i < arraysize:  # this is the main calculation. Numerical Int is done in do_rk
        # ttnd is time to next data point for save array, relative to current time. it's dt

        if i == thisIndex:  # done with this v_step
            v_step_count += 1
            v_lp = v_steps[v_step_count]
            thisIndex = int( indexFor_V_step[v_step_count])
            print( "done with vstep {:d} at i = {:d}, v_lp = {:g}\n".format( v_step_count, i, v_lp))

        if ( v_lp == 0):  # doing a hold; assume that data every 0.1 sec is ok. If not change it above
            ttnd = step_size  # calculate the time to the next data point needed
        else:
            ttnd = write_inc / v_lp  # calculate the time to the next data point needed

        tnow = 0  # time since last data written to arrays

            # if i < 10:
            # print(i, jj-1, calcTime, ttnd, nextStep, hold_time, hold_time_steps[jj-1])
            # if i > 990:
            # print(i, jj-1, calcTime, ttnd, nextStep, hold_time, hold_time_steps[jj-1])

        while tnow < ttnd:
            if (tnow + nextStep) >= ttnd:  # set time to next data (ttnd)
                timeStep = ttnd - tnow
            else:
                timeStep = nextStep

            AtimeStep_err = 1  # this is for adaptive time step control; it gets set to 0 when the time step is ok and we get a new step size (called nextStep)

            while AtimeStep_err == 1:
                    # Adaptive time Step size control
                    # save current values, for full step and errors
                old_mu = fs_mu = mu
                old_theta = fs_theta = theta
                old_vel = fs_vel = vel

                HalfStep = timeStep / 2
                fullStep = timeStep

                if law == "d":
                        # do half steps, first one and then a second one to update parameters
                    mu, theta, vel = do_rk_d(mu, theta, HalfStep, vel)  # 4th order Runge-Kutta
                    mu, theta, vel = do_rk_d(mu, theta, HalfStep, vel)  # 4th order Runge-Kutta
                        # full step
                    fs_mu, fs_theta, fs_vel = do_rk_d( fs_mu, fs_theta, fullStep, fs_vel)  # 4th order Runge-Kutta

                elif law == "r":
                    # do half steps, first one and then a second one to update parameters
                    mu, theta, vel = do_rk_r(mu, theta, HalfStep, vel)  # 4th order Runge-Kutta
                    mu, theta, vel = do_rk_r(mu, theta, HalfStep, vel)  # 4th order Runge-Kutta
                    # full step
                    fs_mu, fs_theta, fs_vel = do_rk_r( fs_mu, fs_theta, fullStep, fs_vel)  # 4th order Runge-Kutta

                elif law == "d_rd":
                        # do half steps, first one and then a second one to update parameters
                    mu, theta, vel = do_rk_d_rd(mu, theta, HalfStep, vel)  # 4th order Runge-Kutta
                    mu, theta, vel = do_rk_d_rd(mu, theta, HalfStep, vel)  # 4th order Runge-Kutta
                        # full step
                    fs_mu, fs_theta, fs_vel = do_rk_d_rd( fs_mu, fs_theta, fullStep, fs_vel)  # 4th order Runge-Kutta

                else:
                    print("hmmm, how could that happen? RSF law must be D or R, or rd version\n")
                    exit()

                # print("fullstep: mu = {}, theta = {}, fullSteptime = {}, vel= {}\n".format(fs_mu, fs_theta, fullStep, fs_vel))

                max_err = 0  # evaluate error

                vel_err = vel - fs_vel
                err = math.fabs(vel_err / vel_err_scale)
                max_err = err if err > max_err else max_err

                mu_err = mu - fs_mu
                err = math.fabs(mu_err / mu_err_scale)
                max_err = err if err > max_err else max_err

                theta_err = theta - fs_theta
                err = math.fabs(theta_err / theta_err_scale)
                max_err = err if err > max_err else max_err

                max_err /= EPSILON  # scale relative to tolerance

                if max_err <= 1.0:  # step succeeded, compute size of next step
                    timeStep = fullStep
                    mu = fs_mu
                    theta = fs_theta
                    vel = fs_vel
                    if max_err > ERRCON:
                        nextStep = SAFETY * timeStep * np.exp(PGROW * np.log(max_err))
                    else:
                        nextStep = 4 * timeStep
                    AtimeStep_err = 0
                    # print( "4here i={}, tnow={}, ttnd={}, timeStep={}, nextStep={}\n".format( i, tnow, ttnd, timeStep, nextStep))
                else:
                    timeStep = (
                        SAFETY * timeStep * np.exp(PSHRINK * np.log(max_err))
                    )  # truncation error too large, reduce step size
                    mu = old_mu
                    theta = old_theta
                    vel = old_vel
                # print("max_err={} \n".format(max_err))
                # print("here i={}, tnow={}, ttnd={}, timeStep={}, nextStep={}\n".format(i, tnow, ttnd, timeStep, nextStep))
                # sys.exit()

                # end of adaptive step size while loop

            # while tnow < ttnd:	now update everything with this time step (timeStep) and move forward until tnow = ttnd

            tnow += timeStep
            calcTime += timeStep

            mu += mu_err * FCOR  # fifth order RK, 5th order bit
            theta += theta_err * FCOR
            vel += vel_err * FCOR
            #vel = v_o * np.exp( (mu - mu_o - (fric_b * np.log(v_oOnDc * theta))) * OneOver_fric_a)

            lpDisp += v_lp * timeStep
            slipDisp += vel * timeStep  # update slider displacement
            # print("here i={}, tnow={}, ttnd={}, timeStep={}, nextStep={}\n".format(i, tnow, ttnd, timeStep, nextStep))
            # sys.exit()
            if mu > muPeak:
                muPeak = mu

        # while tnow < ttnd: is at this level
        # save values to array
        # print(i, mu[i-1], lpDisp, vel[i-1], v_lp, theta[i-1], current_mu, last_mu, current_vel, current_theta)
        x_lp_array[i] = lpDisp
        mu_array[i] = mu
        theta_array[i] = theta
        vel_array[i] = vel
        slipDisp_array[i] = slipDisp
        time_array[i] = calcTime
        i += 1

    return [x_lp_array,mu_array,theta_array,vel_array,slipDisp_array,time_array]



 #Update the plots based on simulation results

# Update critical stiffness dynamically
def update_critical_stiffness(*args):
    try:
        a = float(a_var.get())
        b = float(b_var.get())
        Dc = float(Dc_var.get())
        if Dc != 0:
            critical_stiffness.set(f"{(b - a) / Dc:.6f}")
        else:
            critical_stiffness.set("Undefined (Dc=0)")
    except ValueError:
        critical_stiffness.set("Invalid Input")

# Update the plots based on simulation results
def update_plots():
    global simulation_results
    try:
        # Fetch user inputs
        a = float(a_var.get())
        b = float(b_var.get())
        Dc = float(Dc_var.get())
        k = float(k_var.get())
        law = law_var.get()
        v_steps = [float(v.strip()) for v in v_steps_entry.get().split(",")]
        v_times = [float(t.strip()) for t in v_times_entry.get().split(",")]

        # Run simulation
        x_disp, mu_array, theta_array, vel_array, slip_dist, time = run_simulation(a, b, Dc, k, law, v_steps, v_times)

        # Cache results for later use
        simulation_results["x_disp"] = x_disp
        simulation_results["mu_array"] = mu_array
        simulation_results["theta_array"] = theta_array
        simulation_results["vel_array"] = vel_array
        simulation_results["slip_dist"] = slip_dist
        simulation_results["time"] = time
        simulation_results["parameters"] = {
            "a": a,
            "b": b,
            "Dc": Dc,
            "k": k,
            "law": law,
            "v_steps": v_steps,
            "v_times": v_times,
        }

        # Update Friction vs Displacement
        ax1.clear()
        ax1.plot(x_disp, mu_array, label="Friction")
        ax1.set_title("Friction vs Displacement")
        ax1.set_ylabel("Friction Coefficient")
        ax1.legend()

        # Update Velocity vs Displacement
        ax2.clear()
        ax2.plot(x_disp, vel_array, label="Velocity")
        ax2.set_title("Velocity vs Displacement")
        ax2.set_ylabel("Velocity (µm/s)")
        ax2.legend()

        # Update Theta vs Displacement
        ax3.clear()
        ax3.plot(x_disp, theta_array, label="Theta")
        ax3.set_title("Theta vs Displacement")
        ax3.set_xlabel("Displacement (µm)")
        ax3.set_ylabel("Theta (s)")
        ax3.legend()

        # Redraw the canvases
        canvas.draw()
        print("Simulation and plots updated successfully.")
    except Exception as e:
        print(f"Error in updating plots: {str(e)}")

def export_data():
    try:
        # Open save file dialog
        file_path = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text files", "*.txt")])
        if not file_path:
            return

        # Retrieve cached simulation results
        x_disp = simulation_results["x_disp"]
        mu_array = simulation_results["mu_array"]
        theta_array = simulation_results["theta_array"]
        vel_array = simulation_results["vel_array"]
        slip_dist = simulation_results["slip_dist"]
        time = simulation_results["time"]
        parameters = simulation_results["parameters"]

        # Check if simulation results are available
        if len(x_disp) == 0 or len(mu_array) == 0 or len(time) == 0:
            print("No simulation results available to export.")
            return

        # Write the data to the file
        with open(file_path, "w") as file:
            # Write simulation parameters
            file.write("Simulation Parameters:\n")
            file.write(f"a: {parameters['a']}\n")
            file.write(f"b: {parameters['b']}\n")
            file.write(f"Dc: {parameters['Dc']}\n")
            file.write(f"k: {parameters['k']}\n")
            file.write(f"Friction Law: {parameters['law']}\n")
            file.write(f"Velocity Steps: {', '.join(map(str, parameters['v_steps']))}\n")
            file.write(f"Velocity Times: {', '.join(map(str, parameters['v_times']))}\n\n")

            # Write simulation data
            file.write("Simulation Data:\n")
            file.write("Time (s)\tDisplacement (µm)\tFriction Coefficient\tVelocity (µm/s)\tTheta (s)\tSlip Distance (µm)\n")
            for i in range(len(time)):
                file.write(f"{time[i]:.6f}\t{x_disp[i]:.6f}\t{mu_array[i]:.6f}\t{vel_array[i]:.6f}\t{theta_array[i]:.6f}\t{slip_dist[i]:.6f}\n")

        print(f"Data exported successfully to {file_path}")
    except Exception as e:
        print(f"Error exporting data: {str(e)}")

# Global variables to store simulation results
simulation_results = {
    "x_disp": [],
    "mu_array": [],
    "theta_array": [],
    "vel_array": [],
    "slip_dist": [],
    "time": [],
    "parameters": {},
}

# Main window
root = tk.Tk()
root.title("Rate and State Friction Simulator")

# Frame for user inputs
frame = ttk.Frame(root, padding="10")
frame.grid(row=0, column=0, sticky="nsew")

# Variables for user input
a_var = tk.StringVar(value="0.01")
b_var = tk.StringVar(value="0.012")
Dc_var = tk.StringVar(value="1")
k_var = tk.StringVar(value="0.01")
critical_stiffness = tk.StringVar(value="")
law_var = tk.StringVar(value="d")

# Bind critical stiffness update to changes in a, b, Dc
a_var.trace_add("write", update_critical_stiffness)
b_var.trace_add("write", update_critical_stiffness)
Dc_var.trace_add("write", update_critical_stiffness)


# Input fields for parameters
ttk.Label(frame, text="a (dimensionless):").grid(row=0, column=0, sticky="w")
a_entry = ttk.Entry(frame, textvariable=a_var)
a_entry.grid(row=0, column=1, sticky="ew")

ttk.Label(frame, text="b (dimensionless):").grid(row=1, column=0, sticky="w")
b_entry = ttk.Entry(frame, textvariable=b_var)
b_entry.grid(row=1, column=1, sticky="ew")

ttk.Label(frame, text="Dc (µm):").grid(row=2, column=0, sticky="w")
Dc_entry = ttk.Entry(frame, textvariable=Dc_var)
Dc_entry.grid(row=2, column=1, sticky="ew")

ttk.Label(frame, text="k (N/µm):").grid(row=3, column=0, sticky="w")
k_entry = ttk.Entry(frame, textvariable=k_var)
k_entry.grid(row=3, column=1, sticky="ew")

ttk.Label(frame, text="Critical Stiffness ((b-a)/Dc):").grid(row=4, column=0, sticky="w")
critical_value = ttk.Entry(frame, textvariable=critical_stiffness, state="readonly", width=15)
critical_value.grid(row=4, column=1, sticky="ew")

laws = [("Dieterich", "d"), ("Ruina", "r"), ("Dieterich with RD", "d_rd"), ("Ruina with RD", "r_rd")]
for idx, (label, value) in enumerate(laws):
    tk.Radiobutton(frame, text=label, variable=law_var, value=value).grid(row=5 + idx, column=0, sticky="w")

# Velocity steps and times
ttk.Label(frame, text="Velocity Steps (comma-separated):").grid(row=10, column=0, sticky="w")
v_steps_entry = ttk.Entry(frame)
v_steps_entry.insert(0, "3")
v_steps_entry.grid(row=10, column=1, sticky="ew")

ttk.Label(frame, text="Velocity Times (comma-separated):").grid(row=11, column=0, sticky="w")
v_times_entry = ttk.Entry(frame)
v_times_entry.insert(0, "300")
v_times_entry.grid(row=11, column=1, sticky="ew")

# Create matplotlib figures with shared x-axis
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize=(6, 8))
fig.tight_layout(pad=3)
canvas = FigureCanvasTkAgg(fig, master=root)
canvas_widget = canvas.get_tk_widget()
canvas_widget.grid(row=0, column=1, rowspan=3, sticky="nsew")

# Single toolbar for all plots
toolbar_frame = ttk.Frame(root)
toolbar_frame.grid(row=0, column=1, sticky="n")
toolbar = NavigationToolbar2Tk(canvas, toolbar_frame)
toolbar.update()

# Button to run simulation
run_button = ttk.Button(frame, text="Run Simulation", command=update_plots)
run_button.grid(row=12, column=0, sticky="ew")

export_button = ttk.Button(frame, text="Export Data", command=export_data)
export_button.grid(row=13, column=0, sticky="ew")

# Configure window resizing
root.rowconfigure(0, weight=1)
root.rowconfigure(1, weight=1)
root.rowconfigure(2, weight=1)
root.columnconfigure(1, weight=3)

frame.rowconfigure(list(range(13)), weight=1)
frame.columnconfigure(1, weight=1)

root.mainloop()