"""Task 5 of assignment Inviscid flow over an airfoil (TU Delft's AE4130)

This script produces the plots requested for task 5 of the assignment. An assessment of the effect of the panel density
is carried out.

The interval of number of panels studied can be changed modifying the appropriate parameter.

This script requires that `numpy` and `matplotlib` be installed within the Python environment you are running this
script in.
"""

# 3rd party packages
import numpy as np
import matplotlib.pyplot as plt
# Local source
import panelairfoil as panelairfoil


# Definition of NACA airfoil 4 digits
naca_4_digit_name = "0012"
# Array of number of panels for the panel density investigation
no_panels_array = np.arange(200, 0, -10)
# Initialize arrays for lift and moment coefficients
cl_alpha5 = np.zeros(len(no_panels_array))
cm_alpha5 = np.zeros(len(no_panels_array))
cl_alpha15 = np.zeros(len(no_panels_array))
cm_alpha15 = np.zeros(len(no_panels_array))
# Initialize arrays for relative error of lift and moment coefficients
error_cl_alpha5 = np.zeros(len(no_panels_array))
error_cm_alpha5 = np.zeros(len(no_panels_array))
error_cl_alpha15 = np.zeros(len(no_panels_array))
error_cm_alpha15 = np.zeros(len(no_panels_array))
# Iterate through the different numbers of panels
for i in range(len(no_panels_array)):
    # Generate the Naca4DigitPanelled object for the current number of panels
    naca0012 = panelairfoil.Naca4DigitPanelled(naca_4_digit_name, no_panels_array[i])
    # Calculate lift and moment coefficients for 5 and 15 degree angle of attack
    cl_alpha5[i] = naca0012.calculate_cl_vs_alpha(np.array([np.deg2rad(5)]))[:, 1]
    cm_alpha5[i] = naca0012.calculate_cm_vs_alpha(np.array([np.deg2rad(5)]))[:, 1]
    cl_alpha15[i] = naca0012.calculate_cl_vs_alpha(np.array([np.deg2rad(15)]))[:, 1]
    cm_alpha15[i] = naca0012.calculate_cm_vs_alpha(np.array([np.deg2rad(15)]))[:, 1]
    # Compute the relative error w.r.t. the case with 200 panels (first case of the iteration)
    error_cl_alpha5[i] = cl_alpha5[i] / cl_alpha5[0] - 1
    error_cm_alpha5[i] = cm_alpha5[i] / cm_alpha5[0] - 1
    error_cl_alpha15[i] = cl_alpha15[i] / cl_alpha15[0] - 1
    error_cm_alpha15[i] = cm_alpha15[i] / cm_alpha15[0] - 1
# Generate figure
fig, axs = plt.subplots(2, 1)
# Retrieve standard color sequence of matplotlib library
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
# Assign colors for lift and moment coefficients
cl_color = colors[0]
cm_color = colors[1]
# First subplot
# Set y axis label and color
axs[0].set_ylabel(r'$c_l$', color=cl_color)
axs[0].tick_params(axis='y', labelcolor=cl_color)
# Plot results for lift coefficient
axs[0].plot(np.flip(no_panels_array), np.flip(cl_alpha5), color=cl_color, marker='o', markerfacecolor='none')
axs[0].plot(np.flip(no_panels_array), np.flip(cl_alpha15), color=cl_color, marker='v', markerfacecolor='none')
# Settings for minor grid
axs[0].minorticks_on()
axs[0].grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
# Generate a right y-axis and set label and color
cm_axis = axs[0].twinx()
cm_axis.set_ylabel(r'$c_{m,c/4}$', color=cm_color)
cm_axis.tick_params(axis='y', labelcolor=cm_color)
# Plot results For moment coefficient
cm_axis.plot(np.flip(no_panels_array), np.flip(cm_alpha5), color=cm_color, marker='o', markerfacecolor='none')
cm_axis.plot(np.flip(no_panels_array), np.flip(cm_alpha15), color=cm_color, marker='v', markerfacecolor='none')
# Second subplot
# Set x-axis and y-axis label
axs[1].set_xlabel(r'\# panels')
axs[1].set_ylabel(r'Relative error')
# Plot results for relative error
axs[1].plot(np.flip(no_panels_array), np.flip(error_cl_alpha5), color=cl_color, marker='o', markerfacecolor='none')
axs[1].plot(np.flip(no_panels_array), np.flip(error_cl_alpha15), color=cl_color, marker='v', markerfacecolor='none')
axs[1].plot(np.flip(no_panels_array), np.flip(error_cm_alpha5), color=cm_color, marker='o', markerfacecolor='none')
axs[1].plot(np.flip(no_panels_array), np.flip(error_cm_alpha15), color=cm_color, marker='v', markerfacecolor='none')
# Set minor grid
axs[1].minorticks_on()
axs[1].grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
# Create custom artists
alpha5_artist = plt.Line2D((0, 1), (0, 0), color='k', marker='o', linestyle='-', markerfacecolor='none')
alpha15_artist = plt.Line2D((0, 1), (0, 0), color='k', marker='v', linestyle='-', markerfacecolor='none')
delta_cl_artist = plt.Line2D((0, 1), (0, 0), color=cl_color, linestyle='-')
delta_cm_artist = plt.Line2D((0, 1), (0, 0), color=cm_color, linestyle='-')
# Generate legends
axs[0].legend([alpha5_artist, alpha15_artist], [r'$\alpha=5^{\circ}$', r'$\alpha=15^{\circ}$'])
axs[1].legend([delta_cl_artist, delta_cm_artist], [r'$\delta c_l$', r'$\delta c_{m,c/4}$'])
# Save figure
plt.savefig('PanelDensityEffect.pdf', bbox_inches='tight')
