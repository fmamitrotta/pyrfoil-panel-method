"""Task 3 of assignment Inviscid flow over an airfoil (TU Delft's AE4130)

This script produces the plots requested for task 3 of the assignment. Comparisons between a thin and a thick airfoil
and between a symmetrical and a cambered airfoil are carried out.

The number of panels can be changed modifying the appropriate parameter.

This script requires that `numpy` and `matplotlib` be installed within the Python environment you are running this
script in.
"""

# 3rd party packages
import numpy as np
import matplotlib.pyplot as plt
# Local source
import panelairfoil as panelairfoil
from examples.ae4130assignment1 import naca0012experimentaldata as expdata


# ---------------------- Task 3.a: comparison thin vs thick airfoil ---------------------- #
# Set number of panels to be used (must be an even number)
no_panels = 70
# Generate Naca4DigitPanelled object for the NACA 0012 airfoil
naca0012 = panelairfoil.Naca4DigitPanelled('0012', no_panels)
# Generate Naca4DigitPanelled object for NACA 0024
naca0024 = panelairfoil.Naca4DigitPanelled('0024', no_panels)
# Load experimental data
abbott_force_data = expdata.abbott_force_data()
# Retrieve standard color sequence of matplotlib library
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
# Assign colors for NACA 0012 and NACA 0024
naca0012_color = colors[0]
naca0024_color = colors[2]
# Generate figure for airfoil geometry comparison
geometry_comparison_fig, geometry_comparison_axs = plt.subplots()
# Plot NACA 0012 geometry
geometry_comparison_axs.plot(np.concatenate((np.flip(naca0012.x), naca0012.x)),
                             np.concatenate((np.flip(naca0012.z_lower), naca0012.z_upper)),
                             '-', color=naca0012_color, markerfacecolor='none', label='NACA 0012')
# Plot NACA 0024 geometry
geometry_comparison_axs.plot(np.concatenate((np.flip(naca0024.x), naca0024.x)),
                             np.concatenate((np.flip(naca0024.z_lower), naca0024.z_upper)),
                             '-', color=naca0024_color, markerfacecolor='none', label='NACA 0024')
# Set aspect ratio of axes
geometry_comparison_axs.set_aspect('equal', 'box')
# Set minor grid
geometry_comparison_axs.minorticks_on()
geometry_comparison_axs.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
# Set axes label
plt.xlabel(r'$x/c$')
plt.ylabel(r'$y/c$')
# Generate legend
geometry_comparison_axs.legend(loc='upper right')
# Save figure
plt.savefig('SymmetricalAirfoilGeometryComparison.pdf', bbox_inches='tight')
# Define new list of angles of attack for the comparison of pressure distributions
alpha_array = [0, 5, 10, 15]
# Calculate pressure distributions for NACA 0012
naca0012_cp_vs_xc = dict(zip(['alpha=%d' % x for x in alpha_array],
                             [naca0012.calculate_cp_vs_xc(np.deg2rad(x)) for x in alpha_array]))
# Calculate pressure distributions for NACA 0024
naca0024_cp_vs_xc = dict(zip(['alpha=%d' % x for x in alpha_array],
                             [naca0024.calculate_cp_vs_xc(np.deg2rad(x)) for x in alpha_array]))
# Generate figure for the pressure distribution comparison between NACA 0012 and NACA 0024
cp_fig, cp_axs = plt.subplots(2, 2, sharex='all', sharey='all')
# Invert y-axis
cp_axs[0, 0].invert_yaxis()
# Flatten array of axes in order to iterate in an easier way
cp_axs_flattened = cp_axs.flatten('C')
# Iterate through the axes of the subplot
for i, axs in enumerate(cp_axs_flattened):
    # Set key to access data in dictionaries
    key = 'alpha=' + str(alpha_array[i])
    # Plot NACA 0012 pressure distribution
    axs.plot(naca0012_cp_vs_xc[key][:, 0], naca0012_cp_vs_xc[key][:, 1], '-o', color=naca0012_color,
             markerfacecolor='none', label='NACA 0012')
    # Plot NACA 0024 pressure distribution
    axs.plot(naca0024_cp_vs_xc[key][:, 0], naca0024_cp_vs_xc[key][:, 1], '-o', color=naca0024_color,
             markerfacecolor='none', label='NACA 0024')
    # Set minor grid
    axs.minorticks_on()
    axs.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
# Generate legend in the first subplot
cp_axs[0, 0].legend()
# Add a big axis, hide frame
cp_fig.add_subplot(111, frameon=False)
# Hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
# Hide grid of big axis
plt.grid(False)
# Set label of x-axis and y-axis
plt.xlabel(r'$x/c$')
plt.ylabel(r'$C_p$', labelpad=20)
# Save figure
plt.savefig('CpVsXcSymmetricalAirfoilComparison.pdf', bbox_inches='tight')
# Calculate lift curve of NACA 0012 using the same angles of attack of abbott data
naca0012_cl_vs_alpha = naca0012.calculate_cl_vs_alpha(abbott_force_data[:, 0])
# Calculate lift curve of NACA 0024 using the same angles of attack of abbott data
naca0024_cl_vs_alpha = naca0024.calculate_cl_vs_alpha(abbott_force_data[:, 0])
# Generate figure for lift curve comparison between NACA 0012 and NACA 0024
cl_fig, cl_axs = plt.subplots()
# Plot NACA 0012 lift curve
cl_axs.plot(naca0012_cl_vs_alpha[:, 0], naca0012_cl_vs_alpha[:, 1], '-o', color=naca0012_color, markerfacecolor='none',
            label='NACA 0012')
# Plot NACA 0024 lift curve
cl_axs.plot(naca0024_cl_vs_alpha[:, 0], naca0024_cl_vs_alpha[:, 1], '-o', color=naca0024_color, markerfacecolor='none',
            label='NACA 0024')
# Set minor grid
cl_axs.minorticks_on()
cl_axs.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
# Set label of x-axis and y-axis
plt.xlabel(r'$\alpha [^\circ]$')
plt.ylabel(r'$c_l$')
# Generate legend
cl_axs.legend()
# Save figure
plt.savefig('ClVsAlphaSymmetricalAirfoilComparison.pdf', bbox_inches='tight')
# ---------------------------------------------------------------------------------------- #

# --------------- Task 3.b: comparison symmetrical vs cambered airfoil --------------- #
# Generate Naca4DigitPanelled object for NACA 2412
naca2412 = panelairfoil.Naca4DigitPanelled('2412', no_panels)
# Assign color for NACA 2412
naca2412_color = colors[3]
# Generate figure for geometry comparison between NACA 0012 and NACA 2412
geometry_comparison_fig, geometry_comparison_axs = plt.subplots()
# Plot geometry of NACA 0012
geometry_comparison_axs.plot(np.concatenate((np.flip(naca0012.x), naca0012.x)),
                             np.concatenate((np.flip(naca0012.z_lower), naca0012.z_upper)),
                             '-', color=naca0012_color, markerfacecolor='none', label='NACA 0012')
# Plot geometry of NACA 2412
geometry_comparison_axs.plot(np.concatenate((np.flip(naca2412.x), naca2412.x)),
                             np.concatenate((np.flip(naca2412.z_lower), naca2412.z_upper)),
                             '-', color=naca2412_color, markerfacecolor='none', label='NACA 2412')
# Set aspect ratio of axes
geometry_comparison_axs.set_aspect('equal', 'box')
# Set minor grid
geometry_comparison_axs.minorticks_on()
geometry_comparison_axs.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
# Set label of x-axis and y-axis
plt.xlabel(r'$x/c$')
plt.ylabel(r'$y/c$')
# Generate legend
geometry_comparison_axs.legend(loc='upper right')
# Save figure
plt.savefig('CamberedAirfoilGeometryComparison.pdf', bbox_inches='tight')
# Calculate pressure distributions for NACA 2412
naca2412_cp_vs_xc = dict(zip(['alpha=%d' % x for x in alpha_array],
                             [naca2412.calculate_cp_vs_xc(np.deg2rad(x)) for x in alpha_array]))
# Generate figure for pressure distribution comparison between NACA 0012 and NACA 2412
cp_fig, cp_axs = plt.subplots(2, 2, sharex='all', sharey='all')
# Invert y-axis
cp_axs[0, 0].invert_yaxis()
# Flatten array of axes in order to iterate in an easier way
cp_axs_flattened = cp_axs.flatten('C')
# Iterate through the axes of the subplot
for i, axs in enumerate(cp_axs_flattened):
    # Set key to access data in dictionaries
    key = 'alpha=' + str(alpha_array[i])
    # Plot pressure distribution of NACA 0012
    axs.plot(naca0012_cp_vs_xc[key][:, 0], naca0012_cp_vs_xc[key][:, 1], '-o', color=naca0012_color,
             markerfacecolor='none', label='NACA 0012')
    # Plot pressure distribution of NACA 0024
    axs.plot(naca2412_cp_vs_xc[key][:, 0], naca2412_cp_vs_xc[key][:, 1], '-o', color=naca2412_color,
             markerfacecolor='none', label='NACA 2412')
    # Set minor grid
    axs.minorticks_on()
    axs.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
# Generate legend in first subplot
cp_axs[0, 0].legend()
# Add a big axis, hide frame
cp_fig.add_subplot(111, frameon=False)
# Hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
# Hide grid of big axis
plt.grid(False)
# Set label of x-axis and y-axis
plt.xlabel(r'$x/c$')
plt.ylabel(r'$C_p$', labelpad=20)
# Save figure
plt.savefig('CpVsXcCamberedAirfoilComparison.pdf', bbox_inches='tight')
# Calculate lift curve for NACA 2412
naca2412_cl_vs_alpha = naca2412.calculate_cl_vs_alpha(abbott_force_data[:, 0])
# Generate figure for lift curve comparison between NACA 0012 and NACA 2412
cl_fig, cl_axs = plt.subplots()
# Plot lift curve of NACA 0012
cl_axs.plot(naca0012_cl_vs_alpha[:, 0], naca0012_cl_vs_alpha[:, 1], '-o', color=naca0012_color, markerfacecolor='none',
            label='NACA 0012')
# Plot lift curve of NACA 2412
cl_axs.plot(naca2412_cl_vs_alpha[:, 0], naca2412_cl_vs_alpha[:, 1], '-o', color=naca2412_color, markerfacecolor='none',
            label='NACA 2412')
# Set minor grid
cl_axs.minorticks_on()
cl_axs.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
# Set label of x-axis and y-axis
plt.xlabel(r'$\alpha [^\circ]$')
plt.ylabel(r'$c_l$')
# Generate legend
cl_axs.legend()
# Save figure
plt.savefig('ClVsAlphaCamberedAirfoilComparison.pdf', bbox_inches='tight')
# ------------------------------------------------------------------------------------ #
