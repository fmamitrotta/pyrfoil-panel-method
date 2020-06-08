"""Task 2 of assignment Inviscid flow over an airfoil (TU Delft's AE4130)

This script produces the plots requested for task 2 of the assignment. A comparison of results for NACA 0012 with
experimental data from literature is carried out.

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


# --------------------- Task 2: comparison of results for NACA 0012 with literature --------------------- #
# Set number of panels to be used (must be an even number)
no_panels = 70
# Generate Naca4DigitPanelled object for the NACA 0012 airfoil
naca0012 = panelairfoil.Naca4DigitPanelled('0012', no_panels)
# Load experimental data
gregory_pressure_data = expdata.gregory_pressure_data()
ladson_pressure_data = expdata.ladson_pressure_data()
abbott_force_data = expdata.abbott_force_data()
# Calculate pressure distribution for alpha = 0, 10, 15
alpha_array = [0, 10, 15]
naca0012_cp_vs_xc = dict(zip(['alpha=%d' % x for x in alpha_array],
                             [naca0012.calculate_cp_vs_xc(np.deg2rad(x)) for x in alpha_array]))
# Distinguish between upper and lower pressures
naca0012_cp_vs_xc_grouped = {
    'upper': {key: value[int(np.size(value, 0) / 2):, :] for key, value in naca0012_cp_vs_xc.items()},
    'lower': {key: np.flipud(value[0:int(np.size(value, 0) / 2), :])
              for key, value in naca0012_cp_vs_xc.items()}}
# Distinguish between upper and lower pressures of Ladson experimental data
ladson_data_grouped = {'upper': {key: value[int(np.size(value, 0) / 2):, :]
                                 for key, value in ladson_pressure_data.items()},
                       'lower': {key: np.flipud(value[0:int(np.size(value, 0) / 2) + 1, :])
                                 for key, value in ladson_pressure_data.items()}}
# Calculate difference w.r.t. ladson and gregory data
difference_vs_ladson_data = {key1: {key2: np.interp(ladson_data_grouped[key2][key1 + ', fixed transition'][:, 0],
                                    naca0012_cp_vs_xc_grouped[key2][key1][:, 0],
                                    naca0012_cp_vs_xc_grouped[key2][key1][:, 1], left=np.NAN,
                                    right=np.NAN) - ladson_data_grouped[key2][
                                                        key1 + ', fixed transition'][:, 1] for key2 in
                                    naca0012_cp_vs_xc_grouped} for key1 in naca0012_cp_vs_xc}
difference_vs_gregory_data = {
    key: np.interp(gregory_pressure_data[key][:, 0], naca0012_cp_vs_xc_grouped['upper'][key][:, 0],
                   naca0012_cp_vs_xc_grouped['upper'][key][:, 1], left=np.NAN, right=np.NAN) - gregory_pressure_data[
                                                            key][:, 1] for key in
    naca0012_cp_vs_xc}
# Calculate bias and random error w.r.t. ladson and gregory data
error_vs_ladson_data = {
    key1: {key2: [np.nanmean(difference_vs_ladson_data[key1][key2]), np.nanstd(difference_vs_ladson_data[key1][key2])]
           for key2 in naca0012_cp_vs_xc_grouped} for key1 in naca0012_cp_vs_xc}
error_vs_gregory_data = {
    key: [np.nanmean(difference_vs_gregory_data[key]), np.nanstd(difference_vs_gregory_data[key])]
    for key in naca0012_cp_vs_xc}
# Generate figure for pressure curves
cp_fig, cp_axs = plt.subplots(1, len(alpha_array), sharex='all', sharey='all')
cp_fig.set_size_inches(cp_fig.get_size_inches()[0] * 1.5, cp_fig.get_size_inches()[1] / 1.5)
# Invert y-axis
cp_axs[0].invert_yaxis()
# Iterate through the different axes of the subplot (so also through the angles of attack)
for i in range(len(cp_axs)):
    # Set key to access data in dictionaries
    key = 'alpha=' + str(alpha_array[i])
    # Plot gregory data
    cp_axs[i].plot(gregory_pressure_data[key][:, 0], gregory_pressure_data[key][:, 1],
                   'kv', markerfacecolor='none', label='Gregory and O\'Reilly')
    # Plot ladson data
    cp_axs[i].plot(ladson_pressure_data[key][:, 0], ladson_pressure_data[key][:, 1],
                   'k^', markerfacecolor='none', label='Ladson et al.')
    # Plot numerical data
    cp_axs[i].plot(naca0012_cp_vs_xc[key][:, 0], naca0012_cp_vs_xc[key][:, 1], '-o', markerfacecolor='none',
                   label='Present method')
    # Set minor grid for current
    cp_axs[i].minorticks_on()
    cp_axs[i].grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
# Generate legend on first subplot
cp_axs[0].legend()
# Add a big axis, hide frame
cp_fig.add_subplot(111, frameon=False)
# Hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
# Don't show grid of big axis
plt.grid(False)
# Set label of x-axis and y-axis
plt.xlabel(r'$x/c$')
plt.ylabel(r'$C_p$', labelpad=20)
# Save figure
plt.savefig('CpVsXcLiteratureComparison.pdf', bbox_inches='tight')
# Generate figure for cl vs alpha curve
cl_fig, cl_axs = plt.subplots()
# Plot abbott data
cl_axs.plot(abbott_force_data[:, 0], abbott_force_data[:, 1], 'k<', markerfacecolor='none',
            label='Abbott and von Doenhoff')
# Calculate numerical cl vs alpha curve for the same angles of attack of abbott data
naca0012_cl_vs_alpha = naca0012.calculate_cl_vs_alpha(abbott_force_data[:, 0])
# Plot cl vs alpha numerical curve
cl_axs.plot(naca0012_cl_vs_alpha[:, 0], naca0012_cl_vs_alpha[:, 1], '-o', markerfacecolor='none',
            label='Present method')
# Set minor grid
cl_axs.minorticks_on()
cl_axs.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
# Set label for x-axis and y-axis
plt.xlabel(r'$\alpha [^\circ]$')
plt.ylabel(r'$c_l$')
# Generate legend
cl_axs.legend()
# Generate right y-axis
cl_error_axis = cl_axs.twinx()
# Retrieve standard color sequence of matplotlib library
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
# Assign color to lift coefficient relative error
cl_error_color = colors[1]
# Set right y-axis label and color
cl_error_axis.set_ylabel(r'$\delta c_l$', color=cl_error_color)
cl_error_axis.tick_params(axis='y', labelcolor=cl_error_color)
# Plot relative error in lift coefficient
cl_error_axis.plot(naca0012_cl_vs_alpha[:, 0], np.divide(naca0012_cl_vs_alpha[:, 1], abbott_force_data[:, 1]) - 1,
                   color=cl_error_color,
                   marker='o', markerfacecolor='none')
# Save cl vs alpha figure
plt.savefig('ClVsAlphaLiteratureComparison.pdf', bbox_inches='tight')
# Calculate approximated slope of numerical and experimental lift curves
alpha_limit = 10
abbott_data_cl_slope_fit = np.polyfit(
    abbott_force_data[(abbott_force_data[:, 0] <= alpha_limit) & (abbott_force_data[:, 0] >= -alpha_limit), 0],
    abbott_force_data[(abbott_force_data[:, 0] <= alpha_limit) & (abbott_force_data[:, 0] >= -alpha_limit), 1],
    1)
panel_method_cl_slope_fit = np.polyfit(
    naca0012_cl_vs_alpha[(naca0012_cl_vs_alpha[:, 0] <= alpha_limit) & (naca0012_cl_vs_alpha[:, 0] >= -alpha_limit), 0],
    naca0012_cl_vs_alpha[(naca0012_cl_vs_alpha[:, 0] <= alpha_limit) & (naca0012_cl_vs_alpha[:, 0] >= -alpha_limit), 1],
    1)
# Calculate relative error between the two slopes
slope_relative_error = panel_method_cl_slope_fit[0] / abbott_data_cl_slope_fit[0] - 1
# ------------------------------------------------------------------------------------------------------- #
