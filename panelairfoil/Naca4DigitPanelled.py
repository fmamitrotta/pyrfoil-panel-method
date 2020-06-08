"""Naca4DigitPanelled module

This module allows the user to work with Naca4DigitPanelled objects.

This module requires `numpy` and `matplotlib` to be installed within the Python environment you are using this module
in.

This module contains the following classes:

    * Naca4DigitPanelled - class representing a NACA 4-digit airfoil discretized according to the panel method
"""

# 3rd party packages
import numpy as np
import matplotlib.pyplot as plt
# Local source
from .cosspace import cosspace
from .LinearStrengthVortexPanel import LinearStrengthVortexPanel


class Naca4DigitPanelled:
    """
    A class used to represent a NACA 4-digit airfoil discretized according to the panel method

    ...

    Attributes
    ----------
    name : str
        string containing the 4 digits that identify the airfoil
    no_panels : int
        number of panels for the discretization of the airfoil (must be an even number)
    no_x_points : int
        number of chord-wise points used to calculate the airfoil coordinates
    panel_array : LinearStrengthVortexPanel
        array of LinearStrengthVortexPanel objects, representing the panels used to discretize the airfoil (in the
        future this may be substituted by an array of more general Panel objects)
    normal_velocity_influence_coefficients: numpy array
        array of the influence coefficients for the normal velocity component
    tangential_velocity_influence_coefficients: numpy array
        array of the influence coefficients for the tangential velocity component

    Methods
    -------
    generate_panels()
        Discretizes the airfoil generating an array of LinearStrengthVortexPanel objects
    calculate_influence_coefficients()
        Calculates the influence coefficients of the airfoil
    determine_vortices_strength(qinf=1, alpha=0)
        Solves the linear system of equations and finds the strengths of the vortices
    plot_panels()
        Plots the geometry of the panels discretizing the airfoil
    calculate_cp_vs_xc(alpha)
        Calculates the pressure coefficients along the chord for a given angle of attack
    calculate_cl_vs_alpha(alpha_array)
        Calculates the lift coefficient curve for a given interval of angles of attack
    calculate_cm_vs_alpha(alpha_array)
        Calculates the moment coefficient curve for a given interval of angles of attack
    """

    def __init__(self, name=None, no_panels=None, no_x_points=1001):
        """
        Parameters
        ----------
        name : str
            The 4 digit identifying the airfoil
        no_panels : int
            Number of panels used for the discretization of the airfoil
        no_x_points : int, optional
            Number of chord-wise points used to evaluate the airfoil coordinates
        """

        # Name of the airfoil, it must be a string containing the 4 digit
        self.name = name
        # Number of panels, it must be an even number
        self.no_panels = no_panels
        # Number of points for the calculation of the airfoil coordinates
        self.no_x_points = no_x_points
        # Call the generate_panels method to discretize the geometry of the airfoil into no_panels panels
        self.panel_array = self.generate_panels()
        # Calculate the normal and tangential velocity influence coefficients with the appropriate method
        self.normal_velocity_influence_coefficients, self.tangential_velocity_influence_coefficients = \
            self.calculate_influence_coefficients()
        # Initialize the array of the unknown vortex strengths
        self.__gamma_array = None

    @property
    def name(self):
        return self.__name

    @name.setter
    def name(self, value):
        self.__name = value
        if hasattr(self, 'panel_array'):
            # If the array of panels has been generated once, recalculate it with the new airfoil
            self.panel_array = self.generate_panels()
            self.influence_coefficients_array = self.calculate_influence_coefficients()

    @property
    def no_panels(self):
        return self.__no_panels

    @no_panels.setter
    def no_panels(self, value):
        self.__no_panels = value
        if hasattr(self, 'panel_array'):
            # If the array of panels has been generated once, recalculate it with the new number of panels
            self.panel_array = self.generate_panels()
            self.influence_coefficients_array = self.calculate_influence_coefficients()

    @property
    def no_x_points(self):
        return self.__no_x_points

    @no_x_points.setter
    def no_x_points(self, value):
        self.__no_x_points = value
        if hasattr(self, 'panel_array'):
            # If the array of panels has been generated once, recalculate it with the new number x points
            self.panel_array = self.generate_panels()
            self.influence_coefficients_array = self.calculate_influence_coefficients()

    @property
    def t(self):
        # Maximum thickness as a fraction of the chord
        return float(self.name[2:]) / 100

    @property
    def m(self):
        # Maximum camber as a fraction of the chord
        return float(self.name[0]) / 100

    @property
    def p(self):
        # Location of maximum camber
        return float(self.name[1]) / 10

    @property
    def x(self):
        # Definition of the non-dimensional x coordinate of the airfoil
        return cosspace(num=self.no_x_points)

    @property
    def yt(self):
        # Distribution of the non-dimensional half thickness over the chord
        return 5 * self.t * (0.2969 * np.sqrt(self.x) -
                             0.1260 * self.x -
                             0.3516 * self.x ** 2 +
                             0.2843 * self.x ** 3 -
                             0.1036 * self.x ** 4)

    @property
    def yc(self):
        # Distribution of the mean camber line over the chord
        if self.p != 0:
            # Use formula only if the location of maximum camber is different from
            # zero in order to avoid division by zero
            yc = np.concatenate(
                (self.m / self.p ** 2 * (2 * self.p * self.x[self.x <= self.p] - self.x[self.x <= self.p] ** 2),
                 self.m / (1 - self.p) ** 2 * ((1 - 2 * self.p) + 2 * self.p * self.x[self.x > self.p] -
                                               self.x[self.x > self.p] ** 2)),
                axis=0)
        else:
            # If location of maximum camber is zero, then assume symmetrical airfoil
            # and set mean camber line to zero
            yc = np.zeros(np.size(self.x))
        # Return value
        return yc

    @property
    def z_upper(self):
        return self.yc + self.yt

    @property
    def z_lower(self):
        return self.yc - self.yt

    # Method for the generation of an array of panels
    def generate_panels(self):
        # Define coordinates of the panels edges
        x_panel_edge = cosspace(num=self.no_panels / 2 + 1)
        z_panel_edge_upper = np.interp(x_panel_edge, self.x, self.z_upper)
        z_panel_edge_lower = np.interp(x_panel_edge, self.x, self.z_lower)
        # Assemble arrays containing the ordered coordinates of the panels edges
        x1_array = np.concatenate((np.flip(x_panel_edge[1:]), x_panel_edge[0:-1]))
        x2_array = np.concatenate((np.flip(x_panel_edge[0:-1]), x_panel_edge[1:]))
        z1_array = np.concatenate((np.flip(z_panel_edge_lower[1:]),
                                   z_panel_edge_upper[0:-1]))
        z2_array = np.concatenate((np.flip(z_panel_edge_lower[0:-1]),
                                   z_panel_edge_upper[1:]))
        # Generate an array of LinearStrengthVortexPanel objects
        panel_array = np.ndarray((self.no_panels,), dtype=np.object)
        for i in range(self.no_panels):
            panel_array[i] = LinearStrengthVortexPanel(
                x1_array[i], x2_array[i], z1_array[i], z2_array[i],
                parent_airfoil=self)
        # Return array
        return panel_array

    # Method for the calculation of the panels' influence coefficients
    def calculate_influence_coefficients(self):
        # Initialize influence coefficients matrix for normal velocity
        a = np.zeros((self.no_panels + 1, self.no_panels + 1))
        # Initialize influence coefficients matrix for tangential velocity
        b = np.zeros((self.no_panels, self.no_panels + 1))
        # Fill last row of influence coefficient matrix for normal velocity (known a priori)
        a[self.no_panels, 0] = 1
        a[self.no_panels, self.no_panels] = 1
        # Iterate through the collocation points
        for i in range(self.no_panels):
            # Calculate first and last influence coefficient of current collocation point,
            # corresponding to the velocity induced by the first and last vortex strengths
            a[i, 0] = np.dot(
                self.panel_array[0].vor2dl(self.panel_array[i].control_point_global_coordinate)[1, :],
                self.panel_array[i].n_i)
            a[i, self.no_panels] = np.dot(
                self.panel_array[self.no_panels - 1].vor2dl(self.panel_array[i].control_point_global_coordinate)[2, :],
                self.panel_array[i].n_i)
            b[i, 0] = np.dot(
                self.panel_array[0].vor2dl(self.panel_array[i].control_point_global_coordinate)[1, :],
                self.panel_array[i].t_i)  # *np.sign(np.cos(self.panel_array[i].alpha_i)))
            b[i, self.no_panels] = np.dot(
                self.panel_array[self.no_panels - 1].vor2dl(self.panel_array[i].control_point_global_coordinate)[2, :],
                self.panel_array[i].t_i)  # *np.sign(np.cos(self.panel_array[i].alpha_i)))
            # Iterate through the strength of the remaining vortices
            for j in range(1, self.no_panels):
                # Calculate influence coefficient
                a[i, j] = np.dot(
                    self.panel_array[j - 1].vor2dl(self.panel_array[i].control_point_global_coordinate)[2, :] +
                    self.panel_array[j].vor2dl(self.panel_array[i].control_point_global_coordinate)[1, :],
                    self.panel_array[i].n_i)
                b[i, j] = np.dot(
                    self.panel_array[j - 1].vor2dl(self.panel_array[i].control_point_global_coordinate)[2, :] +
                    self.panel_array[j].vor2dl(self.panel_array[i].control_point_global_coordinate)[1, :],
                    self.panel_array[i].t_i)
        # Return arrays
        return a, b

    # Method for the solution of the linear system
    def determine_vortices_strength(self, qinf=1, alpha=0):
        # Initialize RHS vector
        rhs = np.zeros((self.no_panels + 1, 1))
        # Iterate through the collocation points
        for i in range(self.no_panels):
            # Establish boundary condition (RHS) for current collocation point
            rhs[i] = np.dot(-np.array([qinf * np.cos(alpha), qinf * np.sin(alpha)]), self.panel_array[i].n_i)
        # Solve linear system
        self.__gamma_array = np.linalg.solve(
            self.normal_velocity_influence_coefficients, rhs)
        # Assign the calculated strength of the vortices
        for i in range(self.no_panels):
            self.panel_array[i].gamma1 = self.__gamma_array[i, 0]
            self.panel_array[i].gamma2 = self.__gamma_array[i + 1, 0]

    # Method for the visualization of the panels
    def plot_panels(self):
        # Initialize the arrays containing the coordinates of the panels' edges
        x_panel_edge = np.zeros(self.no_panels + 1)
        z_panel_edge = np.zeros(self.no_panels + 1)
        # Initialize the arrays containing the coordinates of the panels' control point
        x_control_point = np.zeros(self.no_panels)
        z_control_point = np.zeros(self.no_panels)
        # Initialize the arrays containing the components of the panels' normal vector
        normal_x_component = np.zeros(self.no_panels)
        normal_z_component = np.zeros(self.no_panels)
        # Iterate through the panels to retrieve the various coordinates and components
        for i in range(self.no_panels):
            x_panel_edge[i] = self.panel_array[i].x1
            z_panel_edge[i] = self.panel_array[i].z1
            x_control_point[i] = self.panel_array[i].control_point_global_coordinate[0, 0]
            z_control_point[i] = self.panel_array[i].control_point_global_coordinate[1, 0]
            normal_x_component[i] = self.panel_array[i].n_i[0]
            normal_z_component[i] = self.panel_array[i].n_i[1]
        # Assign the coordinates of the last edge
        x_panel_edge[-1] = self.panel_array[-1].x2
        z_panel_edge[-1] = self.panel_array[-1].z2
        # Plot panels
        plt.plot(x_panel_edge, z_panel_edge, '-|')
        # Plot control points
        plt.plot(x_control_point, z_control_point, 'x')
        # Plot normal vectors
        plt.gca().quiver(x_control_point, z_control_point,
                         normal_x_component, normal_z_component,
                         angles='xy', scale_units='xy', scale=10, color='r')
        # Adjust axes
        plt.gca().axis('equal')

    # Method for the calculation of the pressure coefficient distribution at a given angle of attack
    def calculate_cp_vs_xc(self, alpha):
        # Determine the strength of the vortices for the considered condition
        self.determine_vortices_strength(alpha=alpha)
        # Initialize x and cp vectors
        x = np.zeros((self.no_panels, 1))
        cp = np.zeros((self.no_panels, 1))
        # Determine x and cp by iteration through the panels
        for i in range(self.no_panels):
            x[i] = self.panel_array[i].control_point_global_coordinate[0, 0]
            # Perturbation velocity at the control point of the current panel
            qt_j = self.panel_array[i].qtinf_j(alpha=alpha) + np.dot(
                self.tangential_velocity_influence_coefficients[i, :], self.__gamma_array.T[0, :])
            # Pressure coefficient for the current panel
            cp[i] = 1 - qt_j ** 2
        # Return assembled array
        return np.concatenate((x, cp), axis=1)

    # Method for the calculation of the cl vs alpha curve
    def calculate_cl_vs_alpha(self, alpha_array):
        # Initialize final array
        if isinstance(alpha_array, (int, float)):
            alpha_array = [alpha_array]
        cl_vs_alpha = np.zeros((len(alpha_array), 2))
        # Iterate through the input array of angles of attack
        for i in range(np.size(cl_vs_alpha, 0)):
            self.determine_vortices_strength(alpha=alpha_array[i])
            # Calculate cp distribution for the current angle of attack
            cp_vs_x = self.calculate_cp_vs_xc(np.deg2rad(alpha_array[i]))
            # Initialize array of delta lift coefficient
            delta_cl = np.zeros(self.no_panels)
            # Determine vector normal to free-stream velocity
            n_qinf = np.array([-np.cos(np.pi / 2 - np.deg2rad(alpha_array[i])),
                               np.sin(np.pi / 2 - np.deg2rad(alpha_array[i]))])
            # Iterate through the panels to obtain each portion of lift coefficient
            for j in range(self.no_panels):
                delta_cl[j] = np.dot(-cp_vs_x[j] * self.panel_array[j].x2p * self.panel_array[j].n_i, n_qinf)
            # Store cl as the summation of all delta lift
            cl_vs_alpha[i, :] = np.array([[alpha_array[i], np.sum(delta_cl)]])
        # Return array
        return cl_vs_alpha

    # Method for the calculation of the cm vs alpha curve
    def calculate_cm_vs_alpha(self, alpha_array):
        # Initialize final array
        cm_vs_alpha = np.zeros((len(alpha_array), 2))
        for i in range(np.size(cm_vs_alpha, 0)):
            # Calculate cp distribution for the current angle of attack
            cp_vs_xc = self.calculate_cp_vs_xc(np.deg2rad(alpha_array[i]))
            # Distinguish between pressure on lower and upper surface
            cp_lower = np.flip(cp_vs_xc[0:int(np.size(cp_vs_xc, 0) / 2), 1])
            cp_upper = cp_vs_xc[int(np.size(cp_vs_xc, 0) / 2):, 1]
            # Calculate upper and lower surface derivatives
            dzu_dx = np.interp(cp_vs_xc[int(np.size(cp_vs_xc, 0) / 2):, 0], self.x, np.gradient(self.z_upper, self.x))
            dzl_dx = np.interp(cp_vs_xc[int(np.size(cp_vs_xc, 0) / 2):, 0], self.x, np.gradient(self.z_lower, self.x))
            # Interpolate upper and lower surface coordinates to panel points
            zu = np.interp(cp_vs_xc[int(np.size(cp_vs_xc, 0) / 2):, 0], self.x, self.z_upper)
            zl = np.interp(cp_vs_xc[int(np.size(cp_vs_xc, 0) / 2):, 0], self.x, self.z_lower)
            # Calculate the leading edge moment coefficient from the integration of the pressure distribution
            cm_le = np.trapz((cp_upper - cp_lower) * cp_vs_xc[int(np.size(cp_vs_xc, 0) / 2):, 0],
                             cp_vs_xc[int(np.size(cp_vs_xc, 0) / 2):, 0]) + np.trapz(
                cp_upper * dzu_dx * zu - cp_lower * dzl_dx * zl, cp_vs_xc[int(np.size(cp_vs_xc, 0) / 2):, 0])
            # Calculate the lift coefficient from the integration of the pressure distribution
            cl = self.calculate_cl_vs_alpha(np.deg2rad(alpha_array[i]))[0, 1]
            # Fill the output array with the current angle of attack and the corresponding quarter chord moment
            # coefficient
            cm_vs_alpha[i, :] = np.array([[alpha_array[i], cm_le + cl / 4]])
            # Return assembled array
        return cm_vs_alpha
