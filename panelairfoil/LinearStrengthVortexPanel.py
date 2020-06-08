"""LinearStrengthVortexPanel module

This module allows the user to work with LinearStrengthVortexPanel objects.

This module requires `numpy` to be installed within the Python environment you are using this module in.

This module contains the following classes:

    * LinearStrengthVortexPanel - class for the definition of panels with linear-strength vortex
"""

# 3rd party packages
import numpy as np


# Class for the definition of panels with linear-strength vortex
class LinearStrengthVortexPanel:
    """
    A class for the definition of panels with linear-strength vortex

    ...

    Attributes
    ----------
    x1 : float
        x-coordinate of panel's first point in global reference frame
    z1 : float
        z-coordinate of panel's first point in global reference frame
    x2 : float
        x-coordinate of panel's second point in global reference frame
    z2 : float
        z-coordinate of panel's second point in global reference frame
    gamma1: float
        vortex strength at panel's first point
    gamma2: float
        vortex strength at panel's second point
    parent_airfoil: Naca4DigitPanelled
        parent airfoil object where the panel belongs to
    x1p: float
        x-coordinate of panel's first point in local reference frame
    z1p: float
        z-coordinate of panel's first point in local reference frame
    alpha_i: float
        panel angle
    local_2_global_rotation_matrix: numpy array
        matrix to rotate from local to global reference frame
    global_2_local_rotation_matrix: numpy array
        matrix to rotate from gloabl to local reference frame

    Methods
    -------
    calculate_alpha_i()
        Calculates the panel angle
    calculate_local_2_global_rotation_matrix()
        Calculates the rotation matrix from local to global reference frame
    calculate_global_2_local_rotation_matrix()
        Calculates the rotation matrix from global to local reference frame
    point_2_local_coordinates(point_global_coordinates)
        Transforms a point's coordinates from global to local reference frame
    velocity_2_global_frame(velocity_local_frame)
        Transforms a velocity vector from local to global reference frame
    vor2dl(point_global_coordinates)
        Calculates the velocity components due to the linearly varying strength vortex panel
    qtinf_j(qinf=1, alpha=0)
        Calculates the panel-tangential component of the free-stream
    qt_j(qinf=1, alpha=0)
        Calculates the perturbation velocity at the collocation point of the panel
    delta_lift_j(qinf=1, alpha=0)
        Calculates the lift of the panel
    """

    def __init__(self, x1=None, x2=None, z1=None, z2=None, gamma1=1, gamma2=1, parent_airfoil=None):
        """
        Parameters
        ----------
        x1 : float
            x-coordinate of panel's first point in global reference frame
        x2 : float
            x-coordinate of panel's second point in global reference frame
        z1 : float
            z-coordinate of panel's first point in global reference frame
        z2 : int, optional
            z-coordinate of panel's second point in global reference frame
        gamma1 : float, optional
            vortex strength at the first point of the panel
        gamma2 : float, optional
            vortex strength at the second point of the panel
        parent_airfoil: Naca4DigitPanelled
            parent airfoil object where the panel belongs to
        """

        # First point of panel in global coordinates
        self.x1 = x1
        self.z1 = z1
        # Second point of panel in global coordinates
        self.x2 = x2
        self.z2 = z2
        # Strength of vortex at first point
        self.gamma1 = gamma1
        # Strength of vortex at second point
        self.gamma2 = gamma2
        # Airfoil object which the panel belongs to
        self.parent_airfoil = parent_airfoil
        # First point of panel in local coordinates
        self.x1p = 0
        self.z1p = 0
        # Panel angle
        self.alpha_i = self.calculate_alpha_i()
        # Matrix for rotation from local to global reference frame
        self.local_2_global_rotation_matrix = self.calculate_local_2_global_rotation_matrix()
        # Matrix for rotation from global to local reference frame
        self.global_2_local_rotation_matrix = self.calculate_global_2_local_rotation_matrix()

    @property
    def x1(self):
        return self.__x1

    @x1.setter
    def x1(self, value):
        self.__x1 = value
        if hasattr(self, 'alpha_i'):
            # If panel angle has been set once, recalculate it (and the rotation matrices) with the new value of x1
            self.alpha_i = self.calculate_alpha_i()
            self.local_2_global_rotation_matrix = self.calculate_local_2_global_rotation_matrix()
            self.global_2_local_rotation_matrix = self.calculate_global_2_local_rotation_matrix()

    @property
    def z1(self):
        return self.__z1

    @z1.setter
    def z1(self, value):
        self.__z1 = value
        if hasattr(self, 'alpha_i'):
            # If panel angle has been set once, recalculate it (and the rotation matrices) with the new value of z1
            self.alpha_i = self.calculate_alpha_i()
            self.local_2_global_rotation_matrix = self.calculate_local_2_global_rotation_matrix()
            self.global_2_local_rotation_matrix = self.calculate_global_2_local_rotation_matrix()

    @property
    def x2(self):
        return self.__x2

    @x2.setter
    def x2(self, value):
        self.__x2 = value
        if hasattr(self, 'alpha_i'):
            # If panel angle has been set once, recalculate it (and the rotation matrices) with the new value of x2
            self.alpha_i = self.calculate_alpha_i()
            self.local_2_global_rotation_matrix = self.calculate_local_2_global_rotation_matrix()
            self.global_2_local_rotation_matrix = self.calculate_global_2_local_rotation_matrix()

    @property
    def z2(self):
        return self.__z2

    @z2.setter
    def z2(self, value):
        self.__z2 = value
        if hasattr(self, 'alpha_i'):
            # If panel angle has been set once, recalculate it (and the rotation matrices) with the new value of z2
            self.alpha_i = self.calculate_alpha_i()
            self.local_2_global_rotation_matrix = self.calculate_local_2_global_rotation_matrix()
            self.global_2_local_rotation_matrix = self.calculate_global_2_local_rotation_matrix()

    @property
    # Hidden property defining the local coordinates of the second point of the panel
    def __point2_local_coordinates(self):
        return self.point_2_local_coordinates(np.array([[self.x2],
                                                        [self.z2]]))

    @property
    def x2p(self):
        return self.__point2_local_coordinates[0, 0]

    @property
    def z2p(self):
        return self.__point2_local_coordinates[1, 0]

    @property
    # Property defining the global coordinates of the control point of the panel
    def control_point_global_coordinate(self):
        return 0.5 * np.array([[(self.x1 + self.x2)],
                               [(self.z1 + self.z2)]])

    @property
    # Property defining the normal vector of the panel
    def n_i(self):
        return np.array([np.sin(self.alpha_i), np.cos(self.alpha_i)])

    @property
    # Property defining the tangential vector of the panel
    def t_i(self):
        return np.array([np.cos(self.alpha_i), -np.sin(self.alpha_i)])

    # Method for the calculation of the panel angle
    def calculate_alpha_i(self):
        # x-axis of global reference frame
        global_x_axis = np.array([1, 0])
        # x-axis of panel (in global coordinates)
        local_x_axis = np.array([self.x2 - self.x1, self.z2 - self.z1])
        # Calculate angle with cosine formula
        angle = np.arccos(np.dot(global_x_axis, local_x_axis) /
                          (np.linalg.norm(global_x_axis) * np.linalg.norm(local_x_axis)))
        # Correct angle in order to follow convention
        if self.z2 >= self.z1:
            angle = 2 * np.pi - angle
        return angle

    # Method for the calculation of the rotation matrix from local to global reference frame
    def calculate_local_2_global_rotation_matrix(self):
        return np.array([[np.cos(self.alpha_i), np.sin(self.alpha_i)],
                         [-np.sin(self.alpha_i), np.cos(self.alpha_i)]])

    # Method for the calculation of the rotation matrix from global to local reference frame
    def calculate_global_2_local_rotation_matrix(self):
        return np.array([[np.cos(self.alpha_i), -np.sin(self.alpha_i)],
                         [np.sin(self.alpha_i), np.cos(self.alpha_i)]])

    # Method for the transformation of a point's coordinates from global to local reference frame
    def point_2_local_coordinates(self, point_global_coordinates):
        multiplication_vector = point_global_coordinates - np.array([[self.x1], [self.z1]])
        local_coordinates = np.matmul(self.global_2_local_rotation_matrix,
                                      multiplication_vector)
        return local_coordinates

    # Method for the transformation of a velocity vector from local to global reference frame
    def velocity_2_global_frame(self, velocity_local_frame):
        global_frame = np.matmul(self.local_2_global_rotation_matrix,
                                 velocity_local_frame)
        return global_frame

    # Method for the calculation of the velocity components due to the linearly varying strength vortex panel
    def vor2dl(self, point_global_coordinates):
        # Transform the coordinates of the input point from global to local reference frame
        point_local_coordinates = self.point_2_local_coordinates(point_global_coordinates)
        xp = point_local_coordinates[0, 0]
        zp = point_local_coordinates[1, 0]
        # Calculate the distance between the query point and the points of the panel
        r1 = np.sqrt(xp ** 2 + zp ** 2)
        r2 = np.sqrt((xp - self.x2p) ** 2 + (zp - self.z2p) ** 2)
        # Calculate the orientation of the segments connecting the query point to the two points of the panel
        if np.linalg.norm(point_local_coordinates - self.point_2_local_coordinates(
                self.control_point_global_coordinate)) < np.finfo(float).eps:
            # If query point coincides with control point of the panel, then explicitly set the orientation in order to
            # avoid numerical issues with the sign of the angles
            theta1 = 0
            theta2 = np.pi
        else:
            theta1 = np.arctan2(zp - self.z1p, xp - self.x1p)
            theta2 = np.arctan2(zp - self.z2p, xp - self.x2p)
        # Calculate velocity components coming from the contribution of the leading singularity strength
        uap = self.gamma1 / (2 * np.pi * (self.x2p - self.x1p)) * (
                (self.x2p - xp) * (theta2 - theta1) - zp * np.log(r2 / r1))
        wap = self.gamma1 / (2 * np.pi * (self.x2p - self.x1p)) * (
                (xp - self.x2p) * np.log(r1 / r2) - (self.x2p - self.x1p) + zp * (theta2 - theta1))
        # Calculate velocity components coming from the contribution of the trailing singularity strength
        ubp = self.gamma2 / (2 * np.pi * (self.x2p - self.x1p)) * (
                (xp - self.x1p) * (theta2 - theta1) + zp * np.log(r2 / r1))
        wbp = self.gamma2 / (2 * np.pi * (self.x2p - self.x1p)) * (
                (self.x1p - xp) * np.log(r1 / r2) + (self.x2p - self.x1p) - zp * (theta2 - theta1))
        # Calculate total velocity components
        velocity_local_frame = np.array([[uap + ubp],
                                         [wap + wbp]])
        velocity_global_frame = self.velocity_2_global_frame(velocity_local_frame)
        # Return assembled array
        return np.concatenate((velocity_global_frame.T,
                               np.transpose(self.velocity_2_global_frame(np.array([[uap], [wap]]))),
                               np.transpose(self.velocity_2_global_frame(np.array([[ubp], [wbp]])))), axis=0)

    # Method for the calculation of the panel-tangential component of the free-stream
    def qtinf_j(self, qinf=1, alpha=0):
        return np.dot(np.array([qinf * np.cos(alpha), qinf * np.sin(alpha)]), self.t_i)

    # Method for the calculation of the perturbation velocity at the collocation point of the panel
    def qt_j(self, qinf=1, alpha=0):
        return self.qtinf_j(qinf, alpha) + (self.gamma1 + self.gamma2) / 4 * np.sign(np.cos(self.alpha_i))

    # Method for the calculation of the lift of the panel
    def delta_lift_j(self, rho=1, qinf=1):
        return rho * qinf * (self.gamma1 + self.gamma2) / 2 * self.x2p
