"""panelairfoil package

This package allows the user to perform inviscid aerodynamic calculation for 2D NACA airfoil according to the panel
method.

This package requires `numpy` and `matplotlib` to be installed within the Python environment you are using this package
in.

This package contains the following modules:

    * cosspace - defines function to obtain cosine spaced numbers over a specified interval
    * LinearStrengthVortexPanel - defines LinearStrengthVortexPanel class allowing user to work with panels having
    vortex with linearly varying strength
    * Naca4DigitPanelled - defines Naca4DigitPanelled class allowing user to work with NACA 4-digit airfoils
     discretized according to the panel method
"""

# Local source
from .cosspace import cosspace
from .LinearStrengthVortexPanel import LinearStrengthVortexPanel
from .Naca4DigitPanelled import Naca4DigitPanelled
