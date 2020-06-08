"""cosspace module

This module allows the user to obtain a cosine spaced array of numbers.

This module requires `numpy` to be installed within the Python environment you are using this module in.

This module contains the following functions:

    * cosspace - cosine spaced numbers over a specified interval
"""

# 3rd party packages
import numpy as np


def cosspace(start: float = 0, stop: float = 1, num: int = 50) -> np.array:
    """Return cosine spaced numbers over a specified interval."""
    # Calculate the angles corresponding to a semicircle divided into num parts
    beta = np.linspace(0, np.pi, num)
    # Find the station corresponding to the angles according to the formula in Section 11.2.1 of Katz, J., & Plotkin, A.
    # (2001). Low-Speed Aerodynamics (Cambridge Aerospace Series). Cambridge: Cambridge University Press.
    # doi:10.1017/CBO9780511810329
    return (stop - start) * 0.5 * (1 - np.cos(beta)) + start
