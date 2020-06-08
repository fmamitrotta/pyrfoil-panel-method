# pyrfoil-panel-method
Object-oriented python implementation of the panel method for airfoil calculations.

## Motivation
This project was created as a follow-up of my work for the first assignment of the course Aircraft Aerodynamics of TU Delft (AE4130), dealing with the calculation of inviscid flow over an airfoil by means of the famous panel method. The motivation behind the creation and maintenance of this project stems from the desire to have an object-oriented implementation of the panel method as opposed to the more classical functional or procedural implementation. Such object-oriented implementation is considered to give advantages in terms of modularity and flexibility, together with the convenient opportunity to link mathematical formulations to specific methods of classes representing meaningful modelling elements (e.g. panels employing a certain kind of singularity).

## Features
The current version is limited to NACA airfoils and to calculations employing a vortex with a linearly varying strength over each panel. User-defined airfoils and more types of singularity are planned advancements.

## Installation
The `panelairfoil` package, which contains the definition of all classes and functions of the framework, has not yet been uploaded to PyPi. For the moment you have to clone the repository to your system in order to use the calculation framework.

## Usage
Generate your discretized airfoil setting the 4 digit of the NACA designation and the number of panels to be used (must be an even number):
```python
import panelairfoil as panelairfoil # import panelairfoil package
naca0012 = panelairfoil.Naca4DigitPanelled('0012', 70)  # define a NACA 0012 airfoil discretized in 70 panels
```
Once the `Naca4DigitPanelled` object has been created, you can execute several operations, such as plotting the discretized geometry via the method `plot_panels` and calculating the pressure coefficients along the chord for a given angle of attack via the method `calculate_cp_vs_xc`.

The airfoil is discretized by default with panels employing a vortex with linearly varying strength. In the future the user will be free to choose among different singularities.

The examples folder provides sample scripts demonstrating the use of the object-oriented calculation framework.

## Contributing
Please don't hesistate to throw feedback and suggestions. Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License
[GPL-3.0](https://choosealicense.com/licenses/gpl-3.0/)
