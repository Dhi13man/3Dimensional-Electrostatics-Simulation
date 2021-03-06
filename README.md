# 3Dimensional-Electrostatics-Simulation
This is a MATLAB program to simulate Electrostatic phenomenon, particularly Coulomb's Law. We maintain a simulated 'Charge space', assumed to be an infinite homogenous insulating medium. The assumption is that every charged body in it remains stationary.

![](https://raw.githubusercontent.com/Dhi13man/3Dimensional-Electrostatics-Simulation/master/examples/100_charge_space.png)

----
Simulation Properties:
--
1. Let `charge_space` be a Structure Array with `x_coord`, `y_coord`, `z_coord` coordinates and `mag` magnitude as properties of the Structure, representing the simulated medium obeying the necessary assumptions of Coulomb’s Law.
2. Then, `length(charged_space)` gives the number of Charged particles in charge_space.
3. Let `charge_space_permittivity` represents the permittivity of the medium.
4. The `charge_space` is visualized as a 3 Dimensional Scatter plot for the Simulation.
----

# Prerequisites and Instructions:
1. Requires MATLAB for execution of CoulombsLaw.m script.
2. Run CoulombsLaw.m in MATLAB and call any of the Client functions `menu()`, `calc()`, `plot_graphs()`, `plot_space()` to use them as required.
3. Further instructions and detailed information about the functions are mentioned below.

----

# Main script CoulombsLaw.m Functions and Usage:

Client Functions that wrap the Back-end Functions (for usage):
--
A. The `plot_space()` client function **plots the entirety of the charge_space and all the charges present in it in 3 Dimensional Space as a Scatter plot**. 

B. The `menu()` client function creates the **"Charge Space Menu"** and can be called to 
1. **Add charged bodies** into `charge_space`. Positive Charges are Red and Negative are Blue.
2. **Remove charged bodies** from `charge_space`. Or 6. **Clear entirety** of the `charge_space`.
3. **Generate N random charged bodies** in `charge_space`.
4. **Display information** i.e. Location and Magnitude of any charged body.
5. **Plot the charge_space** as a 3-Dimensional scatter plot and **Exit** out of **Charge Space Menu**.

C. The `calc()` client function creates the **"Calculate Values with Coulomb's Law Menu"** and can be called to
1. **Change permittivity**, i.e `charge_space_permittivity` of `charge_space`.
2. **Calculate net force on any charge** in `charge_space`.
3 **Calculate force between any two charges** in `charge_space`.
4. **Calculate force between any two temporary charges** defined by parameters.
5. **Calculate net Electric Field strength** in any point in `charge_space`.
6. **Exit** out of **Calculate Values with Coulomb's Law Menu**.

D. The `plot_graphs()` client function can be used to **graphically test the relationship of the parameters with the force**.

Back-end Functions:
--
**A. Charge Space Manipulation**
1. `clear_charge_space()`: Removes all charges from Charge Space.
2. `charge_place()`: Creates a charge at coordinates pointed by parameters `x, y, z` in charge space of magnitude given by parameter `mag` and plots it if parameter `show` is 1 or [].
3. `display_charge()`: Display location and magnitude of charge given by parameter `obj` in charge space. `obj` can be numeric or Structure.
4. `charge_destroy()`: Removes charge given by parameter `obj` from charge space. `obj` can be numeric or Structure.
5. `n_random_charges()`: Generate number of random charges as given by parameter `number`.

**B. Locating Charges**
1. `get_chargenum_bycoord()`: Get the number of the Charge pointed by given coordinate parameters `x, y, z` in charge_space.
2. `get_charge_bynum()`: Get Charge object pointed by parameter `charge_number` in charge_space. 

**C. Finding physical Force and Electric Field Strength values**
1. `charge_pair_force()`: Calculates the three Dimensional components of force on charge parameter `obj` by charge parameter `obj2`.
2. `charge_net_force()`: Superposition Theorem: Calculates the three dimensional components of the net force on charge parameter `obj` by all other Charges.
3. `net_field_at()`: Superposition Theorem: Calculates the three dimensional components of the net Electric Field strength at coordinate parameters `x, y, z` in Charge Space.

----

# Repository Folder Contents:
**1. CoulombsLaw.m:** The main MATLAB script containing all the utility and backend functions to make the Simulation work. Functions are explained above in detail.\
**2. CoulombsLaw.pptx** Presentation prepared for the purpose of this project. Theoretical concepts explained from Ground Zero.
**3. EMFT_Project.prj:** Contains the project structure information. Can be ignored entirely during application.\
**4. examples folder:** Contains example images of the things the script can accomplish.\
**5. resources folder:** Contains the Project resources. Do not modify without expertise.\
**6. README.md:** This Github Repository information file.\
**7. LICENSE:** BSD 3 License information that governs the Software. Do not modify.

----

Examples:
--
![](https://github.com/Dhi13man/3Dimensional-Electrostatics-Simulation/blob/master/examples/3_charge_space.png)\
**Fig:** Simulation of a 3D Charge Space in MATLAB, where 3 charges were generated in
1. (10, 0, 20); Magnitude = 100 C
2. (-10, 0, 70); Magnitude = 1000 C
3. (5, 12, -7); Magnitude = -20 C

----
![](https://raw.githubusercontent.com/Dhi13man/3Dimensional-Electrostatics-Simulation/master/examples/net_field_calc_3_charge_space.png)\
**Fig:** Net Electric Field strength calculated in Charge Space at (0, 0, 0), and verified to be correct.

----
![](https://raw.githubusercontent.com/Dhi13man/3Dimensional-Electrostatics-Simulation/master/examples/net_force_calc_3_charge_space.png)\
**Fig:** Net Force calculated on Charge 3 (5, 12, -7) by the other two Charges, and verified to be correct.

----
