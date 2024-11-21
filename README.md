# Interface-Maker
A python3 code to create slabs and interfaces for first-principles calculations. 

This code is based on the A. Zur et al. paper: "Lattice match: An application to heteroepitaxy, Journal of applied physics 55(2) (1984) 378-386". 

By setting the tolerance for the misfit of lattice vectors and angles, the code can automatically find the interfaces with the desired lattice match, also ensuring the matching area / supercell size to be small enough. The code can also generate the slabs with the desired thickness and vacuum, and the interfaces with the desired gap.

⭐️ Star the repository if you find it useful :) ⭐️

## Workflow
![Workflow](res/image.png)

## Usage
1. Prepare the input POSCAR file and put it in the `input` folder. Note that the conventional cell should be used.
2. Set the parameters in the `interface_maker.py` file:

    ```python
    # Input bulk structures, need the conventional cell
    LOWER_CONV = 'input/POSCAR_LCO_MP_R_3c_Conv.vasp'
    UPPER_CONV = 'input/POSCAR_LNO_MP_I4mmm_Conv.vasp'

    # Option 1: Set maximum Miller indices of h, k, l for lower and upper slabs
    MAX_H, MAX_K, MAX_L = 1, 1, 1

    # # Option 2: Assign the specific Miller indices for lower and upper slabs
    # LOWER_HKL, UPPER_HKL = (0, 0, 1), (0, 0, 1)

    # Minimum thickness of the slab, without vacuum, in Angstrom
    MIN_SLAB_THICKNESS = 20

    # Slab vacuum and interface gap, in Angstrom
    SLAB_VACUUM, INTERFACE_GAP = 10, 2

    # Area range for the matched interfaces, in A^2
    MIN_AREA, MAX_AREA = 250, 2500

    # Tolerance for the misfit of lattice vectors and angles
    UV_TOL, ANGLE_TOL = 0.05, 5

    # Run the shape filter or not, which will only keep the near-diamond shape interfaces
    SHAPE_FILTER = True
    ```
3. Run the following command, and the output files will be saved in the `output` folder:
    ```bash
    python3 interface_maker.py
    ```
4. Check the `intf_profiles.txt` or `intf_profiles.csv` file to see the summary of the generated interfaces. Modify the parameters if needed and rerun the code to get the desired interfaces. A example of the `intf_profiles.txt` file is shown below:
    ```
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                        Interface Maker v1.0                    
                By Guangchen Liu, gliu4@wpi.edu               

    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        Miller indices considered for lower and upper slabs:    

            7 non-equivalent planes for lower slab:           

                        001 010 011 100 101                     
                            110 111                           

            5 non-equivalent planes for upper slab:           

                        001 010 011 110 111                     

    ------------------------------------------------------------

    Search results for matched interfaces with area within 2500 A^2: 

    Lower hkl           Upper hkl           Area (A^2)          
    (0, 0, 1)           (0, 0, 1)           262.5354            
    (0, 0, 1)           (0, 1, 0)           495.9002            
    (0, 0, 1)           (0, 1, 1)           466.7296            
    (0, 0, 1)           (1, 1, 0)           350.0472            
    (0, 0, 1)           (1, 1, 1)           350.0472            
    (0, 1, 0)           (0, 0, 1)           571.8600            
    (0, 1, 0)           (0, 1, 0)           1919.8158           
    (0, 1, 0)           (0, 1, 1)           776.0958            
    (0, 1, 0)           (1, 1, 0)           490.1657            
    (0, 1, 0)           (1, 1, 1)           285.9300            
    (0, 1, 1)           (0, 0, 1)           752.7283            
    (0, 1, 1)           (0, 1, 0)           1204.3653           
    (0, 1, 1)           (0, 1, 1)           1003.6377           
    (0, 1, 1)           (1, 1, 0)           351.2732            
    (0, 1, 1)           (1, 1, 1)           351.2732            
    (1, 0, 0)           (0, 0, 1)           290.2569            
    (1, 0, 0)           (0, 1, 0)           787.8402            
    (1, 0, 0)           (0, 1, 1)           870.7707            
    (1, 0, 0)           (1, 1, 0)           1326.8888           
    (1, 0, 0)           (1, 1, 1)           290.2569            
    (1, 0, 1)           (0, 0, 1)           753.1087            
    (1, 0, 1)           (0, 1, 0)           1154.7667           
    (1, 0, 1)           (0, 1, 1)           1054.3522           
    (1, 0, 1)           (1, 1, 0)           351.4507            
    (1, 0, 1)           (1, 1, 1)           351.4507            
    (1, 1, 0)           (0, 0, 1)           522.8883            
    (1, 1, 0)           (0, 1, 0)           290.4935            
    (1, 1, 0)           (0, 1, 1)           290.4935            
    (1, 1, 0)           (1, 1, 0)           1801.0596           
    (1, 1, 0)           (1, 1, 1)           290.4935            
    (1, 1, 1)           (0, 0, 1)           323.0970            
    (1, 1, 1)           (0, 1, 0)           646.1939            
    (1, 1, 1)           (0, 1, 1)           258.4776            
    (1, 1, 1)           (1, 1, 0)           581.5745            
    (1, 1, 1)           (1, 1, 1)           516.9551            
            
    Total number of interfaces found: 35            

    ---------------------- Interface 1-1 -----------------------
    Total atoms:                            1016
    Lower / Upper hkl:                      (001) / (001)
    Lower / Upper area (A^2):               262.54 / 261.77

    U misfit (%):                           0.596748
    V misfit (%):                           0.882101
    Angle misfit (°):                       0.210358
    Area misfit (%):                        0.289945

    Transformed matrix for lower slab:
    3.000000  0.000000
    0.000000  3.000000

    Transformed matrix for upper slab:
    1.000000  4.000000
    -4.000000  1.000000


    ---------------------- Interface 1-2 -----------------------
    Total atoms:                            1016
    Lower / Upper hkl:                      (001) / (001)
    Lower / Upper area (A^2):               262.54 / 261.77

    U misfit (%):                           0.596748
    V misfit (%):                           0.882101
    Angle misfit (°):                       0.210358
    Area misfit (%):                        0.289945

    Transformed matrix for lower slab:
    3.000000  0.000000
    0.000000  3.000000

    Transformed matrix for upper slab:
    1.000000  4.000000
    -4.000000  1.000000


    ---------------------- Interface 1-3 -----------------------
    Total atoms:                            1016
    Lower / Upper hkl:                      (001) / (001)
    Lower / Upper area (A^2):               262.54 / 261.77

    U misfit (%):                           0.596748
    V misfit (%):                           0.882101
    Angle misfit (°):                       0.210358
    Area misfit (%):                        0.289945

    Transformed matrix for lower slab:
    3.000000  0.000000
    0.000000  3.000000

    Transformed matrix for upper slab:
    1.000000  4.000000
    -4.000000  1.000000


    ---------------------- Interface 1-4 -----------------------
    Total atoms:                            1016
    Lower / Upper hkl:                      (001) / (001)
    Lower / Upper area (A^2):               262.54 / 261.77

    U misfit (%):                           0.596748
    V misfit (%):                           0.882101
    Angle misfit (°):                       0.210358
    Area misfit (%):                        0.289945

    Transformed matrix for lower slab:
    3.000000  0.000000
    0.000000  3.000000

    Transformed matrix for upper slab:
    1.000000  4.000000
    -4.000000  1.000000


    ---------------------- Interface 2-1 -----------------------
    Total atoms:                            1860
    Lower / Upper hkl:                      (001) / (010)
    Lower / Upper area (A^2):               495.90 / 491.31

    U misfit (%):                           3.570834
    V misfit (%):                           2.801181
    Angle misfit (°):                       2.525753
    Area misfit (%):                        0.925514

    Transformed matrix for lower slab:
    -3.000000  2.000000
    -4.000000  -3.000000

    Transformed matrix for upper slab:
    1.000000  4.000000
    2.000000  -2.000000


    ---------------------- Interface 2-2 -----------------------
    ...
    
    ```