
import os
import shutil
import numpy as np
from itertools import product
from ase.io import read, write
from ase.build import surface, make_supercell

def reduce(a, b):
    a_list, b_list = [a], [b]
    T_list = [np.array([[1, 0], [0, 1]])]
    while True:
        if np.dot(a, b) < 0:
            b = -b
            a_list.append(a)
            b_list.append(b)
            T_ = np.array([[1, 0], [0, -1]])
            T = np.dot(T_, T_list[-1])
            T_list.append(T)

        if np.linalg.norm(a) > np.linalg.norm(b):
            a, b = b, a
            a_list.append(a)
            b_list.append(b)
            T_ = np.array([[0, 1], [1, 0]])
            T = np.dot(T_, T_list[-1])
            T_list.append(T)
            continue

        if np.linalg.norm(b) > np.linalg.norm(b + a):
            b = b + a
            a_list.append(a)
            b_list.append(b)
            T_ = np.array([[1, 0], [1, 1]])
            T = np.dot(T_, T_list[-1])
            T_list.append(T)
            continue

        if np.linalg.norm(b) > np.linalg.norm(b - a):
            b = b - a
            a_list.append(a)
            b_list.append(b)
            T_ = np.array([[1, 0], [-1, 1]])
            T = np.dot(T_, T_list[-1])
            T_list.append(T)
            continue

        break

    return a_list, b_list, T_list


def find_int(area_0, area_1, max_area):
    ratio = area_0 / area_1
    N_0 = int(max_area / area_1) + 1
    N_1 = int(max_area / area_0) + 1

    # Initialize the variables
    n_0 = 1
    n_1 = 1
    min_diff = abs(ratio - n_0 / n_1)

    # Loop over all the possible integers and store all the possible i, j
    int_list = [[1, 1, 1]]
    
    for n_0 in range(1, N_0):
        for n_1 in range(1, N_1):
            diff = abs(ratio - n_0 / n_1)
            if diff < min_diff:
                min_diff = diff
                int_list.append([n_0, n_1, n_0 / n_1])

    return int_list, ratio

def find_ijm(N):
    ijm_list = []
    for i in range(1, N + 1):
        if N % i == 0:
            m = N // i
            for j in range(m):
                ijm_list.append([i, j, m])
    return ijm_list

def cal_mis(u_length_0, v_length_0, uv_angle_0, u_length_1, v_length_1, uv_angle_1):
    u_mis = np.abs(u_length_0 - u_length_1) / u_length_0
    v_mis = np.abs(v_length_0 - v_length_1) / v_length_0
    angle_mis = np.abs(uv_angle_0 - uv_angle_1)

    return u_mis, v_mis, angle_mis

def gen_intf(i, profile):
    hkl_0, hkl_1 = profile[0], profile[1]
    T_0 = np.array(profile[7:11]).reshape(2, 2)
    T_1 = np.array(profile[11:15]).reshape(2, 2)
    
    # Transform the 2x2 T matrix to 3x3 T matrix
    T_0 = np.vstack([T_0, [0, 0]])
    T_0 = np.hstack([T_0, [[0], [0], [1]]])
    T_1 = np.vstack([T_1, [0, 0]])
    T_1 = np.hstack([T_1, [[0], [0], [1]]])

    # Read the slab data
    cell_name_0 = f'{LOWER_CONV.split("/")[-1].split(".")[0]}'
    cell_name_1 = f'{UPPER_CONV.split("/")[-1].split(".")[0]}'
    slab_0 = read(f'output/slabs/slab_{hkl_0}_{cell_name_0}.vasp')
    slab_1 = read(f'output/slabs/slab_{hkl_1}_{cell_name_1}.vasp')

    # Transform the slab data
    slab_0 = make_supercell(slab_0, T_0, order='atom-major')
    slab_1 = make_supercell(slab_1, T_1, order='atom-major')

    # Write the transformed slab data
    save_path = f'output/interfaces/intf_{i+1}_{hkl_0}_{hkl_1}'
    os.makedirs(save_path)
    write(f'{save_path}/intf_{i+1}_slab_0_{hkl_0}.vasp', slab_0, format='vasp', direct=False)
    write(f'{save_path}/intf_{i+1}_slab_1_{hkl_1}.vasp', slab_1, format='vasp', direct=False)

    # Compare the areas of slab_0 and slab_1, let slab_0 be the lower slab with the larger area
    area_0, area_1 = profile[5], profile[6]
    reversed = False
    if area_0 < area_1:
        slab_0, slab_1 = slab_1, slab_0
        reversed = True

    # Get the thickness of slab_0, slab_1, and interface
    z_top_0 = max([atom.position[2] for atom in slab_0])
    z_bottom_0 = min([atom.position[2] for atom in slab_0])
    z_top_1 = max([atom.position[2] for atom in slab_1])
    z_bottom_1 = min([atom.position[2] for atom in slab_1])
    slab_0_thickness = z_top_0 - z_bottom_0
    slab_1_thickness = z_top_1 - z_bottom_1
    interface_thickness = slab_0_thickness + slab_1_thickness + SLAB_VACUUM * 2 + INTERFACE_GAP

    # Create the interface 
    interface = slab_0.copy()
    interface.set_cell([slab_0.cell[0], slab_0.cell[1], [0, 0, interface_thickness]])

    # Reverse the z-axis
    interface.positions[:, 2] = interface_thickness-interface.positions[:, 2]

    # Get the global coordinates of the atoms in interface and slab_1
    cell = interface.cell
    slab_1_cell = slab_1.cell
    
    cell_ = cell.copy()
    cell_[2] = slab_1_cell[2]
    slab_1_positions = slab_1.positions

    # Transform the slab_1 to the global coordinates of interface
    slab_1_positions_frac = np.dot(np.linalg.inv(slab_1_cell.T), slab_1_positions.T).T
    slab_1_positions_global = np.dot(cell_.T, slab_1_positions_frac.T).T
    slab_1.cell = cell
    slab_1.positions = slab_1_positions_global

    # Shift the interface lower
    z_bottom_interface = min([atom.position[2] for atom in interface])
    z_disp_interface = z_bottom_interface - SLAB_VACUUM
    interface.translate([0, 0, -z_disp_interface])

    # Shift the slab_1 upper
    z_top_slab_1 = max([atom.position[2] for atom in slab_1])
    z_disp_slab_1 = interface_thickness - SLAB_VACUUM - z_top_slab_1
    slab_1.translate([0, 0, z_disp_slab_1])

    # Create the interface
    interface.extend(slab_1)

    # Reverse the z-axis if needed
    if reversed:
        interface.positions[:, 2] = interface_thickness-interface.positions[:, 2]
    
    # Write the interface
    write(f'output/interfaces/intf_{i+1}_{hkl_0}_{hkl_1}.vasp', interface, format='vasp', direct=False)

    # Store the lattice matching data
    with open('output/interfaces/intf_profile.txt', 'a') as f:
        f.write(f' Interface {i+1} '.center(60, '-') + '\n')

        f.write('Total atoms:'.ljust(40) + f'{len(interface)}\n')
        f.write('Lower / Upper hkl:'.ljust(40) + f'({hkl_0}) / ({hkl_1})\n')
        f.write('Lower / Upper area (A^2):'.ljust(40) + f'{area_0:.2f} / {area_1:.2f}\n')
        f.write('\n')

        f.write('U misfit (%):'.ljust(40) + f'{profile[2] * 100:.6f}\n')
        f.write('V misfit (%):'.ljust(40) + f'{profile[3] * 100:.6f}\n')
        f.write('Angle misfit (°):'.ljust(40) + f'{profile[4]:.6f}\n')
        f.write('Area misfit (%):'.ljust(40) + f'{np.abs(area_0 - area_1) / area_0 * 100:.6f}\n')
        f.write('\n')

        f.write('Transformed matrix for lower slab:\n')
        f.write(f'{T_0[0][0]:.6f}  {T_0[0][1]:.6f}\n')
        f.write(f'{T_0[1][0]:.6f}  {T_0[1][1]:.6f}\n')
        f.write('\n')

        f.write('Transformed matrix for upper slab:\n')
        f.write(f'{T_1[0][0]:.6f}  {T_1[0][1]:.6f}\n')
        f.write(f'{T_1[1][0]:.6f}  {T_1[1][1]:.6f}\n')
        f.write('\n\n')
    
    # Store the lattice matching data in a csv file
    with open('output/interfaces/intf_profile.csv', 'a') as f:
        f.write(f'{i+1},{len(interface)},{str(hkl_0)},{str(hkl_1)},{area_0:.6f},{area_1:.6f},{profile[2]*100:.6f},{profile[3]*100:.6f},{profile[4]:.6f},{np.abs(area_0-area_1)/area_0*100:.6f},{T_0[0][0]:.6f},{T_0[0][1]:.6f},{T_0[1][0]:.6f},{T_0[1][1]:.6f},{T_1[0][0]:.6f},{T_1[0][1]:.6f},{T_1[1][0]:.6f},{T_1[1][1]:.6f}\n')

def trim(data):
    data = np.array(data)
    # Get the lattice parameters
    lattice_params = data[:, -3:]

    # Compare lattice parameters and get the unique indices
    same_idx = []
    for i in range(len(lattice_params)):
        for j in range(i+1, len(lattice_params)):
            if np.allclose(lattice_params[i], lattice_params[j]):
                same_idx.append(j)
                
    # Trim the data
    data = np.delete(data, same_idx, axis=0)
    data = data.tolist()

    return data

def slab_maker(cell_conv, miller_indices, vacuum):
    cell_name = f'{cell_conv.split("/")[-1].split(".")[0]}'

    data = []

    for h, k, l in miller_indices:
        atom = read(cell_conv)
        slab = surface(lattice=atom, indices=(h, k, l), layers=1, vacuum=vacuum, tol=1e-10, periodic=True)

        # Increase the number of layers if the slab is too short
        layers = 1
        while True:
            z_top = max([atom.position[2] for atom in slab])
            z_bottom = min([atom.position[2] for atom in slab])
            slab_length = z_top - z_bottom
            if slab_length < MIN_SLAB_LENGTH:
                layers += 1
                slab = surface(lattice=atom, indices=(h, k, l), layers=layers, vacuum=vacuum, tol=1e-10, periodic=True)
            else:
                break
        
        # Get cell parameters
        cell = slab.cell
        a, b = cell[0][:2], cell[1][:2]
        S = np.linalg.norm(np.cross(a, b))

        # Reduce the cell vectors
        a_list, b_list, T_list = reduce(a, b)
        a_r, b_r, T = a_list[-1], b_list[-1], T_list[-1]
        a_length_r, b_length_r = np.linalg.norm(a_r), np.linalg.norm(b_r)
        ab_angle_r = np.arccos(np.dot(a_r, b_r) / (a_length_r * b_length_r)) * 180 / np.pi
        
        # Store the slab data for each Miller index
        data.append([h, k, l, S, *a, *b, a_length_r, b_length_r, ab_angle_r])
        write(f'output/slabs/slab_{h}{k}{l}_{cell_name}.vasp', slab, format='vasp', direct=False)

    # Compare lattice parameters and delete the same ones
    data = trim(data)

    ''' Data format:
    0 - Miller index: hkl
    1 - Area of the slab
    2 - Cell vector ax
    3 - Cell vector ay
    4 - Cell vector bx
    5 - Cell vector by
    6 - Length of reduced cell vector a
    7 - Length of reduced cell vector b
    8 - Angle between reduced cell vectors a and b
    '''
    data = [[f'{int(i[0])}{int(i[1])}{int(i[2])}', *i[3:]] for i in data]

    return data

def pair_slabs(data_lower, data_upper, max_area):
    # Get miller indices and areas
    hkl_0 = [i[0] for i in data_lower]
    areas_0 = [i[1] for i in data_lower]
    hkl_1 = [i[0] for i in data_upper]
    areas_1 = [i[1] for i in data_upper]

    # Create the target dictionary {hkl: area}
    target_0 = dict(zip(hkl_0, areas_0))
    target_1 = dict(zip(hkl_1, areas_1))

    pairs = []
    for i, j in product(hkl_0, hkl_1):
        int_list, ratio = find_int(target_0[i], target_1[j], max_area)

        ''' Data format:
        0 - Miller index hkl of lower slab
        1 - Miller index hkl of upper slab
        2 - Area of lower slab
        3 - Area of upper slab
        4 - Ratio of the areas = area_lower / area_upper
        5 - Integer n_0
        6 - Integer n_1
        7 - Integer ratio n_0 / n_1
        '''
        pairs.append([i, j, target_0[i], target_1[j], ratio, int_list])

    return pairs

def cal_uv(data_slab, hkl, n):
    # Find the possible i, j, m
    ijm_list = find_ijm(n)

    # Get the vectors a, b
    data_hkl = [i for i in data_slab if i[0] == hkl][0]
    s = data_hkl[1]
    a = data_hkl[2:4]
    b = data_hkl[4:6]
    ab = np.array([a, b])
    
    # Calculate the supercell vectors u, v
    data = []
    for ijm in ijm_list:
        ijm_matrix = np.array([[ijm[0], ijm[1]], [0, ijm[2]]])
        supercell = np.dot(ijm_matrix, ab)
        u, v = supercell[0], supercell[1]

        S = np.linalg.norm(np.cross(u, v))
        u_r_list, v_r_list, T_list = reduce(u, v)

        # Only save the final reduced cell vectors
        u_r, v_r, T_r = u_r_list[-1], v_r_list[-1], T_list[-1]
        u_length_r, v_length_r = np.linalg.norm(u_r), np.linalg.norm(v_r)
        uv_angle_r = np.arccos(np.dot(u_r, v_r) / (u_length_r * v_length_r)) * 180 / np.pi

        # Store the data for each i, j, m
        data.append([S, n, ijm[0], ijm[1], ijm[2], *T_r.flatten(), u_length_r, v_length_r, uv_angle_r])
    
    # Compare lattice parameters and delete the same ones
    data = trim(data)

    ''' Data format:
    0 - Miller index: hkl
    1 - Area of the supercell
    2 - Scaling factor n
    3 - Integer i
    4 - Integer j
    5 - Integer m
    6 - Reduced super cell matrix T_r_1
    7 - Reduced super cell matrix T_r_2
    8 - Reduced super cell matrix T_r_3
    9 - Reduced super cell matrix T_r_4
    10 - Length of reduced super cell vector u
    11 - Length of reduced super cell vector v
    12 - Angle between reduced super cell vectors u and v
    '''
    data = [[hkl, *i] for i in data]

    return data

def lattice_match(data_pairs, data_ab_lower, data_ab_upper):
    data_matched = []
    for i, row in enumerate(data_pairs):
        hkl_0 = row[0]
        hkl_1 = row[1]

        int_list = row[5]
        for int_ in int_list:
            n_0 = int(int_[0])
            n_1 = int(int_[1])

            data_uv_lower = cal_uv(data_ab_lower, hkl_0, n_1)
            data_uv_upper = cal_uv(data_ab_upper, hkl_1, n_0)

            data_mis = []
            for j, uv_lower in enumerate(data_uv_lower):
                for k, uv_upper in enumerate(data_uv_upper):
                    mis = cal_mis(uv_lower[10], uv_lower[11], uv_lower[12], uv_upper[10], uv_upper[11], uv_upper[12])
                    data_mis.append([j, k, *mis])

            data_mis = np.array(data_mis)
            data_mis = data_mis[(data_mis[:, 2] < UV_TOL) & (data_mis[:, 3] < UV_TOL) & (data_mis[:, 4] < ANGLE_TOL)]

            for mis in data_mis:
                i, j = int(mis[0]), int(mis[1])
                S_lower = data_uv_lower[i][1]
                S_upper = data_uv_upper[j][1]
                ijm_lower = data_uv_lower[i][3:6]
                ijm_upper = data_uv_upper[j][3:6]
                t_lower = data_uv_lower[i][6:10]
                t_upper = data_uv_upper[j][6:10]

                ijm_lower = np.array([[ijm_lower[0], ijm_lower[1]], [0, ijm_lower[2]]])
                ijm_upper = np.array([[ijm_upper[0], ijm_upper[1]], [0, ijm_upper[2]]])
                t_lower = np.array([[t_lower[0], t_lower[1]], [t_lower[2], t_lower[3]]])
                t_upper = np.array([[t_upper[0], t_upper[1]], [t_upper[2], t_upper[3]]])
                T_lower = np.dot(t_lower, ijm_lower)
                T_upper = np.dot(t_upper, ijm_upper)
                
                ''' Data format:
                0 - Miller index hkl of lower slab
                1 - Miller index hkl of upper slab
                2 - u_mis
                3 - v_mis
                4 - angle_mis
                5 - Area of lower slab
                6 - Area of upper slab
                7 - Transformed matrix T1 of lower slab
                8 - Transformed matrix T2 of lower slab
                9 - Transformed matrix T3 of lower slab
                10 - Transformed matrix T4 of lower slab
                11 - Transformed matrix T1 of upper slab
                12 - Transformed matrix T2 of upper slab
                13 - Transformed matrix T3 of upper slab
                14 - Transformed matrix T4 of upper slab
                '''
                data_matched.append([hkl_0, hkl_1, *mis[2:], S_lower, S_upper, *T_lower.flatten(), *T_upper.flatten()])
    return data_matched

def find_hkl(h_max, k_max, l_max):
    hkl_list = []
    for h in range(h_max+1):
        for k in range(k_max+1):
            for l in range(l_max+1):
                if h == 0 and k == 0 and l == 0:
                    continue
                # Reduce the h, k, l using the greatest common divisor
                gcd = np.gcd.reduce([h, k, l])
                if gcd == 1:
                    hkl_list.append((h, k, l))
    return hkl_list

def main():
    LOWER_HKL, UPPER_HKL = find_hkl(H_MAX, K_MAX, L_MAX), find_hkl(H_MAX, K_MAX, L_MAX)
    print(f'\nTotal number of Miller indices: {len(LOWER_HKL)}')

    # Create slabs folder
    if not os.path.exists('output/slabs'):
        os.makedirs('output/slabs')
    else:
        shutil.rmtree('output/slabs')
        os.makedirs('output/slabs')

    # Create interfaces folder
    if not os.path.exists(f'output/interfaces'):
        os.makedirs(f'output/interfaces')
    else:
        shutil.rmtree(f'output/interfaces')
        os.makedirs(f'output/interfaces')

    with open('output/interfaces/intf_profile.txt', 'a') as f:
        f.write('+'.center(60, '+') + '\n\n')
        f.write('Interface Maker v1.0'.center(60) + '\n')
        f.write('By Guangchen Liu, gliu4@wpi.edu'.center(60) + '\n\n')
        f.write('+'.center(60, '+') + '\n\n')
        f.write(f'Miller indices considered: {len(LOWER_HKL)}'.center(60) + '\n\n')
        for i in range(0, len(LOWER_HKL), 3):
            f.write(' '.join([str(i) for i in LOWER_HKL[i:i+3]]).center(60) + '\n')
        f.write('\n')
        f.write(f'Search results for matched interfaces with area within {SEARCH_AREA[0]} - {SEARCH_AREA[1]} A^2: \n\n')
        f.write(f'{"Lower hkl":<20}{"Upper hkl":<20}{"Area (A^2)":<20}{"Found Number":<20}\n')
    
    with open('output/interfaces/intf_profile.csv', 'a') as f:
        f.write('Interface ID,Total atoms,Lower hkl,Upper hkl,Lower area,Upper area,U misfit (%),V misfit (%),Angle misfit (°),Area misfit (%),T_0_1,T_0_2,T_0_3,T_0_4,T_1_1,T_1_2,T_1_3,T_1_4\n')

    # Get the product of the Miller indices
    print('\nFinding matched interfaces...')
    data_matched_all = []
    for i, j in product(LOWER_HKL, UPPER_HKL):
        LOWER_HKL = [i]
        UPPER_HKL = [j]

        # Create slabs for lower and upper materials
        data_ab_lower = slab_maker(cell_conv=LOWER_CONV, miller_indices=LOWER_HKL, vacuum=SLAB_VACUUM)
        data_ab_upper = slab_maker(cell_conv=UPPER_CONV, miller_indices=UPPER_HKL, vacuum=SLAB_VACUUM)

        # Match the lattices of the lower and upper slabs
        max_area, upper_limit = SEARCH_AREA
        while True:
            # print(f'\nTrying to find matched interfaces for {LOWER_HKL[0]} and {UPPER_HKL[0]} with area within {max_area} A^2...')
            data_pairs = pair_slabs(data_ab_lower, data_ab_upper, max_area)
            data_matched = lattice_match(data_pairs, data_ab_lower, data_ab_upper)
            if len(data_matched) == 0:
                max_area += 10
                if max_area > upper_limit:
                    with open('output/interfaces/intf_profile.txt', 'a') as f:
                        f.write(f'{str(LOWER_HKL[0]):<20}{str(UPPER_HKL[0]):<20}{"-":<20}{"0":<20}\n')
                        print(f'\nNo matched interfaces found for {LOWER_HKL[0]} and {UPPER_HKL[0]} with area within {upper_limit} A^2')
                    break
            else:
                data_matched_all.extend(data_matched)
                with open('output/interfaces/intf_profile.txt', 'a') as f:
                    f.write(f'{str(LOWER_HKL[0]):<20}{str(UPPER_HKL[0]):<20}{max_area:<20}{len(data_matched):<20}\n')
                    print(f'\nFound {len(data_matched)} matched interfaces for {LOWER_HKL[0]} and {UPPER_HKL[0]} with area within {max_area} A^2')
                break

    with open('output/interfaces/intf_profile.txt', 'a') as f:
        f.write(f'\nTotal number of interfaces found: {len(data_matched_all)}'.center(60) + '\n\n')

    # # Write the matched interfaces
    # if not len(data_matched) == 0:
    #     print('\n'.ljust(4) + '---> Creating interfaces...')
    #     for i, profile in enumerate(data_matched):
    #         gen_intf(i, profile)
    #     print('\n'.ljust(4) + '---> All interfaces are created successfully!\n')

if __name__ == '__main__':
    # Input bulk structures, need the conventional cell
    LOWER_CONV = 'input/POSCAR_LCO_MP_R_3c_Conv.vasp'
    UPPER_CONV = 'input/POSCAR_LNO_MP_I4mmm_Conv.vasp'

    # Maximum Miller indices for h, k, l or assign the specific Miller indices for lower and upper slabs
    H_MAX, K_MAX, L_MAX = 1, 1, 1

    # Minimum length of the slab, without vacuum, in Angstrom
    MIN_SLAB_LENGTH = 10

    # Slab vacuum and interface gap, in Angstrom
    SLAB_VACUUM, INTERFACE_GAP = 10, 2

    # Search area of the interface, in A^2
    SEARCH_AREA = [100, 2000]

    # Tolerance for the misfit of lattice vectors and angles
    UV_TOL, ANGLE_TOL = 0.02, 1

    main()

