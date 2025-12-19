

import streamlit as st
import io
import zipfile
import datetime
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


def find_int(area_0, area_1, area):
    ratio = area_0 / area_1
    N_0 = int(area / area_1) + 1
    N_1 = int(area / area_0) + 1

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


def gen_intf(i, profile, slabs_lower, slabs_upper, SLAB_VACUUM, INTERFACE_GAP, profile_txt, profile_csv, zip_buffer):
    hkl_0, hkl_1 = profile[0], profile[1]
    T_0 = np.array(profile[7:11]).reshape(2, 2)
    T_1 = np.array(profile[11:15]).reshape(2, 2)
    
    # Transform the 2x2 T matrix to 3x3 T matrix
    T_0 = np.vstack([T_0, [0, 0]])
    T_0 = np.hstack([T_0, [[0], [0], [1]]])
    T_1 = np.vstack([T_1, [0, 0]])
    T_1 = np.hstack([T_1, [[0], [0], [1]]])

    # Read the slab data
    slab_0 = slabs_lower[hkl_0]
    slab_1 = slabs_upper[hkl_1]

    # Transform the slab data
    slab_0 = make_supercell(slab_0, T_0, order='atom-major')
    slab_1 = make_supercell(slab_1, T_1, order='atom-major')

    # Create slab_0_reverse and slab_1_reverse for the reversed z-axis
    slab_0_reverse = slab_0.copy()
    slab_1_reverse = slab_1.copy()
    slab_0_reverse.positions[:, 2] = slab_0.cell[2, 2] - slab_0.positions[:, 2]
    slab_1_reverse.positions[:, 2] = slab_1.cell[2, 2] - slab_1.positions[:, 2]

    # Helper to write ASE atoms object to the zip buffer
    def add_to_zip(structure, filename):
        vasp_buffer = io.StringIO()
        write(vasp_buffer, structure, format='vasp', direct=False, sort=True)
        zip_file.writestr(filename, vasp_buffer.getvalue())
    
    with zipfile.ZipFile(zip_buffer, 'a', zipfile.ZIP_DEFLATED) as zip_file:
        base_name = f'intf_{i+1}_{hkl_0}_{hkl_1}'
        add_to_zip(slab_0, f'{base_name}/slab_0_{hkl_0}.vasp')
        add_to_zip(slab_1, f'{base_name}/slab_1_{hkl_1}.vasp')

    comb = [(slab_0, slab_1), (slab_0_reverse, slab_1), (slab_0, slab_1_reverse), (slab_0_reverse, slab_1_reverse)]
    for j, (slab_0, slab_1) in enumerate(comb):

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
        
        # Write the interface to the zip buffer
        with zipfile.ZipFile(zip_buffer, 'a', zipfile.ZIP_DEFLATED) as zip_file:
            add_to_zip(interface, f'intf_{i+1}_{j+1}_{hkl_0}_{hkl_1}.vasp')

        # Store the lattice matching data
        profile_txt.write(f' Interface {i+1}-{j+1} '.center(60, '-') + '\n')

        profile_txt.write('Total atoms:'.ljust(40) + f'{len(interface)}\n')
        profile_txt.write('Lower / Upper hkl:'.ljust(40) + f'({hkl_0}) / ({hkl_1})\n')
        profile_txt.write('Lower / Upper area (A^2):'.ljust(40) + f'{area_0:.2f} / {area_1:.2f}\n')
        profile_txt.write('\n')

        profile_txt.write('U misfit (%):'.ljust(40) + f'{profile[2] * 100:.6f}\n')
        profile_txt.write('V misfit (%):'.ljust(40) + f'{profile[3] * 100:.6f}\n')
        profile_txt.write('Angle misfit (°):'.ljust(40) + f'{profile[4]:.6f}\n')
        profile_txt.write('Area misfit (%):'.ljust(40) + f'{np.abs(area_0 - area_1) / area_0 * 100:.6f}\n')
        profile_txt.write('\n')

        profile_txt.write('Transformed matrix for lower slab:\n')
        profile_txt.write(f'{T_0[0][0]:.6f}  {T_0[0][1]:.6f}\n')
        profile_txt.write(f'{T_0[1][0]:.6f}  {T_0[1][1]:.6f}\n')
        profile_txt.write('\n')

        profile_txt.write('Transformed matrix for upper slab:\n')
        profile_txt.write(f'{T_1[0][0]:.6f}  {T_1[0][1]:.6f}\n')
        profile_txt.write(f'{T_1[1][0]:.6f}  {T_1[1][1]:.6f}\n')
        profile_txt.write('\n\n')
        
        # Store the lattice matching data in a csv file
        profile_csv.write(f'{i+1},{j+1},{len(interface)},{str(hkl_0)},{str(hkl_1)},{area_0:.6f},{area_1:.6f},{profile[2]*100:.6f},{profile[3]*100:.6f},{profile[4]:.6f},{np.abs(area_0-area_1)/area_0*100:.6f},{T_0[0][0]:.6f},{T_0[0][1]:.6f},{T_0[1][0]:.6f},{T_0[1][1]:.6f},{T_1[0][0]:.6f},{T_1[0][1]:.6f},{T_1[1][0]:.6f},{T_1[1][1]:.6f}\n')

    return zip_buffer, profile_txt, profile_csv


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

    return data, same_idx


def slab_maker(cell_conv, miller_indices, vacuum, MIN_THICKNESS):
    # cell_name = f'{cell_conv.split("/")[-1].split(".")[0]}'

    data = []
    slabs = []
    slabs_store = {}

    for h, k, l in miller_indices:
        # atom = read(cell_conv)
        atom = cell_conv
        slab = surface(lattice=atom, indices=(h, k, l), layers=1, vacuum=vacuum, tol=1e-10, periodic=True)

        # Increase the number of layers if the slab is too short
        layers = 1
        while True:
            z_top = max([atom.position[2] for atom in slab])
            z_bottom = min([atom.position[2] for atom in slab])
            slab_length = z_top - z_bottom
            if slab_length < MIN_THICKNESS:
                layers += 1
                slab = surface(lattice=atom, indices=(h, k, l), layers=layers, vacuum=vacuum, tol=1e-10, periodic=True)
            else:
                break
        
        # Get cell parameters
        cell = slab.cell
        a, b = cell[0][:2], cell[1][:2]
        a_ = np.array([cell[0][0], cell[0][1], 0.0])
        b_ = np.array([cell[1][0], cell[1][1], 0.0])
        S = np.linalg.norm(np.cross(a_, b_))

        # Reduce the cell vectors
        a_list, b_list, T_list = reduce(a, b)
        a_r, b_r, T = a_list[-1], b_list[-1], T_list[-1]
        a_length_r, b_length_r = np.linalg.norm(a_r), np.linalg.norm(b_r)
        ab_angle_r = np.arccos(np.dot(a_r, b_r) / (a_length_r * b_length_r)) * 180 / np.pi
        
        # Store the slab data for each Miller index
        data.append([h, k, l, S, *a, *b, a_length_r, b_length_r, ab_angle_r])
        slabs.append(slab)

    # Compare lattice parameters and delete the same ones
    data, same_idx = trim(data)
    
    # Write the slabs
    for i, slab in enumerate(slabs):
        h, k, l = miller_indices[i]
        if i not in same_idx:
            slabs_store[f'{int(h)}{int(k)}{int(l)}'] = slab

    data = [[f'{int(i[0])}{int(i[1])}{int(i[2])}', *i[3:]] for i in data]

    return data, slabs_store


def pair_slabs(data_lower, data_upper, area):
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
        int_list, ratio = find_int(target_0[i], target_1[j], area)

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

        u_ = np.array([supercell[0][0], supercell[0][1], 0.0])
        v_ = np.array([supercell[1][0], supercell[1][1], 0.0])

        S = np.linalg.norm(np.cross(u_, v_))
        u_r_list, v_r_list, T_list = reduce(u, v)

        # Only save the final reduced cell vectors
        u_r, v_r, T_r = u_r_list[-1], v_r_list[-1], T_list[-1]
        u_length_r, v_length_r = np.linalg.norm(u_r), np.linalg.norm(v_r)
        uv_angle_r = np.arccos(np.dot(u_r, v_r) / (u_length_r * v_length_r)) * 180 / np.pi

        # Store the data for each i, j, m
        data.append([S, n, ijm[0], ijm[1], ijm[2], *T_r.flatten(), u_length_r, v_length_r, uv_angle_r])
    
    # Compare lattice parameters and delete the same ones
    data, same_idx = trim(data)

    data = [[hkl, *i] for i in data]

    return data


def lattice_match(data_pairs, data_ab_lower, data_ab_upper, UV_TOL, ANGLE_TOL):
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
            data_mis = data_mis[(data_mis[:, 2] < (UV_TOL / 100)) & (data_mis[:, 3] < (UV_TOL / 100)) & (data_mis[:, 4] < ANGLE_TOL)]

            for mis in data_mis:
                i, j = int(mis[0]), int(mis[1])
                S_lower = data_uv_lower[i][1]
                S_upper = data_uv_upper[j][1]
                ijm_lower = data_uv_lower[i][3:6]
                ijm_upper = data_uv_upper[j][3:6]
                t_lower = data_uv_lower[i][6:10]
                t_upper = data_uv_upper[j][6:10]
                u_lower = data_uv_lower[i][10]
                v_lower = data_uv_lower[i][11]
                angle_lower = data_uv_lower[i][12]
                u_upper = data_uv_upper[j][10]
                v_upper = data_uv_upper[j][11]
                angle_upper = data_uv_upper[j][12]

                ijm_lower = np.array([[ijm_lower[0], ijm_lower[1]], [0, ijm_lower[2]]])
                ijm_upper = np.array([[ijm_upper[0], ijm_upper[1]], [0, ijm_upper[2]]])
                t_lower = np.array([[t_lower[0], t_lower[1]], [t_lower[2], t_lower[3]]])
                t_upper = np.array([[t_upper[0], t_upper[1]], [t_upper[2], t_upper[3]]])
                T_lower = np.dot(t_lower, ijm_lower)
                T_upper = np.dot(t_upper, ijm_upper)
                
                data_matched.append([hkl_0, hkl_1, *mis[2:], S_lower, S_upper, *T_lower.flatten(), *T_upper.flatten(), u_lower, v_lower, angle_lower, u_upper, v_upper, angle_upper])
    
    return data_matched


def filter_data(data_matched, MIN_AREA, SHAPE_FILTER):
    if len(data_matched) > 1:
        data_matched = np.array(data_matched)

        # # Compare the areas and get the closest area to the MIN_AREA
        # area_diff = data_matched[:, 5].astype(float) - MIN_AREA
        # area_diff_sorted = np.sort(area_diff)
        # # Get the idx of the first positive value in area_diff_sorted
        # idx_0 = np.where(area_diff_sorted > 0)[0][0]
        # idxes = np.where(area_diff == area_diff_sorted[idx_0])[0]
        # data_matched = data_matched[idxes]

        if SHAPE_FILTER:
            # Compare the u, v lengths and calculate the uv_ratio = u / v
            u, v = data_matched[:, 15], data_matched[:, 16]
            # Change the u, v to float type and calculate the uv_ratio
            u, v = u.astype(float), v.astype(float)
            uv_ratio = np.abs(u / v - 1)
            # Filter the data using the min uv_ratio
            min_idx = np.argmin(uv_ratio)
            data_matched = [data_matched[min_idx].tolist()]
            # Change items in data_matched to float type
            for i in range(2, len(data_matched[0])):
                data_matched[0][i] = float(data_matched[0][i])
        else:
            data_matched = data_matched.tolist()
            # Change items in data_matched to float type
            for data in data_matched:
                for i in range(2, len(data)):
                    data[i] = float(data[i])
        
        min_area = data_matched[0][5]
        return data_matched, min_area
    else:
        min_area = data_matched[0][5]
        return data_matched, min_area


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


def interface_maker(session_state):
    log_container = st.empty()
    log_buffer = ''
    if session_state.assign_specific:
        LOWER_HKL, UPPER_HKL = session_state.LOWER_HKL, session_state.UPPER_HKL
        log_buffer += f'Assigned Miller indices:\nLower slab: {LOWER_HKL} \nupper slab: {UPPER_HKL}'
    else:
        MAX_H, MAX_K, MAX_L = session_state.MAX_H, session_state.MAX_K, session_state.MAX_L
        LOWER_HKL, UPPER_HKL = find_hkl(MAX_H, MAX_K, MAX_L), find_hkl(MAX_H, MAX_K, MAX_L)
        log_buffer += f'Assigned Miller indices:\nLower slab: {LOWER_HKL} \nupper slab: {UPPER_HKL}'
    # Display the log
    with log_container:
        st.text_area(
            label='Interface Generation Details', 
            value=log_buffer,
            height=600,
            )

    # Unpack the session state variables
    timestamp = session_state.timestamp
    LOWER_CONV = session_state.LOWER_CONV
    UPPER_CONV = session_state.UPPER_CONV
    UV_TOL = session_state.UV_TOL
    ANGLE_TOL = session_state.ANGLE_TOL
    MIN_THICKNESS = session_state.MIN_THICKNESS
    INTERFACE_VACUUM = session_state.INTERFACE_VACUUM
    INTERFACE_GAP = session_state.INTERFACE_GAP
    MIN_AREA = session_state.MIN_AREA
    MAX_AREA = session_state.MAX_AREA
    SHAPE_FILTER = session_state.SHAPE_FILTER

    # Create slabs for lower and upper materials
    data_ab_lower, slabs_lower = slab_maker(cell_conv=LOWER_CONV, miller_indices=LOWER_HKL, vacuum=INTERFACE_VACUUM, MIN_THICKNESS=MIN_THICKNESS)
    data_ab_upper, slabs_upper = slab_maker(cell_conv=UPPER_CONV, miller_indices=UPPER_HKL, vacuum=INTERFACE_VACUUM, MIN_THICKNESS=MIN_THICKNESS)
    lower_khl, upper_khl = [i[0] for i in data_ab_lower], [i[0] for i in data_ab_upper]
    area_lower, area_upper = [i[1] for i in data_ab_lower], [i[1] for i in data_ab_upper]

    profile_txt = io.StringIO()
    profile_txt.write('+'.center(60, '+') + '\n\n')
    profile_txt.write('Interface Maker v1.0'.center(60) + '\n')
    profile_txt.write('By Guangchen Liu, gliu4@wpi.edu'.center(60) + '\n\n')
    profile_txt.write('+'.center(60, '+') + '\n\n')
    if not session_state.assign_specific:
        profile_txt.write('Miller indices considered for lower and upper slabs:'.center(60) + '\n\n')
        profile_txt.write(f'{len(lower_khl)} non-equivalent planes for lower slab:'.center(60) + '\n\n')
        for i in range(0, len(lower_khl), 5):
            profile_txt.write(' '.join([f'{str(j)}' for j in lower_khl[i:i+5]]).center(60) + '\n')
        profile_txt.write('\n')
        profile_txt.write(f'{len(upper_khl)} non-equivalent planes for upper slab:'.center(60) + '\n\n')
        for i in range(0, len(upper_khl), 5):
            profile_txt.write(' '.join([f'{str(j)}' for j in upper_khl[i:i+5]]).center(60) + '\n')
    else:
        profile_txt.write(f'Aassigned Miller indices:'.center(60) + '\n')
        profile_txt.write(f'Lower slab: {LOWER_HKL[0]}'.center(60) + '\n')
        profile_txt.write(f'Upper slab: {UPPER_HKL[0]}'.center(60) + '\n')
    profile_txt.write('\n')
    profile_txt.write('-'.center(60, '-') + '\n\n')
    if SHAPE_FILTER:
        profile_txt.write('Warning: Shape filter is ON! '.center(60) + '\n')
        profile_txt.write('Only the most diamond-like interface will be kept!'.center(60) + '\n')
    else:
        profile_txt.write('Warning: Shape filter is OFF! '.center(60) + '\n')
        profile_txt.write('All matched interfaces will be kept!'.center(60) + '\n')
    profile_txt.write('\n')
    profile_txt.write('-'.center(60, '-') + '\n\n')
    profile_txt.write(f'Search results for matched interfaces with area within {MAX_AREA} A^2: \n\n')
    profile_txt.write(f'{"Lower hkl":<20}{"Upper hkl":<20}{"Area (A^2)":<20}\n')

    profile_csv = io.StringIO()
    profile_csv.write('Interface ID,Surface ID,Total atoms,Lower hkl,Upper hkl,Lower area,Upper area,U misfit (%),V misfit (%),Angle misfit (°),Area misfit (%),T_0_1,T_0_2,T_0_3,T_0_4,T_1_1,T_1_2,T_1_3,T_1_4\n')

    # Get the product of the Miller indices
    log_buffer += '\n\nFinding matched interfaces...'
    # Display the log
    with log_container:
        st.text_area(
            label='Interface Generation Details', 
            value=log_buffer,
            height=600,
            )
    data_matched_all = []
    for i, j in product(LOWER_HKL, UPPER_HKL):
        LOWER_HKL = [i]
        UPPER_HKL = [j]

        # Get the slab data in data_ab_lower and data_ab_upper using the Miller indices
        lower_hkl = f'{LOWER_HKL[0][0]}{LOWER_HKL[0][1]}{LOWER_HKL[0][2]}'
        upper_hkl = f'{UPPER_HKL[0][0]}{UPPER_HKL[0][1]}{UPPER_HKL[0][2]}'
        data_ab_lower_hkl = [i for i in data_ab_lower if i[0] == lower_hkl]
        data_ab_upper_hkl = [i for i in data_ab_upper if i[0] == upper_hkl]

        if len(data_ab_lower_hkl) == 1 and len(data_ab_upper_hkl) == 1:
            # Match the lattices of the lower and upper slabs
            while True:
                data_pairs = pair_slabs(data_ab_lower_hkl, data_ab_upper_hkl, MAX_AREA)
                data_matched = lattice_match(data_pairs, data_ab_lower_hkl, data_ab_upper_hkl , UV_TOL, ANGLE_TOL)
                if len(data_matched) == 0:
                    profile_txt.write(f'{str(LOWER_HKL[0]):<20}{str(UPPER_HKL[0]):<20}{"-":<20}{"0":<20}\n')
                    log_buffer += '\n'.ljust(4) + f'---> No matched interfaces found for {LOWER_HKL[0]} and {UPPER_HKL[0]} within {MAX_AREA} Å^2'
                    # Display the log
                    with log_container:
                        st.text_area(
                            label='Interface Generation Details',
                            value=log_buffer,
                            height=600,
                            )
                    break
                else:
                    data_matched, min_area = filter_data(data_matched, MIN_AREA, SHAPE_FILTER)
                    data_matched_all.extend(data_matched)
                    profile_txt.write(f'{str(LOWER_HKL[0]):<20}{str(UPPER_HKL[0]):<20}{min_area:<20.4f}\n')
                    log_buffer += '\n'.ljust(4) + f'---> Found {len(data_matched)} matched interfaces for {LOWER_HKL[0]} and {UPPER_HKL[0]} within {min_area:.4f} Å^2'
                    # Display the log
                    with log_container:
                        st.text_area(
                            label='Interface Generation Details', 
                            value=log_buffer,
                            height=600,
                            )
                    break

    profile_txt.write(f'\nTotal number of interfaces found: {len(data_matched_all)}'.center(60) + '\n\n')
    log_buffer += '\n'.ljust(4) + f'---> Total number of interfaces found: {len(data_matched_all)}'
    # Display the log
    with log_container:
        st.text_area(
            label='Interface Generation Details', 
            value=log_buffer,
            height=600,
            )

    # Write the matched interfaces
    if not len(data_matched_all) == 0:
        log_buffer += '\n\nGenerating interfaces...'
        # Display the log
        with log_container:
            st.text_area(
                label='Interface Generation Details', 
                value=log_buffer,
                height=600,
                )
        # Create in-memory zip buffer
        zip_buffer = io.BytesIO()
        for i, profile in enumerate(data_matched_all):
            zip_buffer, profile_txt, profile_csv = gen_intf(i, profile, slabs_lower, slabs_upper, INTERFACE_VACUUM, INTERFACE_GAP, profile_txt, profile_csv, zip_buffer)
            
            log_buffer += '\n'.ljust(4) + f'---> Generated interface {i+1} for ({profile[0][0]}, {profile[0][1]}, {profile[0][2]}) and ({profile[1][0]},  {profile[1][1]}, {profile[1][2]})...'
            # Display the log
            with log_container:
                st.text_area(
                    label='Interface Generation Details', 
                    value=log_buffer,
                    height=600,
                    )

        # Write the profile txt and csv to the zip buffer
        with zipfile.ZipFile(zip_buffer, 'a', zipfile.ZIP_DEFLATED) as zip_file:
            zip_file.writestr(f'intf_profile.txt', profile_txt.getvalue())
            zip_file.writestr(f'intf_profile.csv', profile_csv.getvalue())
        
        # Close the profile txt and csv
        profile_txt.close()
        profile_csv.close()

    return zip_buffer


def main():
    # Set the page config
    st.set_page_config(
        page_title='Interface-Maker', 
        layout='wide', 
        page_icon=':material/apps:', 
        menu_items={
        'Get Help': 'https://github.com/aguang5241/Interface-Maker',
        'Report a bug': 'mailto:gliu4@wpi.edu',
        'About': '# Interface-Maker  '
        '\n**Design Customize & Simulate**  '
        '\n\n*Developed by Guangchen Liu*  '
        '\n*IMPD Group, Worcester Polytechnic Institute, MA USA*',
    })
    
    # Add app logo to the sidebar
    st.sidebar.image('res/logo.png', width='stretch')
    # Set the sidebar title
    st.sidebar.title('Interface-Maker  [![GitHub stars](https://img.shields.io/github/stars/aguang5241/Interface-Maker?style=social)](https://github.com/aguang5241/Interface-Maker)')
    # Add a description to the sidebar
    st.sidebar.markdown('An application for generating customizable slabs and interfaces for first-principles simulations.', unsafe_allow_html=True)
    # Add a citation link to the sidebar
    st.sidebar.divider()
    st.sidebar.markdown('If you find this application useful, please consider citing our publication:')
    st.sidebar.markdown('[![DOI](https://img.shields.io/badge/DOI-10.1016/j.mtphys.2025.101940-blue)](https://doi.org/10.1016/j.mtphys.2025.101940)')
    # Add contact information: gliu4@wpi.edu
    # st.sidebar.divider()
    st.sidebar.markdown('For any questions or suggestions, please contact:')
    st.sidebar.markdown('[![Email](https://img.shields.io/badge/Email-yzhong@wpi.edu-white?logo=mail.ru&logoColor=white)](mailto:yzhong@wpi.edu)')
    st.sidebar.markdown('[![Email](https://img.shields.io/badge/Email-gliu4@wpi.edu-white?logo=mail.ru&logoColor=white)](mailto:gliu4@wpi.edu)')
    st.sidebar.markdown('[![LinkedIn](https://img.shields.io/badge/LinkedIn-Guangchen%20Liu-white?logo=linkedin&logoColor=white)](https://www.linkedin.com/in/aguang5241)')
    
    # Add a title to the main page
    st.title('Interface-Maker  [![GitHub stars](https://img.shields.io/github/stars/aguang5241/Interface-Maker?style=social)](https://github.com/aguang5241/Interface-Maker)')
    # Add a description to the main page
    st.markdown('An application for generating customizable slabs and interfaces for first-principles simulations.', unsafe_allow_html=True)

    # Add a section for uploading files of lower and upper systems
    st.divider()
    st.subheader('Upload Your Structures')
    st.caption('Upload the conventional cell of the lower and upper systems in VASP format (.vasp/.poscar).')
    col1, col2 = st.columns(2, gap='medium', border=True)
    with col1:
        lower_input = st.file_uploader('Lower system', type=['vasp', 'POSCAR'])
    with col2:
        upper_input = st.file_uploader('Upper system', type=['vasp', 'POSCAR'])

    # Check if the files are uploaded
    st.session_state.files_uploaded = False
    if lower_input is not None and upper_input is not None:
        try:
            # Read the uploaded files
            lower_content = lower_input.getvalue().decode("utf-8")
            lower_content = io.StringIO(lower_content)
            LOWER_CONV = read(lower_content, format='vasp')
            upper_content = upper_input.getvalue().decode("utf-8")
            upper_content = io.StringIO(upper_content)
            UPPER_CONV = read(upper_content, format='vasp')
            # Set the files_uploaded flag to True
            st.session_state.files_uploaded = True
            st.success('Files uploaded successfully!')
            st.write('Lower System:', LOWER_CONV)
            st.write('Upper System:', UPPER_CONV)
            # Set session_state variables
            st.session_state.LOWER_CONV = LOWER_CONV
            st.session_state.UPPER_CONV = UPPER_CONV
        except Exception as e:
            st.error(f'Sorry, we cannot read the uploaded files. Please check the file format and try again.')
            st.stop()
    
    if st.session_state.files_uploaded:
        # Add a section for Miller Indices
        st.divider()
        st.subheader('Define Miller Indices')
        st.caption('Define the maximum Miller indices of h, k, l for lower and upper slabs. If you are interested in specific Miller indices, please assign them by checking the box below.')
        assign_specific = st.checkbox(
            'Assign Specific Miller Indices for Lower and Upper Systems',
            value=False,
        )
        st.session_state.assign_specific = assign_specific
        if not assign_specific:
            with st.container(border=True):
                col3, col4, col5 = st.columns(3, gap='medium', border=False)
                with col3:
                    MAX_H = st.number_input('Maximum H', 0, 10, 1, step=1)
                with col4:
                    MAX_K = st.number_input('Maximum K', 0, 10, 1, step=1)
                with col5:
                    MAX_L = st.number_input('Maximum L', 0, 10, 1, step=1)
                # Ensure that at least one of the Miller indices is greater than 0
                if MAX_H == 0 and MAX_K == 0 and MAX_L == 0:
                    st.error('Please set at least one of the Miller indices greater than 0.')
                    st.stop()
                # Set session_state variables
                st.session_state.MAX_H = MAX_H
                st.session_state.MAX_K = MAX_K
                st.session_state.MAX_L = MAX_L
        else:
            st.markdown('**Upper Slab:**')
            with st.container(border=True):
                col6, col7, col8 = st.columns(3, gap='medium', border=False)
                with col6:
                    h_upper = st.number_input('H', 0, 10, 0, step=1, key='h_upper')
                with col7:
                    k_upper = st.number_input('K', 0, 10, 0, step=1, key='k_upper')
                with col8:
                    l_upper = st.number_input('L', 0, 10, 1, step=1, key='l_upper')
                # Ensure that at least one of the Miller indices is greater than 0
                if h_upper == 0 and k_upper == 0 and l_upper == 0:
                    st.error('Please set at least one of the Miller indices greater than 0.')
                    st.stop()
            st.markdown('**Lower Slab:**')
            with st.container(border=True):
                col3, col4, col5 = st.columns(3, gap='medium', border=False)
                with col3:
                    h_lower = st.number_input('H', 0, 10, 0, step=1, key='h_lower')
                with col4:
                    k_lower = st.number_input('K', 0, 10, 0, step=1, key='k_lower')
                with col5:
                    l_lower = st.number_input('L', 0, 10, 1, step=1, key='l_lower')
                # Ensure that at least one of the Miller indices is greater than 0
                if h_lower == 0 and k_lower == 0 and l_lower == 0:
                    st.error('Please set at least one of the Miller indices greater than 0.')
                    st.stop()
            # LOWER_HKL = (h_lower, k_lower, l_lower)
            # UPPER_HKL = (h_upper, k_upper, l_upper)
            LOWER_HKL = [(h_lower, k_lower, l_lower)]
            UPPER_HKL = [(h_upper, k_upper, l_upper)]
            # Set session_state variables
            st.session_state.LOWER_HKL = LOWER_HKL
            st.session_state.UPPER_HKL = UPPER_HKL

        # Add a section for Lattice Matching
        st.divider()
        st.subheader('Set Lattice Matching Parameters')
        st.caption('Set the tolerance for the misfit of lattice vectors (in %) and angle (in degree).')
        col9, col10 = st.columns(2, gap='medium', border=True)
        with col9:
            UV_TOL = st.slider(
                label='Lattice Vector Tolerance (%)',
                min_value=0,
                max_value=20,
                value=5,
                step=1,
            )
        with col10:
            ANGLE_TOL = st.slider(
                label='Angle Tolerance (°)',
                min_value=0,
                max_value=20,
                value=5,
                step=1,
            )
        # Set session_state variable
        st.session_state.UV_TOL = UV_TOL
        st.session_state.ANGLE_TOL = ANGLE_TOL

        # Add a section for Interface Geometry
        st.divider()
        st.subheader('Customize Interface Geometry')
        st.caption('Customize the interface geometry including the minimum thickness of the slab, slab vacuum, interface gap, area range for the matched interfaces, and the shape filter option.')
        MIN_THICKNESS = st.slider(
            label='Minimum Thickness of the Interface ($Å$)',
            min_value=1,
            max_value=1000,
            value=20,
            step=1,
        )
        INTERFACE_VACUUM = st.slider(
            label='Vacuum Layer Thickness ($Å$)',
            min_value=0,
            max_value=100,
            value=10,
            step=1,
        )
        INTERFACE_GAP = st.slider(
            label='Interface Gap Thickness ($Å$)',
            min_value=0,
            max_value=100,
            value=2,
            step=1,
        )
        area = st.slider(
            label='Area Range for the Matched Interfaces ($Å^2$)',
            min_value=1,
            max_value=5000,
            value=(50, 500),
            step=1,
        )
        MIN_AREA, MAX_AREA = area[0], area[1]
        SHAPE_FILTER = st.checkbox(
            'Shape Filter',
            value=True,
            help='Only keep the square-like interfaces',
        )
        # Set session_state variable
        st.session_state.MIN_THICKNESS = MIN_THICKNESS
        st.session_state.INTERFACE_VACUUM = INTERFACE_VACUUM
        st.session_state.INTERFACE_GAP = INTERFACE_GAP
        st.session_state.MIN_AREA = MIN_AREA
        st.session_state.MAX_AREA = MAX_AREA
        st.session_state.SHAPE_FILTER = SHAPE_FILTER

    # Add a button to run the interface maker
    st.divider()
    if st.button('Generate Interfaces', type='primary', width='stretch'):
        timestamp = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
        st.session_state.timestamp = timestamp
        if not st.session_state.files_uploaded:
            st.error('Please upload the conventional cell of the lower and upper systems in VASP format (.vasp/.poscar).')
            st.stop()
        with st.spinner('Generating interfaces...'):
            try:
                zip_buffer = interface_maker(st.session_state)
                st.success('Interfaces Generated Successfully!')
                st.balloons()
                # Finalize zip in memory
                zip_buffer.seek(0)
                # Streamlit download button
                st.download_button(
                    label='Download Generated Interfaces',
                    data=zip_buffer,
                    file_name=f'Interface-Maker_{timestamp}.zip',
                    mime='application/zip',
                    icon=':material/download:',
                    width='stretch',
                )
            except Exception as e:
                st.error(f'Sorry, we cannot generate the interfaces. Please check the input parameters and try again.')
                st.stop()

    # Add a footer
    st.divider()
    st.markdown(
    '''
    <div style='color: rgba(0, 0, 0, 0.4); font-weight: bold;'>
        Copyright © 2025 Interface-Maker | Developed by Guangchen Liu. All rights reserved.
    </div>
    ''',
    unsafe_allow_html=True
    )


if __name__ == '__main__':
    main()