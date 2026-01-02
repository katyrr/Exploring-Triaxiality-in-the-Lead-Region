import pytest
import subprocess
import sys

def test_point():

    point_res = subprocess.run(
        [sys.executable, "src/main.py", "test_point"],
        capture_output=True,
        text=True,
        check=True,
    )
    point_out = point_res.stdout
    # check the correct config file was used for the tets:
    assert not "not found" in point_out

    assert "eps = [0.260, 0.260]" in point_out
    assert "gamma = [24.0, 24.0]" in point_out

    # check that all steps were started (suggesting no errors in prior stages)
    assert "Started running gampn" in point_out
    assert "Started running asyrmo" in point_out
    assert "Started running probamo" in point_out

    # check that the end of the script was reached (suggesting no runtime errors)
    assert "total runtime" in point_out

def test_gamma_range():

    gamma_res = subprocess.run(
        [sys.executable, "src/main.py", "test_gamma_range"],
        capture_output=True,
        text=True,
        check=True,
    )
    gamma_out = gamma_res.stdout

    # check the correct config file was used for the tets:
    assert not "not found" in gamma_out
    assert "eps = [0.280, 0.280]" in gamma_out
    assert "gamma = [16.0, 34.0]" in gamma_out

    # check that all steps were started (suggesting no errors in prior stages)
    assert "Started running gampn" in gamma_out
    assert "Started running asyrmo" in gamma_out
    assert "Started running probamo" in gamma_out

    # check that the end of the script was reached (suggesting no runtime errors)
    assert "total runtime" in gamma_out

def test_eps_range():

    eps_res = subprocess.run(
        [sys.executable, "src/main.py", "test_eps_range"],
        capture_output=True,
        text=True,
        check=True,
    )
    eps_out = eps_res.stdout

    # check the correct config file was used for the tets:
    assert not "not found" in eps_out
    assert "eps = [0.250, 0.350]" in eps_out
    assert "gamma = [22.0, 22.0]" in eps_out

    # check that all steps were started (suggesting no errors in prior stages)
    assert "Started running gampn" in eps_out
    assert "Started running asyrmo" in eps_out
    assert "Started running probamo" in eps_out

    # check that the end of the script was reached (suggesting no runtime errors)
    assert "total runtime" in eps_out

def test_mesh():

    mesh_res = subprocess.run(
        [sys.executable, "src/main.py", "test_mesh"],
        capture_output=True,
        text=True,
        check=True,
    )
    mesh_out = mesh_res.stdout

    # check the correct config file was used for the tets:
    assert not "not found" in mesh_out
    assert "eps = [0.280, 0.280]" in mesh_out
    assert "gamma = [16.0, 34.0]" in mesh_out

    # check that all steps were started (suggesting no errors in prior stages)
    assert "Started running gampn" in mesh_out
    assert "Started running asyrmo" in mesh_out
    assert "Started running probamo" in mesh_out

    # check that the end of the script was reached (suggesting no runtime errors)
    assert "total runtime" in mesh_out


'''
DEBEUG 44: {
    'spin_1/2_energies': [507.2, 625.3], 'spin_1/2_mag_moments': [0.33, 0.7171], 'spin_1/2_quad_moments': [0.0, 0.0], 
    'spin_3/2_energies': [379.1, 596.4, 764.2], 'spin_3/2_mag_moments': [-0.1609, 0.2866, 0.2933], 'spin_3/2_quad_moments': [-1.825, 1.9423, 1.8873], 
    'spin_5/2_energies': [0.0, 410.4, 716.3, 862.2], 'spin_5/2_mag_moments': [0.5394, 0.7278, 0.9458, 1.243], 'spin_5/2_quad_moments': [-3.7555, 1.4857, 2.2263, 2.8666], 
    'gs_spin_strings': '5/2', 'gs_spin_floats': 2.5, 'gs_mag_moments': 0.5394, 'gs_quad_moments': -3.7555, 
    'spin_7/2_energies': [27.3, 165.6, 599.1, 827.4], 'spin_7/2_mag_moments': [0.3687, 0.5826, 0.5382, 1.0668], 'spin_7/2_quad_moments': [-4.6146, -0.933, 2.2884, 3.2978], 
    'spin_9/2_energies': [237.5, 381.7, 657.2], 'spin_9/2_mag_moments': [1.2364, 0.721, 1.7448], 'spin_9/2_quad_moments': [-1.3429, 0.5785, 3.2512], 
    'spin_11/2_energies': [495.0, 639.3], 'spin_11/2_mag_moments': [1.9462, 0.8239], 'spin_11/2_quad_moments': [0.563, 1.6165], 
    'spin_13/2_energies': [786.4, 962.9], 'spin_13/2_mag_moments': [2.5636, 1.1852], 'spin_13/2_quad_moments': [1.904, 2.2354]}

DEBEUG 45: {
    'spin_1/2_energies': [0.0], 'spin_1/2_mag_moments': [0.4727], 'spin_1/2_quad_moments': [0.0], 
    'gs_spin_strings': '1/2', 'gs_spin_floats': 0.5, 'gs_mag_moments': 0.4727, 'gs_quad_moments': 0.0, 
    'spin_3/2_energies': [138.2, 612.6, 735.1], 'spin_3/2_mag_moments': [0.5932, -0.2367, -1.6387], 'spin_3/2_quad_moments': [-3.9829, 3.9938, 0.8415], 
    'spin_5/2_energies': [157.5, 284.2, 744.0, 858.3], 'spin_5/2_mag_moments': [0.9969, 0.9464, 0.411, -0.4545], 'spin_5/2_quad_moments': [-5.6316, 6.9897, -1.2774, -1.6312], 
    'spin_7/2_energies': [428.8, 483.5, 573.1, 934.9], 'spin_7/2_mag_moments': [1.3149, 1.0853, -1.3768, 1.1193], 'spin_7/2_quad_moments': [0.3011, -5.6763, -5.2585, -3.8344], 
    'spin_9/2_energies': [510.2, 571.6, 644.7], 'spin_9/2_mag_moments': [1.6705, -0.7851, 1.6713], 'spin_9/2_quad_moments': [-6.8406, 10.5839, -2.0198], 
    'spin_11/2_energies': [615.9, 799.2, 872.1], 'spin_11/2_mag_moments': [-0.6962, -0.1659, 2.0344], 'spin_11/2_quad_moments': [-7.0042, 5.0622, -3.4351], 
    'spin_13/2_energies': [nan], 'spin_13/2_mag_moments': [nan]}
'''