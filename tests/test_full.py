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
    assert "Returned from gampn" in point_out
    assert "Returned from asyrmo" in point_out
    assert "Returned from probamo" in point_out

    # check that the end of the script was reached (suggesting no runtime errors)
    assert "Total runtime" in point_out

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
    assert "Returned from gampn" in gamma_out
    assert "Returned from asyrmo" in gamma_out
    assert "Returned from probamo" in gamma_out

    # check that the end of the script was reached (suggesting no runtime errors)
    assert "Total runtime" in gamma_out

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
    assert "eps = [0.001, 0.301]" in eps_out
    assert "gamma = [35.0, 35.0]" in eps_out

    # check that all steps were started (suggesting no errors in prior stages)
    assert "Returned from gampn" in eps_out
    assert "Returned from asyrmo" in eps_out
    assert "Returned from probamo" in eps_out

    # check that the end of the script was reached (suggesting no runtime errors)
    assert "Total runtime" in eps_out

    # check some numerical outputs for consistency with previous runs
    assert "Energies of Spin 3/2 States / keV                           19.4  ± 1.2" in eps_out
    assert "Ground State Magnetic Dipole Moment / $μ_{N}$               2.0   ± 0.3" in eps_out

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
    assert "eps = [0.001, 0.500]" in mesh_out
    assert "gamma = [0.0, 60.0]" in mesh_out

    # check that all steps were started (suggesting no errors in prior stages)
    assert "Returned from gampn" in mesh_out
    assert "Returned from asyrmo" in mesh_out
    assert "Returned from probamo" in mesh_out

    # check that the end of the script was reached (suggesting no runtime errors)
    assert "Total runtime" in mesh_out
