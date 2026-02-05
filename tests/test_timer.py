import pytest
import time
from src.classes.timer import Timer

def test_use():

    test_timer = Timer()

    test_timer.start()
    time.sleep(1)
    test_timer.stop()

    lapsed_time = test_timer.get_lapsed_time()

    assert f"{lapsed_time:.1f}" == "1.0"

def test_misuse(capsys):

    test_timer = Timer()

    test_timer.get_lapsed_time()
    out, _ = capsys.readouterr()
    assert "WARNING: Could not get lapsed time (recorded start and/or end time are None)." in out

    test_timer.stop()
    out, _ = capsys.readouterr()
    assert "WARNING: Could not stop timer (not started yet)." in out

    test_timer.start()
    test_timer.get_lapsed_time()
    out, _ = capsys.readouterr()
    assert "WARNING: Could not get lapsed time (the timer is still running)." in out

    test_timer.start()
    out, _ = capsys.readouterr()
    assert "WARNING: Could not start timer (already running)." in out

    test_timer.stop()
    test_timer.stop()
    out, _ = capsys.readouterr()
    assert "WARNING: Could not stop timer (not started yet)." in out

    test_timer.start()
    test_timer.get_lapsed_time()
    out, _ = capsys.readouterr()
    assert "WARNING: Could not get lapsed time (the timer is still running)." in out






