import matplotlib.pyplot as plt
import numpy as np
import pytest
from implementation.beeler_reuter_0d import BeelerReuter0D, Stimulation


def prepare_model(model_class, dt, curr_dur, curr_value, t_prebeats):
    """
    Prepares a 2D cardiac model with a stimulation protocol.

    Parameters
    ----------
    model_class : Callable
        The cardiac model class to be instantiated.
    dt : float
        Time step for the simulation (ms or model units).
    curr_value : float
        Amplitude of the stimulus current (μA/cm² or model units).
    curr_dur : float
        Duration of each stimulus pulse (ms or model units).
    t_prebeats : float
        Interval between preconditioning stimuli (ms or model units).

    Returns
    -------
    model : Model
        Configured and initialized model ready for simulation.
    """

    stimulations = [
        Stimulation(t_start=0.0, duration=curr_dur, amplitude=curr_value),
        Stimulation(t_start=t_prebeats, duration=curr_dur, amplitude=curr_value),
        Stimulation(t_start=2 * t_prebeats, duration=curr_dur, amplitude=curr_value),
        Stimulation(t_start=3 * t_prebeats, duration=curr_dur, amplitude=curr_value),
    ]

    model = model_class(dt=dt, stimulations=stimulations)

    return model


def calculate_apd(u, dt, beat_index=3):
    """
    Calculates the action potential duration (APD) for a single beat (third by default).

    Parameters
    ----------
    u : np.ndarray
        Membrane potential time series.
    dt : float
        Time step of the simulation (ms).
    beat_index : int, optional
        Index of the beat to analyze (default is 3).

    Returns
    -------
    apd : float or None
        Duration of the action potential (ms or model units),
        or None if no complete AP was found.
    """
    threshold = u.max() - np.ptp(u) * 0.9
    print(threshold)

    up_idx = np.where((u[:-1] < threshold) & (u[1:] >= threshold))[0]
    down_idx = np.where((u[:-1] > threshold) & (u[1:] <= threshold))[0]

    if len(up_idx) <= beat_index or len(down_idx) == 0:
        return None

    ap_start = up_idx[beat_index]
    ap_end_candidates = down_idx[down_idx > ap_start]
    if len(ap_end_candidates) == 0:
        return None

    ap_end = ap_end_candidates[0]
    return (ap_end - ap_start) * dt


def test_model_attributes():
    """
    Test that the model has the expected attributes.
    Checks for the presence of key variables and parameters in the 0D Model.
    """
    model = BeelerReuter0D(dt=0.01, stimulations=[])

    assert "u" in model.variables, "Model should have variable 'u'"
    assert "x1" in model.variables, "Model should have variable 'x1'"
    assert "m" in model.variables, "Model should have variable 'm'"
    assert "h" in model.variables, "Model should have variable 'h'"
    assert "j" in model.variables, "Model should have variable 'j'"
    assert "d" in model.variables, "Model should have variable 'd'"
    assert "f" in model.variables, "Model should have variable 'f'"
    assert "cai" in model.variables, "Model should have variable 'cai'"


def test_model_step():
    """Test if all state variables are updated for each step."""

    t_prebeats = 1000.0  # interval between preconditioning stimuli (ms or model units).
    t_calc = 1000.0  # time after the last preconditioning beat to continue recording (ms or model units).
    t_max = 3 * t_prebeats + t_calc
    model = prepare_model(
        BeelerReuter0D, dt=0.01, curr_dur=0.25, curr_value=1.0, t_prebeats=t_prebeats
    )
    model.run(t_max=t_max)

    for variable in model.variables:
        values = np.array(model.history[variable])
        assert np.all(~np.isnan(values))
        assert len(values) == int(t_max / model.dt)

@pytest.mark.skip()
def test_model_run():
    """
    Test the model run for a short simulation.
    Runs the 0D Model with a predefined stimulation protocol and checks
    that the membrane potential 'u' stays within expected physiological ranges.
    """
    t_prebeats = 1000.0  # interval between preconditioning stimuli (ms or model units).
    t_calc = 1000.0  # time after the last preconditioning beat to continue recording (ms or model units).
    t_max = 3 * t_prebeats + t_calc
    model = prepare_model(
        BeelerReuter0D, dt=0.01, curr_dur=0.25, curr_value=1.0, t_prebeats=t_prebeats
    )
    model.run(t_max=t_max)
    u = np.array(model.history["u"])

    assert np.max(u) == pytest.approx(28, abs=1)
    assert np.min(u) == pytest.approx(-84, abs=1)

    apd = calculate_apd(u, model.dt, beat_index=1)
    assert apd == pytest.approx(285, abs=1)


def compare_figure():
    """
    Comparison between the model and figure data.
    The data is extracted using WebPlotDigitizer (https://apps.automeris.io/wpd4/).
    The stimulation offset is choosen in order to align with the data.
    """
    stimulations = [Stimulation(t_start=45.0, duration=0.25, amplitude=1.0)]
    t_max = 500.0

    model = BeelerReuter0D(dt=0.01, stimulations=stimulations)
    model.run(t_max=t_max)
    u = np.array(model.history["u"])

    t_ref, u_ref = np.loadtxt("tests/action_potential.csv", delimiter=",", skiprows=2).T

    plt.plot(np.arange(len(u)) * model.dt, u)
    plt.scatter(t_ref, u_ref, c="r")
    plt.show()
