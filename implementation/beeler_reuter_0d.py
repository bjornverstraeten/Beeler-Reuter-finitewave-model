"""
This module provides a simple interface to run the model in a 0D setting,
i.e., without spatial dimensions. It includes class for defining stimulation protocols
and a class for the 0D model itself.

"""

from beeler_reuter import ops


class Stimulation:
    """
    Stimulus protocol for the 0D model.

    Parameters
    ----------
    t_start : float
        Start time (ms) of the first stimulus window.
    duration : float
        Duration (ms) of a single pulse.
    amplitude : float
        Pulse amplitude in the same units as du/dt contribution (typically "units/ms").

    Method
    ------
    stim(t: float) -> float
        Returns the instantaneous stimulus value at time t.

    """

    def __init__(self, t_start: float, duration: float, amplitude: float):
        self.t_start = t_start
        self.duration = duration
        self.amplitude = amplitude

    def stim(self, t: float) -> float:
        return (
            self.amplitude if self.t_start <= t < self.t_start + self.duration else 0.0
        )


class BeelerReuter0D:
    """
    Beeler-Reuter OD implementation.

    Parameters
    ----------

    dt : float
        Time step size (ms).
    stimulations : list[Stimulation]
        List of stimulation protocols to apply during the simulation.

    Attributes
    ----------
    variables : dict[str, float]
        Current state variables of the model.
    parameters : dict[str, float]
        Model parameters.
    history : dict[str, list[float]]
        Time history of state variables for post-processing.

    Methods
    -------
    step(i: int)
        Perform a single time step update.
    run(t_max: float)
        Run the simulation up to time t_max.
    """

    def __init__(self, dt: float, stimulations: list[Stimulation]):
        self.dt = dt
        self.stimulations = stimulations
        self.variables = ops.get_variables()
        self.parameters = ops.get_parameters()
        self.history: dict[str, list] = {s: [] for s in self.variables}

    def step(self, i: int):
        """
        Perform a single time step update.

        Parameters
        ----------
        i : int
            Current time step index.
        """

        du = self.dt * ops.calc_rhs(
            self.variables["u"],
            self.variables["x1"],
            self.variables["m"],
            self.variables["h"],
            self.variables["j"],
            self.variables["d"],
            self.variables["f"],
            self.variables["cai"],
            self.parameters["gna"],
            self.parameters["gnac"],
            self.parameters["gs"],
            self.parameters["E_Na"],
            self.parameters["cm"],
        )
        dcai = self.dt * ops.calc_dcai(
            self.variables["u"],
            self.variables["d"],
            self.variables["f"],
            self.variables["cai"],
            self.parameters["gs"],
        )
        dx1 = self.dt * ops.calc_dx1(
            self.variables["u"],
            self.variables["x1"],
            self.parameters["a_x1"],
            self.parameters["b_x1"],
        )
        dm = self.dt * ops.calc_dm(
            self.variables["u"],
            self.variables["m"],
            self.parameters["a_m"],
            self.parameters["b_m"],
        )
        dh = self.dt * ops.calc_dh(
            self.variables["u"],
            self.variables["h"],
            self.parameters["a_h"],
            self.parameters["b_h"],
        )
        dj = self.dt * ops.calc_dj(
            self.variables["u"],
            self.variables["j"],
            self.parameters["a_j"],
            self.parameters["b_j"],
        )
        dd = self.dt * ops.calc_dd(
            self.variables["u"],
            self.variables["d"],
            self.parameters["a_d"],
            self.parameters["b_d"],
        )
        df = self.dt * ops.calc_df(
            self.variables["u"],
            self.variables["f"],
            self.parameters["a_f"],
            self.parameters["b_f"],
        )

        stim = sum(stim.stim(t=self.dt * i) for stim in self.stimulations)

        self.variables["u"] += du + stim
        self.variables["cai"] += dcai
        self.variables["x1"] += dx1
        self.variables["m"] += dm
        self.variables["h"] += dh
        self.variables["j"] += dj
        self.variables["d"] += dd
        self.variables["f"] += df

    def run(self, t_max: float):
        """
        Run the simulation up to time t_max.

        Parameters
        ----------
        t_max : float
            Maximum simulation time.
        """
        n_steps = int(round(t_max / self.dt))
        for i in range(n_steps):
            self.step(i)
            for s in self.variables:
                self.history[s].append(self.variables[s])
