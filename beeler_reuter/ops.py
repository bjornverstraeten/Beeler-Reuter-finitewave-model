"""
ops.py — mathematical core of the model.

This module provides functions to compute the model equations,
as well as functions to retrieve default parameters and initial
values for the state variables.

The Beeler-Reuter model is a cardiac ionic model of mammalian ventricular myocardial fibres using the Hodgkin-Huxley formalism.
The model is based on ionic currents that have been measured by the voltage-clamp method.
The ionic current contain four components:

References:
- Beeler GW, Reuter H. Reconstruction of the action potential of ventricular myocardial fibres.
   J Physiol. 1977 Jun;268(1):177-210.

DOI: <https://doi.org/10.1113/jphysiol.1977.sp011853>
"""

__all__ = (
    "get_variables",
    "get_parameters",
    "calc_rhs",
    "calc_dx1",
    "calc_dm",
    "calc_dh",
    "calc_dj",
    "calc_df",
)

import math


def get_variables() -> dict[str, float]:
    """
    Returns default initial values for state variables.

    State Variables:
    u   =  -84       # Membrane potential (mV)
    cai =  1e-7      # Intracellular calcium concentration (mole)
    x1  =  0.00595   # Gating variable for the time-dependent outward current ix1
    m   =  0.0118    # Activation variable for the fast inward sodium current ina
    h   =  0.985     # Inactivation variable for the fast inward sodium current ina
    j   =  0.972     # Slow inactivation variable for the fast inward sodium current ina
    d   =  0.000527  # Activation gate for the slow inward current isi
    f   =  0.997     # Inactivation gate for the slow inward current isi
    """
    return {
        "u": -84,
        "cai": 1e-7,
        "x1": 0.00595,
        "m": 0.0118,
        "h": 0.985,
        "j": 0.972,
        "d": 0.000527,
        "f": 0.997,
    }


def get_parameters() -> dict[str, float | list[float]]:
    """
    Returns default parameter values for the model.

    Ion Channel Conductances (mmho/cm2):
    gna = 4.0     # Fast sodium (Na+) conductance
    gnac = 0.003  # Fast sodium (Na+) conductance
    gs = 0.09    # Slow inward calcium (Ca2+) conductance

    Membrane Capacity:
    cm = 1.0 (uF/cm2)

    Equilibrium potential (mV):
    E_Na = 50.0

    Rate constants:
    | Rate (ms⁻¹) | C₁ (ms⁻¹) | C₂ (mV⁻¹) | C₃ (mV) | C₄ ((mV·ms)⁻¹) | C₅ (mV) | C₆ (mV⁻¹) | C₇ |
    |-------------|-----------|-----------|---------|----------------|---------|-----------|----|
    | αₓ₁         | 0.0005    | 0.083     | 50      | 0              | 0       | 0.057     | 1  |
    | βₓ₁         | 0.0013    | -0.06     | 20      | 0              | 0       | -0.04     | 1  |
    | αₘ          | 0         | 0         | 47      | -1             | 47      | -0.1      | -1 |
    | βₘ          | 40        | -0.056    | 72      | 0              | 0       | 0         | 0  |
    | αₕ          | 0.126     | -0.25     | 77      | 0              | 0       | 0         | 0  |
    | βₕ          | 1.7       | 0         | 22.5    | 0              | 0       | -0.082    | 1  |
    | αⱼ          | 0.055     | -0.25     | 78      | 0              | 0       | -0.2      | 1  |
    | βⱼ          | 0.3       | 0         | 32      | 0              | 0       | -0.1      | 1  |
    | α_d         | 0.095     | 0.01      | -5      | 0              | 0       | -0.072    | 1  |
    | β_d         | 0.07      | -0.017    | 44      | 0              | 0       | 0.05      | 1  |
    | α_f         | 0.012     | 0.08      | 28      | 0              | 0       | 0.15      | 1  |
    | β_f         | 0.0065    | -0.02     | 30      | 0              | 0       | -0.2      | 1  |
    """

    return {
        "gna": 4.0,
        "gnac": 0.003,
        "gs": 0.09,
        "cm": 1.0,
        "E_Na": 50.0,
        "a_x1": [0.0005, 0.083, 50, 0, 0, 0.057, 1],
        "b_x1": [0.0013, -0.06, 20, 0, 0, -0.04, 1],
        "a_m": [0, 0, 47, -1, 47, -0.1, -1],
        "b_m": [40, -0.056, 72, 0, 0, 0, 0],
        "a_h": [0.126, -0.25, 77, 0, 0, 0, 0],
        "b_h": [1.7, 0, 22.5, 0, 0, -0.082, 1],
        "a_j": [0.055, -0.25, 78, 0, 0, -0.2, 1],
        "b_j": [0.3, 0, 32, 0, 0, -0.1, 1],
        "a_d": [0.095, 0.01, -5, 0, 0, -0.072, 1],
        "b_d": [0.07, -0.017, 44, 0, 0, 0.05, 1],
        "a_f": [0.012, 0.08, 28, 0, 0, 0.15, 1],
        "b_f": [0.0065, -0.02, 30, 0, 0, -0.2, 1],
    }


def calc_rhs(u, x1, m, h, j, d, f, cai, gna, gnac, gs, E_Na, cm, exp=math.exp) -> float:
    """
    Computes the right-hand side of the Beeler-Reuter model for the transmembrane potential u.
    This function implements the ordinary differential equation governing the
    evolution of the transmembrane potential `u`, which represents the electrical
    activity of cardiac cells. The equation includes contributions from various ionic
    currents that influence the membrane potential.

    Parameters
    ----------
    u : float
        Membrane potential [mV].
    x1 : float
        Gating variable for the time-dependent outward current ix1.
    m : float
        Activation variable for the fast inward sodium current ina.
    h : float
        Inactivation variable for the fast inward sodium current ina.
    j : float
        Slow inactivation variable for the fast inward sodium current ina.
    d : float
        Activation gate for the slow inward current isi.
    f : float
        Inactivation gate for the slow inward current isi.
    cai : float
        Intracellular calcium concentration [mole].
    gna : float
        Fast sodium (Na+) conductance [mmho/cm2].
    gnac : float
        Fast sodium (Na+) conductance [mmho/cm2].
    gs : float
        Slow inward calcium (Ca2+) conductance [mmho/cm2].
    E_Na : float
        Equilibrium potential for sodium [mV].
    cm : float
        Membrane capacity [uF/cm2].
    exp : callable, optional
        Exponential function to use (default: math.exp).

    Returns
    -------
    float
        The computed right-hand side of the differential equation for `u`.
    """
    ik1 = calc_ik1(u, exp=exp)
    ix1 = calc_ix1(u, x1, exp=exp)
    ina = calc_ina(u, m, h, j, gna, gnac, E_Na)
    isi = calc_isi(u, d, f, cai, gs)
    return -ik1 - ix1 - ina - isi


def calc_ik1(u, exp=math.exp) -> float:
    """
    Computes the time-independent inward rectifier potassium current (ik1).
    This current is responsible for maintaining the resting membrane potential.

    Parameters
    ----------
    u : float
        Membrane potential [mV].
    exp : callable, optional
        Exponential function to use (default: math.exp).

    Returns
    -------
    float
        The computed ik1 current.
    """
    a1 = 4 * (exp(0.04 * (u + 85)) - 1) / (exp(0.08 * (u + 53)) + exp(0.04 * (u + 53)))
    if u == -23:
        return 0.35 * (a1 + 1)

    return 0.35 * (a1 + 0.2 * (u + 23) / (1 - exp(-0.04 * (u + 23))))


def calc_ix1(u, x1, exp=math.exp) -> float:
    """
    Computes the time-dependent outward potassium current (ix1).
    This current is activated by depolarization and contributes to repolarization.

    Parameters
    ----------
    u : float
        Membrane potential [mV].
    x1 : float
        Gating variable for the time-dependent outward current ix1.
    exp : callable, optional
        Exponential function to use (default: math.exp).

    Returns
    -------
    float
        The computed ix1 current.
    """
    return x1 * (0.8 * (exp(0.04 * (u + 77)) - 1) / exp(0.04 * (u + 35)))


def calc_ina(u, m, h, j, gna, gnac, E_Na) -> float:
    """
    Computes the fast inward sodium current (ina).
    This current is responsible for the rapid upstroke of the action potential.

    Parameters
    ----------
    u : float
        Membrane potential [mV].
    m : float
        Activation variable for the fast inward sodium current ina.
    h : float
        Inactivation variable for the fast inward sodium current ina.
    j : float
        Slow inactivation variable for the fast inward sodium current ina.
    gna : float
        Fast sodium (Na+) conductance [mmho/cm2].
    gnac : float
        Fast sodium (Na+) conductance [mmho/cm2].
    E_Na : float
        Equilibrium potential for sodium [mV].

    Returns
    -------
    float
        The computed ina current.
    """
    return (gna * m**3 * h * j + gnac) * (u - E_Na)


def calc_isi(u, d, f, cai, gs, ln=math.log) -> float:
    """
    Computes the slow inward calcium current (isi).
    This current is responsible for the plateau phase of the action potential.

    Parameters
    ----------
    u : float
        Membrane potential [mV].
    d : float
        Activation gate for the slow inward current isi.
    f : float
        Inactivation gate for the slow inward current isi.
    cai : float
        Intracellular calcium concentration [mole].
    gs : float
        Slow inward calcium (Ca2+) conductance [mmho/cm2].
    ln : callable, optional
        Natural logarithm function to use (default: math.log).

    Returns
    -------
    float
        The computed isi current.
    """
    E_s = -82.3 - 13.0287 * ln(cai)
    return gs * d * f * (u - E_s)


def calc_dcai(u, d, f, cai, gs) -> float:
    """
    Computes the rate of change of intracellular calcium concentration (dcai/dt).
    This function models the dynamics of intracellular calcium concentration.

    Parameters
    ----------
    cai : float
        Intracellular calcium concentration [mole].
    isi : float
        Slow inward calcium current [μA/μF].

    Returns
    -------
    float
        The computed rate of change of intracellular calcium concentration.
    """
    isi = calc_isi(u, d, f, cai, gs)
    return -(10 ** (-7)) * isi + 0.07 * (10 ** (-7) - cai)


def calc_dy(u, y, a_co, b_co) -> float:
    """
    Computes the time derivative of a generic gating variable (dy/dt).
    This function is used to update the state of gating variables in the model.

    Parameters
    ----------
    u : float
        Membrane potential [mV].
    y : float
        Current value of the gating variable.
    a_co : list[float]
        Parameters for the alpha rate function.
    b_co : list[float]
        Parameters for the beta rate function.

    Returns
    -------
    float
        The computed time derivative of the gating variable.
    """
    a = calc_rate(u, a_co)
    b = calc_rate(u, b_co)
    ty = 1 / (a + b)
    y_inf = a / (a + b)
    return (y_inf - y) / ty


def calc_rate(u, consts, exp=math.exp) -> float:
    """
    Computes the rate constant for a gating variable.
    This function calculates the voltage-dependent rate constants for gating variables.

    Parameters
    ----------
    u : float
        Membrane potential [mV].
    consts : list[float]
        Parameters for the rate function.
    exp : callable, optional
        Exponential function to use (default: math.exp).

    Returns
    -------
    float
        The computed rate constant.
    """
    return (
        consts[0] * exp(consts[1] * (u + consts[2])) + consts[3] * (u + consts[4])
    ) / (exp(consts[5] * (u + consts[2])) + consts[6])


def calc_dx1(u, x1, a_x1, b_x1) -> float:
    """
    Computes the time derivative of the gating variable x1 (dx1/dt).
    This function updates the state of the gating variable x1.

    Parameters
    ----------
    u : float
        Membrane potential [mV].
    x1 : float
        Current value of the gating variable x1.
    a_x1 : list[float]
        Parameters for the alpha rate function for x1.
    b_x1 : list[float]
        Parameters for the beta rate function for x1.

    Returns
    -------
    float
        The computed time derivative of the gating variable x1.
    """
    return calc_dy(u, x1, a_x1, b_x1)


def calc_dm(u, m, a_m, b_m) -> float:
    """
    Computes the time derivative of the gating variable m (dm/dt).
    This function updates the state of the gating variable m.

    Parameters
    ----------
    u : float
        Membrane potential [mV].
    m : float
        Current value of the gating variable m.
    a_m : list[float]
        Parameters for the alpha rate function for m.
    b_m : list[float]
        Parameters for the beta rate function for m.

    Returns
    -------
    float
        The computed time derivative of the gating variable m.
    """
    return calc_dy(u, m, a_m, b_m)


def calc_dh(u, h, a_h, b_h) -> float:
    """
    Computes the time derivative of the gating variable h (dh/dt).
    This function updates the state of the gating variable h.

    Parameters
    ----------
    u : float
        Membrane potential [mV].
    h : float
        Current value of the gating variable h.
    a_h : list[float]
        Parameters for the alpha rate function for h.
    b_h : list[float]
        Parameters for the beta rate function for h.

    Returns
    -------
    float
        The computed time derivative of the gating variable h.
    """
    return calc_dy(u, h, a_h, b_h)


def calc_dj(u, j, a_j, b_j) -> float:
    """
    Computes the time derivative of the gating variable j (dj/dt).
    This function updates the state of the gating variable j.

    Parameters
    ----------
    u : float
        Membrane potential [mV].
    j : float
        Current value of the gating variable j.
    a_j : list[float]
        Parameters for the alpha rate function for j.
    b_j : list[float]
        Parameters for the beta rate function for j.

    Returns
    -------
    float
        The computed time derivative of the gating variable j.
    """
    return calc_dy(u, j, a_j, b_j)


def calc_dd(u, d, a_d, b_d) -> float:
    """
    Computes the time derivative of the gating variable d (df/dt).
    This function updates the state of the gating variable d.

    Parameters
    ----------
    u : float
        Membrane potential [mV].
    d : float
        Current value of the gating variable d.
    a_d : list[float]
        Parameters for the alpha rate function for d.
    b_d : list[float]
        Parameters for the beta rate function for d.

    Returns
    -------
    float
        The computed time derivative of the gating variable d.
    """
    return calc_dy(u, d, a_d, b_d)


def calc_df(u, f, a_f, b_f) -> float:
    """
    Computes the time derivative of the gating variable f (df/dt).
    This function updates the state of the gating variable f.

    Parameters
    ----------
    u : float
        Membrane potential [mV].
    f : float
        Current value of the gating variable f.
    a_f : list[float]
        Parameters for the alpha rate function for f.
    b_f : list[float]
        Parameters for the beta rate function for f.

    Returns
    -------
    float
        The computed time derivative of the gating variable f.
    """
    return calc_dy(u, f, a_f, b_f)
