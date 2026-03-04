## Finitewave model template (replace with the model name)

The Beeler-Reuter model is a cardiac ionic model of mammalian ventricular myocardial fibres using the Hodgkin-Huxley formalism.
The model is based on ionic currents that have been measured by the voltage-clamp method.
The ionic current contain four components:
A fast inward sodium current,
a slow inward current,
and a time-dependent and independent outward potassium current.
Additionally, the intracellular calcium concentration is included as a dynamical variable.

This model implementation can be used separately from the Finitewave, allowing for standalone simulations and testing of the model dynamics without the need for the entire framework.

### Reference

Beeler GW, Reuter H.
Reconstruction of the action potential of ventricular myocardial fibres.
J Physiol. 1977 Jun;268(1):177-210.

DOI: <https://doi.org/10.1113/jphysiol.1977.sp011853>

### How to use

```bash
python -m examples.beeler_reuter_example
```

### How to test

```bash
python -m pytest -q
```

### Repository structure

```text
.
├── beeler_reuter/                   # equations package (ops.py)
│   ├── __init__.py
│   └── ops.py                       # model equations (pure functions)
├── implementation/                  # 0D model implementation
│   ├── __init__.py
│   └── beeler_reuter_0d.py
├── example/
│   └── beeler_reuter_example.py     # minimal script to run a short trace
├── tests/
│   └── test.py                      # smoke test; reproducibility checks
├── .gitignore
├── LICENSE                          # MIT
├── pyproject.toml                   # configuration file
└── README.md                        # this file
```

### Variables

Model state variables: description, units and ranges (optional)

- `u = -84` - Membrane potential (mV)
- `cai = 1e-7` - Intracellular calcium concentration (mole)
- `x1 = 0.00595` - Gating variable for the time-dependent outward current $I_{x1}$
- `m = 0.0118` - Activation variable for the fast inward sodium current $I_{Na}$
- `h = 0.985` - Inactivation variable for the fast inward sodium current $I_{Na}$
- `j = 0.972` - Slow inactivation variable for the fast inward sodium current $I_{Na}$
- `d = 0.000527` - Activation gate for the slow inward current $I_s$
- `f = 0.997` - Inactivation gate for the slow inward current $I_s$

### Parameters

Parameters and their default values

Ion Channel Conductances ($mmho/cm2$):

- `gna = 4.0` - Fast sodium ($Na^+$) conductance
- `gnac = 0.003` - Fast sodium ($Na^+$) conductance
- `gs = 0.09` - Slow inward calcium ($Ca^{2+}$) conductance

Membrane Capacity:

- cm = 1.0 ($\mu F/cm^2$)

Equilibrium potential (mV):

- $E_{Na}$ = 50.0

Rate constants:

| Rate (ms⁻¹)   | C₁ (ms⁻¹) | C₂ (mV⁻¹) | C₃ (mV) | C₄ ((mV·ms)⁻¹) | C₅ (mV) | C₆ (mV⁻¹) | C₇  |
| ------------- | --------- | --------- | ------- | -------------- | ------- | --------- | --- |
| $\alpha_{x1}$ | 0.0005    | 0.083     | 50      | 0              | 0       | 0.057     | 1   |
| $\beta_{x1}$  | 0.0013    | -0.06     | 20      | 0              | 0       | -0.04     | 1   |
| $\alpha_m$    | 0         | 0         | 47      | -1             | 47      | -0.1      | -1  |
| $\beta_m$     | 40        | -0.056    | 72      | 0              | 0       | 0         | 0   |
| $\alpha_h$    | 0.126     | -0.25     | 77      | 0              | 0       | 0         | 0   |
| $\beta_h$     | 1.7       | 0         | 22.5    | 0              | 0       | -0.082    | 1   |
| $\alpha_j$    | 0.055     | -0.25     | 78      | 0              | 0       | -0.2      | 1   |
| $\beta_j$     | 0.3       | 0         | 32      | 0              | 0       | -0.1      | 1   |
| $\alpha_d$    | 0.095     | 0.01      | -5      | 0              | 0       | -0.072    | 1   |
| $\beta_d$     | 0.07      | -0.017    | 44      | 0              | 0       | 0.05      | 1   |
| $\alpha_f$    | 0.012     | 0.08      | 28      | 0              | 0       | 0.15      | 1   |
| $\beta_f$     | 0.0065    | -0.02     | 30      | 0              | 0       | -0.2      | 1   |
