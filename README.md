# Numerical Solvers: ODEs, the Heat Equation and the 2D Poisson Problem

A small Python project that implements classical numerical methods from scratch and checks them the way a numerical analyst would: against known solutions, with measured convergence orders, stability limits and solver timings.

![Heat diffusion](figures/heat_diffusion.gif)

## What is in here

| Part | Problem | Methods | What is demonstrated |
|---|---|---|---|
| 1. ODEs | $y' = f(t, y)$ | Explicit Euler, Heun (RK2), RK4, backward Euler (Newton) | Observed order of convergence, conservation properties, stiffness |
| 2. Heat equation | $u_t = \alpha\,u_{xx}$ on $(0,1)$, Dirichlet BCs | Explicit (FTCS), backward Euler, Crank-Nicolson | Order in time and space, the stability condition $r \le 1/2$ |
| 3. Poisson equation | $-\Delta u = f$ on the unit square | Jacobi (Python loops and NumPy), conjugate gradients, sparse direct | Iteration counts and run time as the grid grows |

Everything in `solvers/` is written by hand. SciPy is used only for sparse linear algebra (`splu`, `spsolve`, `cg`) and, in one experiment, `solve_ivp` as an independent reference solution.

## Methods in brief

**ODE integrators.** With step size $h$:

- Euler: $y_{n+1} = y_n + h f(t_n, y_n)$, order 1
- Heun: average of the slopes at both ends of the step, order 2
- RK4: the classical four-stage method, order 4
- Backward Euler: $y_{n+1} = y_n + h f(t_{n+1}, y_{n+1})$, solved with Newton's method, order 1 and A-stable

**Heat equation.** With mesh ratio $r = \alpha\,\Delta t / \Delta x^2$ and $A = \mathrm{tridiag}(1,-2,1)$ on the interior points:

- Explicit: $u^{n+1} = (I + rA)\,u^n$, stable only if $r \le 1/2$
- Backward Euler: $(I - rA)\,u^{n+1} = u^n$, unconditionally stable
- Crank-Nicolson: $(I - \tfrac{r}{2}A)\,u^{n+1} = (I + \tfrac{r}{2}A)\,u^n$, unconditionally stable, second order in time

All three use the second-order central difference in space.

**Poisson problem.** The 5-point stencil on an $n \times n$ grid of interior points gives a sparse linear system $Au = h^2 f$. It is solved with Jacobi iteration, conjugate gradients and a sparse direct solver. The test problem uses the manufactured solution $u(x,y) = p(x)\,p(y)$ with $p(s) = s(1-s)e^{s}$, so the exact error can be measured.

## Results

All numbers below come from `python run_all.py`.

### 1. ODE solvers

Test problem $y' = -y + \sin t$, $y(0)=1$, error measured at $t = 5$.

| Method | Theoretical order | Observed order |
|---|---|---|
| Explicit Euler | 1 | 1.03 |
| Backward Euler | 1 | 0.97 |
| Heun (RK2) | 2 | 2.02 |
| RK4 | 4 | 4.00 |

![ODE convergence](figures/ode_convergence.png)

**Nonlinear systems.** For the Lotka-Volterra model the quantity $H = dx - c\ln x + by - a\ln y$ is conserved. After $t = 40$ with $h = 0.01$ its drift is $3.7\times10^{-1}$ for Euler, $1.4\times10^{-5}$ for Heun and $6.6\times10^{-10}$ for RK4. For a damped pendulum over $t \in [0, 10]$ with $h = 0.025$, the angle error is $1.7$ rad for Euler and $2.9\times10^{-6}$ rad for RK4.

![Lotka-Volterra](figures/ode_lotka_volterra.png)

**Stiffness.** For $y' = -\lambda(y - \cos t) - \sin t$ with $\lambda = 1000$, explicit Euler is stable only for $h < 2/\lambda = 0.002$. With $h = 0.0125$ it blows up (values of order $10^{79}$), while backward Euler follows the exact solution $\cos t$.

![Stiff problem](figures/ode_stiff.png)

### 2. Heat equation

Test problem $u_0 = \sin(\pi x)$ with exact solution $e^{-\pi^2 t}\sin(\pi x)$, $T = 0.1$.

| Study | Scheme | Observed order |
|---|---|---|
| Time (dx = 1/2000) | Backward Euler | 0.98 |
| Time (dx = 1/2000) | Crank-Nicolson | 2.02 |
| Space (4000 time steps) | Crank-Nicolson | 2.00 |

![Heat convergence](figures/heat_convergence.png)

**Stability.** With a high-frequency component in the initial data, the explicit scheme stays bounded for $r = 0.4$ and $r = 0.5$ but grows by a factor of about $10^{6}$ in 60 steps for $r = 0.6$. Backward Euler and Crank-Nicolson stay bounded even for $r = 10$.

![Heat stability](figures/heat_stability.png)

### 3. Poisson solvers

Iterations and run time to reach a relative residual below $10^{-8}$. The iteration counts are deterministic; times are from one run and will differ on another machine.

| n (grid is n x n) | Jacobi iterations | CG iterations | Jacobi, Python loops | Jacobi, NumPy | CG | Sparse direct |
|---|---|---|---|---|---|---|
| 15 | 940 | 45 | 0.15 s | 0.007 s | 0.002 s | 0.002 s |
| 31 | 3,760 | 92 | 2.6 s | 0.04 s | 0.004 s | 0.004 s |
| 63 | 15,020 | 187 | not run | 0.41 s | 0.009 s | 0.012 s |
| 127 | about 60,000 | 379 | not run | 4.4 s | 0.05 s | 0.07 s |

Jacobi needs about 4 times more iterations each time the grid spacing is halved, while CG needs about 2 times more. Vectorising the Jacobi loop with NumPy is roughly 60 times faster than plain Python loops at $n = 31$. The maximum error against the exact solution falls by a factor of about 4 per refinement, as expected for a second-order scheme.

![Poisson performance](figures/poisson_performance.png)

## Getting started

Requires Python 3.9 or newer.

```bash
git clone <this-repository-url>
cd numerical-solvers
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

Run the tests:

```bash
pytest
```

Run one experiment, or all of them (this regenerates everything in `figures/`):

```bash
python -m experiments.ode_convergence
python -m experiments.ode_systems
python -m experiments.ode_stiff
python -m experiments.heat_convergence
python -m experiments.heat_stability
python -m experiments.poisson_performance
python -m experiments.make_animation
python run_all.py
```

Use the library directly:

```python
import numpy as np
from solvers import ode, heat

# ODE: y' = -y, y(0) = 1
t, y = ode.integrate("rk4", lambda t, y: -y, 1.0, 0.0, 2.0, 50)

# Heat equation with Crank-Nicolson
x, u = heat.solve_heat("crank_nicolson", lambda x: np.sin(np.pi * x), nx=100, nt=200, T=0.1)
```

## Repository layout

```
numerical-solvers/
├── solvers/
│   ├── ode.py          # Euler, Heun, RK4, backward Euler, order estimation
│   ├── heat.py         # explicit, backward Euler, Crank-Nicolson for u_t = alpha u_xx
│   └── poisson.py      # Jacobi (loops, NumPy), CG, sparse direct for -laplace(u) = f
├── experiments/        # scripts that produce the figures and tables above
├── tests/              # pytest checks of orders, stability and solver agreement
├── figures/            # generated plots and animation
├── run_all.py
└── requirements.txt
```

## Tests

The test suite checks that:

- the observed convergence orders match theory (1, 1, 2, 4 for the ODE solvers; 1 and 2 in time and 2 in space for the heat schemes; 2 for the Poisson discretisation),
- the explicit heat scheme is stable for $r \le 1/2$ and unstable for $r = 0.6$, while the implicit schemes stay bounded for $r = 10$,
- explicit Euler is unstable and backward Euler stable on the stiff problem,
- all four Poisson solvers return the same solution, and CG needs far fewer iterations than Jacobi.

## Limitations and possible extensions

- Fixed step sizes only; adaptive step-size control (for example Runge-Kutta-Fehlberg) would be a natural next step.
- The heat equation is solved in 1D with homogeneous Dirichlet boundary conditions; Neumann or time-dependent boundary conditions and a 2D version would extend it.
- The Poisson solvers do not include multigrid or preconditioned CG, which are the standard tools for large problems.
- Backward Euler uses a finite-difference Jacobian; supplying an analytic Jacobian would be faster for large systems.
- No parallel or GPU computing is used.

## References

- R. J. LeVeque, *Finite Difference Methods for Ordinary and Partial Differential Equations*, SIAM, 2007.
- E. Süli and D. Mayers, *An Introduction to Numerical Analysis*, Cambridge University Press, 2003.

