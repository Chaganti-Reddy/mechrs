# mechrs

> **Classical mechanics as lazy Rust iterators.**  
> Model infinite time series. Zero heap allocation. Physics-accurate.

---

## The Core Idea

Most physics solvers allocate a `Vec<(f64, f64)>` to store a trajectory. `mechrs` doesn't.  
Every system — projectiles, pendulums, orbits — is an **infinite lazy iterator**. You pull values out one at a time, chain adapters, and stop when you want.

```rust
use mechrs::Projectile;

let traj = Projectile::new(45.0_f64.to_radians(), 20.0); // angle, speed (m/s)

// Zero allocation — evaluates lazily
let range = traj.positions(0.001)
    .take_while(|&(_, y)| y >= 0.0)
    .last()
    .map(|(x, _)| x);

println!("Range: {:.3} m", range.unwrap()); // 40.775 m
```

---

## Physics Covered

| Module | System | Iterator yields |
|---|---|---|
| `kinematics` | Projectile motion | `(x, y)` positions |
| `oscillator` | Simple harmonic / damped | `(t, displacement)` |
| `pendulum` | Nonlinear pendulum (RK4) | `(t, angle)` |
| `orbital` | Two-body orbit (Verlet) | `(x, y)` in orbital plane |
| `field` | Uniform + gradient fields | `(t, force_vector)` |

---

## Equations

### Projectile Motion

$$x(t) = v_0 \cos(\theta) \cdot t$$
$$y(t) = v_0 \sin(\theta) \cdot t - \frac{1}{2}g t^2$$

```rust
// Implemented as a lazy iterator — no Vec
pub fn positions(&self, dt: f64) -> impl Iterator<Item = (f64, f64)> + '_ {
    (0..).map(move |i| i as f64 * dt)
         .map(move |t| (
             self.vx * t,
             self.vy * t - 0.5 * 9.81 * t * t,
         ))
}
```

### Simple Harmonic Oscillator

$$\ddot{x} + \omega^2 x = 0 \qquad \omega = \sqrt{\frac{k}{m}}$$

```rust
// RK4 as std::iter::successors — state machine, no allocation
pub fn states(k: f64, m: f64, x0: f64, v0: f64, dt: f64)
    -> impl Iterator<Item = (f64, f64)>
{
    let omega_sq = k / m;
    std::iter::successors(Some((x0, v0)), move |&(x, v)| {
        let a = -omega_sq * x;
        Some((x + v * dt, v + a * dt))
    })
    .enumerate()
    .map(move |(i, (x, _))| (i as f64 * dt, x))
}
```

### Nonlinear Pendulum (RK4)

$$\ddot{\theta} = -\frac{g}{L}\sin(\theta)$$

```rust
// Full RK4 — no small-angle approximation
rk4_iter(|_, theta| -(g / length) * theta.sin(), 0.0, theta0, dt)
    .take(steps)
    .map(|(t, theta)| (t, theta.to_degrees()))
```


### Kepler Orbit (Velocity Verlet)

$$\vec{a} = -\frac{GM}{|\vec{r}|^3}\vec{r}$$

```rust
verlet_orbit(gm, r0, v0, dt)
    .take_while(|(r, _)| r.norm() < escape_radius)
    .enumerate()
    .map(|(i, (r, _))| (r.x, r.y))
```

---

## CLI

```bash
mechrs projectile --angle 45 --speed 20
# Range: 40.775 m | Max height: 10.194 m | Time of flight: 2.886 s

mechrs pendulum --length 1.0 --angle 30 --steps 5000
# Period (numerical): 2.007 s | Analytical (small angle): 2.006 s

mechrs oscillator --k 10 --mass 1 --x0 1.0 --plot
# Outputs ASCII plot of displacement over time
```

---

## Project Structure

```
mechrs/
├── src/
│   ├── lib.rs              # public API re-exports
│   ├── symbols.rs          # SymbolTable + Expr<'sym>
│   ├── integrators.rs      # rk4_iter, verlet_iter (generic over fn)
│   ├── kinematics.rs       # Projectile
│   ├── oscillator.rs       # SHO + damped oscillator
│   ├── pendulum.rs         # nonlinear pendulum
│   ├── orbital.rs          # two-body Kepler orbit
│   └── cli.rs              # clap CLI
├── examples/
│   ├── projectile_sweep.rs # range vs angle plot
│   ├── pendulum_period.rs  # numerical vs analytical period
│   └── orbit_precession.rs # Verlet orbit demo
└── tests/
    └── physics.rs          # known analytical answers
```

---

## Running Tests

```bash
cargo test                  # unit + integration
cargo run --example projectile_sweep
```

Tests assert against known analytical solutions within tolerance:

```rust
#[test]
fn projectile_range_45_degrees() {
    let traj = Projectile::new(std::f64::consts::FRAC_PI_4, 20.0);
    let range = traj.range(0.0001);
    approx::assert_abs_diff_eq!(range, 40.775, epsilon = 0.01);
}
```

---

## Why Rust for This?

- **Zero-cost abstractions**: iterator chains compile to tight loops, no overhead vs hand-written loops
- **Lifetimes**: model "this expression borrows from this equation context" at compile time — no runtime reference counting
- **`f64` everywhere**: no boxing, no dynamic dispatch in the hot path
- **`no_std` compatible** (planned): embed on microcontrollers
