//! Simple harmonic and damped harmonic oscillator via RK4.

/// A spring-mass oscillator, optionally with damping.
pub struct Oscillator {
    pub k: f64,
    pub mass: f64,
    pub damping: f64,
}

impl Oscillator {
    /// Undamped oscillator (damping = 0).
    pub fn new(k: f64, mass: f64) -> Self {
        Self { k, mass, damping: 0.0 }
    }

    /// Damped oscillator with damping coefficient `b`.
    pub fn damped(k: f64, mass: f64, b: f64) -> Self {
        Self { k, mass, damping: b }
    }

    /// Lazy iterator of `(time, displacement)` pairs.
    ///
    /// Uses RK4 on the 2D state `(x, v)` by reducing to two first-order ODEs:
    /// - `dx/dt = v`
    /// - `dv/dt = -(k/m)x - (b/m)v`
    pub fn states(&self, x0: f64, v0: f64, dt: f64) -> impl Iterator<Item = (f64, f64)> {
        let omega_sq = self.k / self.mass;
        let gamma = self.damping / self.mass;

        // Pack (x, v) into a single f64 via a 2-element array inside rk4 over x,
        // but rk4_iter only handles scalar ODEs. We implement the state machine directly.
        std::iter::successors(Some((0.0_f64, x0, v0)), move |&(t, x, v)| {
            // RK4 over 2D state (x, v)
            let f = |_t: f64, x: f64, v: f64| -> (f64, f64) {
                (v, -omega_sq * x - gamma * v)
            };
            let (dx1, dv1) = f(t, x, v);
            let (dx2, dv2) = f(t + dt / 2.0, x + dt / 2.0 * dx1, v + dt / 2.0 * dv1);
            let (dx3, dv3) = f(t + dt / 2.0, x + dt / 2.0 * dx2, v + dt / 2.0 * dv2);
            let (dx4, dv4) = f(t + dt, x + dt * dx3, v + dt * dv3);
            let x_next = x + (dt / 6.0) * (dx1 + 2.0 * dx2 + 2.0 * dx3 + dx4);
            let v_next = v + (dt / 6.0) * (dv1 + 2.0 * dv2 + 2.0 * dv3 + dv4);
            Some((t + dt, x_next, v_next))
        })
        .map(|(t, x, _v)| (t, x))
    }

    /// Analytical angular frequency `ω = √(k/m)`.
    pub fn omega(&self) -> f64 {
        (self.k / self.mass).sqrt()
    }

    /// Analytical period `T = 2π√(m/k)` (undamped only).
    pub fn analytical_period(&self) -> f64 {
        2.0 * std::f64::consts::PI / self.omega()
    }

    /// Numerically detect the period via zero-crossing of displacement.
    pub fn numerical_period(&self, x0: f64, dt: f64) -> Option<f64> {
        let mut states = self.states(x0, 0.0, dt);
        let mut prev = states.next()?;
        for curr in states {
            if prev.1 > 0.0 && curr.1 <= 0.0 {
                // Starting at x0 > 0, v = 0: first downward zero crossing is at T/4
                return Some(prev.0 * 4.0);
            }
            prev = curr;
        }
        None
    }

    /// Lazy stream of mechanical energy `E = ½kx² + ½mv²`.
    pub fn energy_stream(&self, x0: f64, v0: f64, dt: f64) -> impl Iterator<Item = f64> {
        let k = self.k;
        let m = self.mass;
        std::iter::successors(Some((0.0_f64, x0, v0)), move |&(t, x, v)| {
            let omega_sq = k / m;
            let (dx1, dv1) = (v, -omega_sq * x);
            let (dx2, dv2) = (v + dt/2.0*dv1, -omega_sq*(x + dt/2.0*dx1));
            let (dx3, dv3) = (v + dt/2.0*dv2, -omega_sq*(x + dt/2.0*dx2));
            let (dx4, dv4) = (v + dt*dv3, -omega_sq*(x + dt*dx3));
            let x_next = x + (dt/6.0)*(dx1 + 2.0*dx2 + 2.0*dx3 + dx4);
            let v_next = v + (dt/6.0)*(dv1 + 2.0*dv2 + 2.0*dv3 + dv4);
            Some((t + dt, x_next, v_next))
        })
        .map(move |(_, x, v)| 0.5 * k * x * x + 0.5 * m * v * v)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn undamped_period_matches_analytical() {
        let osc = Oscillator::new(10.0, 1.0);
        let numerical = osc.numerical_period(1.0, 0.0001).unwrap();
        let analytical = osc.analytical_period();
        approx::assert_abs_diff_eq!(numerical, analytical, epsilon = 0.001);
    }

    #[test]
    fn energy_conserved_undamped() {
        let osc = Oscillator::new(10.0, 1.0);
        let energies: Vec<f64> = osc.energy_stream(1.0, 0.0, 0.0001).take(10000).collect();
        let e0 = energies[0];
        for e in &energies {
            approx::assert_abs_diff_eq!(e, &e0, epsilon = 1e-4);
        }
    }

    #[test]
    fn damping_reduces_amplitude() {
        let osc = Oscillator::damped(10.0, 1.0, 0.5);
        let early: f64 = osc.states(1.0, 0.0, 0.001).nth(100).map(|(_, x)| x.abs()).unwrap();
        let late: f64 = osc.states(1.0, 0.0, 0.001).nth(5000).map(|(_, x)| x.abs()).unwrap();
        assert!(late < early, "Damping should reduce amplitude over time");
    }
}
