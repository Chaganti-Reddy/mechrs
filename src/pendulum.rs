//! Nonlinear pendulum solved with RK4.
//!
//! Compares exact `θ'' = -(g/L)sin(θ)` against the small-angle approximation `θ'' = -(g/L)θ`.

const G: f64 = 9.81;

/// A simple pendulum of given length.
pub struct Pendulum {
    pub length: f64,
    pub g: f64,
}

impl Pendulum {
    /// Create a pendulum with `g = 9.81 m/s²`.
    pub fn new(length: f64) -> Self {
        Self { length, g: G }
    }

    /// Full nonlinear RK4 iterator: `θ'' = -(g/L)sin(θ)`.
    /// Yields `(time, angle_rad)`.
    pub fn states(&self, theta0_rad: f64, dt: f64) -> impl Iterator<Item = (f64, f64)> {
        let gl = self.g / self.length;
        std::iter::successors(Some((0.0_f64, theta0_rad, 0.0_f64)), move |&(t, th, w)| {
            let f = |_t: f64, th: f64, w: f64| -> (f64, f64) { (w, -gl * th.sin()) };
            let (dth1, dw1) = f(t, th, w);
            let (dth2, dw2) = f(t + dt/2.0, th + dt/2.0*dth1, w + dt/2.0*dw1);
            let (dth3, dw3) = f(t + dt/2.0, th + dt/2.0*dth2, w + dt/2.0*dw2);
            let (dth4, dw4) = f(t + dt, th + dt*dth3, w + dt*dw3);
            let th_next = th + (dt/6.0)*(dth1 + 2.0*dth2 + 2.0*dth3 + dth4);
            let w_next  = w  + (dt/6.0)*(dw1  + 2.0*dw2  + 2.0*dw3  + dw4);
            Some((t + dt, th_next, w_next))
        })
        .map(|(t, th, _)| (t, th))
    }

    /// Small-angle approximation RK4 iterator: `θ'' = -(g/L)θ`.
    pub fn small_angle_states(&self, theta0_rad: f64, dt: f64) -> impl Iterator<Item = (f64, f64)> {
        let gl = self.g / self.length;
        std::iter::successors(Some((0.0_f64, theta0_rad, 0.0_f64)), move |&(t, th, w)| {
            let f = |_t: f64, th: f64, w: f64| -> (f64, f64) { (w, -gl * th) };
            let (dth1, dw1) = f(t, th, w);
            let (dth2, dw2) = f(t + dt/2.0, th + dt/2.0*dth1, w + dt/2.0*dw1);
            let (dth3, dw3) = f(t + dt/2.0, th + dt/2.0*dth2, w + dt/2.0*dw2);
            let (dth4, dw4) = f(t + dt, th + dt*dth3, w + dt*dw3);
            let th_next = th + (dt/6.0)*(dth1 + 2.0*dth2 + 2.0*dth3 + dth4);
            let w_next  = w  + (dt/6.0)*(dw1  + 2.0*dw2  + 2.0*dw3  + dw4);
            Some((t + dt, th_next, w_next))
        })
        .map(|(t, th, _)| (t, th))
    }

    /// Analytical period (small-angle only): `T = 2π√(L/g)`.
    pub fn analytical_period(&self) -> f64 {
        2.0 * std::f64::consts::PI * (self.length / self.g).sqrt()
    }

    /// Numerical period via first downward zero-crossing (time × 2).
    pub fn numerical_period(&self, theta0_rad: f64, dt: f64) -> f64 {
        let mut states = self.states(theta0_rad, dt);
        let mut prev = states.next().unwrap_or((0.0, 0.0));
        for curr in states {
            if prev.1 > 0.0 && curr.1 <= 0.0 {
                // Starting at theta0 > 0, omega = 0: first downward zero crossing is at T/4
                return prev.0 * 4.0;
            }
            prev = curr;
        }
        f64::NAN
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_angle_period_matches_analytical() {
        let p = Pendulum::new(1.0);
        let numerical = p.numerical_period(0.05_f64.to_radians(), 0.0001);
        approx::assert_abs_diff_eq!(numerical, p.analytical_period(), epsilon = 0.001);
    }

    #[test]
    fn large_angle_period_longer_than_small_angle() {
        let p = Pendulum::new(1.0);
        let large = p.numerical_period(80.0_f64.to_radians(), 0.0001);
        assert!(large > p.analytical_period(),
            "Large angle period ({:.4}) should exceed small-angle ({:.4})", large, p.analytical_period());
    }

    #[test]
    fn nonlinear_vs_small_angle_diverge_at_large_amplitude() {
        let p = Pendulum::new(1.0);
        let theta0 = 60.0_f64.to_radians();
        let dt = 0.001;
        let steps = 500;

        let nonlinear: Vec<f64> = p.states(theta0, dt).take(steps).map(|(_, th)| th).collect();
        let small_angle: Vec<f64> = p.small_angle_states(theta0, dt).take(steps).map(|(_, th)| th).collect();

        let max_diff = nonlinear.iter().zip(small_angle.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);

        assert!(max_diff > 0.01, "60° amplitude should show measurable divergence: {max_diff}");
    }
}
