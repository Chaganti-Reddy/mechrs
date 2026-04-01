//! Generic numerical integrators as lazy iterators.
//!
//! Both [`euler_iter`] and [`rk4_iter`] return infinite iterators of `(t, y)` pairs
//! using `std::iter::successors` — no Vec, no heap allocation in the hot path.

/// Euler method iterator for `dy/dt = f(t, y)`.
///
/// Simple first-order method. Fast but accumulates error quickly.
///
/// # Example
/// ```
/// use mechrs::integrators::euler_iter;
/// let (t, y) = euler_iter(|_, y| y, 0.0, 1.0, 0.001).nth(1000).unwrap();
/// // y ≈ e^1.0 ≈ 2.718 (with Euler error)
/// ```
pub fn euler_iter(
    f: impl Fn(f64, f64) -> f64 + 'static,
    t0: f64,
    y0: f64,
    dt: f64,
) -> impl Iterator<Item = (f64, f64)> {
    std::iter::successors(Some((t0, y0)), move |&(t, y)| {
        let y_next = y + dt * f(t, y);
        Some((t + dt, y_next))
    })
}

/// Fourth-order Runge-Kutta iterator for `dy/dt = f(t, y)`.
///
/// Much more accurate than Euler at the same step size.
///
/// # Example
/// ```
/// use mechrs::integrators::rk4_iter;
/// use approx::assert_abs_diff_eq;
/// let (_, y) = rk4_iter(|_, y| y, 0.0, 1.0, 0.001).nth(1000).unwrap();
/// assert_abs_diff_eq!(y, 1.0_f64.exp(), epsilon = 1e-6);
/// ```
pub fn rk4_iter(
    f: impl Fn(f64, f64) -> f64 + 'static,
    t0: f64,
    y0: f64,
    dt: f64,
) -> impl Iterator<Item = (f64, f64)> {
    std::iter::successors(Some((t0, y0)), move |&(t, y)| {
        let k1 = f(t, y);
        let k2 = f(t + dt / 2.0, y + dt / 2.0 * k1);
        let k3 = f(t + dt / 2.0, y + dt / 2.0 * k2);
        let k4 = f(t + dt, y + dt * k3);
        let y_next = y + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        Some((t + dt, y_next))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rk4_exponential_growth() {
        // dy/dt = y, y(0) = 1 => y(t) = e^t
        let (_, y) = rk4_iter(|_, y| y, 0.0, 1.0, 0.001).nth(1000).unwrap();
        approx::assert_abs_diff_eq!(y, 1.0_f64.exp(), epsilon = 1e-6);
    }

    #[test]
    fn rk4_exponential_decay() {
        // dy/dt = -y, y(0) = 1 => y(1) ≈ 0.3679
        let (_, y) = rk4_iter(|_, y| -y, 0.0, 1.0, 0.001).nth(1000).unwrap();
        approx::assert_abs_diff_eq!(y, (-1.0_f64).exp(), epsilon = 1e-6);
    }

    #[test]
    fn euler_worse_than_rk4() {
        let euler_err = {
            let (_, y) = euler_iter(|_, y| y, 0.0, 1.0, 0.01).nth(100).unwrap();
            (y - 1.0_f64.exp()).abs()
        };
        let rk4_err = {
            let (_, y) = rk4_iter(|_, y| y, 0.0, 1.0, 0.01).nth(100).unwrap();
            (y - 1.0_f64.exp()).abs()
        };
        assert!(rk4_err < euler_err, "RK4 should be more accurate than Euler");
    }
}
