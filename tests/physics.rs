//! Integration tests using known analytical solutions.

use mechrs::{
    integrators::rk4_iter,
    kinematics::Projectile,
    orbital::{self, Vec2},
    oscillator::Oscillator,
    pendulum::Pendulum,
};

// ── Integrators ──────────────────────────────────────────────────────────────

#[test]
fn rk4_matches_exact_exponential() {
    let (_, y) = rk4_iter(|_, y| y, 0.0, 1.0, 0.001).nth(1000).unwrap();
    approx::assert_abs_diff_eq!(y, 1.0_f64.exp(), epsilon = 1e-6);
}

// ── Projectile ───────────────────────────────────────────────────────────────

#[test]
fn projectile_range_45_degrees() {
    let p = Projectile::new(std::f64::consts::FRAC_PI_4, 20.0);
    approx::assert_abs_diff_eq!(p.range(0.0001), 40.775, epsilon = 0.01);
}

#[test]
fn projectile_zero_range_at_zero_angle() {
    let p = Projectile::new(0.0, 20.0);
    // Horizontal launch from y=0 immediately returns (dt = small step)
    // range is effectively 0 (first step already at y < 0 from floating point)
    assert!(p.range(0.01) < 1.0);
}

// ── Oscillator ───────────────────────────────────────────────────────────────

#[test]
fn oscillator_period_k1_m1() {
    // T = 2π ≈ 6.283
    let osc = Oscillator::new(1.0, 1.0);
    let num = osc.numerical_period(1.0, 0.0001).unwrap();
    approx::assert_abs_diff_eq!(num, std::f64::consts::TAU, epsilon = 0.001);
}

// ── Pendulum ─────────────────────────────────────────────────────────────────

#[test]
fn pendulum_small_angle_matches_analytical() {
    let p = Pendulum::new(1.0);
    let num = p.numerical_period(0.05_f64.to_radians(), 0.0001);
    approx::assert_abs_diff_eq!(num, p.analytical_period(), epsilon = 0.001);
}

#[test]
fn pendulum_large_angle_period_increases() {
    let p = Pendulum::new(1.0);
    let small = p.numerical_period(5.0_f64.to_radians(), 0.0001);
    let large = p.numerical_period(80.0_f64.to_radians(), 0.0001);
    assert!(large > small * 1.1, "80° period should be significantly longer than 5° period");
}

// ── Orbital ──────────────────────────────────────────────────────────────────

#[test]
fn circular_orbit_conserves_radius() {
    let gm = 1.327e20_f64;
    let r = 1.496e11_f64;
    let v = orbital::circular_speed(gm, r);
    let r0 = Vec2::new(r, 0.0);
    let v0 = Vec2::new(0.0, v);
    let dt = 3600.0 * 24.0;

    let max_drift = orbital::verlet_orbit(gm, r0, v0, dt)
        .take(365)
        .map(|(pos, _)| (pos.norm() / r - 1.0).abs())
        .fold(0.0_f64, f64::max);

    assert!(max_drift < 1e-3, "Radius drifts by {:.2e} (> 0.1%)", max_drift);
}
