//! Two-body Kepler orbit via Velocity Verlet integration.

use std::ops::{Add, Mul, Sub};

/// 2D vector for orbital mechanics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vec2 {
    pub x: f64,
    pub y: f64,
}

impl Vec2 {
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    pub fn zero() -> Self {
        Self { x: 0.0, y: 0.0 }
    }

    /// Euclidean norm (magnitude).
    pub fn norm(&self) -> f64 {
        (self.x * self.x + self.y * self.y).sqrt()
    }

    /// Unit vector in the same direction.
    pub fn normalized(&self) -> Self {
        let n = self.norm();
        Self {
            x: self.x / n,
            y: self.y / n,
        }
    }

    /// Dot product.
    pub fn dot(&self, other: Self) -> f64 {
        self.x * other.x + self.y * other.y
    }
}

impl Add for Vec2 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

impl Sub for Vec2 {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
        }
    }
}

impl Mul<f64> for Vec2 {
    type Output = Self;
    fn mul(self, s: f64) -> Self {
        Self {
            x: self.x * s,
            y: self.y * s,
        }
    }
}

/// Gravitational acceleration at position `r` for central mass `gm = G*M`.
fn gravity(gm: f64, r: Vec2) -> Vec2 {
    let r3 = r.norm().powi(3);
    r * (-gm / r3)
}

/// Velocity Verlet orbit iterator.
///
/// Yields `(position, velocity)` at each time step.
/// Velocity Verlet conserves energy much better than Euler for orbital dynamics.
///
/// # Arguments
/// - `gm`: G × central mass (m³/s²)
/// - `r0`: initial position (m)
/// - `v0`: initial velocity (m/s)
/// - `dt`: time step (s)
pub fn verlet_orbit(gm: f64, r0: Vec2, v0: Vec2, dt: f64) -> impl Iterator<Item = (Vec2, Vec2)> {
    // State: (position, velocity, acceleration)
    let a0 = gravity(gm, r0);
    std::iter::successors(Some((r0, v0, a0)), move |&(r, v, a)| {
        let r_next = r + v * dt + a * (0.5 * dt * dt);
        let a_next = gravity(gm, r_next);
        let v_next = v + (a + a_next) * (0.5 * dt);
        Some((r_next, v_next, a_next))
    })
    .map(|(r, v, _)| (r, v))
}

/// Specific orbital energy: `E = v²/2 - GM/r`.
pub fn specific_energy(gm: f64, r: Vec2, v: Vec2) -> f64 {
    0.5 * v.dot(v) - gm / r.norm()
}

/// Circular orbit speed at radius `r`: `v = √(GM/r)`.
pub fn circular_speed(gm: f64, r: f64) -> f64 {
    (gm / r).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    const GM_SUN: f64 = 1.327e20; // m³/s²
    const AU: f64 = 1.496e11; // m

    #[test]
    fn circular_orbit_constant_radius() {
        let r = AU;
        let v = circular_speed(GM_SUN, r);
        let r0 = Vec2::new(r, 0.0);
        let v0 = Vec2::new(0.0, v);
        let dt = 3600.0 * 24.0; // 1 day in seconds

        let radii: Vec<f64> = verlet_orbit(GM_SUN, r0, v0, dt)
            .take(365)
            .map(|(pos, _)| pos.norm())
            .collect();

        let initial = radii[0];
        for radius in &radii {
            // Allow 0.1% relative drift (Verlet is symplectic but not exact)
            approx::assert_abs_diff_eq!(radius / initial, 1.0, epsilon = 1e-3);
        }
    }

    #[test]
    fn energy_conserved_circular() {
        let r = AU;
        let v = circular_speed(GM_SUN, r);
        let r0 = Vec2::new(r, 0.0);
        let v0 = Vec2::new(0.0, v);
        let dt = 3600.0 * 24.0;
        let e0 = specific_energy(GM_SUN, r0, v0);

        let energies: Vec<f64> = verlet_orbit(GM_SUN, r0, v0, dt)
            .take(365)
            .map(|(pos, vel)| specific_energy(GM_SUN, pos, vel))
            .collect();

        for e in &energies {
            approx::assert_abs_diff_eq!(e, &e0, epsilon = e0.abs() * 1e-6);
        }
    }

    #[test]
    fn vec2_ops() {
        let a = Vec2::new(3.0, 4.0);
        assert_eq!(a.norm(), 5.0);
        let b = Vec2::new(1.0, 0.0);
        let c = Vec2::new(0.0, 1.0);
        assert_eq!((b + c), Vec2::new(1.0, 1.0));
    }
}
