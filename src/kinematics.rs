//! Projectile motion under constant gravity.
//!
//! All methods return lazy iterators — no allocation, no precomputed Vec.

const G: f64 = 9.81;

/// A projectile launched at a given angle and speed.
pub struct Projectile {
    vx: f64,
    vy: f64,
}

impl Projectile {
    /// Create a projectile from launch angle (radians) and initial speed (m/s).
    pub fn new(angle_rad: f64, speed: f64) -> Self {
        Self {
            vx: speed * angle_rad.cos(),
            vy: speed * angle_rad.sin(),
        }
    }

    /// Infinite lazy iterator of `(x, y)` positions at each time step.
    pub fn positions(&self, dt: f64) -> impl Iterator<Item = (f64, f64)> + '_ {
        (0..).map(move |i| {
            let t = i as f64 * dt;
            (self.vx * t, self.vy * t - 0.5 * G * t * t)
        })
    }

    /// Horizontal range when the projectile returns to y = 0 (metres).
    pub fn range(&self, dt: f64) -> f64 {
        self.positions(dt)
            .take_while(|&(_, y)| y >= 0.0)
            .last()
            .map(|(x, _)| x)
            .unwrap_or(0.0)
    }

    /// Maximum height reached (metres).
    pub fn max_height(&self, dt: f64) -> f64 {
        self.positions(dt)
            .take_while(|&(_, y)| y >= 0.0)
            .map(|(_, y)| y)
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Total time of flight until y returns to 0 (seconds).
    pub fn time_of_flight(&self, dt: f64) -> f64 {
        self.positions(dt)
            .take_while(|&(_, y)| y >= 0.0)
            .count() as f64 * dt
    }

    /// Lazy running-maximum-height stream via `.scan()`.
    pub fn max_height_stream(&self, dt: f64) -> impl Iterator<Item = f64> + '_ {
        self.positions(dt)
            .take_while(|&(_, y)| y >= 0.0)
            .scan(f64::NEG_INFINITY, |max, (_, y)| {
                *max = max.max(y);
                Some(*max)
            })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn range_45_degrees() {
        // Analytical: v²/g = 400/9.81 ≈ 40.775 m
        let p = Projectile::new(std::f64::consts::FRAC_PI_4, 20.0);
        approx::assert_abs_diff_eq!(p.range(0.0001), 40.775, epsilon = 0.01);
    }

    #[test]
    fn max_height_90_degrees() {
        // Analytical: v²/(2g) = 400/19.62 ≈ 20.387 m
        let p = Projectile::new(std::f64::consts::FRAC_PI_2, 20.0);
        approx::assert_abs_diff_eq!(p.max_height(0.0001), 20.0_f64.powi(2) / (2.0 * G), epsilon = 0.01);
    }

    #[test]
    fn complementary_angles_same_range() {
        // 30° and 60° give the same range
        let p30 = Projectile::new(30.0_f64.to_radians(), 20.0);
        let p60 = Projectile::new(60.0_f64.to_radians(), 20.0);
        approx::assert_abs_diff_eq!(p30.range(0.0001), p60.range(0.0001), epsilon = 0.01);
    }

    #[test]
    fn max_range_at_45_degrees() {
        let best = (0..90_u32)
            .map(|deg| {
                let p = Projectile::new((deg as f64).to_radians(), 20.0);
                (deg, p.range(0.001))
            })
            .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
            .unwrap();
        assert_eq!(best.0, 45);
    }
}
