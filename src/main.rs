//! mechrs CLI — run physics simulations from the command line.

use clap::{Parser, Subcommand};
use mechrs::{kinematics::Projectile, oscillator::Oscillator, orbital, pendulum::Pendulum};

#[derive(Parser)]
#[command(name = "mechrs", about = "Classical mechanics as lazy Rust iterators", version)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Projectile motion under gravity.
    Projectile {
        /// Launch angle in degrees.
        #[arg(long, default_value_t = 45.0)]
        angle: f64,
        /// Initial speed in m/s.
        #[arg(long, default_value_t = 20.0)]
        speed: f64,
    },
    /// Nonlinear pendulum (RK4).
    Pendulum {
        /// Pendulum length in metres.
        #[arg(long, default_value_t = 1.0)]
        length: f64,
        /// Initial angle in degrees.
        #[arg(long, default_value_t = 10.0)]
        angle: f64,
    },
    /// Spring-mass oscillator.
    Oscillator {
        /// Spring constant k (N/m).
        #[arg(long, default_value_t = 10.0)]
        k: f64,
        /// Mass in kg.
        #[arg(long, default_value_t = 1.0)]
        mass: f64,
        /// Initial displacement in metres.
        #[arg(long, default_value_t = 1.0)]
        x0: f64,
        /// Damping coefficient (0 = undamped).
        #[arg(long, default_value_t = 0.0)]
        damping: f64,
    },
    /// Two-body Kepler orbit (Velocity Verlet).
    Orbit {
        /// Orbital radius in AU.
        #[arg(long, default_value_t = 1.0)]
        radius: f64,
        /// GM in units of 10^20 m³/s².
        #[arg(long, default_value_t = 1.327)]
        gm: f64,
    },
}

fn main() {
    let cli = Cli::parse();

    match cli.command {
        Command::Projectile { angle, speed } => {
            let p = Projectile::new(angle.to_radians(), speed);
            let range = p.range(0.0001);
            let height = p.max_height(0.0001);
            let tof = p.time_of_flight(0.0001);
            println!("Projectile: {angle}° at {speed} m/s");
            println!("  Range:           {range:.3} m");
            println!("  Max height:      {height:.3} m");
            println!("  Time of flight:  {tof:.3} s");
        }

        Command::Pendulum { length, angle } => {
            let p = Pendulum::new(length);
            let theta0 = angle.to_radians();
            let analytical = p.analytical_period();
            let numerical = p.numerical_period(theta0, 0.0001);
            let error_pct = (numerical - analytical).abs() / analytical * 100.0;
            println!("Pendulum: L={length} m, θ₀={angle}°");
            println!("  Analytical period (small-angle): {analytical:.4} s");
            println!("  Numerical period:                {numerical:.4} s");
            println!("  Small-angle error:               {error_pct:.4}%");
        }

        Command::Oscillator { k, mass, x0, damping } => {
            let osc = if damping > 0.0 {
                Oscillator::damped(k, mass, damping)
            } else {
                Oscillator::new(k, mass)
            };
            let analytical = osc.analytical_period();
            let numerical = osc.numerical_period(x0, 0.0001);
            println!("Oscillator: k={k} N/m, m={mass} kg, x0={x0} m, b={damping}");
            println!("  Analytical period: {analytical:.4} s");
            if let Some(num) = numerical {
                println!("  Numerical period:  {num:.4} s");
            } else {
                println!("  Numerical period:  (overdamped, no oscillation)");
            }

            // ASCII displacement plot
            println!("\n  Displacement (first 2 periods, normalized):");
            let steps = (analytical * 2.0 / 0.01) as usize;
            let max_w = 40;
            for (t, x) in osc.states(x0, 0.0, 0.01).take(steps) {
                let bars = ((x / x0) * max_w as f64).round() as i32;
                let center = max_w as i32;
                let line: String = (0..(2 * max_w + 1) as i32)
                    .map(|i| {
                        if i == center { '|' }
                        else if bars >= 0 && i > center && i <= center + bars { '█' }
                        else if bars < 0 && i >= center + bars && i < center { '█' }
                        else { ' ' }
                    })
                    .collect();
                println!("  {t:5.2}s {line}");
            }
        }

        Command::Orbit { radius, gm } => {
            let r = radius * 1.496e11; // AU to metres
            let gm_si = gm * 1e20;
            let v = orbital::circular_speed(gm_si, r);
            let r0 = orbital::Vec2::new(r, 0.0);
            let v0 = orbital::Vec2::new(0.0, v);
            let dt = 3600.0 * 24.0; // 1 day

            let e0 = orbital::specific_energy(gm_si, r0, v0);
            let (r_min, r_max, _) = orbital::verlet_orbit(gm_si, r0, v0, dt)
                .take(365)
                .map(|(pos, _)| pos.norm())
                .fold((f64::INFINITY, f64::NEG_INFINITY, 0usize), |(mn, mx, n), r| (mn.min(r), mx.max(r), n + 1));

            println!("Orbit: r={radius} AU, GM={gm}×10²⁰ m³/s²");
            println!("  Circular speed:  {:.3} km/s", v / 1000.0);
            println!("  Specific energy: {e0:.3e} J/kg");
            println!("  Min radius: {:.4} AU", r_min / 1.496e11);
            println!("  Max radius: {:.4} AU", r_max / 1.496e11);
            println!("  (365-day simulation)");
        }
    }
}
