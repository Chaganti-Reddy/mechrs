//! Print a range table for angles 0°–90° and show the sweep as ASCII bars.
use mechrs::kinematics::Projectile;

fn main() {
    let speed = 20.0;
    println!("Projectile range sweep (v = {speed} m/s)\n");
    println!("{:>5}  {:>10}  {}", "Angle", "Range (m)", "");
    println!("{:-<5}  {:-<10}  {:-<40}", "", "", "");

    let max_range = speed * speed / 9.81;
    for deg in (0..=90).step_by(5) {
        let p = Projectile::new((deg as f64).to_radians(), speed);
        let range = p.range(0.001);
        let bars = ((range / max_range) * 40.0).round() as usize;
        println!("{deg:>4}°  {range:>10.3}  {}", "█".repeat(bars));
    }
}
