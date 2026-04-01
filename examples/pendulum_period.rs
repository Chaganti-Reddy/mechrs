//! Compare analytical vs numerical period for 5 amplitudes.
use mechrs::pendulum::Pendulum;

fn main() {
    let p = Pendulum::new(1.0);
    let analytical = p.analytical_period();

    println!("Pendulum period comparison (L = 1.0 m)\n");
    println!(
        "{:>8}  {:>12}  {:>12}  {:>10}",
        "Angle", "Analytical", "Numerical", "Error %"
    );
    println!("{:-<8}  {:-<12}  {:-<12}  {:-<10}", "", "", "", "");

    for &deg in &[5.0_f64, 15.0, 30.0, 60.0, 80.0] {
        let theta0 = deg.to_radians();
        let numerical = p.numerical_period(theta0, 0.0001);
        let error = (numerical - analytical).abs() / analytical * 100.0;
        println!("{deg:>7.1}°  {analytical:>12.6}  {numerical:>12.6}  {error:>9.4}%");
    }
}
