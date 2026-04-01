//! ASCII orbital trace — Verlet orbit demo.
use mechrs::orbital::{circular_speed, verlet_orbit, Vec2};

fn main() {
    let gm = 1.327e20_f64;
    let r = 1.496e11_f64;
    let v = circular_speed(gm, r);
    let r0 = Vec2::new(r, 0.0);

    // Elliptical orbit: reduce speed by 20%
    let v0 = Vec2::new(0.0, v * 0.8);
    let dt = 3600.0 * 24.0;

    // Collect one orbit (approx)
    let points: Vec<(f64, f64)> = verlet_orbit(gm, r0, v0, dt)
        .take(500)
        .map(|(pos, _)| (pos.x / r, pos.y / r))
        .collect();

    // ASCII render on a 60×30 grid
    let width = 60usize;
    let height = 30usize;
    let mut grid = vec![vec![' '; width]; height];

    for (x, y) in &points {
        let col = ((*x + 1.6) / 3.2 * (width as f64 - 1.0)).round() as isize;
        let row = ((-*y + 1.6) / 3.2 * (height as f64 - 1.0)).round() as isize;
        if col >= 0 && col < width as isize && row >= 0 && row < height as isize {
            grid[row as usize][col as usize] = '*';
        }
    }

    // Mark sun at origin
    let sc = ((1.6 / 3.2) * (width as f64 - 1.0)).round() as usize;
    let sr = ((1.6 / 3.2) * (height as f64 - 1.0)).round() as usize;
    if sr < height && sc < width {
        grid[sr][sc] = '☀';
    }

    println!("Elliptical orbit (v = 0.8 × circular speed)\n");
    for row in &grid {
        println!("{}", row.iter().collect::<String>());
    }
}
