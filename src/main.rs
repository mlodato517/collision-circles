extern crate glutin_window;
extern crate graphics;
extern crate opengl_graphics;
extern crate piston;

use glutin_window::GlutinWindow as Window;
use opengl_graphics::{GlGraphics, OpenGL};
use piston::event_loop::{EventSettings, Events};
use piston::input::{RenderArgs, RenderEvent, UpdateArgs, UpdateEvent};
use piston::window::WindowSettings;
use rand::prelude::*;
use std::fmt;

#[derive(Debug)]
struct Vector2D {
    x: f64,
    y: f64,
}

impl Vector2D {
    fn calc_len(&self) -> f64 {
        self.dot(self).sqrt()
    }
    fn dot(&self, other: &Self) -> f64 {
        self.x * other.x + self.y * other.y
    }
}

impl fmt::Display for Vector2D {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({:.2},{:.2})", self.x, self.y)
    }
}

pub struct App {
    gl: GlGraphics, // OpenGL drawing backend.
    circles: Vec<Vector2D>,
    velocities: Vec<Vector2D>,
    updates: u64,
}

const CIRCLE_RADIUS: f64 = 10.0;
const SCREEN_SIZE: f64 = 200.0;

impl App {
    fn render(&mut self, args: &RenderArgs) {
        use graphics::*;

        const WHITE: [f32; 4] = [1.0, 1.0, 1.0, 1.0];
        const RED: [f32; 4] = [1.0, 0.0, 0.0, 1.0];

        let circles = &self.circles;
        let square = rectangle::square(0.0, 0.0, CIRCLE_RADIUS * 2.0);

        self.gl.draw(args.viewport(), |c, gl| {
            clear(WHITE, gl);

            for circle in circles {
                let transform = c
                    .transform
                    .trans(circle.x, circle.y)
                    .trans(-CIRCLE_RADIUS, -CIRCLE_RADIUS);

                ellipse(RED, square, transform, gl);
            }
        });
    }

    fn update<R>(&mut self, args: &UpdateArgs, rng: &mut R)
    where
        R: Rng,
    {
        // Update each circle's position
        for (c, v) in self.circles.iter_mut().zip(&self.velocities) {
            c.x += v.x * args.dt;
            c.y += v.y * args.dt;
        }

        // If we collide with the boundary, change direction
        for (c, v) in self.circles.iter_mut().zip(self.velocities.iter_mut()) {
            if (c.x + CIRCLE_RADIUS >= SCREEN_SIZE && v.x > 0.0)
                || (c.x - CIRCLE_RADIUS <= 0.0 && v.x < 0.0)
            {
                v.x *= -1.0;
            }
            if (c.y + CIRCLE_RADIUS >= SCREEN_SIZE && v.y > 0.0)
                || (c.y - CIRCLE_RADIUS <= 0.0 && v.y < 0.0)
            {
                v.y *= -1.0;
            }
        }

        // Check for collisions
        let len = self.circles.len();
        for i in 1..len {
            let (first, second) = self.circles.split_at_mut(i);
            let c = &mut first[first.len() - 1];
            for (j, other_c) in second.iter_mut().enumerate() {
                // Fix the indices to match our wackiness
                let j = i + j;
                let i = i - 1;

                let sum_rad_sqr = (2.0 * CIRCLE_RADIUS).powi(2);
                let dist = (c.x - other_c.x).powi(2) + (c.y - other_c.y).powi(2);
                if dist <= sum_rad_sqr {
                    let (first, second) = self.velocities.split_at_mut(i + 1);
                    let v = &mut first[first.len() - 1];
                    let other_v = &mut second[j - i - 1];
                    handle_collision(c, other_c, v, other_v);
                }
            }
        }

        // Make a new circle every so often.
        if self.updates % 100 == 0 {
            self.add_ball(rng);
        }

        self.updates += 1;
    }

    fn add_ball<R>(&mut self, rng: &mut R)
    where
        R: Rng,
    {
        loop {
            let (c, v) = gen_circle(rng);
            let circle_already_there = self.circles.iter().any(|existing_c| {
                let sum_rad = 2.0 * CIRCLE_RADIUS;
                let dist = ((c.x - existing_c.x).powi(2) + (c.y - existing_c.y).powi(2)).sqrt();
                dist <= sum_rad
            });
            if !circle_already_there {
                self.circles.push(c);
                self.velocities.push(v);
                break;
            }
        }
    }
}

fn main() {
    // Change this to OpenGL::V2_1 if not working.
    let opengl = OpenGL::V3_2;

    // Create an Glutin window.
    let mut window: Window = WindowSettings::new("colliding circles", [SCREEN_SIZE, SCREEN_SIZE])
        .graphics_api(opengl)
        .exit_on_esc(true)
        .build()
        .unwrap();

    // Create a new game and run it.
    let mut app = App {
        gl: GlGraphics::new(opengl),
        circles: vec![],
        velocities: vec![],
        updates: 0,
    };

    let mut rng = thread_rng();
    let mut events = Events::new(EventSettings::new());
    while let Some(e) = events.next(&mut window) {
        if let Some(args) = e.render_args() {
            app.render(&args);
        }

        if let Some(args) = e.update_args() {
            app.update(&args, &mut rng);
        }
    }
}

fn gen_circle<R>(rng: &mut R) -> (Vector2D, Vector2D)
where
    R: Rng,
{
    (
        Vector2D {
            x: rng.gen_range(CIRCLE_RADIUS..SCREEN_SIZE - CIRCLE_RADIUS),
            y: rng.gen_range(CIRCLE_RADIUS..SCREEN_SIZE - CIRCLE_RADIUS),
        },
        Vector2D {
            x: rng.gen_range(-100.0..100.0),
            y: rng.gen_range(-100.0..100.0),
        },
    )
}

/// Updates the velocities of two colliding circles.
/// Credit for this algorithm goes to Brian Smith.
///
/// The general idea is that we:
///
/// 1. move to a reference frame where one circle is stationary
/// 2. project the velocity of the moving circle onto the vector between the circles' centers
/// 3. subtract the transferred momentum from the moving circle and add it to the other
///
/// This currently assumes masses are equal.
fn handle_collision(c: &Vector2D, other_c: &Vector2D, v: &mut Vector2D, other_v: &mut Vector2D) {
    // Subtracting the velocities is like moving to an inertial reference frame where
    // one circle is stationary. This avoids having to calculate the transferred
    // momentum in both directions.
    let delta_v = Vector2D {
        x: v.x - other_v.x,
        y: v.y - other_v.y,
    };

    // Get the vector between the centers. Momentum is transferred along this.
    let collision_vec = Vector2D {
        x: c.x - other_c.x,
        y: c.y - other_c.y,
    };
    let collision_len = collision_vec.calc_len();

    // Using the dot product we find the angle between the velocity
    // and the collision vector. Then we project the x and y values
    // onto this axis to find out how much is transferred.
    //
    // Here:
    //   - `c` is the "collision vector" - it is the difference of the centers of the circles.
    //   - `v` is the "delta velocity vector" - it is the difference of the velocities.
    //   - `t` is the "transferred vector" - it is how much of `v` is projected along `c`
    //
    // c . v = c.len * v.len * cos(theta) // "definition" of dot product
    // t.len = v.len * cos(theta)         // definition of cosine
    // t.len = c . v / c.len              // previous two lines simplified
    let transferred_mag = collision_vec.dot(&delta_v) / collision_len;

    // The resulting vector is the collision vector scaled down to the unit vector
    // and then scaled up by the magnitude of the transferred vector.
    let scale = transferred_mag / collision_len;
    let transferred = Vector2D {
        x: collision_vec.x * scale,
        y: collision_vec.y * scale,
    };

    // We subtract it from the initial velocity and add it to the other.
    v.x -= transferred.x;
    v.y -= transferred.y;

    other_v.x += transferred.x;
    other_v.y += transferred.y;
}

#[cfg(test)]
mod tests {
    use super::*;

    macro_rules! assert_close {
        ($f1:expr, $f2:literal) => {
            assert!(
                ($f1 - $f2).abs() <= 1.0e-10,
                "{}({:.2}) - {}",
                stringify!($f1),
                $f1,
                $f2
            );
        };
    }

    #[test]
    fn test_horizontal_collision_one_stopped() {
        let c = Vector2D { x: 0.0, y: 0.0 };
        let mut v = Vector2D { x: 10.0, y: 0.0 };
        let other_c = Vector2D { x: 1.0, y: 0.0 };
        let mut other_v = Vector2D { x: 0.0, y: 0.0 };

        handle_collision(&c, &other_c, &mut v, &mut other_v);
        assert_close!(v.x, 0.0);
        assert_close!(v.y, 0.0);

        assert_close!(other_v.x, 10.0);
        assert_close!(other_v.y, 0.0);
    }

    #[test]
    fn test_horizontal_collision() {
        let c = Vector2D { x: 0.0, y: 0.0 };
        let mut v = Vector2D { x: 10.0, y: 0.0 };
        let other_c = Vector2D { x: 1.0, y: 0.0 };
        let mut other_v = Vector2D { x: -10.0, y: 0.0 };

        handle_collision(&c, &other_c, &mut v, &mut other_v);
        assert_close!(v.x, -10.0);
        assert_close!(v.y, 0.0);

        assert_close!(other_v.x, 10.0);
        assert_close!(other_v.y, 0.0);
    }

    #[test]
    fn test_horizontal_offset_collision() {
        let c = Vector2D { x: 0.0, y: 0.0 };
        let mut v = Vector2D { x: 10.0, y: 0.0 };
        let other_c = Vector2D { x: 1.0, y: 1.0 };
        let mut other_v = Vector2D { x: -10.0, y: 0.0 };

        handle_collision(&c, &other_c, &mut v, &mut other_v);
        assert_close!(v.x, 0.0);
        assert_close!(v.y, -10.0);

        assert_close!(other_v.x, 0.0);
        assert_close!(other_v.y, 10.0);
    }

    #[test]
    fn test_vertical_collision_one_stopped() {
        let c = Vector2D { x: 0.0, y: 0.0 };
        let mut v = Vector2D { x: 0.0, y: 10.0 };
        let other_c = Vector2D { x: 0.0, y: 10.0 };
        let mut other_v = Vector2D { x: 0.0, y: 0.0 };

        handle_collision(&c, &other_c, &mut v, &mut other_v);
        assert_close!(v.x, 0.0);
        assert_close!(v.y, 0.0);

        assert_close!(other_v.x, 0.0);
        assert_close!(other_v.y, 10.0);
    }

    #[test]
    fn test_vertical_collision() {
        let c = Vector2D { x: 0.0, y: 0.0 };
        let mut v = Vector2D { x: 0.0, y: 10.0 };
        let other_c = Vector2D { x: 0.0, y: 10.0 };
        let mut other_v = Vector2D { x: 0.0, y: -10.0 };

        handle_collision(&c, &other_c, &mut v, &mut other_v);
        assert_close!(v.x, 0.0);
        assert_close!(v.y, -10.0);

        assert_close!(other_v.x, 0.0);
        assert_close!(other_v.y, 10.0);
    }

    #[test]
    fn test_45_down_collision() {
        let c = Vector2D { x: 0.0, y: 0.0 };
        let mut v = Vector2D { x: 10.0, y: 10.0 };
        let other_c = Vector2D { x: 2.0, y: 0.0 };
        let mut other_v = Vector2D { x: -10.0, y: 10.0 };

        handle_collision(&c, &other_c, &mut v, &mut other_v);
        assert_close!(v.x, -10.0);
        assert_close!(v.y, 10.0);

        assert_close!(other_v.x, 10.0);
        assert_close!(other_v.y, 10.0);
    }

    #[test]
    fn test_45_up_collision() {
        let c = Vector2D { x: 0.0, y: 0.0 };
        let mut v = Vector2D { x: 10.0, y: -10.0 };
        let other_c = Vector2D { x: 2.0, y: 0.0 };
        let mut other_v = Vector2D { x: -10.0, y: -10.0 };

        handle_collision(&c, &other_c, &mut v, &mut other_v);
        assert_close!(v.x, -10.0);
        assert_close!(v.y, -10.0);

        assert_close!(other_v.x, 10.0);
        assert_close!(other_v.y, -10.0);
    }
}
