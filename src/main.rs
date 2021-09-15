extern crate glutin_window;
extern crate graphics;
extern crate opengl_graphics;
extern crate piston;

mod collision;

use glutin_window::GlutinWindow as Window;
use opengl_graphics::{GlGraphics, OpenGL};
use piston::event_loop::{EventSettings, Events};
use piston::input::{RenderArgs, RenderEvent, UpdateArgs, UpdateEvent};
use piston::window::WindowSettings;
use rand::prelude::*;
use std::fmt;

#[derive(Debug)]
pub struct Circle {
    origin: Vector2D,
    radius: f64,
}

#[derive(Debug)]
pub struct Vector2D {
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
    circles: Vec<Circle>,
    velocities: Vec<Vector2D>,
    updates: u64,
}

const SCREEN_SIZE: f64 = 200.0;

impl App {
    fn render(&mut self, args: &RenderArgs) {
        use graphics::*;

        const WHITE: [f32; 4] = [1.0, 1.0, 1.0, 1.0];
        const RED: [f32; 4] = [1.0, 0.0, 0.0, 1.0];

        let circles = &self.circles;

        self.gl.draw(args.viewport(), |c, gl| {
            clear(WHITE, gl);

            for circle in circles {
                let transform = c
                    .transform
                    .trans(circle.origin.x, circle.origin.y)
                    .trans(-circle.radius, -circle.radius);

                let circle_box = rectangle::square(0.0, 0.0, circle.radius * 2.0);
                ellipse(RED, circle_box, transform, gl);
            }
        });
    }

    fn update<R>(&mut self, args: &UpdateArgs, rng: &mut R)
    where
        R: Rng,
    {
        // Update each circle's position
        for (c, v) in self.circles.iter_mut().zip(&self.velocities) {
            c.origin.x += v.x * args.dt;
            c.origin.y += v.y * args.dt;
        }

        // If we collide with the boundary, change direction
        for (c, v) in self.circles.iter_mut().zip(self.velocities.iter_mut()) {
            if (c.origin.x + c.radius >= SCREEN_SIZE && v.x > 0.0)
                || (c.origin.x - c.radius <= 0.0 && v.x < 0.0)
            {
                v.x *= -1.0;
            }
            if (c.origin.y + c.radius >= SCREEN_SIZE && v.y > 0.0)
                || (c.origin.y - c.radius <= 0.0 && v.y < 0.0)
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

                // TODO Consider extracting this entire block or passing more information
                // (`diff`, `delta_v`, `delta_v`, etc.) into more functions to reduce duplication.
                let sum_rad_sqr = (c.radius + other_c.radius).powi(2);
                let dist = (c.origin.x - other_c.origin.x).powi(2)
                    + (c.origin.y - other_c.origin.y).powi(2);
                let diff = dist - sum_rad_sqr;
                if diff < 0.0 {
                    let (first, second) = self.velocities.split_at_mut(i + 1);
                    let v = &mut first[first.len() - 1];
                    let other_v = &mut second[j - i - 1];

                    let backup_time = collision::handle_pre_collision(c, other_c, v, other_v);
                    collision::handle_collision(c, other_c, v, other_v);
                    c.origin.x += v.x * backup_time;
                    c.origin.y += v.y * backup_time;
                    other_c.origin.x += other_v.x * backup_time;
                    other_c.origin.y += other_v.y * backup_time;
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
            let c = gen_circle(rng);
            let v = gen_velocity(rng);
            let circle_already_there = self.circles.iter().any(|existing_c| {
                let sum_rad = c.radius + existing_c.radius;
                let dist = ((c.origin.x - existing_c.origin.x).powi(2)
                    + (c.origin.y - existing_c.origin.y).powi(2))
                .sqrt();
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

/// Generates a random circle with a random radius somewhere on the screen.
fn gen_circle<R>(rng: &mut R) -> Circle
where
    R: Rng,
{
    let radius = rng.gen_range(10.0..20.0);
    Circle {
        radius,
        origin: Vector2D {
            x: rng.gen_range(radius..SCREEN_SIZE - radius),
            y: rng.gen_range(radius..SCREEN_SIZE - radius),
        },
    }
}

/// Generates a random velocity vector.
fn gen_velocity<R>(rng: &mut R) -> Vector2D
where
    R: Rng,
{
    Vector2D {
        x: rng.gen_range(-100.0..100.0),
        y: rng.gen_range(-100.0..100.0),
    }
}
