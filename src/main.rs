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

#[derive(Debug)]
struct Point {
    x: f64,
    y: f64,
}
#[derive(Debug)]
struct Velocity {
    x: f64,
    y: f64,
}
pub struct App {
    gl: GlGraphics, // OpenGL drawing backend.
    circles: Vec<Point>,
    velocities: Vec<Velocity>,
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

fn gen_circle<R>(rng: &mut R) -> (Point, Velocity)
where
    R: Rng,
{
    (
        Point {
            x: rng.gen_range(CIRCLE_RADIUS..SCREEN_SIZE - CIRCLE_RADIUS),
            y: rng.gen_range(CIRCLE_RADIUS..SCREEN_SIZE - CIRCLE_RADIUS),
        },
        Velocity {
            x: rng.gen_range(-100.0..100.0),
            y: rng.gen_range(-100.0..100.0),
        },
    )
}
