use crate::{Circle, Vector2D};

/// Deals with reversing the two collided circles until they are _just_ touching (instead
/// of, perhaps, overlapping a good big). This allows the velocity calculation to have the
/// proper angle while also ensuring the balls don't get "stuck" together.
///
/// When reversing the circles we need to keep track of how much time we've reversed so that,
/// after the collision takes place, we can play that time back. This ensures the circles are
/// updated properly each `dt`.
pub fn handle_pre_collision(
    c: &mut Circle,
    other_c: &mut Circle,
    v: &Vector2D,
    other_v: &Vector2D,
) -> f64 {
    // Move to reference frame where one circle is static.
    let delta_v = Vector2D {
        x: v.x - other_v.x,
        y: v.y - other_v.y,
    };
    let delta_v_len = delta_v.calc_len();

    // We imagine a triangle. One leg is `c1 - c2` and one leg is
    // r1 + r2. A third leg of unknown size points in the
    // direction of `v1 - v2` called `x`. If theta is the angle between `x`
    // and `c1 - c2` then the law of cosines says:
    //
    // (r1+r2)^2 = c_len^2 + x_len^2 - 2*c_len*x_len*cos(theta)
    //
    // Where `c_len` is the length of `c1-c2`. Solving for x using the quadratic formula gives:
    //
    // x = c_len*cos(theta) - sqrt((c_len*cos(theta))^2 - (c_len^2 - (r1+r2)^2))
    //
    // We'll affectionally call c_len*cos(theta) "c_theta"
    // and x "backup_len" because it's how far we have to reverse.
    let delta_c = Vector2D {
        x: c.origin.x - other_c.origin.x,
        y: c.origin.y - other_c.origin.y,
    };
    let diff = delta_c.x.powi(2) + delta_c.y.powi(2) - (c.radius + other_c.radius).powi(2);
    let c_theta = delta_v.dot(&delta_c) / delta_v_len;

    let backup_len = c_theta + (c_theta.powi(2) - diff).sqrt();
    let backup_time = backup_len / delta_v_len;

    c.origin.x -= v.x * backup_time;
    c.origin.y -= v.y * backup_time;

    other_c.origin.x -= other_v.x * backup_time;
    other_c.origin.y -= other_v.y * backup_time;

    // In debug mode assert the radiuses _just_ touch.
    debug_assert!(
        ((c.origin.x - other_c.origin.x).powi(2) + (c.origin.y - other_c.origin.y).powi(2)
            - (c.radius + other_c.radius).powi(2))
        .abs()
            < 1.0e-10
    );

    backup_time
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
pub fn handle_collision(c: &Circle, other_c: &Circle, v: &mut Vector2D, other_v: &mut Vector2D) {
    // Subtracting the velocities is like moving to an inertial reference frame where
    // one circle is stationary. This avoids having to calculate the transferred
    // momentum in both directions.
    let delta_v = Vector2D {
        x: v.x - other_v.x,
        y: v.y - other_v.y,
    };

    // Get the vector between the centers. Momentum is transferred along this.
    let collision_vec = Vector2D {
        x: c.origin.x - other_c.origin.x,
        y: c.origin.y - other_c.origin.y,
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
        ($f1:expr, $f2:expr) => {
            assert!(
                ($f1 - $f2).abs() <= 1.0e-10,
                "{}({:.2}) - {}",
                stringify!($f1),
                $f1,
                $f2
            );
        };
    }

    mod handle_collision {
        use super::*;
        #[test]
        fn test_horizontal_collision_one_stopped() {
            let c = Circle {
                radius: 10.0,
                origin: Vector2D { x: 0.0, y: 0.0 },
            };
            let mut v = Vector2D { x: 10.0, y: 0.0 };
            let other_c = Circle {
                radius: 10.0,
                origin: Vector2D { x: 1.0, y: 0.0 },
            };
            let mut other_v = Vector2D { x: 0.0, y: 0.0 };

            handle_collision(&c, &other_c, &mut v, &mut other_v);
            assert_close!(v.x, 0.0);
            assert_close!(v.y, 0.0);

            assert_close!(other_v.x, 10.0);
            assert_close!(other_v.y, 0.0);
        }

        #[test]
        fn test_horizontal_collision() {
            let c = Circle {
                radius: 10.0,
                origin: Vector2D { x: 0.0, y: 0.0 },
            };
            let mut v = Vector2D { x: 10.0, y: 0.0 };
            let other_c = Circle {
                radius: 10.0,
                origin: Vector2D { x: 1.0, y: 0.0 },
            };
            let mut other_v = Vector2D { x: -10.0, y: 0.0 };

            handle_collision(&c, &other_c, &mut v, &mut other_v);
            assert_close!(v.x, -10.0);
            assert_close!(v.y, 0.0);

            assert_close!(other_v.x, 10.0);
            assert_close!(other_v.y, 0.0);
        }

        #[test]
        fn test_horizontal_offset_collision() {
            let c = Circle {
                radius: 10.0,
                origin: Vector2D { x: 0.0, y: 0.0 },
            };
            let mut v = Vector2D { x: 10.0, y: 0.0 };
            let other_c = Circle {
                radius: 10.0,
                origin: Vector2D { x: 1.0, y: 1.0 },
            };
            let mut other_v = Vector2D { x: -10.0, y: 0.0 };

            handle_collision(&c, &other_c, &mut v, &mut other_v);
            assert_close!(v.x, 0.0);
            assert_close!(v.y, -10.0);

            assert_close!(other_v.x, 0.0);
            assert_close!(other_v.y, 10.0);
        }

        #[test]
        fn test_vertical_collision_one_stopped() {
            let c = Circle {
                radius: 10.0,
                origin: Vector2D { x: 0.0, y: 0.0 },
            };
            let mut v = Vector2D { x: 0.0, y: 10.0 };
            let other_c = Circle {
                radius: 10.0,
                origin: Vector2D { x: 0.0, y: 10.0 },
            };
            let mut other_v = Vector2D { x: 0.0, y: 0.0 };

            handle_collision(&c, &other_c, &mut v, &mut other_v);
            assert_close!(v.x, 0.0);
            assert_close!(v.y, 0.0);

            assert_close!(other_v.x, 0.0);
            assert_close!(other_v.y, 10.0);
        }

        #[test]
        fn test_vertical_collision() {
            let c = Circle {
                radius: 10.0,
                origin: Vector2D { x: 0.0, y: 0.0 },
            };
            let mut v = Vector2D { x: 0.0, y: 10.0 };
            let other_c = Circle {
                radius: 10.0,
                origin: Vector2D { x: 0.0, y: 10.0 },
            };
            let mut other_v = Vector2D { x: 0.0, y: -10.0 };

            handle_collision(&c, &other_c, &mut v, &mut other_v);
            assert_close!(v.x, 0.0);
            assert_close!(v.y, -10.0);

            assert_close!(other_v.x, 0.0);
            assert_close!(other_v.y, 10.0);
        }

        #[test]
        fn test_45_down_collision() {
            let c = Circle {
                radius: 10.0,
                origin: Vector2D { x: 0.0, y: 0.0 },
            };
            let mut v = Vector2D { x: 10.0, y: 10.0 };
            let other_c = Circle {
                radius: 10.0,
                origin: Vector2D { x: 2.0, y: 0.0 },
            };
            let mut other_v = Vector2D { x: -10.0, y: 10.0 };

            handle_collision(&c, &other_c, &mut v, &mut other_v);
            assert_close!(v.x, -10.0);
            assert_close!(v.y, 10.0);

            assert_close!(other_v.x, 10.0);
            assert_close!(other_v.y, 10.0);
        }

        #[test]
        fn test_45_up_collision() {
            let c = Circle {
                radius: 10.0,
                origin: Vector2D { x: 0.0, y: 0.0 },
            };
            let mut v = Vector2D { x: 10.0, y: -10.0 };
            let other_c = Circle {
                radius: 10.0,
                origin: Vector2D { x: 2.0, y: 0.0 },
            };
            let mut other_v = Vector2D { x: -10.0, y: -10.0 };

            handle_collision(&c, &other_c, &mut v, &mut other_v);
            assert_close!(v.x, -10.0);
            assert_close!(v.y, -10.0);

            assert_close!(other_v.x, 10.0);
            assert_close!(other_v.y, -10.0);
        }
    }

    mod handle_pre_collision {
        use super::*;

        #[test]
        fn test_horizontal_pre_collision_one_stopped() {
            let mut c = Circle {
                radius: 10.0,
                origin: Vector2D { x: 0.0, y: 0.0 },
            };
            let v = Vector2D { x: 10.0, y: 0.0 };
            let mut other_c = Circle {
                radius: 10.0,
                origin: Vector2D { x: 10.0, y: 0.0 },
            };
            let other_v = Vector2D { x: 0.0, y: 0.0 };

            handle_pre_collision(&mut c, &mut other_c, &v, &other_v);
            assert_close!(c.origin.x, -10.0);
            assert_close!(c.origin.y, 0.0);

            assert_close!(other_c.origin.x, 10.0);
            assert_close!(other_c.origin.y, 0.0);
        }

        #[test]
        fn test_horizontal_pre_collision() {
            let mut c = Circle {
                radius: 10.0,
                origin: Vector2D { x: 0.0, y: 0.0 },
            };
            let v = Vector2D { x: 10.0, y: 0.0 };
            let mut other_c = Circle {
                radius: 10.0,
                origin: Vector2D { x: 10.0, y: 0.0 },
            };
            let other_v = Vector2D { x: -10.0, y: 0.0 };

            handle_pre_collision(&mut c, &mut other_c, &v, &other_v);
            assert_close!(c.origin.x, -5.0);
            assert_close!(c.origin.y, 0.0);

            assert_close!(other_c.origin.x, 15.0);
            assert_close!(other_c.origin.y, 0.0);
        }

        #[test]
        fn test_vertical_pre_collision_one_stopped() {
            let mut c = Circle {
                radius: 10.0,
                origin: Vector2D { x: 0.0, y: 0.0 },
            };
            let v = Vector2D { x: 0.0, y: 10.0 };
            let mut other_c = Circle {
                radius: 10.0,
                origin: Vector2D { x: 0.0, y: 10.0 },
            };
            let other_v = Vector2D { x: 0.0, y: 0.0 };

            handle_pre_collision(&mut c, &mut other_c, &v, &other_v);
            assert_close!(c.origin.x, 0.0);
            assert_close!(c.origin.y, -10.0);

            assert_close!(other_c.origin.x, 0.0);
            assert_close!(other_c.origin.y, 10.0);
        }

        #[test]
        fn test_horizontal_pre_collision_with_vertical_velocity() {
            let mut c = Circle {
                radius: 10.0,
                origin: Vector2D { x: 0.0, y: 0.0 },
            };
            let v = Vector2D { x: 0.0, y: 10.0 };
            let mut other_c = Circle {
                radius: 10.0,
                origin: Vector2D { x: 10.0, y: 0.0 },
            };
            let other_v = Vector2D { x: 0.0, y: 0.0 };

            handle_pre_collision(&mut c, &mut other_c, &v, &other_v);
            assert_close!(c.origin.x, 0.0);
            assert_close!(c.origin.y, -10.0 * 3.0f64.sqrt());

            assert_close!(other_c.origin.x, 10.0);
            assert_close!(other_c.origin.y, 0.0);
        }

        #[test]
        fn test_horizontal_pre_collision_one_stopped_past_center() {
            let mut c = Circle {
                radius: 10.0,
                origin: Vector2D { x: 11.0, y: 0.0 },
            };
            let v = Vector2D { x: 10.0, y: 0.0 };
            let mut other_c = Circle {
                radius: 10.0,
                origin: Vector2D { x: 10.0, y: 0.0 },
            };
            let other_v = Vector2D { x: 0.0, y: 0.0 };

            handle_pre_collision(&mut c, &mut other_c, &v, &other_v);
            assert_close!(c.origin.x, -10.0);
            assert_close!(c.origin.y, 0.0);

            assert_close!(other_c.origin.x, 10.0);
            assert_close!(other_c.origin.y, 0.0);
        }

        #[test]
        fn test_horizontal_pre_collision_one_stopped_past_center_backward() {
            let mut c = Circle {
                radius: 10.0,
                origin: Vector2D { x: 11.0, y: 0.0 },
            };
            let v = Vector2D { x: -10.0, y: 0.0 };
            let mut other_c = Circle {
                radius: 10.0,
                origin: Vector2D { x: 10.0, y: 0.0 },
            };
            let other_v = Vector2D { x: 0.0, y: 0.0 };

            handle_pre_collision(&mut c, &mut other_c, &v, &other_v);
            assert_close!(c.origin.x, 30.0);
            assert_close!(c.origin.y, 0.0);

            assert_close!(other_c.origin.x, 10.0);
            assert_close!(other_c.origin.y, 0.0);
        }
    }
}
