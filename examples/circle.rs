use nannou::prelude::*;
use nannou_raycast::ray2d::Ray2d;

const IOR_GLASS: f32 = 1.500;

fn main() {
    nannou::app(model).update(update).run();
}

struct Model {
    _window: window::Id,
    circle_pos: Vec2,
    circle_radius: f32,
}

fn model(app: &App) -> Model {
    let _window = app.new_window().view(view).build().unwrap();
    let circle_pos = vec2(0.0, 0.0);
    let circle_radius = 100.0;
    Model {
        _window,
        circle_pos,
        circle_radius,
    }
}

fn update(_app: &App, _model: &mut Model, _update: Update) {}

fn view(app: &App, model: &Model, frame: Frame) {
    let draw = app.draw();
    draw.background().color(BLACK);
    draw.ellipse()
        .x_y(model.circle_pos.x, model.circle_pos.y)
        .radius(model.circle_radius)
        .color(STEELBLUE);

    for i in (0..360).step_by(6) {
        let mut ray = Ray2d::default();
        ray.orig = app.mouse.position();
        ray.set_dir_from_angle(deg_to_rad(i as f32));

        if let Some((dist_to_collision, surface_normal)) =
            ray.intersect_circle(&model.circle_pos, &model.circle_radius)
        {
            let collision_point = ray.orig + ray.dir * dist_to_collision;

            // draw collisions
            draw.arrow()
                .color(WHITE)
                .weight(1.0)
                .start(ray.orig)
                .end(collision_point);

            // do not draw reflections and refractions if the mouse is inside the circle, just to avoid confusion
            if app.mouse.position().distance(model.circle_pos) > model.circle_radius {
                // draw refractions
                let refraction = ray.refract(surface_normal, IOR_GLASS);
                draw.arrow()
                    .color(BLUEVIOLET)
                    .weight(1.0)
                    .start(collision_point)
                    .end(collision_point + refraction * 2000.0);

                //draw reflections
                let reflection = ray.reflect(surface_normal);
                draw.arrow()
                    .color(STEELBLUE)
                    .weight(1.0)
                    .start(collision_point)
                    .end(collision_point + reflection * 2000.0);
            }
        } else {
            ray.draw(&draw, 2000.0, 1.0, rgb(1.0, 1.0, 1.0));
        }
    }

    draw.to_frame(app, &frame).unwrap();
}
