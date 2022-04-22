use nannou::prelude::Vec2;

#[derive(Debug, Copy, Clone)]
pub enum BoundingVolume {
    Circle { position: Vec2, radius: f32 },
    Aabb { min: Vec2, max: Vec2 },
}