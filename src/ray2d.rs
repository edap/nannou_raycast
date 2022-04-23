#![allow(dead_code)]
use crate::bounding_volume;
use nannou::prelude::*;

#[derive(Debug, Clone, Copy)]
pub struct Ray2d {
    pub orig: Vec2,
    pub dir: Vec2,
}

impl Default for Ray2d {
    fn default() -> Self {
        Self {
            orig: vec2(0.0, 0.0),
            dir: vec2(1.0, 0.0),
        }
    }
}

impl Ray2d {
    pub fn new(orig: Vec2, dir: Vec2) -> Self {
        Self {
            orig,
            dir: dir.normalize(),
        }
    }

    pub fn reflect(&self, surface_normal: Vec2) -> Vec2 {
        //I - 2.0 * dot(N, I) * N
        // https://www.khronos.org/registry/OpenGL-Refpages/gl4/html/reflect.xhtml
        //
        self.dir - surface_normal.normalize() * (2.0 * surface_normal.dot(self.dir))
    }

    pub fn refract(&self, surface_normal: Vec2, ior: f32) -> Vec2 {
        // https://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-to-shading/reflection-refraction-fresnel

        let mut cosi = clamp(-1.0, 1.0, self.dir.dot(surface_normal));
        let (mut etai, mut etat) = (1.0, ior);
        let mut n = surface_normal.normalize();
        if cosi < 0.0 {
            cosi = -cosi;
        } else {
            std::mem::swap(&mut etai, &mut etat);
            n = -surface_normal;
        }
        let eta = etai / etat;
        let k = 1.0 - eta * eta * (1.0 - cosi * cosi);
        if k < f32::zero() {
            self.dir.normalize() * 0.0
        } else {
            self.dir.normalize() * eta + n.normalize() * (eta * cosi - k.sqrt())
        }
    }

    // in case of material like glass, that are both refractive and reflective, fresnel equation find out how much
    // light is refracted and how much light is reflected
    // reference https://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-to-shading/reflection-refraction-fresnel
    pub fn fresnel(&self, surface_normal: Vec2, ior: f32) -> f32 {
        let i_dot_n = self.dir.dot(surface_normal) as f32;
        let mut eta_i = 1.0;
        let mut eta_t = ior;
        if i_dot_n > 0.0 {
            eta_i = eta_t;
            eta_t = 1.0;
        }

        let sin_t = eta_i / eta_t * (1.0 - i_dot_n * i_dot_n).max(0.0).sqrt();
        if sin_t > 1.0 {
            //Total internal reflection
            1.0
        } else {
            let cos_t = (1.0 - sin_t * sin_t).max(0.0).sqrt();
            let cos_i = cos_t.abs();
            let r_s = ((eta_t * cos_i) - (eta_i * cos_t)) / ((eta_t * cos_i) + (eta_i * cos_t));
            let r_p = ((eta_i * cos_i) - (eta_t * cos_t)) / ((eta_i * cos_i) + (eta_t * cos_t));
            (r_s * r_s + r_p * r_p) / 2.0
        }
    }

    pub fn draw(&self, draw: &Draw, mag: f32, weight: f32, col: Rgb) {
        draw.arrow()
            .color(col)
            .weight(weight)
            .start(self.orig)
            .end(self.dir.normalize() * mag);
    }

    pub fn look_at(&mut self, x: f32, y: f32) {
        self.dir.x = x - self.orig.x;
        self.dir.y = y - self.orig.y;
        self.dir = self.dir.normalize();
    }

    pub fn set_dir_from_angle(&mut self, angle_in_radians: f32) {
        self.dir.x = angle_in_radians.cos();
        self.dir.y = angle_in_radians.sin();
    }

    pub fn intersect_segment(&self, x1: &f32, y1: &f32, x2: &f32, y2: &f32) -> Option<f32> {
        let x3 = self.orig.x;
        let y3 = self.orig.y;
        let x4 = self.orig.x + self.dir.x;
        let y4 = self.orig.y + self.dir.y;
        let den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);

        let t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / den;
        let u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / den;

        if den != 0.0 && t > 0.0 && t < 1.0 && u > 0.0 {
            Some(u)
        } else {
            None
        }
    }

    pub fn intersect_polyline(&self, points: &[Vec2]) -> Option<(f32, Vec2)> {
        if points.len() <= 1 {
            return None;
        }
        let mut distance: f32 = Float::infinity();
        let mut surface_normal: Vec2 = vec2(0.0, 0.0);
        for index in 0..points.len() - 1 {
            if let Some(collision_distance) = self.intersect_segment(
                &points[index].x,
                &points[index].y,
                &points[index + 1].x,
                &points[index + 1].y,
            ) {
                if collision_distance < distance {
                    let segment_dir = (points[index] - points[index + 1]).normalize();
                    surface_normal = vec2(segment_dir.y, -segment_dir.x);
                    distance = collision_distance;
                }
            }
        }
        if distance < Float::infinity() {
            Some((distance, surface_normal))
        } else {
            None
        }
    }

    pub fn intersect_circle(&self, center: &Vec2, radius: &f32) -> Option<(f32, Vec2)> {
        let l = *center - self.orig;
        let adj = l.dot(self.dir);
        let d2 = l.dot(l) - (adj * adj);
        let radius2 = radius * radius;
        if d2 > radius2 {
            return None;
        }
        let thc = (radius2 - d2).sqrt();
        let t0 = adj - thc;
        let t1 = adj + thc;
        if t0 < f32::zero() && t1 < f32::zero() {
            return None;
        }
        let inside = self.orig.distance(*center) <= *radius;
        let distance = if t0 < t1 && !inside { t0 } else { t1 };
        // TODO, test this surface_normal
        let surface_normal = ((self.orig * distance) - *center).normalize();
        Some((distance, surface_normal))
    }

    // credits to https://github.com/rustgd/collision-rs/blob/master/src/volume/aabb/aabb2.rs
    pub fn intersect_aabb(&self, min: &Vec2, max: &Vec2) -> Option<f32> {
        let mut tmax: f32 = Float::infinity();
        let mut tmin: f32 = Float::neg_infinity();
        if self.dir.x != f32::zero() {
            let tx1 = (min.x - self.orig.x) / self.dir.x;
            let tx2 = (max.x - self.orig.x) / self.dir.x;
            tmin = tmin.max(tx1.min(tx2));
            tmax = tmax.min(tx1.max(tx2));
        } else if self.orig.x <= min.x || self.orig.x >= max.x {
            return None;
        }

        if self.dir.y != 0.0 {
            let ty1 = (min.y - self.orig.y) / self.dir.y;
            let ty2 = (max.y - self.orig.y) / self.dir.y;
            tmin = tmin.max(ty1.min(ty2));
            tmax = tmax.min(ty1.max(ty2));
        } else if self.orig.y <= min.y || self.orig.y >= max.y {
            return None;
        }

        if (tmin < f32::zero() && tmax < f32::zero()) || tmax < tmin {
            None
        } else {
            let t = if tmin >= f32::zero() { tmin } else { tmax };
            Some(t)
        }
    }

    pub fn intersect_bounding_volume(
        &self,
        volume: &bounding_volume::BoundingVolume,
    ) -> Option<f32> {
        match volume {
            bounding_volume::BoundingVolume::Circle { position, radius } => {
                match self.intersect_circle(position, radius) {
                    Some((dist, _surface_normal)) => Some(dist),
                    _ => None,
                }
            }
            bounding_volume::BoundingVolume::Aabb { min, max } => self.intersect_aabb(min, max),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const IOR_GLASS: f32 = 1.500;
    const IOR_GOLD: f32 = 0.470;
    const IOR_DIAMOND: f32 = 2.418;

    #[test]
    fn take_parametesrs_for_initialization() {
        let r = Ray2d::new(vec2(2.0, 0.0), vec2(0.0, 1.0));
        assert_eq!(r.dir, vec2(0.0, 1.0));
        assert_eq!(r.orig, vec2(2.0, 0.0));
    }

    #[test]
    fn have_a_default_direction() {
        let r = Ray2d::default();
        assert_eq!(r.dir, vec2(1.0, 0.0));
    }

    #[test]
    fn rotate_the_vector_towards_the_look_at_point() {
        let mut r = Ray2d::default();
        r.look_at(-1.0, -1.0);
        assert_eq!(r.dir, vec2(-0.70710677, -0.70710677));
    }

    #[test]
    fn calculate_the_refraction_vector_for_different_ior_coeficient() {
        let ray = Ray2d::new(vec2(0.0, 0.0), vec2(1.0, 1.0));
        let intersection = ray.intersect_circle(&vec2(2.0, 1.7), &1.0).unwrap();
        // GLASS
        assert_eq!(
            ray.refract(intersection.1, IOR_GLASS),
            vec2(0.72593915, 0.687759)
        );
        // GOLD
        assert_eq!(
            ray.refract(intersection.1, IOR_GOLD),
            vec2(0.63922864, 0.76901674)
        );
        // DIAMOND
        assert_eq!(
            ray.refract(intersection.1, IOR_DIAMOND),
            vec2(0.7398675, 0.6727526)
        );
    }

    #[test]
    fn calculate_the_reflection_vector() {
        let ray = Ray2d::new(vec2(0.0, 0.0), vec2(1.0, 1.0));
        let intersection = ray.intersect_circle(&vec2(2.0, 2.0), &1.0).unwrap();

        assert_eq!(ray.reflect(intersection.1), vec2(-0.70710677, -0.70710677));
    }

    #[test]
    fn calculate_the_fresnel_coeficient_for_different_ior_coeficient() {
        let ray = Ray2d::new(vec2(0.0, 0.0), vec2(1.0, 1.0));
        let intersection = ray.intersect_circle(&vec2(2.0, 1.7), &1.0).unwrap();
        // GLASS
        assert_eq!(ray.fresnel(intersection.1, IOR_GLASS), 0.040000003);
        // GOLD
        assert_eq!(ray.fresnel(intersection.1, IOR_GOLD), 0.12999213);
        // DIAMOND
        assert_eq!(ray.fresnel(intersection.1, IOR_DIAMOND), 0.17211089);
    }

    #[test]
    fn calculate_the_intersection_with_a_segment() {
        let mut r = Ray2d::default();
        r.dir = vec2(0.0, 1.0);
        let start_segment = vec2(-1.0, 2.0);
        let end_segment = vec2(1.0, 2.0);
        let distance_to_intersection = r.intersect_segment(
            &start_segment.x,
            &start_segment.y,
            &end_segment.x,
            &end_segment.y,
        );
        assert_eq!(distance_to_intersection.unwrap(), 2.0);
    }

    #[test]
    fn return_none_if_there_isnt_an_intersection_with_a_segment() {
        let mut r = Ray2d::default();
        r.dir = vec2(0.0, -1.0);
        let start_segment = vec2(-1.0, 2.0);
        let end_segment = vec2(1.0, 2.0);
        let distance_to_intersection = r.intersect_segment(
            &start_segment.x,
            &start_segment.y,
            &end_segment.x,
            &end_segment.y,
        );
        assert_eq!(distance_to_intersection, None);
    }

    #[test]
    fn calculate_the_intersection_with_a_circle() {
        let ray = Ray2d::new(vec2(0.0, 0.0), vec2(1.0, 1.0));
        let intersection = ray.intersect_circle(&vec2(2.0, 2.0), &1.0);
        assert_eq!(
            intersection.unwrap(),
            (1.8284273, vec2(-0.70710677, -0.70710677))
        )
    }

    #[test]
    fn return_none_if_there_is_no_intersection_with_a_circle() {
        let ray = Ray2d::new(vec2(0.0, 0.0), vec2(0.0, 1.0));
        let intersection = ray.intersect_circle(&vec2(1.0, 1.0), &0.5);
        assert_eq!(intersection, None)
    }

    #[test]
    fn calculate_the_intersection_with_a_polyline() {
        let ray = Ray2d::new(vec2(0.0, 0.0), vec2(1.0, 1.0));
        let radius = 150.0;
        let mut points: Vec<Vec2> = Vec::new();
        for i in 0..=180 {
            let radian = deg_to_rad(i as f32);
            let x = radian.sin() * radius;
            let y = radian.cos() * radius;
            points.push(vec2(x, y));
        }
        let intersection = ray.intersect_polyline(&points);
        assert_eq!(intersection.unwrap(), (150.0, vec2(0.70090896, 0.7132508)))
    }

    #[test]
    fn return_none_if_there_is_no_intersection_with_a_polyline() {
        let ray = Ray2d::new(vec2(0.0, 0.0), vec2(-1.0, -1.0));
        let radius = 150.0;
        let mut points: Vec<Vec2> = Vec::new();
        for i in 0..=180 {
            let radian = deg_to_rad(i as f32);
            let x = radian.sin() * radius;
            let y = radian.cos() * radius;
            points.push(vec2(x, y));
        }
        let intersection = ray.intersect_polyline(&points);
        assert_eq!(intersection, None)
    }

    #[test]
    fn calculate_the_intersection_with_a_bounding_volume_circle() {
        let ray = Ray2d::new(vec2(0.0, 0.0), vec2(1.0, 1.0));
        let bv = bounding_volume::BoundingVolume::Circle {
            position: vec2(2.0, 2.0),
            radius: 1.0,
        };
        let intersection = ray.intersect_bounding_volume(&bv);
        assert_eq!(intersection.unwrap(), 1.8284273);
    }

    #[test]
    fn return_none_if_there_is_no_intersection_with_a_bounding_volume_circle() {
        let ray = Ray2d::new(vec2(0.0, 0.0), vec2(-1.0, -1.0));
        let bv = bounding_volume::BoundingVolume::Circle {
            position: vec2(2.0, 2.0),
            radius: 1.0,
        };
        let intersection = ray.intersect_bounding_volume(&bv);
        assert_eq!(intersection, None);
    }
    #[test]
    fn calculate_the_intersection_with_a_bounding_volume_aabb() {
        let ray = Ray2d::new(vec2(0.0, 0.0), vec2(1.0, 1.0));
        let bv = bounding_volume::BoundingVolume::Aabb {
            min: vec2(1.0, 1.0),
            max: vec2(2.0, 2.0),
        };
        let intersection = ray.intersect_bounding_volume(&bv);
        assert_eq!(intersection.unwrap(), 1.4142135);
    }

    #[test]
    fn return_none_if_there_is_no_intersection_with_a_bounding_volume_aabb() {
        let ray = Ray2d::new(vec2(-0.1, -0.1), vec2(-1.0, -1.0));
        let bv = bounding_volume::BoundingVolume::Aabb {
            min: vec2(0.0, 0.0),
            max: vec2(1.0, 1.0),
        };
        let intersection = ray.intersect_bounding_volume(&bv);
        assert_eq!(intersection, None);
    }
}
