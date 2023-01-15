use approx::{AbsDiffEq, RelativeEq, UlpsEq, ulps_eq};
use bevy::prelude::*;
use crate::{Curve, ApproxVec};


/// Represents a circular segment
#[derive(Debug, PartialEq)]
pub struct Arc {
	center: Vec3,
	axis1: Vec3,
	axis2: Vec3,
	radius: f32,
	angle: f32,
	arc_len: f32,
}

impl Arc {
	pub fn new(center: Vec3, axis1: Vec3, axis2: Vec3, radius: f32, angle: f32, arc_len: f32) -> Self {
		Self{ center, axis1, axis2, radius, angle, arc_len }
	}

	pub fn from_points(start: Vec3, tan: Vec3, end: Vec3) -> Self {
		let center;
		let axis1;
		let axis2;
		let radius;
		let angle;
		let arc_len;

		assert!(tan.is_normalized());
		let start_to_mid = end - start;
		let normal = start_to_mid.cross(tan);
		let perp_axis = tan.cross(normal);
		let denominator = 2.0 * perp_axis.dot(start_to_mid);

		if ulps_eq!(0.0, denominator) {
			// The radius is infinite, so use a straight line. Place the center in the middle of the line
			center = start + start_to_mid * 0.5;
			axis1 = start_to_mid.normalize();
			axis2 = axis1;
			radius = f32::INFINITY;
			angle = 0.0;
			arc_len = start_to_mid.length();
		} else {
			let center_dist = start_to_mid.dot(start_to_mid) / denominator;
			center = start + perp_axis * center_dist;
			axis2 = tan;
			axis1 = normal.cross(axis2).normalize();
			let perp_axis_mag = perp_axis.length();
			radius = (center_dist * perp_axis_mag).abs();
			if ulps_eq!(0.0, radius) {
				angle = 0.0;
				arc_len = 0.0;
			} else {
				let inv_radius = 1.0 / radius;
				let center_to_mid_dir = (start - center) + start_to_mid;
				let center_to_end_dir = (start - center) * inv_radius;
				let twist = perp_axis.dot(start_to_mid);
				angle = center_to_end_dir.dot(center_to_mid_dir).acos() * twist.signum();
				arc_len = radius * angle;
			}
		}

		Self{ center, axis1, axis2, radius, angle, arc_len }
	}

	pub fn arc_len(&self) -> f32 {
		self.arc_len
	}

	pub fn angle(&self) -> f32 {
		self.angle
	}

	pub fn set_angle(&mut self, angle: f32) {
		self.angle = angle;
	}
}

impl Curve for Arc {
    fn sample(&self, t: f32) -> (Vec3, Vec3) {
        let angle = self.angle * t;
		let sin_rot = angle.sin();
		let cos_rot = angle.cos();
		let pos = (self.center + self.axis1 * self.radius * cos_rot) + (self.axis2 * self.radius * sin_rot);
		let tan = (cos_rot * self.axis2) - (sin_rot * self.axis1);
		(pos, tan)
    }

    fn real_length(&self) -> f32 {
        self.arc_len
    }
}

impl AbsDiffEq for Arc {

    type Epsilon = f32;

    fn default_epsilon() -> Self::Epsilon {
        Self::Epsilon::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.center.abs_diff_eq(other.center, epsilon)
		&& self.axis1.abs_diff_eq(other.axis1, epsilon)
		&& self.axis2.abs_diff_eq(other.axis2, epsilon)
		&& self.radius.abs_diff_eq(&other.radius, epsilon)
		&& self.angle.abs_diff_eq(&other.angle, epsilon)
		&& self.arc_len.abs_diff_eq(&other.arc_len, epsilon)
    }
}

impl RelativeEq for Arc {
    fn default_max_relative() -> Self::Epsilon {
        Self::Epsilon::default_max_relative()
    }

    fn relative_eq(&self, other: &Self, epsilon: Self::Epsilon, max_relative: Self::Epsilon) -> bool {
        ApproxVec::from(self.center).relative_eq(&ApproxVec::from(other.center), epsilon, max_relative)
		&& ApproxVec::from(self.axis1).relative_eq(&ApproxVec::from(other.axis1), epsilon, max_relative)
		&& ApproxVec::from(self.axis2).relative_eq(&ApproxVec::from(other.axis2), epsilon, max_relative)
		&& self.radius.relative_eq(&other.radius, epsilon, max_relative)
		&& self.angle.relative_eq(&other.angle, epsilon, max_relative)
		&& self.arc_len.relative_eq(&other.arc_len, epsilon, max_relative)
    }
}

impl UlpsEq for Arc {
    fn default_max_ulps() -> u32 {
        Self::Epsilon::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        ApproxVec::from(self.center).ulps_eq(&ApproxVec::from(other.center), epsilon, max_ulps)
		&& ApproxVec::from(self.axis1).ulps_eq(&ApproxVec::from(other.axis1), epsilon, max_ulps)
		&& ApproxVec::from(self.axis2).ulps_eq(&ApproxVec::from(other.axis2), epsilon, max_ulps)
		&& self.radius.ulps_eq(&other.radius, epsilon, max_ulps)
		&& self.angle.ulps_eq(&other.angle, epsilon, max_ulps)
		&& self.arc_len.ulps_eq(&other.arc_len, epsilon, max_ulps)
    }
}


#[cfg(test)]
mod methods {

	use super::*;
	
	mod from_points {

		use super::*;
		use std::f32::consts::PI;
		use approx::assert_ulps_eq;
		use crate::vec::ApproxVec;

		#[test]
		fn straight_line() {
			let arc = Arc::from_points(-Vec3::X, Vec3::X, Vec3::X);
			assert_ulps_eq!(ApproxVec::ZERO, arc.center.into());
			assert_ulps_eq!(f32::INFINITY, arc.radius);
			assert_ulps_eq!(0.0, arc.angle);
			assert_ulps_eq!(2.0, arc.arc_len);
		}

		#[test]
		fn right_unit_quarter_turn() {
			let arc = Arc::from_points(Vec3::ZERO, -Vec3::Z, Vec3::new(1.0, 0.0, -1.0));
			assert_ulps_eq!(ApproxVec::X, arc.center.into());
			assert_ulps_eq!(-ApproxVec::X, arc.axis1.into());
			assert_ulps_eq!(-ApproxVec::Z, arc.axis2.into());
			assert_ulps_eq!(1.0, arc.radius);
			assert_ulps_eq!(PI / 2.0, arc.angle);
			assert_ulps_eq!(PI / 2.0, arc.arc_len);
		}

		#[test]
		fn upwards_unit_half_turn() {
			let arc = Arc::from_points(Vec3::ZERO, -Vec3::Z, Vec3::new(0.0, 1.0, -1.0));
			assert_ulps_eq!(ApproxVec::Y, arc.center.into());
			assert_ulps_eq!(-ApproxVec::Y, arc.axis1.into());
			assert_ulps_eq!(-ApproxVec::Z, arc.axis2.into());
			assert_ulps_eq!(1.0, arc.radius);
			assert_ulps_eq!(PI / 2.0, arc.angle);
			assert_ulps_eq!(PI / 2.0, arc.arc_len);
		}

		#[test]
		fn half_unit_turn_around_center() {
			let arc = Arc::from_points(-Vec3::X, -Vec3::Z, Vec3::X);
			assert_ulps_eq!(ApproxVec::ZERO, arc.center.into());
			assert_ulps_eq!(-ApproxVec::X, arc.axis1.into());
			assert_ulps_eq!(-ApproxVec::Z, arc.axis2.into());
			assert_ulps_eq!(1.0, arc.radius);
			assert_ulps_eq!(PI, arc.angle);
			assert_ulps_eq!(PI, arc.arc_len);
		}

		#[test]
		fn large_quarter_turn() {
			let arc = Arc::from_points(Vec3::new(5.0, 0.0, -20.0), -Vec3::Z, Vec3::new(-10.0, 0.0, -35.0));
			assert_ulps_eq!(ApproxVec::new(-10.0, 0.0, -20.0), arc.center.into());
			assert_ulps_eq!(ApproxVec::X, arc.axis1.into());
			assert_ulps_eq!(-ApproxVec::Z, arc.axis2.into());
			assert_ulps_eq!(15.0, arc.radius);
			assert_ulps_eq!(PI / 2.0, arc.angle);
			assert_ulps_eq!(PI / 2.0 * 15.0, arc.arc_len);
		}

		#[test]
		fn off_center_2() {
			let arc = Arc::from_points(Vec3::new(-2.0, 0.0, 2.0), -Vec3::Z, Vec3::new(-0.5, 0.0, 0.5));
			assert_ulps_eq!(ApproxVec::new(-0.5, 0.0, 2.0), arc.center.into());
			assert_ulps_eq!(-ApproxVec::X, arc.axis1.into());
			assert_ulps_eq!(-ApproxVec::Z, arc.axis2.into());
			assert_ulps_eq!(1.5, arc.radius);
			assert_ulps_eq!(PI / 2.0, arc.angle);
			assert_ulps_eq!(PI / 2.0 * 1.5, arc.arc_len);

			let arc = Arc::from_points(Vec3::new(-0.5, 0.0, 0.5), -Vec3::X, Vec3::new(-2.0, 0.0, 2.0));
			assert_ulps_eq!(ApproxVec::new(-0.5, 0.0, 2.0), arc.center.into());
			assert_ulps_eq!(-ApproxVec::Z, arc.axis1.into());
			assert_ulps_eq!(-ApproxVec::X, arc.axis2.into());
			assert_ulps_eq!(1.5, arc.radius);
			assert_ulps_eq!(PI / 2.0, arc.angle);
			assert_ulps_eq!(PI / 2.0 * 1.5, arc.arc_len);
		}

		#[test]
		fn off_center_1() {
			let arc = Arc::from_points(Vec3::new(1.0, 0.0, -1.0), Vec3::Z, Vec3::new(-0.5, 0.0, 0.5));
			assert_ulps_eq!(ApproxVec::new(-0.5, 0.0, -1.0), arc.center.into());
			assert_ulps_eq!(ApproxVec::X, arc.axis1.into());
			assert_ulps_eq!(ApproxVec::Z, arc.axis2.into());
			assert_ulps_eq!(1.5, arc.radius);
			assert_ulps_eq!(PI / 2.0, arc.angle);
			assert_ulps_eq!(PI / 2.0 * 1.5, arc.arc_len);

			let arc = Arc::from_points(Vec3::new(-0.5, 0.0, 0.5), Vec3::X, Vec3::new(1.0, 0.0, -1.0));
			assert_ulps_eq!(ApproxVec::new(-0.5, 0.0, -1.0), arc.center.into());
			assert_ulps_eq!(ApproxVec::Z, arc.axis1.into());
			assert_ulps_eq!(ApproxVec::X, arc.axis2.into());
			assert_ulps_eq!(1.5, arc.radius);
			assert_ulps_eq!(PI / 2.0, arc.angle);
			assert_ulps_eq!(PI / 2.0 * 1.5, arc.arc_len);
		}

		#[test]
		#[should_panic]
		fn tangent_not_normalized() {
			Arc::from_points(Vec3::ZERO, Vec3::ONE, Vec3::ZERO);
		}
	}
}

#[cfg(test)]
mod traits {

	use super::*;
	use std::f32::consts::PI;
	use approx::{assert_abs_diff_eq, assert_abs_diff_ne, assert_relative_eq, assert_relative_ne, assert_ulps_eq, assert_ulps_ne};

	mod curve {

		use super::*;

		mod sample {

			use super::*;
			use std::f32::consts::SQRT_2;
			use approx::assert_ulps_eq;
			use crate::vec::ApproxVec;

			#[test]
			fn basic() {
				let arc = Arc::from_points(-Vec3::X, -Vec3::Z, -Vec3::Z);
				let (pos, tan) = arc.sample(0.0);
				assert_ulps_eq!(ApproxVec::new(-1.0, 0.0, 0.0), pos.into());
				assert_ulps_eq!(-ApproxVec::Z, tan.into());
				let (pos, tan) = arc.sample(0.5);
				assert_ulps_eq!(ApproxVec::new(-SQRT_2 / 2.0, 0.0, -SQRT_2 / 2.0), pos.into());
				assert_ulps_eq!(ApproxVec::new(SQRT_2 / 2.0, 0.0, -SQRT_2 / 2.0), tan.into());
				let (pos, tan) = arc.sample(1.0);
				assert_ulps_eq!(ApproxVec::new(0.0, 0.0, -1.0), pos.into());
				assert_ulps_eq!(ApproxVec::X, tan.into());
			}

			#[test]
			fn basic_radius_2() {
				let arc = Arc::from_points(Vec3::new(-2.0, 0.0, 0.0), -Vec3::Z, Vec3::new(0.0, 0.0, -2.0));
				let (pos, tan) = arc.sample(0.0);
				assert_ulps_eq!(ApproxVec::new(-2.0, 0.0, 0.0), pos.into());
				assert_ulps_eq!(-ApproxVec::Z, tan.into());
				let (pos, tan) = arc.sample(0.5);
				assert_ulps_eq!(ApproxVec::new(-SQRT_2, 0.0, -SQRT_2), pos.into());
				assert_ulps_eq!(ApproxVec::new(SQRT_2 / 2.0, 0.0, -SQRT_2 / 2.0), tan.into());
				let (pos, tan) = arc.sample(1.0);
				assert_ulps_eq!(ApproxVec::new(0.0, 0.0, -2.0), pos.into());
				assert_ulps_eq!(ApproxVec::X, tan.into());
			}
		}
	}

	#[test]
	fn abs_diff_eq() {
		let a = Arc::from_points(Vec3::new(-2.0, 0.0, 0.0), -Vec3::Z, Vec3::new(0.0, 0.0, -2.0));
		let b = Arc::from_points(Vec3::new(2.0, 0.0, 0.0), -Vec3::Z, Vec3::new(0.0, 0.0, -2.0));
		let c = Arc::new(Vec3::ZERO, Vec3::X, -Vec3::Z, 2.0, PI / 2.0, PI);
		assert_abs_diff_ne!(a, b);
		assert_abs_diff_ne!(a, c);
		assert_abs_diff_eq!(b, c);
	}

	#[test]
	fn relative_eq() {
		let a = Arc::from_points(Vec3::new(-2.0, 0.0, 0.0), -Vec3::Z, Vec3::new(0.0, 0.0, -2.0));
		let b = Arc::from_points(Vec3::new(2.0, 0.0, 0.0), -Vec3::Z, Vec3::new(0.0, 0.0, -2.0));
		let c = Arc::new(Vec3::ZERO, Vec3::X, -Vec3::Z, 2.0, PI / 2.0, PI);
		assert_relative_ne!(a, b);
		assert_relative_ne!(a, c);
		assert_relative_eq!(b, c);
	}

	#[test]
	fn ulps_eq() {
		let a = Arc::from_points(Vec3::new(-2.0, 0.0, 0.0), -Vec3::Z, Vec3::new(0.0, 0.0, -2.0));
		let b = Arc::from_points(Vec3::new(2.0, 0.0, 0.0), -Vec3::Z, Vec3::new(0.0, 0.0, -2.0));
		let c = Arc::new(Vec3::ZERO, Vec3::X, -Vec3::Z, 2.0, PI / 2.0, PI);
		assert_ulps_ne!(a, b);
		assert_ulps_ne!(a, c);
		assert_ulps_eq!(b, c);
	}
}