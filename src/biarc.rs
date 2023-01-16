use std::f32::consts::{PI, TAU};
use approx::{AbsDiffEq};
use bevy::prelude::*;
use crate::{Curve, Arc};


#[derive(Debug, PartialEq)]
pub struct BiArc {
	arc_1: Arc,
	arc_2: Arc
}

impl BiArc {
	pub fn new(arc_1: Arc, arc_2: Arc) -> Self {
		Self{ arc_1, arc_2 }
	}

	pub fn from_endpoints(p1: Vec3, t1: Vec3, p2: Vec3, t2: Vec3) -> Self {
		let (arc_1, arc_2, _center) = Self::compute_arcs(p1, t1, p2, t2);
		Self{ arc_1, arc_2 }
	}

	pub fn compute_arcs(p1: Vec3, t1: Vec3, p2: Vec3, t2: Vec3) -> (Arc, Arc, Option<Vec3>) {
		assert!( t1.is_normalized() );
		assert!( t2.is_normalized() );

		let epsilon = 0.001;
		let v = p2 - p1;
		let v_dot_v = v.dot(v);

		// if the control points are equal, we don't need to interpolate
		if v_dot_v < epsilon
		{
			return (
				Arc::new(p1, v, v, 0.0, 0.0, 0.0),
				Arc::new(p1, v, v, 0.0, 0.0, 0.0),
				None,
			);
		}

		let mut p_arc_1: Arc;
		let mut p_arc_2: Arc;
		// compute the denominator for the quadratic formula
		let t = t1 + t2;
		let v_dot_t = v.dot(t);
		let t1_dot_t2 = t1.dot(t2);
		let denominator = 2.0 * (1.0 - t1_dot_t2);

		let d;

		if denominator < epsilon {
			let v_dot_t2 = v.dot(t2);
			if v_dot_t2.abs() < epsilon {
				let v_mag = v_dot_v.sqrt();
				let inv_v_mag_sqr = 1.0 / v_dot_v;

				// compute the normal to the plane containing the arcs
				// (this has length vMag)
				let plane_normal = v.cross(t2);

				// compute the axis perpendicular to the tangent direction and
				// aligned with the circles (this has length vMag*vMag)
				let perp_axis = plane_normal.cross(v);

				let radius = v_mag * 0.25;

				let center_to_p1 = v * -0.25;

				// interpolate across two semicircles
				p_arc_1 = Arc::new(
					p1 - center_to_p1,
					center_to_p1,
					perp_axis * radius * inv_v_mag_sqr,
					radius,
					PI,
					radius * PI,
				);

				p_arc_2 = Arc::new(
					p2 + center_to_p1,
					center_to_p1 * -1.0,
					perp_axis * -radius * inv_v_mag_sqr,
					radius,
					PI,
					radius * PI,
				);

				return (p_arc_1, p_arc_2, None);
			} else {
				d = v_dot_v / (4.0 * v_dot_t2);
			}
		} else {
			// use the positive result of the quadratic formula
			let discriminant = v_dot_t * v_dot_t + denominator * v_dot_v;
			d = (-v_dot_t + discriminant.sqrt()) / denominator;
		}

		// compute the connection point (i.e. the mid point)
		let mut pm = t1 - t2;
		pm = p2 + pm * d;
		pm = pm + p1;
		pm = pm * 0.5;

		// compute vectors from the end points to the mid point
		let p1_to_pm = pm - p1;
		let p2_to_pm = pm - p2;

		// compute the arcs
		p_arc_1 = Arc::from_points(p1, t1, p1 + p1_to_pm);
		p_arc_2 = Arc::from_points(p2, -t2, p2 + p2_to_pm);

		// use the longer path around the circle if d is negative
		if d < 0.0
		{
			p_arc_1.set_angle(p_arc_1.angle().signum() * TAU - p_arc_1.angle());
			p_arc_2.set_angle(p_arc_2.angle().signum() * TAU - p_arc_2.angle());
		}

		(p_arc_1, p_arc_2, Some(pm))
	}
}

impl Curve for BiArc {
	fn sample(&self, t: f32) -> (Vec3, Vec3) {
		let total_len = self.real_length();
		let arc_1_frac_len = self.arc_1.real_length() / total_len;
		let arc_2_frac_len = self.arc_2.real_length() / total_len;
		if t >= arc_1_frac_len {
			return self.arc_2.sample((t - arc_1_frac_len) / arc_2_frac_len);
		} else {
			return self.arc_1.sample(t / arc_1_frac_len);
		}
	}

	fn real_length(&self) -> f32 {
		self.arc_1.arc_len() + self.arc_2.arc_len()
	}
}

impl AbsDiffEq for BiArc {

	type Epsilon = f32;

	fn default_epsilon() -> Self::Epsilon {
		Self::Epsilon::default_epsilon()
	}

	fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
		self.arc_1.abs_diff_eq(&other.arc_1, epsilon)
		&& self.arc_2.abs_diff_eq(&other.arc_2, epsilon)
	}
}


#[cfg(test)]
mod tests {

	use super::*;
	use std::f32::consts::{FRAC_1_SQRT_2, SQRT_2};
	use approx::{abs_diff_eq, abs_diff_ne};
	use crate::ApproxVec;

	mod methods {

		use super::*;

		mod from_endpoints {

			use super::*;

			#[test]
			fn basic() {
				let expected_arc_1 = Arc::from_points(Vec3::new(-1.0, 0.0, -1.0), Vec3::Z, Vec3::ZERO);
				let expected_arc_2 = Arc::from_points(Vec3::new(1.0, 0.0, 1.0), -Vec3::Z, Vec3::ZERO);
				let expected_biarc = BiArc::new(expected_arc_1, expected_arc_2);
				let p1 = Vec3::new(-1.0, 0.0, -1.0);
				let t1 = Vec3::Z;
				let p2 = Vec3::new(1.0, 0.0, 1.0);
				let t2 = Vec3::Z;
				let result_biarc = BiArc::from_endpoints(p1, t1, p2, t2);
				assert!(
					abs_diff_eq!(expected_biarc, result_biarc),
					"Expected endpoints ({}, {}) and ({}, {}) to result in biarc {:?}, but {:?} was returned",
					p1, t1, p2, t2, expected_biarc, result_biarc
				)
			}

			#[test]
			fn asymmetric() {
				let expected_arc_1 = Arc::from_points(Vec3::new(-2.0, 0.0, 2.0), -Vec3::Z, Vec3::new(-0.5, 0.0, 0.5));
				let expected_arc_2 = Arc::from_points(Vec3::new(1.0, 0.0, -1.0), Vec3::Z, Vec3::new(-0.5, 0.0, 0.5));
				let expected_biarc = BiArc::new(expected_arc_1, expected_arc_2);
				let p1 = Vec3::new(-2.0, 0.0, 2.0);
				let t1 = -Vec3::Z;
				let p2 = Vec3::new(1.0, 0.0, -1.0);
				let t2 = -Vec3::Z;
				let result_biarc = BiArc::from_endpoints(p1, t1, p2, t2);
				assert!(
					abs_diff_eq!(expected_biarc, result_biarc),
					"Expected endpoints ({}, {}) and ({}, {}) to result in biarc {:?}, but {:?} was returned",
					p1, t1, p2, t2, expected_biarc, result_biarc
				)
			}
		}

		mod compute_arcs {

			use super::*;

			#[test]
			fn basic() {
				let expected_arc_1 = Arc::from_points(Vec3::new(-1.0, 0.0, -1.0), Vec3::Z, Vec3::ZERO);
				let expected_arc_2 = Arc::from_points(Vec3::new(1.0, 0.0, 1.0), -Vec3::Z, Vec3::ZERO);
				let expected_center = Vec3::ZERO;
				let p1 = Vec3::new(-1.0, 0.0, -1.0);
				let t1 = Vec3::Z;
				let p2 = Vec3::new(1.0, 0.0, 1.0);
				let t2 = Vec3::Z;
				let result = BiArc::compute_arcs(p1, t1, p2, t2);
				let (result_arc_1, result_arc_2, result_center_opt) = result;
				assert!(
					abs_diff_eq!(expected_arc_1, result_arc_1),
					"Expected endpoints ({}, {}) and ({}, {}) to result in arc 1 parameters {:?}, but {:?} was returned",
					p1, t1, p2, t2, expected_arc_1, result_arc_1,
				);
				assert!(
					abs_diff_eq!(expected_arc_2, result_arc_2),
					"Expected endpoints ({}, {}) and ({}, {}) to result in arc 2 parameters {:?}, but {:?} was returned",
					p1, t1, p2, t2, expected_arc_2, result_arc_2,
				);
				assert!(result_center_opt.is_some());
				let result_center = result_center_opt.unwrap();
				assert!(
					abs_diff_eq!(ApproxVec::from(expected_center), ApproxVec::from(result_center)),
					"Expected endpoints ({} {}) and ({} {}) to result in an arc join point at {}, but {} was returned",
					p1, t1, p2, t2, expected_center, result_center,
				);
			}

			#[test]
			fn asymmetric() {
				let expected_arc_1 = Arc::from_points(Vec3::new(-2.0, 0.0, 2.0), -Vec3::Z, Vec3::new(-0.5, 0.0, 0.5));
				let expected_arc_2 = Arc::from_points(Vec3::new(1.0, 0.0, -1.0), Vec3::Z, Vec3::new(-0.5, 0.0, 0.5));
				let expected_center = Vec3::new(-0.5, 0.0, 0.5);
				let p1 = Vec3::new(-2.0, 0.0, 2.0);
				let t1 = -Vec3::Z;
				let p2 = Vec3::new(1.0, 0.0, -1.0);
				let t2 = -Vec3::Z;
				let result = BiArc::compute_arcs(p1, t1, p2, t2);
				let (result_arc_1, result_arc_2, result_center_opt) = result;
				assert!(
					abs_diff_eq!(expected_arc_1, result_arc_1),
					"Expected endpoints ({}, {}) and ({}, {}) to result in arc 1 parameters {:?}, but {:?} was returned",
					p1, t1, p2, t2, expected_arc_1, result_arc_1,
				);
				assert!(
					abs_diff_eq!(expected_arc_2, result_arc_2),
					"Expected endpoints ({}, {}) and ({}, {}) to result in arc 2 parameters {:?}, but {:?} was returned",
					p1, t1, p2, t2, expected_arc_2, result_arc_2,
				);
				assert!(result_center_opt.is_some());
				let result_center = result_center_opt.unwrap();
				assert!(
					abs_diff_eq!(ApproxVec::from(expected_center), ApproxVec::from(result_center)),
					"Expected endpoints ({} {}) and ({} {}) to result in an arc join point at {}, but {} was returned",
					p1, t1, p2, t2, expected_center, result_center,
				);
			}
		}
	}

	mod traits {

		use super::*;

		mod curve {

			use super::*;
			
			mod sample {

				use super::*;

				fn sample_test_case(arc: BiArc, test_cases: Vec<(f32, (Vec3, Vec3))>) -> Result<(), String> {
					for (t, case) in test_cases.iter() {
						let (expected_pos, expected_tan) = case;
						let (result_pos, result_tan) = arc.sample(*t);
						if abs_diff_ne!(ApproxVec::from(expected_pos), ApproxVec::from(result_pos)) {
							return Err(format!(
								"Expected sample at t={} to return position {}, but {} was returned",
								t, expected_pos, result_pos
							));
						}
						if abs_diff_ne!(ApproxVec::from(expected_tan), ApproxVec::from(result_tan)) {
							return Err(format!(
								"Expected sample at t={} to return tangent vector {}, but {} was returned",
								t, expected_tan, result_tan
							));
						}
					}
					Ok(())
				}

				#[test]
				fn basic() {
					let arc_1 = Arc::from_points(Vec3::new(-1.0, 0.0, -1.0), Vec3::Z, Vec3::ZERO);
					let arc_2 = Arc::from_points(Vec3::ZERO, Vec3::X, Vec3::new(1.0, 0.0, 1.0));
					let biarc = BiArc::new(arc_1, arc_2);
					let test_cases = vec![
						(0.0, (Vec3::new(-1.0, 0.0, -1.0), Vec3::Z)),
						(0.25, (Vec3::new(-FRAC_1_SQRT_2, 0.0, FRAC_1_SQRT_2 - 1.0), Vec3::new(FRAC_1_SQRT_2, 0.0, FRAC_1_SQRT_2))),
						(0.5, (Vec3::ZERO, Vec3::X)),
						(0.75, (Vec3::new(FRAC_1_SQRT_2, 0.0, 1.0 - FRAC_1_SQRT_2), Vec3::new(FRAC_1_SQRT_2, 0.0, FRAC_1_SQRT_2))),
						(1.0, (Vec3::new(1.0, 0.0, 1.0), Vec3::Z)),
					];
					sample_test_case(biarc, test_cases).unwrap();
				}

				#[test]
				fn asymmetric() {
					let arc_1 = Arc::from_points(Vec3::new(-2.0, 0.0, 2.0), -Vec3::Z, Vec3::ZERO);
					let arc_2 = Arc::from_points(Vec3::ZERO, Vec3::X, Vec3::new(1.0, 0.0, -1.0));
					let biarc = BiArc::new(arc_1, arc_2);
					let test_cases = vec![
						(0.0, (Vec3::new(-2.0, 0.0, 2.0), -Vec3::Z)),
						(1.0/3.0, (Vec3::new(-SQRT_2, 0.0, 2.0 - SQRT_2), Vec3::new(FRAC_1_SQRT_2, 0.0, -FRAC_1_SQRT_2))),
						(2.0/3.0, (Vec3::ZERO, Vec3::X)),
						(2.5/3.0, (Vec3::new(FRAC_1_SQRT_2, 0.0, FRAC_1_SQRT_2 - 1.0), Vec3::new(FRAC_1_SQRT_2, 0.0, -FRAC_1_SQRT_2))),
						(1.0, (Vec3::new(1.0, 0.0, -1.0), -Vec3::Z)),
					];
					sample_test_case(biarc, test_cases).unwrap();
				}
			}
		}

		mod abs_diff_eq {

			use super::*;

			#[test]
			#[ignore]
			fn basic() {
				let arc_1 = Arc::from_points(Vec3::NEG_ONE, Vec3::Z, Vec3::ZERO);
				let arc_2 = Arc::from_points(Vec3::ZERO, Vec3::X, Vec3::ONE);
				let _biarc = BiArc::new(arc_1, arc_2);
			}
		}
	}
}
