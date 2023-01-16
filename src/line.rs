use bevy::prelude::Vec3;
use crate::Curve;



/// Straight line between two points, interpolated linearly
pub struct LinearCurve {
	p1: Vec3,
	p2: Vec3,
}

impl LinearCurve {

	/// Creates a new linear curve with the given endpoints
	pub fn new(p1: Vec3, p2: Vec3) -> Self {
		Self{ p1, p2 }
	}
}

impl Curve for LinearCurve {

    fn sample(&self, t: f32) -> (Vec3, Vec3) {
        let pos = self.p1.lerp(self.p2, t);
		let tan = (self.p2 - self.p1).normalize();
		(pos, tan)
    }

    fn real_length(&self) -> f32 {
        self.p1.distance(self.p2)
    }
}


#[cfg(test)]
mod traits {

	mod curve {

		use approx::ulps_eq;
		use bevy::prelude::Vec3;
		use crate::{Curve, line::LinearCurve, vec::ApproxVec};

		#[test]
		fn sample() {
			let line = LinearCurve::new(Vec3::NEG_ONE, Vec3::ONE);
			let expected_tan = Vec3::ONE.normalize();
			let test_cases = vec![
				(0.0, Vec3::new(-1.0, -1.0, -1.0)),
				(0.25, Vec3::new(-0.5, -0.5, -0.5)),
				(0.5, Vec3::new(0.0, 0.0, 0.0)),
				(0.75, Vec3::new(0.5, 0.5, 0.5)),
				(1.0, Vec3::new(1.0, 1.0, 1.0)),
			];
			for (t, expected_position) in test_cases.iter() {
				let (pos, tan) = line.sample(*t);
				assert!(
					ulps_eq!(ApproxVec::from(pos), ApproxVec::from(expected_position)),
					"Expected position of sample {} to be {}, but {} was returned",
					t, expected_position, pos,
				);
				assert!(
					ulps_eq!(ApproxVec::from(expected_tan), ApproxVec::from(tan)),
					"Expected tangent of sample at {} to be {}, but {} was returned",
					t, expected_tan, tan,
				)
			}
		}
	}
}
