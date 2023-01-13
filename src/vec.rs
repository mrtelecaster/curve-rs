use std::ops::Neg;
use bevy::prelude::Vec3;
use approx::{AbsDiffEq, RelativeEq, UlpsEq};


/// Wrapper type around a [Bevy engine `Vec3`](Vec3) that supports approximate comparison using
/// [the `approx` crate](approx)
/// 
/// Mostly used for unit testing
#[derive(Debug, PartialEq)]
pub struct ApproxVec { vec: Vec3 }

impl ApproxVec {

	pub const ZERO: Self = Self::splat(0.0);

	pub const X: Self = Self::new(1.0, 0.0, 0.0);

	pub const Y: Self = Self::new(0.0, 1.0, 0.0);

	pub const Z: Self = Self::new(0.0, 0.0, 1.0);

	/// Creates a new vector with the given X, Y, and Z values
	#[inline]
	pub const fn new(x: f32, y: f32, z: f32) -> Self {
		Self{ vec: Vec3::new(x, y, z) }
	}

	#[inline]
	pub const fn splat(val: f32) -> Self {
		Self{ vec: Vec3::splat(val) }
	}
}

impl From<Vec3> for ApproxVec {
    fn from(vec: Vec3) -> Self {
        Self{ vec }
    }
}

impl From<&Vec3> for ApproxVec {
    fn from(v: &Vec3) -> Self {
        Self{ vec: v.clone() }
    }
}

impl Neg for ApproxVec {

    type Output = Self;

    fn neg(self) -> Self::Output {
        Self{ vec: -self.vec }
    }
}

impl AbsDiffEq for ApproxVec {

    type Epsilon = f32;

    fn default_epsilon() -> Self::Epsilon {
        f32::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.vec.x.abs_diff_eq(&other.vec.x, epsilon)
		&& self.vec.y.abs_diff_eq(&other.vec.y, epsilon)
		&& self.vec.z.abs_diff_eq(&other.vec.z, epsilon)
    }
}

impl RelativeEq for ApproxVec {

    fn default_max_relative() -> Self::Epsilon {
        f32::default_max_relative()
    }

    fn relative_eq(&self, other: &Self, epsilon: Self::Epsilon, max_relative: Self::Epsilon) -> bool {
        self.vec.x.relative_eq(&other.vec.x, epsilon, max_relative)
		&& self.vec.y.relative_eq(&other.vec.y, epsilon, max_relative)
		&& self.vec.z.relative_eq(&other.vec.z, epsilon, max_relative)
    }
}

impl UlpsEq for ApproxVec {

    fn default_max_ulps() -> u32 {
        f32::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        self.vec.x.ulps_eq(&other.vec.x, epsilon, max_ulps)
		&& self.vec.y.ulps_eq(&other.vec.y, epsilon, max_ulps)
		&& self.vec.z.ulps_eq(&other.vec.z, epsilon, max_ulps)
    }
}


mod tests {

	mod traits {

		mod abs_diff_eq {

			use crate::vec::ApproxVec;
			use approx::{assert_abs_diff_eq, assert_abs_diff_ne};

			/// Basic test case demonstrating intended usage
			#[test]
			fn basic() {
				assert_abs_diff_eq!(ApproxVec::new(0.0, 0.0, 0.0), ApproxVec::new(0.0, 0.0, 0.0));
				assert_abs_diff_eq!(ApproxVec::new(1.0, 2.0, 3.0), ApproxVec::new(1.0, 2.0, 3.0));
				assert_abs_diff_eq!(ApproxVec::new(-1.0, 2.0, -3.0), ApproxVec::new(-1.0, 2.0, -3.0));
				assert_abs_diff_eq!(ApproxVec::new(1.0, -2.0, 3.0), ApproxVec::new(1.0, -2.0, 3.0));

				assert_abs_diff_ne!(ApproxVec::new(1.0, 2.0, 3.0), ApproxVec::new(2.0, 3.0, 4.0));
				assert_abs_diff_ne!(ApproxVec::new(1.0, 2.0, 3.0), ApproxVec::new(0.0, 2.0, 3.0));
				assert_abs_diff_ne!(ApproxVec::new(1.0, 2.0, 3.0), ApproxVec::new(1.0, 0.0, 3.0));
				assert_abs_diff_ne!(ApproxVec::new(1.0, 2.0, 3.0), ApproxVec::new(1.0, 2.0, 0.0));
			}

			/// Check equality of values that are close to zero without being approximately equal
			#[test]
			fn zero() {
				assert_abs_diff_eq!(ApproxVec::new(0.0, 0.0, 0.0), ApproxVec::new(0.0, 0.0, 0.0));
				assert_abs_diff_eq!(ApproxVec::new(-0.0, 0.0, -0.0), ApproxVec::new(0.0, -0.0, -0.0));

				assert_abs_diff_ne!(ApproxVec::new(0.000001, 0.0, 0.0), ApproxVec::new(0.0, 0.0, 0.0));
				assert_abs_diff_ne!(ApproxVec::new(0.0, 0.0, 0.0), ApproxVec::new(0.0, 0.000001, 0.0));
				assert_abs_diff_ne!(ApproxVec::new(0.0, 0.0, -0.000001), ApproxVec::new(0.0, 0.0, 0.0));
				assert_abs_diff_ne!(ApproxVec::new(0.0, 0.0, 0.0), ApproxVec::new(-0.000001, 0.0, 0.0));
			}

			/// Test approximate equality of large numbers
			#[test]
			fn big() {
				assert_abs_diff_eq!(ApproxVec::new(100000000.0, 0.0, 0.0), ApproxVec::new(100000001.0, 0.0, 0.0));
				assert_abs_diff_eq!(ApproxVec::new(0.0, 100000001.0, 0.0), ApproxVec::new(0.0, 100000000.0, 0.0));
				assert_abs_diff_eq!(ApproxVec::new(0.0, 0.0, -100000000.0), ApproxVec::new(0.0, 0.0, -100000001.0));
				assert_abs_diff_eq!(ApproxVec::new(-100000001.0, 0.0, 0.0), ApproxVec::new(-100000000.0, 0.0, 0.0));

				assert_abs_diff_ne!(ApproxVec::new(0.0, 10000.0, 0.0), ApproxVec::new(0.0, 10001.0, 0.0));
				assert_abs_diff_ne!(ApproxVec::new(0.0, 0.0, 10001.0), ApproxVec::new(0.0, 0.0, 10000.0));
				assert_abs_diff_ne!(ApproxVec::new(-10000.0, 0.0, 0.0), ApproxVec::new(-10001.0, 0.0, 0.0));
				assert_abs_diff_ne!(ApproxVec::new(0.0, -10001.0, 0.0), ApproxVec::new(0.0, -10000.0, 0.0));
			}

			/// Test approximate equality of numbers around 1
			#[test]
			fn med() {
				assert_abs_diff_eq!(ApproxVec::new(1.0000001, 0.0, 0.0), ApproxVec::new(1.0000002, 0.0, 0.0));
				assert_abs_diff_eq!(ApproxVec::new(0.0, 1.0000002, 0.0), ApproxVec::new(0.0, 1.0000001, 0.0));
				assert_abs_diff_eq!(ApproxVec::new(0.0, 0.0, -1.0000001), ApproxVec::new(0.0, 0.0, -1.0000002));
				assert_abs_diff_eq!(ApproxVec::new(-1.0000002, 0.0, 0.0), ApproxVec::new(-1.0000001, 0.0, 0.0));

				assert_abs_diff_ne!(ApproxVec::new(0.0, 1.000001, 0.0), ApproxVec::new(0.0, 1.000002, 0.0));
				assert_abs_diff_ne!(ApproxVec::new(0.0, 0.0, 1.000002), ApproxVec::new(0.0, 0.0, 1.000001));
				assert_abs_diff_ne!(ApproxVec::new(-1.000001, 0.0, 0.0), ApproxVec::new(-1.000002, 0.0, 0.0));
				assert_abs_diff_ne!(ApproxVec::new(0.0, -1.000002, 0.0), ApproxVec::new(0.0, -1.000001, 0.0));
			}

			/// Test approximate equality of small numbers close to zero
			#[test]
			fn small() {
				assert_abs_diff_eq!(ApproxVec::new(0.0, 0.0,0.000010001), ApproxVec::new(0.0, 0.0, 0.000010002));
				assert_abs_diff_eq!(ApproxVec::new(0.000010002, 0.0, 0.0), ApproxVec::new(0.000010001, 0.0, 0.0));
				assert_abs_diff_eq!(ApproxVec::new(0.0, -0.000010001, 0.0), ApproxVec::new(0.0, -0.000010002, 0.0));
				assert_abs_diff_eq!(ApproxVec::new(0.0, 0.0, -0.000010002), ApproxVec::new(0.0, 0.0, -0.000010001));

				assert_abs_diff_ne!(ApproxVec::new(0.000001002, 0.0, 0.0), ApproxVec::new(0.0000001001, 0.0, 0.0));
				assert_abs_diff_ne!(ApproxVec::new(0.0, 0.000001001, 0.0), ApproxVec::new(0.0, 0.0000001002, 0.0));
				assert_abs_diff_ne!(ApproxVec::new(0.0, 0.0, -0.000001002), ApproxVec::new(0.0, 0.0, -0.0000001001));
				assert_abs_diff_ne!(ApproxVec::new(-0.000001001, 0.0, 0.0), ApproxVec::new(-0.0000001002, 0.0, 0.0));
			}
		}

		mod relative_eq {

			use approx::{assert_relative_eq, assert_relative_ne};
			use crate::vec::ApproxVec;

			/// Demonstrate typical usage of [`ApproxVec::relative_eq`]
			#[test]
			fn basic() {
				assert_relative_eq!(ApproxVec::new(1.0, 2.0, 3.0), ApproxVec::new(1.0, 2.0, 3.0));
				assert_relative_eq!(ApproxVec::new(1.0, 2.0, 3.0), ApproxVec::new(1.0000001, 2.0000001, 3.0000001));

				assert_relative_ne!(ApproxVec::new(1.0, 2.0, 3.0), ApproxVec::new(1.000001, 2.000001, 3.000001));
			}

			/// Test values near zero
			#[test]
			fn zero() {
				assert_relative_eq!(ApproxVec::new(0.0, 0.0, 0.0), ApproxVec::new(0.0, 0.0, 0.0));
				assert_relative_eq!(ApproxVec::new(0.0, 0.0, 0.0), ApproxVec::new(0.0, -0.0, 0.0));
				assert_relative_eq!(ApproxVec::new(0.0, 0.0, -0.0), ApproxVec::new(0.0, 0.0, 0.0));
				assert_relative_eq!(ApproxVec::new(-0.0, 0.0, 0.0), ApproxVec::new(-0.0, 0.0, 0.0));

				assert_relative_ne!(ApproxVec::new(0.0, 0.000001, 0.0), ApproxVec::new(0.0, 0.0, 0.0));
				assert_relative_ne!(ApproxVec::new(0.0, 0.0, 0.0), ApproxVec::new(0.0, 0.0, 0.000001));
				assert_relative_ne!(ApproxVec::new(-0.000001, 0.0, 0.0), ApproxVec::new(0.0, 0.0, 0.0));
				assert_relative_ne!(ApproxVec::new(0.0, 0.0, 0.0), ApproxVec::new(0.0, -0.000001, 0.0));
			}
		}

		mod ulps_eq {

			use approx::{assert_ulps_eq, assert_ulps_ne};
			use crate::vec::ApproxVec;

			/// Demonstrate typical usage of [`ApproxVec::ulps_eq`] in unit test assertions
			#[test]
			fn basic() {
				assert_ulps_eq!(ApproxVec::new(1.0, 2.0, 3.0), ApproxVec::new(1.0, 2.0, 3.0));
				assert_ulps_eq!(ApproxVec::new(1.0, 2.0, 3.0), ApproxVec::new(1.0000001, 2.0000001, 3.0000001));

				assert_ulps_ne!(ApproxVec::new(1.0, 2.0, 3.0), ApproxVec::new(1.000001, 2.000001, 3.000001));
			}

			/// Test values near zero
			#[test]
			fn zero() {
				assert_ulps_eq!(ApproxVec::new(0.0, 0.0, 0.0), ApproxVec::new(0.0, 0.0, 0.0));
				assert_ulps_eq!(ApproxVec::new(0.0, 0.0, 0.0), ApproxVec::new(0.0, -0.0, 0.0));
				assert_ulps_eq!(ApproxVec::new(0.0, 0.0, -0.0), ApproxVec::new(0.0, 0.0, 0.0));
				assert_ulps_eq!(ApproxVec::new(-0.0, 0.0, 0.0), ApproxVec::new(-0.0, 0.0, 0.0));

				assert_ulps_ne!(ApproxVec::new(0.0, 0.000001, 0.0), ApproxVec::new(0.0, 0.0, 0.0));
				assert_ulps_ne!(ApproxVec::new(0.0, 0.0, 0.0), ApproxVec::new(0.0, 0.0, 0.000001));
				assert_ulps_ne!(ApproxVec::new(-0.000001, 0.0, 0.0), ApproxVec::new(0.0, 0.0, 0.0));
				assert_ulps_ne!(ApproxVec::new(0.0, 0.0, 0.0), ApproxVec::new(0.0, -0.000001, 0.0));
			}

			#[test]
			fn infinity() {
				assert_ulps_eq!(ApproxVec::new(f32::INFINITY, 0.0, 0.0), ApproxVec::new(f32::INFINITY, 0.0, 0.0));
				assert_ulps_eq!(ApproxVec::new(0.0, f32::NEG_INFINITY, 0.0), ApproxVec::new(0.0, f32::NEG_INFINITY, 0.0));
				assert_ulps_eq!(ApproxVec::new(0.0, 0.0, f32::INFINITY), ApproxVec::new(0.0, 0.0, f32::MAX));
				assert_ulps_eq!(ApproxVec::new(f32::NEG_INFINITY, 0.0, 0.0), ApproxVec::new(-f32::MAX, 0.0, 0.0));

				assert_ulps_ne!(ApproxVec::new(0.0, f32::NEG_INFINITY, 0.0), ApproxVec::new(0.0, f32::INFINITY, 0.0));
			}
		}
	}
}
