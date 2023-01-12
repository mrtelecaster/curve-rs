use approx::{AbsDiffEq, RelativeEq, UlpsEq};



/// A position in 3D space
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vec3<T> where T: PartialEq {
	/// X value of this vector
	pub x: T,
	/// Y value of this vector
	pub y: T,
	/// Z value of this vector
	pub z: T,
}

impl<T> Vec3<T> where T: PartialEq {
	/// 
	pub fn new(x: T, y: T, z: T) -> Self {
		Self{ x, y, z }
	}
}

impl<T> AbsDiffEq for Vec3<T> where T: AbsDiffEq<Epsilon=T> + Clone + Copy + PartialEq {

    type Epsilon = T;

    fn default_epsilon() -> Self::Epsilon {
        T::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.x.abs_diff_eq(&other.x, epsilon)
		&& self.y.abs_diff_eq(&other.y, epsilon)
		&& self.z.abs_diff_eq(&other.z, epsilon)
    }
}

impl<T> RelativeEq for Vec3<T> where T: AbsDiffEq<Epsilon=T> + Clone + Copy + RelativeEq {
    fn default_max_relative() -> Self::Epsilon {
        T::default_max_relative()
    }

    fn relative_eq(&self, other: &Self, epsilon: Self::Epsilon, max_relative: Self::Epsilon) -> bool {
        self.x.relative_eq(&other.x, epsilon, max_relative)
		&& self.y.relative_eq(&other.y, epsilon, max_relative)
		&& self.z.relative_eq(&other.z, epsilon, max_relative)
    }
}

impl<T> UlpsEq for Vec3<T> where T: AbsDiffEq<Epsilon=T> + Clone + Copy + UlpsEq {
    fn default_max_ulps() -> u32 {
        T::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        self.x.ulps_eq(&other.x, epsilon, max_ulps)
		&& self.y.ulps_eq(&other.y, epsilon, max_ulps)
		&& self.z.ulps_eq(&other.z, epsilon, max_ulps)
    }
}



#[cfg(test)]
mod tests {

	use super::*;

	mod approx_impl {

		use super::*;
		use approx::{assert_abs_diff_eq, assert_relative_eq, assert_ulps_eq};

		#[test]
		pub fn abs_diff_eq() {
			let a = Vec3::new(1.0, -2.0, 3.0);
			let b = Vec3::new(1.0, -2.0, 3.0);
			assert_abs_diff_eq!(a, b);
		}

		#[test]
		pub fn relative_eq() {
			let a = Vec3::new(1.0, -2.0, 3.0);
			let b = Vec3::new(1.0, -2.0, 3.0);
			assert_relative_eq!(a, b);
		}

		#[test]
		pub fn ulps_eq() {
			let a = Vec3::new(1.0, -2.0, 3.0);
			let b = Vec3::new(1.0, -2.0, 3.0);
			assert_ulps_eq!(a, b);
		}
	}
}
