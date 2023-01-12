//! Used to define a curve that you can interpolate smoothly from the start to end position


use std::ops::Div;

mod vec;


/// A curve that can be smoothly interpolated over
pub trait Curve<T> where T: Div<Output=T> + PartialEq {

	/// Samples the curve at the given fractional position, with `0.0` being the start of the curve
	/// and `1.0` being the end. Returns a tuple with the position and tangent vector of the 
	fn sample(&self, t: T) -> (vec::Vec3<T>, vec::Vec3<T>);

	/// Gets the length of the path traced by this curve in world units
	fn real_length(&self) -> T;

	/// Samples the curve at the given point by real distance along the curve
	///
	/// The returned value contains the position and tangent vector
	fn sample_dist(&self, dist: T) -> (vec::Vec3<T>, vec::Vec3<T>) {
		self.sample(dist / self.real_length())
	}
}



