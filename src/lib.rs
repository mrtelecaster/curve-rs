//! Used to define a curve that you can interpolate smoothly from the start to end position


use bevy::prelude::Vec3;

mod arc; pub use arc::Arc;
mod line; pub use line::LinearCurve;
mod vec; pub use vec::ApproxVec;


/// A curve that can be smoothly interpolated over
pub trait Curve {

	/// Samples the curve at the given fractional position, with `0.0` being the start of the curve
	/// and `1.0` being the end. Returns a tuple with the position and tangent vector of the 
	fn sample(&self, t: f32) -> (Vec3, Vec3);

	/// Gets the length of the path traced by this curve in world units
	fn real_length(&self) -> f32;

	/// Samples the curve at the given point by real distance along the curve
	///
	/// The returned value contains the position and tangent vector
	fn sample_dist(&self, dist: f32) -> (Vec3, Vec3) {
		self.sample(dist / self.real_length())
	}
}



