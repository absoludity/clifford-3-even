/// Represents a rotor in the even sub-algebra of Cl(3).
///
/// A rotor is composed of a scalar part and three bivector components:
/// r = b₀ + b₁e₂₃ + b₂e₃₁ + b₃e₁₂
/// where b₀ is the scalar and b₁, b₂, b₃ are the bivector components.
use std::fmt;
use std::ops::{Add, Div, Mul, Sub};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rotor {
    /// Scalar component
    pub scalar: f64,
    /// e₂₃ bivector component
    pub e23: f64,
    /// e₃₁ bivector component
    pub e31: f64,
    /// e₁₂ bivector component
    pub e12: f64,
}

impl fmt::Display for Rotor {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Start with the scalar part
        write!(f, "{}", self.scalar)?;

        // Add each bivector component with its corresponding label
        // Only show components that aren't zero (or very close to zero)
        if self.e23.abs() > 1e-10 {
            write!(
                f,
                "{}{}Ix",
                if self.e23 >= 0.0 { " + " } else { " - " },
                self.e23.abs()
            )?;
        }

        if self.e31.abs() > 1e-10 {
            write!(
                f,
                "{}{}Iy",
                if self.e31 >= 0.0 { " + " } else { " - " },
                self.e31.abs()
            )?;
        }

        if self.e12.abs() > 1e-10 {
            write!(
                f,
                "{}{}Iz",
                if self.e12 >= 0.0 { " + " } else { " - " },
                self.e12.abs()
            )?;
        }

        Ok(())
    }
}

impl Add for Rotor {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        // Addition is component-wise
        Self::new(
            self.scalar + rhs.scalar,
            self.e23 + rhs.e23,
            self.e31 + rhs.e31,
            self.e12 + rhs.e12,
        )
    }
}

impl Sub for Rotor {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        // Subtraction is component-wise
        Self::new(
            self.scalar - rhs.scalar,
            self.e23 - rhs.e23,
            self.e31 - rhs.e31,
            self.e12 - rhs.e12,
        )
    }
}

impl Mul for Rotor {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        // Rotor multiplication formula:
        // (a + b₁e₂₃ + b₂e₃₁ + b₃e₁₂) * (c + d₁e₂₃ + d₂e₃₁ + d₃e₁₂)
        
        // Scalar part: a*c - b₁*d₁ - b₂*d₂ - b₃*d₃
        let scalar = self.e12.mul_add(-rhs.e12, self.e31.mul_add(-rhs.e31, self.scalar.mul_add(rhs.scalar, -(self.e23 * rhs.e23))));
        
        // e₂₃ part: a*d₁ + b₁*c + b₂*d₃ - b₃*d₂
        let e23 = self.e12.mul_add(-rhs.e31, self.e31.mul_add(rhs.e12, self.scalar.mul_add(rhs.e23, self.e23 * rhs.scalar)));
        
        // e₃₁ part: a*d₂ + b₂*c + b₃*d₁ - b₁*d₃
        let e31 = self.e23.mul_add(-rhs.e12, self.e12.mul_add(rhs.e23, self.scalar.mul_add(rhs.e31, self.e31 * rhs.scalar)));
        
        // e₁₂ part: a*d₃ + b₃*c + b₁*d₂ - b₂*d₁
        let e12 = self.e31.mul_add(-rhs.e23, self.e23.mul_add(rhs.e31, self.scalar.mul_add(rhs.e12, self.e12 * rhs.scalar)));
        
        Self::new(scalar, e23, e31, e12)
    }
}

impl Div for Rotor {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        // Division by a rotor r is multiplication by its inverse r⁻¹
        // For a normalized rotor, the inverse is the reverse
        // For a general rotor, we need to divide by the norm squared
        
        // Calculate the norm squared
        let norm_squared = rhs.e12.mul_add(rhs.e12, rhs.e31.mul_add(rhs.e31, rhs.scalar.mul_add(rhs.scalar, rhs.e23 * rhs.e23)));
        
        // Calculate the inverse (reverse divided by norm squared)
        let inverse = Self::new(
            rhs.scalar / norm_squared,
            -rhs.e23 / norm_squared,
            -rhs.e31 / norm_squared,
            -rhs.e12 / norm_squared,
        );
        
        // Multiply by the inverse
        self * inverse
    }
}

impl Rotor {
    /// Creates a new rotor with the given components.
    #[must_use] pub const fn new(scalar: f64, e23: f64, e31: f64, e12: f64) -> Self {
        Self {
            scalar,
            e23,
            e31,
            e12,
        }
    }

    /// Creates an identity rotor (scalar = 1, bivectors = 0).
    #[must_use] pub const fn identity() -> Self {
        Self::new(1.0, 0.0, 0.0, 0.0)
    }

    /// Creates a rotor representing a rotation around the specified axis by the given angle.
    ///
    /// # Arguments
    ///
    /// * `axis` - A 3D vector [x, y, z] representing the axis of rotation
    /// * `angle` - The angle of rotation in radians
    ///            Note: When applying this rotor using the sandwich product (r * v * `r.reverse()`),
    ///            the resulting rotation will be by 2*angle.
    ///
    /// # Returns
    ///
    /// A rotor representing the specified rotation
    #[must_use] pub fn from_axis_angle(axis: [f64; 3], angle: f64) -> Self {
        // Normalize the axis vector
        let magnitude = axis[2].mul_add(axis[2], axis[0].mul_add(axis[0], axis[1] * axis[1])).sqrt();
        
        // Handle zero magnitude case
        if magnitude < 1e-10 {
            return Self::identity();
        }
        
        let normalized_axis = [
            axis[0] / magnitude,
            axis[1] / magnitude,
            axis[2] / magnitude,
        ];
        
        // Use the full angle (not half) for the rotor
        let cos_angle = angle.cos();
        let sin_angle = angle.sin();
        
        // Calculate bivector components (scaled by sin_angle)
        let e23 = normalized_axis[0] * sin_angle;
        let e31 = normalized_axis[1] * sin_angle;
        let e12 = normalized_axis[2] * sin_angle;
        
        Self::new(cos_angle, e23, e31, e12)
    }

    /// Normalizes the rotor to ensure it represents a proper rotation.
    #[must_use] pub fn normalize(&self) -> Self {
        let norm_squared = self.e12.mul_add(self.e12, self.e31.mul_add(self.e31, self.scalar.mul_add(self.scalar, self.e23 * self.e23)));

        let norm = norm_squared.sqrt();

        Self {
            scalar: self.scalar / norm,
            e23: self.e23 / norm,
            e31: self.e31 / norm,
            e12: self.e12 / norm,
        }
    }

    /// Returns the reverse of this rotor.
    ///
    /// The reverse of a rotor flips the sign of all bivector components
    /// while keeping the scalar component unchanged.
    #[must_use] pub fn reverse(&self) -> Self {
        Self {
            scalar: self.scalar,
            e23: -self.e23,
            e31: -self.e31,
            e12: -self.e12,
        }
    }
    
    /// Checks if this rotor is approximately equal to another rotor.
    ///
    /// # Arguments
    ///
    /// * `other` - The rotor to compare with
    /// * `epsilon` - The maximum allowed difference between components
    ///
    /// # Returns
    ///
    /// `true` if all components are within `epsilon` of each other
    #[must_use] pub fn approx_eq(&self, other: &Self, epsilon: f64) -> bool {
        (self.scalar - other.scalar).abs() < epsilon &&
        (self.e23 - other.e23).abs() < epsilon &&
        (self.e31 - other.e31).abs() < epsilon &&
        (self.e12 - other.e12).abs() < epsilon
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case(Rotor::new(1.0, 3.0, 2.0, 1.0), "1 + 3Ix + 2Iy + 1Iz")]
    #[case(Rotor::new(1.0, -3.0, 2.0, -1.0), "1 - 3Ix + 2Iy - 1Iz")]
    #[case(Rotor::new(1.0, 0.0, 0.0, 0.0), "1")]
    fn test_display(#[case] rotor: Rotor, #[case] expected: &str) {
        assert_eq!(format!("{}", rotor), expected);
    }

    #[rstest]
    #[case(Rotor::new(1.0, 1.0, 2.0, 3.0), Rotor::new(1.0, -1.0, -2.0, -3.0))]
    #[case(Rotor::new(2.0, -1.0, 0.0, 3.0), Rotor::new(2.0, 1.0, 0.0, -3.0))]
    #[case(Rotor::new(1.0, 0.0, 0.0, 0.0), Rotor::new(1.0, 0.0, 0.0, 0.0))]
    fn test_reverse(#[case] rotor: Rotor, #[case] expected: Rotor) {
        assert_eq!(rotor.reverse(), expected);
    }
    
    #[rstest]
    // Rotation around x-axis by π/4 (will rotate by π/2 when applied)
    #[case([1.0, 0.0, 0.0], std::f64::consts::FRAC_PI_4, 
           Rotor::new(std::f64::consts::FRAC_1_SQRT_2, std::f64::consts::FRAC_1_SQRT_2, 0.0, 0.0))]
    // Rotation around y-axis by π/4 (will rotate by π/2 when applied)
    #[case([0.0, 1.0, 0.0], std::f64::consts::FRAC_PI_4, 
           Rotor::new(std::f64::consts::FRAC_1_SQRT_2, 0.0, std::f64::consts::FRAC_1_SQRT_2, 0.0))]
    // Rotation around z-axis by π/4 (will rotate by π/2 when applied)
    #[case([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_4, 
           Rotor::new(std::f64::consts::FRAC_1_SQRT_2, 0.0, 0.0, std::f64::consts::FRAC_1_SQRT_2))]
    // Zero vector should return identity
    #[case([0.0, 0.0, 0.0], std::f64::consts::PI, Rotor::identity())]
    fn test_from_axis_angle(#[case] axis: [f64; 3], #[case] angle: f64, #[case] expected: Rotor) {
        let result = Rotor::from_axis_angle(axis, angle);
        
        // Use the approx_eq method
        assert!(result.approx_eq(&expected, 1e-10), 
                "Expected {:?} but got {:?}", expected, result);
    }

    #[test]
    fn test_normalized_display() {
        // Test with normalized rotor
        let r4 = Rotor::new(1.0, 3.0, 2.0, 1.0).normalize();
        // We don't test the exact string since normalization will produce floating
        // point values. Just make sure it formats without errors
        let formatted = format!("{}", r4);
        assert!(formatted.contains("Ix"));
        assert!(formatted.contains("Iy"));
        assert!(formatted.contains("Iz"));
    }
    
    #[rstest]
    // Identity * Identity = Identity
    #[case(
        Rotor::identity(),
        Rotor::identity(),
        Rotor::identity()
    )]
    // Identity * R = R
    #[case(
        Rotor::identity(),
        Rotor::new(2.0, 3.0, 4.0, 5.0),
        Rotor::new(2.0, 3.0, 4.0, 5.0)
    )]
    // R * Identity = R
    #[case(
        Rotor::new(2.0, 3.0, 4.0, 5.0),
        Rotor::identity(),
        Rotor::new(2.0, 3.0, 4.0, 5.0)
    )]
    // Specific example: (1 + 1Ix) * (1 + 1Iy) = 1 + 1Ix + 1Iy + 1Iz
    #[case(
        Rotor::new(1.0, 1.0, 0.0, 0.0),
        Rotor::new(1.0, 0.0, 1.0, 0.0),
        Rotor::new(1.0, 1.0, 1.0, 1.0)
    )]
    fn test_rotor_multiplication(#[case] a: Rotor, #[case] b: Rotor, #[case] expected: Rotor) {
        let result = a * b;
        
        // Use the approx_eq method
        assert!(result.approx_eq(&expected, 1e-10), 
                "Expected {:?} but got {:?}", expected, result);
    }
    
    #[rstest]
    // Rotation around x-axis by π/4 should flip y to z when applied (which is a π/2 rotation)
    #[case(
        Rotor::from_axis_angle([1.0, 0.0, 0.0], std::f64::consts::FRAC_PI_4),
        Rotor::new(0.0, 0.0, 1.0, 0.0),  // Pure y-axis bivector (e31)
        Rotor::new(0.0, 0.0, 0.0, 1.0)   // Pure z-axis bivector (e12)
    )]
    // Rotation around y-axis by π/4 should flip z to x when applied (which is a π/2 rotation)
    #[case(
        Rotor::from_axis_angle([0.0, 1.0, 0.0], std::f64::consts::FRAC_PI_4),
        Rotor::new(0.0, 0.0, 0.0, 1.0),  // Pure z-axis bivector (e12)
        Rotor::new(0.0, 1.0, 0.0, 0.0)   // Pure x-axis bivector (e23)
    )]
    // Rotation around z-axis by π/4 should flip x to y when applied (which is a π/2 rotation)
    #[case(
        Rotor::from_axis_angle([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_4),
        Rotor::new(0.0, 1.0, 0.0, 0.0),  // Pure x-axis bivector (e23)
        Rotor::new(0.0, 0.0, 1.0, 0.0)   // Pure y-axis bivector (e31)
    )]
    fn test_rotor_rotation(#[case] r: Rotor, #[case] v: Rotor, #[case] expected: Rotor) {
        // Apply the rotation using the sandwich product: r * v * r.reverse()
        let result = r * v * r.reverse();
        
        // Use the approx_eq method
        assert!(result.approx_eq(&expected, 1e-10), 
                "Expected {:?} but got {:?}", expected, result);
    }
    
    #[rstest]
    // Basic addition
    #[case(
        Rotor::new(1.0, 2.0, 3.0, 4.0),
        Rotor::new(5.0, 6.0, 7.0, 8.0),
        Rotor::new(6.0, 8.0, 10.0, 12.0)
    )]
    // Addition with identity
    #[case(
        Rotor::identity(),
        Rotor::new(5.0, 6.0, 7.0, 8.0),
        Rotor::new(6.0, 6.0, 7.0, 8.0)
    )]
    // Addition with negative components
    #[case(
        Rotor::new(1.0, -2.0, 3.0, -4.0),
        Rotor::new(5.0, 6.0, -7.0, 8.0),
        Rotor::new(6.0, 4.0, -4.0, 4.0)
    )]
    fn test_rotor_addition(#[case] a: Rotor, #[case] b: Rotor, #[case] expected: Rotor) {
        let result = a + b;
        
        assert!(result.approx_eq(&expected, 1e-10), 
                "Expected {:?} but got {:?}", expected, result);
    }
    
    #[rstest]
    // Basic subtraction
    #[case(
        Rotor::new(5.0, 6.0, 7.0, 8.0),
        Rotor::new(1.0, 2.0, 3.0, 4.0),
        Rotor::new(4.0, 4.0, 4.0, 4.0)
    )]
    // Subtraction with identity
    #[case(
        Rotor::new(5.0, 6.0, 7.0, 8.0),
        Rotor::identity(),
        Rotor::new(4.0, 6.0, 7.0, 8.0)
    )]
    // Subtraction with negative components
    #[case(
        Rotor::new(5.0, 6.0, -7.0, 8.0),
        Rotor::new(1.0, -2.0, 3.0, -4.0),
        Rotor::new(4.0, 8.0, -10.0, 12.0)
    )]
    fn test_rotor_subtraction(#[case] a: Rotor, #[case] b: Rotor, #[case] expected: Rotor) {
        let result = a - b;
        
        assert!(result.approx_eq(&expected, 1e-10), 
                "Expected {:?} but got {:?}", expected, result);
    }
    
    #[rstest]
    // Division by identity
    #[case(
        Rotor::new(2.0, 3.0, 4.0, 5.0),
        Rotor::identity(),
        Rotor::new(2.0, 3.0, 4.0, 5.0)
    )]
    // Division of identity by a rotor
    #[case(
        Rotor::identity(),
        Rotor::new(2.0, 0.0, 0.0, 0.0),
        Rotor::new(0.5, 0.0, 0.0, 0.0)
    )]
    // Division by a normalized rotor should be multiplication by its reverse
    #[case(
        Rotor::new(2.0, 3.0, 4.0, 5.0),
        Rotor::from_axis_angle([1.0, 0.0, 0.0], std::f64::consts::FRAC_PI_4),
        Rotor::new(2.0, 3.0, 4.0, 5.0) * Rotor::from_axis_angle([1.0, 0.0, 0.0], std::f64::consts::FRAC_PI_4).reverse()
    )]
    // Self-division should yield identity
    #[case(
        Rotor::new(2.0, 3.0, 4.0, 5.0),
        Rotor::new(2.0, 3.0, 4.0, 5.0),
        Rotor::identity()
    )]
    fn test_rotor_division(#[case] a: Rotor, #[case] b: Rotor, #[case] expected: Rotor) {
        let result = a / b;
        
        assert!(result.approx_eq(&expected, 1e-10), 
                "Expected {:?} but got {:?}", expected, result);
    }
}
