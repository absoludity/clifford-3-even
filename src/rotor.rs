/// Represents a rotor in the even sub-algebra of Cl(3).
///
/// A rotor is composed of a scalar part and three bivector components:
/// r = b₀ + b₁Ix + b₂Iy + b₃Iz
/// where b₀ is the scalar and b₁, b₂, b₃ are the bivector components.
use std::fmt;
use std::ops::{Add, Div, Mul, Sub};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rotor {
    /// Scalar component
    pub scalar: f64,
    /// Ix bivector component
    pub ix: f64,
    /// Iy bivector component
    pub iy: f64,
    /// Iz bivector component
    pub iz: f64,
}

impl fmt::Display for Rotor {
    /// Formats the rotor with configurable decimal precision.
    ///
    /// # Examples
    ///
    /// ```
    /// use clifford_3_even::Rotor;
    ///
    /// let r = Rotor::new(1.23456789, 2.34567891, 3.45678912, 4.56789123);
    ///
    /// // Default formatting (uses f64 default)
    /// assert_eq!(format!("{}", r), "1.23456789 + 2.34567891Ix + 3.45678912Iy + 4.56789123Iz");
    ///
    /// // With 2 decimal places
    /// assert_eq!(format!("{:.2}", r), "1.23 + 2.35Ix + 3.46Iy + 4.57Iz");
    ///
    /// // With 0 decimal places
    /// assert_eq!(format!("{:.0}", r), "1 + 2Ix + 3Iy + 5Iz");
    /// ```
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Check which components are non-zero
        let scalar_nonzero = self.scalar.abs() > 1e-10;
        let ix_nonzero = self.ix.abs() > 1e-10;
        let iy_nonzero = self.iy.abs() > 1e-10;
        let iz_nonzero = self.iz.abs() > 1e-10;

        // Check if we have any non-zero bivector components
        let has_bivector_terms = ix_nonzero || iy_nonzero || iz_nonzero;

        // Only show scalar if it's non-zero OR if there are no other terms
        let show_scalar = scalar_nonzero || !has_bivector_terms;

        let mut first_term = true;

        // Use precision from formatter if specified, otherwise use default f64 formatting
        if let Some(precision) = f.precision() {
            // Show scalar part if needed
            if show_scalar {
                write!(f, "{:.prec$}", self.scalar, prec = precision)?;
                first_term = false;
            }

            // Add each bivector component with its corresponding label
            if ix_nonzero {
                if first_term {
                    write!(f, "{:.prec$}Ix", self.ix, prec = precision)?;
                    first_term = false;
                } else {
                    write!(
                        f,
                        "{}{:.prec$}Ix",
                        if self.ix >= 0.0 { " + " } else { " - " },
                        self.ix.abs(),
                        prec = precision
                    )?;
                }
            }

            if iy_nonzero {
                if first_term {
                    write!(f, "{:.prec$}Iy", self.iy, prec = precision)?;
                    first_term = false;
                } else {
                    write!(
                        f,
                        "{}{:.prec$}Iy",
                        if self.iy >= 0.0 { " + " } else { " - " },
                        self.iy.abs(),
                        prec = precision
                    )?;
                }
            }

            if iz_nonzero {
                if first_term {
                    write!(f, "{:.prec$}Iz", self.iz, prec = precision)?;
                } else {
                    write!(
                        f,
                        "{}{:.prec$}Iz",
                        if self.iz >= 0.0 { " + " } else { " - " },
                        self.iz.abs(),
                        prec = precision
                    )?;
                }
            }
        } else {
            // Use default f64 formatting when no precision is specified
            // Show scalar part if needed
            if show_scalar {
                write!(f, "{}", self.scalar)?;
                first_term = false;
            }

            // Add each bivector component with its corresponding label
            if ix_nonzero {
                if first_term {
                    write!(f, "{}Ix", self.ix)?;
                    first_term = false;
                } else {
                    write!(
                        f,
                        "{}{}Ix",
                        if self.ix >= 0.0 { " + " } else { " - " },
                        self.ix.abs()
                    )?;
                }
            }

            if iy_nonzero {
                if first_term {
                    write!(f, "{}Iy", self.iy)?;
                    first_term = false;
                } else {
                    write!(
                        f,
                        "{}{}Iy",
                        if self.iy >= 0.0 { " + " } else { " - " },
                        self.iy.abs()
                    )?;
                }
            }

            if iz_nonzero {
                if first_term {
                    write!(f, "{}Iz", self.iz)?;
                } else {
                    write!(
                        f,
                        "{}{}Iz",
                        if self.iz >= 0.0 { " + " } else { " - " },
                        self.iz.abs()
                    )?;
                }
            }
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
            self.ix + rhs.ix,
            self.iy + rhs.iy,
            self.iz + rhs.iz,
        )
    }
}

impl Sub for Rotor {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        // Subtraction is component-wise
        Self::new(
            self.scalar - rhs.scalar,
            self.ix - rhs.ix,
            self.iy - rhs.iy,
            self.iz - rhs.iz,
        )
    }
}

impl Mul for Rotor {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        // Rotor multiplication formula for:
        // (a_0 + a_x*Ix + a_y*Iy + a_z*Iz) * (b_0 + b_x*Ix + b_y*Iy + b_z*Iz)

        // Scalar part: a_0*b_0 - a_x*b_x - a_y*b_y - a_z*b_z
        let scalar = self.iz.mul_add(
            -rhs.iz,
            self.iy.mul_add(
                -rhs.iy,
                self.scalar.mul_add(rhs.scalar, -(self.ix * rhs.ix)),
            ),
        );

        // Ix part: a_0*b_x + a_x*b_0 - a_y*b_z + a_z*b_y
        let ix = self.iz.mul_add(
            rhs.iy,
            self.iy
                .mul_add(-rhs.iz, self.scalar.mul_add(rhs.ix, self.ix * rhs.scalar)),
        );

        // Iy part: a_0*b_y + a_y*b_0 - a_z*b_x + a_x*b_z
        let iy = self.ix.mul_add(
            rhs.iz,
            self.iz
                .mul_add(-rhs.ix, self.scalar.mul_add(rhs.iy, self.iy * rhs.scalar)),
        );

        // Iz part: a_0*b_z + a_z*b_0 - a_x*b_y + a_y*b_x
        let iz = self.iy.mul_add(
            rhs.ix,
            self.ix
                .mul_add(-rhs.iy, self.scalar.mul_add(rhs.iz, self.iz * rhs.scalar)),
        );

        Self::new(scalar, ix, iy, iz)
    }
}

impl Div for Rotor {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        // Division by a rotor r is multiplication by its inverse r⁻¹
        // For a normalized rotor, the inverse is the reverse
        // For a general rotor, we need to divide by the norm squared

        // Calculate the norm squared
        let norm_squared = rhs.iz.mul_add(
            rhs.iz,
            rhs.iy
                .mul_add(rhs.iy, rhs.scalar.mul_add(rhs.scalar, rhs.ix * rhs.ix)),
        );

        // Calculate the inverse (reverse divided by norm squared)
        let inverse = Self::new(
            rhs.scalar / norm_squared,
            -rhs.ix / norm_squared,
            -rhs.iy / norm_squared,
            -rhs.iz / norm_squared,
        );

        // Multiply by the inverse
        self * inverse
    }
}

impl Rotor {
    /// Creates a new rotor with the given components.
    #[must_use]
    pub const fn new(scalar: f64, ix: f64, iy: f64, iz: f64) -> Self {
        Self { scalar, ix, iy, iz }
    }

    /// Creates an identity rotor (scalar = 1, bivectors = 0).
    #[must_use]
    pub const fn identity() -> Self {
        Self::new(1.0, 0.0, 0.0, 0.0)
    }

    /// Creates a rotor representing a rotation around the specified axis by the given angle.
    ///
    /// # Arguments
    ///
    /// * `axis` - A 3D vector [x, y, z] representing the axis of rotation
    /// * `angle` - The angle of rotation in radians
    ///   Note: When applying this rotor using the sandwich product (r * v * `r.reverse()`),
    ///   the resulting rotation will be by 2*angle.
    ///
    /// # Returns
    ///
    /// A rotor representing the specified rotation
    #[must_use]
    pub fn from_axis_angle(axis: [f64; 3], angle: f64) -> Self {
        // Normalize the axis vector
        let magnitude = axis[2]
            .mul_add(axis[2], axis[0].mul_add(axis[0], axis[1] * axis[1]))
            .sqrt();

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
        let ix = normalized_axis[0] * sin_angle;
        let iy = normalized_axis[1] * sin_angle;
        let iz = normalized_axis[2] * sin_angle;

        Self::new(cos_angle, ix, iy, iz)
    }

    /// Returns the squared magnitude of the rotor.
    #[must_use]
    pub fn magnitude_squared(&self) -> f64 {
        self.iz.mul_add(
            self.iz,
            self.iy
                .mul_add(self.iy, self.scalar.mul_add(self.scalar, self.ix * self.ix)),
        )
    }

    /// Normalizes the rotor to ensure it represents a proper rotation.
    /// Returns the norm (the factor by which the rotor was divided).
    pub fn normalize(&mut self) -> f64 {
        let norm_squared = self.magnitude_squared();
        let norm = norm_squared.sqrt();

        self.scalar /= norm;
        self.ix /= norm;
        self.iy /= norm;
        self.iz /= norm;

        norm
    }

    /// Returns the reverse of this rotor.
    ///
    /// The reverse of a rotor flips the sign of all bivector components
    /// while keeping the scalar component unchanged.
    #[must_use]
    pub fn reverse(&self) -> Self {
        Self {
            scalar: self.scalar,
            ix: -self.ix,
            iy: -self.iy,
            iz: -self.iz,
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
    #[must_use]
    pub fn approx_eq(&self, other: &Self, epsilon: f64) -> bool {
        (self.scalar - other.scalar).abs() < epsilon
            && (self.ix - other.ix).abs() < epsilon
            && (self.iy - other.iy).abs() < epsilon
            && (self.iz - other.iz).abs() < epsilon
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Ix, Iz};
    use rstest::rstest;

    #[rstest]
    #[case(Rotor::new(1.0, 3.0, 2.0, 1.0), "1 + 3Ix + 2Iy + 1Iz")]
    #[case(Rotor::new(1.0, -3.0, 2.0, -1.0), "1 - 3Ix + 2Iy - 1Iz")]
    #[case(Rotor::new(1.0, 0.0, 0.0, 0.0), "1")]
    #[case(Rotor::new(0.0, 2.0, 3.0, 1.0), "2Ix + 3Iy + 1Iz")]
    #[case(Rotor::new(0.0, 0.0, 0.0, 0.0), "0")]
    // Test single bivector component (no scalar)
    #[case(Rotor::new(0.0, 3.0, 0.0, 0.0), "3Ix")]
    #[case(Rotor::new(0.0, 0.0, -2.5, 0.0), "-2.5Iy")]
    #[case(Rotor::new(0.0, 0.0, 0.0, 1.0), "1Iz")]
    // Test negative first term
    #[case(Rotor::new(0.0, -1.0, 2.0, 0.0), "-1Ix + 2Iy")]
    fn test_display(#[case] rotor: Rotor, #[case] expected: &str) {
        assert_eq!(format!("{}", rotor), expected);
    }

    #[rstest]
    #[case(Rotor::new(1.0, 1.0, 2.0, 3.0), Rotor::new(1.0, -1.0, -2.0, -3.0))]
    #[case(Rotor::new(2.0, -1.0, 0.0, 3.0), Rotor::new(2.0, 1.0, 0.0, -3.0))]
    #[case(Rotor::new(1.0, 0.0, 0.0, 0.0), Rotor::new(1.0, 0.0, 0.0, 0.0))]
    // Test that Ix constant gets changed to -Ix under reversion
    #[case(Ix, Rotor::new(0.0, -1.0, 0.0, 0.0))]
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

        // Use the approx_eq method with Display format
        assert!(
            result.approx_eq(&expected, 1e-10),
            "Expected {} but got {}",
            expected,
            result
        );
    }

    #[test]
    fn test_normalized_display() {
        // Test with normalized rotor
        let mut r4 = Rotor::new(1.0, 3.0, 2.0, 1.0);
        let norm = r4.normalize();

        // Test the actual formatted output with 2 decimal places
        let formatted = format!("{:.2}", r4);
        // The original rotor (1, 3, 2, 1) has magnitude sqrt(15) ≈ 3.872983
        // After normalization: (1/3.872983, 3/3.872983, 2/3.872983, 1/3.872983)
        // Which gives approximately: (0.258199, 0.774597, 0.516398, 0.258199)
        assert_eq!(formatted, "0.26 + 0.77Ix + 0.52Iy + 0.26Iz");

        // Also verify that normalize returned the expected value.
        assert!(
            norm.eq(&15.0_f64.sqrt()),
            "Expected {} but got {}",
            15.0_f64.sqrt(),
            norm
        );
    }

    #[rstest]
    // Identity rotor should have magnitude_squared = 1
    #[case(Rotor::identity(), 1.0)]
    // Unit rotor with all components = 0.5 should have magnitude_squared = 1.0
    #[case(Rotor::new(0.5, 0.5, 0.5, 0.5), 1.0)]
    // Rotor with components (1, 2, 3, 4) should have magnitude_squared = 30
    #[case(Rotor::new(1.0, 2.0, 3.0, 4.0), 30.0)]
    // Rotor with components (3, 4, 0, 0) should have magnitude_squared = 25
    #[case(Rotor::new(3.0, 4.0, 0.0, 0.0), 25.0)]
    fn test_magnitude_squared(#[case] rotor: Rotor, #[case] expected: f64) {
        let result = rotor.magnitude_squared();
        assert!(
            (result - expected).abs() < 1e-10,
            "Expected {} but got {}",
            expected,
            result
        );

        // Also check that rotor * rotor.reverse() gives the same scalar value
        let product = rotor * rotor.reverse();
        assert!(
            (product.scalar - expected).abs() < 1e-10,
            "Expected rotor * rotor.reverse() scalar to be {} but got {}",
            expected,
            product.scalar
        );
    }

    #[test]
    fn test_precision_formatting() {
        let rotor = Rotor::new(1.23456789, 2.34567891, 3.45678912, 4.56789123);

        // Test default formatting (should use f64 default)
        let default_format = format!("{}", rotor);
        assert_eq!(
            default_format,
            "1.23456789 + 2.34567891Ix + 3.45678912Iy + 4.56789123Iz"
        );

        // Test with 2 decimal places
        let two_decimal = format!("{:.2}", rotor);
        assert_eq!(two_decimal, "1.23 + 2.35Ix + 3.46Iy + 4.57Iz");

        // Test with 0 decimal places (explicit)
        let no_decimal = format!("{:.0}", rotor);
        assert_eq!(no_decimal, "1 + 2Ix + 3Iy + 5Iz");

        // Test with high precision
        let high_precision = format!("{:.8}", rotor);
        assert_eq!(
            high_precision,
            "1.23456789 + 2.34567891Ix + 3.45678912Iy + 4.56789123Iz"
        );

        // Test precision formatting with zero scalar
        let zero_scalar = Rotor::new(0.0, 1.23456789, 2.34567891, 0.0);
        let zero_scalar_2dp = format!("{:.2}", zero_scalar);
        assert_eq!(zero_scalar_2dp, "1.23Ix + 2.35Iy");
    }

    #[test]
    fn test_original_formatting_preserved() {
        // Verify that original test cases still work with default f64 formatting
        let integer_values = Rotor::new(1.0, 3.0, 2.0, 1.0);
        assert_eq!(format!("{}", integer_values), "1 + 3Ix + 2Iy + 1Iz");

        let mixed_signs = Rotor::new(1.0, -3.0, 2.0, -1.0);
        assert_eq!(format!("{}", mixed_signs), "1 - 3Ix + 2Iy - 1Iz");

        let identity = Rotor::new(1.0, 0.0, 0.0, 0.0);
        assert_eq!(format!("{}", identity), "1");
    }

    #[rstest]
    // Identity * Identity = Identity
    #[case(Rotor::identity(), Rotor::identity(), Rotor::identity())]
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
    // Specific example: (1 + 1Ix) * (1 + 1Iy) = 1 + 1Ix + 1Iy - 1Iz
    #[case(
        Rotor::new(1.0, 1.0, 0.0, 0.0),
        Rotor::new(1.0, 0.0, 1.0, 0.0),
        Rotor::new(1.0, 1.0, 1.0, -1.0)
    )]
    fn test_rotor_multiplication(#[case] a: Rotor, #[case] b: Rotor, #[case] expected: Rotor) {
        let result = a * b;

        // Use the approx_eq method
        assert!(
            result.approx_eq(&expected, 1e-10),
            "Expected {:?} but got {:?}",
            expected,
            result
        );
    }

    #[rstest]
    // Rotation around x-axis by π/4 should flip y to z when applied (which is a π/2 rotation)
    #[case(
        Rotor::from_axis_angle([1.0, 0.0, 0.0], std::f64::consts::FRAC_PI_4),
        Rotor::new(0.0, 0.0, 1.0, 0.0),  // Pure y-axis bivector (e31)
        Rotor::new(0.0, 0.0, 0.0, -1.0)  // Pure z-axis bivector (e12) with negative sign
    )]
    // Rotation around y-axis by π/4 should flip z to x when applied (which is a π/2 rotation)
    #[case(
        Rotor::from_axis_angle([0.0, 1.0, 0.0], std::f64::consts::FRAC_PI_4),
        Rotor::new(0.0, 0.0, 0.0, 1.0),  // Pure z-axis bivector (e12)
        Rotor::new(0.0, -1.0, 0.0, 0.0)  // Pure x-axis bivector (e23) with negative sign
    )]
    // Rotation around z-axis by π/4 should flip x to y when applied (which is a π/2 rotation)
    #[case(
        Rotor::from_axis_angle([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_4),
        Rotor::new(0.0, 1.0, 0.0, 0.0),  // Pure x-axis bivector (e23)
        Rotor::new(0.0, 0.0, -1.0, 0.0)  // Pure y-axis bivector (e31) with negative sign
    )]
    fn test_rotor_rotation(#[case] r: Rotor, #[case] v: Rotor, #[case] expected: Rotor) {
        // Apply the rotation using the sandwich product: r * v * r.reverse()
        let result = r * v * r.reverse();

        // Use the approx_eq method
        assert!(
            result.approx_eq(&expected, 1e-10),
            "Expected {:?} but got {:?}",
            expected,
            result
        );
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

        assert!(
            result.approx_eq(&expected, 1e-10),
            "Expected {:?} but got {:?}",
            expected,
            result
        );
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

        assert!(
            result.approx_eq(&expected, 1e-10),
            "Expected {:?} but got {:?}",
            expected,
            result
        );
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

        assert!(
            result.approx_eq(&expected, 1e-10),
            "Expected {:?} but got {:?}",
            expected,
            result
        );
    }

    #[rstest]
    // Test Pauli X operation as a quantum gate using Ix.reverse() * ψ * Iz
    // In the computational basis:
    // |0⟩ → |1⟩ (represented as 1 → -Iy)
    #[case::zero_to_one(
        Rotor::new(1.0, 0.0, 0.0, 0.0),  // Input: |0⟩ state (scalar 1)
        Rotor::new(0.0, 0.0, -1.0, 0.0)  // Expected: |1⟩ state (-Iy)
    )]
    // |1⟩ → |0⟩ (represented as -Iy → 1)
    #[case::one_to_zero(
        Rotor::new(0.0, 0.0, -1.0, 0.0),  // Input: |1⟩ state (-Iy)
        Rotor::new(1.0, 0.0, 0.0, 0.0)    // Expected: |0⟩ state (scalar 1)
    )]
    // In the Pauli X eigenbasis:
    // |+⟩ → |+⟩ (represented as normalized (1-Iy))
    #[case::plus_eigenstate(
        Rotor::new(std::f64::consts::FRAC_1_SQRT_2, 0.0, -std::f64::consts::FRAC_1_SQRT_2, 0.0),  // Input: |+⟩ state
        Rotor::new(std::f64::consts::FRAC_1_SQRT_2, 0.0, -std::f64::consts::FRAC_1_SQRT_2, 0.0)   // Expected: |+⟩ state unchanged
    )]
    // |-⟩ → -|-⟩ (represented as normalized (1+Iy))
    #[case::minus_eigenstate(
        Rotor::new(std::f64::consts::FRAC_1_SQRT_2, 0.0, std::f64::consts::FRAC_1_SQRT_2, 0.0),   // Input: |-⟩ state
        Rotor::new(-std::f64::consts::FRAC_1_SQRT_2, 0.0, -std::f64::consts::FRAC_1_SQRT_2, 0.0)  // Expected: |-⟩ state with sign flipped
    )]
    fn test_pauli_x_quantum_gate(#[case] state: Rotor, #[case] expected: Rotor) {
        // Apply the Pauli X operation using: Ix.reverse() * state * Iz
        let result = Ix.reverse() * state * Iz;

        // Use the approx_eq method with Display format
        assert!(
            result.approx_eq(&expected, 1e-10),
            "Expected {} but got {}",
            expected,
            result
        );
    }
}
