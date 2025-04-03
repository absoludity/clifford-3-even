# Rotors - the even sub-algebra of the Clifford 3 dimensional algebra

A Rust implementation of simple rotors in the even sub-algebra of Cl(3).

## What are Rotors?

Rotors are mathematical objects used to represent rotations in 3D space. Unlike
matrices or quaternions, rotors arise naturally from the geometric product in
Clifford Algebra.

A rotor is composed of a scalar part and three bivector components:
```
r = b_0 + b_x yz + b_y zx + b_z xy = b_0 + b_x Ix + b_y Iy + bz Iz
```

where:
- b_0, b_x, b_y and b_z are scalar numbers,
- x, y and z are orthogonal unit vectors (xx = yy = zz = 1), which anti-commute with each other (xy = -yx),
- I is the vector product xyz, a pseudo-scalar (it commutes like a scalar). So Ix = xyzx = yz

## Features

- Create rotors from axis-angle representation
- Perform rotor multiplication, addition, subtraction, and division
- Normalize rotors
- Calculate the reverse of a rotor
- Apply rotations using the sandwich product (r * v * r.reverse())
- The library exports 3 consts, Ix, Iy and Iz, for convenience.

## Usage

```rust
use clifford_3_even::Rotor;

// Create a rotor representing a rotation around the x-axis by π/4 radians
let r = Rotor::from_axis_angle([1.0, 0.0, 0.0], std::f64::consts::FRAC_PI_4);

// When applied using the sandwich product, this will rotate by π/2 radians
// For example, to rotate a bivector representing the y-axis:
let y_axis = Rotor::new(0.0, 0.0, 1.0, 0.0);  // Pure y-axis bivector (e31)
let rotated = r * y_axis * r.reverse();
// rotated will be approximately equal to Rotor::new(0.0, 0.0, 0.0, 1.0) (z-axis)
```

## License

This project is licensed under either of:

* MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)
* Apache License, Version 2.0, ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)

at your option.
