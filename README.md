# Rotors - the even sub-algebra of the Clifford 3 dimensional algebra

A Rust implementation of simple rotors in the even sub-algebra of Cl(3).

## What are Rotors?

Rotors are mathematical objects used to represent rotations in 3D space. Unlike
matrices or quaternions, rotors arise naturally from the geometric product in
Clifford Algebra.

A rotor is composed of a scalar part and three bivector components:
```
r = b₀ + b₁e₂₃ + b₂e₃₁ + b₃e₁₂
```

where:
- b₀ is the scalar component
- b₁, b₂, b₃ are the bivector components corresponding to the three planes of rotation
- e₂₃, e₃₁, e₁₂ are the basis bivectors (often denoted as Ix, Iy, Iz)

## Features

- Create rotors from axis-angle representation
- Perform rotor multiplication, addition, subtraction, and division
- Normalize rotors
- Calculate the reverse of a rotor
- Apply rotations using the sandwich product (r * v * r.reverse())

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
