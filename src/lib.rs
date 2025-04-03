pub mod rotor;
pub use rotor::Rotor;

/// The X Pauli operator bivector (e₂₃)
#[allow(non_upper_case_globals)]
pub const Ix: Rotor = Rotor::new(0.0, 1.0, 0.0, 0.0);

/// The Y Pauli operator bivector (e₃₁)
#[allow(non_upper_case_globals)]
pub const Iy: Rotor = Rotor::new(0.0, 0.0, 1.0, 0.0);

/// The Z Pauli operator bivector (e₁₂)
#[allow(non_upper_case_globals)]
pub const Iz: Rotor = Rotor::new(0.0, 0.0, 0.0, 1.0);
