# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Nothing yet

### Changed
- Nothing yet

### Fixed
- Nothing yet

## [0.1.1] - 2024-01-XX

### Added
- Configurable precision formatting for `Display` trait implementation
- `magnitude_squared()` method to calculate the squared magnitude of a rotor
- Comprehensive test suite with 45+ test cases covering:
  - Display formatting with various precision levels
  - Magnitude calculations with verification via `rotor * rotor.reverse()`
  - Edge cases for zero components and normalization
- GitHub Actions CI/CD pipeline with:
  - Multi-platform testing (Ubuntu, Windows, macOS)
  - Multiple Rust version support (MSRV 1.85.0, stable, beta, nightly)
  - Code formatting, linting, and documentation checks
  - Automated publishing workflow

### Changed
- `normalize()` method now modifies the rotor in-place and returns the norm factor
- Display formatting now omits zero scalar components when bivector terms are present
- Improved documentation with examples for precision formatting

### Fixed
- Documentation indentation issues that caused clippy warnings
- Recursive alias definition in cargo config that prevented formatting checks

## [0.1.0] - 2024-01-XX

### Added
- Initial implementation of `Rotor` struct for 3D Clifford algebra even sub-algebra
- Basic arithmetic operations: addition, subtraction, multiplication, division
- Rotor creation from axis-angle representation
- Normalization and reverse operations
- Display formatting with mathematical notation
- Comprehensive test suite for quantum gate operations (Pauli X gate)
- Support for rotation operations via sandwich product
- MIT/Apache-2.0 dual licensing

[Unreleased]: https://github.com/absoludity/clifford-3-even/compare/v0.1.1...HEAD
[0.1.1]: https://github.com/absoludity/clifford-3-even/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/absoludity/clifford-3-even/releases/tag/v0.1.0