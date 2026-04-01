//! # mechrs
//!
//! Classical mechanics as lazy Rust iterators.
//! Every physical system yields an infinite iterator — pull values on demand, zero heap allocation.

pub mod integrators;
pub mod kinematics;
pub mod orbital;
pub mod oscillator;
pub mod pendulum;
pub mod symbols;
