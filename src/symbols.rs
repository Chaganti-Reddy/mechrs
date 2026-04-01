//! Symbol table with lifetime-bounded expression references.
//!
//! `Expr<'sym>` borrows its name from the owning `SymbolTable` — it cannot outlive it.
//! This is enforced at compile time with no runtime cost.

use std::collections::HashMap;

/// A table of named physical quantities.
pub struct SymbolTable {
    vars: HashMap<String, f64>,
}

/// An expression borrowed from a [`SymbolTable`].
///
/// The lifetime `'sym` ensures this cannot outlive the table it was created from.
pub struct Expr<'sym> {
    /// Name borrowed (zero-copy) from the symbol table.
    pub name: &'sym str,
    pub value: f64,
    pub unit: &'static str,
}

impl SymbolTable {
    /// Create an empty symbol table.
    pub fn new() -> Self {
        Self {
            vars: HashMap::new(),
        }
    }

    /// Insert or update a named value.
    pub fn insert(&mut self, name: impl Into<String>, value: f64) {
        self.vars.insert(name.into(), value);
    }

    /// Borrow a single expression by name.
    ///
    /// The returned `Expr<'_>` borrows from `self` — it cannot outlive this table.
    pub fn bind(&self, name: &str) -> Option<Expr<'_>> {
        self.vars.get_key_value(name).map(|(k, &v)| Expr {
            name: k.as_str(),
            value: v,
            unit: infer_unit(k),
        })
    }

    /// Iterate over all bound expressions.
    pub fn bind_all(&self) -> impl Iterator<Item = Expr<'_>> {
        self.vars.iter().map(|(k, &v)| Expr {
            name: k.as_str(),
            value: v,
            unit: infer_unit(k),
        })
    }
}

impl Default for SymbolTable {
    fn default() -> Self {
        Self::new()
    }
}

/// Infer a SI unit string from a common variable name.
fn infer_unit(name: &str) -> &'static str {
    match name {
        "velocity" | "speed" | "vx" | "vy" | "v0" => "m/s",
        "mass" | "m" => "kg",
        "length" | "l" | "radius" | "r" => "m",
        "time" | "t" => "s",
        "angle" | "theta" => "rad",
        "energy" | "e" => "J",
        "force" | "f" => "N",
        "acceleration" | "a" | "g" => "m/s²",
        _ => "",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn expr_borrows_name_not_cloned() {
        let mut table = SymbolTable::new();
        table.insert("velocity", 9.8);
        let expr = table.bind("velocity").unwrap();
        assert_eq!(expr.name, "velocity");
        assert_eq!(expr.value, 9.8);
        assert_eq!(expr.unit, "m/s");
    }

    #[test]
    fn bind_all_yields_all_entries() {
        let mut table = SymbolTable::new();
        table.insert("mass", 1.0);
        table.insert("velocity", 2.0);
        let count = table.bind_all().count();
        assert_eq!(count, 2);
    }

    #[test]
    fn bind_missing_returns_none() {
        let table = SymbolTable::new();
        assert!(table.bind("missing").is_none());
    }
}
