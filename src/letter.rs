/// Letter type is just a character.
pub type Letter = char;

/// Constant definition of a gap "letter".
pub const GAP: Letter = '-';

/// Extension trait over primitive letter types.
pub trait NormalizeLetter {
    /// This method normalizes `Self` into a value of `Letter` type.
    ///
    /// E.g. `None` becomes `'-'` (gap).
    fn normalize_letter(self) -> Letter;
}

// reflexive implementation
impl NormalizeLetter for Letter {
    fn normalize_letter(self) -> Letter {
        self
    }
}

// generic reference auto-implementation
impl<'a, L> NormalizeLetter for &'a L
where
    L: NormalizeLetter + Copy,
{
    fn normalize_letter(self) -> Letter {
        (*self).normalize_letter()
    }
}

// normalizes optional letter into non-optional
impl<L> NormalizeLetter for Option<L>
where
    L: NormalizeLetter,
{
    fn normalize_letter(self) -> Letter {
        self.map_or(GAP, L::normalize_letter)
    }
}
