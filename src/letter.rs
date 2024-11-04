pub type Letter = char;

pub const GAP: Letter = '-';

pub trait NormalizeLetter {
    fn normalize_letter(self) -> Letter;
}

impl NormalizeLetter for Letter {
    fn normalize_letter(self) -> Letter {
        self
    }
}

impl<'a, L> NormalizeLetter for &'a L
where
    L: NormalizeLetter + Copy,
{
    fn normalize_letter(self) -> Letter {
        (*self).normalize_letter()
    }
}

impl<L> NormalizeLetter for Option<L>
where
    L: NormalizeLetter,
{
    fn normalize_letter(self) -> Letter {
        self.map_or(GAP, L::normalize_letter)
    }
}
