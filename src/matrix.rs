use std::ops::{Index, IndexMut};

use crate::values::Score;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct AlignmentMatrix {
    buf: Vec<Score>,
    width: usize,
}

impl AlignmentMatrix {
    pub fn new(height: usize, width: usize) -> Self {
        Self { buf: vec![0; height * width], width }
    }

    pub fn height(&self) -> usize {
        self.buf.len() / self.width()
    }

    pub fn width(&self) -> usize {
        self.width
    }

    pub fn get_ref(&self, i: usize, j: usize) -> Option<&Score> {
        if j >= self.width {
            None
        } else {
            self.buf.get(i * self.width + j)
        }
    }

    pub fn get_mut(&mut self, i: usize, j: usize) -> Option<&mut Score> {
        if j >= self.width {
            None
        } else {
            self.buf.get_mut(i * self.width + j)
        }
    }

    pub fn get(&self, i: usize, j: usize) -> Option<Score> {
        self.get_ref(i, j).copied()
    }

    #[must_use]
    pub fn set(&mut self, i: usize, j: usize, score: Score) -> bool {
        if let Some(ref_mut) = self.get_mut(i, j) {
            *ref_mut = score;
            true
        } else {
            false
        }
    }
}

impl Index<(usize, usize)> for AlignmentMatrix {
    type Output = Score;

    fn index(&self, (i, j): (usize, usize)) -> &Self::Output {
        self.get_ref(i, j).expect("invalid index")
    }
}

impl IndexMut<(usize, usize)> for AlignmentMatrix {
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut Self::Output {
        self.get_mut(i, j).expect("invalid index")
    }
}

impl Index<[usize; 2]> for AlignmentMatrix {
    type Output = Score;

    fn index(&self, index: [usize; 2]) -> &Self::Output {
        &self[(index[0], index[1])]
    }
}

impl IndexMut<[usize; 2]> for AlignmentMatrix {
    fn index_mut(&mut self, index: [usize; 2]) -> &mut Self::Output {
        &mut self[(index[0], index[1])]
    }
}
