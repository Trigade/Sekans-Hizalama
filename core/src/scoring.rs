use serde::{Serialize, Serializer};
use std::fmt;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Scoring {
    pub match_score: i32,
    pub mismatch: i32,
    pub gap: i32,
}

impl Default for Scoring {
    fn default() -> Self {
        Self { match_score: 1, mismatch: -1, gap: -1 }
    }
}

impl Scoring {
    pub const fn new(match_score: i32, mismatch: i32, gap: i32) -> Self {
        Self { match_score, mismatch, gap }
    }

    pub fn validate(&self) -> Result<(), AlignError> {
        if self.match_score <= 0 {
            return Err(AlignError::InvalidScoring("match_score pozitif olmalı"));
        }
        if self.mismatch >= self.match_score {
            return Err(AlignError::InvalidScoring("mismatch, match_score'dan küçük olmalı"));
        }
        if self.gap > 0 {
            return Err(AlignError::InvalidScoring("gap sıfır veya negatif olmalı"));
        }
        Ok(())
    }

    #[inline(always)]
    pub fn substitution(&self, x: u8, y: u8) -> i32 {
        if x == y { self.match_score } else { self.mismatch }
    }
}

#[repr(u8)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum Dir {
    #[default]
    None = 0,
    Diag = 1,
    Up = 2,
    Left = 3,
}

impl Serialize for Dir {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_u8(*self as u8)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AlignError {
    InvalidScoring(&'static str),
    InvalidChar { seq: char, pos: usize, byte: u8 },
    TooLong { len: usize, max: usize },
}

impl fmt::Display for AlignError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidScoring(msg) => write!(f, "Geçersiz skorlama: {msg}"),
            Self::InvalidChar { seq, pos, byte } => write!(
                f,
                "Dizi {seq}, konum {pos}: geçersiz karakter (0x{byte:02X}). \
                 Yalnızca A-Z harfleri kabul edilir."
            ),
            Self::TooLong { len, max } => {
                write!(f, "Dizi çok uzun: {len} karakter (üst sınır {max})")
            }
        }
    }
}

impl std::error::Error for AlignError {}

pub fn validate_sequence(seq: &[u8], label: char) -> Result<(), AlignError> {
    match seq.iter().position(|b| !b.is_ascii_uppercase()) {
        Some(pos) => Err(AlignError::InvalidChar { seq: label, pos, byte: seq[pos] }),
        None => Ok(()),
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct AlignmentResult {
    pub score: i32,
    pub aligned_a: String,
    pub aligned_b: String,
    pub rows: usize,
    pub cols: usize,
    pub matrix: Vec<i32>,
    pub dirs: Vec<Dir>,
    pub path: Vec<(usize, usize)>,
    pub start: (usize, usize),
    pub end: (usize, usize),
    pub stats: AlignmentStats,
}

impl AlignmentResult {
    #[inline]
    pub fn idx(&self, i: usize, j: usize) -> usize {
        i * self.cols + j
    }

    #[inline]
    pub fn cell(&self, i: usize, j: usize) -> i32 {
        self.matrix[self.idx(i, j)]
    }

    #[inline]
    pub fn dir(&self, i: usize, j: usize) -> Dir {
        self.dirs[self.idx(i, j)]
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct AlignmentStats {
    pub length: usize,
    pub matches: usize,
    pub mismatches: usize,
    pub gaps: usize,
    pub identity: f64,
}

impl AlignmentStats {
    pub fn from_aligned(a: &[u8], b: &[u8]) -> Self {
        debug_assert_eq!(a.len(), b.len());
        let (mut matches, mut mismatches, mut gaps) = (0, 0, 0);
        for (&x, &y) in a.iter().zip(b) {
            if x == b'-' || y == b'-' {
                gaps += 1;
            } else if x == y {
                matches += 1;
            } else {
                mismatches += 1;
            }
        }
        let length = a.len();
        let identity = if length == 0 { 0.0 } else { matches as f64 * 100.0 / length as f64 };
        Self { length, matches, mismatches, gaps, identity }
    }
}