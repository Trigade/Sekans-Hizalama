pub mod blosum;
pub mod gotoh;
pub mod hirschberg;
pub mod nw;
pub mod scoring;
pub mod sw;

use scoring::{AlignError, AlignmentResult, Scoring, Substitution};
use serde::Deserialize;
use wasm_bindgen::prelude::*;

pub const MAX_SEQ_LEN: usize = 50_000;
pub const MAX_SEQ_LEN_SCORE_ONLY: usize = 100_000;
pub const VIS_CELLS: usize = 120 * 120;
pub const MAX_FULL_CELLS: usize = 25_000_000;
pub const MAX_HIRSCHBERG_CELLS: usize = 100_000_000;

#[derive(Deserialize, Clone, Copy, Debug, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
pub enum Mode {
    Global,
    Local,
}

#[derive(Deserialize, Clone, Copy, Debug, PartialEq, Eq, Default)]
#[serde(rename_all = "lowercase")]
pub enum SubstKind {
    #[default]
    Simple,
    Blosum62,
}

#[derive(Deserialize, Clone, Copy, Debug)]
#[serde(rename_all = "camelCase")]
pub struct AlignOptions {
    pub mode: Mode,
    #[serde(default)]
    pub substitution: SubstKind,
    #[serde(rename = "match", default = "default_match")]
    pub match_score: i32,
    #[serde(default = "default_mismatch")]
    pub mismatch: i32,
    #[serde(default)]
    pub gap_open: i32,
    pub gap_extend: i32,
}

fn default_match() -> i32 {
    1
}

fn default_mismatch() -> i32 {
    -1
}

impl AlignOptions {
    pub fn scoring(&self) -> Scoring {
        let subst = match self.substitution {
            SubstKind::Simple => Substitution::Simple {
                match_score: self.match_score,
                mismatch: self.mismatch,
            },
            SubstKind::Blosum62 => Substitution::Blosum62,
        };
        Scoring::with_gaps(subst, self.gap_open, self.gap_extend)
    }
}

fn normalize(input: &str) -> Vec<u8> {
    input
        .bytes()
        .filter(|b| !b.is_ascii_whitespace())
        .map(|b| b.to_ascii_uppercase())
        .collect()
}

fn check_len(seq: &[u8], max: usize) -> Result<(), AlignError> {
    if seq.len() > max {
        Err(AlignError::TooLong { len: seq.len(), max })
    } else {
        Ok(())
    }
}

fn check_cells(cells: usize, max: usize, what: &'static str) -> Result<(), AlignError> {
    if cells > max {
        Err(AlignError::TooLarge { cells, max, what })
    } else {
        Ok(())
    }
}

pub fn run_alignment(a: &[u8], b: &[u8], opts: &AlignOptions) -> Result<AlignmentResult, AlignError> {
    check_len(a, MAX_SEQ_LEN)?;
    check_len(b, MAX_SEQ_LEN)?;

    let s = opts.scoring();
    s.validate()?;

    let cells = (a.len() + 1) * (b.len() + 1);
    let visual = cells <= VIS_CELLS;

    let mut result = match (opts.mode, s.is_affine()) {
        (Mode::Global, false) if visual => nw::needleman_wunsch(a, b, s)?,
        (Mode::Global, false) => {
            check_cells(cells, MAX_HIRSCHBERG_CELLS, "Hirschberg")?;
            hirschberg::hirschberg(a, b, s)?
        }
        (Mode::Local, false) => {
            check_cells(cells, MAX_FULL_CELLS, "Smith–Waterman")?;
            sw::smith_waterman(a, b, s)?
        }
        (Mode::Global, true) => {
            check_cells(cells, MAX_FULL_CELLS, "Gotoh")?;
            gotoh::gotoh_global(a, b, s, visual)?
        }
        (Mode::Local, true) => {
            check_cells(cells, MAX_FULL_CELLS, "Gotoh")?;
            gotoh::gotoh_local(a, b, s, visual)?
        }
    };

    if !visual {
        result.strip_matrices();
    }
    Ok(result)
}

fn to_js_err(e: AlignError) -> JsError {
    JsError::new(&e.to_string())
}

fn parse_options(opts: JsValue) -> Result<AlignOptions, JsError> {
    serde_wasm_bindgen::from_value(opts).map_err(|e| JsError::new(&format!("Geçersiz seçenekler: {e}")))
}

#[wasm_bindgen]
pub fn align(a: &str, b: &str, opts: JsValue) -> Result<JsValue, JsError> {
    let opts = parse_options(opts)?;
    let (a, b) = (normalize(a), normalize(b));
    let result = run_alignment(&a, &b, &opts).map_err(to_js_err)?;
    serde_wasm_bindgen::to_value(&result).map_err(|e| JsError::new(&e.to_string()))
}

#[wasm_bindgen(js_name = substitutionScore)]
pub fn substitution_score(x: &str, y: &str, opts: JsValue) -> Result<i32, JsError> {
    let opts = parse_options(opts)?;
    let first = |s: &str| s.bytes().next().map(|c| c.to_ascii_uppercase()).unwrap_or(b'X');
    Ok(opts.scoring().substitution(first(x), first(y)))
}

#[wasm_bindgen(js_name = alignGlobalScore)]
pub fn align_global_score(
    a: &str,
    b: &str,
    match_score: i32,
    mismatch: i32,
    gap: i32,
) -> Result<i32, JsError> {
    let (a, b) = (normalize(a), normalize(b));
    check_len(&a, MAX_SEQ_LEN_SCORE_ONLY).map_err(to_js_err)?;
    check_len(&b, MAX_SEQ_LEN_SCORE_ONLY).map_err(to_js_err)?;

    nw::needleman_wunsch_score(&a, &b, Scoring::new(match_score, mismatch, gap)).map_err(to_js_err)
}

#[wasm_bindgen(js_name = alignLocalScore)]
pub fn align_local_score(
    a: &str,
    b: &str,
    match_score: i32,
    mismatch: i32,
    gap: i32,
) -> Result<i32, JsError> {
    let (a, b) = (normalize(a), normalize(b));
    check_len(&a, MAX_SEQ_LEN_SCORE_ONLY).map_err(to_js_err)?;
    check_len(&b, MAX_SEQ_LEN_SCORE_ONLY).map_err(to_js_err)?;

    sw::smith_waterman_score(&a, &b, Scoring::new(match_score, mismatch, gap)).map_err(to_js_err)
}

#[wasm_bindgen]
pub fn version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn opts(mode: Mode, substitution: SubstKind, gap_open: i32, gap_extend: i32) -> AlignOptions {
        AlignOptions { mode, substitution, match_score: 1, mismatch: -1, gap_open, gap_extend }
    }

    #[test]
    fn normalize_bosluk_ve_buyuk_harf() {
        assert_eq!(normalize("ac gt\r\nAc\tgT"), b"ACGTACGT".to_vec());
    }

    #[test]
    fn uzunluk_limiti() {
        let long = vec![b'A'; MAX_SEQ_LEN + 1];
        assert!(check_len(&long, MAX_SEQ_LEN).is_err());
        assert!(check_len(&long[..MAX_SEQ_LEN], MAX_SEQ_LEN).is_ok());
    }

    #[test]
    fn algoritma_secimi() {
        let small = b"GATTACA";
        let big = vec![b'A'; 500];
        let cases = [
            (opts(Mode::Global, SubstKind::Simple, 0, -1), &small[..], "needleman-wunsch"),
            (opts(Mode::Global, SubstKind::Simple, 0, -1), &big[..], "hirschberg"),
            (opts(Mode::Local, SubstKind::Simple, 0, -1), &small[..], "smith-waterman"),
            (opts(Mode::Global, SubstKind::Blosum62, -11, -1), &small[..], "gotoh-global"),
            (opts(Mode::Local, SubstKind::Blosum62, -11, -1), &small[..], "gotoh-local"),
        ];
        for (o, seq, expected) in cases {
            let r = run_alignment(seq, seq, &o).unwrap();
            assert_eq!(r.algorithm, expected);
        }
    }

    #[test]
    fn buyuk_matriste_matris_gonderilmez() {
        let big = vec![b'A'; 500];
        let r = run_alignment(&big, &big, &opts(Mode::Local, SubstKind::Simple, 0, -1)).unwrap();
        assert!(r.matrix.is_empty() && r.dirs.is_empty() && r.layers.is_none());
        assert_eq!(r.score, 500);
    }

    #[test]
    fn kucuk_afin_katmanlari_icerir() {
        let r = run_alignment(b"ACGT", b"AGT", &opts(Mode::Global, SubstKind::Simple, -3, -1)).unwrap();
        assert!(r.layers.is_some());
        assert_eq!(r.matrix.len(), 5 * 4);
    }
}