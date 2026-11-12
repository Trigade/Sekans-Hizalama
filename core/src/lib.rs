pub mod nw;
pub mod scoring;
pub mod sw;

use scoring::{AlignError, Scoring};
use wasm_bindgen::prelude::*;

pub const MAX_SEQ_LEN: usize = 5_000;

pub const MAX_SEQ_LEN_SCORE_ONLY: usize = 100_000;

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

fn to_js_err(e: AlignError) -> JsError {
    JsError::new(&e.to_string())
}

#[wasm_bindgen(js_name = alignGlobal)]
pub fn align_global(
    a: &str,
    b: &str,
    match_score: i32,
    mismatch: i32,
    gap: i32,
) -> Result<JsValue, JsError> {
    let (a, b) = (normalize(a), normalize(b));
    check_len(&a, MAX_SEQ_LEN).map_err(to_js_err)?;
    check_len(&b, MAX_SEQ_LEN).map_err(to_js_err)?;

    let result = nw::needleman_wunsch(&a, &b, Scoring::new(match_score, mismatch, gap))
        .map_err(to_js_err)?;

    serde_wasm_bindgen::to_value(&result).map_err(|e| JsError::new(&e.to_string()))
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

    nw::needleman_wunsch_score(&a, &b, Scoring::new(match_score, mismatch, gap))
        .map_err(to_js_err)
}

#[wasm_bindgen(js_name = alignLocal)]
pub fn align_local(
    a: &str,
    b: &str,
    match_score: i32,
    mismatch: i32,
    gap: i32,
) -> Result<JsValue, JsError> {
    let (a, b) = (normalize(a), normalize(b));
    check_len(&a, MAX_SEQ_LEN).map_err(to_js_err)?;
    check_len(&b, MAX_SEQ_LEN).map_err(to_js_err)?;

    let result = sw::smith_waterman(&a, &b, Scoring::new(match_score, mismatch, gap))
        .map_err(to_js_err)?;

    serde_wasm_bindgen::to_value(&result).map_err(|e| JsError::new(&e.to_string()))
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

    sw::smith_waterman_score(&a, &b, Scoring::new(match_score, mismatch, gap))
        .map_err(to_js_err)
}

#[wasm_bindgen]
pub fn version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

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
}