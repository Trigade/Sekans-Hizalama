use crate::scoring::{
    validate_sequence, AlignError, AlignmentResult, AlignmentStats, Dir, Layers, Scoring,
};

pub const NEG: i32 = i32::MIN / 4;

const ST_M: u8 = 0;
const ST_X: u8 = 1;
const ST_Y: u8 = 2;
const ST_START: u8 = 3;

#[inline(always)]
fn best3(m: i32, x: i32, y: i32) -> (i32, u8) {
    if m >= x && m >= y {
        (m, ST_M)
    } else if x >= y {
        (x, ST_X)
    } else {
        (y, ST_Y)
    }
}

#[inline(always)]
fn dir_of(state: u8) -> Dir {
    match state {
        ST_M => Dir::Diag,
        ST_X => Dir::Up,
        ST_Y => Dir::Left,
        _ => Dir::None,
    }
}

pub fn gotoh_global(
    a: &[u8],
    b: &[u8],
    s: Scoring,
    keep_layers: bool,
) -> Result<AlignmentResult, AlignError> {
    gotoh(a, b, s, false, keep_layers)
}

pub fn gotoh_local(
    a: &[u8],
    b: &[u8],
    s: Scoring,
    keep_layers: bool,
) -> Result<AlignmentResult, AlignError> {
    gotoh(a, b, s, true, keep_layers)
}

struct Grid {
    rows: usize,
    cols: usize,
    tb: Vec<u8>,
    matrix: Vec<i32>,
    dirs: Vec<Dir>,
    layers: Option<Layers>,
}

impl Grid {
    #[inline(always)]
    fn store(&mut self, k: usize, local: bool, m: i32, x: i32, y: i32) {
        let (v, st) = best3(m, x, y);
        if local && v <= 0 {
            self.matrix[k] = 0;
            self.dirs[k] = Dir::None;
        } else {
            self.matrix[k] = v;
            self.dirs[k] = dir_of(st);
        }
        if let Some(l) = self.layers.as_mut() {
            l.m[k] = m;
            l.x[k] = x;
            l.y[k] = y;
        }
    }
}

#[allow(clippy::needless_range_loop)]
fn gotoh(
    a: &[u8],
    b: &[u8],
    s: Scoring,
    local: bool,
    keep_layers: bool,
) -> Result<AlignmentResult, AlignError> {
    s.validate()?;
    validate_sequence(a, 'A')?;
    validate_sequence(b, 'B')?;

    let rows = a.len() + 1;
    let cols = b.len() + 1;
    let n = rows * cols;
    let e = s.gap_extend;
    let oe = s.gap_open + s.gap_extend;

    let mut g = Grid {
        rows,
        cols,
        tb: vec![0u8; n],
        matrix: vec![0i32; n],
        dirs: vec![Dir::None; n],
        layers: keep_layers.then(|| Layers { m: vec![NEG; n], x: vec![NEG; n], y: vec![NEG; n] }),
    };

    let mut pm = vec![NEG; cols];
    let mut px = vec![NEG; cols];
    let mut py = vec![NEG; cols];
    let mut cm = vec![NEG; cols];
    let mut cx = vec![NEG; cols];
    let mut cy = vec![NEG; cols];

    if !local {
        pm[0] = 0;
        for j in 1..cols {
            py[j] = s.gap_open + e * j as i32;
            g.tb[j] = (if j == 1 { ST_M } else { ST_Y }) << 4;
        }
    }
    for j in 0..cols {
        g.store(j, local, pm[j], px[j], py[j]);
    }

    let mut best = 0i32;
    let mut best_cell = (0usize, 0usize);

    for i in 1..rows {
        let ca = a[i - 1];
        let row = i * cols;

        cm[0] = NEG;
        cy[0] = NEG;
        cx[0] = if local { NEG } else { s.gap_open + e * i as i32 };
        if !local {
            g.tb[row] = (if i == 1 { ST_M } else { ST_X }) << 2;
        }
        g.store(row, local, cm[0], cx[0], cy[0]);

        for j in 1..cols {
            let (diag, dst) = best3(pm[j - 1], px[j - 1], py[j - 1]);
            let (mprev, mst) = if local && diag <= 0 { (0, ST_START) } else { (diag, dst) };
            let mv = mprev + s.substitution(ca, b[j - 1]);

            let (xv, xst) = best3(pm[j] + oe, px[j] + e, py[j] + oe);
            let (yv, yst) = best3(cm[j - 1] + oe, cx[j - 1] + oe, cy[j - 1] + e);

            cm[j] = mv;
            cx[j] = xv;
            cy[j] = yv;

            let k = row + j;
            g.tb[k] = mst | (xst << 2) | (yst << 4);
            g.store(k, local, mv, xv, yv);

            if local && mv > best {
                best = mv;
                best_cell = (i, j);
            }
        }

        std::mem::swap(&mut pm, &mut cm);
        std::mem::swap(&mut px, &mut cx);
        std::mem::swap(&mut py, &mut cy);
    }

    let (score, start_state, end_cell) = if local {
        (best, ST_M, best_cell)
    } else {
        let (v, st) = best3(pm[cols - 1], px[cols - 1], py[cols - 1]);
        (v, st, (rows - 1, cols - 1))
    };

    let cap = rows + cols;
    let mut out_a: Vec<u8> = Vec::with_capacity(cap);
    let mut out_b: Vec<u8> = Vec::with_capacity(cap);
    let mut path: Vec<(usize, usize)> = Vec::with_capacity(cap);

    if !local || score > 0 {
        let (mut i, mut j) = end_cell;
        let mut state = start_state;
        path.push((i, j));

        while i > 0 || j > 0 {
            let k = i * g.cols + j;
            let t = g.tb[k];
            g.dirs[k] = dir_of(state);
            let prev = match state {
                ST_M => {
                    out_a.push(a[i - 1]);
                    out_b.push(b[j - 1]);
                    i -= 1;
                    j -= 1;
                    t & 3
                }
                ST_X => {
                    out_a.push(a[i - 1]);
                    out_b.push(b'-');
                    i -= 1;
                    (t >> 2) & 3
                }
                _ => {
                    out_a.push(b'-');
                    out_b.push(b[j - 1]);
                    j -= 1;
                    (t >> 4) & 3
                }
            };
            path.push((i, j));
            if prev == ST_START {
                break;
            }
            state = prev;
        }
    }

    out_a.reverse();
    out_b.reverse();
    path.reverse();

    let stats = AlignmentStats::from_aligned(&out_a, &out_b);
    let start = path.first().copied().unwrap_or(end_cell);

    Ok(AlignmentResult {
        score,
        aligned_a: out_a.into_iter().map(char::from).collect(),
        aligned_b: out_b.into_iter().map(char::from).collect(),
        rows: g.rows,
        cols: g.cols,
        matrix: g.matrix,
        dirs: g.dirs,
        path,
        start,
        end: end_cell,
        stats,
        algorithm: if local { "gotoh-local" } else { "gotoh-global" },
        layers: g.layers,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scoring::Substitution;

    type Case = (&'static [u8], &'static [u8], i32, i32, bool, i32);

    fn rescore(r: &AlignmentResult, s: Scoring) -> i32 {
        let a = r.aligned_a.as_bytes();
        let b = r.aligned_b.as_bytes();
        let mut total = 0;
        let mut k = 0;
        while k < a.len() {
            if a[k] == b'-' || b[k] == b'-' {
                let in_a = a[k] == b'-';
                let mut len = 0;
                while k < a.len() && (if in_a { a[k] == b'-' } else { b[k] == b'-' }) {
                    len += 1;
                    k += 1;
                }
                total += s.gap_cost(len);
            } else {
                total += s.substitution(a[k], b[k]);
                k += 1;
            }
        }
        total
    }

    #[test]
    fn afin_global_tek_uzun_bosluk_tercih_eder() {
        let s = Scoring::with_gaps(Substitution::Simple { match_score: 2, mismatch: -1 }, -4, -1);
        let r = gotoh_global(b"ACGTTTTACG", b"ACGACG", s, false).unwrap();
        assert_eq!(r.aligned_b.matches('-').count(), 4);
        assert!(r.aligned_b.contains("----"));
        assert_eq!(rescore(&r, s), r.score);
    }

    #[test]
    fn acma_sifirsa_lineer_ile_ayni_skor() {
        let lin = Scoring::new(1, -1, -1);
        let aff = Scoring::with_gaps(Substitution::Simple { match_score: 1, mismatch: -1 }, 0, -1);
        for (a, b) in [(&b"GATTACA"[..], &b"GCATGCU"[..]), (b"HEAGAWGHEE", b"PAWHEAE"), (b"", b"ACG")] {
            let nw = crate::nw::needleman_wunsch(a, b, lin).unwrap().score;
            let g = gotoh_global(a, b, aff, false).unwrap().score;
            assert_eq!(nw, g);
            let sw = crate::sw::smith_waterman(a, b, lin).unwrap().score;
            let gl = gotoh_local(a, b, aff, false).unwrap().score;
            assert_eq!(sw, gl);
        }
    }

    #[test]
    fn global_yol_koseden_koseye_ve_skor_tutarli() {
        let s = Scoring::blosum62(-11, -1);
        let r = gotoh_global(b"HEAGAWGHEE", b"PAWHEAE", s, true).unwrap();
        assert_eq!(r.path.first(), Some(&(0, 0)));
        assert_eq!(r.path.last(), Some(&(r.rows - 1, r.cols - 1)));
        assert_eq!(rescore(&r, s), r.score);
        assert!(r.layers.is_some());
    }

    #[test]
    fn lokal_skor_tutarli_ve_negatif_yok() {
        let s = Scoring::blosum62(-11, -1);
        let r = gotoh_local(b"MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMF", b"MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLV", s, false).unwrap();
        assert!(r.score > 0);
        assert_eq!(rescore(&r, s), r.score);
        assert!(r.matrix.iter().all(|&v| v >= 0));
        assert_eq!(r.path.last(), Some(&r.end));
    }

    #[test]
    fn lokal_benzerlik_yoksa_bos() {
        let s = Scoring::with_gaps(Substitution::Simple { match_score: 1, mismatch: -1 }, -3, -1);
        let r = gotoh_local(b"AAAA", b"TTTT", s, false).unwrap();
        assert_eq!(r.score, 0);
        assert!(r.path.is_empty());
    }

    #[test]
    fn biopython_referans_skorlari() {
        let cases: &[Case] = &REFERENCE;
        for &(a, b, open, ext, local, expected) in cases {
            let s = Scoring::blosum62(open, ext);
            let r = if local { gotoh_local(a, b, s, false) } else { gotoh_global(a, b, s, false) }.unwrap();
            assert_eq!(r.score, expected, "{:?} {:?} local={local}", std::str::from_utf8(a), std::str::from_utf8(b));
            assert_eq!(rescore(&r, s), expected);
        }
    }

    const REFERENCE: [Case; 8] = [
        (b"HEAGAWGHEE", b"PAWHEAE", -11, -1, false, 1),
        (b"HEAGAWGHEE", b"PAWHEAE", -11, -1, true, 17),
        (b"MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTK", b"MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQ", -11, -1, false, 74),
        (b"MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTK", b"MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQ", -11, -1, true, 77),
        (b"MKTAYIAKQRQISFVKSHFSRQ", b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ", -10, -2, false, 77),
        (b"MKTAYIAKQRQISFVKSHFSRQ", b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ", -10, -2, true, 109),
        (b"WWWWCCCC", b"CCCCWWWW", -5, -1, false, 26),
        (b"WWWWCCCC", b"CCCCWWWW", -5, -1, true, 44),
    ];
}