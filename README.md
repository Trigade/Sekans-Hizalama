# Sekans Hizalama Görselleştirici · Rust + WebAssembly

Needleman–Wunsch (global) ve Smith–Waterman (lokal) sekans hizalama algoritmalarını **tarayıcıda, adım adım** görselleştiren interaktif bir araç. Hesaplama çekirdeği Rust ile yazılmış, WebAssembly'ye derlenmiştir. Kurulum gerektirmez; tüm hesaplama kullanıcının tarayıcısında yapılır.

**🔗 Canlı demo: https://KULLANICI_ADI.github.io/seq-align-wasm/**

<!-- Buraya demo GIF'i ekle: ![Demo](docs/demo.gif) -->

> Biyoenformatik Algoritmalar dersi kapsamında geliştirilmiştir.

---

## Özellikler

- **İki algoritma:** Needleman–Wunsch (global) ve Smith–Waterman (lokal) hizalama
- **Animasyonlu DP matrisi:** Hücreler anti-diyagonal dalga cephesi sırasıyla dolar, ardından traceback yolu vurgulanır
- **Hücre açıklaması:** Herhangi bir hücrenin üzerine gelindiğinde değerinin hangi komşudan, hangi formülle hesaplandığı gösterilir
- **Ayarlanabilir skorlama:** Eşleşme, uyuşmazlık ve boşluk cezası; değişiklikler anında yeniden hesaplanır
- **FASTA desteği:** `>` ile başlayan başlık satırları ve boşluklar otomatik temizlenir
- **Hizalama istatistikleri:** Skor, özdeşlik yüzdesi, eşleşme / uyuşmazlık / boşluk sayıları
- **Benchmark:** Aynı algoritmanın saf JavaScript ve Rust/WASM sürümlerini MCUPS cinsinden karşılaştırır
- **Hazır örnekler:** Wikipedia örnekleri, Durbin protein örneği, gömülü EcoRI motifi (`GAATTC`), insan hemoglobin α/β zincirleri

## Algoritmalar

Her iki algoritma da dinamik programlama ile `(n+1) × (m+1)` boyutlu bir skor matrisi doldurur:

```
H(i, j) = max( H(i-1, j-1) + s(aᵢ, bⱼ),   ↖ eşleşme / uyuşmazlık
               H(i-1, j)   + gap,          ↑ A'da harf, B'de boşluk
               H(i, j-1)   + gap )         ← A'da boşluk, B'de harf
```

| | Needleman–Wunsch (global) | Smith–Waterman (lokal) |
|---|---|---|
| Amaç | İki dizinin tamamını hizalar | En benzer alt bölgeyi bulur |
| Kenar değerleri | `i × gap`, `j × gap` | `0` |
| Hücre formülü | yukarıdaki `max` | `max(0, …)` — negatif değer olmaz |
| Traceback başlangıcı | Sağ alt köşe | Matrisin en yüksek hücresi |
| Traceback bitişi | `(0, 0)` | İlk `0` değerli hücre |

Eşitlik durumunda yön önceliği `↖ > ↑ > ←` olarak sabittir; bu, sonuçların deterministik ve test edilebilir olmasını sağlar.

**Karmaşıklık:** Tam hizalama O(n·m) zaman ve bellek kullanır. Yalnızca skor hesaplayan sürümler (`*_score`) iki satır tutarak **O(m) bellek** ile çalışır.

## Mimari

```
┌──────────────────────────┐      ┌────────────────────┐      ┌──────────────────────┐
│  Rust çekirdek (core)    │      │  WASM köprüsü      │      │  Arayüz (web)        │
│  scoring.rs · nw.rs ·    │ ───▶ │  lib.rs            │ ───▶ │  main.js · animate.js│
│  sw.rs                   │      │  wasm-bindgen      │      │  bench.js            │
│  saf Rust, WASM'dan      │      │  normalizasyon,    │      │  vanilla JS,         │
│  bağımsız                │      │  limitler, hatalar │      │  framework yok       │
└──────────────────────────┘      └────────────────────┘      └──────────────────────┘
```

Tasarım kararları:

- **Katmanlı çekirdek:** Algoritmalar `wasm_bindgen`'e bağımlı değildir. Aynı kod `cargo test` ile native olarak test edilir ve ileride bir masaüstü uygulamasında (Tauri) değiştirilmeden kullanılabilir.
- **Düz matris:** `Vec<Vec<i32>>` yerine tek blok `Vec<i32>` (`index = i * cols + j`). Satır başına ayrı bellek tahsisi yoktur, önbellek dostudur ve JS tarafına tek dizi olarak aktarılır.
- **1 byte'lık yön kodu:** Traceback yönleri `#[repr(u8)]` enum olarak tutulur ve JS'e string yerine `0–3` sayısı olarak gönderilir.
- **Panik yerine `Result`:** Geçersiz karakter, mantıksız skorlama veya aşırı uzun dizi durumunda JS tarafına anlamlı bir `Error` fırlatılır.
- **Performanslı animasyon:** DOM bir kez oluşturulur, animasyon yalnızca CSS sınıfı ekler (layout tetiklemeyen `opacity` / `transform` geçişleri). Tüm animasyonlar `AbortController` ile iptal edilebilir.

## Proje yapısı

```
seq-align-wasm/
├── core/                    # Rust kütüphanesi
│   ├── Cargo.toml
│   └── src/
│       ├── lib.rs           # WASM'a açılan fonksiyonlar
│       ├── scoring.rs       # Skorlama, hata tipleri, sonuç yapıları
│       ├── nw.rs            # Needleman–Wunsch + birim testleri
│       └── sw.rs            # Smith–Waterman + birim testleri
├── web/
│   ├── index.html
│   ├── style.css
│   ├── main.js              # Uygulama akışı, WASM çağrıları
│   ├── animate.js           # Matris çizimi ve animasyon
│   └── bench.js             # JS referans uygulaması + benchmark
└── .github/workflows/
    └── deploy.yml           # Test → WASM derleme → GitHub Pages
```

## Yerelde çalıştırma

Gereksinimler: [Rust](https://rustup.rs), `wasm32-unknown-unknown` hedefi ve [wasm-pack](https://github.com/drager/wasm-pack).

```bash
rustup target add wasm32-unknown-unknown
cargo install wasm-pack

cd core
cargo test
wasm-pack build --target web --release --out-dir ../web/pkg

cd ../web
python -m http.server 8000
```

Ardından tarayıcıda `http://localhost:8000` adresini aç.

> Tarayıcılar `.wasm` dosyalarını `file://` ile açılan sayfalara yüklemez; bu yüzden basit bir HTTP sunucusu gerekir.

## Testler

```bash
cd core
cargo test
```

Testler; Wikipedia'daki bilinen örnekleri, kenar durumlarını (boş diziler, ortak harf olmayan diziler), hizalamadan yeniden hesaplanan skorun matris skoruyla eşitliğini, lineer bellekli sürümlerin tam sürümle aynı skoru vermesini ve geçersiz girdilerin reddedilmesini doğrular.

## Benchmark

Rastgele DNA dizileri (sabit tohumlu PRNG), her boyut için 1 ısınma + 5 ölçüm, medyan süre. Her ölçümde iki tarafın skorunun eşit olduğu doğrulanır.

<!-- Kendi ölçümlerini buraya yaz -->

| Boyut | Hücre | JavaScript | Rust/WASM | Hızlanma |
|---|---|---|---|---|
| 500 × 500 | 0,25 M | – ms | – ms | –× |
| 1.000 × 1.000 | 1 M | – ms | – ms | –× |
| 2.000 × 2.000 | 4 M | – ms | – ms | –× |
| 4.000 × 4.000 | 16 M | – ms | – ms | –× |

## Yol haritası

- [ ] Rust tarafında sınır kontrolsüz döngüler ve WASM SIMD ile optimizasyon
- [ ] Affine gap cezası (Gotoh algoritması)
- [ ] BLOSUM62 ikame matrisi (protein hizalama)
- [ ] Hirschberg algoritması ile lineer bellekli tam hizalama
- [ ] k-mer indeksli BLAST benzeri veri tabanı araması
- [ ] Tauri ile masaüstü sürümü

## Kaynaklar

- Needleman, S. B. & Wunsch, C. D. (1970). *A general method applicable to the search for similarities in the amino acid sequence of two proteins.* Journal of Molecular Biology, 48(3), 443–453.
- Smith, T. F. & Waterman, M. S. (1981). *Identification of common molecular subsequences.* Journal of Molecular Biology, 147(1), 195–197.
- Durbin, R., Eddy, S., Krogh, A. & Mitchison, G. (1998). *Biological Sequence Analysis.* Cambridge University Press.