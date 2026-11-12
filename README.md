# Sekans Hizalama Görselleştirici · Rust + WebAssembly

Needleman–Wunsch, Smith–Waterman, Gotoh ve Hirschberg sekans hizalama algoritmalarını **tarayıcıda, adım adım** görselleştiren interaktif bir araç. Hesaplama çekirdeği Rust ile yazılmış, WebAssembly'ye derlenmiştir. Kurulum gerektirmez; tüm hesaplama kullanıcının tarayıcısında yapılır.

**🔗 Canlı demo: https://trigade.github.io/Sekans-Hizalama/**

<!-- Buraya demo GIF'i ekle: ![Demo](docs/demo.gif) -->

> Biyoenformatik Algoritmalar dersi kapsamında geliştirilmiştir.

---

## Özellikler

- **Dört algoritma:** Needleman–Wunsch (global), Smith–Waterman (lokal), Gotoh (afin boşluk, global ve lokal), Hirschberg (lineer bellekli global)
- **BLOSUM62:** Protein dizileri için standart ikame matrisi
- **Otomatik algoritma seçimi:** Seçilen ayarlara ve dizi boyutuna göre doğru algoritma kullanılır
- **Animasyonlu DP matrisi:** Hücreler anti-diyagonal dalga cephesi sırasıyla dolar, ardından traceback yolu vurgulanır
- **Hücre açıklaması:** Herhangi bir hücrenin üzerine gelindiğinde değerinin hangi komşudan, hangi formülle hesaplandığı gösterilir
- **Ayarlanabilir skorlama:** Basit eşleşme/uyuşmazlık ya da BLOSUM62; lineer ya da afin (açma + uzatma) boşluk cezası; değişiklikler anında yeniden hesaplanır
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

### Gotoh: afin boşluk cezası

Gerçek biyolojik dizilerde tek bir uzun boşluk, birçok kısa boşluktan daha olasıdır. Afin modelde `k` uzunluklu bir boşluğun maliyeti `açma + k × uzatma` olur. Gotoh (1982) bunu üç matrisle O(n·m) sürede çözer:

```
M(i, j) = s(aᵢ, bⱼ) + max( M, X, Y )(i-1, j-1)                    harf–harf
X(i, j) = max( M(i-1, j) + açma + uzatma,  X(i-1, j) + uzatma,
               Y(i-1, j) + açma + uzatma )                         B'de boşluk
Y(i, j) = max( M(i, j-1) + açma + uzatma,  Y(i, j-1) + uzatma,
               X(i, j-1) + açma + uzatma )                         A'da boşluk
```

Arayüzde bir hücrenin üzerine gelindiğinde üç matrisin değerleri ve hangisinin kazandığı gösterilir.

### Hirschberg: lineer bellek

Tam DP matrisi O(n·m) bellek ister: 10.000 harflik iki dizi için 400 MB. Hirschberg (1975), A dizisini ortadan bölüp ileri ve geri yönde yalnızca son satırları hesaplayarak B'deki en iyi kesme noktasını bulur ve iki yarıyı özyinelemeli çözer. Zaman yine O(n·m), bellek ise **O(n + m)**. Görselleştirilemeyecek kadar büyük global hizalamalarda otomatik olarak kullanılır.

| Koşul | Kullanılan algoritma |
|---|---|
| Global, lineer boşluk, matris ≤ 120×120 | Needleman–Wunsch (matris çizilir) |
| Global, lineer boşluk, daha büyük | Hirschberg |
| Lokal, lineer boşluk | Smith–Waterman |
| Afin boşluk | Gotoh (global veya lokal) |

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
Sekans-Hizalama/
├── core/                    # Rust kütüphanesi
│   ├── Cargo.toml
│   └── src/
│       ├── lib.rs           # WASM'a açılan fonksiyonlar
│       ├── scoring.rs       # Skorlama, hata tipleri, sonuç yapıları
│       ├── nw.rs            # Needleman–Wunsch
│       ├── sw.rs            # Smith–Waterman
│       ├── gotoh.rs         # Gotoh (afin boşluk, global + lokal)
│       ├── hirschberg.rs    # Hirschberg (lineer bellek)
│       └── blosum.rs        # BLOSUM62 tablosu
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

32 birim testi; Wikipedia'daki bilinen örnekleri, Biopython `PairwiseAligner` ile hesaplanmış BLOSUM62 + afin referans skorlarını, Hirschberg'in rastgele dizilerde Needleman–Wunsch ile aynı skoru vermesini, kenar durumlarını (boş diziler, ortak harf olmayan diziler), hizalamadan yeniden hesaplanan skorun matris skoruyla eşitliğini, lineer bellekli sürümlerin tam sürümle aynı skoru vermesini ve geçersiz girdilerin reddedilmesini doğrular.

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
- [x] Afin boşluk cezası (Gotoh algoritması)
- [x] BLOSUM62 ikame matrisi (protein hizalama)
- [x] Hirschberg algoritması ile lineer bellekli tam hizalama
- [ ] k-mer indeksli BLAST benzeri veri tabanı araması
- [ ] Tauri ile masaüstü sürümü

## Kaynaklar

- Needleman, S. B. & Wunsch, C. D. (1970). *A general method applicable to the search for similarities in the amino acid sequence of two proteins.* Journal of Molecular Biology, 48(3), 443–453.
- Smith, T. F. & Waterman, M. S. (1981). *Identification of common molecular subsequences.* Journal of Molecular Biology, 147(1), 195–197.
- Gotoh, O. (1982). *An improved algorithm for matching biological sequences.* Journal of Molecular Biology, 162(3), 705–708.
- Hirschberg, D. S. (1975). *A linear space algorithm for computing maximal common subsequences.* Communications of the ACM, 18(6), 341–343.
- Henikoff, S. & Henikoff, J. G. (1992). *Amino acid substitution matrices from protein blocks.* PNAS, 89(22), 10915–10919.
- Durbin, R., Eddy, S., Krogh, A. & Mitchison, G. (1998). *Biological Sequence Analysis.* Cambridge University Press.