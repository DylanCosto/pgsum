//! SHA-256 digests of input files, written as `sha256:<hex>`.
//!
//! Hashing a 3 GB reference takes over a second, so `cached_file_sha256` remembers digests in
//! `digests.tsv` in the cache directory (`$PGSUM_CACHE_DIR`, else `$XDG_CACHE_HOME/pgsum`, else
//! `~/.cache/pgsum`), keyed by the file's canonical path, size, modification time and (on Unix) device and
//! inode. A file that changes in any of these is hashed again. `PGSUM_NO_DIGEST_CACHE=1` disables the cache.

use std::fs::File;
use std::io::{self, Read};
use std::path::{Path, PathBuf};

use sha2::{Digest, Sha256};

use crate::{Error, Result};

pub fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

pub fn file_sha256(path: &Path) -> Result<String> {
    let mut reader = HashingReader::new(File::open(path).map_err(Error::io(path))?);
    io::copy(&mut reader, &mut io::sink()).map_err(Error::io(path))?;
    Ok(reader.finish())
}

/// Hashes every byte read through it.
pub struct HashingReader<R> {
    inner: R,
    hasher: Sha256,
    pub bytes: u64,
}

impl<R: Read> HashingReader<R> {
    pub fn new(inner: R) -> Self {
        HashingReader {
            inner,
            hasher: Sha256::new(),
            bytes: 0,
        }
    }

    pub fn finish(self) -> String {
        format!("sha256:{}", hex(&self.hasher.finalize()))
    }
}

impl<R: Read> Read for HashingReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.inner.read(buf)?;
        self.hasher.update(&buf[..n]);
        self.bytes += n as u64;
        Ok(n)
    }
}

fn cache_dir() -> Option<PathBuf> {
    if std::env::var_os("PGSUM_NO_DIGEST_CACHE").is_some_and(|v| !v.is_empty() && v != "0") {
        return None;
    }
    let from = |var: &str| std::env::var_os(var).filter(|v| !v.is_empty()).map(PathBuf::from);
    from("PGSUM_CACHE_DIR")
        .or_else(|| from("XDG_CACHE_HOME").map(|d| d.join("pgsum")))
        .or_else(|| from("HOME").map(|d| d.join(".cache/pgsum")))
}

/// The fingerprint a cached digest is valid for: canonical path, size, modification time in nanoseconds
/// and, on Unix, device and inode.
fn fingerprint(path: &Path) -> Option<String> {
    let canonical = std::fs::canonicalize(path).ok()?;
    let meta = std::fs::metadata(&canonical).ok()?;
    let modified = meta
        .modified()
        .ok()?
        .duration_since(std::time::UNIX_EPOCH)
        .ok()?
        .as_nanos();
    #[cfg(unix)]
    let node = {
        use std::os::unix::fs::MetadataExt;
        format!("{}:{}", meta.dev(), meta.ino())
    };
    #[cfg(not(unix))]
    let node = String::from("-");
    let path = canonical.to_str()?;
    (!path.contains(['\t', '\n'])).then(|| format!("{path}\t{}\t{modified}\t{node}", meta.len()))
}

/// `file_sha256`, remembered across runs while the file's fingerprint is unchanged (see the module notes).
/// Any problem with the cache falls back to hashing.
pub fn cached_file_sha256(path: &Path) -> Result<String> {
    let (Some(dir), Some(key)) = (cache_dir(), fingerprint(path)) else {
        return file_sha256(path);
    };
    let text = std::fs::read_to_string(dir.join("digests.tsv")).unwrap_or_default();
    let cached = text.lines().find_map(|l| {
        l.rsplit_once('\t')
            .filter(|(k, _)| *k == key)
            .map(|(_, v)| v.to_owned())
    });
    if let Some(digest) = cached.filter(|d| d.starts_with("sha256:")) {
        return Ok(digest);
    }
    let digest = file_sha256(path)?;
    store(&dir, &key, &text, &digest);
    Ok(digest)
}

/// Record a digest computed while reading `path` some other way (e.g. while scanning it), so a later
/// `cached_file_sha256` needn't read the file again.
pub fn remember_file_sha256(path: &Path, digest: &str) {
    if let (Some(dir), Some(key)) = (cache_dir(), fingerprint(path)) {
        let text = std::fs::read_to_string(dir.join("digests.tsv")).unwrap_or_default();
        store(&dir, &key, &text, digest);
    }
}

fn store(dir: &Path, key: &str, text: &str, digest: &str) {
    let cache = dir.join("digests.tsv");
    // Keep other files' entries, dropping older ones for this path and any for files that no longer exist,
    // and replace the cache atomically.
    let path_prefix = format!("{}\t", key.split('\t').next().unwrap_or_default());
    let mut kept: String = text
        .lines()
        .filter(|l| !l.starts_with(&path_prefix))
        .filter(|l| l.split('\t').next().is_some_and(|p| Path::new(p).exists()))
        .map(|l| format!("{l}\n"))
        .collect();
    kept.push_str(&format!("{key}\t{digest}\n"));
    let tmp = dir.join(format!("digests.tsv.{}.tmp", std::process::id()));
    let _ = std::fs::create_dir_all(dir)
        .and_then(|_| std::fs::write(&tmp, kept))
        .and_then(|_| std::fs::rename(&tmp, &cache))
        .inspect_err(|_| {
            let _ = std::fs::remove_file(&tmp);
        });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cached_digest_follows_the_file() {
        let dir = std::env::temp_dir().join(format!("pgsum-digest-test-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        // SAFETY: tests in this module are the only readers of PGSUM_CACHE_DIR.
        unsafe { std::env::set_var("PGSUM_CACHE_DIR", dir.join("cache")) };
        let file = dir.join("data.txt");
        std::fs::write(&file, "one").unwrap();
        let first = cached_file_sha256(&file).unwrap();
        assert_eq!(first, file_sha256(&file).unwrap());
        assert_eq!(cached_file_sha256(&file).unwrap(), first);
        assert!(
            std::fs::read_to_string(dir.join("cache/digests.tsv"))
                .unwrap()
                .contains(&first)
        );
        std::thread::sleep(std::time::Duration::from_millis(20));
        std::fs::write(&file, "two, longer").unwrap();
        let second = cached_file_sha256(&file).unwrap();
        assert_eq!(second, file_sha256(&file).unwrap());
        assert_ne!(first, second);
        assert_eq!(
            std::fs::read_to_string(dir.join("cache/digests.tsv"))
                .unwrap()
                .lines()
                .count(),
            1
        );
        // Entries for deleted files are dropped when the cache is next written.
        let other = dir.join("other.txt");
        std::fs::write(&other, "three").unwrap();
        cached_file_sha256(&other).unwrap();
        std::fs::remove_file(&file).unwrap();
        let third = dir.join("third.txt");
        std::fs::write(&third, "four").unwrap();
        cached_file_sha256(&third).unwrap();
        let cache = std::fs::read_to_string(dir.join("cache/digests.tsv")).unwrap();
        assert_eq!(cache.lines().count(), 2, "{cache}");
        std::fs::remove_dir_all(&dir).unwrap();
    }
}

/// MD5 (RFC 1321), only to compare files with published MD5 manifests; never used as an identity on its own.
#[derive(Clone)]
pub struct Md5 {
    state: [u32; 4],
    buffer: Vec<u8>,
    length: u64,
}

impl Default for Md5 {
    fn default() -> Self {
        Md5 {
            state: [0x67452301, 0xefcdab89, 0x98badcfe, 0x10325476],
            buffer: Vec::with_capacity(64),
            length: 0,
        }
    }
}

const MD5_SHIFTS: [u32; 64] = [
    7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 5, 9, 14, 20, 5, 9, 14, 20, 5, 9, 14, 20, 5, 9, 14, 20,
    4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23, 6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15,
    21,
];

impl Md5 {
    pub fn update(&mut self, mut data: &[u8]) {
        self.length = self.length.wrapping_add(data.len() as u64);
        if !self.buffer.is_empty() {
            let take = (64 - self.buffer.len()).min(data.len());
            self.buffer.extend_from_slice(&data[..take]);
            data = &data[take..];
            if self.buffer.len() == 64 {
                let block: [u8; 64] = self.buffer[..].try_into().expect("64 bytes");
                self.block(&block);
                self.buffer.clear();
            }
        }
        let mut chunks = data.chunks_exact(64);
        for chunk in &mut chunks {
            self.block(chunk.try_into().expect("64 bytes"));
        }
        self.buffer.extend_from_slice(chunks.remainder());
    }

    fn block(&mut self, block: &[u8; 64]) {
        let m: Vec<u32> = block
            .chunks_exact(4)
            .map(|c| u32::from_le_bytes(c.try_into().expect("4")))
            .collect();
        let [mut a, mut b, mut c, mut d] = self.state;
        for i in 0..64 {
            let (f, g) = match i / 16 {
                0 => ((b & c) | (!b & d), i),
                1 => ((d & b) | (!d & c), (5 * i + 1) % 16),
                2 => (b ^ c ^ d, (3 * i + 5) % 16),
                _ => (c ^ (b | !d), (7 * i) % 16),
            };
            let k = ((i as f64 + 1.0).sin().abs() * 4294967296.0) as u32;
            let rotated = a
                .wrapping_add(f)
                .wrapping_add(k)
                .wrapping_add(m[g])
                .rotate_left(MD5_SHIFTS[i]);
            (a, d, c) = (d, c, b);
            b = b.wrapping_add(rotated);
        }
        for (s, v) in self.state.iter_mut().zip([a, b, c, d]) {
            *s = s.wrapping_add(v);
        }
    }

    /// The digest as lower-case hex.
    pub fn finish(mut self) -> String {
        let bits = self.length.wrapping_mul(8);
        let mut tail = vec![0x80u8];
        while (self.buffer.len() + tail.len()) % 64 != 56 {
            tail.push(0);
        }
        tail.extend_from_slice(&bits.to_le_bytes());
        let length = self.length;
        self.update(&tail);
        self.length = length;
        hex(&self.state.iter().flat_map(|s| s.to_le_bytes()).collect::<Vec<u8>>())
    }
}

#[cfg(test)]
mod md5_tests {
    use super::Md5;

    #[test]
    fn rfc_1321_test_suite() {
        for (input, expected) in [
            ("", "d41d8cd98f00b204e9800998ecf8427e"),
            ("a", "0cc175b9c0f1b6a831c399e269772661"),
            ("abc", "900150983cd24fb0d6963f7d28e17f72"),
            ("message digest", "f96b697d7cb7938d525a2f31aaf161d0"),
            ("abcdefghijklmnopqrstuvwxyz", "c3fcd3d76192e4007dfb496cca67e13b"),
            (
                "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789",
                "d174ab98d277d9f5a5611c2c9f419d9f",
            ),
            (
                "12345678901234567890123456789012345678901234567890123456789012345678901234567890",
                "57edf4a22be3c955ac49da2e2107b67a",
            ),
        ] {
            let mut whole = Md5::default();
            whole.update(input.as_bytes());
            assert_eq!(whole.finish(), expected, "{input:?}");
            // The same digest when fed in uneven pieces.
            let mut pieces = Md5::default();
            for chunk in input.as_bytes().chunks(7) {
                pieces.update(chunk);
            }
            assert_eq!(pieces.finish(), expected);
        }
    }
}
