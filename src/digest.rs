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
