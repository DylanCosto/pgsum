//! Native BCF decoding into bounded in-memory VCF records for the shared genotype rules.
//! No temporary converted file is written; hashing stays on the original compressed input.
use crate::{Error, Result};
use noodles_vcf::variant::io::Write as _;
use std::io::{self, BufRead, Read};
use std::path::Path;

/// Detect the decompressed magic rather than relying on a filename extension.
pub fn is_bcf(path: &Path) -> Result<bool> {
    let mut input = std::io::BufReader::new(std::fs::File::open(path).map_err(Error::io(path))?);
    let gzip = input.fill_buf().map_err(Error::io(path))?.starts_with(&[0x1f, 0x8b]);
    let mut source: Box<dyn Read> = if gzip {
        Box::new(flate2::read::MultiGzDecoder::new(input))
    } else {
        Box::new(input)
    };
    let mut magic = [0; 3];
    match source.read_exact(&mut magic) {
        Ok(()) => Ok(magic == *b"BCF"),
        Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => Ok(false),
        Err(e) => Err(Error::io(path)(e)),
    }
}

pub(crate) struct TextReader<R> {
    reader: noodles_bcf::io::Reader<R>,
    header: noodles_vcf::Header,
    record: noodles_bcf::Record,
    buffer: Vec<u8>,
    position: usize,
    eof: bool,
}

impl<R: Read> TextReader<R> {
    /// `input` is already decompressed; the same upstream reader hashes the original bytes.
    pub fn new(input: R) -> io::Result<Self> {
        let mut reader = noodles_bcf::io::Reader::from(input);
        let header = reader.read_header()?;
        let mut buffer = Vec::new();
        noodles_vcf::io::Writer::new(&mut buffer).write_header(&header)?;
        Ok(Self {
            reader,
            header,
            record: Default::default(),
            buffer,
            position: 0,
            eof: false,
        })
    }
}

impl<R: Read> BufRead for TextReader<R> {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        if self.position == self.buffer.len() && !self.eof {
            self.buffer.clear();
            self.position = 0;
            if self.reader.read_record(&mut self.record)? == 0 {
                self.eof = true;
            } else {
                let samples = self.record.samples()?;
                for series in samples.series() {
                    let series = series?;
                    let name = series.name(&self.header)?;
                    if !self.header.formats().contains_key(name) {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("BCF FORMAT/{name} lacks a type definition"),
                        ));
                    }
                }
                // The codec has unimplemented branches for uncommon encodings. Surface those as
                // an input error, never continue with an incomplete or substituted genotype.
                std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    noodles_vcf::io::Writer::new(&mut self.buffer).write_variant_record(&self.header, &self.record)
                }))
                .map_err(|_| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        "unsupported or malformed BCF field encoding",
                    )
                })??;
            }
        }
        Ok(&self.buffer[self.position..])
    }
    fn consume(&mut self, amount: usize) {
        self.position = (self.position + amount).min(self.buffer.len());
    }
}

impl<R: Read> Read for TextReader<R> {
    fn read(&mut self, output: &mut [u8]) -> io::Result<usize> {
        let available = self.fill_buf()?;
        let n = output.len().min(available.len());
        output[..n].copy_from_slice(&available[..n]);
        self.consume(n);
        Ok(n)
    }
}
