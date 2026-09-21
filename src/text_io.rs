use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;

fn is_gzipped(path: &str) -> bool {
    Path::new(path)
        .extension()
        .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"))
}

pub fn open_input(path: Option<&str>) -> Result<Box<dyn Read + Send>> {
    match path {
        Some(path) => {
            let file =
                File::open(path).with_context(|| format!("failed to open input file: {path}"))?;
            if is_gzipped(path) {
                Ok(Box::new(MultiGzDecoder::new(file)))
            } else {
                Ok(Box::new(file))
            }
        }
        None => Ok(Box::new(io::stdin())),
    }
}

pub fn open_optional_input(path: Option<&str>) -> Result<Option<Box<dyn Read + Send>>> {
    path.map(|path| open_input(Some(path))).transpose()
}

/// Keeps compression ownership until the command explicitly checks finalization.
pub enum Output {
    Plain(BufWriter<Box<dyn Write>>),
    Gzip(Box<GzEncoder<BufWriter<Box<dyn Write>>>>),
}

impl Output {
    pub fn borrowed(&mut self) -> Box<dyn Write + '_> {
        Box::new(self)
    }

    fn finish(self) -> io::Result<()> {
        let mut writer = match self {
            Self::Plain(writer) => writer,
            Self::Gzip(writer) => writer.finish()?,
        };
        writer.flush()
    }
}

impl Write for Output {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        match self {
            Self::Plain(writer) => writer.write(bytes),
            Self::Gzip(writer) => writer.write(bytes),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            Self::Plain(writer) => writer.flush(),
            Self::Gzip(writer) => writer.flush(),
        }
    }
}

pub fn open_output(path: Option<&str>, compresslevel: u32) -> Result<Output> {
    let destination: Box<dyn Write> = match path {
        Some(path) => Box::new(
            File::create(path).with_context(|| format!("failed to create output file: {path}"))?,
        ),
        None => Box::new(io::stdout().lock()),
    };
    let writer = BufWriter::new(destination);
    Ok(if path.is_some_and(is_gzipped) {
        Output::Gzip(Box::new(GzEncoder::new(
            writer,
            Compression::new(compresslevel),
        )))
    } else {
        Output::Plain(writer)
    })
}

pub fn finish_outputs(outputs: impl IntoIterator<Item = Option<Output>>) -> io::Result<()> {
    for output in outputs.into_iter().flatten() {
        output.finish()?;
    }
    Ok(())
}
