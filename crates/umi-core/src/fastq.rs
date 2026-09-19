use std::io::Write;

use crate::error::ExtractError;

pub fn write_fastq_record<W: Write>(
    writer: &mut W,
    id: &[u8],
    seq: &[u8],
    qual: &[u8],
) -> Result<(), ExtractError> {
    writer.write_all(b"@")?;
    writer.write_all(id)?;
    writer.write_all(b"\n")?;
    writer.write_all(seq)?;
    writer.write_all(b"\n+\n")?;
    writer.write_all(qual)?;
    writer.write_all(b"\n")?;
    Ok(())
}
