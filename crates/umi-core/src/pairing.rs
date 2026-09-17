//! Paired-end read handling before grouping, mirroring the start of
//! `umi_tools.sam_methods.get_bundles.__call__`.

use rust_htslib::bam::Record;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairPolicy {
    Discard,
    Use,
    /// Pass the read through ungrouped; only `group` can do this.
    Output,
}

impl PairPolicy {
    #[must_use]
    pub fn parse(name: &str) -> Option<Self> {
        match name {
            "discard" => Some(Self::Discard),
            "use" => Some(Self::Use),
            "output" => Some(Self::Output),
            _ => None,
        }
    }
}

/// `--paired` and how unmapped mates, unpaired reads and chimeric pairs are treated.
#[derive(Debug, Clone, Copy)]
pub struct PairingOptions {
    pub paired: bool,
    pub unmapped_reads: PairPolicy,
    pub chimeric_pairs: PairPolicy,
    pub unpaired_reads: PairPolicy,
}

impl Default for PairingOptions {
    fn default() -> Self {
        Self {
            paired: false,
            unmapped_reads: PairPolicy::Discard,
            chimeric_pairs: PairPolicy::Use,
            unpaired_reads: PairPolicy::Use,
        }
    }
}

#[derive(Debug, thiserror::Error)]
pub enum PairingError {
    #[error(
        "Cannot use --{0}=output. If you want to retain these reads without deduplicating them, use the group command"
    )]
    OutputNotAllowed(&'static str),
}

/// What `get_bundles` does with one read before grouping.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Triage {
    /// Whether the read is a read2, which never counts as input.
    pub is_read2: bool,
    /// Ungrouped copies to write. Two are possible: `umi_tools` yields an
    /// unpaired read for `--unpaired-reads=output` and then keeps processing it.
    pub copies: u8,
    /// Whether the read goes on to grouping.
    pub grouped: bool,
}

impl PairingOptions {
    /// `validateSamOptions` for commands that cannot pass reads through ungrouped.
    ///
    /// # Errors
    ///
    /// Returns the option set to `output`.
    pub fn validate(&self, can_output: bool) -> Result<(), PairingError> {
        if can_output {
            return Ok(());
        }
        for (name, policy) in [
            ("unmapped-reads", self.unmapped_reads),
            ("chimeric-pairs", self.chimeric_pairs),
            ("unpaired-reads", self.unpaired_reads),
        ] {
            if policy == PairPolicy::Output {
                return Err(PairingError::OutputNotAllowed(name));
            }
        }
        Ok(())
    }

    /// Whether unmapped reads (and unmapped read2s) are passed through.
    #[must_use]
    pub const fn outputs_unmapped(&self) -> bool {
        matches!(self.unmapped_reads, PairPolicy::Use | PairPolicy::Output)
    }

    /// `passthrough` is true for `group`, which writes ungrouped reads.
    #[must_use]
    pub fn triage(&self, record: &Record, passthrough: bool) -> Triage {
        let outputs_unmapped = passthrough && self.outputs_unmapped();
        let skipped = |copies| Triage {
            is_read2: false,
            copies,
            grouped: false,
        };
        let mut copies = 0;

        if record.is_last_in_template() {
            let written = passthrough && (!record.is_unmapped() || outputs_unmapped);
            return Triage {
                is_read2: true,
                copies: u8::from(written),
                grouped: false,
            };
        }

        if self.paired && !record.is_paired() {
            match self.unpaired_reads {
                PairPolicy::Discard => return skipped(0),
                PairPolicy::Output => copies += 1,
                PairPolicy::Use => {}
            }
        }

        if record.is_unmapped() {
            return skipped(copies + u8::from(outputs_unmapped));
        }

        if self.paired && record.is_mate_unmapped() && self.unmapped_reads != PairPolicy::Use {
            return skipped(copies + u8::from(outputs_unmapped));
        }

        if record.is_paired() && record.tid() != record.mtid() {
            match self.chimeric_pairs {
                PairPolicy::Discard => return skipped(copies),
                PairPolicy::Output => return skipped(copies + 1),
                PairPolicy::Use => {}
            }
        }

        Triage {
            is_read2: false,
            copies,
            grouped: true,
        }
    }
}

#[cfg(test)]
mod tests {
    use rust_htslib::bam::HeaderView;

    use super::*;

    fn read(flag: u16, mate_contig: &str) -> Record {
        let header = HeaderView::from_bytes(b"@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1000\n");
        let line = format!("r\t{flag}\tchr1\t101\t60\t10M\t{mate_contig}\t201\t100\t*\t*");
        Record::from_sam(&header, line.as_bytes()).unwrap()
    }

    const PAIRED_R1: u16 = 0x1 | 0x40;
    const PAIRED_R2: u16 = 0x1 | 0x80;

    fn paired(unmapped: PairPolicy, chimeric: PairPolicy, unpaired: PairPolicy) -> PairingOptions {
        PairingOptions {
            paired: true,
            unmapped_reads: unmapped,
            chimeric_pairs: chimeric,
            unpaired_reads: unpaired,
        }
    }

    #[test]
    fn read2_is_written_only_by_group() {
        let options = PairingOptions::default();
        let mapped_r2 = read(PAIRED_R2, "=");
        assert_eq!(options.triage(&mapped_r2, true).copies, 1);
        assert_eq!(options.triage(&mapped_r2, false).copies, 0);
        assert!(options.triage(&mapped_r2, false).is_read2);

        let unmapped_r2 = read(PAIRED_R2 | 0x4, "=");
        assert_eq!(options.triage(&unmapped_r2, true).copies, 0);
        let keep_unmapped = paired(PairPolicy::Use, PairPolicy::Use, PairPolicy::Use);
        assert_eq!(keep_unmapped.triage(&unmapped_r2, true).copies, 1);
    }

    #[test]
    fn mate_unmapped_follows_unmapped_reads() {
        let r1 = read(PAIRED_R1 | 0x8, "=");
        let discard = PairingOptions {
            paired: true,
            ..PairingOptions::default()
        };
        assert!(!discard.triage(&r1, false).grouped);
        let use_it = paired(PairPolicy::Use, PairPolicy::Use, PairPolicy::Use);
        assert!(use_it.triage(&r1, false).grouped);
        let output = paired(PairPolicy::Output, PairPolicy::Use, PairPolicy::Use);
        let triage = output.triage(&r1, true);
        assert_eq!((triage.copies, triage.grouped), (1, false));
    }

    #[test]
    fn chimeric_pairs_are_detected_without_the_paired_option() {
        let r1 = read(PAIRED_R1, "chr2");
        let discard = PairingOptions {
            chimeric_pairs: PairPolicy::Discard,
            ..PairingOptions::default()
        };
        assert!(!discard.triage(&r1, false).grouped);
        assert!(PairingOptions::default().triage(&r1, false).grouped);
    }

    #[test]
    fn unpaired_output_is_written_and_still_grouped() {
        let single = read(0, "*");
        let output = paired(PairPolicy::Discard, PairPolicy::Use, PairPolicy::Output);
        let triage = output.triage(&single, true);
        assert_eq!((triage.copies, triage.grouped), (1, true));
        let discard = paired(PairPolicy::Discard, PairPolicy::Use, PairPolicy::Discard);
        assert!(!discard.triage(&single, true).grouped);
        assert!(PairingOptions::default().triage(&single, true).grouped);
    }

    #[test]
    fn output_is_only_for_group() {
        let output = paired(PairPolicy::Output, PairPolicy::Use, PairPolicy::Use);
        assert!(output.validate(true).is_ok());
        assert!(matches!(
            output.validate(false),
            Err(PairingError::OutputNotAllowed("unmapped-reads"))
        ));
    }
}
