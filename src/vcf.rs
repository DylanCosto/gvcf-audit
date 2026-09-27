use crate::{Args, Result, State, fail, reference::Reference};
use serde::Serialize;
use std::collections::{BTreeMap, HashSet};

#[derive(Default, Serialize)]
pub struct Header {
    pub sample: String,
    pub caller: String,
    pub sample_count: usize,
    pub declared_contigs: usize,
    pub unmapped_header_contigs: Vec<String>,
    pub reblocked: bool,
    #[serde(skip)]
    pub column: usize,
    #[serde(skip)]
    pub columns: usize,
    #[serde(skip)]
    pub refcall: bool,
    #[serde(skip)]
    pub callers: HashSet<String>,
    #[serde(skip)]
    contigs: HashSet<usize>,
    #[serde(skip)]
    contig_names: HashSet<String>,
}

impl Header {
    pub fn line(&mut self, line: &str, args: &Args, reference: &Reference) -> Result<bool> {
        let lower = line.to_ascii_lowercase();
        if lower.starts_with("##deepvariant_version=") || lower.starts_with("##source=deepvariant")
        {
            self.callers.insert("deepvariant".into());
        }
        if lower.starts_with("##source=dragen")
            || lower.starts_with("##dragenversion=")
            || lower.starts_with("##dragencommandline=")
        {
            self.callers.insert("dragen".into());
        }
        if lower.starts_with("##gatkcommandline=<id=haplotypecaller,")
            || lower.starts_with("##gatkcommandline.haplotypecaller=")
        {
            self.callers.insert("gatk".into());
        }
        if lower.starts_with("##gatkcommandline=<id=reblockgvcf,") {
            self.reblocked = true;
            self.callers.insert("gatk".into());
        }
        if line.starts_with("##FILTER=<ID=RefCall,")
            && line.contains("Genotyping model thinks this site is reference")
        {
            self.refcall = true;
            self.callers.insert("deepvariant".into());
        }
        if let Some(text) = line
            .strip_prefix("##contig=<")
            .and_then(|s| s.strip_suffix('>'))
        {
            let value = |key: &str| {
                text.split(',')
                    .find_map(|f| f.strip_prefix(key))
                    .map(|v| v.trim_matches('"'))
            };
            let name = value("ID=").ok_or("Contig header has no ID")?;
            if name.is_empty() || !self.contig_names.insert(name.into()) {
                return fail("Empty or duplicate contig header ID");
            }
            let length = value("length=")
                .or_else(|| value("Length="))
                .map(str::parse::<u64>)
                .transpose()?;
            if let Some(id) = reference.resolve_optional(name)? {
                if !self.contigs.insert(id) {
                    return fail("Duplicate or aliased duplicate contig header");
                }
                if length.is_some_and(|n| n != reference.contigs[id].length) {
                    return fail(format!("Reference length mismatch for {name}"));
                }
            } else {
                self.unmapped_header_contigs.push(name.into());
            }
            self.declared_contigs += 1;
        }
        if line.starts_with("#CHROM\t") {
            let f: Vec<_> = line.split('\t').collect();
            if f.len() < 10
                || f[..9]
                    != [
                        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                    ]
            {
                return fail("Expected a VCF header with FORMAT and sample columns");
            }
            let samples = &f[9..];
            let unique: HashSet<_> = samples.iter().collect();
            if unique.len() != samples.len() || samples.contains(&"") {
                return fail("Duplicate or empty sample names");
            }
            let selected = match args.sample.as_deref() {
                Some(s) => samples
                    .iter()
                    .position(|&n| n == s)
                    .ok_or("Requested sample is not in the VCF")?,
                None if samples.len() == 1 => 0,
                None => return fail("Multiple samples: choose one with --sample NAME"),
            };
            self.sample = samples[selected].into();
            self.sample_count = samples.len();
            self.column = 9 + selected;
            self.columns = f.len();
            self.caller = if args.caller != "auto" {
                args.caller.clone()
            } else if self.callers.len() == 1 {
                self.callers.iter().next().unwrap().clone()
            } else {
                "generic".into()
            };
            Ok(true)
        } else if line.starts_with("##") {
            Ok(false)
        } else {
            fail("Unexpected line before #CHROM header")
        }
    }
}

pub struct Record<'a> {
    pub chrom: &'a str,
    pub contig: usize,
    pub start: u64,
    pub end: u64,
    pub state: State,
}

fn number(value: Option<&str>, field: &str) -> Result<Option<f64>> {
    match value {
        None | Some(".") => Ok(None),
        Some(s) => {
            let v: f64 = s
                .parse()
                .map_err(|_| format!("Invalid FORMAT/{field}: {s:?}"))?;
            if !v.is_finite() || v < 0.0 {
                return fail(format!("Invalid FORMAT/{field}: {s:?}"));
            }
            Ok(Some(v))
        }
    }
}

pub fn parse<'a>(
    line: &'a str,
    h: &Header,
    args: &Args,
    reference: &Reference,
    adaptations: &mut BTreeMap<String, u64>,
) -> Result<Record<'a>> {
    let f: Vec<_> = line.split('\t').collect();
    if f.len() != h.columns {
        return fail("Record column count differs from #CHROM header");
    }
    let contig = reference.resolve(f[0])?;
    let start = f[1]
        .parse::<u64>()?
        .checked_sub(1)
        .ok_or("VCF POS must be at least 1")?;
    if f[3].is_empty() || !f[3].bytes().all(|b| b.is_ascii_alphabetic()) {
        return fail("Invalid REF allele");
    }
    let ref_end = start
        .checked_add(f[3].len() as u64)
        .ok_or("REF coordinate overflow")?;
    let mut end = None;
    for item in f[7].split(';') {
        if item == "END" {
            return fail("INFO/END has no value");
        }
        if let Some(v) = item.strip_prefix("END=") {
            if end.is_some() {
                return fail("Duplicate INFO/END");
            }
            end = Some(v.parse::<u64>()?);
        }
    }
    let end = end.unwrap_or(ref_end);
    if end < ref_end || end > reference.contigs[contig].length {
        return fail("Record span is outside the reference or shorter than REF");
    }
    let alts: Vec<_> = f[4].split(',').collect();
    if alts.iter().any(|a| a.is_empty()) || (alts.len() > 1 && alts.contains(&".")) {
        return fail("Empty ALT allele or missing ALT mixed with other alleles");
    }
    let block = alts.iter().all(|a| matches!(*a, "." | "<NON_REF>" | "<*>"));
    let formats: Vec<_> = f[8].split(':').collect();
    let values: Vec<_> = f[h.column].split(':').collect();
    if values.len() > formats.len() {
        return fail("More sample fields than FORMAT keys");
    }
    let mut unique = HashSet::new();
    if formats
        .iter()
        .any(|key| key.is_empty() || !unique.insert(*key))
    {
        return fail("Empty or duplicate FORMAT key");
    }
    let value = |key: &str| {
        formats
            .iter()
            .position(|&k| k == key)
            .and_then(|i| values.get(i).copied())
    };
    let gt_text = value("GT").unwrap_or(".");
    if gt_text.contains('/') && gt_text.contains('|') {
        return fail("Mixed phasing separators in GT");
    }
    let gt: Vec<Option<usize>> = gt_text
        .split(['/', '|'])
        .map(|a| {
            if a == "." {
                Ok(None)
            } else {
                let n: usize = a.parse().map_err(|_| format!("Invalid GT allele {a:?}"))?;
                if n > alts.len() || (n > 0 && alts[n - 1] == ".") {
                    return fail("GT allele index is not in ALT");
                }
                Ok(Some(n))
            }
        })
        .collect::<Result<_>>()?;
    let dp = number(value("DP"), "DP")?;
    let min_dp = number(value("MIN_DP"), "MIN_DP")?;
    for key in ["DP", "MIN_DP"] {
        if let Some(v) = value(key)
            && v != "."
            && (v.is_empty() || !v.bytes().all(|b| b.is_ascii_digit()))
        {
            return fail(format!("FORMAT/{key} must be a nonnegative integer or ."));
        }
    }
    let gq = number(value("GQ"), "GQ")?;
    if block && min_dp.is_none() {
        *adaptations
            .entry("reference_blocks_missing_MIN_DP".into())
            .or_default() += 1;
    }
    let depth = if block {
        if min_dp.is_none() && args.allow_block_dp && dp.is_some() {
            *adaptations
                .entry("reference_blocks_using_DP".into())
                .or_default() += 1;
            dp
        } else {
            min_dp
        }
    } else {
        dp
    };
    let all_ref = gt.iter().all(|a| *a == Some(0));
    let refcall = h.caller == "deepvariant"
        && f[6] == "RefCall"
        && all_ref
        && (h.refcall || args.caller == "deepvariant");
    if refcall {
        *adaptations
            .entry("DeepVariant_RefCall_accepted".into())
            .or_default() += 1;
    }
    let confidence_block = block
        && all_ref
        && f[6] == "."
        && matches!(h.caller.as_str(), "deepvariant" | "gatk" | "dragen");
    if confidence_block {
        *adaptations
            .entry("reference_blocks_with_unset_site_filter".into())
            .or_default() += 1;
    }
    let filters_pass =
        f[6] == "PASS" || refcall || confidence_block || (f[6] == "." && args.allow_unfiltered);
    let ft_pass = value("FT").is_none_or(|s| matches!(s, "PASS" | "."));
    let known_alleles = f[3]
        .bytes()
        .all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'));
    let called_literal = gt.iter().flatten().all(|&a| {
        a == 0
            || alts[a - 1]
                .bytes()
                .all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
    });
    let complex = !block
        && (end != ref_end
            || alts.iter().any(|a| {
                !matches!(*a, "<NON_REF>" | "<*>")
                    && (a.len() != f[3].len()
                        || !a
                            .bytes()
                            .all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')))
            }));
    let state = if !reference.matches(contig, start, f[3].as_bytes())? {
        State::ReferenceMismatch
    } else if !known_alleles {
        State::ReferenceAmbiguous
    } else if gt.iter().all(Option::is_none) {
        State::NoCall
    } else if gt.iter().any(Option::is_none) {
        State::PartialNoCall
    } else if gt.len() > 2 {
        State::UnsupportedPloidy
    } else if !ft_pass {
        State::GenotypeFiltered
    } else if !filters_pass {
        if f[6] == "." {
            State::FilterNotAssessed
        } else {
            State::SiteFiltered
        }
    } else if depth.is_none() || gq.is_none() {
        State::QualityMissing
    } else if depth.unwrap() < args.min_dp {
        State::LowDepth
    } else if gq.unwrap() < args.min_gq {
        State::LowGq
    } else if !called_literal || (block && !all_ref) {
        State::UnsupportedAllele
    } else if complex {
        State::ComplexVariant
    } else if all_ref {
        State::CallableReference
    } else {
        State::CallableVariant
    };
    Ok(Record {
        chrom: f[0],
        contig,
        start,
        end,
        state,
    })
}
