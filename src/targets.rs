use crate::{Result, fail, reader};
use clap::{Parser, ValueEnum};
use std::{
    collections::{BTreeSet, HashMap, HashSet},
    fs::{self, OpenOptions},
    io::{BufRead, BufWriter, Write},
    path::{Path, PathBuf},
    time::{SystemTime, UNIX_EPOCH},
};

#[derive(Parser)]
#[command(
    about = "Convert exon annotations to gene targets",
    after_help = "Reads GTF or GFF3, plain text or gzip. Only exon features become targets.\nCoordinates change from 1-based inclusive to 0-based half-open on both strands.\nUse annotations from the same assembly as your audit reference; no liftover or contig renaming is performed."
)]
pub struct Args {
    /// GTF or GFF3 annotation file
    #[arg(long)]
    annotation: PathBuf,
    /// Detect from the filename or GFF3 header; use explicitly for other filenames
    #[arg(long, value_enum, default_value = "auto")]
    format: Format,
    /// New targets TSV file; existing paths are never overwritten
    #[arg(long)]
    out: PathBuf,
    /// Exact annotation contig name to include (repeat for multiple contigs)
    #[arg(long)]
    contig: Vec<String>,
    /// Exact gene_id (GTF) or gene ID (GFF3) to include; repeat as needed
    #[arg(long)]
    gene: Vec<String>,
    /// Suppress the completion summary
    #[arg(long)]
    quiet: bool,
}

#[derive(Clone, Copy, PartialEq, Eq, ValueEnum)]
enum Format {
    Auto,
    Gtf,
    Gff3,
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord)]
struct Target {
    contig: String,
    start: u64,
    end: u64,
    gene: String,
    exon: String,
}

#[derive(PartialEq, Eq)]
struct Node {
    contig: String,
    kind: String,
    strand: String,
    parents: Vec<String>,
}

struct Exon {
    contig: String,
    start: u64,
    end: u64,
    strand: String,
    id: String,
    parents: Vec<String>,
    line: usize,
}

const GENERATED: &str = "gvcf-audit:exon:";

fn identifier(s: &str) -> Result<String> {
    if s.trim().is_empty() || s == "." || s.chars().any(char::is_control) {
        return fail("Empty, missing or control-containing identifier");
    }
    Ok(s.to_string())
}

fn decode(s: &str) -> Result<String> {
    let mut bytes = Vec::with_capacity(s.len());
    let mut i = 0;
    while i < s.len() {
        if s.as_bytes()[i] == b'%' {
            let hex = s.get(i + 1..i + 3).ok_or("Invalid percent escape")?;
            bytes.push(u8::from_str_radix(hex, 16).map_err(|_| "Invalid percent escape")?);
            i += 3;
        } else {
            bytes.push(s.as_bytes()[i]);
            i += 1;
        }
    }
    Ok(String::from_utf8(bytes)?)
}

fn keep_attribute(attrs: &mut HashMap<String, String>, key: &str, value: String) -> Result<()> {
    if matches!(key, "gene_id" | "exon_id" | "ID" | "Parent")
        && attrs.insert(key.to_string(), value).is_some()
    {
        return fail(format!("Duplicate {key} attribute"));
    }
    Ok(())
}

fn gtf_without_comment(text: &str) -> &str {
    if !text.contains('#') {
        return text;
    }
    let mut quoted = false;
    let mut escaped = false;
    for (i, c) in text.char_indices() {
        if escaped {
            escaped = false;
        } else if quoted && c == '\\' {
            escaped = true;
        } else if c == '"' {
            quoted = !quoted;
        } else if c == '#' && !quoted {
            return text[..i].trim_end();
        }
    }
    text
}

fn gtf_attributes(mut text: &str) -> Result<HashMap<String, String>> {
    let mut attrs = HashMap::new();
    while !text.trim().is_empty() {
        text = text.trim_start();
        let split = text
            .find(char::is_whitespace)
            .ok_or("Expected GTF key and value")?;
        let key = &text[..split];
        text = text[split..].trim_start();
        let value;
        if let Some(rest) = text.strip_prefix('"') {
            let mut parsed = String::new();
            let mut chars = rest.char_indices();
            let mut end = None;
            while let Some((i, c)) = chars.next() {
                match c {
                    '"' => {
                        end = Some(i + 1);
                        break;
                    }
                    '\\' => {
                        let (_, escaped) = chars.next().ok_or("Unfinished GTF escape")?;
                        if !matches!(escaped, '"' | '\\') {
                            return fail("Unsupported GTF escape");
                        }
                        parsed.push(escaped);
                    }
                    _ => parsed.push(c),
                }
            }
            text = &rest[end.ok_or("Unclosed GTF quoted value")?..];
            value = parsed;
        } else {
            let end = text.find(';').unwrap_or(text.len());
            value = text[..end].trim().to_string();
            if value.split_whitespace().count() != 1 {
                return fail("Expected one GTF value");
            }
            text = &text[end..];
        }
        text = text.trim_start();
        if !text.is_empty() {
            text = text
                .strip_prefix(';')
                .ok_or("Expected semicolon after GTF value")?;
        }
        keep_attribute(&mut attrs, key, value)?;
    }
    Ok(attrs)
}

fn gff_attributes(text: &str) -> Result<HashMap<String, String>> {
    let mut attrs = HashMap::new();
    if text == "." {
        return Ok(attrs);
    }
    for field in text.split(';').filter(|s| !s.is_empty()) {
        let (key, value) = field
            .split_once('=')
            .ok_or("Expected GFF3 key=value attribute")?;
        // Split Parent lists before decoding, since an encoded comma is part of an ID.
        keep_attribute(&mut attrs, &decode(key)?, value.to_string())?;
    }
    Ok(attrs)
}

fn exon_id(value: Option<String>, start: u64, end: u64, strand: &str) -> Result<String> {
    if let Some(value) = value {
        let id = identifier(&value)?;
        if id.starts_with(GENERATED) {
            return fail(format!(
                "Exon IDs starting with {GENERATED} are reserved for generated IDs"
            ));
        }
        Ok(id)
    } else {
        Ok(format!("{GENERATED}{start}-{end}:{strand}"))
    }
}

fn is_gene(kind: &str) -> bool {
    matches!(
        kind,
        "gene"
            | "pseudogene"
            | "ncRNA_gene"
            | "transposable_element_gene"
            | "SO:0000704"
            | "SO:0000336"
            | "SO:0001263"
    )
}

fn directive<'a>(text: &'a str, name: &str) -> Option<&'a str> {
    text.strip_prefix(name)
        .filter(|rest| rest.is_empty() || rest.starts_with(char::is_whitespace))
        .map(str::trim)
}

fn ancestors(
    id: &str,
    nodes: &HashMap<String, Node>,
    cache: &mut HashMap<String, BTreeSet<String>>,
    visiting: &mut HashSet<String>,
) -> Result<BTreeSet<String>> {
    if let Some(genes) = cache.get(id) {
        return Ok(genes.clone());
    }
    if visiting.len() >= 256 {
        return fail("GFF3 parent chain exceeds 256 levels");
    }
    if !visiting.insert(id.to_string()) {
        return fail(format!("GFF3 parent cycle at {id}"));
    }
    let node = nodes
        .get(id)
        .ok_or_else(|| format!("Missing GFF3 parent {id}"))?;
    let mut genes = BTreeSet::new();
    if is_gene(&node.kind) {
        genes.insert(id.to_string());
    } else {
        if node.parents.is_empty() {
            return fail(format!("GFF3 feature {id} has no gene ancestor"));
        }
        for parent in &node.parents {
            let other = nodes
                .get(parent)
                .ok_or_else(|| format!("Missing GFF3 parent {parent}"))?;
            same_locus(&node.contig, &node.strand, other)?;
            genes.extend(ancestors(parent, nodes, cache, visiting)?);
        }
    }
    visiting.remove(id);
    cache.insert(id.to_string(), genes.clone());
    Ok(genes)
}

fn same_locus(contig: &str, strand: &str, parent: &Node) -> Result<()> {
    if contig != parent.contig
        || (matches!(strand, "+" | "-")
            && matches!(parent.strand.as_str(), "+" | "-")
            && strand != parent.strand)
    {
        return fail("GFF3 child and parent have incompatible contigs or strands");
    }
    Ok(())
}

fn write_targets(out: &Path, rows: &BTreeSet<Target>) -> Result<()> {
    let parent = out
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    fs::create_dir_all(parent)?;
    let temp = parent.join(format!(
        ".gvcf-targets-{}-{}",
        std::process::id(),
        SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos()
    ));
    let file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&temp)?;
    let result = (|| -> Result<()> {
        let mut f = BufWriter::new(file);
        writeln!(f, "contig\tstart\tend\tgene\texon")?;
        for r in rows {
            writeln!(
                f,
                "{}\t{}\t{}\t{}\t{}",
                r.contig, r.start, r.end, r.gene, r.exon
            )?;
        }
        f.flush()?;
        drop(f);
        // Publishing a hard link fails if any destination already exists, including a symlink.
        fs::hard_link(&temp, out)?;
        Ok(())
    })();
    let _ = fs::remove_file(temp);
    result
}

pub fn run(args: &Args) -> Result<bool> {
    if fs::symlink_metadata(&args.out).is_ok() {
        return fail("Output already exists; choose a new --out file");
    }
    let mut format = args.format;
    if format == Format::Auto {
        let name = args.annotation.to_string_lossy().to_ascii_lowercase();
        let name = name
            .strip_suffix(".gz")
            .or_else(|| name.strip_suffix(".bgz"))
            .unwrap_or(&name);
        format = if name.ends_with(".gtf") {
            Format::Gtf
        } else if name.ends_with(".gff3") || name.ends_with(".gff") {
            Format::Gff3
        } else {
            Format::Auto
        };
    }
    let mut rows = BTreeSet::new();
    let mut nodes = HashMap::new();
    let mut exons = Vec::new();
    let mut regions: HashMap<String, (u64, u64)> = HashMap::new();
    let mut exon_records = 0u64;
    let mut fallback_records = 0u64;
    for (line_no, line) in reader(&args.annotation)?.lines().enumerate() {
        let line = line.map_err(|e| format!("Annotation line {}: {e}", line_no + 1))?;
        let text = line.trim_end_matches('\r');
        let result = (|| -> Result<bool> {
            if text.trim_end() == "##FASTA" {
                if format != Format::Gff3 {
                    return fail("FASTA directive requires GFF3");
                }
                return Ok(false);
            }
            if let Some(version) = directive(text, "##gff-version") {
                if !matches!(version.split('.').next(), Some("3"))
                    || version
                        .split('.')
                        .any(|part| part.is_empty() || !part.bytes().all(|b| b.is_ascii_digit()))
                    || format == Format::Gtf
                {
                    return fail("Expected GFF3 annotation; GFF1/GFF2 are not supported");
                }
                format = Format::Gff3;
                return Ok(true);
            }
            if text.trim().is_empty() || text.starts_with('#') {
                if let Some(region) = directive(text, "##sequence-region") {
                    let f: Vec<_> = region.split_whitespace().collect();
                    if f.len() != 3 {
                        return fail("Invalid sequence-region directive");
                    }
                    let contig = identifier(&decode(f[0])?)?;
                    let start: u64 = f[1].parse()?;
                    let end: u64 = f[2].parse()?;
                    if start == 0 || end < start {
                        return fail("Invalid sequence-region bounds");
                    }
                    if regions.insert(contig, (start - 1, end)).is_some() {
                        return fail("Repeated sequence-region directive");
                    }
                }
                return Ok(true);
            }
            if format == Format::Auto {
                return fail("Cannot detect annotation format; use --format gtf or --format gff3");
            }
            let text = if format == Format::Gtf {
                gtf_without_comment(text)
            } else {
                text
            };
            let f: Vec<_> = text.split('\t').collect();
            if f.len() != 9 {
                return fail("Expected nine tab-separated annotation columns");
            }
            let contig = identifier(&if format == Format::Gff3 {
                decode(f[0])?
            } else {
                f[0].to_string()
            })?;
            if contig.starts_with('#') || contig.chars().any(char::is_whitespace) {
                return fail(
                    "Contig name cannot contain whitespace or start with # in targets TSV",
                );
            }
            let start: u64 = f[3].parse()?;
            let end: u64 = f[4].parse()?;
            if start == 0 || end < start {
                return fail("Expected positive, increasing 1-based coordinates");
            }
            if !matches!(f[6], "+" | "-" | "." | "?") {
                return fail("Invalid strand");
            }
            if !args.contig.is_empty() && !args.contig.contains(&contig) {
                return Ok(true);
            }
            let start = start - 1;
            let is_exon = matches!(f[2], "exon" | "SO:0000147");
            if format == Format::Gtf {
                if !is_exon {
                    return Ok(true);
                }
                let mut attrs = gtf_attributes(f[8])?;
                let gene = identifier(attrs.get("gene_id").ok_or("GTF exon is missing gene_id")?)?;
                exon_records += 1;
                fallback_records += u64::from(!attrs.contains_key("exon_id"));
                let exon = exon_id(attrs.remove("exon_id"), start, end, f[6])?;
                if !args.gene.is_empty() && !args.gene.contains(&gene) {
                    return Ok(true);
                }
                rows.insert(Target {
                    contig,
                    start,
                    end,
                    gene,
                    exon,
                });
            } else {
                let attrs = gff_attributes(f[8])?;
                let id = attrs
                    .get("ID")
                    .map(|id| {
                        if id.contains(',') {
                            return fail("GFF3 ID must be a single value");
                        }
                        identifier(&decode(id)?)
                    })
                    .transpose()?;
                let parents = attrs
                    .get("Parent")
                    .map(|s| {
                        s.split(',')
                            .map(|p| identifier(&decode(p)?))
                            .collect::<Result<BTreeSet<_>>>()
                    })
                    .transpose()?
                    .unwrap_or_default()
                    .into_iter()
                    .collect::<Vec<_>>();
                if let Some(id) = &id {
                    let node = Node {
                        contig: contig.clone(),
                        kind: f[2].to_string(),
                        strand: f[6].to_string(),
                        parents: parents.clone(),
                    };
                    if let Some(previous) = nodes.get(id) {
                        if previous != &node {
                            return fail(format!("Conflicting GFF3 definitions for ID {id}"));
                        }
                    } else {
                        nodes.insert(id.clone(), node);
                    }
                }
                if is_exon {
                    if parents.is_empty() {
                        return fail("GFF3 exon is missing Parent");
                    }
                    exon_records += 1;
                    fallback_records += u64::from(id.is_none());
                    exons.push(Exon {
                        contig,
                        start,
                        end,
                        strand: f[6].to_string(),
                        id: exon_id(id, start, end, f[6])?,
                        parents,
                        line: line_no + 1,
                    });
                }
            }
            Ok(true)
        })();
        match result {
            Ok(false) => break,
            Ok(true) => {}
            Err(e) => return fail(format!("Annotation line {}: {e}", line_no + 1)),
        }
    }
    let mut cache = HashMap::new();
    for exon in exons {
        let result = (|| -> Result<()> {
            let mut genes = BTreeSet::new();
            for parent in &exon.parents {
                let node = nodes
                    .get(parent)
                    .ok_or_else(|| format!("Missing GFF3 parent {parent}"))?;
                same_locus(&exon.contig, &exon.strand, node)?;
                genes.extend(ancestors(parent, &nodes, &mut cache, &mut HashSet::new())?);
            }
            for gene in genes {
                same_locus(&exon.contig, &exon.strand, &nodes[&gene])?;
                if args.gene.is_empty() || args.gene.contains(&gene) {
                    rows.insert(Target {
                        contig: exon.contig.clone(),
                        start: exon.start,
                        end: exon.end,
                        gene,
                        exon: exon.id.clone(),
                    });
                }
            }
            Ok(())
        })();
        result.map_err(|e| format!("Annotation line {}: {e}", exon.line))?;
    }
    for row in &rows {
        if regions
            .get(&row.contig)
            .is_some_and(|&(s, e)| row.start < s || row.end > e)
        {
            return fail(format!(
                "Exon {} exceeds sequence-region bounds; circular wrapping is not supported",
                row.exon
            ));
        }
    }
    for contig in &args.contig {
        if !rows.iter().any(|r| &r.contig == contig) {
            return fail(format!("No exon targets match requested contig {contig}"));
        }
    }
    for gene in &args.gene {
        if !rows.iter().any(|r| &r.gene == gene) {
            return fail(format!("No exon targets match requested gene ID {gene}"));
        }
    }
    if rows.is_empty() {
        return fail("No exon targets found; CDS/UTR features are not used to infer exons");
    }
    write_targets(&args.out, &rows)?;
    if !args.quiet {
        let genes: BTreeSet<_> = rows.iter().map(|r| (&r.contig, &r.gene)).collect();
        eprintln!(
            "Read {exon_records} exon records ({fallback_records} without exon IDs); wrote {} unique targets for {} genes to {}",
            rows.len(),
            genes.len(),
            args.out.display()
        );
    }
    Ok(true)
}
