use crate::{Result, fail};
use memmap2::{Mmap, MmapOptions};
use std::{
    cell::{Ref, RefCell},
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

pub struct Contig {
    pub name: String,
    pub length: u64,
    offset: u64,
    bases: u64,
    width: u64,
}

pub struct Reference {
    pub contigs: Vec<Contig>,
    names: HashMap<String, usize>,
    file: File,
    cache: RefCell<Option<(usize, Mmap)>>,
}

impl Reference {
    pub fn open(path: &Path) -> Result<Self> {
        let index = std::fs::read_to_string(format!("{}.fai", path.display()))
            .map_err(|e| format!("Cannot read FASTA index: {e}. Create it with samtools faidx."))?;
        let file = File::open(path)?;
        let file_size = file.metadata()?.len();
        let mut names = HashMap::new();
        let mut contigs = Vec::new();
        for line in index.lines() {
            let f: Vec<_> = line.split('\t').collect();
            if f.len() < 5 || f[0].is_empty() {
                return fail("Malformed FASTA index");
            }
            let c = Contig {
                name: f[0].into(),
                length: f[1].parse()?,
                offset: f[2].parse()?,
                bases: f[3].parse()?,
                width: f[4].parse()?,
            };
            if c.length == 0 || c.bases == 0 || c.width < c.bases {
                return fail("Invalid FASTA index dimensions");
            }
            let last = c
                .offset
                .checked_add(
                    ((c.length - 1) / c.bases)
                        .checked_mul(c.width)
                        .ok_or("FASTA index overflow")?,
                )
                .and_then(|p| p.checked_add((c.length - 1) % c.bases))
                .ok_or("FASTA index overflow")?;
            if last >= file_size || names.insert(c.name.clone(), contigs.len()).is_some() {
                return fail("FASTA index is out of bounds or contains duplicate contigs");
            }
            contigs.push(c);
        }
        if contigs.is_empty() {
            return fail("Empty FASTA index");
        }
        Self::validate_layout(&file, &contigs)?;
        Ok(Self {
            contigs,
            names,
            file,
            cache: RefCell::new(None),
        })
    }

    fn validate_layout(file: &File, contigs: &[Contig]) -> Result<()> {
        let mut current: Option<usize> = None;
        let mut bases = 0u64;
        let mut start = 0u64;
        let mut input = BufReader::with_capacity(1024 * 1024, file);
        if input.fill_buf()?.starts_with(&[31, 139]) {
            return fail("Use an uncompressed FASTA with its .fai index");
        }
        let mut buffer = Vec::new();
        loop {
            buffer.clear();
            let width = input.read_until(b'\n', &mut buffer)? as u64;
            if width == 0 {
                break;
            }
            let line = buffer.strip_suffix(b"\n").unwrap_or(&buffer);
            let line = line.strip_suffix(b"\r").unwrap_or(line);
            if line.is_empty() && current.is_none_or(|id| bases == contigs[id].length) {
                start += width;
                continue;
            }
            if let Some(header) = line.strip_prefix(b">") {
                if let Some(id) = current
                    && bases != contigs[id].length
                {
                    return fail("FASTA sequence length differs from its index; rebuild the .fai");
                }
                let id = current.map_or(0, |id| id + 1);
                let c = contigs
                    .get(id)
                    .ok_or("FASTA has contigs missing from its index")?;
                let name = header
                    .split(|b| b.is_ascii_whitespace())
                    .next()
                    .unwrap_or_default();
                if name != c.name.as_bytes() || c.offset != start + width {
                    return fail(
                        "FASTA headers or offsets differ from its index; rebuild the .fai",
                    );
                }
                current = Some(id);
                bases = 0;
            } else {
                let id = current.ok_or("FASTA must begin with a > header")?;
                let c = &contigs[id];
                let take = c.bases.min(c.length - bases);
                let expected = c.offset + bases / c.bases * c.width;
                if take == 0
                    || line.len() as u64 != take
                    || start != expected
                    || (bases + take < c.length && width != c.width)
                {
                    return fail("FASTA sequence layout differs from its index; rebuild the .fai");
                }
                bases += take;
            }
            start += width;
        }
        if current != Some(contigs.len() - 1) || bases != contigs.last().unwrap().length {
            return fail(
                "FASTA contigs or sequence lengths differ from its index; rebuild the .fai",
            );
        }
        Ok(())
    }

    pub fn resolve(&self, name: &str) -> Result<usize> {
        self.resolve_optional(name)?
            .ok_or_else(|| format!("Contig {name:?} is not in the reference").into())
    }

    pub fn resolve_optional(&self, name: &str) -> Result<Option<usize>> {
        if let Some(&id) = self.names.get(name) {
            return Ok(Some(id));
        }
        let stripped = name.strip_prefix("chr").unwrap_or(name);
        let aliases: Vec<_> = self
            .contigs
            .iter()
            .enumerate()
            .filter(|(_, c)| {
                let other = c.name.strip_prefix("chr").unwrap_or(&c.name);
                other == stripped || (matches!(other, "M" | "MT") && matches!(stripped, "M" | "MT"))
            })
            .map(|(i, _)| i)
            .collect();
        match aliases.as_slice() {
            [id] => Ok(Some(*id)),
            [] => Ok(None),
            _ => fail(format!("Contig {name:?} has ambiguous reference aliases")),
        }
    }

    fn mapped(&self, id: usize) -> Result<Ref<'_, Mmap>> {
        if self
            .cache
            .borrow()
            .as_ref()
            .is_none_or(|(cached, _)| *cached != id)
        {
            let c = &self.contigs[id];
            let span = (c.length - 1) / c.bases * c.width + (c.length - 1) % c.bases + 1;
            let length = usize::try_from(span)
                .map_err(|_| "Reference contig is too large for this platform")?;
            // The reference must remain unchanged while its read-only mapping is in use.
            let map = unsafe {
                MmapOptions::new()
                    .offset(c.offset)
                    .len(length)
                    .map(&self.file)?
            };
            *self.cache.borrow_mut() = Some((id, map));
        }
        Ok(Ref::map(self.cache.borrow(), |cache| {
            &cache.as_ref().unwrap().1
        }))
    }

    pub fn matches(&self, id: usize, start: u64, text: &[u8]) -> Result<bool> {
        let c = &self.contigs[id];
        if start
            .checked_add(text.len() as u64)
            .is_none_or(|e| e > c.length)
        {
            return fail("REF allele extends past the reference contig");
        }
        let data = self.mapped(id)?;
        for (i, &b) in text.iter().enumerate() {
            let p = start + i as u64;
            let offset = p / c.bases * c.width + p % c.bases;
            if !data[offset as usize].eq_ignore_ascii_case(&b) {
                return Ok(false);
            }
        }
        Ok(true)
    }

    pub fn first_mismatch(
        &self,
        id: usize,
        start: u64,
        text: &str,
    ) -> Result<Option<(u64, char, char)>> {
        let c = &self.contigs[id];
        let data = self.mapped(id)?;
        for (i, b) in text.bytes().enumerate() {
            let p = start + i as u64;
            let actual = data[(p / c.bases * c.width + p % c.bases) as usize];
            if !actual.eq_ignore_ascii_case(&b) {
                return Ok(Some((p, b as char, actual as char)));
            }
        }
        Ok(None)
    }

    pub fn ambiguous_runs(&self, id: usize) -> Result<Vec<(u64, u64)>> {
        let c = &self.contigs[id];
        let data = self.mapped(id)?;
        let mut runs = Vec::new();
        let mut open = None;
        let mut p = 0;
        while p < c.length {
            let take = (c.length - p).min(c.bases);
            let offset = (p / c.bases * c.width) as usize;
            for (i, &b) in data[offset..offset + take as usize].iter().enumerate() {
                if !b.is_ascii_alphabetic() {
                    return fail("FASTA index does not match sequence layout");
                }
                let pos = p + i as u64;
                if matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T') {
                    if let Some(start) = open.take() {
                        runs.push((start, pos));
                    }
                } else if open.is_none() {
                    open = Some(pos);
                }
            }
            p += take;
        }
        if let Some(start) = open {
            runs.push((start, c.length));
        }
        Ok(runs)
    }
}
