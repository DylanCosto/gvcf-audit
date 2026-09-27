use crate::{Result, STATES, State, fail, reader, reference::Reference};
use serde::Serialize;
use std::{
    cmp::Reverse,
    collections::BinaryHeap,
    fs::File,
    io::{BufRead, BufWriter, Write},
    path::Path,
};

#[derive(Serialize)]
pub struct Region {
    pub contig: String,
    pub start: u64,
    pub end: u64,
    pub name: String,
    pub bases_by_state: [u64; STATES.len()],
}

pub struct Scope {
    pub regions: Vec<Region>,
    pub by_contig: Vec<Vec<usize>>,
    pub union: Vec<Vec<(u64, u64)>>,
}

impl Scope {
    pub fn load(bed: Option<&Path>, reference: &Reference) -> Result<Self> {
        let n = reference.contigs.len();
        let mut out = Self {
            regions: Vec::new(),
            by_contig: vec![Vec::new(); n],
            union: vec![Vec::new(); n],
        };
        if let Some(path) = bed {
            for (line_no, line) in reader(path)?.lines().enumerate() {
                let line = line?;
                if line.trim().is_empty()
                    || line.starts_with('#')
                    || line.starts_with("track ")
                    || line.starts_with("browser ")
                {
                    continue;
                }
                let f: Vec<_> = line.split('\t').collect();
                if f.len() < 3 {
                    return fail(format!(
                        "BED line {}: expected at least 3 tab-separated columns",
                        line_no + 1
                    ));
                }
                let id = reference.resolve(f[0])?;
                let start: u64 = f[1].parse()?;
                let end: u64 = f[2].parse()?;
                if start >= end || end > reference.contigs[id].length {
                    return fail(format!(
                        "BED line {}: empty or out-of-bounds interval",
                        line_no + 1
                    ));
                }
                let name = f
                    .get(3)
                    .filter(|s| !s.is_empty())
                    .map(|s| s.to_string())
                    .unwrap_or_else(|| format!("{}:{start}-{end}", f[0]));
                out.add(id, reference.contigs[id].name.clone(), start, end, name);
            }
            if out.regions.is_empty() {
                return fail("The BED contains no intervals");
            }
        } else {
            for (id, c) in reference.contigs.iter().enumerate() {
                out.add(id, c.name.clone(), 0, c.length, c.name.clone());
            }
        }
        for id in 0..n {
            out.by_contig[id].sort_by_key(|&i| (out.regions[i].start, out.regions[i].end));
            for &i in &out.by_contig[id] {
                let r = &out.regions[i];
                if let Some(last) = out.union[id].last_mut().filter(|v| r.start <= v.1) {
                    last.1 = last.1.max(r.end);
                } else {
                    out.union[id].push((r.start, r.end));
                }
            }
        }
        Ok(out)
    }

    fn add(&mut self, id: usize, contig: String, start: u64, end: u64, name: String) {
        self.by_contig[id].push(self.regions.len());
        self.regions.push(Region {
            contig,
            start,
            end,
            name,
            bases_by_state: [0; STATES.len()],
        });
    }

    pub fn bases(&self) -> u64 {
        self.union.iter().flatten().map(|(s, e)| e - s).sum()
    }
}

struct BedSink {
    output: BufWriter<File>,
    pending: Option<(String, u64, u64, State)>,
}

impl BedSink {
    fn new(path: &Path) -> Result<Self> {
        Ok(Self {
            output: BufWriter::new(File::create(path)?),
            pending: None,
        })
    }
    fn push(&mut self, c: &str, start: u64, end: u64, state: State) -> Result<()> {
        if let Some((old_c, _, old_end, old_state)) = self.pending.as_mut()
            && old_c == c
            && *old_end == start
            && *old_state == state
        {
            *old_end = end;
            return Ok(());
        }
        self.flush_pending()?;
        self.pending = Some((c.into(), start, end, state));
        Ok(())
    }
    fn flush_pending(&mut self) -> Result<()> {
        if let Some((c, s, e, state)) = self.pending.take() {
            writeln!(self.output, "{c}\t{s}\t{e}\t{}", state.name())?;
        }
        Ok(())
    }
    fn finish(&mut self) -> Result<()> {
        self.flush_pending()?;
        self.output.flush()?;
        Ok(())
    }
}

pub struct Output {
    pub scope: Scope,
    pub counts: [u64; STATES.len()],
    pub contig_counts: Vec<[u64; STATES.len()]>,
    callable: BedSink,
    unresolved: BedSink,
}

impl Output {
    pub fn new(scope: Scope, dir: &Path) -> Result<Self> {
        Ok(Self {
            contig_counts: vec![[0; STATES.len()]; scope.union.len()],
            scope,
            counts: [0; STATES.len()],
            callable: BedSink::new(&dir.join("callable.bed"))?,
            unresolved: BedSink::new(&dir.join("unresolved.bed"))?,
        })
    }
    pub fn finish(&mut self, dir: &Path) -> Result<()> {
        self.callable.finish()?;
        self.unresolved.finish()?;
        if self.counts.iter().sum::<u64>() != self.scope.bases() {
            return fail("Internal interval accounting error");
        }
        let mut f = BufWriter::new(File::create(dir.join("regions.tsv"))?);
        write!(
            f,
            "contig\tstart\tend\tname\tbases\tcallable_bases\tcallable_percent"
        )?;
        for state in STATES {
            write!(f, "\t{}", state.name())?;
        }
        writeln!(f)?;
        for r in &self.scope.regions {
            let size = r.end - r.start;
            if r.bases_by_state.iter().sum::<u64>() != size {
                return fail("Internal region accounting error");
            }
            let called = r.bases_by_state[0] + r.bases_by_state[1];
            write!(
                f,
                "{}\t{}\t{}\t{}\t{size}\t{called}\t{:.4}",
                r.contig,
                r.start,
                r.end,
                r.name,
                100.0 * called as f64 / size as f64
            )?;
            for n in r.bases_by_state {
                write!(f, "\t{n}")?;
            }
            writeln!(f)?;
        }
        f.flush()?;
        Ok(())
    }
}

pub struct Emitter {
    contig: usize,
    masks: Vec<(u64, u64)>,
    mask: usize,
    scope: Vec<(u64, u64)>,
    interval: usize,
    regions: Vec<usize>,
    next_region: usize,
    active_regions: Vec<usize>,
}

impl Emitter {
    pub fn new(contig: usize, output: &Output, reference: &Reference) -> Result<Self> {
        let scope = output.scope.union[contig].clone();
        let masks = if scope.is_empty() {
            Vec::new()
        } else {
            reference.ambiguous_runs(contig)?
        };
        Ok(Self {
            contig,
            masks,
            mask: 0,
            scope,
            interval: 0,
            regions: output.scope.by_contig[contig].clone(),
            next_region: 0,
            active_regions: Vec::new(),
        })
    }

    fn emit(&mut self, start: u64, end: u64, state: State, out: &mut Output) -> Result<()> {
        while self.interval < self.scope.len() && self.scope[self.interval].1 <= start {
            self.interval += 1;
        }
        let mut i = self.interval;
        while i < self.scope.len() && self.scope[i].0 < end {
            let mut p = start.max(self.scope[i].0);
            let stop = end.min(self.scope[i].1);
            while p < stop {
                while self.mask < self.masks.len() && self.masks[self.mask].1 <= p {
                    self.mask += 1;
                }
                let (next, label) = match self.masks.get(self.mask) {
                    Some(&(s, e)) if s <= p => (e.min(stop), State::ReferenceAmbiguous),
                    Some(&(s, _)) => (s.min(stop), state),
                    None => (stop, state),
                };
                self.piece(p, next, label, out)?;
                p = next;
            }
            if self.scope[i].1 >= end {
                break;
            }
            i += 1;
        }
        Ok(())
    }

    fn piece(&mut self, start: u64, end: u64, state: State, out: &mut Output) -> Result<()> {
        if start == end {
            return Ok(());
        }
        out.counts[state as usize] += end - start;
        out.contig_counts[self.contig][state as usize] += end - start;
        while self.next_region < self.regions.len()
            && out.scope.regions[self.regions[self.next_region]].start < end
        {
            self.active_regions.push(self.regions[self.next_region]);
            self.next_region += 1;
        }
        self.active_regions
            .retain(|&i| out.scope.regions[i].end > start);
        for &i in &self.active_regions {
            let r = &mut out.scope.regions[i];
            r.bases_by_state[state as usize] += r.end.min(end) - r.start.max(start);
        }
        let name = &out.scope.regions[out.scope.by_contig[self.contig][0]].contig;
        if state.callable() {
            out.callable.push(name, start, end, state)?;
        } else {
            out.unresolved.push(name, start, end, state)?;
        }
        Ok(())
    }
}

pub struct Sweep {
    pub contig: usize,
    pub raw_name: String,
    pub last_start: u64,
    cursor: u64,
    ends: BinaryHeap<Reverse<(u64, usize)>>,
    counts: [u64; STATES.len()],
    active: u64,
    emitter: Emitter,
}

impl Sweep {
    pub fn new(
        contig: usize,
        raw_name: String,
        out: &Output,
        reference: &Reference,
    ) -> Result<Self> {
        Ok(Self {
            contig,
            raw_name,
            last_start: 0,
            cursor: 0,
            ends: BinaryHeap::new(),
            counts: [0; STATES.len()],
            active: 0,
            emitter: Emitter::new(contig, out, reference)?,
        })
    }

    pub fn advance(&mut self, stop: u64, out: &mut Output) -> Result<()> {
        while self.cursor < stop {
            let next = self.ends.peek().map(|e| e.0.0.min(stop)).unwrap_or(stop);
            let state = match self.active {
                0 => State::NoRecord,
                1 => STATES[self.counts.iter().position(|&n| n != 0).unwrap()],
                _ => State::OverlappingRecords,
            };
            self.emitter.emit(self.cursor, next, state, out)?;
            self.cursor = next;
            while self.ends.peek().is_some_and(|e| e.0.0 == self.cursor) {
                let Reverse((_, state)) = self.ends.pop().unwrap();
                self.counts[state] -= 1;
                self.active -= 1;
            }
        }
        Ok(())
    }

    pub fn add(&mut self, start: u64, end: u64, state: State, out: &mut Output) -> Result<()> {
        if start < self.last_start {
            return fail("Records are not sorted by position");
        }
        self.advance(start, out)?;
        self.last_start = start;
        self.ends.push(Reverse((end, state as usize)));
        self.counts[state as usize] += 1;
        self.active += 1;
        Ok(())
    }
}
