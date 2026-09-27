use crate::{Result, STATES, fail};
use std::{
    collections::BTreeMap,
    fs::File,
    io::{BufRead, BufReader, Seek, SeekFrom},
    path::{Path, PathBuf},
};

#[derive(Clone)]
pub struct Segment {
    pub start: u64,
    pub end: u64,
    pub state: usize,
}

fn parse(line: &str, callable: bool) -> Result<(String, Segment)> {
    let f: Vec<_> = line.trim_end_matches(['\n', '\r']).split('\t').collect();
    if f.len() != 4 || f[0].is_empty() {
        return fail("Expected four columns in an audit BED");
    }
    let start = f[1].parse()?;
    let end = f[2].parse()?;
    let state = STATES
        .iter()
        .position(|s| s.name() == f[3])
        .ok_or("Unknown BED state")?;
    if start >= end || STATES[state].callable() != callable {
        return fail("Invalid interval or state in an audit BED");
    }
    Ok((f[0].into(), Segment { start, end, state }))
}

// Index contig byte ranges so audits can use different contig orders without loading every interval.
pub struct BedIndex {
    path: PathBuf,
    callable: bool,
    ranges: BTreeMap<String, (u64, u64)>,
}

impl BedIndex {
    pub fn open(
        path: &Path,
        callable: bool,
        scope: &BTreeMap<String, Vec<(u64, u64)>>,
    ) -> Result<Self> {
        let mut reader = BufReader::new(File::open(path)?);
        let mut ranges = BTreeMap::new();
        let mut line = String::new();
        let mut offset = 0;
        let mut previous = String::new();
        let mut previous_end = 0;
        loop {
            line.clear();
            let n = reader.read_line(&mut line)? as u64;
            if n == 0 {
                break;
            }
            let (contig, segment) = parse(&line, callable)
                .map_err(|e| format!("{} at byte {offset}: {e}", path.display()))?;
            if !scope.contains_key(&contig) {
                return fail("BED contig is outside the reported scope");
            }
            if contig != previous {
                if ranges.contains_key(&contig) {
                    return fail("BED contig appears in multiple blocks");
                }
                ranges.insert(contig.clone(), (offset, offset + n));
                previous = contig.clone();
                previous_end = 0;
            }
            if segment.start < previous_end {
                return fail("Audit BED intervals overlap or are out of order");
            }
            previous_end = segment.end;
            ranges.get_mut(&contig).unwrap().1 = offset + n;
            offset += n;
        }
        Ok(Self {
            path: path.into(),
            callable,
            ranges,
        })
    }

    fn reader(&self, contig: &str) -> Result<ContigReader> {
        let (start, end) = self.ranges.get(contig).copied().unwrap_or((0, 0));
        let mut file = File::open(&self.path)?;
        file.seek(SeekFrom::Start(start))?;
        Ok(ContigReader {
            reader: BufReader::new(file),
            remaining: end - start,
            contig: contig.into(),
            callable: self.callable,
            line: String::new(),
        })
    }
}

struct ContigReader {
    reader: BufReader<File>,
    remaining: u64,
    contig: String,
    callable: bool,
    line: String,
}
impl ContigReader {
    fn next(&mut self) -> Result<Option<Segment>> {
        if self.remaining == 0 {
            return Ok(None);
        }
        self.line.clear();
        let n = self.reader.read_line(&mut self.line)? as u64;
        if n == 0 || n > self.remaining {
            return fail("Audit BED changed while being compared");
        }
        self.remaining -= n;
        let (contig, segment) = parse(&self.line, self.callable)?;
        if contig != self.contig {
            return fail("Audit BED changed while being compared");
        }
        Ok(Some(segment))
    }
}

pub struct Partition<'a> {
    readers: [ContigReader; 2],
    pending: [Option<Segment>; 2],
    scope: &'a [(u64, u64)],
    interval: usize,
    cursor: u64,
}
impl<'a> Partition<'a> {
    pub fn new(indices: &[BedIndex; 2], contig: &str, scope: &'a [(u64, u64)]) -> Result<Self> {
        let mut readers = [indices[0].reader(contig)?, indices[1].reader(contig)?];
        let pending = [readers[0].next()?, readers[1].next()?];
        Ok(Self {
            readers,
            pending,
            scope,
            interval: 0,
            cursor: scope[0].0,
        })
    }
    pub fn next(&mut self) -> Result<Option<Segment>> {
        let which = match (&self.pending[0], &self.pending[1]) {
            (Some(a), Some(b)) => usize::from(a.start > b.start),
            (Some(_), None) => 0,
            (None, Some(_)) => 1,
            (None, None) => {
                if self.interval != self.scope.len() {
                    return fail("Audit BEDs do not cover the reported scope");
                }
                return Ok(None);
            }
        };
        let segment = self.pending[which].take().unwrap();
        if self.interval == self.scope.len()
            || segment.start != self.cursor
            || segment.end > self.scope[self.interval].1
        {
            return fail(
                "Audit BEDs contain gaps, overlaps or intervals outside the reported scope",
            );
        }
        self.cursor = segment.end;
        if self.cursor == self.scope[self.interval].1 {
            self.interval += 1;
            if let Some(&(start, _)) = self.scope.get(self.interval) {
                self.cursor = start;
            }
        }
        self.pending[which] = self.readers[which].next()?;
        Ok(Some(segment))
    }
}
