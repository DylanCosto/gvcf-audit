import gzip
import json
import random
import subprocess

from support import random_seed, workspace

repo, study, exe, temporary = workspace()
cases = study / "cases"
cases.mkdir()
checks = []


def feature(kind, start=1, end=10, attrs=".", contig="demo", strand="+"):
    return "\t".join(map(str, [contig, ".", kind, start, end, ".", strand, ".", attrs]))


def run(name, text, expected=None, fmt="gtf", extra=(), error=None):
    src = cases / (name + "." + fmt)
    out = cases / (name + ".tsv")
    src.write_text(text)
    p = subprocess.run(
        [str(exe), "targets", "--annotation", str(src), "--out", str(out), *extra],
        text=True,
        capture_output=True,
    )
    assert p.returncode == (1 if error else 0), (name, p.returncode, p.stderr)
    if error:
        assert error in p.stderr and not out.exists(), (name, p.stderr)
    else:
        lines = out.read_text().splitlines()
        assert lines[0] == "contig\tstart\tend\tgene\texon"
        rows = [
            tuple([f[0], int(f[1]), int(f[2]), f[3], f[4]])
            for l in lines[1:]
            if (f := l.split("\t"))
        ]
        if expected is not None:
            assert rows == sorted(set(expected)), (name, rows, expected)
    assert not list(cases.glob(".gvcf-targets-*"))
    checks.append(name)
    return src, out, p


base = feature("exon", 1, 10, 'gene_id "G"; transcript_id "T"; exon_id "E";')
run("gtf-basic", base, [("demo", 0, 10, "G", "E")])
run(
    "negative-single-base",
    feature("exon", 7, 7, 'gene_id "G"; exon_id "E";', strand="-"),
    [("demo", 6, 7, "G", "E")],
)
run("duplicate-rows", base + "\n" + base, [("demo", 0, 10, "G", "E")])
run(
    "gtf-attributes",
    feature(
        "exon",
        1,
        10,
        'tag "one"; gene_name "a;b"; gene_id "G;1"; tag "two"; exon_id "E=1"; exon_number 1;',
    ),
    [("demo", 0, 10, "G;1", "E=1")],
)
run(
    "gtf-escaped",
    feature("exon", 1, 10, 'gene_id "G\\"1"; exon_id "E\\\\1";'),
    [("demo", 0, 10, 'G"1', "E\\1")],
)
run(
    "gtf-fallback",
    feature("exon", 1, 10, 'gene_id "G"; exon_number "1";', strand="-"),
    [("demo", 0, 10, "G", "gvcf-audit:exon:0-10:-")],
)
run(
    "gtf-unquoted",
    feature("exon", 1, 10, "gene_id G; exon_id E"),
    [("demo", 0, 10, "G", "E")],
)
run("inline-comment", base + " # explanatory comment", [("demo", 0, 10, "G", "E")])
run(
    "inline-comment-tab", base + "\t# explanatory\tcomment", [("demo", 0, 10, "G", "E")]
)
run(
    "quoted-hash",
    feature("exon", attrs='gene_id "G#1"; exon_id "E#1"; # comment'),
    [("demo", 0, 10, "G#1", "E#1")],
)
run("crlf", base + "\r\n", [("demo", 0, 10, "G", "E")])
run(
    "nonexons-ignored",
    feature("gene", 1, 10, 'gene_id "G";')
    + "\n"
    + feature("CDS", 1, 10, 'gene_id "G";')
    + "\n"
    + base,
    [("demo", 0, 10, "G", "E")],
)
for name, text, error in [
    ("missing-gene", feature("exon", attrs='transcript_id "T";'), "missing gene_id"),
    ("empty-gene", feature("exon", attrs='gene_id "";'), "Empty"),
    ("dot-exon", feature("exon", attrs='gene_id "G"; exon_id ".";'), "Empty"),
    ("duplicate-gene", feature("exon", attrs='gene_id "G"; gene_id "H";'), "Duplicate"),
    ("quote", feature("exon", attrs='gene_id "G;'), "Unclosed"),
    ("separator", feature("exon", attrs='gene_id "G" exon_id "E";'), "semicolon"),
    ("escape", feature("exon", attrs='gene_id "G\\n";'), "Unsupported GTF escape"),
    ("zero", feature("exon", 0, 10, 'gene_id "G";'), "positive"),
    ("reversed", feature("exon", 20, 10, 'gene_id "G";'), "positive"),
    ("negative", feature("exon", -1, 10, 'gene_id "G";'), "invalid digit"),
    ("overflow", feature("exon", 1, 2**64, 'gene_id "G";'), "too large"),
    ("strand", feature("exon", attrs='gene_id "G";', strand="x"), "strand"),
    ("columns", "demo\texon", "nine"),
    (
        "reserved",
        feature("exon", attrs='gene_id "G"; exon_id "gvcf-audit:exon:0-10:+";'),
        "reserved",
    ),
    ("noexons", feature("CDS", attrs='gene_id "G";'), "No exon"),
    ("trailing-bad", base + "\nnot a feature", "nine"),
]:
    run(name, text, fmt="gtf", error=error)
gene = feature("gene", 1, 100, "ID=G")
tx = feature("mRNA", 1, 100, "ID=T;Parent=G")
ex = feature("exon", 1, 10, "ID=E;Parent=T")
gff = "##gff-version 3\n" + gene + "\n" + tx + "\n" + ex
run("gff-basic", gff, [("demo", 0, 10, "G", "E")], fmt="gff3")
run(
    "gff-forward-parent",
    ex + "\n" + tx + "\n" + gene,
    [("demo", 0, 10, "G", "E")],
    fmt="gff3",
)
run(
    "gff-multipart",
    gff + "\n" + feature("exon", 21, 30, "ID=E;Parent=T"),
    [("demo", 0, 10, "G", "E"), ("demo", 20, 30, "G", "E")],
    fmt="gff3",
)
run(
    "gff-direct-gene",
    gene + "\n" + feature("exon", 1, 10, "Parent=G"),
    [("demo", 0, 10, "G", "gvcf-audit:exon:0-10:+")],
    fmt="gff3",
)
run(
    "gff-fasta", gff + "\n##FASTA\n>demo\nAAAA", [("demo", 0, 10, "G", "E")], fmt="gff3"
)
run(
    "gff-encoding",
    feature("gene", 1, 10, "ID=G%2C1")
    + "\n"
    + feature("exon", 1, 10, "ID=E%3B%3D%25+;Parent=G%2C1"),
    [("demo", 0, 10, "G,1", "E;=%+")],
    fmt="gff3",
)
run(
    "gff-encoded-contig",
    feature("gene", 1, 10, "ID=G", contig="chr%3A1")
    + "\n"
    + feature("exon", 1, 10, "Parent=G", contig="chr%3A1"),
    [("chr:1", 0, 10, "G", "gvcf-audit:exon:0-10:+")],
    fmt="gff3",
)
run(
    "gff-unicode",
    feature("gene", 1, 10, "ID=G%C3%A9")
    + "\n"
    + feature("exon", 1, 10, "ID=E;Parent=G%C3%A9"),
    [("demo", 0, 10, "Gé", "E")],
    fmt="gff3",
)
run(
    "gff-so",
    feature("SO:0000704", 1, 10, "ID=G")
    + "\n"
    + feature("SO:0000147", 1, 10, "ID=E;Parent=G"),
    [("demo", 0, 10, "G", "E")],
    fmt="gff3",
)
for kind in ["pseudogene", "ncRNA_gene", "transposable_element_gene"]:
    run(
        "gff-" + kind,
        feature(kind, 1, 10, "ID=G") + "\n" + feature("exon", 1, 10, "ID=E;Parent=G"),
        [("demo", 0, 10, "G", "E")],
        fmt="gff3",
    )
for name, text, error in [
    ("missing-parent", gene + "\n" + ex, "Missing GFF3 parent T"),
    ("parentless-exon", gene + "\n" + feature("exon", attrs="ID=E"), "missing Parent"),
    (
        "parentless-transcript",
        feature("mRNA", attrs="ID=T") + "\n" + ex,
        "no gene ancestor",
    ),
    (
        "cycle",
        feature("mRNA", attrs="ID=T;Parent=U")
        + "\n"
        + feature("mRNA", attrs="ID=U;Parent=T")
        + "\n"
        + ex,
        "cycle",
    ),
    (
        "conflicting-id",
        gff + "\n" + feature("mRNA", attrs="ID=T;Parent=H"),
        "Conflicting",
    ),
    (
        "duplicate-id",
        gene + "\n" + feature("exon", attrs="ID=E;ID=F;Parent=G"),
        "Duplicate",
    ),
    (
        "percent-escape",
        gene + "\n" + feature("exon", attrs="ID=E%ZZ;Parent=G"),
        "percent",
    ),
    ("utf8", gene + "\n" + feature("exon", attrs="ID=E%FF;Parent=G"), "utf-8"),
    ("control", gene + "\n" + feature("exon", attrs="ID=E%09F;Parent=G"), "control"),
    ("contig-comment", feature("gene", attrs="ID=G", contig="%23chr"), "Contig name"),
    (
        "strand-mismatch",
        gene + "\n" + feature("exon", attrs="Parent=G", strand="-"),
        "incompatible",
    ),
    (
        "strand-bridge",
        gene
        + "\n"
        + feature("mRNA", attrs="ID=T;Parent=G", strand=".")
        + "\n"
        + feature("exon", attrs="Parent=T", strand="-"),
        "incompatible",
    ),
    (
        "cross-contig",
        gene + "\n" + feature("exon", attrs="Parent=G", contig="other"),
        "incompatible",
    ),
    ("region-bounds", "##sequence-region demo 1 5\n" + gff, "sequence-region bounds"),
    (
        "late-region-bounds",
        gff + "\n##sequence-region demo 1 5",
        "sequence-region bounds",
    ),
    ("gff2", "##gff-version 2\n" + gff, "GFF1/GFF2"),
    (
        "multi-id",
        gene + "\n" + feature("exon", attrs="ID=E,F;Parent=G"),
        "single value",
    ),
]:
    run("gff-" + name, text, fmt="gff3", error=error)
run(
    "contig-filter",
    base + "\n" + feature("exon", 1, 20, 'gene_id "H"; exon_id "F";', contig="other"),
    [("demo", 0, 10, "G", "E")],
    extra=["--contig", "demo"],
)
run(
    "gene-filter",
    base + "\n" + feature("exon", 1, 20, 'gene_id "H"; exon_id "F";'),
    [("demo", 0, 10, "G", "E")],
    extra=["--gene", "G"],
)
run("unknown-gene", base, extra=["--gene", "missing"], error="No exon targets match")
run(
    "unknown-contig", base, extra=["--contig", "missing"], error="No exon targets match"
)
run(
    "some-missing-gene",
    base,
    extra=["--gene", "G", "--gene", "missing"],
    error="No exon targets match",
)
run("unknown-format", base, fmt="txt", error="Cannot detect")
run(
    "explicit-format",
    base,
    [("demo", 0, 10, "G", "E")],
    fmt="txt",
    extra=["--format", "gtf"],
)
run("header-detect", gff, [("demo", 0, 10, "G", "E")], fmt="txt")
run("header-conflict", gff, fmt="gtf", error="Expected GFF3")
run(
    "header-tabs",
    gff.replace("##gff-version 3", "##gff-version\t3"),
    [("demo", 0, 10, "G", "E")],
    fmt="txt",
)
run(
    "region-tabs",
    "##sequence-region\tdemo 1 5\n" + gff,
    fmt="gff3",
    error="sequence-region bounds",
)
run("bad-header", "##gff-version 3.bad\n" + gff, fmt="gff3", error="Expected GFF3")
run(
    "gff-gene-filter",
    gene
    + "\n"
    + feature("gene", 1, 100, "ID=H")
    + "\n"
    + feature("exon", 1, 10, "ID=E;Parent=G,H"),
    [("demo", 0, 10, "H", "E")],
    fmt="gff3",
    extra=["--gene", "H"],
)
# Generate equivalent GTF and GFF3 representations from a separate interval model.
seed = random_seed(5021)
rng = random.Random(seed)
for trial in range(30):
    expected = []
    gtf = []
    gff = []
    for gi in range(8):
        g = f"G{gi}"
        c = rng.choice(["demo", "other"])
        strand = rng.choice(["+", "-"])
        t = g + "T"
        u = g + "U"
        gff += [
            feature("gene", 1, 100, "ID=" + g, contig=c, strand=strand),
            feature("mRNA", 1, 100, f"ID={t};Parent={g}", contig=c, strand=strand),
            feature("mRNA", 1, 100, f"ID={u};Parent={g}", contig=c, strand=strand),
        ]
        for ei in range(6):
            s = rng.randrange(100)
            e = rng.randrange(s + 1, 101)
            x = g + f"E{ei}"
            expected.append((c, s, e, g, x))
            gff.append(
                feature(
                    "exon", s + 1, e, f"ID={x};Parent={t},{u}", contig=c, strand=strand
                )
            )
            for transcript in [t, u]:
                gtf.append(
                    feature(
                        "exon",
                        s + 1,
                        e,
                        f'gene_id "{g}"; transcript_id "{transcript}"; exon_id "{x}";',
                        contig=c,
                        strand=strand,
                    )
                )
    rng.shuffle(gtf)
    rng.shuffle(gff)
    run(f"random-gtf-{trial}", "\n".join(gtf), expected)
    run(f"random-gff-{trial}", "\n".join(gff), expected, fmt="gff3")
# One exon shared across genes produces a target for each gene.
run(
    "multi-gene",
    gene
    + "\n"
    + feature("gene", 1, 100, "ID=H")
    + "\n"
    + feature("exon", 1, 10, "ID=E;Parent=G,H"),
    [("demo", 0, 10, "G", "E"), ("demo", 0, 10, "H", "E")],
    fmt="gff3",
)
# Existing destinations (including broken symlinks) are preserved.
for name in ["existing", "dangling"]:
    out = cases / (name + ".tsv")
    if name == "existing":
        out.write_text("preserve me")
    else:
        out.symlink_to(cases / "absent")
    p = subprocess.run(
        [
            str(exe),
            "targets",
            "--annotation",
            str(cases / "gtf-basic.gtf"),
            "--out",
            str(out),
        ],
        capture_output=True,
        text=True,
    )
    assert p.returncode == 1 and "Output already exists" in p.stderr
    assert out.read_text() == "preserve me" if name == "existing" else out.is_symlink()
    checks.append(name)
for name, data, expectedcode in [
    ("gzip", gzip.compress(base.encode()), 0),
    (
        "concat-gzip",
        gzip.compress((base + "\n").encode()) + gzip.compress(base.encode()),
        0,
    ),
    ("truncated", gzip.compress(base.encode())[:-8], 1),
]:
    src = cases / (name + ".gtf.gz")
    src.write_bytes(data)
    out = cases / (name + ".tsv")
    p = subprocess.run(
        [str(exe), "targets", "--annotation", str(src), "--out", str(out)],
        capture_output=True,
        text=True,
    )
    assert p.returncode == expectedcode, (name, p.stderr)
    assert out.exists() == (expectedcode == 0)
    checks.append(name)
(study / "controlled-results.json").write_text(
    json.dumps(
        {
            "checks_passed": len(checks),
            "randomized_equivalent_models": 30,
            "checks": checks,
        },
        indent=2,
    )
    + "\n"
)
print(f"{len(checks)} checks passed")
