const data = JSON.parse(document.getElementById('comparison-data').textContent);
for (const kind of ['genes', 'exons']) {
  const search = document.getElementById(`${kind}-search`);
  const sort = document.getElementById(`${kind}-sort`);
  const body = document.getElementById(`${kind}-rows`);
  const previous = document.getElementById(`${kind}-prev`);
  const next = document.getElementById(`${kind}-next`);
  let selected = [];
  let page = 0;
  const pageSize = 50;
  function render() {
    body.replaceChildren();
    const first = page * pageSize;
    for (const row of selected.slice(first, first + pageSize)) {
      const tr = document.createElement('tr');
      const cells = [row.exon ? `${row.gene} / ${row.exon}` : row.gene, row.contig,
        `${(100 * row.before_callable_fraction).toFixed(2)}%`,
        `${(100 * row.after_callable_fraction).toFixed(2)}%`,
        row.gained_callable_bases, row.lost_callable_bases, row.net_callable_bases];
      for (const value of cells) {
        const td = document.createElement('td');
        td.textContent = value;
        tr.append(td);
      }
      body.append(tr);
    }
    document.getElementById(`${kind}-count`).textContent = selected.length
      ? `${first + 1}–${Math.min(first + pageSize, selected.length)} of ${selected.length} ${kind}.`
      : data[kind].length ? `No matching ${kind}.`
      : `No ${kind} summaries. Create both audits with --gene-targets to compare them.`;
    previous.disabled = page === 0;
    next.disabled = first + pageSize >= selected.length;
  }
  function update() {
    const query = search.value.toLowerCase().trim();
    selected = data[kind].filter(row => `${row.gene} ${row.exon ?? ''} ${row.contig}`.toLowerCase().includes(query));
    selected.sort((a, b) => b[sort.value] - a[sort.value]);
    page = 0;
    render();
  }
  search.addEventListener('input', update);
  sort.addEventListener('change', update);
  previous.addEventListener('click', () => { page--; render(); });
  next.addEventListener('click', () => { page++; render(); });
  update();
}
