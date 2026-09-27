(() => {
  const rows = JSON.parse(document.getElementById('region-data').textContent);
  const search = document.getElementById('region-search');
  const sort = document.getElementById('region-sort');
  const body = document.getElementById('region-rows');
  const previous = document.getElementById('region-prev');
  const next = document.getElementById('region-next');
  const fraction = row => row[4] / (row[2] - row[1]);
  let selected = rows;
  let page = 0;
  const pageSize = 50;

  function render() {
    body.replaceChildren();
    const first = page * pageSize;
    for (const row of selected.slice(first, first + pageSize)) {
      const tr = document.createElement('tr');
      const cells = [row[3], `${row[0]}:${row[1]}–${row[2]}`,
        row[2] - row[1], `${(100 * fraction(row)).toFixed(2)}%`];
      cells.forEach((value, index) => {
        const td = document.createElement('td');
        td.textContent = value;
        if (index >= 2) td.className = 'number';
        tr.append(td);
      });
      body.append(tr);
    }
    document.getElementById('region-count').textContent = selected.length
      ? `${first + 1}–${Math.min(first + pageSize, selected.length)} of ${selected.length} regions.`
      : 'No matching regions.';
    previous.disabled = page === 0;
    next.disabled = first + pageSize >= selected.length;
  }

  function update() {
    const query = search.value.toLowerCase().trim();
    selected = rows.filter(row => `${row[3]} ${row[0]}`.toLowerCase().includes(query));
    if (sort.value !== 'input') {
      const direction = sort.value === 'lowest' ? 1 : -1;
      selected.sort((a, b) => direction * (fraction(a) - fraction(b)));
    }
    page = 0;
    render();
  }

  search.addEventListener('input', update);
  sort.addEventListener('change', update);
  previous.addEventListener('click', () => { page--; render(); });
  next.addEventListener('click', () => { page++; render(); });
  render();
})();

for (const kind of ['gene', 'exon']) {
  const search = document.getElementById(`${kind}-search`);
  if (!search) continue;
  const rows = JSON.parse(document.getElementById(`${kind}-data`).textContent);
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
      const cells = kind === 'gene'
        ? [row.gene, row.contig, row.targeted_exons, row.unresolved_exons]
        : [row.gene, row.exon, row.intervals.map(([s, e]) => `${row.contig}:${s}–${e}`).join(', ')];
      cells.push(row.bases, row.unresolved_bases,
        `${(100 * row.callable_fraction).toFixed(2)}%`, row.main_reason?.replaceAll('_', ' ') ?? 'None');
      const tr = document.createElement('tr');
      for (const value of cells) {
        const td = document.createElement('td');
        td.textContent = value;
        if (typeof value === 'number') td.className = 'number';
        tr.append(td);
      }
      body.append(tr);
    }
    document.getElementById(`${kind}-count`).textContent = selected.length
      ? `${first + 1}–${Math.min(first + pageSize, selected.length)} of ${selected.length} ${kind}s.`
      : `No matching ${kind}s.`;
    previous.disabled = page === 0;
    next.disabled = first + pageSize >= selected.length;
  }

  function update() {
    const query = search.value.toLowerCase().trim();
    selected = rows.filter(row => `${row.gene} ${row.exon ?? ''} ${row.contig}`.toLowerCase().includes(query));
    if (sort.value !== 'input') {
      const direction = sort.value === 'lowest' ? 1 : -1;
      selected.sort((a, b) => direction * (a.callable_fraction - b.callable_fraction));
    }
    page = 0;
    render();
  }
  search.addEventListener('input', update);
  sort.addEventListener('change', update);
  previous.addEventListener('click', () => { page--; render(); });
  next.addEventListener('click', () => { page++; render(); });
  update();
}
