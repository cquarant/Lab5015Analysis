const variables={tL:'L',tR:'R',tAve:'DUT mean',tAveRef:'REF mean',tAveRefCor:'corrected REF mean',tLRef:'L REF',tRRef:'R REF',phaseMean:'mean phase',phaseAve:'DUT mean phase',phaseAveRef:'REF mean phase',phaseRefAve:'REF mean phase',phaseL:'L phase',phaseR:'R phase',eRatio:'E(L) / E(R)',eL:'L energy',eR:'R energy',eLRef:'REF L energy',eRRef:'REF R energy',eAveRef:'mean REF energy'};
function describe(name) {
  const [kind,first,...rest]=name.split('_');
  if (first!=='deltaT') return {group:'Energy and phase',label:variables[first]||first};
  const [left,right,...terms]=rest;
  const stages={eCor:'time walk',phaseCor:'phase',eRatioCor:'energy ratio',phaseMeanCor:'mean phase'};
  let label=variables[left]+' − '+variables[right];
  const corrections=terms.filter(t=>stages[t]).map(t=>stages[t]);
  if (corrections.length) label+=' · '+corrections.join(' + ');
  if (terms.includes('vs')) label+=' vs '+(variables[terms[terms.length-1]]||terms[terms.length-1]);
  const group={h1:'Time distributions',p1:'Calibration profiles',h2:'2D correlations'}[kind];
  return {group,label};
}
const defaultQ='h1_deltaT_tL_tR_eRatioCor_phaseMeanCor';
const byId = id => document.getElementById(id);
const ids = ['run','ref','bar','th','method','stage','q'];
const defaults = {run:'4194',ref:'7',bar:'7',th:'10',method:'gaussian',stage:'full',q:defaultQ};
const methods = {direct_halfdiff:'Direct',avg_mixed:'Mixed',matched_halfdiff:'Corrected difference',catr:'CatR'};
const two = n => String(n).padStart(2,'0');
const packCache = {}, summaryCache = {};
let version = 0, currentURL;
async function read(url, type) {
  const response = await fetch(url);
  if (!response.ok) throw new Error('Could not load ' + url);
  return response[type]();
}
function cached(cache, url, type) {
  if (cache.url !== url) {
    cache.url = url;
    cache.data = read(url, type).catch(err => {
      if (cache.url === url) cache.url = null;
      throw err;
    });
  }
  return cache.data;
}
async function render() {
  const n = ++version;
  const s = Object.fromEntries(ids.map(id => [id,byId(id).value]));
  const ref = 'ref'+two(s.ref);
  const group = s.run+'/'+ref;
  const img = byId('diagImg'), download = byId('download');
  img.hidden = download.hidden = true;
  byId('diagStatus').textContent = 'Loading plot…';
  byId('widths').replaceChildren();
  byId('entries').textContent = '';
  const imageTask = async () => {
    const data = await cached(packCache,'diagnostic-images/'+group+'/bar'+two(s.bar)+'.bin','arrayBuffer');
    if (n !== version) return;
    const length = new DataView(data).getUint32(0,true);
    const offsets = JSON.parse(new TextDecoder().decode(new Uint8Array(data,4,length)));
    const slice = offsets[s.th+'/'+s.q];
    if (!slice) throw new Error('Plot unavailable for this selection.');
    if (currentURL) URL.revokeObjectURL(currentURL);
    const start = 4+length+slice[0];
    currentURL = URL.createObjectURL(new Blob([data.slice(start,start+slice[1])],{type:'image/png'}));
    img.src = download.href = currentURL;
    img.alt = describe(s.q).label+', run '+s.run+', '+ref+', DUT '+s.bar+', T1 '+s.th;
    download.download = s.run+'_'+ref+'_'+s.q+'_bar'+two(s.bar)+'_th'+two(s.th)+'.png';
    img.hidden = download.hidden = false;
    byId('diagStatus').textContent = s.q;
  };
  const tableTask = async () => {
    const summary = await cached(summaryCache,'data/'+s.run+'_'+ref+'_summary.json','json');
    if (n !== version) return;
    const g = summary.groups.find(g => g.bar === Number(s.bar) && g.threshold === Number(s.th));
    if (!g) throw new Error('Closure replay unavailable for this selection.');
    const result = g.stages[s.stage][s.method];
    for (const [key,name] of Object.entries(methods)) {
      const row = byId('widths').insertRow();
      row.insertCell().textContent = name;
      const variance = result?.[key];
      row.insertCell().textContent = Number.isFinite(variance) && variance >= 0 ? Math.sqrt(variance).toFixed(2) : '—';
    }
    byId('entries').textContent = g.counts.accepted.toLocaleString('en-GB')+' accepted coincidences. A dash indicates an unavailable width.';
  };
  await Promise.all([[imageTask,'diagStatus'],[tableTask,'entries']].map(async ([task,id]) => {
    try { await task(); } catch (err) { if (n === version) byId(id).textContent = err.message; }
  }));
}
function restore() {
  const query = new URLSearchParams(location.hash.split('?')[1] || location.search);
  for (const id of ids) {
    const field = byId(id);
    field.value = query.get(id) || defaults[id];
    if (!field.value) field.value = defaults[id];
  }
  render();
}
async function init() {
  const options = {run:[4194,4257,4281,4291],ref:[4,5,7,8,9,10,11,12],bar:Array.from({length:16},(_,i)=>i),th:[3,5,7,10,15,20]};
  for (const [id,values] of Object.entries(options)) byId(id).replaceChildren(...values.map(v=>new Option(v,v)));
  const quantities = await read('diagnostic-images/quantities.json','json');
  for (const group of ['Time distributions','Calibration profiles','2D correlations','Energy and phase']) {
    const choices = document.createElement('optgroup');
    choices.label = group;
    for (const name of quantities) {
      const d = describe(name);
      if (d.group === group) choices.append(new Option(d.label,name));
    }
    byId('q').append(choices);
  }
  for (const id of ids) byId(id).addEventListener('change',() => {
    const params = new URLSearchParams(ids.map(k => [k,byId(k).value]));
    history.replaceState(null,'','#diagnostics?'+params);
    render();
  });
  window.addEventListener('hashchange',restore);
  restore();
}
init().catch(err => {byId('diagStatus').textContent = err.message;});
