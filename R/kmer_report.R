#' Write a self-contained interactive HTML report
#'
#' Bundles the admitted primers, their genome positions per virus, the recommended
#' pool and the saturation curve into one standalone HTML file: a ranked primer
#' table (enrichment bars, multi-select) linked to an interactive genome track
#' (hover a row or a mark for details).
#'
#' @param targets Named list of [Biostrings::DNAStringSet] (one per virus).
#' @param admit Admitted candidates (data.frame with `pattern` and `enrichment`).
#' @param cov A [coverage_matrix()] result — the report solves the balanced optimal
#'   pool ([optimize_pool_balanced()]) for every size N from it.
#' @param read_len Effective read length used for genome positions.
#' @param file Output HTML path.
#' @param min_sites Drop primers binding fewer than this many sites across the panel
#'   (default 2 — single-site primers cover a lone locus and add little).
#' @param floor Per-target coverage threshold: any coverage below it is flagged red.
#' @param k_max Largest pool size N offered by the selector.
#' @return The output path, invisibly.
#' @export
kmer_report <- function(targets, admit, cov, read_len = 200, file = "kmer_report.html",
                        min_sites = 2L, floor = 0.5, k_max = 6L) {
  admit <- admit[order(-admit$enrichment), ]

  vinfo <- lapply(names(targets), function(v)
    list(name = v, length = as.integer(Biostrings::width(targets[[v]][1])[1])))

  pos_by_virus <- lapply(names(targets), function(v) {
    mc <- map_coverage(admit$pattern, targets[[v]][1], read_len)
    split(mc[, c("start", "cov_start", "cov_end")], mc$pattern)
  })
  names(pos_by_virus) <- names(targets)

  primers <- lapply(seq_len(nrow(admit)), function(i) {
    pat <- admit$pattern[i]
    pos <- lapply(names(targets), function(v) {
      d <- pos_by_virus[[v]][[pat]]
      if (is.null(d)) list()
      else lapply(seq_len(nrow(d)), function(j) list(b = d$start[j], s = d$cov_start[j], e = d$cov_end[j]))
    })
    names(pos) <- names(targets)
    list(id = i, pattern = pat,
         order = as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(pat))),
         enrichment = round(admit$enrichment[i], 2), pos = pos)
  })

  sites <- vapply(primers, function(p) sum(vapply(p$pos, length, integer(1))), integer(1))
  primers <- primers[sites >= min_sites]

  # restrict the coverage matrix to the displayed primers so the pools stay consistent
  keep <- match(vapply(primers, function(p) p$pattern, character(1)), cov$patterns)
  keep <- keep[!is.na(keep)]
  cov$matrix <- cov$matrix[, keep, drop = FALSE]
  cov$patterns <- cov$patterns[keep]

  # balanced optimal pool for every size N (always feasible, both genomes in focus)
  sb <- saturation_balanced(cov, k_max)
  pools <- lapply(sb$pools, function(r)
    list(k = r$k, patterns = r$patterns, coverage = round(r$coverage, 3),
         per_target = as.list(round(r$per_target, 3)), min_cov = round(r$min_cov, 3)))

  data <- list(
    viruses = vinfo, primers = primers, pools = pools, floor = floor,
    pool = list(n_star = sb$n_star),
    saturation = lapply(seq_len(nrow(sb$curve)), function(i)
      list(k = sb$curve$k[i], coverage = sb$curve$coverage[i], min_cov = sb$curve$min_cov[i])))

  json <- jsonlite::toJSON(data, auto_unbox = TRUE, null = "null", na = "null", digits = 6)
  parts <- strsplit(report_template(), "__DATA_JSON__", fixed = TRUE)[[1]]
  writeLines(paste0(parts[1], json, parts[2]), file)
  invisible(file)
}

report_template <- function() {
r"---(<!doctype html>
<html lang="fr"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>KMerEnrich — rapport de design</title>
<style>
:root{--ink:#1a2027;--muted:#6b7784;--line:#e2e6ea;--bg:#f6f8fa;--card:#fff;--accent:#2a78d6}
*{box-sizing:border-box}
body{margin:0;font:14px/1.5 -apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;color:var(--ink);background:var(--bg)}
.wrap{max-width:1180px;margin:0 auto;padding:24px 18px 60px}
h1{font-size:22px;margin:0 0 4px}h2{font-size:15px;margin:0 0 10px}
.sub{color:var(--muted);margin:0 0 20px}
.cards{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin-bottom:18px}
.card{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:12px 14px}
.card .k{color:var(--muted);font-size:12px}.card .v{font-size:20px;font-weight:600;margin-top:2px}
.card .v.small{font-size:13px;font-weight:500}
.panel{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:16px;margin-bottom:16px}
.row{display:flex;align-items:center;gap:10px}.between{justify-content:space-between}
.btns button,.tools button,select,input{font:inherit;border:1px solid var(--line);background:#fff;border-radius:7px;padding:5px 10px;cursor:pointer}
.btns button:hover,.tools button:hover{border-color:var(--accent);color:var(--accent)}
.tools{display:flex;gap:8px;margin:-6px 0 18px}
select,input{cursor:auto}
.trackbox{margin-top:12px;width:100%;overflow:hidden}
svg.track{width:100%;display:block}
.tick{font-size:9px;fill:var(--muted)}
.glabel{font-size:11px;fill:#37474f;font-family:ui-monospace,Menlo,monospace;font-weight:600}
.legend{display:flex;flex-wrap:wrap;gap:12px;margin-top:8px}
.lg{display:inline-flex;align-items:center;gap:6px;font-family:ui-monospace,Menlo,monospace;font-size:12px}
.sw{width:12px;height:12px;border-radius:3px;display:inline-block}
.grid2{display:grid;grid-template-columns:1.4fr 1fr;gap:16px}
table{width:100%;border-collapse:collapse;font-size:13px}
th,td{text-align:left;padding:6px 8px;border-bottom:1px solid var(--line)}
th{color:var(--muted);font-weight:500;font-size:12px;position:sticky;top:0;background:#fff}
.tblbox{max-height:440px;overflow:auto}
.mono{font-family:ui-monospace,Menlo,monospace}.muted{color:var(--muted)}.num{text-align:right}
tr[data-id]:hover{background:#f0f6ff}
.enr{min-width:150px}.enr .bar{display:inline-block;height:9px;border-radius:3px;vertical-align:middle;margin-right:6px;min-width:2px}
.enr .val{font-variant-numeric:tabular-nums;font-size:12px}
.badge{background:#e7f0ff;color:var(--accent);border-radius:5px;padding:1px 7px;font-size:11px;font-weight:600}
svg.satsvg{width:100%;max-width:460px;height:auto}
.satpt{cursor:pointer}
.tip{position:fixed;display:none;pointer-events:none;background:#263238;color:#fff;padding:6px 9px;border-radius:6px;font-size:12px;z-index:9;box-shadow:0 4px 14px rgba(0,0,0,.2)}
.mark{cursor:pointer}
.lo{color:#e34948;font-weight:600}
.satleg{display:flex;gap:14px;margin-top:6px;font-size:11px;color:var(--muted)}
.satleg i{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:5px;vertical-align:middle}
@media(max-width:820px){.cards{grid-template-columns:repeat(2,1fr)}.grid2{grid-template-columns:1fr}}
</style></head>
<body><div class="wrap">
<h1>KMerEnrich — rapport de design des primers</h1>
<p class="sub">Sélectionne des primers dans la table pour les voir sur le génome. Survole une ligne ou un marqueur pour les détails.</p>
<div class="tools"><button id="dlPool">⬇ Pool recommandé (CSV)</button><button id="dlAll">⬇ Tous les primers (CSV)</button></div>
<div class="cards">
 <div class="card"><div class="k">Virus</div><div class="v small" id="sVir"></div></div>
 <div class="card"><div class="k">Taille du pool (N)</div><div class="v" id="sPool"></div></div>
 <div class="card"><div class="k">Couverture globale</div><div class="v" id="sCov"></div></div>
 <div class="card"><div class="k">Couverture par virus</div><div class="v small" id="sPer"></div></div>
</div>
<div class="panel">
 <div class="row between">
  <div class="muted">Tous les génomes de référence (un par virus)</div>
  <div class="btns"><button id="bPool">Pool recommandé (N*)</button><button id="bAll">Tout</button><button id="bNone">Aucun</button></div>
 </div>
 <div class="row" style="margin-top:12px;gap:14px">
  <label for="nsel">Taille du pool&nbsp;N</label>
  <input type="range" id="nsel" min="1" step="1" style="flex:1;max-width:340px">
  <b id="nval"></b><span class="muted" id="nnote"></span>
 </div>
 <div id="track" class="trackbox"></div>
 <div id="legend" class="legend"></div>
</div>
<div class="grid2">
 <div class="panel">
  <div class="row between"><h2>Meilleurs primers (par enrichissement)</h2><input id="search" placeholder="filtrer une séquence…"></div>
  <div class="tblbox"><table id="tbl"><thead><tr id="thr"></tr></thead><tbody></tbody></table></div>
 </div>
 <div class="panel"><h2>Courbe de saturation</h2><div id="sat"></div>
  <p class="sub" style="margin:8px 0 0">Clique un point pour choisir N. Bleu = couverture globale, orange = virus le moins couvert, rouge = seuil.</p></div>
</div>
</div><div id="tip" class="tip"></div>
<script>
const D = __DATA_JSON__;
const PAL=["#2a78d6","#1baf7a","#eda100","#008300","#4a3aa7","#e34948","#e87ba4","#eb6834"];
const primers=D.primers.slice().sort((a,b)=>b.enrichment-a.enrichment);
primers.forEach((p,i)=>{p.rank=i+1;p.color=PAL[i%PAL.length];});
const byId={};primers.forEach(p=>byId[p.id]=p);
const maxEnr=Math.max.apply(null,primers.map(p=>p.enrichment).concat([1]));
const idByPat={};primers.forEach(p=>idByPat[p.pattern]=p.id);
const pools=D.pools||[];const KMAX=Math.max(1,pools.length);
function poolPatterns(n){return new Set((pools[n-1]||{patterns:[]}).patterns);}
let curN=Math.min(Math.max(D.pool.n_star||1,1),KMAX);
let curPool=poolPatterns(curN);
let sel=new Set([...curPool].map(p=>idByPat[p]).filter(x=>x!=null));
let curV=D.viruses[0].name;
const tip=document.getElementById('tip');
function showTip(h,x,y){tip.innerHTML=h;tip.style.display='block';tip.style.left=(x+12)+'px';tip.style.top=(y+12)+'px';}
function hideTip(){tip.style.display='none';}
function vlen(n){return D.viruses.find(v=>v.name===n).length;}
const W=1000,P=10,LH=86;
function renderTrack(hi){
 const ids=new Set(sel);if(hi!=null)ids.add(hi);
 const VS=D.viruses,H=VS.length*LH;
 let s='<svg viewBox="0 0 '+W+' '+H+'" preserveAspectRatio="none" class="track" style="height:'+H+'px">';
 VS.forEach((v,i)=>{
  const y0=i*LH,g=v.length,xs=pos=>P+(pos/g)*(W-2*P);
  s+='<text x="'+P+'" y="'+(y0+16)+'" class="glabel">'+v.name+' — '+g+' bp</text>';
  s+='<rect x="'+P+'" y="'+(y0+30)+'" width="'+(W-2*P)+'" height="16" rx="3" fill="#eceff1" stroke="#cfd8dc"/>';
  for(let t=0;t<=5;t++){const pos=Math.round(g*t/5),x=xs(pos);
   s+='<line x1="'+x+'" y1="'+(y0+50)+'" x2="'+x+'" y2="'+(y0+56)+'" stroke="#b0bec5"/><text x="'+x+'" y="'+(y0+68)+'" class="tick" text-anchor="middle">'+pos+'</text>';}
  primers.filter(p=>ids.has(p.id)).forEach(p=>{(p.pos[v.name]||[]).forEach(q=>{
   const a=xs(q.s),b=xs(q.e),xb=xs(q.b),op=hi==null?0.92:(hi===p.id?1:0.15);
   s+='<rect class="mark" data-id="'+p.id+'" data-v="'+v.name+'" x="'+Math.min(a,b)+'" y="'+(y0+32)+'" width="'+Math.max(2,Math.abs(b-a))+'" height="12" rx="2" fill="'+p.color+'" opacity="'+op+'"></rect>';
   s+='<circle class="mark" data-id="'+p.id+'" data-v="'+v.name+'" cx="'+xb+'" cy="'+(y0+38)+'" r="2.4" fill="#263238" opacity="'+op+'"></circle>';});});
 });
 s+='</svg>';
 const el=document.getElementById('track');el.innerHTML=s;
 el.querySelectorAll('.mark').forEach(m=>{const p=byId[+m.dataset.id];
  m.addEventListener('mousemove',e=>showTip('<b>'+p.pattern+'</b><br>enrichissement ×'+p.enrichment+'<br>sur '+m.dataset.v,e.clientX,e.clientY));
  m.addEventListener('mouseleave',hideTip);});
}
function renderLegend(){document.getElementById('legend').innerHTML=
 primers.filter(p=>sel.has(p.id)).map(p=>'<span class="lg"><span class="sw" style="background:'+p.color+'"></span>'+p.pattern+'</span>').join('')||'<span class="muted">aucun primer sélectionné</span>';}
function renderTable(f){
 const q=(f||'').toUpperCase();
 document.querySelector('#tbl tbody').innerHTML=primers.filter(p=>!q||p.pattern.indexOf(q)>=0).map(p=>{
  const siteCells=D.viruses.map(v=>'<td class="num">'+((p.pos[v.name]||[]).length)+'</td>').join('');
  const w=(p.enrichment/maxEnr*100).toFixed(0);
  return '<tr data-id="'+p.id+'"><td><input type="checkbox" data-id="'+p.id+'"'+(sel.has(p.id)?' checked':'')+'></td>'+
   '<td>'+p.rank+'</td><td class="mono">'+p.pattern+'</td><td class="mono muted">'+p.order+'</td>'+
   '<td class="enr"><span class="bar" style="width:'+w+'%;background:'+p.color+'"></span><span class="val">×'+p.enrichment+'</span></td>'+
   siteCells+'<td>'+(curPool.has(p.pattern)?'<span class="badge">pool</span>':'')+'</td></tr>';}).join('');
 const tb=document.querySelector('#tbl tbody');
 tb.querySelectorAll('input[type=checkbox]').forEach(c=>c.addEventListener('change',()=>{
  const id=+c.dataset.id;c.checked?sel.add(id):sel.delete(id);renderTrack();renderLegend();}));
 tb.querySelectorAll('tr').forEach(tr=>{tr.addEventListener('mouseenter',()=>renderTrack(+tr.dataset.id));
  tr.addEventListener('mouseleave',()=>renderTrack());});
}
function renderSat(){
 const c=D.saturation;const box=document.getElementById('sat');
 if(!c.length){box.innerHTML='<p class="muted">non disponible</p>';return;}
 const w=460,h=200,pl=34;
 const mk=Math.max.apply(null,c.map(d=>d.k));
 const mc=Math.max.apply(null,c.map(d=>Math.max(d.coverage||0,d.min_cov||0)).concat([0.01,D.floor]));
 const sx=k=>pl+(mk>1?(k-1)/(mk-1):0)*(w-pl-14),sy=v=>h-28-((v||0)/mc)*(h-28-12);
 let s='<svg viewBox="0 0 '+w+' '+h+'" class="satsvg">';
 s+='<line x1="'+pl+'" y1="'+(h-28)+'" x2="'+(w-6)+'" y2="'+(h-28)+'" stroke="#cfd8dc"/><line x1="'+pl+'" y1="10" x2="'+pl+'" y2="'+(h-28)+'" stroke="#cfd8dc"/>';
 s+='<line x1="'+pl+'" y1="'+sy(D.floor)+'" x2="'+(w-6)+'" y2="'+sy(D.floor)+'" stroke="#e34948" stroke-dasharray="3 3" opacity="0.55"/><text x="'+(w-6)+'" y="'+(sy(D.floor)-3)+'" class="tick" text-anchor="end" fill="#e34948">seuil '+(D.floor*100).toFixed(0)+'%</text>';
 s+='<line x1="'+sx(curN)+'" y1="10" x2="'+sx(curN)+'" y2="'+(h-28)+'" stroke="#90a4ae" stroke-dasharray="4 3"/><text x="'+(sx(curN)+3)+'" y="18" class="tick">N='+curN+'</text>';
 s+='<polyline points="'+c.map(d=>sx(d.k)+','+sy(d.min_cov)).join(' ')+'" fill="none" stroke="#eda100" stroke-width="1.6"/>';
 s+='<polyline points="'+c.map(d=>sx(d.k)+','+sy(d.coverage)).join(' ')+'" fill="none" stroke="#2a78d6" stroke-width="2"/>';
 c.forEach(d=>{const on=d.k===curN;s+='<circle class="satpt" cx="'+sx(d.k)+'" cy="'+sy(d.coverage)+'" r="'+(on?5:3.6)+'" fill="#2a78d6" stroke="#fff" stroke-width="'+(on?2:0)+'" data-k="'+d.k+'" data-c="'+(d.coverage*100).toFixed(0)+'" data-m="'+(d.min_cov*100).toFixed(0)+'"/>';});
 s+='<text x="4" y="14" class="tick">'+(mc*100).toFixed(0)+'%</text><text x="4" y="'+(h-30)+'" class="tick">0%</text>';
 for(let k=1;k<=mk;k++){s+='<text x="'+sx(k)+'" y="'+(h-12)+'" class="tick" text-anchor="middle">'+k+'</text>';}
 s+='</svg><div class="satleg"><span><i style="background:#2a78d6"></i>couverture globale</span><span><i style="background:#eda100"></i>virus le moins couvert</span></div>';
 box.innerHTML=s;
 box.querySelectorAll('.satpt').forEach(pt=>{pt.style.cursor='pointer';
  pt.addEventListener('click',()=>setN(+pt.dataset.k));
  pt.addEventListener('mousemove',e=>showTip('N = '+pt.dataset.k+'<br>global '+pt.dataset.c+'%<br>pire virus '+pt.dataset.m+'%',e.clientX,e.clientY));
  pt.addEventListener('mouseleave',hideTip);});
}
function sync(){renderTrack();renderLegend();renderTable(document.getElementById('search').value);}
document.getElementById('thr').innerHTML='<th></th><th>#</th><th>Primer</th><th>À commander</th><th>Enrichissement ×N</th>'+D.viruses.map(v=>'<th class="num">Sites '+v.name+'</th>').join('')+'<th></th>';
document.getElementById('bPool').onclick=()=>setN(D.pool.n_star||1);
document.getElementById('bAll').onclick=()=>{sel=new Set(primers.map(p=>p.id));sync();};
document.getElementById('bNone').onclick=()=>{sel=new Set();sync();};
document.getElementById('search').addEventListener('input',e=>renderTable(e.target.value));
function prows(only){return primers.filter(p=>!only||p.inPool).map(p=>({primer:p.pattern,a_commander:p.order,enrichissement:p.enrichment,sites:D.viruses.reduce((n,v)=>n+((p.pos[v.name]||[]).length),0),pool:p.inPool?'oui':'non'}));}
function toCSV(rows){const c=Object.keys(rows[0]||{primer:''});return [c.join(',')].concat(rows.map(r=>c.map(k=>r[k]).join(','))).join('\n');}
function dlCSV(name,rows){const b=new Blob([toCSV(rows)],{type:'text/csv'}),a=document.createElement('a');a.href=URL.createObjectURL(b);a.download=name;document.body.appendChild(a);a.click();a.remove();}
document.getElementById('dlPool').onclick=()=>dlCSV('pool_recommande.csv',prows(true));
document.getElementById('dlAll').onclick=()=>dlCSV('primers.csv',prows(false));
function updateCards(){
 const Pp=pools[curN-1]||{coverage:null,per_target:{}};
 document.getElementById('sVir').textContent=D.viruses.map(v=>v.name).join(', ');
 document.getElementById('sPool').textContent=curN+(curN===D.pool.n_star?' (N*)':'');
 document.getElementById('sCov').textContent=Pp.coverage==null?'—':(Pp.coverage*100).toFixed(0)+'%';
 document.getElementById('sPer').innerHTML=Object.keys(Pp.per_target).map(k=>{const val=Pp.per_target[k],str=(val==null?'—':(val*100).toFixed(0)+'%');return k+': '+(val!=null&&val<D.floor?'<span class="lo">'+str+'</span>':str);}).join(' · ')||'—';
}
function setN(n){curN=Math.min(Math.max(n,1),KMAX);curPool=poolPatterns(curN);
 sel=new Set([...curPool].map(p=>idByPat[p]).filter(x=>x!=null));
 document.getElementById('nsel').value=curN;document.getElementById('nval').textContent='N = '+curN;updateCards();renderSat();sync();}
const nsel=document.getElementById('nsel');nsel.max=KMAX;nsel.value=curN;
nsel.addEventListener('input',e=>setN(+e.target.value));
document.getElementById('nnote').textContent='N* recommandé = '+(D.pool.n_star||'—');
document.getElementById('nval').textContent='N = '+curN;
updateCards();renderTrack();renderLegend();renderTable('');renderSat();
</script></body></html>
)---"
}
