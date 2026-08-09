import { generateOligos, generateStructuredOligos, logoMatrix, logoStackBases, parseSequenceFile, toFasta, toJson, toTsv } from "./core.js";

document.querySelectorAll("[data-common]").forEach(target => target.append(document.querySelector("#common-template").content.cloneNode(true)));
const panels = document.querySelectorAll(".panel");
document.querySelectorAll(".tools button").forEach(button => button.addEventListener("click", () => {
  document.querySelectorAll(".tools button").forEach(item => {
    item.classList.toggle("active", item === button);
    item.setAttribute("aria-pressed", String(item === button));
  });
  panels.forEach(panel => panel.classList.toggle("active", panel.id === button.dataset.panel));
  document.querySelector("#results").classList.add("hidden");
}));

const number = (data, key) => data.get(key) === "" ? null : Number(data.get(key));
function commonOptions(form) {
  const data = new FormData(form);
  return { count:number(data,"count"),seed:data.get("seed"),prefix:data.get("prefix"),fivePrime:data.get("fivePrime").trim(),threePrime:data.get("threePrime").trim(),fiveRandomLength:number(data,"fiveRandomLength"),threeRandomLength:number(data,"threeRandomLength"),sameRandom:data.has("sameRandom"),deduplicate:data.has("deduplicate"),minStem:number(data,"minStem"),minLoop:number(data,"minLoop"),maxLoop:number(data,"maxLoop"),minHpRun:number(data,"minHpRun"),minOverlap:number(data,"minOverlap")};
}

function download(content, filename, type) {
  const url = URL.createObjectURL(new Blob([content], {type})), anchor = document.createElement("a");
  anchor.href = url; anchor.download = filename; anchor.click(); setTimeout(() => URL.revokeObjectURL(url), 1000);
}
function showEntries(entries, form, basename) {
  const result = document.querySelector("#results"), formats = [];
  result.className = "result";
  result.removeAttribute("role");
  if (form.elements.fasta.checked) formats.push(["FASTA", toFasta(entries), `${basename}.fasta`, "text/plain"]);
  if (form.elements.json.checked) formats.push(["JSON", toJson(entries), `${basename}.json`, "application/json"]);
  if (form.elements.tsv.checked) formats.push(["TSV", toTsv(entries), `${basename}.tsv`, "text/tab-separated-values"]);
  const rows = entries.slice(0, 100).map(item => `<tr><td>${escapeHtml(item.name)}</td><td class="sequence" title="${item.sequence}">${item.sequence}</td><td>${item.length}</td><td>${(item.gc_content*100).toFixed(1)}%</td><td>${item.tm == null ? "—" : item.tm.toFixed(1)}</td><td>${item.has_hairpin ? "Yes" : "No"}</td><td>${item.has_homopolymer ? "Yes" : "No"}</td></tr>`).join("");
  result.innerHTML = `<div class="result-head"><div><p class="eyebrow">Generated locally</p><h2>${entries.length} sequences</h2></div><div class="download-links"></div></div><table><thead><tr><th>Name</th><th>Sequence</th><th>Length</th><th>GC</th><th>Tm °C</th><th>Hairpin</th><th>Homopolymer</th></tr></thead><tbody>${rows}</tbody></table>${entries.length>100?"<p>Preview limited to the first 100 sequences. Downloads contain every entry.</p>":""}`;
  for (const [label,content,filename,type] of formats) { const button=document.createElement("button"); button.textContent=`Download ${label}`; button.addEventListener("click",()=>download(content,filename,type)); result.querySelector(".download-links").append(button); }
  result.classList.remove("hidden"); result.scrollIntoView({behavior:window.matchMedia("(prefers-reduced-motion: reduce)").matches?"auto":"smooth"});
}
const escapeHtml = value => String(value).replace(/[&<>"']/g, char => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"})[char]);
function handleError(error) {
  const result = document.querySelector("#results");
  result.className = "result error";
  result.setAttribute("role", "alert");
  result.textContent = error instanceof Error ? error.message : String(error);
  result.scrollIntoView({behavior: window.matchMedia("(prefers-reduced-motion: reduce)").matches ? "auto" : "smooth"});
}

document.querySelector("#random-form").addEventListener("submit", event => { event.preventDefault(); try { const options={...commonOptions(event.currentTarget),length:number(new FormData(event.currentTarget),"length")}; showEntries(generateOligos(options),event.currentTarget,"oligos"); } catch(error){handleError(error);} });
document.querySelector("#structured-form").addEventListener("submit", event => { event.preventDefault(); try { const data=new FormData(event.currentTarget), options={...commonOptions(event.currentTarget),type:data.get("type"),halfLength:number(data,"halfLength"),outerArmLength:number(data,"outerArmLength"),innerHalfLength:number(data,"innerHalfLength"),spacerLength:number(data,"spacerLength")}; showEntries(generateStructuredOligos(options),event.currentTarget,"structured-oligos"); } catch(error){handleError(error);} });

const dropZone=document.querySelector("#drop-zone"), fileInput=document.querySelector("#logo-file");
for (const eventName of ["dragenter","dragover"]) dropZone.addEventListener(eventName,event=>{event.preventDefault();dropZone.classList.add("drag")});
for (const eventName of ["dragleave","drop"]) dropZone.addEventListener(eventName,event=>{event.preventDefault();dropZone.classList.remove("drag")});
dropZone.addEventListener("drop",event=>{if(event.dataTransfer.files.length){fileInput.files=event.dataTransfer.files;dropZone.querySelector("strong").textContent=event.dataTransfer.files[0].name;}});
fileInput.addEventListener("change",()=>{if(fileInput.files[0])dropZone.querySelector("strong").textContent=fileInput.files[0].name;});

function renderLogo(sequences, options) {
  const matrix=logoMatrix(sequences,options.logoType), cssWidth=Math.max(700,Number(options.width)*96), cssHeight=Math.max(240,Number(options.height)*96), margin={top:options.title?48:20,right:20,bottom:45,left:58}, plotWidth=cssWidth-margin.left-margin.right, plotHeight=cssHeight-margin.top-margin.bottom;
  const max=options.logoType==="counts"?sequences.length:options.logoType==="probability"?1:2, colors=options.colorScheme==="base_pairing"?{A:"#e69f00",T:"#e69f00",C:"#0072b2",G:"#0072b2"}:options.colorScheme==="najafabadi"?{A:"#00a08a",C:"#f2ad00",G:"#f98400",T:"#5bbcD6"}:{A:"#2e9f57",C:"#2474b5",G:"#e89b22",T:"#d9473f"};
  let body=`<rect width="100%" height="100%" fill="white"/><line x1="${margin.left}" y1="${margin.top+plotHeight}" x2="${cssWidth-margin.right}" y2="${margin.top+plotHeight}" stroke="#17231f"/><line x1="${margin.left}" y1="${margin.top}" x2="${margin.left}" y2="${margin.top+plotHeight}" stroke="#17231f"/>`;
  if(options.title)body+=`<text x="${cssWidth/2}" y="28" text-anchor="middle" font-family="Georgia" font-size="20">${escapeHtml(options.title)}</text>`;
  matrix.forEach((column,i)=>{const x=margin.left+i*plotWidth/matrix.length, w=Math.max(1,plotWidth/matrix.length*.9); let y=margin.top+plotHeight; for(const base of logoStackBases(column,options.stackOrder)){const h=column[base]/max*plotHeight;if(h>.3){y-=h;body+=`<rect x="${x}" y="${y}" width="${w}" height="${h}" fill="${colors[base]}"/><text x="${x+w/2}" y="${y+h/2}" text-anchor="middle" dominant-baseline="central" fill="white" font-family="Arial" font-weight="700" font-size="${Math.min(w*.7,h*.7)}">${base}</text>`;}} if(matrix.length<=60)body+=`<text x="${x+w/2}" y="${margin.top+plotHeight+18}" text-anchor="middle" font-size="9">${i+1}</text>`;});
  body+=`<text transform="translate(16 ${margin.top+plotHeight/2}) rotate(-90)" text-anchor="middle" font-size="12">${options.logoType==="information"?"Information (bits)":options.logoType==="probability"?"Probability":"Count"}</text><text x="${margin.left+plotWidth/2}" y="${cssHeight-5}" text-anchor="middle" font-size="12">Position</text>`;
  return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${cssWidth} ${cssHeight}" role="img" aria-label="DNA sequence logo">${body}</svg>`;
}
document.querySelector("#logo-form").addEventListener("submit",async event=>{event.preventDefault();try{const file=fileInput.files[0];if(!file)throw new Error("Choose a FASTA or JSON file.");const data=new FormData(event.currentTarget),sequences=parseSequenceFile(await file.text(),data.get("format"),file.name),options=Object.fromEntries(data);const svg=renderLogo(sequences,options),result=document.querySelector("#logo-result");result.innerHTML=`<div class="result-head"><div><p class="eyebrow">Rendered locally</p><h2>${sequences.length} sequences</h2></div><div class="download-links"><button id="download-svg">Download SVG</button></div></div><div class="logo-wrap">${svg}</div>`;result.querySelector("#download-svg").addEventListener("click",()=>download(svg,"sequence-logo.svg","image/svg+xml"));result.classList.remove("hidden");result.scrollIntoView({behavior:window.matchMedia("(prefers-reduced-motion: reduce)").matches?"auto":"smooth"});}catch(error){handleError(error);}});
