const BASES = "ACGT";
const COMPLEMENT = { A: "T", C: "G", G: "C", T: "A", N: "N" };

export function seededRandom(seed) {
  if (seed === "" || seed === null || seed === undefined) return Math.random;
  let state = Number(seed) >>> 0;
  return () => {
    state |= 0;
    state = (state + 0x6d2b79f5) | 0;
    let value = Math.imul(state ^ (state >>> 15), 1 | state);
    value = (value + Math.imul(value ^ (value >>> 7), 61 | value)) ^ value;
    return ((value ^ (value >>> 14)) >>> 0) / 4294967296;
  };
}

const randomSequence = (length, random, alphabet = BASES) =>
  Array.from({ length }, () => alphabet[Math.floor(random() * alphabet.length)]).join("");

export const reverseComplement = sequence =>
  [...sequence.toUpperCase()].reverse().map(base => COMPLEMENT[base] ?? base).join("");

export function entropy(sequence) {
  if (!sequence.length) return 0;
  return [...BASES].reduce((total, base) => {
    const count = [...sequence].filter(value => value === base).length;
    if (!count) return total;
    const p = count / sequence.length;
    return total - p * Math.log2(p);
  }, 0);
}

const longestHomopolymer = sequence => {
  let longest = 0, run = 0, previous = "";
  for (const base of sequence) {
    run = base === previous ? run + 1 : 1;
    longest = Math.max(longest, run);
    previous = base;
  }
  return longest;
};

const hasHairpin = (sequence, minStem, minLoop, maxLoop) => {
  for (let stem = minStem; stem <= Math.floor(sequence.length / 2); stem++) {
    for (let loop = minLoop; loop <= maxLoop; loop++) {
      const required = 2 * stem + loop;
      if (required > sequence.length) break;
      for (let i = 0; i <= sequence.length - required; i++) {
        const left = sequence.slice(i, i + stem);
        const right = sequence.slice(i + stem + loop, i + required);
        if (/^[ACGT]+$/.test(left + right) && right === reverseComplement(left)) return true;
      }
    }
  }
  return false;
};

const hasTandemRepeat = sequence => {
  for (let size = 2; size <= 4; size++) {
    for (let i = 0; i <= sequence.length - size * 3; i++) {
      const unit = sequence.slice(i, i + size);
      if (sequence.slice(i, i + size * 3) === unit.repeat(3)) return true;
    }
  }
  return false;
};

const isLowComplexity = sequence => {
  if (!sequence.length) return false;
  const width = Math.min(10, sequence.length);
  for (let i = 0; i <= sequence.length - width; i++) {
    const window = sequence.slice(i, i + width);
    if (Math.max(...[...BASES].map(base => [...window].filter(v => v === base).length)) / width >= 0.7) return true;
  }
  return false;
};

const NN = {AA:[-7.9,-22.2],TT:[-7.9,-22.2],AT:[-7.2,-20.4],TA:[-7.2,-21.3],CA:[-8.5,-22.7],TG:[-8.5,-22.7],GT:[-8.4,-22.4],AC:[-8.4,-22.4],CT:[-7.8,-21],AG:[-7.8,-21],GA:[-8.2,-22.2],TC:[-8.2,-22.2],CG:[-10.6,-27.2],GC:[-9.8,-24.4],GG:[-8,-19.9],CC:[-8,-19.9]};
export function meltingTemperature(sequence) {
  const clean = sequence.replace(/[^ACGT]/g, "");
  if (clean.length < 2) return null;
  let dh = 0, ds = 0;
  for (let i = 0; i < clean.length - 1; i++) { const pair = NN[clean.slice(i, i + 2)]; dh += pair[0]; ds += pair[1]; }
  for (const end of [clean[0], clean.at(-1)]) { const init = "GC".includes(end) ? [0.1,-2.8] : [2.3,4.1]; dh += init[0]; ds += init[1]; }
  const self = clean === reverseComplement(clean);
  if (self) ds -= 1.4;
  const concentration = self ? 250e-9 : 250e-9 / 4;
  return dh * 1000 / (ds + 1.987 * Math.log(concentration)) - 273.15 + 16.6 * Math.log10(0.05);
}

function analyse(sequence, name, options) {
  const clean = sequence.replace(/N/g, "");
  const composition = Object.fromEntries([...BASES].map(base => [base, [...sequence].filter(v => v === base).length]));
  const longest = longestHomopolymer(clean);
  return { name, sequence, length: sequence.length, gc_content: clean.length ? (composition.G + composition.C) / clean.length : 0,
    entropy: entropy(clean), tm: meltingTemperature(clean), base_composition: composition, longest_homopolymer: longest,
    has_homopolymer: longest >= options.minHpRun, is_low_complexity: isLowComplexity(clean),
    is_palindrome: sequence === reverseComplement(sequence), has_hairpin: hasHairpin(sequence, options.minStem, options.minLoop, options.maxLoop),
    has_tandem_repeat: hasTandemRepeat(clean), complementary_to: [] };
}

function addComplementarity(entries, minOverlap) {
  for (let i = 0; i < entries.length; i++) for (let j = i + 1; j < entries.length; j++) {
    const a = entries[i].sequence.replace(/N/g, ""), rc = reverseComplement(entries[j].sequence.replace(/N/g, ""));
    let matches = false;
    for (let k = 0; k <= a.length - minOverlap && !matches; k++) matches = rc.includes(a.slice(k, k + minOverlap));
    if (matches) { entries[i].complementary_to.push(entries[j].name); entries[j].complementary_to.push(entries[i].name); }
  }
}

function flankGenerator(fixed, length, same, count, random) {
  if (fixed) return () => fixed.toUpperCase();
  if (!length) return () => "";
  if (!same && count > 4 ** length) throw new Error(`Cannot generate ${count} unique flanks of length ${length}.`);
  const shared = same ? randomSequence(length, random) : null, seen = new Set();
  return () => { if (shared !== null) return shared; let value; do value = randomSequence(length, random); while (seen.has(value)); seen.add(value); return value; };
}

function validate(options, structured = false) {
  for (const [label, value] of [["count", options.count], ["minimum stem", options.minStem], ["minimum loop", options.minLoop], ["homopolymer run", options.minHpRun], ["minimum overlap", options.minOverlap]]) if (!Number.isInteger(value) || value < 1) throw new Error(`${label} must be at least 1.`);
  if (options.maxLoop < options.minLoop) throw new Error("Maximum loop must be at least minimum loop.");
  for (const [label, fixed, length] of [["5′",options.fivePrime,options.fiveRandomLength],["3′",options.threePrime,options.threeRandomLength]]) {
    if (fixed && length) throw new Error(`${label} fixed and random flanks are mutually exclusive.`);
    if (fixed && !/^[ACGT]+$/i.test(fixed)) throw new Error(`${label} flank must contain only A, C, G and T.`);
  }
  if (options.sameRandom && !options.fiveRandomLength && !options.threeRandomLength) throw new Error("Shared random flank requires a random flank length.");
  if (structured && ![0,2,3,4,5,6].includes(options.spacerLength)) throw new Error("Spacer length must be 0 or 2–6.");
}

export function generateOligos(options) {
  validate(options); if (options.length < 1) throw new Error("Length must be at least 1.");
  const random = seededRandom(options.seed), five = flankGenerator(options.fivePrime, options.fiveRandomLength, options.sameRandom, options.count, random), three = flankGenerator(options.threePrime, options.threeRandomLength, options.sameRandom, options.count, random);
  const width = String(options.count).length;
  let entries = Array.from({length: options.count}, (_, i) => analyse(five() + randomSequence(options.length, random) + three(), `${options.prefix}${String(i + 1).padStart(width,"0")}`, options));
  if (options.deduplicate || options.sameRandom) entries = entries.filter((item, i, all) => all.findIndex(other => other.sequence === item.sequence) === i);
  addComplementarity(entries, options.minOverlap); return entries;
}

export function generateStructuredOligos(options) {
  validate(options, true);
  for (const [label,value] of [["Half length",options.halfLength],["Outer arm length",options.outerArmLength],["Inner half length",options.innerHalfLength]]) if (value < 1) throw new Error(`${label} must be at least 1.`);
  const random = seededRandom(options.seed), types = options.type === "all" ? ["palindrome","inverted_repeat","at_rich"] : [options.type], total = options.count * types.length;
  const five = flankGenerator(options.fivePrime, options.fiveRandomLength, options.sameRandom, total, random), three = flankGenerator(options.threePrime, options.threeRandomLength, options.sameRandom, total, random), width = String(total).length;
  let index = 0, entries = [];
  for (const type of types) for (let n = 0; n < options.count; n++) {
    let left, right, innerLeft = "", innerRight = "", spacer, sequence, oligoType;
    if (type === "inverted_repeat") { left = randomSequence(options.outerArmLength, random); right = reverseComplement(left); innerLeft = randomSequence(options.innerHalfLength, random); innerRight = reverseComplement(innerLeft); spacer = randomSequence(options.spacerLength, random); sequence = left + innerLeft + spacer + innerRight + right; oligoType = type; }
    else { left = randomSequence(options.halfLength, random, type === "at_rich" ? "AT" : BASES); right = reverseComplement(left); spacer = type === "at_rich" ? "N".repeat(options.spacerLength) : randomSequence(options.spacerLength, random); sequence = left + spacer + right; oligoType = type === "at_rich" ? "at_rich_palindrome" : "palindromic_motif"; }
    sequence = five() + sequence + three(); const item = analyse(sequence, `${options.prefix}${String(++index).padStart(width,"0")}`, options);
    Object.assign(item, {oligo_type:oligoType,left_arm:left,right_arm:right,spacer,inner_left:innerLeft,inner_right:innerRight,is_palindrome:right===reverseComplement(left),inner_is_palindrome:innerLeft ? innerRight===reverseComplement(innerLeft) : false,min_stem:options.minStem,min_loop:options.minLoop,max_loop:options.maxLoop,min_hp_run:options.minHpRun}); entries.push(item);
  }
  if (options.deduplicate || options.sameRandom) entries = entries.filter((item, i, all) => all.findIndex(other => other.sequence === item.sequence) === i);
  addComplementarity(entries, options.minOverlap); return entries;
}

export const toFasta = entries => entries.map(item => `>${item.name}\n${item.sequence}\n`).join("");
export const toJson = entries => JSON.stringify(entries, null, 2) + "\n";
export function toTsv(entries) {
  if (!entries.length) return "";
  const structured = "oligo_type" in entries[0];
  const headers = structured ? ["name","sequence","length","oligo_type","left_arm","right_arm","spacer","inner_left","inner_right","is_palindrome","inner_is_palindrome","gc_content","entropy","tm","has_hairpin","has_homopolymer","has_tandem_repeat","complementary_to"] : ["name","sequence","length","gc_content","entropy","tm","count_A","count_C","count_G","count_T","longest_homopolymer","has_homopolymer","is_low_complexity","is_palindrome","has_hairpin","has_tandem_repeat","complementary_to"];
  const value = (item, key) => key.startsWith("count_") ? item.base_composition[key.at(-1)] : key === "complementary_to" ? item[key].join(",") : ["gc_content","entropy"].includes(key) ? item[key].toFixed(4) : key === "tm" ? (item[key] == null ? "" : item[key].toFixed(2)) : item[key];
  return [headers.join("\t"), ...entries.map(item => headers.map(key => value(item,key)).join("\t"))].join("\n") + "\n";
}

export function parseSequenceFile(text, format = "auto", filename = "") {
  const kind = format === "auto" ? (/\.(fa|fasta|fna|fas)$/i.test(filename) || text.trimStart().startsWith(">") ? "fasta" : "json") : format;
  let sequences;
  if (kind === "json") { const data = JSON.parse(text); if (!Array.isArray(data)) throw new Error("JSON must contain an array of oligos."); sequences = data.map(item => item.sequence); }
  else { sequences = []; let current = ""; for (const raw of text.split(/\r?\n/)) { const line = raw.trim(); if (line.startsWith(">")) { if (current) sequences.push(current); current = ""; } else if (line) current += line.toUpperCase(); } if (current) sequences.push(current); }
  if (!sequences.length) throw new Error("No sequences were found.");
  if (sequences.some(sequence => typeof sequence !== "string" || !/^[ACGTN]+$/i.test(sequence))) throw new Error("Sequences must contain only DNA bases (ACGTN).");
  return sequences.map(sequence => sequence.toUpperCase());
}

export function logoMatrix(sequences, type = "counts") {
  const width = Math.min(...sequences.map(value => value.length));
  return Array.from({length: width}, (_, position) => {
    const counts = Object.fromEntries([...BASES].map(base => [base, sequences.filter(seq => seq[position] === base).length]));
    if (type === "counts") return counts;
    const observed = Object.values(counts).reduce((total, count) => total + count, 0);
    const probability = Object.fromEntries([...BASES].map(base => [base, observed ? counts[base] / observed : 0]));
    if (type === "probability") return probability;
    const h = -Object.values(probability).reduce((sum,p) => sum + (p ? p * Math.log2(p) : 0), 0), information = 2 - h;
    return Object.fromEntries([...BASES].map(base => [base, probability[base] * information]));
  });
}

export function logoStackBases(column, stackOrder) {
  if (stackOrder === "alphabetical") return ["T", "G", "C", "A"];
  return ["A", "C", "G", "T"].sort((a, b) => column[b] - column[a]);
}
