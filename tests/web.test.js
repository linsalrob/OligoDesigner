import test from "node:test";
import assert from "node:assert/strict";
import { generateOligos, generateStructuredOligos, logoMatrix, meltingTemperature, parseSequenceFile, reverseComplement, toFasta, toJson, toTsv } from "../web/core.js";

const defaults={count:3,length:12,seed:42,prefix:"oligo",fivePrime:"",threePrime:"",fiveRandomLength:null,threeRandomLength:null,sameRandom:false,deduplicate:false,minStem:4,minLoop:3,maxLoop:8,minHpRun:4,minOverlap:10};
const structured={...defaults,count:2,prefix:"soligo",type:"all",halfLength:6,outerArmLength:8,innerHalfLength:6,spacerLength:2};

test("reverse complement and Tm match core DNA semantics",()=>{assert.equal(reverseComplement("GAATTC"),"GAATTC");assert.equal(reverseComplement("AACG"),"CGTT");assert.ok(Number.isFinite(meltingTemperature("GCATGCATGCAT")));});
test("random generation is reproducible and analysed",()=>{const a=generateOligos(defaults),b=generateOligos(defaults);assert.deepEqual(a,b);assert.equal(a.length,3);assert.equal(a[0].sequence.length,12);assert.match(a[0].sequence,/^[ACGT]+$/);assert.equal(typeof a[0].has_hairpin,"boolean");});
test("fixed and random flanks are supported",()=>{const fixed=generateOligos({...defaults,fivePrime:"AAA",threePrime:"TT"});assert.ok(fixed.every(item=>item.sequence.startsWith("AAA")&&item.sequence.endsWith("TT")));const random=generateOligos({...defaults,fiveRandomLength:3});assert.equal(new Set(random.map(item=>item.sequence.slice(0,3))).size,3);});
test("invalid mutually exclusive flank options are rejected",()=>assert.throws(()=>generateOligos({...defaults,fivePrime:"AAA",fiveRandomLength:3}),/mutually exclusive/));
test("structured all creates every structure type",()=>{const items=generateStructuredOligos(structured);assert.equal(items.length,6);assert.deepEqual(new Set(items.map(item=>item.oligo_type)),new Set(["palindromic_motif","inverted_repeat","at_rich_palindrome"]));assert.ok(items.every(item=>item.is_palindrome));});
test("structured spacer validates supported values",()=>assert.throws(()=>generateStructuredOligos({...structured,spacerLength:1}),/Spacer length/));
test("all download serializations include every entry",()=>{const entries=generateOligos(defaults);assert.equal((toFasta(entries).match(/^>/gm)||[]).length,3);assert.equal(JSON.parse(toJson(entries)).length,3);assert.equal(toTsv(entries).trim().split("\n").length,4);});
test("FASTA parser handles multiline input and auto detection",()=>assert.deepEqual(parseSequenceFile(">a\nACGT\nAC\n>b\nTTTT\n","auto","input.fa"),["ACGTAC","TTTT"]));
test("JSON parser reads generated output",()=>{const entries=generateOligos(defaults);assert.deepEqual(parseSequenceFile(toJson(entries),"json"),entries.map(item=>item.sequence));});
test("logo matrices support counts, probability, and information",()=>{const seqs=["AA","AC","AG","AT"];assert.deepEqual(logoMatrix(seqs,"counts")[1],{A:1,C:1,G:1,T:1});assert.deepEqual(logoMatrix(seqs,"probability")[1],{A:.25,C:.25,G:.25,T:.25});assert.deepEqual(logoMatrix(seqs,"information")[1],{A:0,C:0,G:0,T:0});});
test("logo truncates unequal sequences to shortest",()=>assert.equal(logoMatrix(["ACGT","AC"],"counts").length,2));
