//! Per-term evidence rows (`pgsum evidence`).
//!
//! For every term of a pack, one JSON object: the term's status, the genotype call it was scored from, its
//! effect dosage and contribution, and every gVCF record whose span covers the term's position. A caller that
//! keeps per-term evidence can store these rows instead of reading the gVCF and assessing each term itself. The
//! genotype table must be extracted with `--term-positions`. See DESIGN.md, "Per-term evidence rows".
//!
//! The call is rebuilt from the records (FORMAT/FT and reference-block anchor checks, the documented RefCall
//! convention, then the diploid DP10/GQ20 site rules), and its state must equal pgsum's own call for the term;
//! any difference stops the run.

use std::collections::BTreeSet;
use std::io::Write;

use serde_json::{Map, Value, json};
use sha2::{Digest, Sha256};

use crate::digest::hex;
use crate::genotype::State;
use crate::genotypes::{GenotypeTable, position_key, unpack_key};
use crate::pack::Pack;
use crate::reference::Reference;
use crate::score::{Options, Plan, outcome, plan};
use crate::term::CONTIGS;
use crate::{Error, Result, invalid};

/// One gVCF record: its fields, INFO flags as `true`, sample fields as strings, `END` from INFO or the REF
/// length, whether it is a reference block, and the digest of the line with its newline.
fn record_value(line: &str, refcall_defined: bool, caller_version: Option<&str>) -> Result<Value> {
    let fields: Vec<&str> = line.split('\t').collect();
    let [chrom, position, id, r, alt, qual, filters, info, fmt, sample] = fields[..] else {
        return invalid!("record does not have ten fields");
    };
    let pos: u64 = position
        .parse()
        .map_err(|_| Error::Invalid(format!("invalid POS {position:?}")))?;
    let mut attributes = Map::new();
    for item in info.split(';').filter(|i| *i != ".") {
        match item.split_once('=') {
            Some((k, v)) => attributes.insert(k.to_owned(), Value::String(v.to_owned())),
            None => attributes.insert(item.to_owned(), Value::Bool(true)),
        };
    }
    let end: u64 = match attributes.get("END") {
        Some(Value::String(v)) => v.parse().map_err(|_| Error::Invalid(format!("invalid END {v:?}")))?,
        _ => pos + r.len() as u64 - 1,
    };
    let keys: Vec<&str> = fmt.split(':').collect();
    let values: Vec<&str> = sample.split(':').collect();
    let mut sample_fields = Map::new();
    for (k, v) in keys.iter().zip(&values) {
        sample_fields.insert((*k).to_owned(), Value::String((*v).to_owned()));
    }
    let gt = sample_fields
        .get("GT")
        .and_then(Value::as_str)
        .unwrap_or(".")
        .to_owned();
    let alts: Vec<&str> = alt.split(',').collect();
    let indices: Vec<Value> = gt
        .split(['/', '|'])
        .map(|a| {
            if a == "." {
                Value::Null
            } else {
                a.parse::<u64>().map(Value::from).unwrap_or(Value::Null)
            }
        })
        .collect();
    let mut hash = Sha256::new();
    hash.update(line.as_bytes());
    hash.update(b"\n");
    Ok(json!({
        "reference_call_filter": if refcall_defined { Value::from("RefCall") } else { Value::Null },
        "source_caller": caller_version.map_or(Value::Null, |v| json!({
            "name": "DeepVariant", "version": v, "refcall_is_reference": refcall_defined
        })),
        "chrom": chrom,
        "pos": pos,
        "end": end,
        "id": id,
        "ref": r,
        "alts": alts,
        "qual": qual,
        "filters": filters.split(';').collect::<Vec<_>>(),
        "info": attributes,
        "phase_set": sample_fields.get("PS").cloned().unwrap_or(Value::Null),
        "sample": sample_fields,
        "raw_gt": gt,
        "allele_indices": indices,
        "reference_block": alts.iter().all(|a| matches!(*a, "<NON_REF>" | "<*>" | ".")),
        "query_record_sha256": format!("sha256:{}", hex(&hash.finalize())),
    }))
}

/// A finite, non-negative number, or nothing.
fn number(v: Option<&Value>) -> Option<f64> {
    let x: f64 = v?.as_str()?.trim().parse().ok()?;
    (x.is_finite() && x >= 0.0).then_some(x)
}

/// The call for a resolved SNV term from its records, or the bare state when a record fails its FORMAT/FT or a
/// reference block's anchor base differs from the reference. `policy` is written into the call.
fn call_value(
    records: &[Value],
    contig: &str,
    pos: u64,
    r: &str,
    a: &str,
    reference: &Reference,
    policy: &str,
) -> (String, Value) {
    let early = |state: &str| (state.to_owned(), json!({ "state": state }));
    let mut adapted: Vec<Value> = Vec::with_capacity(records.len());
    let mut adaptations: Vec<&str> = Vec::new();
    for record in records {
        let ft = record["sample"].get("FT").and_then(Value::as_str);
        if ft.is_some_and(|f| f != "PASS" && f != ".") {
            return early("genotype_filtered");
        }
        let block = record["reference_block"].as_bool().unwrap_or(false);
        if block {
            let rp = record["pos"].as_u64().unwrap_or(0);
            let rr = record["ref"].as_str().unwrap_or("");
            if rp == 0
                || reference
                    .fetch(contig, rp - 1, rp - 1 + rr.len() as u64)
                    .ok()
                    .as_deref()
                    != Some(rr)
            {
                return early("reference_anchor_mismatch");
            }
        }
        let mut record = record.clone();
        if record["reference_call_filter"] == "RefCall"
            && record["filters"] == json!(["RefCall"])
            && record["allele_indices"] == json!([0, 0])
        {
            record["filters"] = json!(["PASS"]);
            adaptations.push("documented_RefCall_reference_convention");
        }
        adapted.push(record);
    }
    let mut result = json!({
        "contig": contig, "pos1": pos, "ref": r, "alts": [a], "policy": policy, "state": "unknown_no_record",
        "gt": null, "dosages": null, "phase_set": null, "clinical_use": false,
    });
    let finish = |mut result: Value, state: &str| {
        result["state"] = Value::from(state);
        result["input_adaptations"] = json!(adaptations);
        (state.to_owned(), result)
    };
    let reference_sequence = reference.fetch(contig, pos - 1, pos - 1 + r.len() as u64).ok();
    if reference_sequence.as_deref() != Some(r) {
        return finish(result, "reference_mismatch");
    }
    let row = match &adapted[..] {
        [] => return finish(result, "unknown_no_record"),
        [row] => row,
        _ => return finish(result, "ambiguous_overlapping_records"),
    };
    let indices: Vec<Option<u64>> = row["allele_indices"]
        .as_array()
        .map(|v| v.iter().map(Value::as_u64).collect())
        .unwrap_or_default();
    if indices.len() != 2 {
        return finish(result, "unsupported_ploidy");
    }
    if indices.iter().any(Option::is_none) {
        let state = if indices.iter().any(Option::is_some) {
            "partial_no_call"
        } else {
            "no_call"
        };
        return finish(result, state);
    }
    if row["filters"]
        .as_array()
        .is_some_and(|f| f.iter().any(|v| v != "PASS" && v != "."))
    {
        return finish(result, "filtered");
    }
    if row["sample"]
        .get("FT")
        .and_then(Value::as_str)
        .is_some_and(|f| f != "." && f != "PASS")
    {
        result["filter_field"] = Value::from("FORMAT/FT");
        return finish(result, "filtered");
    }
    let block = row["reference_block"].as_bool().unwrap_or(false);
    let depth_field = if block { "MIN_DP" } else { "DP" };
    let (depth, gq) = (number(row["sample"].get(depth_field)), number(row["sample"].get("GQ")));
    result["quality"] = json!({ "depth_field": depth_field, "depth": depth, "gq": gq });
    let (Some(depth), Some(gq)) = (depth, gq) else {
        return finish(result, "quality_missing");
    };
    if depth < 10.0 || gq < 20.0 {
        return finish(result, "low_quality");
    }
    let indices: Vec<u64> = indices.into_iter().map(|i| i.expect("no no-calls here")).collect();
    let mapped: Vec<u64> = if block {
        if indices.iter().any(|&i| i != 0) {
            return finish(result, "unsupported_symbolic_genotype");
        }
        if row["end"].as_u64().unwrap_or(0) < pos + r.len() as u64 - 1 {
            return finish(result, "incomplete_reference_span");
        }
        vec![0; indices.len()]
    } else {
        if row["pos"].as_u64() != Some(pos) || row["ref"].as_str() != Some(r) {
            return finish(result, "unsupported_allele_representation");
        }
        let mut mapped = Vec::new();
        for &i in &indices {
            if i == 0 {
                mapped.push(0);
            } else if row["alts"][(i - 1) as usize].as_str() == Some(a) {
                mapped.push(1);
            } else {
                return finish(result, "other_called_allele");
            }
        }
        mapped
    };
    let phase_set = row["phase_set"].as_str();
    let phased = indices.len() == 2
        && row["raw_gt"].as_str().is_some_and(|g| g.contains('|'))
        && phase_set.is_some_and(|p| !p.is_empty() && p.chars().all(|c| c.is_ascii_digit()));
    let separator = if phased { "|" } else { "/" };
    result["gt"] = Value::from(mapped.iter().map(u64::to_string).collect::<Vec<_>>().join(separator));
    result["dosages"] = json!([mapped.iter().filter(|&&m| m == 1).count()]);
    result["phase_set"] = if phased {
        Value::from(phase_set.expect("phased"))
    } else {
        Value::Null
    };
    result["phase_preserved"] = Value::from(phased);
    let state = if mapped.iter().any(|&m| m != 0) {
        "observed_variant"
    } else {
        "observed_reference"
    };
    finish(result, state)
}

/// Write one JSON line per term of `pack`: `ordinal`, `status`, `call`, `effect_dosage`, `contribution` and
/// `source_records`, from the term after ordinal `after` (at most `limit` terms); each call carries `policy` as its policy ID. Terms pgsum
/// reads as sequences (indels) are written with pgsum's status and no call, marked `pgsum_call_not_rebuilt`.
pub fn write_rows(
    pack: &Pack,
    table: &GenotypeTable,
    reference: &Reference,
    policy: &str,
    after: usize,
    limit: Option<usize>,
    out: &mut impl Write,
) -> Result<u64> {
    let extracted = table
        .header
        .packs
        .iter()
        .any(|p| p.pgs_id == pack.header.pgs_id && p.records_sha256 == pack.header.records_sha256);
    if !extracted {
        return invalid!(
            "the genotype table was not extracted with this {} pack",
            pack.header.pgs_id
        );
    }
    if !table.header.states.contains_key(State::RecordsAtPosition.as_str()) {
        return invalid!("the genotype table was not extracted with --term-positions");
    }
    let refcall = table.header.sample.refcall_defined;
    let caller = table.header.sample.deepvariant_version.as_deref();
    let options = Options::default();
    let mut n = 0u64;
    for (i, term) in pack.terms_from(after).take(limit.unwrap_or(usize::MAX)).enumerate() {
        let term = term?;
        let ordinal = (after + i) as u64 + 1;
        // The records at the term's position; a record met twice is kept once.
        let mut records: Vec<Value> = Vec::new();
        if term.contig != 0 {
            let mut seen = BTreeSet::new();
            for line in table.records(position_key(term.contig, term.pos))? {
                let record = record_value(line, refcall, caller)?;
                if seen.insert(record.to_string()) {
                    records.push(record);
                }
            }
        }
        let mine = outcome(&term, table, &options)?;
        let mut row = json!({
            "ordinal": ordinal,
            "status": mine.status,
            "call": null,
            "effect_dosage": mine.effect_dosage,
            "contribution": mine.contribution.as_ref().map(|c| c.to_decimal().to_python_string()),
            "source_records": &records,
        });
        match plan(&term, &options).0 {
            Plan::Unscorable(_) => {}
            Plan::Sequence { .. } => row["pgsum_call_not_rebuilt"] = Value::Bool(true),
            Plan::Snv { key, .. } => {
                let (_, _, r, a) = unpack_key(key);
                let (r, a) = ((r as char).to_string(), (a as char).to_string());
                let contig = CONTIGS[term.contig as usize - 1];
                let (state, call) = call_value(&records, contig, term.pos as u64, &r, &a, reference, policy);
                if state != mine.call_state {
                    return invalid!(
                        "term {ordinal} ({contig}:{}): the rebuilt call is {state}, pgsum's call is {}",
                        term.pos,
                        mine.call_state
                    );
                }
                row["call"] = call;
            }
        }
        serde_json::to_writer(&mut *out, &row).map_err(|e| Error::Invalid(e.to_string()))?;
        out.write_all(b"\n").map_err(Error::io("<evidence rows>"))?;
        n += 1;
    }
    Ok(n)
}
