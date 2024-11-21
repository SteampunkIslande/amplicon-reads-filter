use clap::Parser;
use noodles::sam::alignment::io::Write;
use noodles::sam::alignment::Record;

use std::collections::HashMap;
use std::collections::HashSet;

pub mod cli;
pub mod index;

use noodles::bam;
use noodles::bed;

fn start_spot_on(target: (usize, usize, usize), read: &bam::Record, tolerance: usize) -> bool {
    let (contig_id, start, _) = target;

    match read.reference_sequence_id() {
        Some(Ok(read_contig_id)) if read_contig_id == contig_id => match read.alignment_start() {
            Some(Ok(alignment_start)) => alignment_start.get().abs_diff(start) <= tolerance,
            Some(Err(_)) => false,
            None => false,
        },
        _ => false,
    }
}

fn end_spot_on(target: (usize, usize, usize), read: &bam::Record, tolerance: usize) -> bool {
    let (contig_id, _, end) = target;

    match read.reference_sequence_id() {
        Some(Ok(read_contig_id)) if read_contig_id == contig_id => match read.alignment_end() {
            Some(Ok(alignment_end)) => alignment_end.get().abs_diff(end) <= tolerance,
            Some(Err(_)) => false,
            None => false,
        },
        _ => false,
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = cli::Cli::parse();

    let mut bam_reader = bam::io::reader::Builder::default().build_from_path(args.input())?;
    let mut bed_reader: bed::io::reader::Reader<3, _> =
        bed::io::reader::Builder::default().build_from_path(args.bed())?;

    let on_target_fn = args.on_target().map(|s| s.to_owned()).or(Some(format!(
        "{}.on-target.bam",
        args.input().replace(".bam", "")
    )));
    let mut on_target_writer: Option<bam::io::writer::Writer<_>> = Some(
        bam::io::writer::Builder::default()
            .build_from_path(
                &on_target_fn
                    .as_ref()
                    .expect("Cannot create on-target output file"),
            )
            .expect("Cannot create on-target output file"),
    );

    let off_target_fn = args.off_target().map(|s| s.to_owned()).or(Some(format!(
        "{}.off-target.bam",
        args.input().replace(".bam", "")
    )));
    let mut off_target_writer = if !args.on_target_only() {
        Some(
            bam::io::writer::Builder::default()
                .build_from_path(
                    &off_target_fn
                        .as_ref()
                        .expect("Cannot create off-target output file"),
                )
                .expect("Cannot create off-target output file"),
        )
    } else {
        None
    };

    let start_no_end_fn = args.start_no_end().map(|s| s.to_owned()).or(Some(format!(
        "{}.start-no-end.bam",
        args.input().replace(".bam", "")
    )));
    let mut start_no_end_writer = if !args.on_target_only() {
        Some(
            bam::io::writer::Builder::default()
                .build_from_path(
                    &start_no_end_fn
                        .as_ref()
                        .expect("Cannot create start-no-end output file"),
                )
                .expect("Cannot create start-no-end output file"),
        )
    } else {
        None
    };

    let end_no_start_fn = args.end_no_start().map(|s| s.to_owned()).or(Some(format!(
        "{}.end-no-start.bam",
        args.input().replace(".bam", "")
    )));
    let mut end_no_start_writer = if !args.on_target_only() {
        Some(
            bam::io::writer::Builder::default()
                .build_from_path(
                    &end_no_start_fn
                        .as_ref()
                        .expect("Cannot create end-no-start output file"),
                )
                .expect("Cannot create end-no-start output file"),
        )
    } else {
        None
    };

    let mut bed_intervals: HashSet<(usize, usize, usize)> = HashSet::new();
    let mut bed_record = bed::Record::<3>::default();

    let header = bam_reader.read_header()?;

    let contig_hashmap = header
        .clone()
        .reference_sequences()
        .iter()
        .enumerate()
        .map(|(i, (s, _))| (s.to_owned(), i))
        .collect::<HashMap<bstr::BString, usize>>();

    while bed_reader
        .read_record(&mut bed_record)
        .expect("Could not read BED record")
        > 0
    {
        if let Some(contig_id) = contig_hashmap.get(bed_record.reference_sequence_name()) {
            bed_intervals.insert((
                *contig_id,
                bed_record.feature_start().unwrap().get(),
                bed_record.feature_end().unwrap().unwrap().get(),
            ));
        }
    }

    let bed_intervals = bed_intervals;

    let tolerance = args.tolerance() as usize;

    on_target_writer.as_mut().map(|f| f.write_header(&header));
    off_target_writer.as_mut().map(|f| f.write_header(&header));
    start_no_end_writer
        .as_mut()
        .map(|f| f.write_header(&header));
    end_no_start_writer
        .as_mut()
        .map(|f| f.write_header(&header));
    bam_reader
        .records()
        .filter_map(|r| r.ok())
        .for_each(|record: bam::Record| {
            let start_spot_on = bed_intervals
                .iter()
                .any(|target| start_spot_on(*target, &record, tolerance));
            let end_spot_on = bed_intervals
                .iter()
                .any(|target| end_spot_on(*target, &record, tolerance));
            match (start_spot_on, end_spot_on) {
                (true, true) => {
                    if let Some(writer) = on_target_writer.as_mut() {
                        writer
                            .write_alignment_record(&header, &record)
                            .expect("Could not write on-target record");
                    }
                }
                (false, false) => {
                    if let Some(writer) = off_target_writer.as_mut() {
                        writer
                            .write_alignment_record(&header, &record)
                            .expect("Could not write off-target record");
                    }
                }
                (true, false) => {
                    if let Some(writer) = start_no_end_writer.as_mut() {
                        writer
                            .write_alignment_record(&header, &record)
                            .expect("Could not write incomplete-end record");
                    }
                }
                (false, true) => {
                    if let Some(writer) = end_no_start_writer.as_mut() {
                        writer
                            .write_alignment_record(&header, &record)
                            .expect("Could not write incomplete-start record");
                    }
                }
            }
        });

    // Close all writers

    on_target_writer.as_mut().map(|f| f.finish(&header));
    off_target_writer.as_mut().map(|f| f.finish(&header));
    start_no_end_writer.as_mut().map(|f| f.finish(&header));
    end_no_start_writer.as_mut().map(|f| f.finish(&header));

    // Create index for every output BAM
    if let Some(ref on_target_fn) = on_target_fn {
        let index_fn = format!("{}.bai", on_target_fn);
        let index_writer = std::fs::File::create(index_fn).expect("Cannot create index file");
        let on_target_reader =
            std::fs::File::open(&on_target_fn).expect("Cannot open on-target file");
        index::index(on_target_reader, index_writer).expect("Cannot create index");
    }

    if let Some(ref off_target_fn) = off_target_fn {
        let index_fn = format!("{}.bai", off_target_fn);
        let index_writer = std::fs::File::create(index_fn).expect("Cannot create index file");
        let off_target_reader =
            std::fs::File::open(&off_target_fn).expect("Cannot open off-target file");
        index::index(off_target_reader, index_writer).expect("Cannot create index");
    }

    if let Some(ref start_no_end_fn) = start_no_end_fn {
        let index_fn = format!("{}.bai", start_no_end_fn);
        let index_writer = std::fs::File::create(index_fn).expect("Cannot create index file");
        let start_no_end_reader =
            std::fs::File::open(&start_no_end_fn).expect("Cannot open start-no-end file");
        index::index(start_no_end_reader, index_writer).expect("Cannot create index");
    }

    if let Some(ref end_no_start_fn) = end_no_start_fn {
        let index_fn = format!("{}.bai", end_no_start_fn);
        let index_writer = std::fs::File::create(index_fn).expect("Cannot create index file");
        let end_no_start_reader =
            std::fs::File::open(&end_no_start_fn).expect("Cannot open end-no-start file");
        index::index(end_no_start_reader, index_writer).expect("Cannot create index");
    }

    Ok(())
}
