use noodles::bam;
use noodles::bam::bai;

use std::io;
use std::io::{Read, Write};

use noodles::core::Position;
use noodles::csi::binning_index::{index::reference_sequence::bin::Chunk, Indexer};
use noodles::sam;
use noodles::sam::alignment::Record;

fn is_coordinate_sorted(header: &sam::Header) -> bool {
    use sam::header::record::value::map::header::{sort_order, tag};

    header
        .header()
        .and_then(|hdr| hdr.other_fields().get(&tag::SORT_ORDER))
        .map(|sort_order| sort_order == sort_order::COORDINATE)
        .unwrap_or_default()
}

fn alignment_context(
    record: &bam::Record,
) -> io::Result<(Option<usize>, Option<Position>, Option<Position>)> {
    Ok((
        record.reference_sequence_id().transpose()?,
        record.alignment_start().transpose()?,
        record.alignment_end().transpose()?,
    ))
}

pub fn index(input_bam_stream: impl Read, output_index_stream: impl Write) -> io::Result<()> {
    let mut reader = bam::io::Reader::new(input_bam_stream);
    let header = reader.read_header()?;

    if !is_coordinate_sorted(&header) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "the input BAM must be coordinate-sorted to be indexed",
        ));
    }

    let mut record = bam::Record::default();

    let mut builder = Indexer::default();
    let mut start_position = reader.get_ref().virtual_position();

    while reader.read_record(&mut record)? != 0 {
        let end_position = reader.get_ref().virtual_position();
        let chunk = Chunk::new(start_position, end_position);

        let alignment_context = match alignment_context(&record)? {
            (Some(id), Some(start), Some(end)) => {
                let is_mapped = !record.flags().is_unmapped();
                Some((id, start, end, is_mapped))
            }
            _ => None,
        };

        builder.add_record(alignment_context, chunk)?;

        start_position = end_position;
    }

    let index = builder.build(header.reference_sequences().len());
    println!("Successfully built index");
    let mut writer = bai::io::Writer::new(output_index_stream);
    println!("Writing index");
    writer.write_index(&index)?;

    Ok(())
}
