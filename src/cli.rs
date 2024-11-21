use clap::Parser;

#[derive(Parser, Debug)]
#[clap(version, about="Reads a BAM and a BED, outputs reads from input only if they start and end in one of the specified BED intervals", long_about = None)]
pub struct Cli {
    /// Input BAM file
    #[arg(short, long)]
    input: String,

    /// BED file
    #[arg(short, long)]
    bed: String,

    /// Output on-target reads
    #[arg(long)]
    on_target: Option<String>,

    /// Output off-target reads
    #[arg(long)]
    off_target: Option<String>,

    /// Output reads that start in a BED interval but do not end in one
    #[arg(long)]
    start_no_end: Option<String>,

    /// Output reads that end in a BED interval but do not start in one
    #[arg(long)]
    end_no_start: Option<String>,

    /// Output on-targret only
    #[arg(long)]
    on_target_only: bool,

    /// Tolerance (default: 0)
    /// Number of bases that the read can extend past the BED interval
    #[arg(long, default_value = "0")]
    tolerance: u32,
}

impl Cli {
    pub fn input(&self) -> &str {
        &self.input
    }

    pub fn bed(&self) -> &str {
        &self.bed
    }

    pub fn on_target(&self) -> Option<&str> {
        self.on_target.as_deref()
    }

    pub fn off_target(&self) -> Option<&str> {
        self.off_target.as_deref()
    }

    pub fn tolerance(&self) -> u32 {
        self.tolerance
    }

    pub fn start_no_end(&self) -> Option<&str> {
        self.start_no_end.as_deref()
    }

    pub fn end_no_start(&self) -> Option<&str> {
        self.end_no_start.as_deref()
    }

    pub fn on_target_only(&self) -> bool {
        self.on_target_only
    }
}
