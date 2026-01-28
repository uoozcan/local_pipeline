# Project Structure

```
hla-typing-pipeline/
│
├── main.nf                      # Main Nextflow pipeline script
├── nextflow.config              # Global configuration
├── README.md                    # Main documentation
├── LICENSE                      # MIT License
├── CHANGELOG.md                 # Version history
├── QUICKSTART.sh               # Quick start helper script
│
├── conf/                        # Configuration files
│   ├── params.config           # Parameter definitions
│   └── puhti.config            # CSC Puhti specific config
│
├── modules/                     # Nextflow modules
│   ├── bam_to_fastq.nf        # BAM conversion
│   ├── optitype.nf            # OptiType integration
│   ├── arcashla.nf            # ArcasHLA integration
│   ├── spechla.nf             # SpecHLA integration
│   ├── xhla.nf                # xHLA integration
│   ├── hlahd.nf               # HLA-HD integration
│   ├── hlala.nf               # HLA*LA integration
│   ├── seq2hla.nf             # Seq2HLA integration
│   ├── hlaprofiler.nf         # HLAProfiler integration
│   ├── kourami.nf             # Kourami integration
│   └── majority_voting.nf     # Consensus calling
│
├── bin/                         # Scripts and utilities
│   └── majority_voting.py      # Consensus algorithm
│
├── docs/                        # Documentation
│   ├── INSTALLATION.md         # Installation guide (to be created)
│   └── USAGE.md               # Usage examples
│
├── .github/                     # GitHub specific
│   └── workflows/
│       └── ci.yml              # CI/CD configuration
│
└── .gitignore                   # Git ignore patterns
```

## Key Files

### Pipeline Core
- **main.nf**: Orchestrates the entire workflow
- **nextflow.config**: Defines execution profiles and resource allocation
- **conf/params.config**: All configurable parameters

### Modules
Each tool has its own module with:
- Input/output specifications
- Container definitions
- Resource requirements
- Error handling

### Scripts
- **majority_voting.py**: Implements consensus calling algorithm
- **QUICKSTART.sh**: Interactive setup helper

### Documentation
- **README.md**: Overview and quick reference
- **docs/USAGE.md**: Comprehensive usage examples
- **CHANGELOG.md**: Version history and changes

## File Locations

### Input Data
Place your data in:
- `samples/` - For input files
- `references/` - For reference genomes (not tracked in git)

### Output Data
Results will be generated in:
- `results/` - Default output directory (configurable)
- `work/` - Nextflow temporary files (auto-generated)

### Container Cache
- `singularity_cache/` - Downloaded containers (not tracked in git)

## Configuration Files

### Execution Profiles
- `standard`: Local execution
- `docker`: Docker containers
- `singularity`: Singularity containers
- `slurm`: SLURM scheduler
- `puhti`: CSC Puhti supercomputer

### Custom Configuration
Create custom configs in `conf/` and include them:
```bash
nextflow run main.nf -c conf/my_custom.config
```
