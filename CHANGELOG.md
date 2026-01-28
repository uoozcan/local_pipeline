# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] - 2025-11-21

### Added
- Complete multi-tool HLA typing pipeline with Nextflow DSL2
- Support for three HLA typing tools:
  - OptiType (DNA/RNA-seq)
  - ArcasHLA (RNA-seq optimized, works with DNA)
  - SpecHLA (Exome/WGS with variant calling)
- Majority voting / consensus calling module
- Flexible input support (BAM, CRAM, FASTQ)
- HLA region extraction for space-efficient processing
- Chromosome 6 extraction optimization for SpecHLA
- CSC Puhti HPC profile with SLURM support
- Docker and Singularity container support
- Comprehensive error handling and logging
- Automatic BAM index creation
- Chromosome naming convention detection (chr6 vs 6)
- Exon-only mode for SpecHLA (critical for exome data)
- SLURM submission script with proper parameter handling
- Comprehensive documentation (README, QUICKSTART guide)
- Pipeline execution reports (timeline, trace, DAG)

### Features
- **Smart Input Handling**: Automatically detects and processes different input formats
- **Optimized for Exome Data**: Special exon-only mode for SpecHLA
- **Space Efficient**: Optional HLA region extraction to reduce disk usage
- **Resume Capability**: Nextflow's built-in resume functionality
- **Parallel Processing**: Efficiently processes multiple samples simultaneously
- **Comprehensive Logging**: Detailed logs for each tool and process
- **Resource Management**: Configurable CPU, memory, and time limits
- **Container Ready**: Pre-configured for Docker and Singularity

### Fixed
- SLURM script parameter handling (parameters now properly passed to Nextflow)
- Container paths for CSC Puhti pre-downloaded SIF files
- SpecHLA chromosome naming issues (chr6 vs 6)
- BAM indexing in various scenarios
- Error handling for missing or empty files

### Configuration
- CSC Puhti profile with pre-downloaded containers
- SLURM executor configuration
- Resource labels for different process types
- Automatic retry on failure
- Comprehensive parameter validation

### Documentation
- Complete README with usage examples
- Quick Start guide for new users
- Detailed parameter descriptions
- Troubleshooting section
- CSC Puhti specific instructions
- Performance benchmarks

### Known Issues
- SpecHLA requires exon-only mode (`--spechla_exon_only 1`) for exome data
- Some tools may produce empty results for low-coverage samples
- Chromosome 6 naming must be detected correctly from BAM headers

### Compatibility
- Nextflow: ≥23.04.0
- Docker/Singularity required
- Tested on CSC Puhti (SLURM)
- Compatible with hg19 and hg38 reference builds

## [1.0.0] - 2025-11-XX

### Initial Development
- Basic pipeline structure
- Single tool support
- Initial testing on 1000 Genomes data

---

## Planned for Future Releases

### [2.1.0] - Planned
- [ ] Additional HLA typing tools (HLA-HD, HLA*LA, xHLA)
- [ ] Advanced consensus algorithms
- [ ] Quality score integration
- [ ] Phasing information
- [ ] Interactive HTML reports
- [ ] Database storage support

### [2.2.0] - Planned
- [ ] Cloud deployment (AWS, Google Cloud)
- [ ] GUI interface
- [ ] Real-time monitoring dashboard
- [ ] Automated quality control
- [ ] Population frequency analysis

### [3.0.0] - Planned
- [ ] Machine learning-based consensus
- [ ] Novel allele detection
- [ ] Integration with clinical databases
- [ ] Multi-sample joint calling
