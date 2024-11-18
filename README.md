# SophiaCLI automation scripts

Script designed to automate the workflow for processing Sophia DDM requests. Includes QC metric validation, input file generation, and data upload to Sophia DDM platform.

## IMPORTANT

A regex file is used to determine required metadata for each FASTQ. The naming convention currently approved is:

```
MSK{numbers}_{numbers}_{numbers}[a-z]?_{numbers}[TN]_Pan{numbers}_S{numbers}_R[12]_001.fastq.gz
```

Example valid filenames:
- `MSK24000_01_240000_2440000T_Pan5180_S9_R1_001.fastq.gz`
- `MSK24000_03_240000a_2430000N_Pan5180_S20_R1_001.fastq.gz`

## Prerequisites

- Python 3.x
- Required Python packages:
  - `interop`: For parsing Illumina run metrics
  - Additional dependencies from [Sophia Genetics](https://api-fr.sophiagenetics.com/uploader/cli/docs/):
    - `adegen.py`
    - `sg-upload-v2-wrapper.py`

## Usage

Basic usage:
```bash
python3 sophia.py /path/to/run/folder
```

### Command Line Arguments

- `run_folder`: (Required) Path to the sequencing run folder
- `--skip-upload`: Skip the upload step and only generate JSON
- `--skip-qc`: Skip the QC metrics check
- `--force`: Continue even if QC thresholds are not met
- `--dry-run`: Print commands without executing them

### Examples

Run with all QC checks and upload:
```bash
python3 sophia.py /media/data3/share/241010_A01000_0400_BHM6MHDRX5
```

Generate JSON only:
```bash
python3 pipeline_automation.py /sequencing/runs/240123_run --skip-upload
```

Skip QC checks:
```bash
python3 pipeline_automation.py /sequencing/runs/240123_run --skip-qc
```

Force continue despite QC failures:
```bash
python3 pipeline_automation.py /sequencing/runs/240123_run --force
```

## Workflow Stages

1. **QC Metric Validation**
   - Checks cluster density, percent passing filter, yield, error rate, and Q30 scores
   - Validates against predefined thresholds
   - Can be skipped with `--skip-qc` or forced with `--force`

2. **Input JSON Generation**
   - Runs `adegen.py` to process FASTQ files
   - Automatically handles any required confirmations
   - Creates `{run_name}.inputs.json`

3. **Data Upload**
   - Executes `sg-upload-v2-wrapper.py` to upload data to Sophia DDM
   - Can be skipped with `--skip-upload`

## QC Thresholds

Current thresholds for QC metrics:
- Clusters Passing Filter: ≥ 70%
- Error Rate: ≤ 2%
- Q30 Score: ≥ 80%

## Logging

The script generates detailed logs including:
- Console output for real-time monitoring
- Log files stored in `logs` directory
- Each run creates a separate log file: `logs/{run_name}.log`
- Includes timestamps, log levels, and detailed operation information
- Captures both stdout and stderr from subprocesses

## Error Handling

The script includes comprehensive error handling for:
- Missing or invalid run folders
- QC threshold failures
- Missing required files
- Process timeouts
- Upload failures

All errors are logged to both console and log file.

## Troubleshooting

Common issues and solutions:

1. **QC Threshold Failures**
   - Review logs for specific metrics that failed
   - Use `--force` if metrics are acceptable despite being out of range

2. **Missing Files**
   - Ensure RunParameters.xml exists in run folder
   - Verify FASTQ files are in correct BaseCalls directory
   - Check that all required scripts are in working directory

3. **Upload Failures**
   - Verify network connectivity
   - Check Sophia DDM credentials
   - Review upload wrapper logs for specific errors

## Support

For issues with:
- Run metrics: Check Illumina documentation
- Upload process: Contact Sophia DDM support