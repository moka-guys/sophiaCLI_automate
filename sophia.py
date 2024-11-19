"""
AIO script for automating Sophia CLI processing.

1. Validates sequencing QC metrics (cluster density, PF clusters, yield, error rate, Q30)
2. Generates input JSON files using adegen.py
3. Uploads data to Sophia DDM using sg-upload-v2-wrapper.py

Usage:
    python3 sophia.py /media/data3/share/run_folder [options]

Arguments:
    run_folder (str): Path to the Illumina run folder

Options:
    --skip-upload: Skip the upload step and only generate JSON
    --skip-qc: Skip the QC metrics check
    --force: Continue even if QC thresholds are not met
    --dry-run: Print commands without executing them

QC Thresholds:
    - Clusters Passing Filter: ≥ 70%
    - Error Rate: ≤ 2%
    - Q30 Score: ≥ 80%

File Naming Pattern:
    MSK{numbers}_{numbers}_{numbers}[a-z]?_{numbers}[TN]_Pan{numbers}_S{numbers}_R[12]_001.fastq.gz
"""

#!/usr/bin/env python3

import argparse
import time
import os
import xml.etree.ElementTree as ET
import subprocess
import sys
import logging
from pathlib import Path
from interop import py_interop_run_metrics, py_interop_run, py_interop_summary

def setup_logging(run_name=None):
    """
    Configure logging for the script with both file and console output.
    
    Args:
        run_name (str, optional): Name of the run for the log file. If None, 
                                 only console logging is configured.
    
    Returns:
        logging.Logger: Configured logger instance
    """
    logger = logging.getLogger(__name__)
    
    # Clear any existing handlers
    logger.handlers.clear()
    
    logger.setLevel(logging.INFO)
    
    # Create formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # File handler (if run_name is provided)
    if run_name:
        # Create logs directory if it doesn't exist
        logs_dir = Path("logs")
        logs_dir.mkdir(exist_ok=True)
        
        # Create log file in logs directory
        log_file = logs_dir / f"{run_name}.log"
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.INFO)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        logger.info(f"Logging to file: {log_file}")
    
    return logger

def get_run_metrics(run_folder):
    """
    Get QC metrics for the sequencing run.
    
    Args:
        run_folder (str): Path to the run folder
        
    Returns:
        tuple: QC metrics (cluster_density, percent_pf, yield_g, error_rate, q30_score)
    """
    logger = logging.getLogger(__name__)
    logger.info("===============================================")
    logger.info("CHECKING RUN QC METRICS")
    logger.info("===============================================")
    
    run_metrics = py_interop_run_metrics.run_metrics()
    valid_to_load = py_interop_run.uchar_vector(py_interop_run.MetricCount, 0)
    py_interop_run_metrics.list_summary_metrics_to_load(valid_to_load)
    run_metrics.read(run_folder, valid_to_load)
    summary = py_interop_summary.run_summary()
    py_interop_summary.summarize_run_metrics(run_metrics, summary)

    # Get metrics for first read and lane
    read = 0
    lane = 0

    metrics = {
        'cluster_density': round(float(summary.at(read).at(lane).density().mean()) / 1000, 2),
        'percent_pf': round(float(summary.at(read).at(lane).percent_pf().mean()), 2),
        'yield_g': round(float(summary.at(read).at(lane).yield_g()), 2),
        'error_rate': round(float(summary.at(read).at(lane).error_rate().mean()), 2),
        'q30_score': round(float(summary.total_summary().percent_gt_q30()), 2)
    }

    # Log the metrics
    logger.info(f"Cluster Density: {metrics['cluster_density']} K/mm²")
    logger.info(f"Clusters Passing Filter: {metrics['percent_pf']}%")
    logger.info(f"Yield: {metrics['yield_g']} G")
    logger.info(f"Error Rate: {metrics['error_rate']}%")
    logger.info(f"% >= Q30: {metrics['q30_score']}%")
    
    return metrics

def check_qc_thresholds(metrics):
    """
    Check if QC metrics meet minimum thresholds.
    
    Args:
        metrics (dict): Dictionary of QC metrics
        
    Returns:
        bool: True if all thresholds are met
    """
    logger = logging.getLogger(__name__)
    
    # Define thresholds (adjust these values as needed)
    thresholds = {
        'percent_pf': 70,     # minimum percentage
        'error_rate': 2,      # maximum error rate percentage
        'q30_score': 80       # minimum Q30 percentage
    }
    
    # Check each metric against its threshold
    failed_metrics = []
    
    if metrics['percent_pf'] < thresholds['percent_pf']:
        failed_metrics.append(f"Clusters Passing Filter ({metrics['percent_pf']}%) below {thresholds['percent_pf']}%")
    
    if metrics['error_rate'] > thresholds['error_rate']:
        failed_metrics.append(f"Error Rate ({metrics['error_rate']}%) above {thresholds['error_rate']}%")
    
    if metrics['q30_score'] < thresholds['q30_score']:
        failed_metrics.append(f"Q30 Score ({metrics['q30_score']}%) below {thresholds['q30_score']}%")
    
    if failed_metrics:
        logger.warning("===============================================")
        logger.warning("QC THRESHOLDS NOT MET:")
        logger.warning("===============================================")
        for failure in failed_metrics:
            logger.warning(failure)
        return False
    
    logger.info("All QC thresholds met!")
    return True

def parse_run_parameters(run_folder):
    """
    Extract experiment name from RunParameters.xml file.
    
    Args:
        run_folder (str): Path to the run folder
        
    Returns:
        str: Experiment name from the XML file
    """
    run_params_path = Path(run_folder) / "RunParameters.xml"
    
    if not run_params_path.exists():
        raise FileNotFoundError(f"RunParameters.xml not found in {run_folder}")
    
    tree = ET.parse(run_params_path)
    root = tree.getroot()
    
    experiment_name = root.find(".//ExperimentName")
    if experiment_name is None or not experiment_name.text:
        raise ValueError("ExperimentName tag not found or empty in RunParameters.xml")
        
    return experiment_name.text

def get_basecalls_path(run_folder):
    """
    Construct path to BaseCalls directory.
    
    Args:
        run_folder (str): Path to the run folder
        
    Returns:
        Path: Path to BaseCalls directory
    """
    basecalls_path = Path(run_folder) / "Data" / "Intensities" / "BaseCalls"
    
    if not basecalls_path.exists():
        raise FileNotFoundError(f"BaseCalls directory not found: {basecalls_path}")
        
    return basecalls_path

def get_adegen_version():
    """
    Extract version from adegen.py script.
    
    Returns:
        str: Version number or "Unknown" if not found
    """
    logger = logging.getLogger(__name__)
    try:
        with open("adegen.py", 'r') as f:
            for line in f:
                if line.startswith("VERSION = "):
                    version = line.split("=")[1].strip().strip('"\'')
                    return version
        return "Unknown"
    except Exception as e:
        logger.warning(f"Could not read adegen.py version: {str(e)}")
        return "Unknown"

def run_adegen(fastq_folder, run_name, dry_run=False):
    """
    Run the adegen.py script with specified parameters and automatically handle confirmations.
    
    Args:
        fastq_folder (Path): Path to FASTQ folder
        run_name (str): Name of the run
        dry_run (bool): If True, only print commands without executing
    """
    json_dir = Path("input_json_files")
    json_dir.mkdir(exist_ok=True)
    
    json_file = json_dir / f"{run_name}.inputs.json"
    version = get_adegen_version()

    cmd = [
        "python3",
        "adegen.py",
        str(fastq_folder),
        "-p", "7043",
        "-r", run_name,
        "-o", str(json_file),
        "-x", "msk_regex.txt",
        "-c"
    ]
    
    logger = logging.getLogger(__name__)
    logger.info("===============================================")
    logger.info("STAGE 1/3 - GENERATING INPUTS JSON")
    logger.info("===============================================")
    logger.info(f"Using adegen.py version: {version}")
    logger.info(f"Command to execute: {' '.join(cmd)}")
    
    if dry_run:
        logger.info("[DRY RUN] Command not executed")
        return
    
    try:
        # Create a process with pipes for input and output
        process = subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
            universal_newlines=True
        )
        
        # Function to handle output in real-time
        def print_output(pipe, prefix=''):
            try:
                for line in iter(pipe.readline, ''):
                    line = line.strip()
                    logger.info(f"{prefix}{line}")
                    process.stdin.write('y\n')
                    process.stdin.flush()
                    

            except Exception as e:
                logger.error(f"Error in output handling: {str(e)}")
        
        # Start threads to handle stdout and stderr in real-time
        from threading import Thread
        stdout_thread = Thread(target=print_output, args=(process.stdout,), daemon=True)
        stderr_thread = Thread(target=print_output, args=(process.stderr, 'ERROR: '), daemon=True)
        
        stdout_thread.start()
        stderr_thread.start()
        
        # Wait for completion with a timeout
        try:
            process.wait(timeout=300)  # 5 minute timeout
        except subprocess.TimeoutExpired:
            process.kill()
            raise TimeoutError("adegen.py process timed out after 5 minutes")
        
        # Check if process was successful
        if process.returncode != 0:
            raise RuntimeError(f"adegen.py failed with exit code {process.returncode}")
        
        logger.info("adegen.py completed successfully")
        
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"adegen.py failed with exit code {e.returncode}")
    except Exception as e:
        # Kill the process if it's still running
        try:
            process.kill()
        except:
            pass
        raise RuntimeError(f"Error running adegen.py: {str(e)}")

def run_upload_wrapper(json_file, dry_run=False):
    """
    Run the sg-upload-v2-wrapper.py script for the upload process.
    
    Args:
        json_file (str): Path to the JSON file
        dry_run (bool): If True, only print commands without executing
    """
    cmd_new = ["python3", "sg-upload-v2-wrapper.py", "new", "--json", json_file]
    cmd_upload = ["python3", "sg-upload-v2-wrapper.py", "upload"]
    
    logger = logging.getLogger(__name__)
    logger.info("===============================================")
    logger.info("STAGE 2/3 - CREATING NEW REQUEST ON SOPHIA DDM")
    logger.info("===============================================")
    logger.info(f"Command to execute: {' '.join(cmd_new)}")
    logger.info("===============================================")
    logger.info("STAGE 3/3 - UPLOADING DATA")
    logger.info("===============================================")
    logger.info(f"Command to execute: {' '.join(cmd_upload)}")
    
    if dry_run:
        logger.info("[DRY RUN] Commands not executed")
        return
    
    try:
        result = subprocess.run(cmd_new, check=True, capture_output=True, text=True)
        subprocess.run(cmd_upload, check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Upload wrapper failed with exit code {e.returncode}")

def main():
    # Initialize logger without run_name first
    logger = setup_logging()
    
    parser = argparse.ArgumentParser(description="Automate sequencing pipeline workflow")
    parser.add_argument("run_folder", help="Path to the run folder")
    parser.add_argument("--skip-upload", action="store_true", 
                       help="Skip the upload step and only generate JSON")
    parser.add_argument("--skip-qc", action="store_true",
                       help="Skip the QC metrics check")
    parser.add_argument("--force", action="store_true",
                       help="Continue even if QC thresholds are not met")
    parser.add_argument("--dry-run", action="store_true",
                       help="Print commands without executing them")
    args = parser.parse_args()
    
    try:
        # Get run name first
        run_name = parse_run_parameters(args.run_folder)
        
        # Reinitialize logger with run_name
        logger = setup_logging(run_name)
        
        # Check QC metrics if not skipped
        if not args.skip_qc and not args.dry_run:
            metrics = get_run_metrics(args.run_folder)
            qc_passed = check_qc_thresholds(metrics)
            
            if not qc_passed and not args.force:
                raise ValueError("QC thresholds not met. Use --force to continue anyway.")
        elif args.dry_run:
            logger.info("[DRY RUN] Skipping QC metrics check")
        else:
            logger.info("Skipping QC metrics check")
        
        fastq_folder = get_basecalls_path(args.run_folder)
        
        # Run adegen.py to generate JSON
        run_adegen(fastq_folder, run_name, dry_run=args.dry_run)
        
        if not args.dry_run:
            # Check if JSON file was created
            json_file = Path("input_json_files") / f"{run_name}.inputs.json"
            if not json_file.exists():
                raise FileNotFoundError(f"Expected output file {json_file} was not created")
            logger.info(f"Successfully generated {json_file}")
        
        # Run upload wrapper if not skipped
        if not args.skip_upload:
            json_file = Path("input_json_files") / f"{run_name}.inputs.json"
            run_upload_wrapper(str(json_file), dry_run=args.dry_run)
            if not args.dry_run:
                logger.info("Upload completed!")
        
        status = "[DRY RUN] " if args.dry_run else ""
        logger.info("===============================================")
        logger.info(f"{status}PIPELINE COMPLETED SUCCESSFULLY")
        logger.info("===============================================")
        
    except Exception as e:
        logger.error(f"Pipeline failed: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()