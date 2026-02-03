import time
import watchdog
import argparse
import subprocess
import os
import subprocess
import pysam
from watchdog.observers import Observer
from watchdog.observers.polling import PollingObserver
from watchdog.events import PatternMatchingEventHandler

parser = argparse.ArgumentParser(description='arguments for watch_for_new_proj_dirs')
parser.add_argument('-i','--input_dir', \
                    help='The directory you\'re watching for addition of new BAMs', \
                    required=True)
parser.add_argument('-r','--scarlet_dir', \
                    help='The SCARLET directory downloaded from GitHub (https://github.com/graemefox/SCARLET)', \
                    required=True)
parser.add_argument('-t','--threads', \
                    help='Number threads to use', \
                    required=True)
parser.add_argument('-f','--ref_fasta', \
                    help='Path to genome reference', \
                    required=True)
parser.add_argument('-a','--annotations_gtf', \
                    help='Path to annotations gtf', \
                    required=True)
args = parser.parse_args()

# give args snappy names
watch_dir = args.input_dir
scarlet = args.scarlet_dir
threads = args.threads
reference = args.ref_fasta
annotations = args.annotations_gtf

if __name__ == "__main__":
    patterns = ["*.bam"]
    ignore_patterns = []
    ignore_directories = True
    case_sensitive = False
    my_event_handler = PatternMatchingEventHandler(patterns=patterns, \
                           ignore_patterns=ignore_patterns, \
                           ignore_directories=ignore_directories, \
                           case_sensitive=case_sensitive)

def on_created(event):
    # this is triggered as soon as the new bam is created in the watch_dir
    # check that the file transfer is complete before starting analysis
    historicalSize = -1
    while (historicalSize != os.path.getsize(event.src_path)):
        historicalSize = os.path.getsize(event.src_path)
        time.sleep(5)
    # now can start the analysis
    SAMPLE=os.path.basename(event.src_path)
    SAMPLE=SAMPLE.replace(".bam", "")
    sample_outdir = os.path.join(watch_dir, f"{SAMPLE}_output")
    print("Found a new bam to analyse: " + SAMPLE + ".bam")
    index_command = "samtools index -@" + str(threads) + " " + event.src_path
    subprocess.run(index_command, shell=True, capture_output=True, text=True)
    print("Indexing complete")
    update_wf_human_variation_command = "nextflow pull epi2melabs/wf-human-variation"
    subprocess.run(update_wf_human_variation_command, shell=True, capture_output=True, text=True)
    print("Checked for wf-human-variation updates")

    scarlet_command = (
        f"nextflow run {scarlet}main.nf "
        f"-w {os.path.join(sample_outdir, 'work')} "
        f"-with-report {os.path.join(sample_outdir, f'{SAMPLE}_nextflow_report.html')} "
        f"--sample {SAMPLE} "
        f"--bam {event.src_path} "
        f"--outdir {sample_outdir} "
        f"--reference {reference} "
        f"--annotations {annotations}"
    )

    result = subprocess.run(scarlet_command, shell=True, text=True)
    print(result.stdout)
    tidy_up = "rm -rf " + watch_dir + SAMPLE + "_output/work"
    print("SCARLET analysis complete")
    return()

def on_deleted(event):
   return()

def on_moved(event):
    return()

def on_modified(event):
    return()

my_event_handler.on_created = on_created
my_event_handler.on_deleted = on_deleted
my_event_handler.on_modified = on_modified
my_event_handler.on_moved = on_moved

path = watch_dir
go_recursively = False
my_observer = PollingObserver()
my_observer.schedule(my_event_handler, \
                     path, \
                     recursive=go_recursively)

my_observer.start()
try:
    while True:
        time.sleep(1)
except KeyboardInterrupt:
    my_observer.stop()
    my_observer.join()
