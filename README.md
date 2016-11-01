# FastqAlignmentPipeline
Customized fastq alignment pipeline for Orchestra (Harvard Medical School's LSF based clusters)

### Login to orchestra

ssh {username}@hms.harvard.edu:{path to the folder containing this pipeline scripts}

### Load the conversion scripts

module load seq/rcbio/1.6

### Edit parameters 

vim parameters.txt

### Edit job submission parameters or go with the defaults

vim convertToSubmission.sh

### Run the conversion

bash convertToSubmission.sh

### Run the submission

bash toSubmitPipe.sh parameters.txt 


If any intermediate step fails, then fix it and resubmit. It will ask whether to run the previous steps or not. 
Therefore, you wont need to rerun every step from scratch. 
