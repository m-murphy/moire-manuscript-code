#!/bin/bash                        #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                    #-- the shell for the job
#$ -o output/                      #-- output directory (fill in)
#$ -e err/                         #-- error directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
#$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=.1G                  #-- submits on nodes with enough free memory (required)
#$ -l scratch=1G                   #-- SGE resources (home and scratch disks)
#$ -l h_rt=336:00:00               #-- runtime limit (see above; this requests 336 hours)
#$ -pe smp 20
#$ -t 1-29	                       #-- remove first '#' to specify the number of
                                   #-- tasks if desired (see Tips section on this page)

# Anything under here can be a bash script

# If you used the -t option above, this same script will be run for each task,
# but with $SGE_TASK_ID set to a different value each time (1-10 in this case).
# The commands below are one way to select a different input (PDB codes in
# this example) for each task.  Note that the bash arrays are indexed from 0,
# while task IDs start at 1, so the first entry in the tasks array variable
# is simply a placeholder

#tasks=(0 andara nyangana rundu zambezi)
#input="${tasks[$SGE_TASK_ID]}"

Rscript run_mcmc_namibia.R

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"          # This is useful for debugging and usage purposes,
                                                  # e.g. "did my job exceed its memory request?"
