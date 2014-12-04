# Generate a PBS script for a job, and general utilities for
# waiting for a job to complete.

from shell_command import shellCommand
import sys
from time import sleep
from tempfile import NamedTemporaryFile
import os

# XXX maybe these values could be supplied in a config file?
QSTAT_MAX_TRIES = 5   # number of times to try qstat before failing
QSTAT_ERROR_DELAY = 1 # seconds to sleep while waiting for qstat to recover
QSTAT_DELAY = 10      # seconds to sleep while waiting for job to complete

# this assumes that qstat info for a job will stick around for a while after the job has finished.
def isJobCompleted(jobID, manager='pbs'):
   if manager == 'slurm':
       jobcmd = "sleep 5; sacct -nbXj %s" % jobID  # sleep to let deamon catch up
   else:
       jobcmd = "qstat -f %s" % jobID
   count = 0
   while True:
       (stdout, stderr, exitStatus) = shellCommand(jobcmd)
       # qstat appears to have worked correctly, we can stop trying.
       if exitStatus == 0 or count >= QSTAT_MAX_TRIES:
           break
       count += 1
       sleep(QSTAT_ERROR_DELAY)
   if exitStatus != 0:
       #return (True, exitStatus)
       raise Exception("'%s' returned non-zero exit status %d times, panicking" % (jobcmd, count))
   else:
       # try to fetch the exit status of the job command from the output of qstat.
       jobState = 'R'  #changed default state to 'R'
       exitStatus = None
       for line in stdout.split('\n'):
           ws = line.split()
           if len(ws) == 3:
               if manager == 'slurm':
                   if ws[0] == str(jobID):
                       # all states that cause job termination (ie Completed)
                       if ws[1].strip() in ['CANCELLED', 'COMPLETED', 'FAILED',
                               'NODE_FAIL', 'PREEMPTED', 'SUSPENDED',
                               'TIMEOUT', 'CANCELLED+']:
                           jobState = 'C'
                           exitStatus = int(ws[2].split(':')[0])  # code:signal
                           # check both exit codes
                           exitStatus = int(ws[2].split(':')[1]) if \
                                   exitStatus == 0 else exitStatus
                           # if cancelled, set exitStatus as 1
                           if ws[1].strip() in ['CANCELLED', 'CANCELLED+']:
                               exitStatus = 1
                       else:
                           jobState = 'R'
               else:
                   if ws[0] == 'job_state' and ws[1] == '=':
                       jobState = ws[2]
                   elif ws[0] == 'exit_status' and ws[1] == '=' and ws[2].isdigit():
                       exitStatus = int(ws[2])
       if jobState.upper() == 'C':
           # Job has completed.
           return (True, exitStatus)
       else:
           # Job has not completed.
           return (False, exitStatus)

# returns exit status of job (or None if it can't be determined)
def waitForJobCompletion(jobID, manager='pbs'):
   isFinished, exitCode = isJobCompleted(jobID, manager)
   while(not isFinished):
       sleep(QSTAT_DELAY)
       isFinished, exitCode = isJobCompleted(jobID, manager)
   return exitCode

# returns exit status of job (or None if it can't be determined)
def runJobAndWait(script, stage, logDir='', manager='pbs', verbose=0):
    jobID = script.launch()
    if manager == 'slurm':
        ext = '.slurm'
        prettyJobID = jobID.split()[3]
        jobID = prettyJobID
    else:
        ext = '.pbs'
        prettyJobID = jobID.split('.')[0]
    logFilename = os.path.join(logDir, stage + '.' + prettyJobID + ext)
    with open(logFilename, 'w') as logFile:
        logFile.write(str(script))
    if verbose > 0:
        print('stage = %s, jobID = %s' % (stage, prettyJobID))
    return waitForJobCompletion(jobID, manager)

# Generate a PBS script for a job.
class PBS_Script(object):
    def __init__(self, command, walltime=None, name=None, memInGB=None, queue='batch', moduleList=None, logDir=None):
        self.command = command
        self.queue = queue
        self.name = name
        self.memInGB = memInGB
        self.walltime = walltime
        self.moduleList = moduleList
        self.logDir = logDir

    # render the job script as a string.
    def __str__(self):
        script = ['#!/bin/bash']
        # XXX fixme
        # should include job id in the output name.
        # should use the proper log directory.
        if self.queue == 'terri-smp':
            script.append('#PBS -q terri')
            script.append('#PBS -l procs=8,tpn=8')
        else:
            script.append('#PBS -q %s' % self.queue)
        if self.logDir:
           script.append('#PBS -o %s' % self.logDir)
           script.append('#PBS -e %s' % self.logDir)
        # should put the name of the file in here if possible
        if self.name:
            script.append('#PBS -N %s' % self.name)
        if self.memInGB:
            if self.queue in ['smp', 'terri-smp']:
                script.append('#PBS -l mem=%sgb' % self.memInGB)
            else:
                script.append('#PBS -l pvmem=%sgb' % self.memInGB)
        if self.walltime:
            script.append('#PBS -l walltime=%s' % self.walltime)
        if type(self.moduleList) == list and len(self.moduleList) > 0:
            for item in self.moduleList:
               script.append('module load %s' % item)
        script.append('cd $PBS_O_WORKDIR')
        script.append(self.command)
        return '\n'.join(script) + '\n'

    # create a temporary file to store the job script and then
    # launch it with qsub.
    def launch(self):
        file = NamedTemporaryFile()
        file.write(str(self))
        file.flush()
        command = 'qsub ' + file.name
        (stdout, stderr, returnCode) = shellCommand(command)
        file.close()
        if returnCode == 0:
            return stdout
        else:
            raise(Exception('qsub command failed with exit status: ' + str(returnCode)))

# Generate a SLURM script for a job.
class SLURM_Script(object):
    '''
        NOTE: Current configuration for barcoo testing only!
    '''
    def __init__(self, command, walltime=None, name=None, memInGB=None, queue='main', moduleList=None, logDir=None):
        self.command = command
        self.queue = queue
        self.name = name
        self.memInGB = memInGB
        self.walltime = walltime
        self.moduleList = moduleList
        self.logDir = logDir

    # render the job script as a string.
    def __str__(self):
        script = ['#!/bin/bash']
        # XXX fixme
        # should include job id in the output name.
        # should use the proper log directory.
        
        if self.queue in ['smp','8core']:
            script.append('#SBATCH -p main')
            script.append('#SBATCH --ntasks=1')
            script.append('#SBATCH --cpus-per-task=8')
        else:
            script.append('#SBATCH -p main')
        if self.logDir:
           script.append('#SBATCH -o %s/%s.o%%j' % (self.logDir, self.name))
           script.append('#SBATCH -e %s/%s.e%%j' % (self.logDir, self.name))
        # should put the name of the file in here if possible
        if self.name:
            script.append('#SBATCH -J %s' % self.name)
        if self.memInGB:
            memInMB = self.memInGB*1024
            if self.queue in ['smp','8core']:
                script.append('#SBATCH --mem=%s' % memInMB)
            else:
                script.append('#SBATCH --mem-per-cpu=%s' % memInMB)
        if self.walltime:
            walltime = self.walltime
            if walltime.count(':') > 2:  # convert D:H:M:S to D-H:M:S
                walltime.replace(':', '-', 1)
            script.append('#SBATCH --time=%s' % walltime)
        if type(self.moduleList) == list and len(self.moduleList) > 0:
            for item in self.moduleList:
               script.append('module load %s' % item)
        script.append(self.command)
        return '\n'.join(script) + '\n'

    # create a temporary file to store the job script and then
    # launch it with sbatch.
    def launch(self):
        file = NamedTemporaryFile(delete=False)
        file.write(str(self))
        file.flush()
        command = 'sbatch ' + file.name
        (stdout, stderr, returnCode) = shellCommand(command)
        file.close()
        if returnCode == 0:
            return stdout
        else:
            raise(Exception('sbatch command failed with exit status: ' + str(returnCode)))
