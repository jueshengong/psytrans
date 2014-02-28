#!/usr/bin/env python



import os, sys, traceback
import getpass
from threading import Thread
from subprocess import *

if(sys.hexversion < 0x03000000):
	import Queue
else:
	import queue as Queue


# svmtrain and gnuplot executable

# global parameters and their default values

fold = 5
c_begin, c_end, c_step = -5,  15, 2
g_begin, g_end, g_step =  3, -15, -2
global dataset_pathname, dataset_title, pass_through_string
global out_filename, png_filename

# experimental

nr_local_worker = 1

# process command line options, set global parameters
def process_options(argv=sys.argv):

    global fold
    global c_begin, c_end, c_step
    global g_begin, g_end, g_step
    

def range_f(begin,end,step):
    # like range, but works on non-integer too
    seq = []
    while True:
        if step > 0 and begin > end: break
        if step < 0 and begin < end: break
        seq.append(begin)
        begin = begin + step
    return seq

def permute_sequence(seq):
    n = len(seq)
    if n <= 1: return seq

    mid = int(n/2)
    left = permute_sequence(seq[:mid])
    right = permute_sequence(seq[mid+1:])

    ret = [seq[mid]]
    while left or right:
        if left: ret.append(left.pop(0))
        if right: ret.append(right.pop(0))

    return ret

def calculate_jobs():
    c_seq = permute_sequence(range_f(c_begin,c_end,c_step))
    g_seq = permute_sequence(range_f(g_begin,g_end,g_step))
    nr_c = float(len(c_seq))
    nr_g = float(len(g_seq))
    i = 0
    j = 0
    jobs = []

    while i < nr_c or j < nr_g:
        if i/nr_c < j/nr_g:
            # increase C resolution
            line = []
            for k in range(0,j):
                line.append((c_seq[i],g_seq[k]))
            i = i + 1
            jobs.append(line)
        else:
            # increase g resolution
            line = []
            for k in range(0,i):
                line.append((c_seq[k],g_seq[j]))
            j = j + 1
            jobs.append(line)
    return jobs

class WorkerStopToken:  # used to notify the worker to stop
        pass

class Worker(Thread):

    def __init__(self,name,job_queue,result_queue):
        Thread.__init__(self)
        self.name = name
        self.job_queue = job_queue
        self.result_queue = result_queue

    def run(self):
        while True:
            (cexp,gexp) = self.job_queue.get()
            if cexp is WorkerStopToken:
                self.job_queue.put((cexp,gexp))
                # print('worker {0} stop.'.format(self.name))
                break
            try:
                rate = self.run_one(2.0**cexp,2.0**gexp)
                if rate is None: raise RuntimeError(RuntimeError("get no rate"))
            except:
                # we failed, let others do that and we just quit
            
                traceback.print_exception(sys.exc_info()[0], sys.exc_info()[1], sys.exc_info()[2])
                
                self.job_queue.put((cexp,gexp))
                print('worker {0} quit.'.format(self.name))
                break
            else:
                self.result_queue.put((self.name,cexp,gexp,rate))

    def run_one(self,c,g):
        cmdline = '{0} -c {1} -g {2} -v {3} {4} {5}'.format \
          (svmtrain_exe,c,g,fold,pass_through_string,dataset_pathname)
        result = Popen(cmdline,shell=True,stdout=PIPE).stdout
        for line in result.readlines():
            if str(line).find("Cross") != -1:
                return float(line.split()[-1][0:-1])

def main():
    # put jobs in queue
    jobs = calculate_jobs()
    job_queue = Queue.Queue(0)
    result_queue = Queue.Queue(0)
    for line in jobs:
        for (c, g) in line:
            job_queue.put((c, g))

    job_queue._put = job_queue.queue.appendleft

    # fire local workers
    for i in range(nr_local_worker):
        Worker('local',job_queue,result_queue).start()

    # gather results
    done_jobs = {}

    result_file = open(out_filename, 'w')
    best_rate = -1
    best_c1,best_g1 = None,None

    for line in jobs:
        for (c,g) in line:
            while (c, g) not in done_jobs:
                (worker,c1,g1,rate) = result_queue.get()
                done_jobs[(c1,g1)] = rate
                result_file.write('{0} {1} {2}\n'.format(c1,g1,rate))
                result_file.flush()
                if (rate > best_rate) or (rate==best_rate and g1==best_g1 and c1<best_c1):
                    best_rate = rate
                    best_c1,best_g1=c1,g1
                    best_c = 2.0**c1
                    best_g = 2.0**g1
                print("[{0}] {1} {2} {3} (best c={4}, g={5}, rate={6})".format \
		    (worker,c1,g1,rate, best_c, best_g, best_rate))

    job_queue.put((WorkerStopToken,None))
    print("{0} {1} {2}".format(best_c, best_g, best_rate))

main()
