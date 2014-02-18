import itertools, subprocess, time

runs = range(1,11)
#syncsize = [1, 10, 100, 500, 1000]
syncsize = [1000]
offset = [0]
#offset = [0, 10, 20, 30, 40, 50, 60]

def run_cmd(cmd):
    print cmd
    process = subprocess.Popen(cmd, shell=True,
                               stdout=subprocess.PIPE,
                               stdin=subprocess.PIPE)
    result = process.communicate()
    return result

start = time.time()
configs = list(itertools.product(*[runs,syncsize,offset]))
for config in configs:
    command = "python evolver.py "
    command += str(config[0])+" "+str(config[1])+" "+str(config[2])
    run_cmd(command)

print "time taken:", str(time.time()-start)
