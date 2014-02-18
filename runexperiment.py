# This file is part of Python GRN implementation.
# Architype is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# Architype is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License 
# along with GRN.  If not, see <http://www.gnu.org/licenses/>.
# Author Jonathan Byrne 2014

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
