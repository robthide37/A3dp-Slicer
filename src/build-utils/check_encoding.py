import os
import subprocess
import sys

print("running 'encoding-check.exe' on all files from '"+sys.argv[1]+"'")

for path, subdirs, files in os.walk(sys.argv[1]):
    for name in files:
        # print('encoding-check.exe '+ 'slic3r '+ os.path.join(path, name))
        subprocess.call(['encoding-check.exe', 'slic3r', os.path.join(path, name)])

