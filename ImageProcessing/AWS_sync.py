#!/usr/bin/env python3
from datetime import datetime
import subprocess

now = datetime.now()
print(now)

# Get the date into an fstring
stringNow = now.strftime('%Y-%m-%d-%H:%M:%S')

# Command 
cmd = "raspistill -o ~/Desktop/" + stringNow + ".jpg"
print(cmd.split())

# Take the picture 
subprocess.call(cmd.split())











