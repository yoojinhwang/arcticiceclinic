#!/usr/bin/env python3
from datetime import datetime
import subprocess
import os 

now = datetime.now()
print(now)

# Get the date into an fstring
stringNow = now.strftime('%Y-%m-%d-%H:%M:%S')

# Command 
output_path = os.path.expanduser(f'~/Desktop/ClinicPictures/{stringNow}.jpg')
cmd = f'raspistill -o {output_path}'

# Take the picture 
subprocess.call(cmd.split())

# Sync directory to AWS bin 
sync_cmd = f'aws s3 sync ~/Desktop/ClinicPictures s3://hmc-clinic-2021'
subprocess.call(sync_cmd.split())









