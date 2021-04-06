#!/usr/bin/env python3
from datetime import datetime
import subprocess
import os 
import time

now = datetime.now()

# Get the date into an fstring
stringNow = now.strftime('%Y-%m-%d-%H.%M.%S')

# Command for taking pictures
def take_picture(): 
	output_path = os.path.expanduser(f'~/Desktop/ClinicPictures/{stringNow}.jpg')
	cmd = f'raspistill -o {output_path}'
	return cmd

# Command for recording videos
def record_video(dir): 
	cmd_list = []
	output_path = os.path.expanduser(dir + f'{stringNow}.h264')
	print(output_path)
	cmd1 = 'raspivid -t -15000 -o ' + output_path
	cmd2 = 'MP4Box -add ' + output_path + ' ' + dir + f'{stringNow}.mp4'
	#cmd2 = f'MP4Box -add {stringNow}.h264 {stringNow}.mp4'
	cmd3 = 'rm ' + dir + f'{stringNow}.h264'
	cmd_list.append(cmd1)
	cmd_list.append(cmd2)
	cmd_list.append(cmd3)
	return cmd_list

command_list = record_video(os.path.expanduser('~/Desktop/ClinicPictures/'))

# Take the picture 
#for i in range(len(command_list)):
#	print("index",i)
#	print("command list at i", command_list[i])
#	subprocess.call(command_list[i].split())
#	time.sleep(10)

subprocess.call(command_list[0].split())

# Sync directory to AWS bin 
# sync_cmd = f'aws s3 sync ~/Desktop/ClinicPictures s3://hmc-clinic-2021'
# subprocess.call(sync_cmd.split())









