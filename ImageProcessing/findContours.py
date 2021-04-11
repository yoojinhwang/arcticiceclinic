#!/usr/bin/env python3

import cv2

# Real time testing 
# capture = cv2.VideoCapture(0) 

# Testing from test input file 
capture = cv2.VideoCapture('/home/pi/Desktop/ClinicPictures/test.mp4')

# Background subtraction code
fgbg = cv2.createBackgroundSubtractorMOG2(detectShadows=False)

# Define a codec and create VideoWriter Object 
#fourcc = cv2.VideoWriter_fourcc('M', 'J', 'P', 'G')
#outBackSub = cv2.VideoWriter('backSubtract.mp4', fourcc, 30, (frame_width, frame_height))

while (True): 
	(ret, frame) = capture.read() 
	grayFrame = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
	fgmask = fgbg.apply(grayFrame)
	(thresh, blackAndWhiteFrame) = cv2.threshold(grayFrame,127,255, cv2.THRESH_BINARY)

	# contours find 
	contours, hierarchy = cv2.findContours(blackAndWhiteFrame, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
	print(contours)

	for c in contours: 
	        cv2.drawContours(frame, [c], -1, (0,255,0), 3)

	# Imshow - real time testing
	cv2.imshow('original', frame)
	cv2.imshow('background subtractor', fgmask)

	# TODO: need to fix but should be a key press, keyboard is weird atm
	# for some reason, I need this if statement for the code to run properly
	if cv2.waitKey(1) == 27:
		break 

capture.release()
cv2.destroyAllWindows()

