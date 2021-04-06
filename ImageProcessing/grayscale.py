#!/usr/bin/env python3

import cv2

# Real time testing 
capture = cv2.VideoCapture(0) 

# Testing from test input file 
# capture = cv2.VideoCapture('/home/pi/Desktop/ClinicPictures/test.mp4')

# Background subtraction code
fgbg = cv2.createBackgroundSubtractorMOG2(detectShadows=False)

# Define the frame width and height
#frame_width = int(capture.get(3)) 
#frame_height = int(capture.get(4))

# Define a codec and create VideoWriter Object 
#fourcc = cv2.VideoWriter_fourcc('M', 'J', 'P', 'G')
#outBackSub = cv2.VideoWriter('backSubtract.mp4', fourcc, 30, (frame_width, frame_height))

while (True): 
	(ret, frame) = capture.read() 
	grayFrame = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
	fgmask = fgbg.apply(grayFrame)
	(thresh, blackAndWhiteFrame) = cv2.threshold(grayFrame,127,255, cv2.THRESH_BINARY)
	# testing adaptive thresholding 
	adaptive_thres = cv2.adaptiveThreshold(grayFrame, 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY_INV, 11, 2)
	gaussian_thres = cv2.adaptiveThreshold(grayFrame, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY, 11, 2)
	#outBackSub.write(fgmask)	

	# Imshow - real time testing
	#cv2.imshow('BW', blackAndWhiteFrame)
	cv2.imshow('original', frame)
	cv2.imshow('adaptive threshold', adaptive_thres)
	cv2.imshow('gaussian threshold', gaussian_thres)
	cv2.imshow('background subtractor', fgmask)

	# TODO: need to fix but should be a key press, keyboard is weird atm
	# for some reason, I need this if statement for the code to run properly
	if cv2.waitKey(1) == 27:
		break 

capture.release()
#outBackSub.release()
cv2.destroyAllWindows()
