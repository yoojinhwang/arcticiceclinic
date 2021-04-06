#!/usr/bin/env python3

import cv2

capture = cv2.VideoCapture(0) 
# Later: test with input file 


while (True): 
	(ret, frame) = capture.read() 
	grayFrame = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
	(thresh, blackAndWhiteFrame) = cv2.threshold(grayFrame,127,255, cv2.THRESH_BINARY)
	# testing adaptive thresholding 
	adaptive_thres = cv2.adaptiveThreshold(grayFrame, 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY_INV, 11, 2)
	gaussian_thres = cv2.adaptiveThreshold(grayFrame, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY, 11, 2)
	cv2.imshow('BW', blackAndWhiteFrame)
	cv2.imshow('original', frame)
	cv2.imshow('adaptive threshold', adaptive_thres)
	cv2.imshow('gaussian threshold', gaussian_thres)
	# TODO: need to fix but should be a key press, keyboard is weird atm
	if cv2.waitKey(1) == 27:
		break 

capture.release()
cv2.destroyAllWindows()
