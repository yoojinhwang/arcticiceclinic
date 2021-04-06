import cv2

capture = cv2.VideoCapture(0)
#fgbg = cv2.createBackgroundSubtractorMOG2(detectShadows=False) 
fgbg = cv2.createBackgroundSubtractorKNN()

while True: 
	(ret, frame) = capture.read()
	# grayFrame = cv2.cvtColor(frame, cv2.COLORBGR2GRAY) 
	fgmask = fgbg.apply(frame)
	cv2.imshow('original', frame)
	cv2.imshow('mask', fgmask)

	if cv2.waitKey(1) == 27: 
		break

capture.release() 
cv2.destroyAllWindows()
