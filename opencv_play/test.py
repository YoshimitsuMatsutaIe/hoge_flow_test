import cv2
import sys

camera_id = 0
delay = 1
window_name = 'frame'

cap = cv2.VideoCapture(camera_id)

if not cap.isOpened():
    sys.exit()

while True:
    ret, frame = cap.read()
    cv2.imshow(window_name, frame)
    
    
    
    
    
    
    if cv2.waitKey(delay) & 0xFF == ord('q'):
        break

cv2.destroyWindow(window_name)





import cv2

aruco = cv2.aruco
dictionary = aruco.getPredefinedDictionary(aruco.DICT_4X4_50)

cap = cv2.VideoCapture(0)
while True:
    ret, frame = cap.read()
    cv2.imshow("frame", frame)
    corners, ids, rejectedImgPoints = aruco.detectMarkers(frame, dictionary) 
    aruco.drawDetectedMarkers(frame, corners, ids) 
    
cap.release() 
cv2.destroyAllWindows()