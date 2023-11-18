import mediapipe as mp
import cv2

mp_drawing = mp.solutions.drawing_utils
mp_pose = mp.solutions.pose 
pose = mp_pose.Pose(min_detection_confidence=0.5, min_tracking_confidence=0.5) 
# input video
cap = cv2.VideoCapture("HW13.mp4")
# video width, height, frames per second
frame_width = int(cap.get(3))
frame_height = int(cap.get(4))
fps = int(cap.get(5))
# initialize video writer
fourcc = cv2.VideoWriter_fourcc(*'XVID')
output = cv2.VideoWriter('HW13_pose.mp4', fourcc, fps, (frame_width, frame_height))

with mp_pose.Pose(
    min_detection_confidence=0.5,
    min_tracking_confidence=0.5) as pose:
    
    while cap.isOpened():
        # keep running if there is a frame
        success, image = cap.read()
        # break if there is no more frames
        if not success:
            break
        # prevent modification of the image
        image.flags.writeable = False
        # convert BGR to RGB
        image = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)
        # run mediapipe
        results = pose.process(image)
        # activate writing on image
        image.flags.writeable = True
        # convert from RGB to BGR
        image = cv2.cvtColor(image, cv2.COLOR_RGB2BGR)
        # write landmarks on image
        mp_drawing.draw_landmarks(image, results.pose_landmarks, mp_pose.POSE_CONNECTIONS)
        output.write(image)

        if cv2.waitKey(20) == ord('q'):
            break
cap.release()
cv2.destroyAllWindows()

