# script to start GSAS-II when not installed via pixi
import os
os.environ["GSASII_YOLO_PATH"] = "True"
from GSASII.GSASII import main
main()
