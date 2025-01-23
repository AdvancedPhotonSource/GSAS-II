# Starts GSAS-II when GitHub repo is not installed into current Python
import os
os.environ["GSASII_YOLO_PATH"] = "True"
from GSASII.GSASIIGUI import main
main()
