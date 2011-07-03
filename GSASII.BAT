@echo ==========================================================
@echo                 This file starts GSAS-II
@echo = Run this from the directory where GSAS-II is installed =
@echo ==========================================================
@
@REM Get this script's directory; make sure that the path ends 
@REM    with a single backslash
@set gsasloc=%~dp0\*
@set gsasloc=%gsasloc:\\*=\*%
@set gsasloc=%gsasloc:\*=\%
@
python "%gsasloc%\GSASII.py %1
