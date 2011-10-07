@echo ========================================================================
@echo                General Structure Analysis System-II
@echo      by Robert B. Von Dreele, Argonne National Laboratory(C), 2010
@echo  This product includes software developed by the UChicago Argonne, LLC,
@echo             as Operator of Argonne National Laboratory.'
@echo ========================================================================
@
@REM Get this script's directory; make sure that the path ends
@REM    with a single backslash
@set gsasloc=%~dp0\*
@set gsasloc=%gsasloc:\\*=\*%
@set gsasloc=%gsasloc:\*=\%
@
python %gsasloc%\GSASII.py %1

