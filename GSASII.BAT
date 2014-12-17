@echo ========================================================================
@echo                General Structure Analysis System-II
@echo              by Robert B. Von Dreele and Brian H. Toby
@echo                Argonne National Laboratory(C), 2010
@echo  This product includes software developed by the UChicago Argonne, LLC,
@echo             as Operator of Argonne National Laboratory.
@echo                            Please cite:
@echo      B.H. Toby and R.B. Von Dreele, J. Appl. Cryst. 46, 544-549 (2013)
@echo                   for small angle use also cite:
@echo      R.B. Von Dreele, J. Appl. Cryst. 47, 1784-9 (2014)
@echo ========================================================================
@
@REM Get this script's directory; make sure that the path ends
@REM    with a single backslash
@set gsasloc=%~dp0\*
@set gsasloc=%gsasloc:\\*=\*%
@set gsasloc=%gsasloc:\*=\%
@rem one can use a specific reference to a python via, e.g.:
@set pythonloc=%D:\Python27\python%
%pythonloc% "%gsasloc%\GSASII.py" %1
@REM To keep the window from disappearing with any error messages
pause

