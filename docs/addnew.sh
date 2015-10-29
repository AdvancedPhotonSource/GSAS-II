svn add */*.rst
svn add */*/*/*.html */*/*.html
svn propset svn:mime-type text/html */*/*/*.html */*/*.html
svn propset svn:mime-type application/pdf build/latex/*.pdf
