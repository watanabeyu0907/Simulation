setlocal enabledelayedexpansion

set aname=test2d_003s_stat
set fname=%~n0
set basedir=%~dp0"
pushd %basedir%

REM mkdir %dname%


gnuplot -e "file00='%fname%'; file01='%aname%'" %~n0.plt

platex %fname%.tex
dvips %fname%.dvi -o %fname%.eps
convert2 -define eps:use-cropbox=false -density 600 -units PixelsPerInch -alpha off %fname%.eps %fname%.png
del *.aux *.log *.dvi *-inc.eps
del %fname%.tex

popd

REM del %fname%.eps



