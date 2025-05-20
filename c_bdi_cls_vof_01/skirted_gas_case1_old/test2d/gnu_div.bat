setlocal enabledelayedexpansion

set aname=test2d_001s_stat
set fname=%~n0

REM mkdir %dname%
set basedir=%~dp0"
pushd %basedir%

cd
gnuplot -e "file00='%fname%'; file01='%aname%'" "%fname%.plt"

platex %fname%.tex
dvips %fname%.dvi -o %fname%.eps
convert2 -define eps:use-cropbox=false -density 600 -units PixelsPerInch -alpha off %fname%.eps %fname%.png

REM del *.aux *.log *.dvi *-inc.eps
REM del %fname%.tex
REM del %fname%.eps
popd


pause
