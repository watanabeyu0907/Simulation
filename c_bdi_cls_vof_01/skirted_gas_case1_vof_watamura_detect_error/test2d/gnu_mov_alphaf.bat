setlocal enabledelayedexpansion

set dname=png_alphaf
set aname=test2d_003s_
set fname=%~n0
set basedir=%~dp0"
pushd %basedir%

mkdir %dname%

for /l %%i in (0,1,30) do (

	set /A nn=%%i+0
	set nn=00000!nn!
	set nn=!nn:~-6,6!

	gnuplot -e "file00='%fname%'; file01='%aname%'; nn=%%i" %~n0.plt

	platex %aname%!nn!.tex
	dvips %aname%!nn!.dvi -o %aname%!nn!.eps
	convert2 -define eps:use-cropbox=false -density 400 -units PixelsPerInch -alpha off %aname%!nn!.eps %aname%!nn!.png
	del *.aux *.log *.dvi *-inc.eps
	del %aname%!nn!.tex
	del %aname%!nn!.eps
	REM start "" "%dname%\%aname%!nn!.png"
	REM %aname%!nn!.png
	move %aname%!nn!.png %dname%

)

cd %dname%
ffmpeg -y -framerate 10 -r 10 -i %aname%%%06d.png -vcodec libx264 -pix_fmt yuv420p -vf scale=-1:720 -an %aname%.mp4
%aname%.mp4

popd

REM pause