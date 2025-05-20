@echo off
REM Intel oneAPI の環境設定
call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"

REM コンパイル設定
set ICC=icx
set FLAG=-Wall -Wunknown-pragmas -fast -O3 -fiopenmp -ipo -xCORE-avx512 -qmkl=sequential
set OBJS=main_bdi_cls_vof_01.c
set NAME=main_bdi_cls_vof_01.exe

REM コンパイル実行
%ICC% -o %NAME% %OBJS% %FLAG%

pause
