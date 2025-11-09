rem Build GLPK with Microsoft Visual Studio Community 2015

rem NOTE: Make sure that HOME variable specifies correct path
set HOME="D:\VS2022\VC\Auxiliary\Build"

call %HOME%\vcvarsall.bat x64
copy config_VC config.h
D:\VS2022\VC\Tools\MSVC\14.38.33130\bin\Hostx64\x86\nmake.exe /f Makefile_VC
D:\VS2022\VC\Tools\MSVC\14.38.33130\bin\Hostx64\x86\nmake.exe /f Makefile_VC check

pause
