@echo off

rem Building with MSVC 14.16 (VS 2017)

set commonCompilerFlags=-FR -MTd -nologo -Gm- -GR- -EHa- -fp:strict -Od -Oi -WX -W4 -wd4201 -wd4189 -wd4091 -FC -Z7

set commonLinkerFlags= -incremental:no -opt:ref 
rem For GSL: 
rem set commonCompilerFlags=/LIBPATH:"D:\vcpkg\vcpkg\installed\x64-windows\lib" gslcblas.lib gsl.lib %commonCompilerFlags%

rem For Intel MKL(also need libiomp5md.lib and libiomp5md.dll, these come e.g. with Intel C++ compiler):
set IntelMKLRoot=D:\IntelMKL
set commonCompilerFlags=/openmp /DMKL_ILP64  /I"%IntelMKLRoot%\mkl\2021.1.1\include" %commonCompilerFlags%
set commonLinkerFlags=/nodefaultlib:vcomp /LIBPATH:"%IntelMKLRoot%\mkl\2021.1.1\lib\intel64" /LIBPATH:"%IntelMKLRoot%\compiler\2021.1.2\windows\compiler\lib\intel64_win" mkl_intel_ilp64.lib mkl_intel_thread.lib mkl_core.lib libiomp5md.lib %commonLinkerFlags%

mkdir .\build
pushd .\build
cl /I..\ %commonCompilerFlags% ..\ci_algorithm.cpp /link %commonLinkerFlags%
popd
