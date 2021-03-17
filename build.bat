@echo off

rem MSVC 14.28 (VS 2019)

rem set commonCompilerFlags=-std:c++17 -sanitize=address -EHsc -FR -MTd -nologo -Gm- -GR- -EHa- -fp:precise -Od -Oi -WX -W4 -wd4201 -wd4189 -wd4091 -wd4101 -wd4717 -wd4100 -wd4996 -FC -Z7
set commonCompilerFlags=-std:c++17 -sanitize=address -EHsc -FR -MD -nologo -Gm- -GR- -EHa- -fp:precise -O2 -Oi -WX -W4 -wd4201 -wd4189 -wd4091 -wd4101 -wd4717 -wd4100 -wd4996 -FC -Z7

set commonLinkerFlags=-incremental:no -opt:ref 
rem For GSL: 
rem ## x64 ##
set commonCompilerFlags=/LIBPATH:"D:\vcpkg\vcpkg\installed\x64-windows\lib" "D:\vcpkg\vcpkg\installed\x64-windows\lib\gslcblas.lib" "D:\vcpkg\vcpkg\installed\x64-windows\lib\gsl.lib" %commonCompilerFlags%
rem ## x86 ##
rem set commonCompilerFlags=/LIBPATH:"D:\vcpkg\vcpkg\installed\x86-windows\lib" "D:\vcpkg\vcpkg\installed\x86-windows\lib\gslcblas.lib" "D:\vcpkg\vcpkg\installed\x86-windows\lib\gsl.lib" %commonCompilerFlags%

rem For Intel MKL(also need libiomp5md.lib and libiomp5md.dll, these come e.g. with Intel C++ compiler):
set IntelMKLRoot=D:\IntelMKL
rem ## x64 ##
set commonCompilerFlags=/openmp /DMKL_ILP64  /I"%IntelMKLRoot%\mkl\2021.1.1\include" %commonCompilerFlags%
set commonLinkerFlags=/nodefaultlib:vcomp /LIBPATH:"%IntelMKLRoot%\mkl\2021.1.1\lib\intel64" /LIBPATH:"%IntelMKLRoot%\compiler\2021.1.2\windows\compiler\lib\intel64_win" mkl_intel_ilp64.lib mkl_intel_thread.lib mkl_core.lib libiomp5md.lib %commonLinkerFlags%
rem ## x86 ##
rem set commonCompilerFlags=/openmp /I"%IntelMKLRoot%\mkl\2021.1.1\include" %commonCompilerFlags%
rem set commonLinkerFlags=/nodefaultlib:vcomp /LIBPATH:"%IntelMKLRoot%\mkl\2021.1.1\lib\ia32" /LIBPATH:"%IntelMKLRoot%\compiler\2021.1.2\windows\compiler\lib\ia32_win" mkl_intel_c.lib mkl_intel_thread.lib mkl_core.lib libiomp5md.lib %commonLinkerFlags%

rem For ROOT:
rem set ROOTRoot=D:\root_v6.22.08
rem set commonCompilerFlags= /I"%ROOTRoot%\include" %commonCompilerFlags%
rem set commonLinkerFlags= /LIBPATH:"%ROOTRoot%\lib" libCore.lib libImt.lib libRIO.lib libNet.lib libHist.lib libGraf.lib libGraf3d.lib libGpad.lib libROOTVecOps.lib libTree.lib libTreePlayer.lib libRint.lib libPostscript.lib libMatrix.lib libPhysics.lib libMathCore.lib libThread.lib libROOTDataFrame.lib %commonLinkerFlags%

mkdir .\build
pushd .\build
cl /I..\ %commonCompilerFlags% ..\ci_algorithm.cpp /link %commonLinkerFlags%
popd
