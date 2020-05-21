:: Copyright(c) 2019 Intel Corporation
:: SPDX-License-Identifier: BSD-2-Clause-Patent
@echo off

setlocal
cd /d "%~dp0"

:: Set defaults to prevent inheriting
set "build=y"
:: Default is debug
set "buildtype=Debug"
:: Default is shared
set "shared=ON"
set "GENERATOR="
:: (cmake -G 2>&1 | Select-String -SimpleMatch '*').Line.Split('=')[0].TrimEnd().Replace('* ','')
:: Default is not building unit tests
set "unittest=OFF"
if NOT -%1-==-- call :args %*
if exist CMakeCache.txt del /f /s /q CMakeCache.txt 1>nul
if exist CMakeFiles rmdir /s /q CMakeFiles 1>nul
if NOT "%GENERATOR%"=="" set GENERATOR=-G"%GENERATOR%"

echo Building in %buildtype% configuration

if NOT "%build%"=="y" echo Generating build files

if "%shared%"=="ON" (
    echo Building shared
) else (
    echo Building static
)

if "%unittest%"=="ON" echo Building unit tests

if "%vs%"=="2019" (
    cmake ../.. %GENERATOR% -A x64 -DCMAKE_INSTALL_PREFIX=%SYSTEMDRIVE%\svt-encoders -DBUILD_SHARED_LIBS=%shared% -DBUILD_TESTING=%unittest% %cmake_eflags% || exit 1
) else (
    cmake ../.. %GENERATOR% -DCMAKE_INSTALL_PREFIX=%SYSTEMDRIVE%\svt-encoders -DBUILD_SHARED_LIBS=%shared% -DBUILD_TESTING=%unittest% %cmake_eflags% || exit 1
)

if "%build%"=="y" cmake --build . --config %buildtype%
goto :EOF

:args
if -%1-==-- (
    exit /b
) else if /I "%1"=="/help" (
    call :help
) else if /I "%1"=="help" (
    call :help
) else if /I "%1"=="clean" (
    echo Cleaning build folder
    for %%i in (*) do if not "%%~i" == "build.bat" del "%%~i"
    for /d %%i in (*) do if not "%%~i" == "build.bat" (
        del /f /s /q "%%~i" 1>nul
        rmdir /s /q "%%~i" 1>nul
    )
    exit
) else if /I "%1"=="2019" (
    echo Generating Visual Studio 2019 solution
    set "GENERATOR=Visual Studio 16 2019"
    set vs=2019
    shift
) else if /I "%1"=="2017" (
    echo Generating Visual Studio 2017 solution
    set "GENERATOR=Visual Studio 15 2017 Win64"
    set vs=2017
    shift
) else if /I "%1"=="2015" (
    echo Generating Visual Studio 2015 solution
    set "GENERATOR=Visual Studio 14 2015 Win64"
    set vs=2015
    shift
) else if /I "%1"=="2013" (
    echo Generating Visual Studio 2013 solution
    echo This is currently not officially supported
    set "GENERATOR=Visual Studio 12 2013 Win64"
    set vs=2013
    shift
) else if /I "%1"=="2012" (
    echo Generating Visual Studio 2012 solution
    echo This is currently not officially supported
    set "GENERATOR=Visual Studio 11 2012 Win64"
    set vs=2012
    shift
) else if /I "%1"=="2010" (
    echo Generating Visual Studio 2010 solution
    echo This is currently not officially supported
    set "GENERATOR=Visual Studio 10 2010 Win64"
    set vs=2010
    shift
) else if /I "%1"=="2008" (
    echo Generating Visual Studio 2008 solution
    echo This is currently not officially supported
    set "GENERATOR=Visual Studio 9 2008 Win64"
    set vs=2008
    shift
) else if /I "%1"=="ninja" (
    echo Generating Ninja files
    echo This is currently not officially supported
    set "GENERATOR=Ninja"
    shift
) else if /I "%1"=="msys" (
    echo Generating MSYS Makefiles
    echo This is currently not officially supported
    set "GENERATOR=MSYS Makefiles"
    shift
) else if /I "%1"=="mingw" (
    echo Generating MinGW Makefiles
    echo This is currently not officially supported
    set "GENERATOR=MinGW Makefiles"
    shift
) else if /I "%1"=="unix" (
    echo Generating Unix Makefiles
    echo This is currently not officially supported
    set "GENERATOR=Unix Makefiles"
    shift
) else if /I "%1"=="release" (
    set "buildtype=Release"
    shift
) else if /I "%1"=="debug" (
    set "buildtype=Debug"
    shift
) else if /I "%1"=="test" (
    set "unittest=ON"
    shift
) else if /I "%1"=="static" (
    set "shared=OFF"
    shift
) else if /I "%1"=="shared" (
    set "shared=ON"
    shift
) else if /I "%1"=="nobuild" (
    set "build=n"
    shift
) else (
    echo Unknown argument "%1"
    call :help
)
goto :args

:help
    echo Batch file to build SVT-AV1 on Windows
    echo Usage: build.bat [2019^|2017^|2015^|clean] [release^|debug] [nobuild] [test] [shared^|static]
    exit
goto :EOF
