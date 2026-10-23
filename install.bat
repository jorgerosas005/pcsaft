@echo off
REM Double-click this file to set the project up. It installs the package in
REM editable mode and runs the test suite. Nothing here touches anything
REM outside this folder.

cd /d "%~dp0"
echo.
echo ============================================
echo   pcsaft - setup
echo   folder: %cd%
echo ============================================
echo.

REM Find Python. The "py" launcher ships with the python.org installer;
REM "python" works when Python is on PATH. Try both.
set PY=
where py >nul 2>&1 && set PY=py
if "%PY%"=="" (where python >nul 2>&1 && set PY=python)

if "%PY%"=="" (
    echo   Python was not found.
    echo.
    echo   Install it from https://www.python.org/downloads/
    echo   and TICK "Add python.exe to PATH" on the first screen.
    echo   Then run this file again.
    echo.
    pause
    exit /b 1
)

echo   Using: %PY%
%PY% --version
echo.

echo   [1/3] installing the package (editable)...
%PY% -m pip install -e . --quiet
if errorlevel 1 goto failed

echo   [2/3] installing pytest...
%PY% -m pip install pytest --quiet
if errorlevel 1 goto failed

echo   [3/3] running the tests...
echo.
%PY% -m pytest tests -q
if errorlevel 1 goto failed

echo.
echo   Done. The package is installed.
echo   Try:   %PY% scripts\01_validate_hexane.py
echo.
pause
exit /b 0

:failed
echo.
echo   Something failed above. Copy the error text and ask about it.
echo.
pause
exit /b 1
