@echo off
REM One-time setup + first pipeline run on the Windows laptop.
REM Double-click this file (or run it from a terminal). Leave the window open until it says DONE.
REM Steps: git init + first commit, renv bootstrap, PRISM test profile, PRISM dev profile.
REM Everything printed here is also written to logs\run_dev_windows_console.log

setlocal
cd /d "%~dp0.."
if not exist logs mkdir logs
set CONSOLE=logs\run_dev_windows_console.log

REM ---- find Rscript (PATH first, else newest under Program Files\R) ----
set "RSCRIPT="
for /f "delims=" %%i in ('where Rscript 2^>nul') do if not defined RSCRIPT set "RSCRIPT=%%i"
if not defined RSCRIPT (
  for /f "delims=" %%d in ('dir /b /ad /o-n "C:\Program Files\R\R-*" 2^>nul') do if not defined RSCRIPT set "RSCRIPT=C:\Program Files\R\%%d\bin\Rscript.exe"
)
if not defined RSCRIPT (
  echo ERROR: Rscript not found >> %CONSOLE%
  echo ERROR: Rscript not found. & pause & exit /b 1
)

echo ==== %DATE% %TIME% start ==== > %CONSOLE%
echo Rscript: %RSCRIPT% >> %CONSOLE%
"%RSCRIPT%" --version >> %CONSOLE% 2>&1
git --version >> %CONSOLE% 2>&1

echo [1/4] git init ...
if not exist .git (
  git init -b main >> %CONSOLE% 2>&1
)

echo [2/4] renv bootstrap (installs packages into the project library; a few minutes) ...
"%RSCRIPT%" setup\renv_bootstrap.R >> %CONSOLE% 2>&1
if errorlevel 1 ( echo renv bootstrap FAILED - see %CONSOLE% & echo renv FAILED >> %CONSOLE% & pause & exit /b 1 )

echo Initial git commit ...
git add -A >> %CONSOLE% 2>&1
git commit -q -m "Project scaffold, config profiles, PRISM daily download step" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" -m "Claude-Session: https://claude.ai/code/session_019psHCLDns3ssvMDDJ6pCcm" >> %CONSOLE% 2>&1
git log --oneline -1 >> %CONSOLE% 2>&1

echo [3/4] PRISM test profile (5 days, 2 variables; about a minute) ...
"%RSCRIPT%" scripts\01_prism_daily.R test >> %CONSOLE% 2>&1
if errorlevel 1 ( echo TEST run FAILED - see %CONSOLE% & echo TEST FAILED >> %CONSOLE% & pause & exit /b 1 )

echo [4/4] PRISM dev profile (Arizona, 2023-2025, 4 variables; roughly 3-4 hours, polite downloads) ...
"%RSCRIPT%" scripts\01_prism_daily.R dev >> %CONSOLE% 2>&1
if errorlevel 1 ( echo DEV run finished with errors - see logs & echo DEV ERRORS >> %CONSOLE% & pause & exit /b 1 )

echo ==== %DATE% %TIME% DONE ==== >> %CONSOLE%
echo DONE. Logs are in the logs folder.
pause
