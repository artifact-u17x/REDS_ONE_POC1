@echo off
:: Change to the directory where this script is located
cd /d "%~dp0"

echo Wrapper about to call script
timeout /t 5

:: Get current date and time for the timestamp
for /f "tokens=2 delims==" %%I in ('wmic os get localdatetime /value') do set datetime=%%I
set TIMESTAMP=%datetime:~0,4%-%datetime:~4,2%-%datetime:~6,2%_%datetime:~8,2%%datetime:~10,2%%datetime:~12,2%

:: Create output directory if it doesn't exist
if not exist "consoleOutput" mkdir "consoleOutput"

:: Define output file path (now relative)
set OUTPUT_FILE=consoleOutput\%TIMESTAMP%_consoleOutput.txt

powershell -command "Start-Transcript -Path '%OUTPUT_FILE%' -Append; & '.\updateModel.bat'; Stop-Transcript"