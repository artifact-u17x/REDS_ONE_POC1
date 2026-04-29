@echo off
setlocal enabledelayedexpansion

set URL=earthadvice.org
set OUTPUT_FILE=daily_build.zip
set TOKEN=ZkgZ9N7dI7cI3EnvhkEAleK3mdsJxvaQ

:: Create cellingBatches directory if it doesn't exist
if not exist ".\cellingBatches\" mkdir ".\cellingBatches\"

:: Call the API and store HTTP response code
curl --ssl-no-revoke -X POST  ^
    -F "token=%TOKEN%" ^
    -o "daily_build.zip" ^
    -w "%%{http_code}" ^
    https://%URL%/api/garden-cellings/daily-build > response.txt

:: Read HTTP status code from response.txt
set /p RESPONSE_CODE=<response.txt

echo Response code is:
echo %RESPONSE_CODE%

:: Handle different responses
if "%RESPONSE_CODE%"=="200" (
    echo File downloaded successfully as %OUTPUT_FILE%
        
    :: Get current date and time in format YYYY-MM-DD_HHMMSS
    for /f "tokens=2 delims==" %%I in ('wmic os get localdatetime /format:list') do set datetime=%%I
    set TIMESTAMP=!datetime:~0,4!-!datetime:~4,2!-!datetime:~6,2!_!datetime:~8,2!!datetime:~10,2!!datetime:~12,2!

    :: Extract the zip file
    powershell -command "Expand-Archive -Path '%OUTPUT_FILE%' -DestinationPath '.\temp_extract' -Force"
    
    :: Move and rename the extracted file, with timestamp first, then original filename
    for %%F in (.\temp_extract\*) do (
        set "ORIGINAL_FILENAME=%%~nF"
        move "%%F" ".\cellingBatches\!TIMESTAMP!_!ORIGINAL_FILENAME!.csv"
        echo Extracted file moved to .\cellingBatches\!TIMESTAMP!_!ORIGINAL_FILENAME!.csv
    )
    
    echo Running R model update script...
    "C:\Program Files\R\R-4.5.0\bin\Rscript.exe" .\updateModel.R
    echo R script execution completed.

    REM echo Uploading updated model
    REM copy "newlyUpdatedModel.R" "liveModel10.R"
    REM curl -X POST -F "token=%TOKEN%" -F "script=@.\liveModel10.R" https://%URL%/api/reds-model/1/upload-daily-build
    
    :: Clean up temporary extraction directory
    rmdir /s /q .\temp_extract
    
) else if "%RESPONSE_CODE%"=="403" (
    echo Error: Unauthorized access.
) else if "%RESPONSE_CODE%"=="404" (
    echo Error: File not found.
) else (
    echo Error: Something went wrong. HTTP status: %RESPONSE_CODE%
)

:: Cleanup response file
del response.txt
del %OUTPUT_FILE%

endlocal