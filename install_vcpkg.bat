@echo off
SET VCPKG_DIR=%~dp0\vcpkg
IF NOT EXIST "%VCPKG_DIR%" (
    git clone https://github.com/microsoft/vcpkg.git "%VCPKG_DIR%"
    cd /d "%VCPKG_DIR%"
    bootstrap-vcpkg.bat
)