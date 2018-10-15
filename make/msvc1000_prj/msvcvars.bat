@echo off
REM
REM $Id: msvcvars.bat,v 1.1 2011/02/01 14:28:35 ivanov Exp $
REM

@if not "%VSINSTALLDIR%"=="" goto devenv
@call "%VS100COMNTOOLS%vsvars32.bat"

:devenv

if exist "%VS100COMNTOOLS%..\IDE\VCExpress.*" set DEVENV="%VS100COMNTOOLS%..\IDE\VCExpress"
if exist "%VS100COMNTOOLS%..\IDE\devenv.*" set DEVENV="%VS100COMNTOOLS%..\IDE\devenv"

:end
