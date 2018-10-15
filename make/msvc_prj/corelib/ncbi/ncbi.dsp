# Microsoft Developer Studio Project File - Name="ncbi" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=ncbi - Win32 DebugDLL
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "ncbi.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "ncbi.mak" CFG="ncbi - Win32 DebugDLL"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "ncbi - Win32 DebugDLL" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe
# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "DebugDLL"
# PROP BASE Intermediate_Dir "DebugDLL"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "DebugDLL"
# PROP Intermediate_Dir "DebugDLL"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MDd /W3 /GX /Z7 /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# SUBTRACT BASE CPP /Fr
# ADD CPP /nologo /MDd /W3 /GR /Z7 /Od /I "..\..\..\.." /I "..\..\..\..\corelib" /I "..\..\..\..\connect" /I "..\..\..\..\asnlib" /D "WIN32" /D "_DEBUG" /D "_MBCS" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo
# Begin Target

# Name "ncbi - Win32 DebugDLL"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\..\..\ctools\asn_connection.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\asnlib\asnbufo.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\asnlib\asndebin.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\asnlib\asnenbin.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\asnlib\asngen.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\asnlib\asngenob.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\asnlib\asnio.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\asnlib\asnlex.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\asnlib\asnlext.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\asnlib\asnout.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\asnlib\asnprint.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\asnlib\asntypes.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\binary.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\drwnflpt.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\drwnfspc.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\drwnmfls.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\gifgen.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\matrix.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_ansi_ext.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_ansi_ext.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_base64.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_buffer.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_connection.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_connector.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_connutil.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_core.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_core_c.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_dispd.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_dispd.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_file_connector.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_ftp_connector.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_heapmgr.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_host_info.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_http_connector.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_lb.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_lb.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_lbsmd_stub.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_local.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_local.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_memory_connector.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_priv.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_priv.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_sendmail.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_server_info.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_service.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_service_connector.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_socket.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_socket_connector.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_types.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_util.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\ncbibs.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\ncbienv.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\ncbierr.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\ncbifile.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\ncbilang.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\ncbimath.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\ncbimem.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\ncbimisc.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\ncbimsg.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\ncbiprop.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\ncbisgml.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\ncbistr.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\ncbithr.c
# ADD CPP /D "NCBI_NOTHREADS_AVAIL"
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\ncbitime.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\util\tables\raw_scoremat.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\regex.c
# ADD CPP /D "REGEX_NCBI" /D "REGEX_MALLOC"
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\togif.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\tree.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\tsprintf.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\urlquery.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_version.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\wwwutils.c
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\..\..\ctools\asn_connection.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\asnlib\asndebin.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\asnlib\asnenbin.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\asnlib\asnio.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\asnlib\asnlex.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\asnlib\asnlext.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\asnlib\asnout.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\asnlib\asnprint.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\asnlib\asntypes.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\binary.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\btree.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\gifgen.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\matrix.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\ncbi.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_base64.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_buffer.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_config.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_connection.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_connector.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_connutil.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_core.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_core_c.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_file_connector.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_ftp_connector.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_heapmgr.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_host_info.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_host_infop.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_http_connector.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_memory_connector.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_sendmail.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_server_info.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_server_infop.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_service.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_service_connector.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_service_misc.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_socket.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_socket_connector.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\ncbibs.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\ncbienv.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\ncbierr.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\ncbifile.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\ncbilang.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\ncbilcl.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\ncbimain.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\ncbimath.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\ncbimem.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\ncbimisc.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\ncbimsg.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\ncbiprop.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\ncbisgml.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\ncbistd.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\ncbistr.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\ncbithr.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\ncbitime.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\ncbiwin.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\util\tables\raw_scoremat.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\regex.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\util\tables\tables_export.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\tree.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\corelib\tsprintf.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\api\undefwin.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\urlquery.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\connect\ncbi_version.h
# End Source File
# End Group
# End Target
# End Project
