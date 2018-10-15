/*   mlkludge.h
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*            National Center for Biotechnology Information (NCBI)
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government do not place any restriction on its use or reproduction.
*  We would, however, appreciate having the NCBI and the author cited in
*  any work or product based on this material
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
* ===========================================================================
*
* File Name:  mlkludge.h
*
* Author:  Jonathan Kans
*
* Version Creation Date:   1/30/07
*
* $Revision: 1.1 $
*
* File Description:
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
*
* ==========================================================================
*/


#ifndef _MLKLUDGE_
#define _MLKLUDGE_

#include <objpub.h>
#include <objpubme.h>
#include <objmdrs.h>

/*
#define MlaRequestPtr Mla2RequestPtr
#define TitleMsgPtr TitleMsg2Ptr
#define TitleMsgListPtr TitleMsgList2Ptr
#define MlaBackPtr Mla2BackPtr
*/

#define MlaRequestFree Mla2RequestFree
#define MlaRequestAsnRead Mla22RequestAsnRead
#define MlaRequestAsnWrite Mla2RequestAsnWrite

#define TitleMsgFree TitleMsg2Free
#define TitleMsgNew TitleMsg2New
#define TitleMsgAsnRead TitleMsg2AsnRead
#define TitleMsgAsnWrite TitleMsg2AsnWrite

#define TitleMsgListFree TitleMsgList2Free
#define TitleMsgListNew TitleMsgList2New
#define TitleMsgListAsnRead TitleMsgList2AsnRead
#define TitleMsgListAsnWrite TitleMsgList2AsnWrite

#define MlaBackFree Mla2BackFree
#define MlaBackAsnRead Mla2BackAsnRead
#define MlaBackAsnWrite Mla2BackAsnWrite

#endif /* ndef _MLKLUDGE_ */

