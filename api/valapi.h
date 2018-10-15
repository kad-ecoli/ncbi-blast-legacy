/*   valapi.h
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
* File Name:  valapi.h
*
* Author:  Colleen Bollin
*
* Version Creation Date:   4/8/2009
*
* $Revision: 1.4 $
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

#ifndef _valapi_h_
#define _valapi_h_

#ifdef __cplusplus
extern "C" {
#endif

NLM_EXTERN CommentRulePtr LoadCommentRuleSet (void);

NLM_EXTERN CommentRulePtr GetCommentRuleFromRuleSet (CharPtr prefix);

typedef enum {
  eFieldValid_Valid = 0 ,
  eFieldValid_Invalid,
  eFieldValid_MissingRequiredField,
  eFieldValid_FieldOutOfOrder,
  eFieldValid_DuplicateField,
  eFieldValid_Disallowed
} EFieldValid;


/* error code, field rule violated, value of offending field (if any), extra data provided by user */
typedef void (*StructuredCommentCallback) PROTO ((EFieldValid, FieldRulePtr, UserFieldPtr, UserFieldPtr, Pointer));

NLM_EXTERN EFieldValid 
IsStructuredCommentValidForRule 
(UserObjectPtr uop,
 CommentRulePtr comment_rule,
 StructuredCommentCallback s_callback,
 Pointer s_callback_data);

NLM_EXTERN EFieldValid IsStructuredCommentValid (UserObjectPtr uop, StructuredCommentCallback s_callback, Pointer s_callback_data);

NLM_EXTERN Boolean ReorderStructuredCommentFields (UserObjectPtr uop);


#ifdef __cplusplus 
} 
#endif

#endif
