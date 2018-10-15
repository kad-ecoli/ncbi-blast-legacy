#include <asn.h>

#define NLM_GENERATED_CODE_PROTO

#include <objmacro.h>

static Boolean loaded = FALSE;

#include <asnmacro.h>

#ifndef NLM_EXTERN_LOADS
#define NLM_EXTERN_LOADS {}
#endif

NLM_EXTERN Boolean LIBCALL
objmacroAsnLoad(void)
{

   if ( ! loaded) {
      NLM_EXTERN_LOADS

      if ( ! AsnLoad ())
      return FALSE;
      loaded = TRUE;
   }

   return TRUE;
}



/**************************************************
*    Generated object loaders for Module NCBI-Macro
*    Generated using ASNCODE Revision: 6.18 at Sep 30, 2011  2:02 PM
*
**************************************************/


/**************************************************
*
*    AECRActionNew()
*
**************************************************/
NLM_EXTERN 
AECRActionPtr LIBCALL
AECRActionNew(void)
{
   AECRActionPtr ptr = MemNew((size_t) sizeof(AECRAction));

   ptr -> also_change_mrna = 0;
   return ptr;

}


/**************************************************
*
*    AECRActionFree()
*
**************************************************/
NLM_EXTERN 
AECRActionPtr LIBCALL
AECRActionFree(AECRActionPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   ActionChoiceFree(ptr -> action);
   ConstraintChoiceSetFree(ptr -> constraint);
   return MemFree(ptr);
}


/**************************************************
*
*    AECRActionAsnRead()
*
**************************************************/
NLM_EXTERN 
AECRActionPtr LIBCALL
AECRActionAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   AECRActionPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* AECRAction ::= (self contained) */
      atp = AsnReadId(aip, amp, AECR_ACTION);
   } else {
      atp = AsnLinkType(orig, AECR_ACTION);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = AECRActionNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == AECR_ACTION_action) {
      ptr -> action = ActionChoiceAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == AECR_ACTION_also_change_mrna) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> also_change_mrna = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == AECR_ACTION_constraint) {
      ptr -> constraint = ConstraintChoiceSetAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = AECRActionFree(ptr);
   goto ret;
}



/**************************************************
*
*    AECRActionAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
AECRActionAsnWrite(AECRActionPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, AECR_ACTION);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> action != NULL) {
      if ( ! ActionChoiceAsnWrite(ptr -> action, aip, AECR_ACTION_action)) {
         goto erret;
      }
   }
   av.boolvalue = ptr -> also_change_mrna;
   retval = AsnWrite(aip, AECR_ACTION_also_change_mrna,  &av);
   if (ptr -> constraint != NULL) {
      if ( ! ConstraintChoiceSetAsnWrite(ptr -> constraint, aip, AECR_ACTION_constraint)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ParseActionNew()
*
**************************************************/
NLM_EXTERN 
ParseActionPtr LIBCALL
ParseActionNew(void)
{
   ParseActionPtr ptr = MemNew((size_t) sizeof(ParseAction));

   ptr -> capitalization = 0;
   ptr -> remove_from_parsed = 0;
   return ptr;

}


/**************************************************
*
*    ParseActionFree()
*
**************************************************/
NLM_EXTERN 
ParseActionPtr LIBCALL
ParseActionFree(ParseActionPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   TextPortionFree(ptr -> portion);
   ParseSrcFree(ptr -> src);
   ParseDestFree(ptr -> dest);
   TextTransformSetFree(ptr -> transform);
   return MemFree(ptr);
}


/**************************************************
*
*    ParseActionAsnRead()
*
**************************************************/
NLM_EXTERN 
ParseActionPtr LIBCALL
ParseActionAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ParseActionPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ParseAction ::= (self contained) */
      atp = AsnReadId(aip, amp, PARSE_ACTION);
   } else {
      atp = AsnLinkType(orig, PARSE_ACTION);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ParseActionNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == PARSE_ACTION_portion) {
      ptr -> portion = TextPortionAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == PARSE_ACTION_src) {
      ptr -> src = ParseSrcAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == PARSE_ACTION_dest) {
      ptr -> dest = ParseDestAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == PARSE_ACTION_capitalization) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> capitalization = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == PARSE_ACTION_remove_from_parsed) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> remove_from_parsed = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == PARSE_ACTION_transform) {
      ptr -> transform = TextTransformSetAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == PARSE_ACTION_existing_text) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> existing_text = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ParseActionFree(ptr);
   goto ret;
}



/**************************************************
*
*    ParseActionAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ParseActionAsnWrite(ParseActionPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, PARSE_ACTION);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> portion != NULL) {
      if ( ! TextPortionAsnWrite(ptr -> portion, aip, PARSE_ACTION_portion)) {
         goto erret;
      }
   }
   if (ptr -> src != NULL) {
      if ( ! ParseSrcAsnWrite(ptr -> src, aip, PARSE_ACTION_src)) {
         goto erret;
      }
   }
   if (ptr -> dest != NULL) {
      if ( ! ParseDestAsnWrite(ptr -> dest, aip, PARSE_ACTION_dest)) {
         goto erret;
      }
   }
   av.intvalue = ptr -> capitalization;
   retval = AsnWrite(aip, PARSE_ACTION_capitalization,  &av);
   av.boolvalue = ptr -> remove_from_parsed;
   retval = AsnWrite(aip, PARSE_ACTION_remove_from_parsed,  &av);
   if (ptr -> transform != NULL) {
      if ( ! TextTransformSetAsnWrite(ptr -> transform, aip, PARSE_ACTION_transform)) {
         goto erret;
      }
   }
   av.intvalue = ptr -> existing_text;
   retval = AsnWrite(aip, PARSE_ACTION_existing_text,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    MacroActionListFree()
*
**************************************************/
NLM_EXTERN 
MacroActionListPtr LIBCALL
MacroActionListFree(MacroActionListPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericChoiceSeqOfFree(ptr, (AsnOptFreeFunc) MacroActionChoiceFree);
   return NULL;
}


/**************************************************
*
*    MacroActionListAsnRead()
*
**************************************************/
NLM_EXTERN 
MacroActionListPtr LIBCALL
MacroActionListAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   MacroActionListPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* MacroActionList ::= (self contained) */
      atp = AsnReadId(aip, amp, MACRO_ACTION_LIST);
   } else {
      atp = AsnLinkType(orig, MACRO_ACTION_LIST);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   func = NULL;

   ptr  = AsnGenericChoiceSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) MacroActionChoiceAsnRead, (AsnOptFreeFunc) MacroActionChoiceFree);
   if (isError && ptr  == NULL) {
      goto erret;
   }



ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = MacroActionListFree(ptr);
   goto ret;
}



/**************************************************
*
*    MacroActionListAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
MacroActionListAsnWrite(MacroActionListPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, MACRO_ACTION_LIST);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   retval = AsnGenericChoiceSeqOfAsnWrite(ptr , (AsnWriteFunc) MacroActionChoiceAsnWrite, aip, atp, MACRO_ACTION_LIST_E);
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    SuspectRuleSetFree()
*
**************************************************/
NLM_EXTERN 
SuspectRuleSetPtr LIBCALL
SuspectRuleSetFree(SuspectRuleSetPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericUserSeqOfFree(ptr,  (AsnOptFreeFunc) SuspectRuleFree);
   return NULL;
}


/**************************************************
*
*    SuspectRuleSetAsnRead()
*
**************************************************/
NLM_EXTERN 
SuspectRuleSetPtr LIBCALL
SuspectRuleSetAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   SuspectRuleSetPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* SuspectRuleSet ::= (self contained) */
      atp = AsnReadId(aip, amp, SUSPECT_RULE_SET);
   } else {
      atp = AsnLinkType(orig, SUSPECT_RULE_SET);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   func = NULL;

   ptr  = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) SuspectRuleAsnRead, (AsnOptFreeFunc) SuspectRuleFree);
   if (isError && ptr  == NULL) {
      goto erret;
   }



ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = SuspectRuleSetFree(ptr);
   goto ret;
}



/**************************************************
*
*    SuspectRuleSetAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
SuspectRuleSetAsnWrite(SuspectRuleSetPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, SUSPECT_RULE_SET);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   retval = AsnGenericUserSeqOfAsnWrite(ptr , (AsnWriteFunc) SuspectRuleAsnWrite, aip, atp, SUSPECT_RULE_SET_E);
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    WordSubstitutionNew()
*
**************************************************/
NLM_EXTERN 
WordSubstitutionPtr LIBCALL
WordSubstitutionNew(void)
{
   WordSubstitutionPtr ptr = MemNew((size_t) sizeof(WordSubstitution));

   ptr -> case_sensitive = 0;
   ptr -> whole_word = 0;
   return ptr;

}


/**************************************************
*
*    WordSubstitutionFree()
*
**************************************************/
NLM_EXTERN 
WordSubstitutionPtr LIBCALL
WordSubstitutionFree(WordSubstitutionPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> word);
   AsnGenericBaseSeqOfFree(ptr -> synonyms ,ASNCODE_PTRVAL_SLOT);
   return MemFree(ptr);
}


/**************************************************
*
*    WordSubstitutionAsnRead()
*
**************************************************/
NLM_EXTERN 
WordSubstitutionPtr LIBCALL
WordSubstitutionAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   WordSubstitutionPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* WordSubstitution ::= (self contained) */
      atp = AsnReadId(aip, amp, WORD_SUBSTITUTION);
   } else {
      atp = AsnLinkType(orig, WORD_SUBSTITUTION);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = WordSubstitutionNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == WORD_SUBSTITUTION_word) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> word = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == WORD_SUBSTITUTION_synonyms) {
      ptr -> synonyms = AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_PTRVAL_SLOT, &isError);
      if (isError && ptr -> synonyms == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == WORD_SUBSTITUTION_case_sensitive) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> case_sensitive = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == WORD_SUBSTITUTION_whole_word) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> whole_word = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = WordSubstitutionFree(ptr);
   goto ret;
}



/**************************************************
*
*    WordSubstitutionAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
WordSubstitutionAsnWrite(WordSubstitutionPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, WORD_SUBSTITUTION);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> word != NULL) {
      av.ptrvalue = ptr -> word;
      retval = AsnWrite(aip, WORD_SUBSTITUTION_word,  &av);
   }
   retval = AsnGenericBaseSeqOfAsnWrite(ptr -> synonyms ,ASNCODE_PTRVAL_SLOT, aip, WORD_SUBSTITUTION_synonyms, WORD_SUBSTITUTION_synonyms_E);
   av.boolvalue = ptr -> case_sensitive;
   retval = AsnWrite(aip, WORD_SUBSTITUTION_case_sensitive,  &av);
   av.boolvalue = ptr -> whole_word;
   retval = AsnWrite(aip, WORD_SUBSTITUTION_whole_word,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    WordSubstitutionSetFree()
*
**************************************************/
NLM_EXTERN 
WordSubstitutionSetPtr LIBCALL
WordSubstitutionSetFree(WordSubstitutionSetPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericUserSeqOfFree(ptr,  (AsnOptFreeFunc) WordSubstitutionFree);
   return NULL;
}


/**************************************************
*
*    WordSubstitutionSetAsnRead()
*
**************************************************/
NLM_EXTERN 
WordSubstitutionSetPtr LIBCALL
WordSubstitutionSetAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   WordSubstitutionSetPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* WordSubstitutionSet ::= (self contained) */
      atp = AsnReadId(aip, amp, WORD_SUBSTITUTION_SET);
   } else {
      atp = AsnLinkType(orig, WORD_SUBSTITUTION_SET);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   func = NULL;

   ptr  = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) WordSubstitutionAsnRead, (AsnOptFreeFunc) WordSubstitutionFree);
   if (isError && ptr  == NULL) {
      goto erret;
   }



ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = WordSubstitutionSetFree(ptr);
   goto ret;
}



/**************************************************
*
*    WordSubstitutionSetAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
WordSubstitutionSetAsnWrite(WordSubstitutionSetPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, WORD_SUBSTITUTION_SET);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   retval = AsnGenericUserSeqOfAsnWrite(ptr , (AsnWriteFunc) WordSubstitutionAsnWrite, aip, atp, WORD_SUBSTITUTION_SET_E);
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    StringConstraintNew()
*
**************************************************/
NLM_EXTERN 
StringConstraintPtr LIBCALL
StringConstraintNew(void)
{
   StringConstraintPtr ptr = MemNew((size_t) sizeof(StringConstraint));

   ptr -> match_location = 1;
   ptr -> case_sensitive = 0;
   ptr -> ignore_space = 0;
   ptr -> ignore_punct = 0;
   ptr -> whole_word = 0;
   ptr -> not_present = 0;
   ptr -> is_all_caps = 0;
   ptr -> is_all_lower = 0;
   ptr -> is_all_punct = 0;
   ptr -> ignore_weasel = 0;
   return ptr;

}


/**************************************************
*
*    StringConstraintFree()
*
**************************************************/
NLM_EXTERN 
StringConstraintPtr LIBCALL
StringConstraintFree(StringConstraintPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> match_text);
   WordSubstitutionSetFree(ptr -> ignore_words);
   return MemFree(ptr);
}


/**************************************************
*
*    StringConstraintAsnRead()
*
**************************************************/
NLM_EXTERN 
StringConstraintPtr LIBCALL
StringConstraintAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   StringConstraintPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* StringConstraint ::= (self contained) */
      atp = AsnReadId(aip, amp, STRING_CONSTRAINT);
   } else {
      atp = AsnLinkType(orig, STRING_CONSTRAINT);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = StringConstraintNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == STRING_CONSTRAINT_match_text) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> match_text = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == STRING_CONSTRAINT_match_location) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> match_location = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == STRING_CONSTRAINT_case_sensitive) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> case_sensitive = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == STRING_CONSTRAINT_ignore_space) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> ignore_space = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == STRING_CONSTRAINT_ignore_punct) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> ignore_punct = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == STRING_CONSTRAINT_ignore_words) {
      ptr -> ignore_words = WordSubstitutionSetAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == STRING_CONSTRAINT_whole_word) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> whole_word = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == STRING_CONSTRAINT_not_present) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> not_present = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == STRING_CONSTRAINT_is_all_caps) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> is_all_caps = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == STRING_CONSTRAINT_is_all_lower) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> is_all_lower = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == STRING_CONSTRAINT_is_all_punct) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> is_all_punct = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == STRING_CONSTRAINT_ignore_weasel) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> ignore_weasel = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = StringConstraintFree(ptr);
   goto ret;
}



/**************************************************
*
*    StringConstraintAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
StringConstraintAsnWrite(StringConstraintPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, STRING_CONSTRAINT);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> match_text != NULL) {
      av.ptrvalue = ptr -> match_text;
      retval = AsnWrite(aip, STRING_CONSTRAINT_match_text,  &av);
   }
   av.intvalue = ptr -> match_location;
   retval = AsnWrite(aip, STRING_CONSTRAINT_match_location,  &av);
   av.boolvalue = ptr -> case_sensitive;
   retval = AsnWrite(aip, STRING_CONSTRAINT_case_sensitive,  &av);
   av.boolvalue = ptr -> ignore_space;
   retval = AsnWrite(aip, STRING_CONSTRAINT_ignore_space,  &av);
   av.boolvalue = ptr -> ignore_punct;
   retval = AsnWrite(aip, STRING_CONSTRAINT_ignore_punct,  &av);
   if (ptr -> ignore_words != NULL) {
      if ( ! WordSubstitutionSetAsnWrite(ptr -> ignore_words, aip, STRING_CONSTRAINT_ignore_words)) {
         goto erret;
      }
   }
   av.boolvalue = ptr -> whole_word;
   retval = AsnWrite(aip, STRING_CONSTRAINT_whole_word,  &av);
   av.boolvalue = ptr -> not_present;
   retval = AsnWrite(aip, STRING_CONSTRAINT_not_present,  &av);
   av.boolvalue = ptr -> is_all_caps;
   retval = AsnWrite(aip, STRING_CONSTRAINT_is_all_caps,  &av);
   av.boolvalue = ptr -> is_all_lower;
   retval = AsnWrite(aip, STRING_CONSTRAINT_is_all_lower,  &av);
   av.boolvalue = ptr -> is_all_punct;
   retval = AsnWrite(aip, STRING_CONSTRAINT_is_all_punct,  &av);
   av.boolvalue = ptr -> ignore_weasel;
   retval = AsnWrite(aip, STRING_CONSTRAINT_ignore_weasel,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    StringConstraintSetFree()
*
**************************************************/
NLM_EXTERN 
StringConstraintSetPtr LIBCALL
StringConstraintSetFree(StringConstraintSetPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericUserSeqOfFree(ptr,  (AsnOptFreeFunc) StringConstraintFree);
   return NULL;
}


/**************************************************
*
*    StringConstraintSetAsnRead()
*
**************************************************/
NLM_EXTERN 
StringConstraintSetPtr LIBCALL
StringConstraintSetAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   StringConstraintSetPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* StringConstraintSet ::= (self contained) */
      atp = AsnReadId(aip, amp, STRING_CONSTRAINT_SET);
   } else {
      atp = AsnLinkType(orig, STRING_CONSTRAINT_SET);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   func = NULL;

   ptr  = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) StringConstraintAsnRead, (AsnOptFreeFunc) StringConstraintFree);
   if (isError && ptr  == NULL) {
      goto erret;
   }



ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = StringConstraintSetFree(ptr);
   goto ret;
}



/**************************************************
*
*    StringConstraintSetAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
StringConstraintSetAsnWrite(StringConstraintSetPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, STRING_CONSTRAINT_SET);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   retval = AsnGenericUserSeqOfAsnWrite(ptr , (AsnWriteFunc) StringConstraintAsnWrite, aip, atp, STRING_CONSTRAINT_SET_E);
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    LocationPosConstraintFree()
*
**************************************************/
NLM_EXTERN 
LocationPosConstraintPtr LIBCALL
LocationPosConstraintFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    LocationPosConstraintAsnRead()
*
**************************************************/
NLM_EXTERN 
LocationPosConstraintPtr LIBCALL
LocationPosConstraintAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* LocationPosConstraint ::= (self contained) */
      atp = AsnReadId(aip, amp, LOCATION_POS_CONSTRAINT);
   } else {
      atp = AsnLinkType(orig, LOCATION_POS_CONSTRAINT);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == LOCATION_POS_CONSTRAINT_dist_from_end) {
      choice = LocationPosConstraint_dist_from_end;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == LOCATION_POS_CONSTRAINT_max_dist_from_end) {
      choice = LocationPosConstraint_max_dist_from_end;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == LOCATION_POS_CONSTRAINT_min_dist_from_end) {
      choice = LocationPosConstraint_min_dist_from_end;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    LocationPosConstraintAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
LocationPosConstraintAsnWrite(LocationPosConstraintPtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, LOCATION_POS_CONSTRAINT);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case LocationPosConstraint_dist_from_end:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, LOCATION_POS_CONSTRAINT_dist_from_end, &av);
      break;
   case LocationPosConstraint_max_dist_from_end:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, LOCATION_POS_CONSTRAINT_max_dist_from_end, &av);
      break;
   case LocationPosConstraint_min_dist_from_end:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, LOCATION_POS_CONSTRAINT_min_dist_from_end, &av);
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    LocationConstraintNew()
*
**************************************************/
NLM_EXTERN 
LocationConstraintPtr LIBCALL
LocationConstraintNew(void)
{
   LocationConstraintPtr ptr = MemNew((size_t) sizeof(LocationConstraint));

   ptr -> strand = 0;
   ptr -> seq_type = 0;
   ptr -> partial5 = 0;
   ptr -> partial3 = 0;
   ptr -> location_type = 0;
   return ptr;

}


/**************************************************
*
*    LocationConstraintFree()
*
**************************************************/
NLM_EXTERN 
LocationConstraintPtr LIBCALL
LocationConstraintFree(LocationConstraintPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   LocationPosConstraintFree(ptr -> end5);
   LocationPosConstraintFree(ptr -> end3);
   return MemFree(ptr);
}


/**************************************************
*
*    LocationConstraintAsnRead()
*
**************************************************/
NLM_EXTERN 
LocationConstraintPtr LIBCALL
LocationConstraintAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   LocationConstraintPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* LocationConstraint ::= (self contained) */
      atp = AsnReadId(aip, amp, LOCATION_CONSTRAINT);
   } else {
      atp = AsnLinkType(orig, LOCATION_CONSTRAINT);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = LocationConstraintNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == LOCATION_CONSTRAINT_strand) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> strand = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == LOCATION_CONSTRAINT_seq_type) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> seq_type = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == LOCATION_CONSTRAINT_partial5) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> partial5 = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == LOCATION_CONSTRAINT_partial3) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> partial3 = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == LOCATION_CONSTRAINT_location_type) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> location_type = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == LOCATION_CONSTRAINT_end5) {
      ptr -> end5 = LocationPosConstraintAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == LOCATION_CONSTRAINT_end3) {
      ptr -> end3 = LocationPosConstraintAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = LocationConstraintFree(ptr);
   goto ret;
}



/**************************************************
*
*    LocationConstraintAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
LocationConstraintAsnWrite(LocationConstraintPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, LOCATION_CONSTRAINT);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> strand;
   retval = AsnWrite(aip, LOCATION_CONSTRAINT_strand,  &av);
   av.intvalue = ptr -> seq_type;
   retval = AsnWrite(aip, LOCATION_CONSTRAINT_seq_type,  &av);
   av.intvalue = ptr -> partial5;
   retval = AsnWrite(aip, LOCATION_CONSTRAINT_partial5,  &av);
   av.intvalue = ptr -> partial3;
   retval = AsnWrite(aip, LOCATION_CONSTRAINT_partial3,  &av);
   av.intvalue = ptr -> location_type;
   retval = AsnWrite(aip, LOCATION_CONSTRAINT_location_type,  &av);
   if (ptr -> end5 != NULL) {
      if ( ! LocationPosConstraintAsnWrite(ptr -> end5, aip, LOCATION_CONSTRAINT_end5)) {
         goto erret;
      }
   }
   if (ptr -> end3 != NULL) {
      if ( ! LocationPosConstraintAsnWrite(ptr -> end3, aip, LOCATION_CONSTRAINT_end3)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    FeatQualLegalValNew()
*
**************************************************/
NLM_EXTERN 
FeatQualLegalValPtr LIBCALL
FeatQualLegalValNew(void)
{
   FeatQualLegalValPtr ptr = MemNew((size_t) sizeof(FeatQualLegalVal));

   return ptr;

}


/**************************************************
*
*    FeatQualLegalValFree()
*
**************************************************/
NLM_EXTERN 
FeatQualLegalValPtr LIBCALL
FeatQualLegalValFree(FeatQualLegalValPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> val);
   return MemFree(ptr);
}


/**************************************************
*
*    FeatQualLegalValAsnRead()
*
**************************************************/
NLM_EXTERN 
FeatQualLegalValPtr LIBCALL
FeatQualLegalValAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   FeatQualLegalValPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* FeatQualLegalVal ::= (self contained) */
      atp = AsnReadId(aip, amp, FEAT_QUAL_LEGAL_VAL);
   } else {
      atp = AsnLinkType(orig, FEAT_QUAL_LEGAL_VAL);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = FeatQualLegalValNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == FEAT_QUAL_LEGAL_VAL_qual) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> qual = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == FEAT_QUAL_LEGAL_VAL_val) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> val = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = FeatQualLegalValFree(ptr);
   goto ret;
}



/**************************************************
*
*    FeatQualLegalValAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
FeatQualLegalValAsnWrite(FeatQualLegalValPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, FEAT_QUAL_LEGAL_VAL);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> qual;
   retval = AsnWrite(aip, FEAT_QUAL_LEGAL_VAL_qual,  &av);
   if (ptr -> val != NULL) {
      av.ptrvalue = ptr -> val;
      retval = AsnWrite(aip, FEAT_QUAL_LEGAL_VAL_val,  &av);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    FeatQualLegalValChoiceFree()
*
**************************************************/
NLM_EXTERN 
FeatQualLegalValChoicePtr LIBCALL
FeatQualLegalValChoiceFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case FeatQualLegalValChoice_qual:
      FeatQualLegalValFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    FeatQualLegalValChoiceAsnRead()
*
**************************************************/
NLM_EXTERN 
FeatQualLegalValChoicePtr LIBCALL
FeatQualLegalValChoiceAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* FeatQualLegalValChoice ::= (self contained) */
      atp = AsnReadId(aip, amp, FEAT_QUAL_LEGAL_VAL_CHOICE);
   } else {
      atp = AsnLinkType(orig, FEAT_QUAL_LEGAL_VAL_CHOICE);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == FEAT_QUAL_LEGAL_VAL_CHOICE_qual) {
      choice = FeatQualLegalValChoice_qual;
      func = (AsnReadFunc) FeatQualLegalValAsnRead;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    FeatQualLegalValChoiceAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
FeatQualLegalValChoiceAsnWrite(FeatQualLegalValChoicePtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, FEAT_QUAL_LEGAL_VAL_CHOICE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case FeatQualLegalValChoice_qual:
      writetype = FEAT_QUAL_LEGAL_VAL_CHOICE_qual;
      func = (AsnWriteFunc) FeatQualLegalValAsnWrite;
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    FeatQualLegalSetFree()
*
**************************************************/
NLM_EXTERN 
FeatQualLegalSetPtr LIBCALL
FeatQualLegalSetFree(FeatQualLegalSetPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericChoiceSeqOfFree(ptr, (AsnOptFreeFunc) FeatQualLegalValChoiceFree);
   return NULL;
}


/**************************************************
*
*    FeatQualLegalSetAsnRead()
*
**************************************************/
NLM_EXTERN 
FeatQualLegalSetPtr LIBCALL
FeatQualLegalSetAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   FeatQualLegalSetPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* FeatQualLegalSet ::= (self contained) */
      atp = AsnReadId(aip, amp, FEAT_QUAL_LEGAL_SET);
   } else {
      atp = AsnLinkType(orig, FEAT_QUAL_LEGAL_SET);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   func = NULL;

   ptr  = AsnGenericChoiceSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) FeatQualLegalValChoiceAsnRead, (AsnOptFreeFunc) FeatQualLegalValChoiceFree);
   if (isError && ptr  == NULL) {
      goto erret;
   }



ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = FeatQualLegalSetFree(ptr);
   goto ret;
}



/**************************************************
*
*    FeatQualLegalSetAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
FeatQualLegalSetAsnWrite(FeatQualLegalSetPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, FEAT_QUAL_LEGAL_SET);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   retval = AsnGenericChoiceSeqOfAsnWrite(ptr , (AsnWriteFunc) FeatQualLegalValChoiceAsnWrite, aip, atp, FEAT_QUAL_LEGAL_SET_E);
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    FeatQualChoiceFree()
*
**************************************************/
NLM_EXTERN 
FeatQualChoicePtr LIBCALL
FeatQualChoiceFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case FeatQualChoice_illegal_qual:
      StringConstraintFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    FeatQualChoiceAsnRead()
*
**************************************************/
NLM_EXTERN 
FeatQualChoicePtr LIBCALL
FeatQualChoiceAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* FeatQualChoice ::= (self contained) */
      atp = AsnReadId(aip, amp, FEAT_QUAL_CHOICE);
   } else {
      atp = AsnLinkType(orig, FEAT_QUAL_CHOICE);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == FEAT_QUAL_CHOICE_legal_qual) {
      choice = FeatQualChoice_legal_qual;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == FEAT_QUAL_CHOICE_illegal_qual) {
      choice = FeatQualChoice_illegal_qual;
      func = (AsnReadFunc) StringConstraintAsnRead;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    FeatQualChoiceAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
FeatQualChoiceAsnWrite(FeatQualChoicePtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, FEAT_QUAL_CHOICE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case FeatQualChoice_legal_qual:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, FEAT_QUAL_CHOICE_legal_qual, &av);
      break;
   case FeatQualChoice_illegal_qual:
      writetype = FEAT_QUAL_CHOICE_illegal_qual;
      func = (AsnWriteFunc) StringConstraintAsnWrite;
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    FeatureFieldNew()
*
**************************************************/
NLM_EXTERN 
FeatureFieldPtr LIBCALL
FeatureFieldNew(void)
{
   FeatureFieldPtr ptr = MemNew((size_t) sizeof(FeatureField));

   return ptr;

}


/**************************************************
*
*    FeatureFieldFree()
*
**************************************************/
NLM_EXTERN 
FeatureFieldPtr LIBCALL
FeatureFieldFree(FeatureFieldPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   FeatQualChoiceFree(ptr -> field);
   return MemFree(ptr);
}


/**************************************************
*
*    FeatureFieldAsnRead()
*
**************************************************/
NLM_EXTERN 
FeatureFieldPtr LIBCALL
FeatureFieldAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   FeatureFieldPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* FeatureField ::= (self contained) */
      atp = AsnReadId(aip, amp, FEATURE_FIELD);
   } else {
      atp = AsnLinkType(orig, FEATURE_FIELD);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = FeatureFieldNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == FEATURE_FIELD_type) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> type = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == FEATURE_FIELD_field) {
      ptr -> field = FeatQualChoiceAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = FeatureFieldFree(ptr);
   goto ret;
}



/**************************************************
*
*    FeatureFieldAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
FeatureFieldAsnWrite(FeatureFieldPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, FEATURE_FIELD);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> type;
   retval = AsnWrite(aip, FEATURE_FIELD_type,  &av);
   if (ptr -> field != NULL) {
      if ( ! FeatQualChoiceAsnWrite(ptr -> field, aip, FEATURE_FIELD_field)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    FeatureFieldLegalNew()
*
**************************************************/
NLM_EXTERN 
FeatureFieldLegalPtr LIBCALL
FeatureFieldLegalNew(void)
{
   FeatureFieldLegalPtr ptr = MemNew((size_t) sizeof(FeatureFieldLegal));

   return ptr;

}


/**************************************************
*
*    FeatureFieldLegalFree()
*
**************************************************/
NLM_EXTERN 
FeatureFieldLegalPtr LIBCALL
FeatureFieldLegalFree(FeatureFieldLegalPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   return MemFree(ptr);
}


/**************************************************
*
*    FeatureFieldLegalAsnRead()
*
**************************************************/
NLM_EXTERN 
FeatureFieldLegalPtr LIBCALL
FeatureFieldLegalAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   FeatureFieldLegalPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* FeatureFieldLegal ::= (self contained) */
      atp = AsnReadId(aip, amp, FEATURE_FIELD_LEGAL);
   } else {
      atp = AsnLinkType(orig, FEATURE_FIELD_LEGAL);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = FeatureFieldLegalNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == FEATURE_FIELD_LEGAL_type) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> type = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == FEATURE_FIELD_LEGAL_field) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> field = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = FeatureFieldLegalFree(ptr);
   goto ret;
}



/**************************************************
*
*    FeatureFieldLegalAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
FeatureFieldLegalAsnWrite(FeatureFieldLegalPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, FEATURE_FIELD_LEGAL);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> type;
   retval = AsnWrite(aip, FEATURE_FIELD_LEGAL_type,  &av);
   av.intvalue = ptr -> field;
   retval = AsnWrite(aip, FEATURE_FIELD_LEGAL_field,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    FeatureFieldPairNew()
*
**************************************************/
NLM_EXTERN 
FeatureFieldPairPtr LIBCALL
FeatureFieldPairNew(void)
{
   FeatureFieldPairPtr ptr = MemNew((size_t) sizeof(FeatureFieldPair));

   return ptr;

}


/**************************************************
*
*    FeatureFieldPairFree()
*
**************************************************/
NLM_EXTERN 
FeatureFieldPairPtr LIBCALL
FeatureFieldPairFree(FeatureFieldPairPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   FeatQualChoiceFree(ptr -> field_from);
   FeatQualChoiceFree(ptr -> field_to);
   return MemFree(ptr);
}


/**************************************************
*
*    FeatureFieldPairAsnRead()
*
**************************************************/
NLM_EXTERN 
FeatureFieldPairPtr LIBCALL
FeatureFieldPairAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   FeatureFieldPairPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* FeatureFieldPair ::= (self contained) */
      atp = AsnReadId(aip, amp, FEATURE_FIELD_PAIR);
   } else {
      atp = AsnLinkType(orig, FEATURE_FIELD_PAIR);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = FeatureFieldPairNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == FEATURE_FIELD_PAIR_type) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> type = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == FEATURE_FIELD_PAIR_field_from) {
      ptr -> field_from = FeatQualChoiceAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == FEATURE_FIELD_PAIR_field_to) {
      ptr -> field_to = FeatQualChoiceAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = FeatureFieldPairFree(ptr);
   goto ret;
}



/**************************************************
*
*    FeatureFieldPairAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
FeatureFieldPairAsnWrite(FeatureFieldPairPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, FEATURE_FIELD_PAIR);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> type;
   retval = AsnWrite(aip, FEATURE_FIELD_PAIR_type,  &av);
   if (ptr -> field_from != NULL) {
      if ( ! FeatQualChoiceAsnWrite(ptr -> field_from, aip, FEATURE_FIELD_PAIR_field_from)) {
         goto erret;
      }
   }
   if (ptr -> field_to != NULL) {
      if ( ! FeatQualChoiceAsnWrite(ptr -> field_to, aip, FEATURE_FIELD_PAIR_field_to)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    RnaFeatTypeFree()
*
**************************************************/
NLM_EXTERN 
RnaFeatTypePtr LIBCALL
RnaFeatTypeFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case RnaFeatType_ncRNA:
      MemFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    RnaFeatTypeAsnRead()
*
**************************************************/
NLM_EXTERN 
RnaFeatTypePtr LIBCALL
RnaFeatTypeAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* RnaFeatType ::= (self contained) */
      atp = AsnReadId(aip, amp, RNA_FEAT_TYPE);
   } else {
      atp = AsnLinkType(orig, RNA_FEAT_TYPE);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == RNA_FEAT_TYPE_any) {
      choice = RnaFeatType_any;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == RNA_FEAT_TYPE_preRNA) {
      choice = RnaFeatType_preRNA;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == RNA_FEAT_TYPE_mRNA) {
      choice = RnaFeatType_mRNA;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == RNA_FEAT_TYPE_tRNA) {
      choice = RnaFeatType_tRNA;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == RNA_FEAT_TYPE_rRNA) {
      choice = RnaFeatType_rRNA;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == RNA_FEAT_TYPE_ncRNA) {
      choice = RnaFeatType_ncRNA;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.ptrvalue = av.ptrvalue;
   }
   else if (atp == RNA_FEAT_TYPE_tmRNA) {
      choice = RnaFeatType_tmRNA;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == RNA_FEAT_TYPE_miscRNA) {
      choice = RnaFeatType_miscRNA;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    RnaFeatTypeAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
RnaFeatTypeAsnWrite(RnaFeatTypePtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, RNA_FEAT_TYPE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case RnaFeatType_any:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, RNA_FEAT_TYPE_any, &av);
      break;
   case RnaFeatType_preRNA:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, RNA_FEAT_TYPE_preRNA, &av);
      break;
   case RnaFeatType_mRNA:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, RNA_FEAT_TYPE_mRNA, &av);
      break;
   case RnaFeatType_tRNA:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, RNA_FEAT_TYPE_tRNA, &av);
      break;
   case RnaFeatType_rRNA:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, RNA_FEAT_TYPE_rRNA, &av);
      break;
   case RnaFeatType_ncRNA:
      av.ptrvalue = anp->data.ptrvalue;
      retval = AsnWrite(aip, RNA_FEAT_TYPE_ncRNA, &av);
      break;
   case RnaFeatType_tmRNA:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, RNA_FEAT_TYPE_tmRNA, &av);
      break;
   case RnaFeatType_miscRNA:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, RNA_FEAT_TYPE_miscRNA, &av);
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    RnaQualNew()
*
**************************************************/
NLM_EXTERN 
RnaQualPtr LIBCALL
RnaQualNew(void)
{
   RnaQualPtr ptr = MemNew((size_t) sizeof(RnaQual));

   return ptr;

}


/**************************************************
*
*    RnaQualFree()
*
**************************************************/
NLM_EXTERN 
RnaQualPtr LIBCALL
RnaQualFree(RnaQualPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   RnaFeatTypeFree(ptr -> type);
   return MemFree(ptr);
}


/**************************************************
*
*    RnaQualAsnRead()
*
**************************************************/
NLM_EXTERN 
RnaQualPtr LIBCALL
RnaQualAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   RnaQualPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* RnaQual ::= (self contained) */
      atp = AsnReadId(aip, amp, RNA_QUAL);
   } else {
      atp = AsnLinkType(orig, RNA_QUAL);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = RnaQualNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == RNA_QUAL_type) {
      ptr -> type = RnaFeatTypeAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == RNA_QUAL_field) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> field = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = RnaQualFree(ptr);
   goto ret;
}



/**************************************************
*
*    RnaQualAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
RnaQualAsnWrite(RnaQualPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, RNA_QUAL);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> type != NULL) {
      if ( ! RnaFeatTypeAsnWrite(ptr -> type, aip, RNA_QUAL_type)) {
         goto erret;
      }
   }
   av.intvalue = ptr -> field;
   retval = AsnWrite(aip, RNA_QUAL_field,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    RnaQualPairNew()
*
**************************************************/
NLM_EXTERN 
RnaQualPairPtr LIBCALL
RnaQualPairNew(void)
{
   RnaQualPairPtr ptr = MemNew((size_t) sizeof(RnaQualPair));

   return ptr;

}


/**************************************************
*
*    RnaQualPairFree()
*
**************************************************/
NLM_EXTERN 
RnaQualPairPtr LIBCALL
RnaQualPairFree(RnaQualPairPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   RnaFeatTypeFree(ptr -> type);
   return MemFree(ptr);
}


/**************************************************
*
*    RnaQualPairAsnRead()
*
**************************************************/
NLM_EXTERN 
RnaQualPairPtr LIBCALL
RnaQualPairAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   RnaQualPairPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* RnaQualPair ::= (self contained) */
      atp = AsnReadId(aip, amp, RNA_QUAL_PAIR);
   } else {
      atp = AsnLinkType(orig, RNA_QUAL_PAIR);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = RnaQualPairNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == RNA_QUAL_PAIR_type) {
      ptr -> type = RnaFeatTypeAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == RNA_QUAL_PAIR_field_from) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> field_from = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == RNA_QUAL_PAIR_field_to) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> field_to = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = RnaQualPairFree(ptr);
   goto ret;
}



/**************************************************
*
*    RnaQualPairAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
RnaQualPairAsnWrite(RnaQualPairPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, RNA_QUAL_PAIR);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> type != NULL) {
      if ( ! RnaFeatTypeAsnWrite(ptr -> type, aip, RNA_QUAL_PAIR_type)) {
         goto erret;
      }
   }
   av.intvalue = ptr -> field_from;
   retval = AsnWrite(aip, RNA_QUAL_PAIR_field_from,  &av);
   av.intvalue = ptr -> field_to;
   retval = AsnWrite(aip, RNA_QUAL_PAIR_field_to,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    SourceQualPairNew()
*
**************************************************/
NLM_EXTERN 
SourceQualPairPtr LIBCALL
SourceQualPairNew(void)
{
   SourceQualPairPtr ptr = MemNew((size_t) sizeof(SourceQualPair));

   return ptr;

}


/**************************************************
*
*    SourceQualPairFree()
*
**************************************************/
NLM_EXTERN 
SourceQualPairPtr LIBCALL
SourceQualPairFree(SourceQualPairPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   return MemFree(ptr);
}


/**************************************************
*
*    SourceQualPairAsnRead()
*
**************************************************/
NLM_EXTERN 
SourceQualPairPtr LIBCALL
SourceQualPairAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   SourceQualPairPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* SourceQualPair ::= (self contained) */
      atp = AsnReadId(aip, amp, SOURCE_QUAL_PAIR);
   } else {
      atp = AsnLinkType(orig, SOURCE_QUAL_PAIR);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = SourceQualPairNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == SOURCE_QUAL_PAIR_field_from) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> field_from = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SOURCE_QUAL_PAIR_field_to) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> field_to = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = SourceQualPairFree(ptr);
   goto ret;
}



/**************************************************
*
*    SourceQualPairAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
SourceQualPairAsnWrite(SourceQualPairPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, SOURCE_QUAL_PAIR);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> field_from;
   retval = AsnWrite(aip, SOURCE_QUAL_PAIR_field_from,  &av);
   av.intvalue = ptr -> field_to;
   retval = AsnWrite(aip, SOURCE_QUAL_PAIR_field_to,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    SourceQualChoiceFree()
*
**************************************************/
NLM_EXTERN 
SourceQualChoicePtr LIBCALL
SourceQualChoiceFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    SourceQualChoiceAsnRead()
*
**************************************************/
NLM_EXTERN 
SourceQualChoicePtr LIBCALL
SourceQualChoiceAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* SourceQualChoice ::= (self contained) */
      atp = AsnReadId(aip, amp, SOURCE_QUAL_CHOICE);
   } else {
      atp = AsnLinkType(orig, SOURCE_QUAL_CHOICE);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == SOURCE_QUAL_CHOICE_textqual) {
      choice = SourceQualChoice_textqual;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == SOURCE_QUAL_CHOICE_location) {
      choice = SourceQualChoice_location;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == SOURCE_QUAL_CHOICE_origin) {
      choice = SourceQualChoice_origin;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == SOURCE_QUAL_CHOICE_gcode) {
      choice = SourceQualChoice_gcode;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == SOURCE_QUAL_CHOICE_mgcode) {
      choice = SourceQualChoice_mgcode;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    SourceQualChoiceAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
SourceQualChoiceAsnWrite(SourceQualChoicePtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, SOURCE_QUAL_CHOICE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case SourceQualChoice_textqual:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, SOURCE_QUAL_CHOICE_textqual, &av);
      break;
   case SourceQualChoice_location:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, SOURCE_QUAL_CHOICE_location, &av);
      break;
   case SourceQualChoice_origin:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, SOURCE_QUAL_CHOICE_origin, &av);
      break;
   case SourceQualChoice_gcode:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, SOURCE_QUAL_CHOICE_gcode, &av);
      break;
   case SourceQualChoice_mgcode:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, SOURCE_QUAL_CHOICE_mgcode, &av);
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    SourceQualTextValNew()
*
**************************************************/
NLM_EXTERN 
SourceQualTextValPtr LIBCALL
SourceQualTextValNew(void)
{
   SourceQualTextValPtr ptr = MemNew((size_t) sizeof(SourceQualTextVal));

   return ptr;

}


/**************************************************
*
*    SourceQualTextValFree()
*
**************************************************/
NLM_EXTERN 
SourceQualTextValPtr LIBCALL
SourceQualTextValFree(SourceQualTextValPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> val);
   return MemFree(ptr);
}


/**************************************************
*
*    SourceQualTextValAsnRead()
*
**************************************************/
NLM_EXTERN 
SourceQualTextValPtr LIBCALL
SourceQualTextValAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   SourceQualTextValPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* SourceQualTextVal ::= (self contained) */
      atp = AsnReadId(aip, amp, SOURCE_QUAL_TEXT_VAL);
   } else {
      atp = AsnLinkType(orig, SOURCE_QUAL_TEXT_VAL);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = SourceQualTextValNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == SOURCE_QUAL_TEXT_VAL_srcqual) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> srcqual = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SOURCE_QUAL_TEXT_VAL_val) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> val = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = SourceQualTextValFree(ptr);
   goto ret;
}



/**************************************************
*
*    SourceQualTextValAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
SourceQualTextValAsnWrite(SourceQualTextValPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, SOURCE_QUAL_TEXT_VAL);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> srcqual;
   retval = AsnWrite(aip, SOURCE_QUAL_TEXT_VAL_srcqual,  &av);
   if (ptr -> val != NULL) {
      av.ptrvalue = ptr -> val;
      retval = AsnWrite(aip, SOURCE_QUAL_TEXT_VAL_val,  &av);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    SourceQualValChoiceFree()
*
**************************************************/
NLM_EXTERN 
SourceQualValChoicePtr LIBCALL
SourceQualValChoiceFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case SourceQualValChoice_textqual:
      SourceQualTextValFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    SourceQualValChoiceAsnRead()
*
**************************************************/
NLM_EXTERN 
SourceQualValChoicePtr LIBCALL
SourceQualValChoiceAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* SourceQualValChoice ::= (self contained) */
      atp = AsnReadId(aip, amp, SOURCE_QUAL_VAL_CHOICE);
   } else {
      atp = AsnLinkType(orig, SOURCE_QUAL_VAL_CHOICE);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == SOURCE_QUAL_VAL_CHOICE_textqual) {
      choice = SourceQualValChoice_textqual;
      func = (AsnReadFunc) SourceQualTextValAsnRead;
   }
   else if (atp == SOURCE_QUAL_VAL_CHOICE_location) {
      choice = SourceQualValChoice_location;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == SOURCE_QUAL_VAL_CHOICE_origin) {
      choice = SourceQualValChoice_origin;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == SOURCE_QUAL_VAL_CHOICE_gcode) {
      choice = SourceQualValChoice_gcode;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == SOURCE_QUAL_VAL_CHOICE_mgcode) {
      choice = SourceQualValChoice_mgcode;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    SourceQualValChoiceAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
SourceQualValChoiceAsnWrite(SourceQualValChoicePtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, SOURCE_QUAL_VAL_CHOICE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case SourceQualValChoice_textqual:
      writetype = SOURCE_QUAL_VAL_CHOICE_textqual;
      func = (AsnWriteFunc) SourceQualTextValAsnWrite;
      break;
   case SourceQualValChoice_location:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, SOURCE_QUAL_VAL_CHOICE_location, &av);
      break;
   case SourceQualValChoice_origin:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, SOURCE_QUAL_VAL_CHOICE_origin, &av);
      break;
   case SourceQualValChoice_gcode:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, SOURCE_QUAL_VAL_CHOICE_gcode, &av);
      break;
   case SourceQualValChoice_mgcode:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, SOURCE_QUAL_VAL_CHOICE_mgcode, &av);
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    SourceQualValSetFree()
*
**************************************************/
NLM_EXTERN 
SourceQualValSetPtr LIBCALL
SourceQualValSetFree(SourceQualValSetPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericChoiceSeqOfFree(ptr, (AsnOptFreeFunc) SourceQualValChoiceFree);
   return NULL;
}


/**************************************************
*
*    SourceQualValSetAsnRead()
*
**************************************************/
NLM_EXTERN 
SourceQualValSetPtr LIBCALL
SourceQualValSetAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   SourceQualValSetPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* SourceQualValSet ::= (self contained) */
      atp = AsnReadId(aip, amp, SOURCE_QUAL_VAL_SET);
   } else {
      atp = AsnLinkType(orig, SOURCE_QUAL_VAL_SET);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   func = NULL;

   ptr  = AsnGenericChoiceSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) SourceQualValChoiceAsnRead, (AsnOptFreeFunc) SourceQualValChoiceFree);
   if (isError && ptr  == NULL) {
      goto erret;
   }



ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = SourceQualValSetFree(ptr);
   goto ret;
}



/**************************************************
*
*    SourceQualValSetAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
SourceQualValSetAsnWrite(SourceQualValSetPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, SOURCE_QUAL_VAL_SET);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   retval = AsnGenericChoiceSeqOfAsnWrite(ptr , (AsnWriteFunc) SourceQualValChoiceAsnWrite, aip, atp, SOURCE_QUAL_VAL_SET_E);
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    CDSGeneProtFieldPairNew()
*
**************************************************/
NLM_EXTERN 
CDSGeneProtFieldPairPtr LIBCALL
CDSGeneProtFieldPairNew(void)
{
   CDSGeneProtFieldPairPtr ptr = MemNew((size_t) sizeof(CDSGeneProtFieldPair));

   return ptr;

}


/**************************************************
*
*    CDSGeneProtFieldPairFree()
*
**************************************************/
NLM_EXTERN 
CDSGeneProtFieldPairPtr LIBCALL
CDSGeneProtFieldPairFree(CDSGeneProtFieldPairPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   return MemFree(ptr);
}


/**************************************************
*
*    CDSGeneProtFieldPairAsnRead()
*
**************************************************/
NLM_EXTERN 
CDSGeneProtFieldPairPtr LIBCALL
CDSGeneProtFieldPairAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   CDSGeneProtFieldPairPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* CDSGeneProtFieldPair ::= (self contained) */
      atp = AsnReadId(aip, amp, CDSGENEPROT_FIELD_PAIR);
   } else {
      atp = AsnLinkType(orig, CDSGENEPROT_FIELD_PAIR);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = CDSGeneProtFieldPairNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == CDSGENEPROT_FIELD_PAIR_field_from) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> field_from = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CDSGENEPROT_FIELD_PAIR_field_to) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> field_to = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = CDSGeneProtFieldPairFree(ptr);
   goto ret;
}



/**************************************************
*
*    CDSGeneProtFieldPairAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
CDSGeneProtFieldPairAsnWrite(CDSGeneProtFieldPairPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, CDSGENEPROT_FIELD_PAIR);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> field_from;
   retval = AsnWrite(aip, CDSGENEPROT_FIELD_PAIR_field_from,  &av);
   av.intvalue = ptr -> field_to;
   retval = AsnWrite(aip, CDSGENEPROT_FIELD_PAIR_field_to,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    MolinfoFieldFree()
*
**************************************************/
NLM_EXTERN 
MolinfoFieldPtr LIBCALL
MolinfoFieldFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    MolinfoFieldAsnRead()
*
**************************************************/
NLM_EXTERN 
MolinfoFieldPtr LIBCALL
MolinfoFieldAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* MolinfoField ::= (self contained) */
      atp = AsnReadId(aip, amp, MOLINFO_FIELD);
   } else {
      atp = AsnLinkType(orig, MOLINFO_FIELD);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == MOLINFO_FIELD_molecule) {
      choice = MolinfoField_molecule;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == MOLINFO_FIELD_technique) {
      choice = MolinfoField_technique;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == MOLINFO_FIELD_completedness) {
      choice = MolinfoField_completedness;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == MOLINFO_FIELD_mol_class) {
      choice = MolinfoField_mol_class;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == MOLINFO_FIELD_topology) {
      choice = MolinfoField_topology;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == MOLINFO_FIELD_strand) {
      choice = MolinfoField_strand;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    MolinfoFieldAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
MolinfoFieldAsnWrite(MolinfoFieldPtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, MOLINFO_FIELD);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case MolinfoField_molecule:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, MOLINFO_FIELD_molecule, &av);
      break;
   case MolinfoField_technique:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, MOLINFO_FIELD_technique, &av);
      break;
   case MolinfoField_completedness:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, MOLINFO_FIELD_completedness, &av);
      break;
   case MolinfoField_mol_class:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, MOLINFO_FIELD_mol_class, &av);
      break;
   case MolinfoField_topology:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, MOLINFO_FIELD_topology, &av);
      break;
   case MolinfoField_strand:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, MOLINFO_FIELD_strand, &av);
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    MolinfoMoleculePairNew()
*
**************************************************/
NLM_EXTERN 
MolinfoMoleculePairPtr LIBCALL
MolinfoMoleculePairNew(void)
{
   MolinfoMoleculePairPtr ptr = MemNew((size_t) sizeof(MolinfoMoleculePair));

   return ptr;

}


/**************************************************
*
*    MolinfoMoleculePairFree()
*
**************************************************/
NLM_EXTERN 
MolinfoMoleculePairPtr LIBCALL
MolinfoMoleculePairFree(MolinfoMoleculePairPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   return MemFree(ptr);
}


/**************************************************
*
*    MolinfoMoleculePairAsnRead()
*
**************************************************/
NLM_EXTERN 
MolinfoMoleculePairPtr LIBCALL
MolinfoMoleculePairAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   MolinfoMoleculePairPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* MolinfoMoleculePair ::= (self contained) */
      atp = AsnReadId(aip, amp, MOLINFO_MOLECULE_PAIR);
   } else {
      atp = AsnLinkType(orig, MOLINFO_MOLECULE_PAIR);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = MolinfoMoleculePairNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == MOLINFO_MOLECULE_PAIR_from) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> from = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == MOLINFO_MOLECULE_PAIR_to) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> to = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = MolinfoMoleculePairFree(ptr);
   goto ret;
}



/**************************************************
*
*    MolinfoMoleculePairAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
MolinfoMoleculePairAsnWrite(MolinfoMoleculePairPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, MOLINFO_MOLECULE_PAIR);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> from;
   retval = AsnWrite(aip, MOLINFO_MOLECULE_PAIR_from,  &av);
   av.intvalue = ptr -> to;
   retval = AsnWrite(aip, MOLINFO_MOLECULE_PAIR_to,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    MolinfoTechniquePairNew()
*
**************************************************/
NLM_EXTERN 
MolinfoTechniquePairPtr LIBCALL
MolinfoTechniquePairNew(void)
{
   MolinfoTechniquePairPtr ptr = MemNew((size_t) sizeof(MolinfoTechniquePair));

   return ptr;

}


/**************************************************
*
*    MolinfoTechniquePairFree()
*
**************************************************/
NLM_EXTERN 
MolinfoTechniquePairPtr LIBCALL
MolinfoTechniquePairFree(MolinfoTechniquePairPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   return MemFree(ptr);
}


/**************************************************
*
*    MolinfoTechniquePairAsnRead()
*
**************************************************/
NLM_EXTERN 
MolinfoTechniquePairPtr LIBCALL
MolinfoTechniquePairAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   MolinfoTechniquePairPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* MolinfoTechniquePair ::= (self contained) */
      atp = AsnReadId(aip, amp, MOLINFO_TECHNIQUE_PAIR);
   } else {
      atp = AsnLinkType(orig, MOLINFO_TECHNIQUE_PAIR);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = MolinfoTechniquePairNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == MOLINFO_TECHNIQUE_PAIR_from) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> from = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == MOLINFO_TECHNIQUE_PAIR_to) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> to = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = MolinfoTechniquePairFree(ptr);
   goto ret;
}



/**************************************************
*
*    MolinfoTechniquePairAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
MolinfoTechniquePairAsnWrite(MolinfoTechniquePairPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, MOLINFO_TECHNIQUE_PAIR);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> from;
   retval = AsnWrite(aip, MOLINFO_TECHNIQUE_PAIR_from,  &av);
   av.intvalue = ptr -> to;
   retval = AsnWrite(aip, MOLINFO_TECHNIQUE_PAIR_to,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    MolinfoCompletednessPairNew()
*
**************************************************/
NLM_EXTERN 
MolinfoCompletednessPairPtr LIBCALL
MolinfoCompletednessPairNew(void)
{
   MolinfoCompletednessPairPtr ptr = MemNew((size_t) sizeof(MolinfoCompletednessPair));

   return ptr;

}


/**************************************************
*
*    MolinfoCompletednessPairFree()
*
**************************************************/
NLM_EXTERN 
MolinfoCompletednessPairPtr LIBCALL
MolinfoCompletednessPairFree(MolinfoCompletednessPairPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   return MemFree(ptr);
}


/**************************************************
*
*    MolinfoCompletednessPairAsnRead()
*
**************************************************/
NLM_EXTERN 
MolinfoCompletednessPairPtr LIBCALL
MolinfoCompletednessPairAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   MolinfoCompletednessPairPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* MolinfoCompletednessPair ::= (self contained) */
      atp = AsnReadId(aip, amp, MOLINFO_COMPLETEDNESS_PAIR);
   } else {
      atp = AsnLinkType(orig, MOLINFO_COMPLETEDNESS_PAIR);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = MolinfoCompletednessPairNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == MOLINFO_COMPLETEDNESS_PAIR_from) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> from = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == MOLINFO_COMPLETEDNESS_PAIR_to) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> to = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = MolinfoCompletednessPairFree(ptr);
   goto ret;
}



/**************************************************
*
*    MolinfoCompletednessPairAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
MolinfoCompletednessPairAsnWrite(MolinfoCompletednessPairPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, MOLINFO_COMPLETEDNESS_PAIR);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> from;
   retval = AsnWrite(aip, MOLINFO_COMPLETEDNESS_PAIR_from,  &av);
   av.intvalue = ptr -> to;
   retval = AsnWrite(aip, MOLINFO_COMPLETEDNESS_PAIR_to,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    MolinfoMolClassPairNew()
*
**************************************************/
NLM_EXTERN 
MolinfoMolClassPairPtr LIBCALL
MolinfoMolClassPairNew(void)
{
   MolinfoMolClassPairPtr ptr = MemNew((size_t) sizeof(MolinfoMolClassPair));

   return ptr;

}


/**************************************************
*
*    MolinfoMolClassPairFree()
*
**************************************************/
NLM_EXTERN 
MolinfoMolClassPairPtr LIBCALL
MolinfoMolClassPairFree(MolinfoMolClassPairPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   return MemFree(ptr);
}


/**************************************************
*
*    MolinfoMolClassPairAsnRead()
*
**************************************************/
NLM_EXTERN 
MolinfoMolClassPairPtr LIBCALL
MolinfoMolClassPairAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   MolinfoMolClassPairPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* MolinfoMolClassPair ::= (self contained) */
      atp = AsnReadId(aip, amp, MOLINFO_MOL_CLASS_PAIR);
   } else {
      atp = AsnLinkType(orig, MOLINFO_MOL_CLASS_PAIR);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = MolinfoMolClassPairNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == MOLINFO_MOL_CLASS_PAIR_from) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> from = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == MOLINFO_MOL_CLASS_PAIR_to) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> to = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = MolinfoMolClassPairFree(ptr);
   goto ret;
}



/**************************************************
*
*    MolinfoMolClassPairAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
MolinfoMolClassPairAsnWrite(MolinfoMolClassPairPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, MOLINFO_MOL_CLASS_PAIR);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> from;
   retval = AsnWrite(aip, MOLINFO_MOL_CLASS_PAIR_from,  &av);
   av.intvalue = ptr -> to;
   retval = AsnWrite(aip, MOLINFO_MOL_CLASS_PAIR_to,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    MolinfoTopologyPairNew()
*
**************************************************/
NLM_EXTERN 
MolinfoTopologyPairPtr LIBCALL
MolinfoTopologyPairNew(void)
{
   MolinfoTopologyPairPtr ptr = MemNew((size_t) sizeof(MolinfoTopologyPair));

   return ptr;

}


/**************************************************
*
*    MolinfoTopologyPairFree()
*
**************************************************/
NLM_EXTERN 
MolinfoTopologyPairPtr LIBCALL
MolinfoTopologyPairFree(MolinfoTopologyPairPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   return MemFree(ptr);
}


/**************************************************
*
*    MolinfoTopologyPairAsnRead()
*
**************************************************/
NLM_EXTERN 
MolinfoTopologyPairPtr LIBCALL
MolinfoTopologyPairAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   MolinfoTopologyPairPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* MolinfoTopologyPair ::= (self contained) */
      atp = AsnReadId(aip, amp, MOLINFO_TOPOLOGY_PAIR);
   } else {
      atp = AsnLinkType(orig, MOLINFO_TOPOLOGY_PAIR);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = MolinfoTopologyPairNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == MOLINFO_TOPOLOGY_PAIR_from) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> from = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == MOLINFO_TOPOLOGY_PAIR_to) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> to = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = MolinfoTopologyPairFree(ptr);
   goto ret;
}



/**************************************************
*
*    MolinfoTopologyPairAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
MolinfoTopologyPairAsnWrite(MolinfoTopologyPairPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, MOLINFO_TOPOLOGY_PAIR);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> from;
   retval = AsnWrite(aip, MOLINFO_TOPOLOGY_PAIR_from,  &av);
   av.intvalue = ptr -> to;
   retval = AsnWrite(aip, MOLINFO_TOPOLOGY_PAIR_to,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    MolinfoStrandPairNew()
*
**************************************************/
NLM_EXTERN 
MolinfoStrandPairPtr LIBCALL
MolinfoStrandPairNew(void)
{
   MolinfoStrandPairPtr ptr = MemNew((size_t) sizeof(MolinfoStrandPair));

   return ptr;

}


/**************************************************
*
*    MolinfoStrandPairFree()
*
**************************************************/
NLM_EXTERN 
MolinfoStrandPairPtr LIBCALL
MolinfoStrandPairFree(MolinfoStrandPairPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   return MemFree(ptr);
}


/**************************************************
*
*    MolinfoStrandPairAsnRead()
*
**************************************************/
NLM_EXTERN 
MolinfoStrandPairPtr LIBCALL
MolinfoStrandPairAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   MolinfoStrandPairPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* MolinfoStrandPair ::= (self contained) */
      atp = AsnReadId(aip, amp, MOLINFO_STRAND_PAIR);
   } else {
      atp = AsnLinkType(orig, MOLINFO_STRAND_PAIR);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = MolinfoStrandPairNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == MOLINFO_STRAND_PAIR_from) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> from = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == MOLINFO_STRAND_PAIR_to) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> to = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = MolinfoStrandPairFree(ptr);
   goto ret;
}



/**************************************************
*
*    MolinfoStrandPairAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
MolinfoStrandPairAsnWrite(MolinfoStrandPairPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, MOLINFO_STRAND_PAIR);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> from;
   retval = AsnWrite(aip, MOLINFO_STRAND_PAIR_from,  &av);
   av.intvalue = ptr -> to;
   retval = AsnWrite(aip, MOLINFO_STRAND_PAIR_to,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    MolinfoFieldPairFree()
*
**************************************************/
NLM_EXTERN 
MolinfoFieldPairPtr LIBCALL
MolinfoFieldPairFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case MolinfoFieldPair_molecule:
      MolinfoMoleculePairFree(anp -> data.ptrvalue);
      break;
   case MolinfoFieldPair_technique:
      MolinfoTechniquePairFree(anp -> data.ptrvalue);
      break;
   case MolinfoFieldPair_completedness:
      MolinfoCompletednessPairFree(anp -> data.ptrvalue);
      break;
   case MolinfoFieldPair_mol_class:
      MolinfoMolClassPairFree(anp -> data.ptrvalue);
      break;
   case MolinfoFieldPair_topology:
      MolinfoTopologyPairFree(anp -> data.ptrvalue);
      break;
   case MolinfoFieldPair_strand:
      MolinfoStrandPairFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    MolinfoFieldPairAsnRead()
*
**************************************************/
NLM_EXTERN 
MolinfoFieldPairPtr LIBCALL
MolinfoFieldPairAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* MolinfoFieldPair ::= (self contained) */
      atp = AsnReadId(aip, amp, MOLINFO_FIELD_PAIR);
   } else {
      atp = AsnLinkType(orig, MOLINFO_FIELD_PAIR);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == MOLINFO_FIELD_PAIR_molecule) {
      choice = MolinfoFieldPair_molecule;
      func = (AsnReadFunc) MolinfoMoleculePairAsnRead;
   }
   else if (atp == MOLINFO_FIELD_PAIR_technique) {
      choice = MolinfoFieldPair_technique;
      func = (AsnReadFunc) MolinfoTechniquePairAsnRead;
   }
   else if (atp == MOLINFO_FIELD_PAIR_completedness) {
      choice = MolinfoFieldPair_completedness;
      func = (AsnReadFunc) MolinfoCompletednessPairAsnRead;
   }
   else if (atp == MOLINFO_FIELD_PAIR_mol_class) {
      choice = MolinfoFieldPair_mol_class;
      func = (AsnReadFunc) MolinfoMolClassPairAsnRead;
   }
   else if (atp == MOLINFO_FIELD_PAIR_topology) {
      choice = MolinfoFieldPair_topology;
      func = (AsnReadFunc) MolinfoTopologyPairAsnRead;
   }
   else if (atp == MOLINFO_FIELD_PAIR_strand) {
      choice = MolinfoFieldPair_strand;
      func = (AsnReadFunc) MolinfoStrandPairAsnRead;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    MolinfoFieldPairAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
MolinfoFieldPairAsnWrite(MolinfoFieldPairPtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, MOLINFO_FIELD_PAIR);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case MolinfoFieldPair_molecule:
      writetype = MOLINFO_FIELD_PAIR_molecule;
      func = (AsnWriteFunc) MolinfoMoleculePairAsnWrite;
      break;
   case MolinfoFieldPair_technique:
      writetype = MOLINFO_FIELD_PAIR_technique;
      func = (AsnWriteFunc) MolinfoTechniquePairAsnWrite;
      break;
   case MolinfoFieldPair_completedness:
      writetype = MOLINFO_FIELD_PAIR_completedness;
      func = (AsnWriteFunc) MolinfoCompletednessPairAsnWrite;
      break;
   case MolinfoFieldPair_mol_class:
      writetype = MOLINFO_FIELD_PAIR_mol_class;
      func = (AsnWriteFunc) MolinfoMolClassPairAsnWrite;
      break;
   case MolinfoFieldPair_topology:
      writetype = MOLINFO_FIELD_PAIR_topology;
      func = (AsnWriteFunc) MolinfoTopologyPairAsnWrite;
      break;
   case MolinfoFieldPair_strand:
      writetype = MOLINFO_FIELD_PAIR_strand;
      func = (AsnWriteFunc) MolinfoStrandPairAsnWrite;
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    MolinfoFieldListFree()
*
**************************************************/
NLM_EXTERN 
MolinfoFieldListPtr LIBCALL
MolinfoFieldListFree(MolinfoFieldListPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericChoiceSeqOfFree(ptr, (AsnOptFreeFunc) MolinfoFieldFree);
   return NULL;
}


/**************************************************
*
*    MolinfoFieldListAsnRead()
*
**************************************************/
NLM_EXTERN 
MolinfoFieldListPtr LIBCALL
MolinfoFieldListAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   MolinfoFieldListPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* MolinfoFieldList ::= (self contained) */
      atp = AsnReadId(aip, amp, MOLINFO_FIELD_LIST);
   } else {
      atp = AsnLinkType(orig, MOLINFO_FIELD_LIST);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   func = NULL;

   ptr  = AsnGenericChoiceSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) MolinfoFieldAsnRead, (AsnOptFreeFunc) MolinfoFieldFree);
   if (isError && ptr  == NULL) {
      goto erret;
   }



ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = MolinfoFieldListFree(ptr);
   goto ret;
}



/**************************************************
*
*    MolinfoFieldListAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
MolinfoFieldListAsnWrite(MolinfoFieldListPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, MOLINFO_FIELD_LIST);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   retval = AsnGenericChoiceSeqOfAsnWrite(ptr , (AsnWriteFunc) MolinfoFieldAsnWrite, aip, atp, MOLINFO_FIELD_LIST_E);
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    MolinfoFieldConstraintNew()
*
**************************************************/
NLM_EXTERN 
MolinfoFieldConstraintPtr LIBCALL
MolinfoFieldConstraintNew(void)
{
   MolinfoFieldConstraintPtr ptr = MemNew((size_t) sizeof(MolinfoFieldConstraint));

   ptr -> is_not = 0;
   return ptr;

}


/**************************************************
*
*    MolinfoFieldConstraintFree()
*
**************************************************/
NLM_EXTERN 
MolinfoFieldConstraintPtr LIBCALL
MolinfoFieldConstraintFree(MolinfoFieldConstraintPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MolinfoFieldFree(ptr -> field);
   return MemFree(ptr);
}


/**************************************************
*
*    MolinfoFieldConstraintAsnRead()
*
**************************************************/
NLM_EXTERN 
MolinfoFieldConstraintPtr LIBCALL
MolinfoFieldConstraintAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   MolinfoFieldConstraintPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* MolinfoFieldConstraint ::= (self contained) */
      atp = AsnReadId(aip, amp, MOLINFO_FIELD_CONSTRAINT);
   } else {
      atp = AsnLinkType(orig, MOLINFO_FIELD_CONSTRAINT);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = MolinfoFieldConstraintNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == MOLINFO_FIELD_CONSTRAINT_field) {
      ptr -> field = MolinfoFieldAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == MOLINFO_FIELD_CONSTRAINT_is_not) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> is_not = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = MolinfoFieldConstraintFree(ptr);
   goto ret;
}



/**************************************************
*
*    MolinfoFieldConstraintAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
MolinfoFieldConstraintAsnWrite(MolinfoFieldConstraintPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, MOLINFO_FIELD_CONSTRAINT);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> field != NULL) {
      if ( ! MolinfoFieldAsnWrite(ptr -> field, aip, MOLINFO_FIELD_CONSTRAINT_field)) {
         goto erret;
      }
   }
   av.boolvalue = ptr -> is_not;
   retval = AsnWrite(aip, MOLINFO_FIELD_CONSTRAINT_is_not,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    StructuredCommentFieldFree()
*
**************************************************/
NLM_EXTERN 
StructuredCommentFieldPtr LIBCALL
StructuredCommentFieldFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case StructuredCommentField_named:
      MemFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    StructuredCommentFieldAsnRead()
*
**************************************************/
NLM_EXTERN 
StructuredCommentFieldPtr LIBCALL
StructuredCommentFieldAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* StructuredCommentField ::= (self contained) */
      atp = AsnReadId(aip, amp, STRUCTURED_COMMENT_FIELD);
   } else {
      atp = AsnLinkType(orig, STRUCTURED_COMMENT_FIELD);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == STRUCTURED_COMMENT_FIELD_database) {
      choice = StructuredCommentField_database;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == STRUCTURED_COMMENT_FIELD_named) {
      choice = StructuredCommentField_named;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.ptrvalue = av.ptrvalue;
   }
   else if (atp == STRUCTURED_COMMENT_FIELD_field_name) {
      choice = StructuredCommentField_field_name;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    StructuredCommentFieldAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
StructuredCommentFieldAsnWrite(StructuredCommentFieldPtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, STRUCTURED_COMMENT_FIELD);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case StructuredCommentField_database:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, STRUCTURED_COMMENT_FIELD_database, &av);
      break;
   case StructuredCommentField_named:
      av.ptrvalue = anp->data.ptrvalue;
      retval = AsnWrite(aip, STRUCTURED_COMMENT_FIELD_named, &av);
      break;
   case StructuredCommentField_field_name:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, STRUCTURED_COMMENT_FIELD_field_name, &av);
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    StructuredCommentFieldPairNew()
*
**************************************************/
NLM_EXTERN 
StructuredCommentFieldPairPtr LIBCALL
StructuredCommentFieldPairNew(void)
{
   StructuredCommentFieldPairPtr ptr = MemNew((size_t) sizeof(StructuredCommentFieldPair));

   return ptr;

}


/**************************************************
*
*    StructuredCommentFieldPairFree()
*
**************************************************/
NLM_EXTERN 
StructuredCommentFieldPairPtr LIBCALL
StructuredCommentFieldPairFree(StructuredCommentFieldPairPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   StructuredCommentFieldFree(ptr -> from);
   StructuredCommentFieldFree(ptr -> to);
   return MemFree(ptr);
}


/**************************************************
*
*    StructuredCommentFieldPairAsnRead()
*
**************************************************/
NLM_EXTERN 
StructuredCommentFieldPairPtr LIBCALL
StructuredCommentFieldPairAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   StructuredCommentFieldPairPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* StructuredCommentFieldPair ::= (self contained) */
      atp = AsnReadId(aip, amp, STRUCTURED_COMMENT_FIELD_PAIR);
   } else {
      atp = AsnLinkType(orig, STRUCTURED_COMMENT_FIELD_PAIR);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = StructuredCommentFieldPairNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == STRUCTURED_COMMENT_FIELD_PAIR_from) {
      ptr -> from = StructuredCommentFieldAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == STRUCTURED_COMMENT_FIELD_PAIR_to) {
      ptr -> to = StructuredCommentFieldAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = StructuredCommentFieldPairFree(ptr);
   goto ret;
}



/**************************************************
*
*    StructuredCommentFieldPairAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
StructuredCommentFieldPairAsnWrite(StructuredCommentFieldPairPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, STRUCTURED_COMMENT_FIELD_PAIR);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> from != NULL) {
      if ( ! StructuredCommentFieldAsnWrite(ptr -> from, aip, STRUCTURED_COMMENT_FIELD_PAIR_from)) {
         goto erret;
      }
   }
   if (ptr -> to != NULL) {
      if ( ! StructuredCommentFieldAsnWrite(ptr -> to, aip, STRUCTURED_COMMENT_FIELD_PAIR_to)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    DBLinkFieldPairNew()
*
**************************************************/
NLM_EXTERN 
DBLinkFieldPairPtr LIBCALL
DBLinkFieldPairNew(void)
{
   DBLinkFieldPairPtr ptr = MemNew((size_t) sizeof(DBLinkFieldPair));

   return ptr;

}


/**************************************************
*
*    DBLinkFieldPairFree()
*
**************************************************/
NLM_EXTERN 
DBLinkFieldPairPtr LIBCALL
DBLinkFieldPairFree(DBLinkFieldPairPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   return MemFree(ptr);
}


/**************************************************
*
*    DBLinkFieldPairAsnRead()
*
**************************************************/
NLM_EXTERN 
DBLinkFieldPairPtr LIBCALL
DBLinkFieldPairAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   DBLinkFieldPairPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* DBLinkFieldPair ::= (self contained) */
      atp = AsnReadId(aip, amp, DBLINK_FIELD_PAIR);
   } else {
      atp = AsnLinkType(orig, DBLINK_FIELD_PAIR);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = DBLinkFieldPairNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == DBLINK_FIELD_PAIR_from) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> from = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == DBLINK_FIELD_PAIR_to) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> to = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = DBLinkFieldPairFree(ptr);
   goto ret;
}



/**************************************************
*
*    DBLinkFieldPairAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
DBLinkFieldPairAsnWrite(DBLinkFieldPairPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, DBLINK_FIELD_PAIR);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> from;
   retval = AsnWrite(aip, DBLINK_FIELD_PAIR_from,  &av);
   av.intvalue = ptr -> to;
   retval = AsnWrite(aip, DBLINK_FIELD_PAIR_to,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    PubFieldConstraintNew()
*
**************************************************/
NLM_EXTERN 
PubFieldConstraintPtr LIBCALL
PubFieldConstraintNew(void)
{
   PubFieldConstraintPtr ptr = MemNew((size_t) sizeof(PubFieldConstraint));

   return ptr;

}


/**************************************************
*
*    PubFieldConstraintFree()
*
**************************************************/
NLM_EXTERN 
PubFieldConstraintPtr LIBCALL
PubFieldConstraintFree(PubFieldConstraintPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   StringConstraintFree(ptr -> constraint);
   return MemFree(ptr);
}


/**************************************************
*
*    PubFieldConstraintAsnRead()
*
**************************************************/
NLM_EXTERN 
PubFieldConstraintPtr LIBCALL
PubFieldConstraintAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   PubFieldConstraintPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* PubFieldConstraint ::= (self contained) */
      atp = AsnReadId(aip, amp, PUB_FIELD_CONSTRAINT);
   } else {
      atp = AsnLinkType(orig, PUB_FIELD_CONSTRAINT);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = PubFieldConstraintNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == PUB_FIELD_CONSTRAINT_field) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> field = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == PUB_FIELD_CONSTRAINT_constraint) {
      ptr -> constraint = StringConstraintAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = PubFieldConstraintFree(ptr);
   goto ret;
}



/**************************************************
*
*    PubFieldConstraintAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
PubFieldConstraintAsnWrite(PubFieldConstraintPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, PUB_FIELD_CONSTRAINT);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> field;
   retval = AsnWrite(aip, PUB_FIELD_CONSTRAINT_field,  &av);
   if (ptr -> constraint != NULL) {
      if ( ! StringConstraintAsnWrite(ptr -> constraint, aip, PUB_FIELD_CONSTRAINT_constraint)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    PubFieldSpecialConstraintTypeFree()
*
**************************************************/
NLM_EXTERN 
PubFieldSpecialConstraintTypePtr LIBCALL
PubFieldSpecialConstraintTypeFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    PubFieldSpecialConstraintTypeAsnRead()
*
**************************************************/
NLM_EXTERN 
PubFieldSpecialConstraintTypePtr LIBCALL
PubFieldSpecialConstraintTypeAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* PubFieldSpecialConstraintType ::= (self contained) */
      atp = AsnReadId(aip, amp, PUB_FIELD_SPECIAL_CONSTRAINT_TYPE);
   } else {
      atp = AsnLinkType(orig, PUB_FIELD_SPECIAL_CONSTRAINT_TYPE);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == PUB_FIELD_SPECIAL_CONSTRAINT_TYPE_is_present) {
      choice = PubFieldSpecialConstraintType_is_present;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == PUB_FIELD_SPECIAL_CONSTRAINT_TYPE_is_not_present) {
      choice = PubFieldSpecialConstraintType_is_not_present;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == PUB_FIELD_SPECIAL_CONSTRAINT_TYPE_is_all_caps) {
      choice = PubFieldSpecialConstraintType_is_all_caps;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == PUB_FIELD_SPECIAL_CONSTRAINT_TYPE_is_all_lower) {
      choice = PubFieldSpecialConstraintType_is_all_lower;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == PUB_FIELD_SPECIAL_CONSTRAINT_TYPE_is_all_punct) {
      choice = PubFieldSpecialConstraintType_is_all_punct;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    PubFieldSpecialConstraintTypeAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
PubFieldSpecialConstraintTypeAsnWrite(PubFieldSpecialConstraintTypePtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, PUB_FIELD_SPECIAL_CONSTRAINT_TYPE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case PubFieldSpecialConstraintType_is_present:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, PUB_FIELD_SPECIAL_CONSTRAINT_TYPE_is_present, &av);
      break;
   case PubFieldSpecialConstraintType_is_not_present:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, PUB_FIELD_SPECIAL_CONSTRAINT_TYPE_is_not_present, &av);
      break;
   case PubFieldSpecialConstraintType_is_all_caps:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, PUB_FIELD_SPECIAL_CONSTRAINT_TYPE_is_all_caps, &av);
      break;
   case PubFieldSpecialConstraintType_is_all_lower:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, PUB_FIELD_SPECIAL_CONSTRAINT_TYPE_is_all_lower, &av);
      break;
   case PubFieldSpecialConstraintType_is_all_punct:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, PUB_FIELD_SPECIAL_CONSTRAINT_TYPE_is_all_punct, &av);
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    PubFieldSpecialConstraintNew()
*
**************************************************/
NLM_EXTERN 
PubFieldSpecialConstraintPtr LIBCALL
PubFieldSpecialConstraintNew(void)
{
   PubFieldSpecialConstraintPtr ptr = MemNew((size_t) sizeof(PubFieldSpecialConstraint));

   return ptr;

}


/**************************************************
*
*    PubFieldSpecialConstraintFree()
*
**************************************************/
NLM_EXTERN 
PubFieldSpecialConstraintPtr LIBCALL
PubFieldSpecialConstraintFree(PubFieldSpecialConstraintPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   PubFieldSpecialConstraintTypeFree(ptr -> constraint);
   return MemFree(ptr);
}


/**************************************************
*
*    PubFieldSpecialConstraintAsnRead()
*
**************************************************/
NLM_EXTERN 
PubFieldSpecialConstraintPtr LIBCALL
PubFieldSpecialConstraintAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   PubFieldSpecialConstraintPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* PubFieldSpecialConstraint ::= (self contained) */
      atp = AsnReadId(aip, amp, PUB_FIELD_SPECIAL_CONSTRAINT);
   } else {
      atp = AsnLinkType(orig, PUB_FIELD_SPECIAL_CONSTRAINT);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = PubFieldSpecialConstraintNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == PUB_FIELD_SPECIAL_CONSTRAINT_field) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> field = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == PUB_FIELD_SPECIAL_CONSTRAINT_constraint) {
      ptr -> constraint = PubFieldSpecialConstraintTypeAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = PubFieldSpecialConstraintFree(ptr);
   goto ret;
}



/**************************************************
*
*    PubFieldSpecialConstraintAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
PubFieldSpecialConstraintAsnWrite(PubFieldSpecialConstraintPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, PUB_FIELD_SPECIAL_CONSTRAINT);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> field;
   retval = AsnWrite(aip, PUB_FIELD_SPECIAL_CONSTRAINT_field,  &av);
   if (ptr -> constraint != NULL) {
      if ( ! PubFieldSpecialConstraintTypeAsnWrite(ptr -> constraint, aip, PUB_FIELD_SPECIAL_CONSTRAINT_constraint)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    PublicationConstraintNew()
*
**************************************************/
NLM_EXTERN 
PublicationConstraintPtr LIBCALL
PublicationConstraintNew(void)
{
   PublicationConstraintPtr ptr = MemNew((size_t) sizeof(PublicationConstraint));

   return ptr;

}


/**************************************************
*
*    PublicationConstraintFree()
*
**************************************************/
NLM_EXTERN 
PublicationConstraintPtr LIBCALL
PublicationConstraintFree(PublicationConstraintPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   PubFieldConstraintFree(ptr -> field);
   PubFieldSpecialConstraintFree(ptr -> special_field);
   return MemFree(ptr);
}


/**************************************************
*
*    PublicationConstraintAsnRead()
*
**************************************************/
NLM_EXTERN 
PublicationConstraintPtr LIBCALL
PublicationConstraintAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   PublicationConstraintPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* PublicationConstraint ::= (self contained) */
      atp = AsnReadId(aip, amp, PUBLICATION_CONSTRAINT);
   } else {
      atp = AsnLinkType(orig, PUBLICATION_CONSTRAINT);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = PublicationConstraintNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == PUBLICATION_CONSTRAINT_type) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> type = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == PUBLICATION_CONSTRAINT_field) {
      ptr -> field = PubFieldConstraintAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == PUBLICATION_CONSTRAINT_special_field) {
      ptr -> special_field = PubFieldSpecialConstraintAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = PublicationConstraintFree(ptr);
   goto ret;
}



/**************************************************
*
*    PublicationConstraintAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
PublicationConstraintAsnWrite(PublicationConstraintPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, PUBLICATION_CONSTRAINT);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> type;
   retval = AsnWrite(aip, PUBLICATION_CONSTRAINT_type,  &av);
   if (ptr -> field != NULL) {
      if ( ! PubFieldConstraintAsnWrite(ptr -> field, aip, PUBLICATION_CONSTRAINT_field)) {
         goto erret;
      }
   }
   if (ptr -> special_field != NULL) {
      if ( ! PubFieldSpecialConstraintAsnWrite(ptr -> special_field, aip, PUBLICATION_CONSTRAINT_special_field)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    SourceConstraintNew()
*
**************************************************/
NLM_EXTERN 
SourceConstraintPtr LIBCALL
SourceConstraintNew(void)
{
   SourceConstraintPtr ptr = MemNew((size_t) sizeof(SourceConstraint));

   return ptr;

}


/**************************************************
*
*    SourceConstraintFree()
*
**************************************************/
NLM_EXTERN 
SourceConstraintPtr LIBCALL
SourceConstraintFree(SourceConstraintPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   SourceQualChoiceFree(ptr -> field1);
   SourceQualChoiceFree(ptr -> field2);
   StringConstraintFree(ptr -> constraint);
   return MemFree(ptr);
}


/**************************************************
*
*    SourceConstraintAsnRead()
*
**************************************************/
NLM_EXTERN 
SourceConstraintPtr LIBCALL
SourceConstraintAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   SourceConstraintPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* SourceConstraint ::= (self contained) */
      atp = AsnReadId(aip, amp, SOURCE_CONSTRAINT);
   } else {
      atp = AsnLinkType(orig, SOURCE_CONSTRAINT);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = SourceConstraintNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == SOURCE_CONSTRAINT_field1) {
      ptr -> field1 = SourceQualChoiceAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SOURCE_CONSTRAINT_field2) {
      ptr -> field2 = SourceQualChoiceAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SOURCE_CONSTRAINT_constraint) {
      ptr -> constraint = StringConstraintAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SOURCE_CONSTRAINT_type_constraint) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> type_constraint = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = SourceConstraintFree(ptr);
   goto ret;
}



/**************************************************
*
*    SourceConstraintAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
SourceConstraintAsnWrite(SourceConstraintPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, SOURCE_CONSTRAINT);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> field1 != NULL) {
      if ( ! SourceQualChoiceAsnWrite(ptr -> field1, aip, SOURCE_CONSTRAINT_field1)) {
         goto erret;
      }
   }
   if (ptr -> field2 != NULL) {
      if ( ! SourceQualChoiceAsnWrite(ptr -> field2, aip, SOURCE_CONSTRAINT_field2)) {
         goto erret;
      }
   }
   if (ptr -> constraint != NULL) {
      if ( ! StringConstraintAsnWrite(ptr -> constraint, aip, SOURCE_CONSTRAINT_constraint)) {
         goto erret;
      }
   }
   av.intvalue = ptr -> type_constraint;
   retval = AsnWrite(aip, SOURCE_CONSTRAINT_type_constraint,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    CDSGeneProtPseudoConstraintNew()
*
**************************************************/
NLM_EXTERN 
CDSGeneProtPseudoConstraintPtr LIBCALL
CDSGeneProtPseudoConstraintNew(void)
{
   CDSGeneProtPseudoConstraintPtr ptr = MemNew((size_t) sizeof(CDSGeneProtPseudoConstraint));

   ptr -> is_pseudo = 1;
   return ptr;

}


/**************************************************
*
*    CDSGeneProtPseudoConstraintFree()
*
**************************************************/
NLM_EXTERN 
CDSGeneProtPseudoConstraintPtr LIBCALL
CDSGeneProtPseudoConstraintFree(CDSGeneProtPseudoConstraintPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   return MemFree(ptr);
}


/**************************************************
*
*    CDSGeneProtPseudoConstraintAsnRead()
*
**************************************************/
NLM_EXTERN 
CDSGeneProtPseudoConstraintPtr LIBCALL
CDSGeneProtPseudoConstraintAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   CDSGeneProtPseudoConstraintPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* CDSGeneProtPseudoConstraint ::= (self contained) */
      atp = AsnReadId(aip, amp, CDSGENEPROT_PSEUDO_CONSTRAINT);
   } else {
      atp = AsnLinkType(orig, CDSGENEPROT_PSEUDO_CONSTRAINT);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = CDSGeneProtPseudoConstraintNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == CDSGENEPROT_PSEUDO_CONSTRAINT_feature) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> feature = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CDSGENEPROT_PSEUDO_CONSTRAINT_is_pseudo) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> is_pseudo = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = CDSGeneProtPseudoConstraintFree(ptr);
   goto ret;
}



/**************************************************
*
*    CDSGeneProtPseudoConstraintAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
CDSGeneProtPseudoConstraintAsnWrite(CDSGeneProtPseudoConstraintPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, CDSGENEPROT_PSEUDO_CONSTRAINT);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> feature;
   retval = AsnWrite(aip, CDSGENEPROT_PSEUDO_CONSTRAINT_feature,  &av);
   av.boolvalue = ptr -> is_pseudo;
   retval = AsnWrite(aip, CDSGENEPROT_PSEUDO_CONSTRAINT_is_pseudo,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    CDSGeneProtConstraintFieldFree()
*
**************************************************/
NLM_EXTERN 
CDSGeneProtConstraintFieldPtr LIBCALL
CDSGeneProtConstraintFieldFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    CDSGeneProtConstraintFieldAsnRead()
*
**************************************************/
NLM_EXTERN 
CDSGeneProtConstraintFieldPtr LIBCALL
CDSGeneProtConstraintFieldAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* CDSGeneProtConstraintField ::= (self contained) */
      atp = AsnReadId(aip, amp, CDSGENEPROT_CONSTRAINT_FIELD);
   } else {
      atp = AsnLinkType(orig, CDSGENEPROT_CONSTRAINT_FIELD);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == CDSGENEPROT_CONSTRAINT_FIELD_field) {
      choice = CDSGeneProtConstraintField_field;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    CDSGeneProtConstraintFieldAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
CDSGeneProtConstraintFieldAsnWrite(CDSGeneProtConstraintFieldPtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, CDSGENEPROT_CONSTRAINT_FIELD);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case CDSGeneProtConstraintField_field:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, CDSGENEPROT_CONSTRAINT_FIELD_field, &av);
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    CDSGeneProtQualConstraintNew()
*
**************************************************/
NLM_EXTERN 
CDSGeneProtQualConstraintPtr LIBCALL
CDSGeneProtQualConstraintNew(void)
{
   CDSGeneProtQualConstraintPtr ptr = MemNew((size_t) sizeof(CDSGeneProtQualConstraint));

   return ptr;

}


/**************************************************
*
*    CDSGeneProtQualConstraintFree()
*
**************************************************/
NLM_EXTERN 
CDSGeneProtQualConstraintPtr LIBCALL
CDSGeneProtQualConstraintFree(CDSGeneProtQualConstraintPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   CDSGeneProtConstraintFieldFree(ptr -> field1);
   CDSGeneProtConstraintFieldFree(ptr -> field2);
   StringConstraintFree(ptr -> constraint);
   return MemFree(ptr);
}


/**************************************************
*
*    CDSGeneProtQualConstraintAsnRead()
*
**************************************************/
NLM_EXTERN 
CDSGeneProtQualConstraintPtr LIBCALL
CDSGeneProtQualConstraintAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   CDSGeneProtQualConstraintPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* CDSGeneProtQualConstraint ::= (self contained) */
      atp = AsnReadId(aip, amp, CDSGENEPROT_QUAL_CONSTRAINT);
   } else {
      atp = AsnLinkType(orig, CDSGENEPROT_QUAL_CONSTRAINT);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = CDSGeneProtQualConstraintNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == CDSGENEPROT_QUAL_CONSTRAINT_field1) {
      ptr -> field1 = CDSGeneProtConstraintFieldAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CDSGENEPROT_QUAL_CONSTRAINT_field2) {
      ptr -> field2 = CDSGeneProtConstraintFieldAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CDSGENEPROT_QUAL_CONSTRAINT_constraint) {
      ptr -> constraint = StringConstraintAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = CDSGeneProtQualConstraintFree(ptr);
   goto ret;
}



/**************************************************
*
*    CDSGeneProtQualConstraintAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
CDSGeneProtQualConstraintAsnWrite(CDSGeneProtQualConstraintPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, CDSGENEPROT_QUAL_CONSTRAINT);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> field1 != NULL) {
      if ( ! CDSGeneProtConstraintFieldAsnWrite(ptr -> field1, aip, CDSGENEPROT_QUAL_CONSTRAINT_field1)) {
         goto erret;
      }
   }
   if (ptr -> field2 != NULL) {
      if ( ! CDSGeneProtConstraintFieldAsnWrite(ptr -> field2, aip, CDSGENEPROT_QUAL_CONSTRAINT_field2)) {
         goto erret;
      }
   }
   if (ptr -> constraint != NULL) {
      if ( ! StringConstraintAsnWrite(ptr -> constraint, aip, CDSGENEPROT_QUAL_CONSTRAINT_constraint)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    FieldConstraintNew()
*
**************************************************/
NLM_EXTERN 
FieldConstraintPtr LIBCALL
FieldConstraintNew(void)
{
   FieldConstraintPtr ptr = MemNew((size_t) sizeof(FieldConstraint));

   return ptr;

}


/**************************************************
*
*    FieldConstraintFree()
*
**************************************************/
NLM_EXTERN 
FieldConstraintPtr LIBCALL
FieldConstraintFree(FieldConstraintPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   FieldTypeFree(ptr -> field);
   StringConstraintFree(ptr -> string_constraint);
   return MemFree(ptr);
}


/**************************************************
*
*    FieldConstraintAsnRead()
*
**************************************************/
NLM_EXTERN 
FieldConstraintPtr LIBCALL
FieldConstraintAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   FieldConstraintPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* FieldConstraint ::= (self contained) */
      atp = AsnReadId(aip, amp, FIELD_CONSTRAINT);
   } else {
      atp = AsnLinkType(orig, FIELD_CONSTRAINT);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = FieldConstraintNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == FIELD_CONSTRAINT_field) {
      ptr -> field = FieldTypeAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == FIELD_CONSTRAINT_string_constraint) {
      ptr -> string_constraint = StringConstraintAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = FieldConstraintFree(ptr);
   goto ret;
}



/**************************************************
*
*    FieldConstraintAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
FieldConstraintAsnWrite(FieldConstraintPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, FIELD_CONSTRAINT);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> field != NULL) {
      if ( ! FieldTypeAsnWrite(ptr -> field, aip, FIELD_CONSTRAINT_field)) {
         goto erret;
      }
   }
   if (ptr -> string_constraint != NULL) {
      if ( ! StringConstraintAsnWrite(ptr -> string_constraint, aip, FIELD_CONSTRAINT_string_constraint)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    FieldTypeFree()
*
**************************************************/
NLM_EXTERN 
FieldTypePtr LIBCALL
FieldTypeFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case FieldType_source_qual:
      SourceQualChoiceFree(anp -> data.ptrvalue);
      break;
   case FieldType_feature_field:
      FeatureFieldFree(anp -> data.ptrvalue);
      break;
   case FieldType_rna_field:
      RnaQualFree(anp -> data.ptrvalue);
      break;
   case FieldType_molinfo_field:
      MolinfoFieldFree(anp -> data.ptrvalue);
      break;
   case FieldType_struc_comment_field:
      StructuredCommentFieldFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    FieldTypeAsnRead()
*
**************************************************/
NLM_EXTERN 
FieldTypePtr LIBCALL
FieldTypeAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* FieldType ::= (self contained) */
      atp = AsnReadId(aip, amp, FIELD_TYPE);
   } else {
      atp = AsnLinkType(orig, FIELD_TYPE);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == FIELD_TYPE_source_qual) {
      choice = FieldType_source_qual;
      func = (AsnReadFunc) SourceQualChoiceAsnRead;
   }
   else if (atp == FIELD_TYPE_feature_field) {
      choice = FieldType_feature_field;
      func = (AsnReadFunc) FeatureFieldAsnRead;
   }
   else if (atp == FIELD_TYPE_rna_field) {
      choice = FieldType_rna_field;
      func = (AsnReadFunc) RnaQualAsnRead;
   }
   else if (atp == FIELD_TYPE_cds_gene_prot) {
      choice = FieldType_cds_gene_prot;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == FIELD_TYPE_molinfo_field) {
      choice = FieldType_molinfo_field;
      func = (AsnReadFunc) MolinfoFieldAsnRead;
   }
   else if (atp == FIELD_TYPE_pub) {
      choice = FieldType_pub;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == FIELD_TYPE_struc_comment_field) {
      choice = FieldType_struc_comment_field;
      func = (AsnReadFunc) StructuredCommentFieldAsnRead;
   }
   else if (atp == FIELD_TYPE_misc) {
      choice = FieldType_misc;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == FIELD_TYPE_dblink) {
      choice = FieldType_dblink;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    FieldTypeAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
FieldTypeAsnWrite(FieldTypePtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, FIELD_TYPE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case FieldType_source_qual:
      writetype = FIELD_TYPE_source_qual;
      func = (AsnWriteFunc) SourceQualChoiceAsnWrite;
      break;
   case FieldType_feature_field:
      writetype = FIELD_TYPE_feature_field;
      func = (AsnWriteFunc) FeatureFieldAsnWrite;
      break;
   case FieldType_rna_field:
      writetype = FIELD_TYPE_rna_field;
      func = (AsnWriteFunc) RnaQualAsnWrite;
      break;
   case FieldType_cds_gene_prot:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, FIELD_TYPE_cds_gene_prot, &av);
      break;
   case FieldType_molinfo_field:
      writetype = FIELD_TYPE_molinfo_field;
      func = (AsnWriteFunc) MolinfoFieldAsnWrite;
      break;
   case FieldType_pub:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, FIELD_TYPE_pub, &av);
      break;
   case FieldType_struc_comment_field:
      writetype = FIELD_TYPE_struc_comment_field;
      func = (AsnWriteFunc) StructuredCommentFieldAsnWrite;
      break;
   case FieldType_misc:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, FIELD_TYPE_misc, &av);
      break;
   case FieldType_dblink:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, FIELD_TYPE_dblink, &av);
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    SequenceConstraintMolTypeConstraintFree()
*
**************************************************/
NLM_EXTERN 
SequenceConstraintMolTypeConstraintPtr LIBCALL
SequenceConstraintMolTypeConstraintFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    SequenceConstraintMolTypeConstraintAsnRead()
*
**************************************************/
NLM_EXTERN 
SequenceConstraintMolTypeConstraintPtr LIBCALL
SequenceConstraintMolTypeConstraintAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* SequenceConstraintMolTypeConstraint ::= (self contained) */
      atp = AsnReadId(aip, amp, SEQUENCE_CONSTRAINT_MOL_TYPE_CONSTRAINT);
   } else {
      atp = AsnLinkType(orig, SEQUENCE_CONSTRAINT_MOL_TYPE_CONSTRAINT);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == SEQUENCE_CONSTRAINT_MOL_TYPE_CONSTRAINT_any) {
      choice = SequenceConstraintMolTypeConstraint_any;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == SEQUENCE_CONSTRAINT_MOL_TYPE_CONSTRAINT_nucleotide) {
      choice = SequenceConstraintMolTypeConstraint_nucleotide;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == SEQUENCE_CONSTRAINT_MOL_TYPE_CONSTRAINT_dna) {
      choice = SequenceConstraintMolTypeConstraint_dna;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == SEQUENCE_CONSTRAINT_MOL_TYPE_CONSTRAINT_rna) {
      choice = SequenceConstraintMolTypeConstraint_rna;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == SEQUENCE_CONSTRAINT_MOL_TYPE_CONSTRAINT_protein) {
      choice = SequenceConstraintMolTypeConstraint_protein;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    SequenceConstraintMolTypeConstraintAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
SequenceConstraintMolTypeConstraintAsnWrite(SequenceConstraintMolTypeConstraintPtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, SEQUENCE_CONSTRAINT_MOL_TYPE_CONSTRAINT);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case SequenceConstraintMolTypeConstraint_any:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, SEQUENCE_CONSTRAINT_MOL_TYPE_CONSTRAINT_any, &av);
      break;
   case SequenceConstraintMolTypeConstraint_nucleotide:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, SEQUENCE_CONSTRAINT_MOL_TYPE_CONSTRAINT_nucleotide, &av);
      break;
   case SequenceConstraintMolTypeConstraint_dna:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, SEQUENCE_CONSTRAINT_MOL_TYPE_CONSTRAINT_dna, &av);
      break;
   case SequenceConstraintMolTypeConstraint_rna:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, SEQUENCE_CONSTRAINT_MOL_TYPE_CONSTRAINT_rna, &av);
      break;
   case SequenceConstraintMolTypeConstraint_protein:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, SEQUENCE_CONSTRAINT_MOL_TYPE_CONSTRAINT_protein, &av);
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    QuantityConstraintFree()
*
**************************************************/
NLM_EXTERN 
QuantityConstraintPtr LIBCALL
QuantityConstraintFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    QuantityConstraintAsnRead()
*
**************************************************/
NLM_EXTERN 
QuantityConstraintPtr LIBCALL
QuantityConstraintAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* QuantityConstraint ::= (self contained) */
      atp = AsnReadId(aip, amp, QUANTITY_CONSTRAINT);
   } else {
      atp = AsnLinkType(orig, QUANTITY_CONSTRAINT);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == QUANTITY_CONSTRAINT_equals) {
      choice = QuantityConstraint_equals;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == QUANTITY_CONSTRAINT_greater_than) {
      choice = QuantityConstraint_greater_than;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == QUANTITY_CONSTRAINT_less_than) {
      choice = QuantityConstraint_less_than;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    QuantityConstraintAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
QuantityConstraintAsnWrite(QuantityConstraintPtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, QUANTITY_CONSTRAINT);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case QuantityConstraint_equals:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, QUANTITY_CONSTRAINT_equals, &av);
      break;
   case QuantityConstraint_greater_than:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, QUANTITY_CONSTRAINT_greater_than, &av);
      break;
   case QuantityConstraint_less_than:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, QUANTITY_CONSTRAINT_less_than, &av);
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    SequenceConstraintNew()
*
**************************************************/
NLM_EXTERN 
SequenceConstraintPtr LIBCALL
SequenceConstraintNew(void)
{
   SequenceConstraintPtr ptr = MemNew((size_t) sizeof(SequenceConstraint));

   ptr -> strandedness = 0;
   return ptr;

}


/**************************************************
*
*    SequenceConstraintFree()
*
**************************************************/
NLM_EXTERN 
SequenceConstraintPtr LIBCALL
SequenceConstraintFree(SequenceConstraintPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   SequenceConstraintMolTypeConstraintFree(ptr -> seqtype);
   StringConstraintFree(ptr -> id);
   QuantityConstraintFree(ptr -> num_type_features);
   QuantityConstraintFree(ptr -> num_features);
   QuantityConstraintFree(ptr -> length);
   return MemFree(ptr);
}


/**************************************************
*
*    SequenceConstraintAsnRead()
*
**************************************************/
NLM_EXTERN 
SequenceConstraintPtr LIBCALL
SequenceConstraintAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   SequenceConstraintPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* SequenceConstraint ::= (self contained) */
      atp = AsnReadId(aip, amp, SEQUENCE_CONSTRAINT);
   } else {
      atp = AsnLinkType(orig, SEQUENCE_CONSTRAINT);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = SequenceConstraintNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == SEQUENCE_CONSTRAINT_seqtype) {
      ptr -> seqtype = SequenceConstraintMolTypeConstraintAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SEQUENCE_CONSTRAINT_id) {
      ptr -> id = StringConstraintAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SEQUENCE_CONSTRAINT_feature) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> feature = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SEQUENCE_CONSTRAINT_num_type_features) {
      ptr -> num_type_features = QuantityConstraintAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SEQUENCE_CONSTRAINT_num_features) {
      ptr -> num_features = QuantityConstraintAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SEQUENCE_CONSTRAINT_length) {
      ptr -> length = QuantityConstraintAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SEQUENCE_CONSTRAINT_strandedness) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> strandedness = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = SequenceConstraintFree(ptr);
   goto ret;
}



/**************************************************
*
*    SequenceConstraintAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
SequenceConstraintAsnWrite(SequenceConstraintPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, SEQUENCE_CONSTRAINT);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> seqtype != NULL) {
      if ( ! SequenceConstraintMolTypeConstraintAsnWrite(ptr -> seqtype, aip, SEQUENCE_CONSTRAINT_seqtype)) {
         goto erret;
      }
   }
   if (ptr -> id != NULL) {
      if ( ! StringConstraintAsnWrite(ptr -> id, aip, SEQUENCE_CONSTRAINT_id)) {
         goto erret;
      }
   }
   av.intvalue = ptr -> feature;
   retval = AsnWrite(aip, SEQUENCE_CONSTRAINT_feature,  &av);
   if (ptr -> num_type_features != NULL) {
      if ( ! QuantityConstraintAsnWrite(ptr -> num_type_features, aip, SEQUENCE_CONSTRAINT_num_type_features)) {
         goto erret;
      }
   }
   if (ptr -> num_features != NULL) {
      if ( ! QuantityConstraintAsnWrite(ptr -> num_features, aip, SEQUENCE_CONSTRAINT_num_features)) {
         goto erret;
      }
   }
   if (ptr -> length != NULL) {
      if ( ! QuantityConstraintAsnWrite(ptr -> length, aip, SEQUENCE_CONSTRAINT_length)) {
         goto erret;
      }
   }
   av.intvalue = ptr -> strandedness;
   retval = AsnWrite(aip, SEQUENCE_CONSTRAINT_strandedness,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    TranslationConstraintNew()
*
**************************************************/
NLM_EXTERN 
TranslationConstraintPtr LIBCALL
TranslationConstraintNew(void)
{
   TranslationConstraintPtr ptr = MemNew((size_t) sizeof(TranslationConstraint));

   ptr -> internal_stops = 0;
   return ptr;

}


/**************************************************
*
*    TranslationConstraintFree()
*
**************************************************/
NLM_EXTERN 
TranslationConstraintPtr LIBCALL
TranslationConstraintFree(TranslationConstraintPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   StringConstraintSetFree(ptr -> actual_strings);
   StringConstraintSetFree(ptr -> transl_strings);
   QuantityConstraintFree(ptr -> num_mismatches);
   return MemFree(ptr);
}


/**************************************************
*
*    TranslationConstraintAsnRead()
*
**************************************************/
NLM_EXTERN 
TranslationConstraintPtr LIBCALL
TranslationConstraintAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   TranslationConstraintPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* TranslationConstraint ::= (self contained) */
      atp = AsnReadId(aip, amp, TRANSLATION_CONSTRAINT);
   } else {
      atp = AsnLinkType(orig, TRANSLATION_CONSTRAINT);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = TranslationConstraintNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == TRANSLATION_CONSTRAINT_actual_strings) {
      ptr -> actual_strings = StringConstraintSetAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TRANSLATION_CONSTRAINT_transl_strings) {
      ptr -> transl_strings = StringConstraintSetAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TRANSLATION_CONSTRAINT_internal_stops) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> internal_stops = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TRANSLATION_CONSTRAINT_num_mismatches) {
      ptr -> num_mismatches = QuantityConstraintAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = TranslationConstraintFree(ptr);
   goto ret;
}



/**************************************************
*
*    TranslationConstraintAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
TranslationConstraintAsnWrite(TranslationConstraintPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, TRANSLATION_CONSTRAINT);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> actual_strings != NULL) {
      if ( ! StringConstraintSetAsnWrite(ptr -> actual_strings, aip, TRANSLATION_CONSTRAINT_actual_strings)) {
         goto erret;
      }
   }
   if (ptr -> transl_strings != NULL) {
      if ( ! StringConstraintSetAsnWrite(ptr -> transl_strings, aip, TRANSLATION_CONSTRAINT_transl_strings)) {
         goto erret;
      }
   }
   av.intvalue = ptr -> internal_stops;
   retval = AsnWrite(aip, TRANSLATION_CONSTRAINT_internal_stops,  &av);
   if (ptr -> num_mismatches != NULL) {
      if ( ! QuantityConstraintAsnWrite(ptr -> num_mismatches, aip, TRANSLATION_CONSTRAINT_num_mismatches)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ConstraintChoiceFree()
*
**************************************************/
NLM_EXTERN 
ConstraintChoicePtr LIBCALL
ConstraintChoiceFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case ConstraintChoice_string:
      StringConstraintFree(anp -> data.ptrvalue);
      break;
   case ConstraintChoice_location:
      LocationConstraintFree(anp -> data.ptrvalue);
      break;
   case ConstraintChoice_field:
      FieldConstraintFree(anp -> data.ptrvalue);
      break;
   case ConstraintChoice_source:
      SourceConstraintFree(anp -> data.ptrvalue);
      break;
   case ConstraintChoice_cdsgeneprot_qual:
      CDSGeneProtQualConstraintFree(anp -> data.ptrvalue);
      break;
   case ConstraintChoice_cdsgeneprot_pseudo:
      CDSGeneProtPseudoConstraintFree(anp -> data.ptrvalue);
      break;
   case ConstraintChoice_sequence:
      SequenceConstraintFree(anp -> data.ptrvalue);
      break;
   case ConstraintChoice_pub:
      PublicationConstraintFree(anp -> data.ptrvalue);
      break;
   case ConstraintChoice_molinfo:
      MolinfoFieldConstraintFree(anp -> data.ptrvalue);
      break;
   case ConstraintChoice_field_missing:
      FieldTypeFree(anp -> data.ptrvalue);
      break;
   case ConstraintChoice_translation:
      TranslationConstraintFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    ConstraintChoiceAsnRead()
*
**************************************************/
NLM_EXTERN 
ConstraintChoicePtr LIBCALL
ConstraintChoiceAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ConstraintChoice ::= (self contained) */
      atp = AsnReadId(aip, amp, CONSTRAINT_CHOICE);
   } else {
      atp = AsnLinkType(orig, CONSTRAINT_CHOICE);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == CONSTRAINT_CHOICE_string) {
      choice = ConstraintChoice_string;
      func = (AsnReadFunc) StringConstraintAsnRead;
   }
   else if (atp == CONSTRAINT_CHOICE_location) {
      choice = ConstraintChoice_location;
      func = (AsnReadFunc) LocationConstraintAsnRead;
   }
   else if (atp == CONSTRAINT_CHOICE_field) {
      choice = ConstraintChoice_field;
      func = (AsnReadFunc) FieldConstraintAsnRead;
   }
   else if (atp == CONSTRAINT_CHOICE_source) {
      choice = ConstraintChoice_source;
      func = (AsnReadFunc) SourceConstraintAsnRead;
   }
   else if (atp == CONSTRAINT_CHOICE_cdsgeneprot_qual) {
      choice = ConstraintChoice_cdsgeneprot_qual;
      func = (AsnReadFunc) CDSGeneProtQualConstraintAsnRead;
   }
   else if (atp == CONSTRAINT_CHOICE_cdsgeneprot_pseudo) {
      choice = ConstraintChoice_cdsgeneprot_pseudo;
      func = (AsnReadFunc) CDSGeneProtPseudoConstraintAsnRead;
   }
   else if (atp == CONSTRAINT_CHOICE_sequence) {
      choice = ConstraintChoice_sequence;
      func = (AsnReadFunc) SequenceConstraintAsnRead;
   }
   else if (atp == CONSTRAINT_CHOICE_pub) {
      choice = ConstraintChoice_pub;
      func = (AsnReadFunc) PublicationConstraintAsnRead;
   }
   else if (atp == CONSTRAINT_CHOICE_molinfo) {
      choice = ConstraintChoice_molinfo;
      func = (AsnReadFunc) MolinfoFieldConstraintAsnRead;
   }
   else if (atp == CONSTRAINT_CHOICE_field_missing) {
      choice = ConstraintChoice_field_missing;
      func = (AsnReadFunc) FieldTypeAsnRead;
   }
   else if (atp == CONSTRAINT_CHOICE_translation) {
      choice = ConstraintChoice_translation;
      func = (AsnReadFunc) TranslationConstraintAsnRead;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    ConstraintChoiceAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ConstraintChoiceAsnWrite(ConstraintChoicePtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, CONSTRAINT_CHOICE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case ConstraintChoice_string:
      writetype = CONSTRAINT_CHOICE_string;
      func = (AsnWriteFunc) StringConstraintAsnWrite;
      break;
   case ConstraintChoice_location:
      writetype = CONSTRAINT_CHOICE_location;
      func = (AsnWriteFunc) LocationConstraintAsnWrite;
      break;
   case ConstraintChoice_field:
      writetype = CONSTRAINT_CHOICE_field;
      func = (AsnWriteFunc) FieldConstraintAsnWrite;
      break;
   case ConstraintChoice_source:
      writetype = CONSTRAINT_CHOICE_source;
      func = (AsnWriteFunc) SourceConstraintAsnWrite;
      break;
   case ConstraintChoice_cdsgeneprot_qual:
      writetype = CONSTRAINT_CHOICE_cdsgeneprot_qual;
      func = (AsnWriteFunc) CDSGeneProtQualConstraintAsnWrite;
      break;
   case ConstraintChoice_cdsgeneprot_pseudo:
      writetype = CONSTRAINT_CHOICE_cdsgeneprot_pseudo;
      func = (AsnWriteFunc) CDSGeneProtPseudoConstraintAsnWrite;
      break;
   case ConstraintChoice_sequence:
      writetype = CONSTRAINT_CHOICE_sequence;
      func = (AsnWriteFunc) SequenceConstraintAsnWrite;
      break;
   case ConstraintChoice_pub:
      writetype = CONSTRAINT_CHOICE_pub;
      func = (AsnWriteFunc) PublicationConstraintAsnWrite;
      break;
   case ConstraintChoice_molinfo:
      writetype = CONSTRAINT_CHOICE_molinfo;
      func = (AsnWriteFunc) MolinfoFieldConstraintAsnWrite;
      break;
   case ConstraintChoice_field_missing:
      writetype = CONSTRAINT_CHOICE_field_missing;
      func = (AsnWriteFunc) FieldTypeAsnWrite;
      break;
   case ConstraintChoice_translation:
      writetype = CONSTRAINT_CHOICE_translation;
      func = (AsnWriteFunc) TranslationConstraintAsnWrite;
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    ConstraintChoiceSetFree()
*
**************************************************/
NLM_EXTERN 
ConstraintChoiceSetPtr LIBCALL
ConstraintChoiceSetFree(ConstraintChoiceSetPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericChoiceSeqOfFree(ptr, (AsnOptFreeFunc) ConstraintChoiceFree);
   return NULL;
}


/**************************************************
*
*    ConstraintChoiceSetAsnRead()
*
**************************************************/
NLM_EXTERN 
ConstraintChoiceSetPtr LIBCALL
ConstraintChoiceSetAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ConstraintChoiceSetPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ConstraintChoiceSet ::= (self contained) */
      atp = AsnReadId(aip, amp, CONSTRAINT_CHOICE_SET);
   } else {
      atp = AsnLinkType(orig, CONSTRAINT_CHOICE_SET);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   func = NULL;

   ptr  = AsnGenericChoiceSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) ConstraintChoiceAsnRead, (AsnOptFreeFunc) ConstraintChoiceFree);
   if (isError && ptr  == NULL) {
      goto erret;
   }



ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ConstraintChoiceSetFree(ptr);
   goto ret;
}



/**************************************************
*
*    ConstraintChoiceSetAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ConstraintChoiceSetAsnWrite(ConstraintChoiceSetPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, CONSTRAINT_CHOICE_SET);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   retval = AsnGenericChoiceSeqOfAsnWrite(ptr , (AsnWriteFunc) ConstraintChoiceAsnWrite, aip, atp, CONSTRAINT_CHOICE_SET_E);
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    TextMarkerFree()
*
**************************************************/
NLM_EXTERN 
TextMarkerPtr LIBCALL
TextMarkerFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case TextMarker_free_text:
      MemFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    TextMarkerAsnRead()
*
**************************************************/
NLM_EXTERN 
TextMarkerPtr LIBCALL
TextMarkerAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* TextMarker ::= (self contained) */
      atp = AsnReadId(aip, amp, TEXT_MARKER);
   } else {
      atp = AsnLinkType(orig, TEXT_MARKER);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == TEXT_MARKER_free_text) {
      choice = TextMarker_free_text;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.ptrvalue = av.ptrvalue;
   }
   else if (atp == TEXT_MARKER_digits) {
      choice = TextMarker_digits;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == TEXT_MARKER_letters) {
      choice = TextMarker_letters;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    TextMarkerAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
TextMarkerAsnWrite(TextMarkerPtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, TEXT_MARKER);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case TextMarker_free_text:
      av.ptrvalue = anp->data.ptrvalue;
      retval = AsnWrite(aip, TEXT_MARKER_free_text, &av);
      break;
   case TextMarker_digits:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, TEXT_MARKER_digits, &av);
      break;
   case TextMarker_letters:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, TEXT_MARKER_letters, &av);
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    TextPortionNew()
*
**************************************************/
NLM_EXTERN 
TextPortionPtr LIBCALL
TextPortionNew(void)
{
   TextPortionPtr ptr = MemNew((size_t) sizeof(TextPortion));

   ptr -> case_sensitive = 0;
   ptr -> whole_word = 0;
   return ptr;

}


/**************************************************
*
*    TextPortionFree()
*
**************************************************/
NLM_EXTERN 
TextPortionPtr LIBCALL
TextPortionFree(TextPortionPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   TextMarkerFree(ptr -> left_marker);
   TextMarkerFree(ptr -> right_marker);
   return MemFree(ptr);
}


/**************************************************
*
*    TextPortionAsnRead()
*
**************************************************/
NLM_EXTERN 
TextPortionPtr LIBCALL
TextPortionAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   TextPortionPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* TextPortion ::= (self contained) */
      atp = AsnReadId(aip, amp, TEXT_PORTION);
   } else {
      atp = AsnLinkType(orig, TEXT_PORTION);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = TextPortionNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == TEXT_PORTION_left_marker) {
      ptr -> left_marker = TextMarkerAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TEXT_PORTION_include_left) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> include_left = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TEXT_PORTION_right_marker) {
      ptr -> right_marker = TextMarkerAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TEXT_PORTION_include_right) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> include_right = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TEXT_PORTION_inside) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> inside = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TEXT_PORTION_case_sensitive) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> case_sensitive = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TEXT_PORTION_whole_word) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> whole_word = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = TextPortionFree(ptr);
   goto ret;
}



/**************************************************
*
*    TextPortionAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
TextPortionAsnWrite(TextPortionPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, TEXT_PORTION);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> left_marker != NULL) {
      if ( ! TextMarkerAsnWrite(ptr -> left_marker, aip, TEXT_PORTION_left_marker)) {
         goto erret;
      }
   }
   av.boolvalue = ptr -> include_left;
   retval = AsnWrite(aip, TEXT_PORTION_include_left,  &av);
   if (ptr -> right_marker != NULL) {
      if ( ! TextMarkerAsnWrite(ptr -> right_marker, aip, TEXT_PORTION_right_marker)) {
         goto erret;
      }
   }
   av.boolvalue = ptr -> include_right;
   retval = AsnWrite(aip, TEXT_PORTION_include_right,  &av);
   av.boolvalue = ptr -> inside;
   retval = AsnWrite(aip, TEXT_PORTION_inside,  &av);
   av.boolvalue = ptr -> case_sensitive;
   retval = AsnWrite(aip, TEXT_PORTION_case_sensitive,  &av);
   av.boolvalue = ptr -> whole_word;
   retval = AsnWrite(aip, TEXT_PORTION_whole_word,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    FieldEditNew()
*
**************************************************/
NLM_EXTERN 
FieldEditPtr LIBCALL
FieldEditNew(void)
{
   FieldEditPtr ptr = MemNew((size_t) sizeof(FieldEdit));

   ptr -> location = 0;
   ptr -> case_insensitive = 0;
   return ptr;

}


/**************************************************
*
*    FieldEditFree()
*
**************************************************/
NLM_EXTERN 
FieldEditPtr LIBCALL
FieldEditFree(FieldEditPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> find_txt);
   MemFree(ptr -> repl_txt);
   return MemFree(ptr);
}


/**************************************************
*
*    FieldEditAsnRead()
*
**************************************************/
NLM_EXTERN 
FieldEditPtr LIBCALL
FieldEditAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   FieldEditPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* FieldEdit ::= (self contained) */
      atp = AsnReadId(aip, amp, FIELD_EDIT);
   } else {
      atp = AsnLinkType(orig, FIELD_EDIT);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = FieldEditNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == FIELD_EDIT_find_txt) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> find_txt = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == FIELD_EDIT_repl_txt) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> repl_txt = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == FIELD_EDIT_location) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> location = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == FIELD_EDIT_case_insensitive) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> case_insensitive = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = FieldEditFree(ptr);
   goto ret;
}



/**************************************************
*
*    FieldEditAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
FieldEditAsnWrite(FieldEditPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, FIELD_EDIT);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> find_txt != NULL) {
      av.ptrvalue = ptr -> find_txt;
      retval = AsnWrite(aip, FIELD_EDIT_find_txt,  &av);
   }
   if (ptr -> repl_txt != NULL) {
      av.ptrvalue = ptr -> repl_txt;
      retval = AsnWrite(aip, FIELD_EDIT_repl_txt,  &av);
   }
   av.intvalue = ptr -> location;
   retval = AsnWrite(aip, FIELD_EDIT_location,  &av);
   av.boolvalue = ptr -> case_insensitive;
   retval = AsnWrite(aip, FIELD_EDIT_case_insensitive,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    FieldPairTypeFree()
*
**************************************************/
NLM_EXTERN 
FieldPairTypePtr LIBCALL
FieldPairTypeFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case FieldPairType_source_qual:
      SourceQualPairFree(anp -> data.ptrvalue);
      break;
   case FieldPairType_feature_field:
      FeatureFieldPairFree(anp -> data.ptrvalue);
      break;
   case FieldPairType_rna_field:
      RnaQualPairFree(anp -> data.ptrvalue);
      break;
   case FieldPairType_cds_gene_prot:
      CDSGeneProtFieldPairFree(anp -> data.ptrvalue);
      break;
   case FieldPairType_molinfo_field:
      MolinfoFieldPairFree(anp -> data.ptrvalue);
      break;
   case FieldPairType_struc_comment_field:
      StructuredCommentFieldPairFree(anp -> data.ptrvalue);
      break;
   case FieldPairType_dblink:
      DBLinkFieldPairFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    FieldPairTypeAsnRead()
*
**************************************************/
NLM_EXTERN 
FieldPairTypePtr LIBCALL
FieldPairTypeAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* FieldPairType ::= (self contained) */
      atp = AsnReadId(aip, amp, FIELD_PAIR_TYPE);
   } else {
      atp = AsnLinkType(orig, FIELD_PAIR_TYPE);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == FIELD_PAIR_TYPE_source_qual) {
      choice = FieldPairType_source_qual;
      func = (AsnReadFunc) SourceQualPairAsnRead;
   }
   else if (atp == FIELD_PAIR_TYPE_feature_field) {
      choice = FieldPairType_feature_field;
      func = (AsnReadFunc) FeatureFieldPairAsnRead;
   }
   else if (atp == FIELD_PAIR_TYPE_rna_field) {
      choice = FieldPairType_rna_field;
      func = (AsnReadFunc) RnaQualPairAsnRead;
   }
   else if (atp == FIELD_PAIR_TYPE_cds_gene_prot) {
      choice = FieldPairType_cds_gene_prot;
      func = (AsnReadFunc) CDSGeneProtFieldPairAsnRead;
   }
   else if (atp == FIELD_PAIR_TYPE_molinfo_field) {
      choice = FieldPairType_molinfo_field;
      func = (AsnReadFunc) MolinfoFieldPairAsnRead;
   }
   else if (atp == FIELD_PAIR_TYPE_struc_comment_field) {
      choice = FieldPairType_struc_comment_field;
      func = (AsnReadFunc) StructuredCommentFieldPairAsnRead;
   }
   else if (atp == FIELD_PAIR_TYPE_dblink) {
      choice = FieldPairType_dblink;
      func = (AsnReadFunc) DBLinkFieldPairAsnRead;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    FieldPairTypeAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
FieldPairTypeAsnWrite(FieldPairTypePtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, FIELD_PAIR_TYPE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case FieldPairType_source_qual:
      writetype = FIELD_PAIR_TYPE_source_qual;
      func = (AsnWriteFunc) SourceQualPairAsnWrite;
      break;
   case FieldPairType_feature_field:
      writetype = FIELD_PAIR_TYPE_feature_field;
      func = (AsnWriteFunc) FeatureFieldPairAsnWrite;
      break;
   case FieldPairType_rna_field:
      writetype = FIELD_PAIR_TYPE_rna_field;
      func = (AsnWriteFunc) RnaQualPairAsnWrite;
      break;
   case FieldPairType_cds_gene_prot:
      writetype = FIELD_PAIR_TYPE_cds_gene_prot;
      func = (AsnWriteFunc) CDSGeneProtFieldPairAsnWrite;
      break;
   case FieldPairType_molinfo_field:
      writetype = FIELD_PAIR_TYPE_molinfo_field;
      func = (AsnWriteFunc) MolinfoFieldPairAsnWrite;
      break;
   case FieldPairType_struc_comment_field:
      writetype = FIELD_PAIR_TYPE_struc_comment_field;
      func = (AsnWriteFunc) StructuredCommentFieldPairAsnWrite;
      break;
   case FieldPairType_dblink:
      writetype = FIELD_PAIR_TYPE_dblink;
      func = (AsnWriteFunc) DBLinkFieldPairAsnWrite;
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    ApplyActionNew()
*
**************************************************/
NLM_EXTERN 
ApplyActionPtr LIBCALL
ApplyActionNew(void)
{
   ApplyActionPtr ptr = MemNew((size_t) sizeof(ApplyAction));

   return ptr;

}


/**************************************************
*
*    ApplyActionFree()
*
**************************************************/
NLM_EXTERN 
ApplyActionPtr LIBCALL
ApplyActionFree(ApplyActionPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   FieldTypeFree(ptr -> field);
   MemFree(ptr -> value);
   return MemFree(ptr);
}


/**************************************************
*
*    ApplyActionAsnRead()
*
**************************************************/
NLM_EXTERN 
ApplyActionPtr LIBCALL
ApplyActionAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ApplyActionPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ApplyAction ::= (self contained) */
      atp = AsnReadId(aip, amp, APPLY_ACTION);
   } else {
      atp = AsnLinkType(orig, APPLY_ACTION);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ApplyActionNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == APPLY_ACTION_field) {
      ptr -> field = FieldTypeAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == APPLY_ACTION_value) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> value = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == APPLY_ACTION_existing_text) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> existing_text = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ApplyActionFree(ptr);
   goto ret;
}



/**************************************************
*
*    ApplyActionAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ApplyActionAsnWrite(ApplyActionPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, APPLY_ACTION);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> field != NULL) {
      if ( ! FieldTypeAsnWrite(ptr -> field, aip, APPLY_ACTION_field)) {
         goto erret;
      }
   }
   if (ptr -> value != NULL) {
      av.ptrvalue = ptr -> value;
      retval = AsnWrite(aip, APPLY_ACTION_value,  &av);
   }
   av.intvalue = ptr -> existing_text;
   retval = AsnWrite(aip, APPLY_ACTION_existing_text,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    EditActionNew()
*
**************************************************/
NLM_EXTERN 
EditActionPtr LIBCALL
EditActionNew(void)
{
   EditActionPtr ptr = MemNew((size_t) sizeof(EditAction));

   return ptr;

}


/**************************************************
*
*    EditActionFree()
*
**************************************************/
NLM_EXTERN 
EditActionPtr LIBCALL
EditActionFree(EditActionPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   FieldEditFree(ptr -> edit);
   FieldTypeFree(ptr -> field);
   return MemFree(ptr);
}


/**************************************************
*
*    EditActionAsnRead()
*
**************************************************/
NLM_EXTERN 
EditActionPtr LIBCALL
EditActionAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   EditActionPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* EditAction ::= (self contained) */
      atp = AsnReadId(aip, amp, EDIT_ACTION);
   } else {
      atp = AsnLinkType(orig, EDIT_ACTION);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = EditActionNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == EDIT_ACTION_edit) {
      ptr -> edit = FieldEditAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == EDIT_ACTION_field) {
      ptr -> field = FieldTypeAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = EditActionFree(ptr);
   goto ret;
}



/**************************************************
*
*    EditActionAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
EditActionAsnWrite(EditActionPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, EDIT_ACTION);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> edit != NULL) {
      if ( ! FieldEditAsnWrite(ptr -> edit, aip, EDIT_ACTION_edit)) {
         goto erret;
      }
   }
   if (ptr -> field != NULL) {
      if ( ! FieldTypeAsnWrite(ptr -> field, aip, EDIT_ACTION_field)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    TextTransformFree()
*
**************************************************/
NLM_EXTERN 
TextTransformPtr LIBCALL
TextTransformFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case TextTransform_edit:
      FieldEditFree(anp -> data.ptrvalue);
      break;
   case TextTransform_remove:
      TextPortionFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    TextTransformAsnRead()
*
**************************************************/
NLM_EXTERN 
TextTransformPtr LIBCALL
TextTransformAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* TextTransform ::= (self contained) */
      atp = AsnReadId(aip, amp, TEXT_TRANSFORM);
   } else {
      atp = AsnLinkType(orig, TEXT_TRANSFORM);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == TEXT_TRANSFORM_edit) {
      choice = TextTransform_edit;
      func = (AsnReadFunc) FieldEditAsnRead;
   }
   else if (atp == TEXT_TRANSFORM_caps) {
      choice = TextTransform_caps;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == TEXT_TRANSFORM_remove) {
      choice = TextTransform_remove;
      func = (AsnReadFunc) TextPortionAsnRead;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    TextTransformAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
TextTransformAsnWrite(TextTransformPtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, TEXT_TRANSFORM);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case TextTransform_edit:
      writetype = TEXT_TRANSFORM_edit;
      func = (AsnWriteFunc) FieldEditAsnWrite;
      break;
   case TextTransform_caps:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, TEXT_TRANSFORM_caps, &av);
      break;
   case TextTransform_remove:
      writetype = TEXT_TRANSFORM_remove;
      func = (AsnWriteFunc) TextPortionAsnWrite;
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    TextTransformSetFree()
*
**************************************************/
NLM_EXTERN 
TextTransformSetPtr LIBCALL
TextTransformSetFree(TextTransformSetPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericChoiceSeqOfFree(ptr, (AsnOptFreeFunc) TextTransformFree);
   return NULL;
}


/**************************************************
*
*    TextTransformSetAsnRead()
*
**************************************************/
NLM_EXTERN 
TextTransformSetPtr LIBCALL
TextTransformSetAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   TextTransformSetPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* TextTransformSet ::= (self contained) */
      atp = AsnReadId(aip, amp, TEXT_TRANSFORM_SET);
   } else {
      atp = AsnLinkType(orig, TEXT_TRANSFORM_SET);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   func = NULL;

   ptr  = AsnGenericChoiceSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) TextTransformAsnRead, (AsnOptFreeFunc) TextTransformFree);
   if (isError && ptr  == NULL) {
      goto erret;
   }



ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = TextTransformSetFree(ptr);
   goto ret;
}



/**************************************************
*
*    TextTransformSetAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
TextTransformSetAsnWrite(TextTransformSetPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, TEXT_TRANSFORM_SET);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   retval = AsnGenericChoiceSeqOfAsnWrite(ptr , (AsnWriteFunc) TextTransformAsnWrite, aip, atp, TEXT_TRANSFORM_SET_E);
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ConvertActionNew()
*
**************************************************/
NLM_EXTERN 
ConvertActionPtr LIBCALL
ConvertActionNew(void)
{
   ConvertActionPtr ptr = MemNew((size_t) sizeof(ConvertAction));

   ptr -> strip_name = 0;
   ptr -> keep_original = 0;
   ptr -> capitalization = 0;
   return ptr;

}


/**************************************************
*
*    ConvertActionFree()
*
**************************************************/
NLM_EXTERN 
ConvertActionPtr LIBCALL
ConvertActionFree(ConvertActionPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   FieldPairTypeFree(ptr -> fields);
   return MemFree(ptr);
}


/**************************************************
*
*    ConvertActionAsnRead()
*
**************************************************/
NLM_EXTERN 
ConvertActionPtr LIBCALL
ConvertActionAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ConvertActionPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ConvertAction ::= (self contained) */
      atp = AsnReadId(aip, amp, CONVERT_ACTION);
   } else {
      atp = AsnLinkType(orig, CONVERT_ACTION);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ConvertActionNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == CONVERT_ACTION_fields) {
      ptr -> fields = FieldPairTypeAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CONVERT_ACTION_strip_name) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> strip_name = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CONVERT_ACTION_keep_original) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> keep_original = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CONVERT_ACTION_capitalization) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> capitalization = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CONVERT_ACTION_existing_text) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> existing_text = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ConvertActionFree(ptr);
   goto ret;
}



/**************************************************
*
*    ConvertActionAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ConvertActionAsnWrite(ConvertActionPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, CONVERT_ACTION);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> fields != NULL) {
      if ( ! FieldPairTypeAsnWrite(ptr -> fields, aip, CONVERT_ACTION_fields)) {
         goto erret;
      }
   }
   av.boolvalue = ptr -> strip_name;
   retval = AsnWrite(aip, CONVERT_ACTION_strip_name,  &av);
   av.boolvalue = ptr -> keep_original;
   retval = AsnWrite(aip, CONVERT_ACTION_keep_original,  &av);
   av.intvalue = ptr -> capitalization;
   retval = AsnWrite(aip, CONVERT_ACTION_capitalization,  &av);
   av.intvalue = ptr -> existing_text;
   retval = AsnWrite(aip, CONVERT_ACTION_existing_text,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    CopyActionNew()
*
**************************************************/
NLM_EXTERN 
CopyActionPtr LIBCALL
CopyActionNew(void)
{
   CopyActionPtr ptr = MemNew((size_t) sizeof(CopyAction));

   return ptr;

}


/**************************************************
*
*    CopyActionFree()
*
**************************************************/
NLM_EXTERN 
CopyActionPtr LIBCALL
CopyActionFree(CopyActionPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   FieldPairTypeFree(ptr -> fields);
   return MemFree(ptr);
}


/**************************************************
*
*    CopyActionAsnRead()
*
**************************************************/
NLM_EXTERN 
CopyActionPtr LIBCALL
CopyActionAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   CopyActionPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* CopyAction ::= (self contained) */
      atp = AsnReadId(aip, amp, COPY_ACTION);
   } else {
      atp = AsnLinkType(orig, COPY_ACTION);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = CopyActionNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == COPY_ACTION_fields) {
      ptr -> fields = FieldPairTypeAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == COPY_ACTION_existing_text) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> existing_text = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = CopyActionFree(ptr);
   goto ret;
}



/**************************************************
*
*    CopyActionAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
CopyActionAsnWrite(CopyActionPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, COPY_ACTION);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> fields != NULL) {
      if ( ! FieldPairTypeAsnWrite(ptr -> fields, aip, COPY_ACTION_fields)) {
         goto erret;
      }
   }
   av.intvalue = ptr -> existing_text;
   retval = AsnWrite(aip, COPY_ACTION_existing_text,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    SwapActionNew()
*
**************************************************/
NLM_EXTERN 
SwapActionPtr LIBCALL
SwapActionNew(void)
{
   SwapActionPtr ptr = MemNew((size_t) sizeof(SwapAction));

   return ptr;

}


/**************************************************
*
*    SwapActionFree()
*
**************************************************/
NLM_EXTERN 
SwapActionPtr LIBCALL
SwapActionFree(SwapActionPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   FieldPairTypeFree(ptr -> fields);
   FieldTypeFree(ptr -> field_to);
   return MemFree(ptr);
}


/**************************************************
*
*    SwapActionAsnRead()
*
**************************************************/
NLM_EXTERN 
SwapActionPtr LIBCALL
SwapActionAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   SwapActionPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* SwapAction ::= (self contained) */
      atp = AsnReadId(aip, amp, SWAP_ACTION);
   } else {
      atp = AsnLinkType(orig, SWAP_ACTION);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = SwapActionNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == SWAP_ACTION_fields) {
      ptr -> fields = FieldPairTypeAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SWAP_ACTION_field_to) {
      ptr -> field_to = FieldTypeAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = SwapActionFree(ptr);
   goto ret;
}



/**************************************************
*
*    SwapActionAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
SwapActionAsnWrite(SwapActionPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, SWAP_ACTION);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> fields != NULL) {
      if ( ! FieldPairTypeAsnWrite(ptr -> fields, aip, SWAP_ACTION_fields)) {
         goto erret;
      }
   }
   if (ptr -> field_to != NULL) {
      if ( ! FieldTypeAsnWrite(ptr -> field_to, aip, SWAP_ACTION_field_to)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    AECRParseActionNew()
*
**************************************************/
NLM_EXTERN 
AECRParseActionPtr LIBCALL
AECRParseActionNew(void)
{
   AECRParseActionPtr ptr = MemNew((size_t) sizeof(AECRParseAction));

   ptr -> remove_from_parsed = 0;
   ptr -> remove_left = 0;
   ptr -> remove_right = 0;
   return ptr;

}


/**************************************************
*
*    AECRParseActionFree()
*
**************************************************/
NLM_EXTERN 
AECRParseActionPtr LIBCALL
AECRParseActionFree(AECRParseActionPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   TextPortionFree(ptr -> portion);
   FieldPairTypeFree(ptr -> fields);
   TextTransformSetFree(ptr -> transform);
   return MemFree(ptr);
}


/**************************************************
*
*    AECRParseActionAsnRead()
*
**************************************************/
NLM_EXTERN 
AECRParseActionPtr LIBCALL
AECRParseActionAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   AECRParseActionPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* AECRParseAction ::= (self contained) */
      atp = AsnReadId(aip, amp, AECRPARSE_ACTION);
   } else {
      atp = AsnLinkType(orig, AECRPARSE_ACTION);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = AECRParseActionNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == AECRPARSE_ACTION_portion) {
      ptr -> portion = TextPortionAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == AECRPARSE_ACTION_fields) {
      ptr -> fields = FieldPairTypeAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == AECRPARSE_ACTION_remove_from_parsed) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> remove_from_parsed = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == AECRPARSE_ACTION_remove_left) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> remove_left = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == AECRPARSE_ACTION_remove_right) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> remove_right = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == AECRPARSE_ACTION_transform) {
      ptr -> transform = TextTransformSetAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == AECRPARSE_ACTION_existing_text) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> existing_text = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = AECRParseActionFree(ptr);
   goto ret;
}



/**************************************************
*
*    AECRParseActionAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
AECRParseActionAsnWrite(AECRParseActionPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, AECRPARSE_ACTION);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> portion != NULL) {
      if ( ! TextPortionAsnWrite(ptr -> portion, aip, AECRPARSE_ACTION_portion)) {
         goto erret;
      }
   }
   if (ptr -> fields != NULL) {
      if ( ! FieldPairTypeAsnWrite(ptr -> fields, aip, AECRPARSE_ACTION_fields)) {
         goto erret;
      }
   }
   av.boolvalue = ptr -> remove_from_parsed;
   retval = AsnWrite(aip, AECRPARSE_ACTION_remove_from_parsed,  &av);
   av.boolvalue = ptr -> remove_left;
   retval = AsnWrite(aip, AECRPARSE_ACTION_remove_left,  &av);
   av.boolvalue = ptr -> remove_right;
   retval = AsnWrite(aip, AECRPARSE_ACTION_remove_right,  &av);
   if (ptr -> transform != NULL) {
      if ( ! TextTransformSetAsnWrite(ptr -> transform, aip, AECRPARSE_ACTION_transform)) {
         goto erret;
      }
   }
   av.intvalue = ptr -> existing_text;
   retval = AsnWrite(aip, AECRPARSE_ACTION_existing_text,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    RemoveActionNew()
*
**************************************************/
NLM_EXTERN 
RemoveActionPtr LIBCALL
RemoveActionNew(void)
{
   RemoveActionPtr ptr = MemNew((size_t) sizeof(RemoveAction));

   return ptr;

}


/**************************************************
*
*    RemoveActionFree()
*
**************************************************/
NLM_EXTERN 
RemoveActionPtr LIBCALL
RemoveActionFree(RemoveActionPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   FieldTypeFree(ptr -> field);
   return MemFree(ptr);
}


/**************************************************
*
*    RemoveActionAsnRead()
*
**************************************************/
NLM_EXTERN 
RemoveActionPtr LIBCALL
RemoveActionAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   RemoveActionPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* RemoveAction ::= (self contained) */
      atp = AsnReadId(aip, amp, REMOVE_ACTION);
   } else {
      atp = AsnLinkType(orig, REMOVE_ACTION);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = RemoveActionNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == REMOVE_ACTION_field) {
      ptr -> field = FieldTypeAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = RemoveActionFree(ptr);
   goto ret;
}



/**************************************************
*
*    RemoveActionAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
RemoveActionAsnWrite(RemoveActionPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, REMOVE_ACTION);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> field != NULL) {
      if ( ! FieldTypeAsnWrite(ptr -> field, aip, REMOVE_ACTION_field)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ActionChoiceFree()
*
**************************************************/
NLM_EXTERN 
ActionChoicePtr LIBCALL
ActionChoiceFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case ActionChoice_apply:
      ApplyActionFree(anp -> data.ptrvalue);
      break;
   case ActionChoice_edit:
      EditActionFree(anp -> data.ptrvalue);
      break;
   case ActionChoice_convert:
      ConvertActionFree(anp -> data.ptrvalue);
      break;
   case ActionChoice_copy:
      CopyActionFree(anp -> data.ptrvalue);
      break;
   case ActionChoice_swap:
      SwapActionFree(anp -> data.ptrvalue);
      break;
   case ActionChoice_remove:
      RemoveActionFree(anp -> data.ptrvalue);
      break;
   case ActionChoice_parse:
      AECRParseActionFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    ActionChoiceAsnRead()
*
**************************************************/
NLM_EXTERN 
ActionChoicePtr LIBCALL
ActionChoiceAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ActionChoice ::= (self contained) */
      atp = AsnReadId(aip, amp, ACTION_CHOICE);
   } else {
      atp = AsnLinkType(orig, ACTION_CHOICE);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == ACTION_CHOICE_apply) {
      choice = ActionChoice_apply;
      func = (AsnReadFunc) ApplyActionAsnRead;
   }
   else if (atp == ACTION_CHOICE_edit) {
      choice = ActionChoice_edit;
      func = (AsnReadFunc) EditActionAsnRead;
   }
   else if (atp == ACTION_CHOICE_convert) {
      choice = ActionChoice_convert;
      func = (AsnReadFunc) ConvertActionAsnRead;
   }
   else if (atp == ACTION_CHOICE_copy) {
      choice = ActionChoice_copy;
      func = (AsnReadFunc) CopyActionAsnRead;
   }
   else if (atp == ACTION_CHOICE_swap) {
      choice = ActionChoice_swap;
      func = (AsnReadFunc) SwapActionAsnRead;
   }
   else if (atp == ACTION_CHOICE_remove) {
      choice = ActionChoice_remove;
      func = (AsnReadFunc) RemoveActionAsnRead;
   }
   else if (atp == ACTION_CHOICE_parse) {
      choice = ActionChoice_parse;
      func = (AsnReadFunc) AECRParseActionAsnRead;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    ActionChoiceAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ActionChoiceAsnWrite(ActionChoicePtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, ACTION_CHOICE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case ActionChoice_apply:
      writetype = ACTION_CHOICE_apply;
      func = (AsnWriteFunc) ApplyActionAsnWrite;
      break;
   case ActionChoice_edit:
      writetype = ACTION_CHOICE_edit;
      func = (AsnWriteFunc) EditActionAsnWrite;
      break;
   case ActionChoice_convert:
      writetype = ACTION_CHOICE_convert;
      func = (AsnWriteFunc) ConvertActionAsnWrite;
      break;
   case ActionChoice_copy:
      writetype = ACTION_CHOICE_copy;
      func = (AsnWriteFunc) CopyActionAsnWrite;
      break;
   case ActionChoice_swap:
      writetype = ACTION_CHOICE_swap;
      func = (AsnWriteFunc) SwapActionAsnWrite;
      break;
   case ActionChoice_remove:
      writetype = ACTION_CHOICE_remove;
      func = (AsnWriteFunc) RemoveActionAsnWrite;
      break;
   case ActionChoice_parse:
      writetype = ACTION_CHOICE_parse;
      func = (AsnWriteFunc) AECRParseActionAsnWrite;
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    ParseSrcOrgChoiceFree()
*
**************************************************/
NLM_EXTERN 
ParseSrcOrgChoicePtr LIBCALL
ParseSrcOrgChoiceFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    ParseSrcOrgChoiceAsnRead()
*
**************************************************/
NLM_EXTERN 
ParseSrcOrgChoicePtr LIBCALL
ParseSrcOrgChoiceAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ParseSrcOrgChoice ::= (self contained) */
      atp = AsnReadId(aip, amp, PARSE_SRC_ORG_CHOICE);
   } else {
      atp = AsnLinkType(orig, PARSE_SRC_ORG_CHOICE);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == PARSE_SRC_ORG_CHOICE_source_qual) {
      choice = ParseSrcOrgChoice_source_qual;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == PARSE_SRC_ORG_CHOICE_taxname_after_binomial) {
      choice = ParseSrcOrgChoice_taxname_after_binomial;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    ParseSrcOrgChoiceAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ParseSrcOrgChoiceAsnWrite(ParseSrcOrgChoicePtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, PARSE_SRC_ORG_CHOICE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case ParseSrcOrgChoice_source_qual:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, PARSE_SRC_ORG_CHOICE_source_qual, &av);
      break;
   case ParseSrcOrgChoice_taxname_after_binomial:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, PARSE_SRC_ORG_CHOICE_taxname_after_binomial, &av);
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    ParseSrcOrgNew()
*
**************************************************/
NLM_EXTERN 
ParseSrcOrgPtr LIBCALL
ParseSrcOrgNew(void)
{
   ParseSrcOrgPtr ptr = MemNew((size_t) sizeof(ParseSrcOrg));

   ptr -> type = 0;
   return ptr;

}


/**************************************************
*
*    ParseSrcOrgFree()
*
**************************************************/
NLM_EXTERN 
ParseSrcOrgPtr LIBCALL
ParseSrcOrgFree(ParseSrcOrgPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   ParseSrcOrgChoiceFree(ptr -> field);
   return MemFree(ptr);
}


/**************************************************
*
*    ParseSrcOrgAsnRead()
*
**************************************************/
NLM_EXTERN 
ParseSrcOrgPtr LIBCALL
ParseSrcOrgAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ParseSrcOrgPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ParseSrcOrg ::= (self contained) */
      atp = AsnReadId(aip, amp, PARSE_SRC_ORG);
   } else {
      atp = AsnLinkType(orig, PARSE_SRC_ORG);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ParseSrcOrgNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == PARSE_SRC_ORG_field) {
      ptr -> field = ParseSrcOrgChoiceAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == PARSE_SRC_ORG_type) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> type = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ParseSrcOrgFree(ptr);
   goto ret;
}



/**************************************************
*
*    ParseSrcOrgAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ParseSrcOrgAsnWrite(ParseSrcOrgPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, PARSE_SRC_ORG);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> field != NULL) {
      if ( ! ParseSrcOrgChoiceAsnWrite(ptr -> field, aip, PARSE_SRC_ORG_field)) {
         goto erret;
      }
   }
   av.intvalue = ptr -> type;
   retval = AsnWrite(aip, PARSE_SRC_ORG_type,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ParseSrcGeneralIdFree()
*
**************************************************/
NLM_EXTERN 
ParseSrcGeneralIdPtr LIBCALL
ParseSrcGeneralIdFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case ParseSrcGeneralId_tag:
      MemFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    ParseSrcGeneralIdAsnRead()
*
**************************************************/
NLM_EXTERN 
ParseSrcGeneralIdPtr LIBCALL
ParseSrcGeneralIdAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ParseSrcGeneralId ::= (self contained) */
      atp = AsnReadId(aip, amp, PARSE_SRC_GENERAL_ID);
   } else {
      atp = AsnLinkType(orig, PARSE_SRC_GENERAL_ID);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == PARSE_SRC_GENERAL_ID_whole_text) {
      choice = ParseSrcGeneralId_whole_text;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == PARSE_SRC_GENERAL_ID_db) {
      choice = ParseSrcGeneralId_db;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == PARSE_SRC_GENERAL_ID_tag) {
      choice = ParseSrcGeneralId_tag;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.ptrvalue = av.ptrvalue;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    ParseSrcGeneralIdAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ParseSrcGeneralIdAsnWrite(ParseSrcGeneralIdPtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, PARSE_SRC_GENERAL_ID);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case ParseSrcGeneralId_whole_text:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, PARSE_SRC_GENERAL_ID_whole_text, &av);
      break;
   case ParseSrcGeneralId_db:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, PARSE_SRC_GENERAL_ID_db, &av);
      break;
   case ParseSrcGeneralId_tag:
      av.ptrvalue = anp->data.ptrvalue;
      retval = AsnWrite(aip, PARSE_SRC_GENERAL_ID_tag, &av);
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    ParseSrcFree()
*
**************************************************/
NLM_EXTERN 
ParseSrcPtr LIBCALL
ParseSrcFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case ParseSrc_org:
      ParseSrcOrgFree(anp -> data.ptrvalue);
      break;
   case ParseSrc_structured_comment:
      MemFree(anp -> data.ptrvalue);
      break;
   case ParseSrc_general_id:
      ParseSrcGeneralIdFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    ParseSrcAsnRead()
*
**************************************************/
NLM_EXTERN 
ParseSrcPtr LIBCALL
ParseSrcAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ParseSrc ::= (self contained) */
      atp = AsnReadId(aip, amp, PARSE_SRC);
   } else {
      atp = AsnLinkType(orig, PARSE_SRC);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == PARSE_SRC_defline) {
      choice = ParseSrc_defline;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == PARSE_SRC_flatfile) {
      choice = ParseSrc_flatfile;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == PARSE_SRC_local_id) {
      choice = ParseSrc_local_id;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == PARSE_SRC_org) {
      choice = ParseSrc_org;
      func = (AsnReadFunc) ParseSrcOrgAsnRead;
   }
   else if (atp == PARSE_SRC_comment) {
      choice = ParseSrc_comment;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == PARSE_SRC_bankit_comment) {
      choice = ParseSrc_bankit_comment;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == PARSE_SRC_structured_comment) {
      choice = ParseSrc_structured_comment;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.ptrvalue = av.ptrvalue;
   }
   else if (atp == PARSE_SRC_file_id) {
      choice = ParseSrc_file_id;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == PARSE_SRC_general_id) {
      choice = ParseSrc_general_id;
      func = (AsnReadFunc) ParseSrcGeneralIdAsnRead;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    ParseSrcAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ParseSrcAsnWrite(ParseSrcPtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, PARSE_SRC);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case ParseSrc_defline:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, PARSE_SRC_defline, &av);
      break;
   case ParseSrc_flatfile:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, PARSE_SRC_flatfile, &av);
      break;
   case ParseSrc_local_id:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, PARSE_SRC_local_id, &av);
      break;
   case ParseSrc_org:
      writetype = PARSE_SRC_org;
      func = (AsnWriteFunc) ParseSrcOrgAsnWrite;
      break;
   case ParseSrc_comment:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, PARSE_SRC_comment, &av);
      break;
   case ParseSrc_bankit_comment:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, PARSE_SRC_bankit_comment, &av);
      break;
   case ParseSrc_structured_comment:
      av.ptrvalue = anp->data.ptrvalue;
      retval = AsnWrite(aip, PARSE_SRC_structured_comment, &av);
      break;
   case ParseSrc_file_id:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, PARSE_SRC_file_id, &av);
      break;
   case ParseSrc_general_id:
      writetype = PARSE_SRC_general_id;
      func = (AsnWriteFunc) ParseSrcGeneralIdAsnWrite;
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    ParseDstOrgNew()
*
**************************************************/
NLM_EXTERN 
ParseDstOrgPtr LIBCALL
ParseDstOrgNew(void)
{
   ParseDstOrgPtr ptr = MemNew((size_t) sizeof(ParseDstOrg));

   ptr -> type = 0;
   return ptr;

}


/**************************************************
*
*    ParseDstOrgFree()
*
**************************************************/
NLM_EXTERN 
ParseDstOrgPtr LIBCALL
ParseDstOrgFree(ParseDstOrgPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   SourceQualChoiceFree(ptr -> field);
   return MemFree(ptr);
}


/**************************************************
*
*    ParseDstOrgAsnRead()
*
**************************************************/
NLM_EXTERN 
ParseDstOrgPtr LIBCALL
ParseDstOrgAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ParseDstOrgPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ParseDstOrg ::= (self contained) */
      atp = AsnReadId(aip, amp, PARSE_DST_ORG);
   } else {
      atp = AsnLinkType(orig, PARSE_DST_ORG);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ParseDstOrgNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == PARSE_DST_ORG_field) {
      ptr -> field = SourceQualChoiceAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == PARSE_DST_ORG_type) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> type = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ParseDstOrgFree(ptr);
   goto ret;
}



/**************************************************
*
*    ParseDstOrgAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ParseDstOrgAsnWrite(ParseDstOrgPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, PARSE_DST_ORG);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> field != NULL) {
      if ( ! SourceQualChoiceAsnWrite(ptr -> field, aip, PARSE_DST_ORG_field)) {
         goto erret;
      }
   }
   av.intvalue = ptr -> type;
   retval = AsnWrite(aip, PARSE_DST_ORG_type,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ParseDestFree()
*
**************************************************/
NLM_EXTERN 
ParseDestPtr LIBCALL
ParseDestFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case ParseDest_org:
      ParseDstOrgFree(anp -> data.ptrvalue);
      break;
   case ParseDest_featqual:
      FeatureFieldLegalFree(anp -> data.ptrvalue);
      break;
   case ParseDest_dbxref:
      MemFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    ParseDestAsnRead()
*
**************************************************/
NLM_EXTERN 
ParseDestPtr LIBCALL
ParseDestAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ParseDest ::= (self contained) */
      atp = AsnReadId(aip, amp, PARSE_DEST);
   } else {
      atp = AsnLinkType(orig, PARSE_DEST);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == PARSE_DEST_defline) {
      choice = ParseDest_defline;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == PARSE_DEST_org) {
      choice = ParseDest_org;
      func = (AsnReadFunc) ParseDstOrgAsnRead;
   }
   else if (atp == PARSE_DEST_featqual) {
      choice = ParseDest_featqual;
      func = (AsnReadFunc) FeatureFieldLegalAsnRead;
   }
   else if (atp == PARSE_DEST_comment_descriptor) {
      choice = ParseDest_comment_descriptor;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == PARSE_DEST_dbxref) {
      choice = ParseDest_dbxref;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.ptrvalue = av.ptrvalue;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    ParseDestAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ParseDestAsnWrite(ParseDestPtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, PARSE_DEST);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case ParseDest_defline:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, PARSE_DEST_defline, &av);
      break;
   case ParseDest_org:
      writetype = PARSE_DEST_org;
      func = (AsnWriteFunc) ParseDstOrgAsnWrite;
      break;
   case ParseDest_featqual:
      writetype = PARSE_DEST_featqual;
      func = (AsnWriteFunc) FeatureFieldLegalAsnWrite;
      break;
   case ParseDest_comment_descriptor:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, PARSE_DEST_comment_descriptor, &av);
      break;
   case ParseDest_dbxref:
      av.ptrvalue = anp->data.ptrvalue;
      retval = AsnWrite(aip, PARSE_DEST_dbxref, &av);
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    LocationIntervalNew()
*
**************************************************/
NLM_EXTERN 
LocationIntervalPtr LIBCALL
LocationIntervalNew(void)
{
   LocationIntervalPtr ptr = MemNew((size_t) sizeof(LocationInterval));

   return ptr;

}


/**************************************************
*
*    LocationIntervalFree()
*
**************************************************/
NLM_EXTERN 
LocationIntervalPtr LIBCALL
LocationIntervalFree(LocationIntervalPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   return MemFree(ptr);
}


/**************************************************
*
*    LocationIntervalAsnRead()
*
**************************************************/
NLM_EXTERN 
LocationIntervalPtr LIBCALL
LocationIntervalAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   LocationIntervalPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* LocationInterval ::= (self contained) */
      atp = AsnReadId(aip, amp, LOCATION_INTERVAL);
   } else {
      atp = AsnLinkType(orig, LOCATION_INTERVAL);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = LocationIntervalNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == LOCATION_INTERVAL_from) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> from = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == LOCATION_INTERVAL_to) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> to = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = LocationIntervalFree(ptr);
   goto ret;
}



/**************************************************
*
*    LocationIntervalAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
LocationIntervalAsnWrite(LocationIntervalPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, LOCATION_INTERVAL);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> from;
   retval = AsnWrite(aip, LOCATION_INTERVAL_from,  &av);
   av.intvalue = ptr -> to;
   retval = AsnWrite(aip, LOCATION_INTERVAL_to,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    LocationChoiceFree()
*
**************************************************/
NLM_EXTERN 
LocationChoicePtr LIBCALL
LocationChoiceFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case LocationChoice_interval:
      LocationIntervalFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    LocationChoiceAsnRead()
*
**************************************************/
NLM_EXTERN 
LocationChoicePtr LIBCALL
LocationChoiceAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* LocationChoice ::= (self contained) */
      atp = AsnReadId(aip, amp, LOCATION_CHOICE);
   } else {
      atp = AsnLinkType(orig, LOCATION_CHOICE);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == LOCATION_CHOICE_interval) {
      choice = LocationChoice_interval;
      func = (AsnReadFunc) LocationIntervalAsnRead;
   }
   else if (atp == LOCATION_CHOICE_whole_sequence) {
      choice = LocationChoice_whole_sequence;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == LOCATION_CHOICE_point) {
      choice = LocationChoice_point;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    LocationChoiceAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
LocationChoiceAsnWrite(LocationChoicePtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, LOCATION_CHOICE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case LocationChoice_interval:
      writetype = LOCATION_CHOICE_interval;
      func = (AsnWriteFunc) LocationIntervalAsnWrite;
      break;
   case LocationChoice_whole_sequence:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, LOCATION_CHOICE_whole_sequence, &av);
      break;
   case LocationChoice_point:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, LOCATION_CHOICE_point, &av);
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    SequenceListFree()
*
**************************************************/
NLM_EXTERN 
SequenceListPtr LIBCALL
SequenceListFree(SequenceListPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericBaseSeqOfFree(ptr,ASNCODE_PTRVAL_SLOT);
   return NULL;
}


/**************************************************
*
*    SequenceListAsnRead()
*
**************************************************/
NLM_EXTERN 
SequenceListPtr LIBCALL
SequenceListAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   SequenceListPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* SequenceList ::= (self contained) */
      atp = AsnReadId(aip, amp, SEQUENCE_LIST);
   } else {
      atp = AsnLinkType(orig, SEQUENCE_LIST);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   func = NULL;

   ptr  = AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_PTRVAL_SLOT, &isError);
   if (isError && ptr  == NULL) {
      goto erret;
   }



ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = SequenceListFree(ptr);
   goto ret;
}



/**************************************************
*
*    SequenceListAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
SequenceListAsnWrite(SequenceListPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, SEQUENCE_LIST);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   retval = AsnGenericBaseSeqOfAsnWrite(ptr, ASNCODE_PTRVAL_SLOT, aip, atp, SEQUENCE_LIST_E);
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    SequenceListChoiceFree()
*
**************************************************/
NLM_EXTERN 
SequenceListChoicePtr LIBCALL
SequenceListChoiceFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case SequenceListChoice_list:
      SequenceListFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    SequenceListChoiceAsnRead()
*
**************************************************/
NLM_EXTERN 
SequenceListChoicePtr LIBCALL
SequenceListChoiceAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* SequenceListChoice ::= (self contained) */
      atp = AsnReadId(aip, amp, SEQUENCE_LIST_CHOICE);
   } else {
      atp = AsnLinkType(orig, SEQUENCE_LIST_CHOICE);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == SEQUENCE_LIST_CHOICE_list) {
      choice = SequenceListChoice_list;
      func = (AsnReadFunc) SequenceListAsnRead;
   }
   else if (atp == SEQUENCE_LIST_CHOICE_all) {
      choice = SequenceListChoice_all;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    SequenceListChoiceAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
SequenceListChoiceAsnWrite(SequenceListChoicePtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, SEQUENCE_LIST_CHOICE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case SequenceListChoice_list:
      writetype = SEQUENCE_LIST_CHOICE_list;
      func = (AsnWriteFunc) SequenceListAsnWrite;
      break;
   case SequenceListChoice_all:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, SEQUENCE_LIST_CHOICE_all, &av);
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    ApplyFeatureActionNew()
*
**************************************************/
NLM_EXTERN 
ApplyFeatureActionPtr LIBCALL
ApplyFeatureActionNew(void)
{
   ApplyFeatureActionPtr ptr = MemNew((size_t) sizeof(ApplyFeatureAction));

   ptr -> partial5 = 0;
   ptr -> partial3 = 0;
   ptr -> plus_strand = 1;
   ptr -> add_redundant = 1;
   ptr -> add_mrna = 0;
   ptr -> apply_to_parts = 0;
   ptr -> only_seg_num = -1;
   return ptr;

}


/**************************************************
*
*    ApplyFeatureActionFree()
*
**************************************************/
NLM_EXTERN 
ApplyFeatureActionPtr LIBCALL
ApplyFeatureActionFree(ApplyFeatureActionPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   LocationChoiceFree(ptr -> location);
   SequenceListChoiceFree(ptr -> seq_list);
   FeatQualLegalSetFree(ptr -> fields);
   SourceQualValSetFree(ptr -> src_fields);
   return MemFree(ptr);
}


/**************************************************
*
*    ApplyFeatureActionAsnRead()
*
**************************************************/
NLM_EXTERN 
ApplyFeatureActionPtr LIBCALL
ApplyFeatureActionAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ApplyFeatureActionPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ApplyFeatureAction ::= (self contained) */
      atp = AsnReadId(aip, amp, APPLY_FEATURE_ACTION);
   } else {
      atp = AsnLinkType(orig, APPLY_FEATURE_ACTION);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ApplyFeatureActionNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == APPLY_FEATURE_ACTION_type) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> type = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == APPLY_FEATURE_ACTION_partial5) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> partial5 = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == APPLY_FEATURE_ACTION_partial3) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> partial3 = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == APPLY_FEATURE_ACTION_plus_strand) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> plus_strand = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == APPLY_FEATURE_ACTION_location) {
      ptr -> location = LocationChoiceAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == APPLY_FEATURE_ACTION_seq_list) {
      ptr -> seq_list = SequenceListChoiceAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == APPLY_FEATURE_ACTION_add_redundant) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> add_redundant = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == APPLY_FEATURE_ACTION_add_mrna) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> add_mrna = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == APPLY_FEATURE_ACTION_apply_to_parts) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> apply_to_parts = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == APPLY_FEATURE_ACTION_only_seg_num) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> only_seg_num = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == APPLY_FEATURE_ACTION_fields) {
      ptr -> fields = FeatQualLegalSetAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == APPLY_FEATURE_ACTION_src_fields) {
      ptr -> src_fields = SourceQualValSetAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ApplyFeatureActionFree(ptr);
   goto ret;
}



/**************************************************
*
*    ApplyFeatureActionAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ApplyFeatureActionAsnWrite(ApplyFeatureActionPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, APPLY_FEATURE_ACTION);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> type;
   retval = AsnWrite(aip, APPLY_FEATURE_ACTION_type,  &av);
   av.boolvalue = ptr -> partial5;
   retval = AsnWrite(aip, APPLY_FEATURE_ACTION_partial5,  &av);
   av.boolvalue = ptr -> partial3;
   retval = AsnWrite(aip, APPLY_FEATURE_ACTION_partial3,  &av);
   av.boolvalue = ptr -> plus_strand;
   retval = AsnWrite(aip, APPLY_FEATURE_ACTION_plus_strand,  &av);
   if (ptr -> location != NULL) {
      if ( ! LocationChoiceAsnWrite(ptr -> location, aip, APPLY_FEATURE_ACTION_location)) {
         goto erret;
      }
   }
   if (ptr -> seq_list != NULL) {
      if ( ! SequenceListChoiceAsnWrite(ptr -> seq_list, aip, APPLY_FEATURE_ACTION_seq_list)) {
         goto erret;
      }
   }
   av.boolvalue = ptr -> add_redundant;
   retval = AsnWrite(aip, APPLY_FEATURE_ACTION_add_redundant,  &av);
   av.boolvalue = ptr -> add_mrna;
   retval = AsnWrite(aip, APPLY_FEATURE_ACTION_add_mrna,  &av);
   av.boolvalue = ptr -> apply_to_parts;
   retval = AsnWrite(aip, APPLY_FEATURE_ACTION_apply_to_parts,  &av);
   av.intvalue = ptr -> only_seg_num;
   retval = AsnWrite(aip, APPLY_FEATURE_ACTION_only_seg_num,  &av);
   if (ptr -> fields != NULL) {
      if ( ! FeatQualLegalSetAsnWrite(ptr -> fields, aip, APPLY_FEATURE_ACTION_fields)) {
         goto erret;
      }
   }
   if (ptr -> src_fields != NULL) {
      if ( ! SourceQualValSetAsnWrite(ptr -> src_fields, aip, APPLY_FEATURE_ACTION_src_fields)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    RemoveFeatureActionNew()
*
**************************************************/
NLM_EXTERN 
RemoveFeatureActionPtr LIBCALL
RemoveFeatureActionNew(void)
{
   RemoveFeatureActionPtr ptr = MemNew((size_t) sizeof(RemoveFeatureAction));

   return ptr;

}


/**************************************************
*
*    RemoveFeatureActionFree()
*
**************************************************/
NLM_EXTERN 
RemoveFeatureActionPtr LIBCALL
RemoveFeatureActionFree(RemoveFeatureActionPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   ConstraintChoiceSetFree(ptr -> constraint);
   return MemFree(ptr);
}


/**************************************************
*
*    RemoveFeatureActionAsnRead()
*
**************************************************/
NLM_EXTERN 
RemoveFeatureActionPtr LIBCALL
RemoveFeatureActionAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   RemoveFeatureActionPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* RemoveFeatureAction ::= (self contained) */
      atp = AsnReadId(aip, amp, REMOVE_FEATURE_ACTION);
   } else {
      atp = AsnLinkType(orig, REMOVE_FEATURE_ACTION);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = RemoveFeatureActionNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == REMOVE_FEATURE_ACTION_type) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> type = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == REMOVE_FEATURE_ACTION_constraint) {
      ptr -> constraint = ConstraintChoiceSetAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = RemoveFeatureActionFree(ptr);
   goto ret;
}



/**************************************************
*
*    RemoveFeatureActionAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
RemoveFeatureActionAsnWrite(RemoveFeatureActionPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, REMOVE_FEATURE_ACTION);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> type;
   retval = AsnWrite(aip, REMOVE_FEATURE_ACTION_type,  &av);
   if (ptr -> constraint != NULL) {
      if ( ! ConstraintChoiceSetAsnWrite(ptr -> constraint, aip, REMOVE_FEATURE_ACTION_constraint)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ConvertFromCDSOptionsNew()
*
**************************************************/
NLM_EXTERN 
ConvertFromCDSOptionsPtr LIBCALL
ConvertFromCDSOptionsNew(void)
{
   ConvertFromCDSOptionsPtr ptr = MemNew((size_t) sizeof(ConvertFromCDSOptions));

   return ptr;

}


/**************************************************
*
*    ConvertFromCDSOptionsFree()
*
**************************************************/
NLM_EXTERN 
ConvertFromCDSOptionsPtr LIBCALL
ConvertFromCDSOptionsFree(ConvertFromCDSOptionsPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   return MemFree(ptr);
}


/**************************************************
*
*    ConvertFromCDSOptionsAsnRead()
*
**************************************************/
NLM_EXTERN 
ConvertFromCDSOptionsPtr LIBCALL
ConvertFromCDSOptionsAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ConvertFromCDSOptionsPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ConvertFromCDSOptions ::= (self contained) */
      atp = AsnReadId(aip, amp, CONVERT_FROM_CDS_OPTIONS);
   } else {
      atp = AsnLinkType(orig, CONVERT_FROM_CDS_OPTIONS);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ConvertFromCDSOptionsNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == CONVERT_FROM_CDS_OPTIONS_remove_mRNA) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> remove_mRNA = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CONVERT_FROM_CDS_OPTIONS_remove_gene) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> remove_gene = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CONVERT_FROM_CDS_OPTIONS_remove_transcript_id) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> remove_transcript_id = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ConvertFromCDSOptionsFree(ptr);
   goto ret;
}



/**************************************************
*
*    ConvertFromCDSOptionsAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ConvertFromCDSOptionsAsnWrite(ConvertFromCDSOptionsPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, CONVERT_FROM_CDS_OPTIONS);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.boolvalue = ptr -> remove_mRNA;
   retval = AsnWrite(aip, CONVERT_FROM_CDS_OPTIONS_remove_mRNA,  &av);
   av.boolvalue = ptr -> remove_gene;
   retval = AsnWrite(aip, CONVERT_FROM_CDS_OPTIONS_remove_gene,  &av);
   av.boolvalue = ptr -> remove_transcript_id;
   retval = AsnWrite(aip, CONVERT_FROM_CDS_OPTIONS_remove_transcript_id,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ConvertFeatureSrcOptionsFree()
*
**************************************************/
NLM_EXTERN 
ConvertFeatureSrcOptionsPtr LIBCALL
ConvertFeatureSrcOptionsFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case ConvertFeatureSrcOptions_cds:
      ConvertFromCDSOptionsFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    ConvertFeatureSrcOptionsAsnRead()
*
**************************************************/
NLM_EXTERN 
ConvertFeatureSrcOptionsPtr LIBCALL
ConvertFeatureSrcOptionsAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ConvertFeatureSrcOptions ::= (self contained) */
      atp = AsnReadId(aip, amp, CONVERT_FEATURE_SRC_OPTIONS);
   } else {
      atp = AsnLinkType(orig, CONVERT_FEATURE_SRC_OPTIONS);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == CONVERT_FEATURE_SRC_OPTIONS_cds) {
      choice = ConvertFeatureSrcOptions_cds;
      func = (AsnReadFunc) ConvertFromCDSOptionsAsnRead;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    ConvertFeatureSrcOptionsAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ConvertFeatureSrcOptionsAsnWrite(ConvertFeatureSrcOptionsPtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, CONVERT_FEATURE_SRC_OPTIONS);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case ConvertFeatureSrcOptions_cds:
      writetype = CONVERT_FEATURE_SRC_OPTIONS_cds;
      func = (AsnWriteFunc) ConvertFromCDSOptionsAsnWrite;
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    RegionTypeNew()
*
**************************************************/
NLM_EXTERN 
RegionTypePtr LIBCALL
RegionTypeNew(void)
{
   RegionTypePtr ptr = MemNew((size_t) sizeof(RegionType));

   return ptr;

}


/**************************************************
*
*    RegionTypeFree()
*
**************************************************/
NLM_EXTERN 
RegionTypePtr LIBCALL
RegionTypeFree(RegionTypePtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   return MemFree(ptr);
}


/**************************************************
*
*    RegionTypeAsnRead()
*
**************************************************/
NLM_EXTERN 
RegionTypePtr LIBCALL
RegionTypeAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   RegionTypePtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* RegionType ::= (self contained) */
      atp = AsnReadId(aip, amp, REGION_TYPE);
   } else {
      atp = AsnLinkType(orig, REGION_TYPE);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = RegionTypeNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == REGION_TYPE_create_nucleotide) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> create_nucleotide = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = RegionTypeFree(ptr);
   goto ret;
}



/**************************************************
*
*    RegionTypeAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
RegionTypeAsnWrite(RegionTypePtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, REGION_TYPE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.boolvalue = ptr -> create_nucleotide;
   retval = AsnWrite(aip, REGION_TYPE_create_nucleotide,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ConvertFeatureDstOptionsFree()
*
**************************************************/
NLM_EXTERN 
ConvertFeatureDstOptionsPtr LIBCALL
ConvertFeatureDstOptionsFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case ConvertFeatureDstOptions_region:
      RegionTypeFree(anp -> data.ptrvalue);
      break;
   case ConvertFeatureDstOptions_ncrna_class:
      MemFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    ConvertFeatureDstOptionsAsnRead()
*
**************************************************/
NLM_EXTERN 
ConvertFeatureDstOptionsPtr LIBCALL
ConvertFeatureDstOptionsAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ConvertFeatureDstOptions ::= (self contained) */
      atp = AsnReadId(aip, amp, CONVERT_FEATURE_DST_OPTIONS);
   } else {
      atp = AsnLinkType(orig, CONVERT_FEATURE_DST_OPTIONS);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == CONVERT_FEATURE_DST_OPTIONS_bond) {
      choice = ConvertFeatureDstOptions_bond;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == CONVERT_FEATURE_DST_OPTIONS_site) {
      choice = ConvertFeatureDstOptions_site;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == CONVERT_FEATURE_DST_OPTIONS_region) {
      choice = ConvertFeatureDstOptions_region;
      func = (AsnReadFunc) RegionTypeAsnRead;
   }
   else if (atp == CONVERT_FEATURE_DST_OPTIONS_ncrna_class) {
      choice = ConvertFeatureDstOptions_ncrna_class;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.ptrvalue = av.ptrvalue;
   }
   else if (atp == CONVERT_FEATURE_DST_OPTIONS_remove_original) {
      choice = ConvertFeatureDstOptions_remove_original;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    ConvertFeatureDstOptionsAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ConvertFeatureDstOptionsAsnWrite(ConvertFeatureDstOptionsPtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, CONVERT_FEATURE_DST_OPTIONS);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case ConvertFeatureDstOptions_bond:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, CONVERT_FEATURE_DST_OPTIONS_bond, &av);
      break;
   case ConvertFeatureDstOptions_site:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, CONVERT_FEATURE_DST_OPTIONS_site, &av);
      break;
   case ConvertFeatureDstOptions_region:
      writetype = CONVERT_FEATURE_DST_OPTIONS_region;
      func = (AsnWriteFunc) RegionTypeAsnWrite;
      break;
   case ConvertFeatureDstOptions_ncrna_class:
      av.ptrvalue = anp->data.ptrvalue;
      retval = AsnWrite(aip, CONVERT_FEATURE_DST_OPTIONS_ncrna_class, &av);
      break;
   case ConvertFeatureDstOptions_remove_original:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, CONVERT_FEATURE_DST_OPTIONS_remove_original, &av);
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    ConvertFeatureActionNew()
*
**************************************************/
NLM_EXTERN 
ConvertFeatureActionPtr LIBCALL
ConvertFeatureActionNew(void)
{
   ConvertFeatureActionPtr ptr = MemNew((size_t) sizeof(ConvertFeatureAction));

   return ptr;

}


/**************************************************
*
*    ConvertFeatureActionFree()
*
**************************************************/
NLM_EXTERN 
ConvertFeatureActionPtr LIBCALL
ConvertFeatureActionFree(ConvertFeatureActionPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   ConvertFeatureSrcOptionsFree(ptr -> src_options);
   ConvertFeatureDstOptionsFree(ptr -> dst_options);
   ConstraintChoiceSetFree(ptr -> src_feat_constraint);
   return MemFree(ptr);
}


/**************************************************
*
*    ConvertFeatureActionAsnRead()
*
**************************************************/
NLM_EXTERN 
ConvertFeatureActionPtr LIBCALL
ConvertFeatureActionAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ConvertFeatureActionPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ConvertFeatureAction ::= (self contained) */
      atp = AsnReadId(aip, amp, CONVERT_FEATURE_ACTION);
   } else {
      atp = AsnLinkType(orig, CONVERT_FEATURE_ACTION);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ConvertFeatureActionNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == CONVERT_FEATURE_ACTION_type_from) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> type_from = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CONVERT_FEATURE_ACTION_type_to) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> type_to = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CONVERT_FEATURE_ACTION_src_options) {
      ptr -> src_options = ConvertFeatureSrcOptionsAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CONVERT_FEATURE_ACTION_dst_options) {
      ptr -> dst_options = ConvertFeatureDstOptionsAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CONVERT_FEATURE_ACTION_leave_original) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> leave_original = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CONVERT_FEATURE_ACTION_src_feat_constraint) {
      ptr -> src_feat_constraint = ConstraintChoiceSetAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ConvertFeatureActionFree(ptr);
   goto ret;
}



/**************************************************
*
*    ConvertFeatureActionAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ConvertFeatureActionAsnWrite(ConvertFeatureActionPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, CONVERT_FEATURE_ACTION);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> type_from;
   retval = AsnWrite(aip, CONVERT_FEATURE_ACTION_type_from,  &av);
   av.intvalue = ptr -> type_to;
   retval = AsnWrite(aip, CONVERT_FEATURE_ACTION_type_to,  &av);
   if (ptr -> src_options != NULL) {
      if ( ! ConvertFeatureSrcOptionsAsnWrite(ptr -> src_options, aip, CONVERT_FEATURE_ACTION_src_options)) {
         goto erret;
      }
   }
   if (ptr -> dst_options != NULL) {
      if ( ! ConvertFeatureDstOptionsAsnWrite(ptr -> dst_options, aip, CONVERT_FEATURE_ACTION_dst_options)) {
         goto erret;
      }
   }
   av.boolvalue = ptr -> leave_original;
   retval = AsnWrite(aip, CONVERT_FEATURE_ACTION_leave_original,  &av);
   if (ptr -> src_feat_constraint != NULL) {
      if ( ! ConstraintChoiceSetAsnWrite(ptr -> src_feat_constraint, aip, CONVERT_FEATURE_ACTION_src_feat_constraint)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    EditLocationStrandNew()
*
**************************************************/
NLM_EXTERN 
EditLocationStrandPtr LIBCALL
EditLocationStrandNew(void)
{
   EditLocationStrandPtr ptr = MemNew((size_t) sizeof(EditLocationStrand));

   return ptr;

}


/**************************************************
*
*    EditLocationStrandFree()
*
**************************************************/
NLM_EXTERN 
EditLocationStrandPtr LIBCALL
EditLocationStrandFree(EditLocationStrandPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   return MemFree(ptr);
}


/**************************************************
*
*    EditLocationStrandAsnRead()
*
**************************************************/
NLM_EXTERN 
EditLocationStrandPtr LIBCALL
EditLocationStrandAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   EditLocationStrandPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* EditLocationStrand ::= (self contained) */
      atp = AsnReadId(aip, amp, EDIT_LOCATION_STRAND);
   } else {
      atp = AsnLinkType(orig, EDIT_LOCATION_STRAND);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = EditLocationStrandNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == EDIT_LOCATION_STRAND_strand_from) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> strand_from = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == EDIT_LOCATION_STRAND_strand_to) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> strand_to = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = EditLocationStrandFree(ptr);
   goto ret;
}



/**************************************************
*
*    EditLocationStrandAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
EditLocationStrandAsnWrite(EditLocationStrandPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, EDIT_LOCATION_STRAND);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> strand_from;
   retval = AsnWrite(aip, EDIT_LOCATION_STRAND_strand_from,  &av);
   av.intvalue = ptr -> strand_to;
   retval = AsnWrite(aip, EDIT_LOCATION_STRAND_strand_to,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    Partial5SetActionNew()
*
**************************************************/
NLM_EXTERN 
Partial5SetActionPtr LIBCALL
Partial5SetActionNew(void)
{
   Partial5SetActionPtr ptr = MemNew((size_t) sizeof(Partial5SetAction));

   return ptr;

}


/**************************************************
*
*    Partial5SetActionFree()
*
**************************************************/
NLM_EXTERN 
Partial5SetActionPtr LIBCALL
Partial5SetActionFree(Partial5SetActionPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   return MemFree(ptr);
}


/**************************************************
*
*    Partial5SetActionAsnRead()
*
**************************************************/
NLM_EXTERN 
Partial5SetActionPtr LIBCALL
Partial5SetActionAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   Partial5SetActionPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* Partial5SetAction ::= (self contained) */
      atp = AsnReadId(aip, amp, PARTIAL_5_SET_ACTION);
   } else {
      atp = AsnLinkType(orig, PARTIAL_5_SET_ACTION);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = Partial5SetActionNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == PARTIAL_5_SET_ACTION_constraint) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> constraint = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == PARTIAL_5_SET_ACTION_extend) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> extend = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = Partial5SetActionFree(ptr);
   goto ret;
}



/**************************************************
*
*    Partial5SetActionAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
Partial5SetActionAsnWrite(Partial5SetActionPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, PARTIAL_5_SET_ACTION);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> constraint;
   retval = AsnWrite(aip, PARTIAL_5_SET_ACTION_constraint,  &av);
   av.boolvalue = ptr -> extend;
   retval = AsnWrite(aip, PARTIAL_5_SET_ACTION_extend,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    Partial3SetActionNew()
*
**************************************************/
NLM_EXTERN 
Partial3SetActionPtr LIBCALL
Partial3SetActionNew(void)
{
   Partial3SetActionPtr ptr = MemNew((size_t) sizeof(Partial3SetAction));

   return ptr;

}


/**************************************************
*
*    Partial3SetActionFree()
*
**************************************************/
NLM_EXTERN 
Partial3SetActionPtr LIBCALL
Partial3SetActionFree(Partial3SetActionPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   return MemFree(ptr);
}


/**************************************************
*
*    Partial3SetActionAsnRead()
*
**************************************************/
NLM_EXTERN 
Partial3SetActionPtr LIBCALL
Partial3SetActionAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   Partial3SetActionPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* Partial3SetAction ::= (self contained) */
      atp = AsnReadId(aip, amp, PARTIAL_3_SET_ACTION);
   } else {
      atp = AsnLinkType(orig, PARTIAL_3_SET_ACTION);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = Partial3SetActionNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == PARTIAL_3_SET_ACTION_constraint) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> constraint = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == PARTIAL_3_SET_ACTION_extend) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> extend = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = Partial3SetActionFree(ptr);
   goto ret;
}



/**************************************************
*
*    Partial3SetActionAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
Partial3SetActionAsnWrite(Partial3SetActionPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, PARTIAL_3_SET_ACTION);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> constraint;
   retval = AsnWrite(aip, PARTIAL_3_SET_ACTION_constraint,  &av);
   av.boolvalue = ptr -> extend;
   retval = AsnWrite(aip, PARTIAL_3_SET_ACTION_extend,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    PartialBothSetActionNew()
*
**************************************************/
NLM_EXTERN 
PartialBothSetActionPtr LIBCALL
PartialBothSetActionNew(void)
{
   PartialBothSetActionPtr ptr = MemNew((size_t) sizeof(PartialBothSetAction));

   return ptr;

}


/**************************************************
*
*    PartialBothSetActionFree()
*
**************************************************/
NLM_EXTERN 
PartialBothSetActionPtr LIBCALL
PartialBothSetActionFree(PartialBothSetActionPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   return MemFree(ptr);
}


/**************************************************
*
*    PartialBothSetActionAsnRead()
*
**************************************************/
NLM_EXTERN 
PartialBothSetActionPtr LIBCALL
PartialBothSetActionAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   PartialBothSetActionPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* PartialBothSetAction ::= (self contained) */
      atp = AsnReadId(aip, amp, PARTIAL_BOTH_SET_ACTION);
   } else {
      atp = AsnLinkType(orig, PARTIAL_BOTH_SET_ACTION);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = PartialBothSetActionNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == PARTIAL_BOTH_SET_ACTION_constraint) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> constraint = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == PARTIAL_BOTH_SET_ACTION_extend) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> extend = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = PartialBothSetActionFree(ptr);
   goto ret;
}



/**************************************************
*
*    PartialBothSetActionAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
PartialBothSetActionAsnWrite(PartialBothSetActionPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, PARTIAL_BOTH_SET_ACTION);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> constraint;
   retval = AsnWrite(aip, PARTIAL_BOTH_SET_ACTION_constraint,  &av);
   av.boolvalue = ptr -> extend;
   retval = AsnWrite(aip, PARTIAL_BOTH_SET_ACTION_extend,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    LocationEditTypeFree()
*
**************************************************/
NLM_EXTERN 
LocationEditTypePtr LIBCALL
LocationEditTypeFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case LocationEditType_strand:
      EditLocationStrandFree(anp -> data.ptrvalue);
      break;
   case LocationEditType_set_5_partial:
      Partial5SetActionFree(anp -> data.ptrvalue);
      break;
   case LocationEditType_set_3_partial:
      Partial3SetActionFree(anp -> data.ptrvalue);
      break;
   case LocationEditType_set_both_partial:
      PartialBothSetActionFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    LocationEditTypeAsnRead()
*
**************************************************/
NLM_EXTERN 
LocationEditTypePtr LIBCALL
LocationEditTypeAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* LocationEditType ::= (self contained) */
      atp = AsnReadId(aip, amp, LOCATION_EDIT_TYPE);
   } else {
      atp = AsnLinkType(orig, LOCATION_EDIT_TYPE);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == LOCATION_EDIT_TYPE_strand) {
      choice = LocationEditType_strand;
      func = (AsnReadFunc) EditLocationStrandAsnRead;
   }
   else if (atp == LOCATION_EDIT_TYPE_set_5_partial) {
      choice = LocationEditType_set_5_partial;
      func = (AsnReadFunc) Partial5SetActionAsnRead;
   }
   else if (atp == LOCATION_EDIT_TYPE_clear_5_partial) {
      choice = LocationEditType_clear_5_partial;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == LOCATION_EDIT_TYPE_set_3_partial) {
      choice = LocationEditType_set_3_partial;
      func = (AsnReadFunc) Partial3SetActionAsnRead;
   }
   else if (atp == LOCATION_EDIT_TYPE_clear_3_partial) {
      choice = LocationEditType_clear_3_partial;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == LOCATION_EDIT_TYPE_set_both_partial) {
      choice = LocationEditType_set_both_partial;
      func = (AsnReadFunc) PartialBothSetActionAsnRead;
   }
   else if (atp == LOCATION_EDIT_TYPE_clear_both_partial) {
      choice = LocationEditType_clear_both_partial;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == LOCATION_EDIT_TYPE_convert) {
      choice = LocationEditType_convert;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == LOCATION_EDIT_TYPE_extend_5) {
      choice = LocationEditType_extend_5;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == LOCATION_EDIT_TYPE_extend_3) {
      choice = LocationEditType_extend_3;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    LocationEditTypeAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
LocationEditTypeAsnWrite(LocationEditTypePtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, LOCATION_EDIT_TYPE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case LocationEditType_strand:
      writetype = LOCATION_EDIT_TYPE_strand;
      func = (AsnWriteFunc) EditLocationStrandAsnWrite;
      break;
   case LocationEditType_set_5_partial:
      writetype = LOCATION_EDIT_TYPE_set_5_partial;
      func = (AsnWriteFunc) Partial5SetActionAsnWrite;
      break;
   case LocationEditType_clear_5_partial:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, LOCATION_EDIT_TYPE_clear_5_partial, &av);
      break;
   case LocationEditType_set_3_partial:
      writetype = LOCATION_EDIT_TYPE_set_3_partial;
      func = (AsnWriteFunc) Partial3SetActionAsnWrite;
      break;
   case LocationEditType_clear_3_partial:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, LOCATION_EDIT_TYPE_clear_3_partial, &av);
      break;
   case LocationEditType_set_both_partial:
      writetype = LOCATION_EDIT_TYPE_set_both_partial;
      func = (AsnWriteFunc) PartialBothSetActionAsnWrite;
      break;
   case LocationEditType_clear_both_partial:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, LOCATION_EDIT_TYPE_clear_both_partial, &av);
      break;
   case LocationEditType_convert:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, LOCATION_EDIT_TYPE_convert, &av);
      break;
   case LocationEditType_extend_5:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, LOCATION_EDIT_TYPE_extend_5, &av);
      break;
   case LocationEditType_extend_3:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, LOCATION_EDIT_TYPE_extend_3, &av);
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    EditFeatureLocationActionNew()
*
**************************************************/
NLM_EXTERN 
EditFeatureLocationActionPtr LIBCALL
EditFeatureLocationActionNew(void)
{
   EditFeatureLocationActionPtr ptr = MemNew((size_t) sizeof(EditFeatureLocationAction));

   return ptr;

}


/**************************************************
*
*    EditFeatureLocationActionFree()
*
**************************************************/
NLM_EXTERN 
EditFeatureLocationActionPtr LIBCALL
EditFeatureLocationActionFree(EditFeatureLocationActionPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   LocationEditTypeFree(ptr -> action);
   ConstraintChoiceSetFree(ptr -> constraint);
   return MemFree(ptr);
}


/**************************************************
*
*    EditFeatureLocationActionAsnRead()
*
**************************************************/
NLM_EXTERN 
EditFeatureLocationActionPtr LIBCALL
EditFeatureLocationActionAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   EditFeatureLocationActionPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* EditFeatureLocationAction ::= (self contained) */
      atp = AsnReadId(aip, amp, EDIT_FEATURE_LOCATION_ACTION);
   } else {
      atp = AsnLinkType(orig, EDIT_FEATURE_LOCATION_ACTION);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = EditFeatureLocationActionNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == EDIT_FEATURE_LOCATION_ACTION_type) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> type = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == EDIT_FEATURE_LOCATION_ACTION_action) {
      ptr -> action = LocationEditTypeAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == EDIT_FEATURE_LOCATION_ACTION_retranslate_cds) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> retranslate_cds = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == EDIT_FEATURE_LOCATION_ACTION_constraint) {
      ptr -> constraint = ConstraintChoiceSetAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = EditFeatureLocationActionFree(ptr);
   goto ret;
}



/**************************************************
*
*    EditFeatureLocationActionAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
EditFeatureLocationActionAsnWrite(EditFeatureLocationActionPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, EDIT_FEATURE_LOCATION_ACTION);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> type;
   retval = AsnWrite(aip, EDIT_FEATURE_LOCATION_ACTION_type,  &av);
   if (ptr -> action != NULL) {
      if ( ! LocationEditTypeAsnWrite(ptr -> action, aip, EDIT_FEATURE_LOCATION_ACTION_action)) {
         goto erret;
      }
   }
   av.boolvalue = ptr -> retranslate_cds;
   retval = AsnWrite(aip, EDIT_FEATURE_LOCATION_ACTION_retranslate_cds,  &av);
   if (ptr -> constraint != NULL) {
      if ( ! ConstraintChoiceSetAsnWrite(ptr -> constraint, aip, EDIT_FEATURE_LOCATION_ACTION_constraint)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    MolinfoBlockNew()
*
**************************************************/
NLM_EXTERN 
MolinfoBlockPtr LIBCALL
MolinfoBlockNew(void)
{
   MolinfoBlockPtr ptr = MemNew((size_t) sizeof(MolinfoBlock));

   return ptr;

}


/**************************************************
*
*    MolinfoBlockFree()
*
**************************************************/
NLM_EXTERN 
MolinfoBlockPtr LIBCALL
MolinfoBlockFree(MolinfoBlockPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MolinfoFieldListFree(ptr -> to_list);
   MolinfoFieldListFree(ptr -> from_list);
   ConstraintChoiceSetFree(ptr -> constraint);
   return MemFree(ptr);
}


/**************************************************
*
*    MolinfoBlockAsnRead()
*
**************************************************/
NLM_EXTERN 
MolinfoBlockPtr LIBCALL
MolinfoBlockAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   MolinfoBlockPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* MolinfoBlock ::= (self contained) */
      atp = AsnReadId(aip, amp, MOLINFO_BLOCK);
   } else {
      atp = AsnLinkType(orig, MOLINFO_BLOCK);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = MolinfoBlockNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == MOLINFO_BLOCK_to_list) {
      ptr -> to_list = MolinfoFieldListAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == MOLINFO_BLOCK_from_list) {
      ptr -> from_list = MolinfoFieldListAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == MOLINFO_BLOCK_constraint) {
      ptr -> constraint = ConstraintChoiceSetAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = MolinfoBlockFree(ptr);
   goto ret;
}



/**************************************************
*
*    MolinfoBlockAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
MolinfoBlockAsnWrite(MolinfoBlockPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, MOLINFO_BLOCK);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> to_list != NULL) {
      if ( ! MolinfoFieldListAsnWrite(ptr -> to_list, aip, MOLINFO_BLOCK_to_list)) {
         goto erret;
      }
   }
   if (ptr -> from_list != NULL) {
      if ( ! MolinfoFieldListAsnWrite(ptr -> from_list, aip, MOLINFO_BLOCK_from_list)) {
         goto erret;
      }
   }
   if (ptr -> constraint != NULL) {
      if ( ! ConstraintChoiceSetAsnWrite(ptr -> constraint, aip, MOLINFO_BLOCK_constraint)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    RemoveDescriptorActionNew()
*
**************************************************/
NLM_EXTERN 
RemoveDescriptorActionPtr LIBCALL
RemoveDescriptorActionNew(void)
{
   RemoveDescriptorActionPtr ptr = MemNew((size_t) sizeof(RemoveDescriptorAction));

   return ptr;

}


/**************************************************
*
*    RemoveDescriptorActionFree()
*
**************************************************/
NLM_EXTERN 
RemoveDescriptorActionPtr LIBCALL
RemoveDescriptorActionFree(RemoveDescriptorActionPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   ConstraintChoiceSetFree(ptr -> constraint);
   return MemFree(ptr);
}


/**************************************************
*
*    RemoveDescriptorActionAsnRead()
*
**************************************************/
NLM_EXTERN 
RemoveDescriptorActionPtr LIBCALL
RemoveDescriptorActionAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   RemoveDescriptorActionPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* RemoveDescriptorAction ::= (self contained) */
      atp = AsnReadId(aip, amp, REMOVE_DESCRIPTOR_ACTION);
   } else {
      atp = AsnLinkType(orig, REMOVE_DESCRIPTOR_ACTION);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = RemoveDescriptorActionNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == REMOVE_DESCRIPTOR_ACTION_type) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> type = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == REMOVE_DESCRIPTOR_ACTION_constraint) {
      ptr -> constraint = ConstraintChoiceSetAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = RemoveDescriptorActionFree(ptr);
   goto ret;
}



/**************************************************
*
*    RemoveDescriptorActionAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
RemoveDescriptorActionAsnWrite(RemoveDescriptorActionPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, REMOVE_DESCRIPTOR_ACTION);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> type;
   retval = AsnWrite(aip, REMOVE_DESCRIPTOR_ACTION_type,  &av);
   if (ptr -> constraint != NULL) {
      if ( ! ConstraintChoiceSetAsnWrite(ptr -> constraint, aip, REMOVE_DESCRIPTOR_ACTION_constraint)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    AutodefActionNew()
*
**************************************************/
NLM_EXTERN 
AutodefActionPtr LIBCALL
AutodefActionNew(void)
{
   AutodefActionPtr ptr = MemNew((size_t) sizeof(AutodefAction));

   return ptr;

}


/**************************************************
*
*    AutodefActionFree()
*
**************************************************/
NLM_EXTERN 
AutodefActionPtr LIBCALL
AutodefActionFree(AutodefActionPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericBaseSeqOfFree(ptr -> modifiers ,ASNCODE_INTVAL_SLOT);
   return MemFree(ptr);
}


/**************************************************
*
*    AutodefActionAsnRead()
*
**************************************************/
NLM_EXTERN 
AutodefActionPtr LIBCALL
AutodefActionAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   AutodefActionPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* AutodefAction ::= (self contained) */
      atp = AsnReadId(aip, amp, AUTODEF_ACTION);
   } else {
      atp = AsnLinkType(orig, AUTODEF_ACTION);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = AutodefActionNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == AUTODEF_ACTION_modifiers) {
      ptr -> modifiers = AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_INTVAL_SLOT, &isError);
      if (isError && ptr -> modifiers == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == AUTODEF_ACTION_clause_list_type) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> clause_list_type = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = AutodefActionFree(ptr);
   goto ret;
}



/**************************************************
*
*    AutodefActionAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
AutodefActionAsnWrite(AutodefActionPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, AUTODEF_ACTION);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   retval = AsnGenericBaseSeqOfAsnWrite(ptr -> modifiers ,ASNCODE_INTVAL_SLOT, aip, AUTODEF_ACTION_modifiers, AUTODEF_ACTION_modifiers_E);
   av.intvalue = ptr -> clause_list_type;
   retval = AsnWrite(aip, AUTODEF_ACTION_clause_list_type,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    FixPubCapsActionNew()
*
**************************************************/
NLM_EXTERN 
FixPubCapsActionPtr LIBCALL
FixPubCapsActionNew(void)
{
   FixPubCapsActionPtr ptr = MemNew((size_t) sizeof(FixPubCapsAction));

   ptr -> punct_only = 0;
   return ptr;

}


/**************************************************
*
*    FixPubCapsActionFree()
*
**************************************************/
NLM_EXTERN 
FixPubCapsActionPtr LIBCALL
FixPubCapsActionFree(FixPubCapsActionPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   ConstraintChoiceSetFree(ptr -> constraint);
   return MemFree(ptr);
}


/**************************************************
*
*    FixPubCapsActionAsnRead()
*
**************************************************/
NLM_EXTERN 
FixPubCapsActionPtr LIBCALL
FixPubCapsActionAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   FixPubCapsActionPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* FixPubCapsAction ::= (self contained) */
      atp = AsnReadId(aip, amp, FIX_PUB_CAPS_ACTION);
   } else {
      atp = AsnLinkType(orig, FIX_PUB_CAPS_ACTION);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = FixPubCapsActionNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == FIX_PUB_CAPS_ACTION_title) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> title = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == FIX_PUB_CAPS_ACTION_authors) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> authors = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == FIX_PUB_CAPS_ACTION_affiliation) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> affiliation = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == FIX_PUB_CAPS_ACTION_affil_country) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> affil_country = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == FIX_PUB_CAPS_ACTION_punct_only) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> punct_only = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == FIX_PUB_CAPS_ACTION_constraint) {
      ptr -> constraint = ConstraintChoiceSetAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = FixPubCapsActionFree(ptr);
   goto ret;
}



/**************************************************
*
*    FixPubCapsActionAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
FixPubCapsActionAsnWrite(FixPubCapsActionPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, FIX_PUB_CAPS_ACTION);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.boolvalue = ptr -> title;
   retval = AsnWrite(aip, FIX_PUB_CAPS_ACTION_title,  &av);
   av.boolvalue = ptr -> authors;
   retval = AsnWrite(aip, FIX_PUB_CAPS_ACTION_authors,  &av);
   av.boolvalue = ptr -> affiliation;
   retval = AsnWrite(aip, FIX_PUB_CAPS_ACTION_affiliation,  &av);
   av.boolvalue = ptr -> affil_country;
   retval = AsnWrite(aip, FIX_PUB_CAPS_ACTION_affil_country,  &av);
   av.boolvalue = ptr -> punct_only;
   retval = AsnWrite(aip, FIX_PUB_CAPS_ACTION_punct_only,  &av);
   if (ptr -> constraint != NULL) {
      if ( ! ConstraintChoiceSetAsnWrite(ptr -> constraint, aip, FIX_PUB_CAPS_ACTION_constraint)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    SortFieldsActionNew()
*
**************************************************/
NLM_EXTERN 
SortFieldsActionPtr LIBCALL
SortFieldsActionNew(void)
{
   SortFieldsActionPtr ptr = MemNew((size_t) sizeof(SortFieldsAction));

   return ptr;

}


/**************************************************
*
*    SortFieldsActionFree()
*
**************************************************/
NLM_EXTERN 
SortFieldsActionPtr LIBCALL
SortFieldsActionFree(SortFieldsActionPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   FieldTypeFree(ptr -> field);
   ConstraintChoiceSetFree(ptr -> constraint);
   return MemFree(ptr);
}


/**************************************************
*
*    SortFieldsActionAsnRead()
*
**************************************************/
NLM_EXTERN 
SortFieldsActionPtr LIBCALL
SortFieldsActionAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   SortFieldsActionPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* SortFieldsAction ::= (self contained) */
      atp = AsnReadId(aip, amp, SORT_FIELDS_ACTION);
   } else {
      atp = AsnLinkType(orig, SORT_FIELDS_ACTION);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = SortFieldsActionNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == SORT_FIELDS_ACTION_field) {
      ptr -> field = FieldTypeAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SORT_FIELDS_ACTION_order) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> order = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SORT_FIELDS_ACTION_constraint) {
      ptr -> constraint = ConstraintChoiceSetAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = SortFieldsActionFree(ptr);
   goto ret;
}



/**************************************************
*
*    SortFieldsActionAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
SortFieldsActionAsnWrite(SortFieldsActionPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, SORT_FIELDS_ACTION);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> field != NULL) {
      if ( ! FieldTypeAsnWrite(ptr -> field, aip, SORT_FIELDS_ACTION_field)) {
         goto erret;
      }
   }
   av.intvalue = ptr -> order;
   retval = AsnWrite(aip, SORT_FIELDS_ACTION_order,  &av);
   if (ptr -> constraint != NULL) {
      if ( ! ConstraintChoiceSetAsnWrite(ptr -> constraint, aip, SORT_FIELDS_ACTION_constraint)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    FixCapsActionFree()
*
**************************************************/
NLM_EXTERN 
FixCapsActionPtr LIBCALL
FixCapsActionFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case FixCapsAction_pub:
      FixPubCapsActionFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    FixCapsActionAsnRead()
*
**************************************************/
NLM_EXTERN 
FixCapsActionPtr LIBCALL
FixCapsActionAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* FixCapsAction ::= (self contained) */
      atp = AsnReadId(aip, amp, FIX_CAPS_ACTION);
   } else {
      atp = AsnLinkType(orig, FIX_CAPS_ACTION);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == FIX_CAPS_ACTION_pub) {
      choice = FixCapsAction_pub;
      func = (AsnReadFunc) FixPubCapsActionAsnRead;
   }
   else if (atp == FIX_CAPS_ACTION_src_country) {
      choice = FixCapsAction_src_country;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == FIX_CAPS_ACTION_mouse_strain) {
      choice = FixCapsAction_mouse_strain;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == FIX_CAPS_ACTION_src_qual) {
      choice = FixCapsAction_src_qual;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    FixCapsActionAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
FixCapsActionAsnWrite(FixCapsActionPtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, FIX_CAPS_ACTION);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case FixCapsAction_pub:
      writetype = FIX_CAPS_ACTION_pub;
      func = (AsnWriteFunc) FixPubCapsActionAsnWrite;
      break;
   case FixCapsAction_src_country:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, FIX_CAPS_ACTION_src_country, &av);
      break;
   case FixCapsAction_mouse_strain:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, FIX_CAPS_ACTION_mouse_strain, &av);
      break;
   case FixCapsAction_src_qual:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, FIX_CAPS_ACTION_src_qual, &av);
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    FixFormatActionFree()
*
**************************************************/
NLM_EXTERN 
FixFormatActionPtr LIBCALL
FixFormatActionFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    FixFormatActionAsnRead()
*
**************************************************/
NLM_EXTERN 
FixFormatActionPtr LIBCALL
FixFormatActionAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* FixFormatAction ::= (self contained) */
      atp = AsnReadId(aip, amp, FIX_FORMAT_ACTION);
   } else {
      atp = AsnLinkType(orig, FIX_FORMAT_ACTION);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == FIX_FORMAT_ACTION_collection_date) {
      choice = FixFormatAction_collection_date;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == FIX_FORMAT_ACTION_lat_lon) {
      choice = FixFormatAction_lat_lon;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == FIX_FORMAT_ACTION_primers) {
      choice = FixFormatAction_primers;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == FIX_FORMAT_ACTION_protein_name) {
      choice = FixFormatAction_protein_name;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    FixFormatActionAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
FixFormatActionAsnWrite(FixFormatActionPtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, FIX_FORMAT_ACTION);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case FixFormatAction_collection_date:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, FIX_FORMAT_ACTION_collection_date, &av);
      break;
   case FixFormatAction_lat_lon:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, FIX_FORMAT_ACTION_lat_lon, &av);
      break;
   case FixFormatAction_primers:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, FIX_FORMAT_ACTION_primers, &av);
      break;
   case FixFormatAction_protein_name:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, FIX_FORMAT_ACTION_protein_name, &av);
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    RemoveDuplicateFeatureActionNew()
*
**************************************************/
NLM_EXTERN 
RemoveDuplicateFeatureActionPtr LIBCALL
RemoveDuplicateFeatureActionNew(void)
{
   RemoveDuplicateFeatureActionPtr ptr = MemNew((size_t) sizeof(RemoveDuplicateFeatureAction));

   return ptr;

}


/**************************************************
*
*    RemoveDuplicateFeatureActionFree()
*
**************************************************/
NLM_EXTERN 
RemoveDuplicateFeatureActionPtr LIBCALL
RemoveDuplicateFeatureActionFree(RemoveDuplicateFeatureActionPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   ConstraintChoiceSetFree(ptr -> rd_constraint);
   return MemFree(ptr);
}


/**************************************************
*
*    RemoveDuplicateFeatureActionAsnRead()
*
**************************************************/
NLM_EXTERN 
RemoveDuplicateFeatureActionPtr LIBCALL
RemoveDuplicateFeatureActionAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   RemoveDuplicateFeatureActionPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* RemoveDuplicateFeatureAction ::= (self contained) */
      atp = AsnReadId(aip, amp, REMOVE_DUPLICATE_FEATURE_ACTION);
   } else {
      atp = AsnLinkType(orig, REMOVE_DUPLICATE_FEATURE_ACTION);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = RemoveDuplicateFeatureActionNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == REMOVE_DUPLICATE_FEATURE_ACTION_type) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> type = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == REMOVE_DUPLICATE_FEATURE_ACTION_ignore_partials) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> ignore_partials = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == REMOVE_DUPLICATE_FEATURE_ACTION_case_sensitive) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> case_sensitive = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == REMOVE_DUPLICATE_FEATURE_ACTION_remove_proteins) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> remove_proteins = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == REMOVE_DUPLICATE_FEATURE_ACTION_rd_constraint) {
      ptr -> rd_constraint = ConstraintChoiceSetAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = RemoveDuplicateFeatureActionFree(ptr);
   goto ret;
}



/**************************************************
*
*    RemoveDuplicateFeatureActionAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
RemoveDuplicateFeatureActionAsnWrite(RemoveDuplicateFeatureActionPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, REMOVE_DUPLICATE_FEATURE_ACTION);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> type;
   retval = AsnWrite(aip, REMOVE_DUPLICATE_FEATURE_ACTION_type,  &av);
   av.boolvalue = ptr -> ignore_partials;
   retval = AsnWrite(aip, REMOVE_DUPLICATE_FEATURE_ACTION_ignore_partials,  &av);
   av.boolvalue = ptr -> case_sensitive;
   retval = AsnWrite(aip, REMOVE_DUPLICATE_FEATURE_ACTION_case_sensitive,  &av);
   av.boolvalue = ptr -> remove_proteins;
   retval = AsnWrite(aip, REMOVE_DUPLICATE_FEATURE_ACTION_remove_proteins,  &av);
   if (ptr -> rd_constraint != NULL) {
      if ( ! ConstraintChoiceSetAsnWrite(ptr -> rd_constraint, aip, REMOVE_DUPLICATE_FEATURE_ACTION_rd_constraint)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    GeneXrefTypeNew()
*
**************************************************/
NLM_EXTERN 
GeneXrefTypePtr LIBCALL
GeneXrefTypeNew(void)
{
   GeneXrefTypePtr ptr = MemNew((size_t) sizeof(GeneXrefType));

   return ptr;

}


/**************************************************
*
*    GeneXrefTypeFree()
*
**************************************************/
NLM_EXTERN 
GeneXrefTypePtr LIBCALL
GeneXrefTypeFree(GeneXrefTypePtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   return MemFree(ptr);
}


/**************************************************
*
*    GeneXrefTypeAsnRead()
*
**************************************************/
NLM_EXTERN 
GeneXrefTypePtr LIBCALL
GeneXrefTypeAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   GeneXrefTypePtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* GeneXrefType ::= (self contained) */
      atp = AsnReadId(aip, amp, GENE_XREF_TYPE);
   } else {
      atp = AsnLinkType(orig, GENE_XREF_TYPE);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = GeneXrefTypeNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == GENE_XREF_TYPE_feature) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> feature = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GENE_XREF_TYPE_suppression) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> suppression = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GENE_XREF_TYPE_necessary) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> necessary = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = GeneXrefTypeFree(ptr);
   goto ret;
}



/**************************************************
*
*    GeneXrefTypeAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
GeneXrefTypeAsnWrite(GeneXrefTypePtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, GENE_XREF_TYPE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> feature;
   retval = AsnWrite(aip, GENE_XREF_TYPE_feature,  &av);
   av.intvalue = ptr -> suppression;
   retval = AsnWrite(aip, GENE_XREF_TYPE_suppression,  &av);
   av.intvalue = ptr -> necessary;
   retval = AsnWrite(aip, GENE_XREF_TYPE_necessary,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    XrefTypeFree()
*
**************************************************/
NLM_EXTERN 
XrefTypePtr LIBCALL
XrefTypeFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case XrefType_gene:
      GeneXrefTypeFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    XrefTypeAsnRead()
*
**************************************************/
NLM_EXTERN 
XrefTypePtr LIBCALL
XrefTypeAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* XrefType ::= (self contained) */
      atp = AsnReadId(aip, amp, XREF_TYPE);
   } else {
      atp = AsnLinkType(orig, XREF_TYPE);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == XREF_TYPE_gene) {
      choice = XrefType_gene;
      func = (AsnReadFunc) GeneXrefTypeAsnRead;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    XrefTypeAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
XrefTypeAsnWrite(XrefTypePtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, XREF_TYPE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case XrefType_gene:
      writetype = XREF_TYPE_gene;
      func = (AsnWriteFunc) GeneXrefTypeAsnWrite;
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    RemoveXrefsActionNew()
*
**************************************************/
NLM_EXTERN 
RemoveXrefsActionPtr LIBCALL
RemoveXrefsActionNew(void)
{
   RemoveXrefsActionPtr ptr = MemNew((size_t) sizeof(RemoveXrefsAction));

   return ptr;

}


/**************************************************
*
*    RemoveXrefsActionFree()
*
**************************************************/
NLM_EXTERN 
RemoveXrefsActionPtr LIBCALL
RemoveXrefsActionFree(RemoveXrefsActionPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   XrefTypeFree(ptr -> xref_type);
   ConstraintChoiceSetFree(ptr -> constraint);
   return MemFree(ptr);
}


/**************************************************
*
*    RemoveXrefsActionAsnRead()
*
**************************************************/
NLM_EXTERN 
RemoveXrefsActionPtr LIBCALL
RemoveXrefsActionAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   RemoveXrefsActionPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* RemoveXrefsAction ::= (self contained) */
      atp = AsnReadId(aip, amp, REMOVE_XREFS_ACTION);
   } else {
      atp = AsnLinkType(orig, REMOVE_XREFS_ACTION);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = RemoveXrefsActionNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == REMOVE_XREFS_ACTION_xref_type) {
      ptr -> xref_type = XrefTypeAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == REMOVE_XREFS_ACTION_constraint) {
      ptr -> constraint = ConstraintChoiceSetAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = RemoveXrefsActionFree(ptr);
   goto ret;
}



/**************************************************
*
*    RemoveXrefsActionAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
RemoveXrefsActionAsnWrite(RemoveXrefsActionPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, REMOVE_XREFS_ACTION);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> xref_type != NULL) {
      if ( ! XrefTypeAsnWrite(ptr -> xref_type, aip, REMOVE_XREFS_ACTION_xref_type)) {
         goto erret;
      }
   }
   if (ptr -> constraint != NULL) {
      if ( ! ConstraintChoiceSetAsnWrite(ptr -> constraint, aip, REMOVE_XREFS_ACTION_constraint)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    MakeGeneXrefActionNew()
*
**************************************************/
NLM_EXTERN 
MakeGeneXrefActionPtr LIBCALL
MakeGeneXrefActionNew(void)
{
   MakeGeneXrefActionPtr ptr = MemNew((size_t) sizeof(MakeGeneXrefAction));

   return ptr;

}


/**************************************************
*
*    MakeGeneXrefActionFree()
*
**************************************************/
NLM_EXTERN 
MakeGeneXrefActionPtr LIBCALL
MakeGeneXrefActionFree(MakeGeneXrefActionPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   ConstraintChoiceSetFree(ptr -> constraint);
   return MemFree(ptr);
}


/**************************************************
*
*    MakeGeneXrefActionAsnRead()
*
**************************************************/
NLM_EXTERN 
MakeGeneXrefActionPtr LIBCALL
MakeGeneXrefActionAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   MakeGeneXrefActionPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* MakeGeneXrefAction ::= (self contained) */
      atp = AsnReadId(aip, amp, MAKE_GENE_XREF_ACTION);
   } else {
      atp = AsnLinkType(orig, MAKE_GENE_XREF_ACTION);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = MakeGeneXrefActionNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == MAKE_GENE_XREF_ACTION_feature) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> feature = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == MAKE_GENE_XREF_ACTION_constraint) {
      ptr -> constraint = ConstraintChoiceSetAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = MakeGeneXrefActionFree(ptr);
   goto ret;
}



/**************************************************
*
*    MakeGeneXrefActionAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
MakeGeneXrefActionAsnWrite(MakeGeneXrefActionPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, MAKE_GENE_XREF_ACTION);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> feature;
   retval = AsnWrite(aip, MAKE_GENE_XREF_ACTION_feature,  &av);
   if (ptr -> constraint != NULL) {
      if ( ! ConstraintChoiceSetAsnWrite(ptr -> constraint, aip, MAKE_GENE_XREF_ACTION_constraint)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    AuthorFixActionNew()
*
**************************************************/
NLM_EXTERN 
AuthorFixActionPtr LIBCALL
AuthorFixActionNew(void)
{
   AuthorFixActionPtr ptr = MemNew((size_t) sizeof(AuthorFixAction));

   return ptr;

}


/**************************************************
*
*    AuthorFixActionFree()
*
**************************************************/
NLM_EXTERN 
AuthorFixActionPtr LIBCALL
AuthorFixActionFree(AuthorFixActionPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   ConstraintChoiceSetFree(ptr -> constraint);
   return MemFree(ptr);
}


/**************************************************
*
*    AuthorFixActionAsnRead()
*
**************************************************/
NLM_EXTERN 
AuthorFixActionPtr LIBCALL
AuthorFixActionAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   AuthorFixActionPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* AuthorFixAction ::= (self contained) */
      atp = AsnReadId(aip, amp, AUTHOR_FIX_ACTION);
   } else {
      atp = AsnLinkType(orig, AUTHOR_FIX_ACTION);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = AuthorFixActionNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == AUTHOR_FIX_ACTION_fix_type) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> fix_type = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == AUTHOR_FIX_ACTION_constraint) {
      ptr -> constraint = ConstraintChoiceSetAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = AuthorFixActionFree(ptr);
   goto ret;
}



/**************************************************
*
*    AuthorFixActionAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
AuthorFixActionAsnWrite(AuthorFixActionPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, AUTHOR_FIX_ACTION);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> fix_type;
   retval = AsnWrite(aip, AUTHOR_FIX_ACTION_fix_type,  &av);
   if (ptr -> constraint != NULL) {
      if ( ! ConstraintChoiceSetAsnWrite(ptr -> constraint, aip, AUTHOR_FIX_ACTION_constraint)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    UpdateSequencesActionNew()
*
**************************************************/
NLM_EXTERN 
UpdateSequencesActionPtr LIBCALL
UpdateSequencesActionNew(void)
{
   UpdateSequencesActionPtr ptr = MemNew((size_t) sizeof(UpdateSequencesAction));

   ptr -> add_cit_subs = 0;
   return ptr;

}


/**************************************************
*
*    UpdateSequencesActionFree()
*
**************************************************/
NLM_EXTERN 
UpdateSequencesActionPtr LIBCALL
UpdateSequencesActionFree(UpdateSequencesActionPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> filename);
   return MemFree(ptr);
}


/**************************************************
*
*    UpdateSequencesActionAsnRead()
*
**************************************************/
NLM_EXTERN 
UpdateSequencesActionPtr LIBCALL
UpdateSequencesActionAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   UpdateSequencesActionPtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* UpdateSequencesAction ::= (self contained) */
      atp = AsnReadId(aip, amp, UPDATE_SEQUENCES_ACTION);
   } else {
      atp = AsnLinkType(orig, UPDATE_SEQUENCES_ACTION);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = UpdateSequencesActionNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == UPDATE_SEQUENCES_ACTION_filename) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> filename = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == UPDATE_SEQUENCES_ACTION_add_cit_subs) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> add_cit_subs = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = UpdateSequencesActionFree(ptr);
   goto ret;
}



/**************************************************
*
*    UpdateSequencesActionAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
UpdateSequencesActionAsnWrite(UpdateSequencesActionPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, UPDATE_SEQUENCES_ACTION);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> filename != NULL) {
      av.ptrvalue = ptr -> filename;
      retval = AsnWrite(aip, UPDATE_SEQUENCES_ACTION_filename,  &av);
   }
   av.boolvalue = ptr -> add_cit_subs;
   retval = AsnWrite(aip, UPDATE_SEQUENCES_ACTION_add_cit_subs,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    MacroActionChoiceFree()
*
**************************************************/
NLM_EXTERN 
MacroActionChoicePtr LIBCALL
MacroActionChoiceFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case MacroActionChoice_aecr:
      AECRActionFree(anp -> data.ptrvalue);
      break;
   case MacroActionChoice_parse:
      ParseActionFree(anp -> data.ptrvalue);
      break;
   case MacroActionChoice_add_feature:
      ApplyFeatureActionFree(anp -> data.ptrvalue);
      break;
   case MacroActionChoice_remove_feature:
      RemoveFeatureActionFree(anp -> data.ptrvalue);
      break;
   case MacroActionChoice_convert_feature:
      ConvertFeatureActionFree(anp -> data.ptrvalue);
      break;
   case MacroActionChoice_edit_location:
      EditFeatureLocationActionFree(anp -> data.ptrvalue);
      break;
   case MacroActionChoice_remove_descriptor:
      RemoveDescriptorActionFree(anp -> data.ptrvalue);
      break;
   case MacroActionChoice_autodef:
      AutodefActionFree(anp -> data.ptrvalue);
      break;
   case MacroActionChoice_fix_pub_caps:
      FixPubCapsActionFree(anp -> data.ptrvalue);
      break;
   case MacroActionChoice_sort_fields:
      SortFieldsActionFree(anp -> data.ptrvalue);
      break;
   case MacroActionChoice_apply_molinfo_block:
      MolinfoBlockFree(anp -> data.ptrvalue);
      break;
   case MacroActionChoice_fix_caps:
      FixCapsActionFree(anp -> data.ptrvalue);
      break;
   case MacroActionChoice_fix_format:
      FixFormatActionFree(anp -> data.ptrvalue);
      break;
   case MacroActionChoice_remove_duplicate_features:
      RemoveDuplicateFeatureActionFree(anp -> data.ptrvalue);
      break;
   case MacroActionChoice_remove_xrefs:
      RemoveXrefsActionFree(anp -> data.ptrvalue);
      break;
   case MacroActionChoice_make_gene_xrefs:
      MakeGeneXrefActionFree(anp -> data.ptrvalue);
      break;
   case MacroActionChoice_fix_author:
      AuthorFixActionFree(anp -> data.ptrvalue);
      break;
   case MacroActionChoice_update_sequences:
      UpdateSequencesActionFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    MacroActionChoiceAsnRead()
*
**************************************************/
NLM_EXTERN 
MacroActionChoicePtr LIBCALL
MacroActionChoiceAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* MacroActionChoice ::= (self contained) */
      atp = AsnReadId(aip, amp, MACRO_ACTION_CHOICE);
   } else {
      atp = AsnLinkType(orig, MACRO_ACTION_CHOICE);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == MACRO_ACTION_CHOICE_aecr) {
      choice = MacroActionChoice_aecr;
      func = (AsnReadFunc) AECRActionAsnRead;
   }
   else if (atp == MACRO_ACTION_CHOICE_parse) {
      choice = MacroActionChoice_parse;
      func = (AsnReadFunc) ParseActionAsnRead;
   }
   else if (atp == MACRO_ACTION_CHOICE_add_feature) {
      choice = MacroActionChoice_add_feature;
      func = (AsnReadFunc) ApplyFeatureActionAsnRead;
   }
   else if (atp == MACRO_ACTION_CHOICE_remove_feature) {
      choice = MacroActionChoice_remove_feature;
      func = (AsnReadFunc) RemoveFeatureActionAsnRead;
   }
   else if (atp == MACRO_ACTION_CHOICE_convert_feature) {
      choice = MacroActionChoice_convert_feature;
      func = (AsnReadFunc) ConvertFeatureActionAsnRead;
   }
   else if (atp == MACRO_ACTION_CHOICE_edit_location) {
      choice = MacroActionChoice_edit_location;
      func = (AsnReadFunc) EditFeatureLocationActionAsnRead;
   }
   else if (atp == MACRO_ACTION_CHOICE_remove_descriptor) {
      choice = MacroActionChoice_remove_descriptor;
      func = (AsnReadFunc) RemoveDescriptorActionAsnRead;
   }
   else if (atp == MACRO_ACTION_CHOICE_autodef) {
      choice = MacroActionChoice_autodef;
      func = (AsnReadFunc) AutodefActionAsnRead;
   }
   else if (atp == MACRO_ACTION_CHOICE_removesets) {
      choice = MacroActionChoice_removesets;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == MACRO_ACTION_CHOICE_trim_junk_from_primer_seq) {
      choice = MacroActionChoice_trim_junk_from_primer_seq;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == MACRO_ACTION_CHOICE_trim_stop_from_complete_cds) {
      choice = MacroActionChoice_trim_stop_from_complete_cds;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == MACRO_ACTION_CHOICE_fix_usa_and_states) {
      choice = MacroActionChoice_fix_usa_and_states;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == MACRO_ACTION_CHOICE_synchronize_cds_partials) {
      choice = MacroActionChoice_synchronize_cds_partials;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == MACRO_ACTION_CHOICE_adjust_for_consensus_splice) {
      choice = MacroActionChoice_adjust_for_consensus_splice;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == MACRO_ACTION_CHOICE_fix_pub_caps) {
      choice = MacroActionChoice_fix_pub_caps;
      func = (AsnReadFunc) FixPubCapsActionAsnRead;
   }
   else if (atp == MACRO_ACTION_CHOICE_remove_seg_gaps) {
      choice = MacroActionChoice_remove_seg_gaps;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == MACRO_ACTION_CHOICE_sort_fields) {
      choice = MacroActionChoice_sort_fields;
      func = (AsnReadFunc) SortFieldsActionAsnRead;
   }
   else if (atp == MACRO_ACTION_CHOICE_apply_molinfo_block) {
      choice = MacroActionChoice_apply_molinfo_block;
      func = (AsnReadFunc) MolinfoBlockAsnRead;
   }
   else if (atp == MACRO_ACTION_CHOICE_fix_caps) {
      choice = MacroActionChoice_fix_caps;
      func = (AsnReadFunc) FixCapsActionAsnRead;
   }
   else if (atp == MACRO_ACTION_CHOICE_fix_format) {
      choice = MacroActionChoice_fix_format;
      func = (AsnReadFunc) FixFormatActionAsnRead;
   }
   else if (atp == MACRO_ACTION_CHOICE_fix_spell) {
      choice = MacroActionChoice_fix_spell;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == MACRO_ACTION_CHOICE_remove_duplicate_features) {
      choice = MacroActionChoice_remove_duplicate_features;
      func = (AsnReadFunc) RemoveDuplicateFeatureActionAsnRead;
   }
   else if (atp == MACRO_ACTION_CHOICE_remove_lineage_notes) {
      choice = MacroActionChoice_remove_lineage_notes;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == MACRO_ACTION_CHOICE_remove_xrefs) {
      choice = MacroActionChoice_remove_xrefs;
      func = (AsnReadFunc) RemoveXrefsActionAsnRead;
   }
   else if (atp == MACRO_ACTION_CHOICE_make_gene_xrefs) {
      choice = MacroActionChoice_make_gene_xrefs;
      func = (AsnReadFunc) MakeGeneXrefActionAsnRead;
   }
   else if (atp == MACRO_ACTION_CHOICE_make_bold_xrefs) {
      choice = MacroActionChoice_make_bold_xrefs;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == MACRO_ACTION_CHOICE_fix_author) {
      choice = MacroActionChoice_fix_author;
      func = (AsnReadFunc) AuthorFixActionAsnRead;
   }
   else if (atp == MACRO_ACTION_CHOICE_update_sequences) {
      choice = MacroActionChoice_update_sequences;
      func = (AsnReadFunc) UpdateSequencesActionAsnRead;
   }
   else if (atp == MACRO_ACTION_CHOICE_add_trans_splicing) {
      choice = MacroActionChoice_add_trans_splicing;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == MACRO_ACTION_CHOICE_remove_invalid_ecnumbers) {
      choice = MacroActionChoice_remove_invalid_ecnumbers;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    MacroActionChoiceAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
MacroActionChoiceAsnWrite(MacroActionChoicePtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, MACRO_ACTION_CHOICE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case MacroActionChoice_aecr:
      writetype = MACRO_ACTION_CHOICE_aecr;
      func = (AsnWriteFunc) AECRActionAsnWrite;
      break;
   case MacroActionChoice_parse:
      writetype = MACRO_ACTION_CHOICE_parse;
      func = (AsnWriteFunc) ParseActionAsnWrite;
      break;
   case MacroActionChoice_add_feature:
      writetype = MACRO_ACTION_CHOICE_add_feature;
      func = (AsnWriteFunc) ApplyFeatureActionAsnWrite;
      break;
   case MacroActionChoice_remove_feature:
      writetype = MACRO_ACTION_CHOICE_remove_feature;
      func = (AsnWriteFunc) RemoveFeatureActionAsnWrite;
      break;
   case MacroActionChoice_convert_feature:
      writetype = MACRO_ACTION_CHOICE_convert_feature;
      func = (AsnWriteFunc) ConvertFeatureActionAsnWrite;
      break;
   case MacroActionChoice_edit_location:
      writetype = MACRO_ACTION_CHOICE_edit_location;
      func = (AsnWriteFunc) EditFeatureLocationActionAsnWrite;
      break;
   case MacroActionChoice_remove_descriptor:
      writetype = MACRO_ACTION_CHOICE_remove_descriptor;
      func = (AsnWriteFunc) RemoveDescriptorActionAsnWrite;
      break;
   case MacroActionChoice_autodef:
      writetype = MACRO_ACTION_CHOICE_autodef;
      func = (AsnWriteFunc) AutodefActionAsnWrite;
      break;
   case MacroActionChoice_removesets:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, MACRO_ACTION_CHOICE_removesets, &av);
      break;
   case MacroActionChoice_trim_junk_from_primer_seq:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, MACRO_ACTION_CHOICE_trim_junk_from_primer_seq, &av);
      break;
   case MacroActionChoice_trim_stop_from_complete_cds:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, MACRO_ACTION_CHOICE_trim_stop_from_complete_cds, &av);
      break;
   case MacroActionChoice_fix_usa_and_states:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, MACRO_ACTION_CHOICE_fix_usa_and_states, &av);
      break;
   case MacroActionChoice_synchronize_cds_partials:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, MACRO_ACTION_CHOICE_synchronize_cds_partials, &av);
      break;
   case MacroActionChoice_adjust_for_consensus_splice:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, MACRO_ACTION_CHOICE_adjust_for_consensus_splice, &av);
      break;
   case MacroActionChoice_fix_pub_caps:
      writetype = MACRO_ACTION_CHOICE_fix_pub_caps;
      func = (AsnWriteFunc) FixPubCapsActionAsnWrite;
      break;
   case MacroActionChoice_remove_seg_gaps:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, MACRO_ACTION_CHOICE_remove_seg_gaps, &av);
      break;
   case MacroActionChoice_sort_fields:
      writetype = MACRO_ACTION_CHOICE_sort_fields;
      func = (AsnWriteFunc) SortFieldsActionAsnWrite;
      break;
   case MacroActionChoice_apply_molinfo_block:
      writetype = MACRO_ACTION_CHOICE_apply_molinfo_block;
      func = (AsnWriteFunc) MolinfoBlockAsnWrite;
      break;
   case MacroActionChoice_fix_caps:
      writetype = MACRO_ACTION_CHOICE_fix_caps;
      func = (AsnWriteFunc) FixCapsActionAsnWrite;
      break;
   case MacroActionChoice_fix_format:
      writetype = MACRO_ACTION_CHOICE_fix_format;
      func = (AsnWriteFunc) FixFormatActionAsnWrite;
      break;
   case MacroActionChoice_fix_spell:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, MACRO_ACTION_CHOICE_fix_spell, &av);
      break;
   case MacroActionChoice_remove_duplicate_features:
      writetype = MACRO_ACTION_CHOICE_remove_duplicate_features;
      func = (AsnWriteFunc) RemoveDuplicateFeatureActionAsnWrite;
      break;
   case MacroActionChoice_remove_lineage_notes:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, MACRO_ACTION_CHOICE_remove_lineage_notes, &av);
      break;
   case MacroActionChoice_remove_xrefs:
      writetype = MACRO_ACTION_CHOICE_remove_xrefs;
      func = (AsnWriteFunc) RemoveXrefsActionAsnWrite;
      break;
   case MacroActionChoice_make_gene_xrefs:
      writetype = MACRO_ACTION_CHOICE_make_gene_xrefs;
      func = (AsnWriteFunc) MakeGeneXrefActionAsnWrite;
      break;
   case MacroActionChoice_make_bold_xrefs:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, MACRO_ACTION_CHOICE_make_bold_xrefs, &av);
      break;
   case MacroActionChoice_fix_author:
      writetype = MACRO_ACTION_CHOICE_fix_author;
      func = (AsnWriteFunc) AuthorFixActionAsnWrite;
      break;
   case MacroActionChoice_update_sequences:
      writetype = MACRO_ACTION_CHOICE_update_sequences;
      func = (AsnWriteFunc) UpdateSequencesActionAsnWrite;
      break;
   case MacroActionChoice_add_trans_splicing:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, MACRO_ACTION_CHOICE_add_trans_splicing, &av);
      break;
   case MacroActionChoice_remove_invalid_ecnumbers:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, MACRO_ACTION_CHOICE_remove_invalid_ecnumbers, &av);
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    SearchFuncFree()
*
**************************************************/
NLM_EXTERN 
SearchFuncPtr LIBCALL
SearchFuncFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case SearchFunc_string_constraint:
      StringConstraintFree(anp -> data.ptrvalue);
      break;
   case SearchFunc_prefix_and_numbers:
      MemFree(anp -> data.ptrvalue);
      break;
   case SearchFunc_has_term:
      MemFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    SearchFuncAsnRead()
*
**************************************************/
NLM_EXTERN 
SearchFuncPtr LIBCALL
SearchFuncAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* SearchFunc ::= (self contained) */
      atp = AsnReadId(aip, amp, SEARCH_FUNC);
   } else {
      atp = AsnLinkType(orig, SEARCH_FUNC);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == SEARCH_FUNC_string_constraint) {
      choice = SearchFunc_string_constraint;
      func = (AsnReadFunc) StringConstraintAsnRead;
   }
   else if (atp == SEARCH_FUNC_contains_plural) {
      choice = SearchFunc_contains_plural;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == SEARCH_FUNC_n_or_more_brackets_or_parentheses) {
      choice = SearchFunc_n_or_more_brackets_or_parentheses;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == SEARCH_FUNC_three_numbers) {
      choice = SearchFunc_three_numbers;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == SEARCH_FUNC_underscore) {
      choice = SearchFunc_underscore;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == SEARCH_FUNC_prefix_and_numbers) {
      choice = SearchFunc_prefix_and_numbers;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.ptrvalue = av.ptrvalue;
   }
   else if (atp == SEARCH_FUNC_all_caps) {
      choice = SearchFunc_all_caps;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == SEARCH_FUNC_unbalanced_paren) {
      choice = SearchFunc_unbalanced_paren;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == SEARCH_FUNC_too_long) {
      choice = SearchFunc_too_long;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == SEARCH_FUNC_has_term) {
      choice = SearchFunc_has_term;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.ptrvalue = av.ptrvalue;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    SearchFuncAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
SearchFuncAsnWrite(SearchFuncPtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, SEARCH_FUNC);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case SearchFunc_string_constraint:
      writetype = SEARCH_FUNC_string_constraint;
      func = (AsnWriteFunc) StringConstraintAsnWrite;
      break;
   case SearchFunc_contains_plural:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, SEARCH_FUNC_contains_plural, &av);
      break;
   case SearchFunc_n_or_more_brackets_or_parentheses:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, SEARCH_FUNC_n_or_more_brackets_or_parentheses, &av);
      break;
   case SearchFunc_three_numbers:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, SEARCH_FUNC_three_numbers, &av);
      break;
   case SearchFunc_underscore:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, SEARCH_FUNC_underscore, &av);
      break;
   case SearchFunc_prefix_and_numbers:
      av.ptrvalue = anp->data.ptrvalue;
      retval = AsnWrite(aip, SEARCH_FUNC_prefix_and_numbers, &av);
      break;
   case SearchFunc_all_caps:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, SEARCH_FUNC_all_caps, &av);
      break;
   case SearchFunc_unbalanced_paren:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, SEARCH_FUNC_unbalanced_paren, &av);
      break;
   case SearchFunc_too_long:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, SEARCH_FUNC_too_long, &av);
      break;
   case SearchFunc_has_term:
      av.ptrvalue = anp->data.ptrvalue;
      retval = AsnWrite(aip, SEARCH_FUNC_has_term, &av);
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    SimpleReplaceNew()
*
**************************************************/
NLM_EXTERN 
SimpleReplacePtr LIBCALL
SimpleReplaceNew(void)
{
   SimpleReplacePtr ptr = MemNew((size_t) sizeof(SimpleReplace));

   ptr -> whole_string = 0;
   ptr -> weasel_to_putative = 0;
   return ptr;

}


/**************************************************
*
*    SimpleReplaceFree()
*
**************************************************/
NLM_EXTERN 
SimpleReplacePtr LIBCALL
SimpleReplaceFree(SimpleReplacePtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> replace);
   return MemFree(ptr);
}


/**************************************************
*
*    SimpleReplaceAsnRead()
*
**************************************************/
NLM_EXTERN 
SimpleReplacePtr LIBCALL
SimpleReplaceAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   SimpleReplacePtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* SimpleReplace ::= (self contained) */
      atp = AsnReadId(aip, amp, SIMPLE_REPLACE);
   } else {
      atp = AsnLinkType(orig, SIMPLE_REPLACE);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = SimpleReplaceNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == SIMPLE_REPLACE_replace) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> replace = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SIMPLE_REPLACE_whole_string) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> whole_string = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SIMPLE_REPLACE_weasel_to_putative) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> weasel_to_putative = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = SimpleReplaceFree(ptr);
   goto ret;
}



/**************************************************
*
*    SimpleReplaceAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
SimpleReplaceAsnWrite(SimpleReplacePtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, SIMPLE_REPLACE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> replace != NULL) {
      av.ptrvalue = ptr -> replace;
      retval = AsnWrite(aip, SIMPLE_REPLACE_replace,  &av);
   }
   av.boolvalue = ptr -> whole_string;
   retval = AsnWrite(aip, SIMPLE_REPLACE_whole_string,  &av);
   av.boolvalue = ptr -> weasel_to_putative;
   retval = AsnWrite(aip, SIMPLE_REPLACE_weasel_to_putative,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ReplaceFuncFree()
*
**************************************************/
NLM_EXTERN 
ReplaceFuncPtr LIBCALL
ReplaceFuncFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case ReplaceFunc_simple_replace:
      SimpleReplaceFree(anp -> data.ptrvalue);
      break;
   case ReplaceFunc_haem_replace:
      MemFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    ReplaceFuncAsnRead()
*
**************************************************/
NLM_EXTERN 
ReplaceFuncPtr LIBCALL
ReplaceFuncAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ReplaceFunc ::= (self contained) */
      atp = AsnReadId(aip, amp, REPLACE_FUNC);
   } else {
      atp = AsnLinkType(orig, REPLACE_FUNC);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == REPLACE_FUNC_simple_replace) {
      choice = ReplaceFunc_simple_replace;
      func = (AsnReadFunc) SimpleReplaceAsnRead;
   }
   else if (atp == REPLACE_FUNC_haem_replace) {
      choice = ReplaceFunc_haem_replace;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.ptrvalue = av.ptrvalue;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    ReplaceFuncAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ReplaceFuncAsnWrite(ReplaceFuncPtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, REPLACE_FUNC);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case ReplaceFunc_simple_replace:
      writetype = REPLACE_FUNC_simple_replace;
      func = (AsnWriteFunc) SimpleReplaceAsnWrite;
      break;
   case ReplaceFunc_haem_replace:
      av.ptrvalue = anp->data.ptrvalue;
      retval = AsnWrite(aip, REPLACE_FUNC_haem_replace, &av);
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    ReplaceRuleNew()
*
**************************************************/
NLM_EXTERN 
ReplaceRulePtr LIBCALL
ReplaceRuleNew(void)
{
   ReplaceRulePtr ptr = MemNew((size_t) sizeof(ReplaceRule));

   ptr -> move_to_note = 0;
   return ptr;

}


/**************************************************
*
*    ReplaceRuleFree()
*
**************************************************/
NLM_EXTERN 
ReplaceRulePtr LIBCALL
ReplaceRuleFree(ReplaceRulePtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   ReplaceFuncFree(ptr -> replace_func);
   return MemFree(ptr);
}


/**************************************************
*
*    ReplaceRuleAsnRead()
*
**************************************************/
NLM_EXTERN 
ReplaceRulePtr LIBCALL
ReplaceRuleAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ReplaceRulePtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ReplaceRule ::= (self contained) */
      atp = AsnReadId(aip, amp, REPLACE_RULE);
   } else {
      atp = AsnLinkType(orig, REPLACE_RULE);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ReplaceRuleNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == REPLACE_RULE_replace_func) {
      ptr -> replace_func = ReplaceFuncAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == REPLACE_RULE_move_to_note) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> move_to_note = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ReplaceRuleFree(ptr);
   goto ret;
}



/**************************************************
*
*    ReplaceRuleAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ReplaceRuleAsnWrite(ReplaceRulePtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, REPLACE_RULE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> replace_func != NULL) {
      if ( ! ReplaceFuncAsnWrite(ptr -> replace_func, aip, REPLACE_RULE_replace_func)) {
         goto erret;
      }
   }
   av.boolvalue = ptr -> move_to_note;
   retval = AsnWrite(aip, REPLACE_RULE_move_to_note,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    SuspectRuleNew()
*
**************************************************/
NLM_EXTERN 
SuspectRulePtr LIBCALL
SuspectRuleNew(void)
{
   SuspectRulePtr ptr = MemNew((size_t) sizeof(SuspectRule));

   ptr -> rule_type = 0;
   return ptr;

}


/**************************************************
*
*    SuspectRuleFree()
*
**************************************************/
NLM_EXTERN 
SuspectRulePtr LIBCALL
SuspectRuleFree(SuspectRulePtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   SearchFuncFree(ptr -> find);
   SearchFuncFree(ptr -> except);
   ConstraintChoiceSetFree(ptr -> feat_constraint);
   ReplaceRuleFree(ptr -> replace);
   MemFree(ptr -> description);
   return MemFree(ptr);
}


/**************************************************
*
*    SuspectRuleAsnRead()
*
**************************************************/
NLM_EXTERN 
SuspectRulePtr LIBCALL
SuspectRuleAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   SuspectRulePtr ptr;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* SuspectRule ::= (self contained) */
      atp = AsnReadId(aip, amp, SUSPECT_RULE);
   } else {
      atp = AsnLinkType(orig, SUSPECT_RULE);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = SuspectRuleNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == SUSPECT_RULE_find) {
      ptr -> find = SearchFuncAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SUSPECT_RULE_except) {
      ptr -> except = SearchFuncAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SUSPECT_RULE_feat_constraint) {
      ptr -> feat_constraint = ConstraintChoiceSetAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SUSPECT_RULE_rule_type) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> rule_type = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SUSPECT_RULE_replace) {
      ptr -> replace = ReplaceRuleAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SUSPECT_RULE_description) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> description = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = SuspectRuleFree(ptr);
   goto ret;
}



/**************************************************
*
*    SuspectRuleAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
SuspectRuleAsnWrite(SuspectRulePtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmacroAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, SUSPECT_RULE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> find != NULL) {
      if ( ! SearchFuncAsnWrite(ptr -> find, aip, SUSPECT_RULE_find)) {
         goto erret;
      }
   }
   if (ptr -> except != NULL) {
      if ( ! SearchFuncAsnWrite(ptr -> except, aip, SUSPECT_RULE_except)) {
         goto erret;
      }
   }
   if (ptr -> feat_constraint != NULL) {
      if ( ! ConstraintChoiceSetAsnWrite(ptr -> feat_constraint, aip, SUSPECT_RULE_feat_constraint)) {
         goto erret;
      }
   }
   av.intvalue = ptr -> rule_type;
   retval = AsnWrite(aip, SUSPECT_RULE_rule_type,  &av);
   if (ptr -> replace != NULL) {
      if ( ! ReplaceRuleAsnWrite(ptr -> replace, aip, SUSPECT_RULE_replace)) {
         goto erret;
      }
   }
   if (ptr -> description != NULL) {
      av.ptrvalue = ptr -> description;
      retval = AsnWrite(aip, SUSPECT_RULE_description,  &av);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}

