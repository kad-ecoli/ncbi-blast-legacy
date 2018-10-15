#include <asn.h>

#define NLM_GENERATED_CODE_PROTO

#include <objvalid.h>

static Boolean loaded = FALSE;

#include <asnvalid.h>

#ifndef NLM_EXTERN_LOADS
#define NLM_EXTERN_LOADS {}
#endif

NLM_EXTERN Boolean LIBCALL
objvalidAsnLoad(void)
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
*    Generated object loaders for Module NCBI-Structured-comment-validation
*    Generated using ASNCODE Revision: 6.16 at Dec 11, 2009  7:16 AM
*
**************************************************/


/**************************************************
*
*    FieldRuleNew()
*
**************************************************/
NLM_EXTERN 
FieldRulePtr LIBCALL
FieldRuleNew(void)
{
   FieldRulePtr ptr = MemNew((size_t) sizeof(FieldRule));

   ptr -> required = 0;
   ptr -> severity = 2;
   return ptr;

}


/**************************************************
*
*    FieldRuleFree()
*
**************************************************/
NLM_EXTERN 
FieldRulePtr LIBCALL
FieldRuleFree(FieldRulePtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> field_name);
   MemFree(ptr -> match_expression);
   return MemFree(ptr);
}


/**************************************************
*
*    FieldRuleAsnRead()
*
**************************************************/
NLM_EXTERN 
FieldRulePtr LIBCALL
FieldRuleAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   FieldRulePtr ptr;

   if (! loaded)
   {
      if (! objvalidAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* FieldRule ::= (self contained) */
      atp = AsnReadId(aip, amp, FIELD_RULE);
   } else {
      atp = AsnLinkType(orig, FIELD_RULE);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = FieldRuleNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == FIELD_RULE_field_name) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> field_name = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == FIELD_RULE_match_expression) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> match_expression = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == FIELD_RULE_required) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> required = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == FIELD_RULE_severity) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> severity = av.intvalue;
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
   ptr = FieldRuleFree(ptr);
   goto ret;
}



/**************************************************
*
*    FieldRuleAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
FieldRuleAsnWrite(FieldRulePtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objvalidAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, FIELD_RULE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> field_name != NULL) {
      av.ptrvalue = ptr -> field_name;
      retval = AsnWrite(aip, FIELD_RULE_field_name,  &av);
   }
   if (ptr -> match_expression != NULL) {
      av.ptrvalue = ptr -> match_expression;
      retval = AsnWrite(aip, FIELD_RULE_match_expression,  &av);
   }
   av.boolvalue = ptr -> required;
   retval = AsnWrite(aip, FIELD_RULE_required,  &av);
   av.intvalue = ptr -> severity;
   retval = AsnWrite(aip, FIELD_RULE_severity,  &av);
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
*    FieldSetFree()
*
**************************************************/
NLM_EXTERN 
FieldSetPtr LIBCALL
FieldSetFree(FieldSetPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericUserSeqOfFree(ptr,  (AsnOptFreeFunc) FieldRuleFree);
   return NULL;
}


/**************************************************
*
*    FieldSetAsnRead()
*
**************************************************/
NLM_EXTERN 
FieldSetPtr LIBCALL
FieldSetAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   FieldSetPtr ptr;

   if (! loaded)
   {
      if (! objvalidAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* FieldSet ::= (self contained) */
      atp = AsnReadId(aip, amp, FIELD_SET);
   } else {
      atp = AsnLinkType(orig, FIELD_SET);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   func = NULL;

   ptr  = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) FieldRuleAsnRead, (AsnOptFreeFunc) FieldRuleFree);
   if (isError && ptr  == NULL) {
      goto erret;
   }



ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = FieldSetFree(ptr);
   goto ret;
}



/**************************************************
*
*    FieldSetAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
FieldSetAsnWrite(FieldSetPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objvalidAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, FIELD_SET);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   retval = AsnGenericUserSeqOfAsnWrite(ptr , (AsnWriteFunc) FieldRuleAsnWrite, aip, atp, FIELD_SET_E);
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    DependentFieldRuleNew()
*
**************************************************/
NLM_EXTERN 
DependentFieldRulePtr LIBCALL
DependentFieldRuleNew(void)
{
   DependentFieldRulePtr ptr = MemNew((size_t) sizeof(DependentFieldRule));

   return ptr;

}


/**************************************************
*
*    DependentFieldRuleFree()
*
**************************************************/
NLM_EXTERN 
DependentFieldRulePtr LIBCALL
DependentFieldRuleFree(DependentFieldRulePtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> match_name);
   MemFree(ptr -> value_constraint);
   FieldSetFree(ptr -> other_fields);
   FieldSetFree(ptr -> disallowed_fields);
   return MemFree(ptr);
}


/**************************************************
*
*    DependentFieldRuleAsnRead()
*
**************************************************/
NLM_EXTERN 
DependentFieldRulePtr LIBCALL
DependentFieldRuleAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   DependentFieldRulePtr ptr;

   if (! loaded)
   {
      if (! objvalidAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* DependentFieldRule ::= (self contained) */
      atp = AsnReadId(aip, amp, DEPENDENT_FIELD_RULE);
   } else {
      atp = AsnLinkType(orig, DEPENDENT_FIELD_RULE);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = DependentFieldRuleNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == DEPENDENT_FIELD_RULE_match_name) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> match_name = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == FIELD_RULE_value_constraint) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> value_constraint = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == FIELD_RULE_other_fields) {
      ptr -> other_fields = FieldSetAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == FIELD_RULE_disallowed_fields) {
      ptr -> disallowed_fields = FieldSetAsnRead(aip, atp);
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
   ptr = DependentFieldRuleFree(ptr);
   goto ret;
}



/**************************************************
*
*    DependentFieldRuleAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
DependentFieldRuleAsnWrite(DependentFieldRulePtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objvalidAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, DEPENDENT_FIELD_RULE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> match_name != NULL) {
      av.ptrvalue = ptr -> match_name;
      retval = AsnWrite(aip, DEPENDENT_FIELD_RULE_match_name,  &av);
   }
   if (ptr -> value_constraint != NULL) {
      av.ptrvalue = ptr -> value_constraint;
      retval = AsnWrite(aip, FIELD_RULE_value_constraint,  &av);
   }
   if (ptr -> other_fields != NULL) {
      if ( ! FieldSetAsnWrite(ptr -> other_fields, aip, FIELD_RULE_other_fields)) {
         goto erret;
      }
   }
   if (ptr -> disallowed_fields != NULL) {
      if ( ! FieldSetAsnWrite(ptr -> disallowed_fields, aip, FIELD_RULE_disallowed_fields)) {
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
*    DependentFieldSetFree()
*
**************************************************/
NLM_EXTERN 
DependentFieldSetPtr LIBCALL
DependentFieldSetFree(DependentFieldSetPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericUserSeqOfFree(ptr,  (AsnOptFreeFunc) DependentFieldRuleFree);
   return NULL;
}


/**************************************************
*
*    DependentFieldSetAsnRead()
*
**************************************************/
NLM_EXTERN 
DependentFieldSetPtr LIBCALL
DependentFieldSetAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   DependentFieldSetPtr ptr;

   if (! loaded)
   {
      if (! objvalidAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* DependentFieldSet ::= (self contained) */
      atp = AsnReadId(aip, amp, DEPENDENT_FIELD_SET);
   } else {
      atp = AsnLinkType(orig, DEPENDENT_FIELD_SET);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   func = NULL;

   ptr  = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) DependentFieldRuleAsnRead, (AsnOptFreeFunc) DependentFieldRuleFree);
   if (isError && ptr  == NULL) {
      goto erret;
   }



ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = DependentFieldSetFree(ptr);
   goto ret;
}



/**************************************************
*
*    DependentFieldSetAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
DependentFieldSetAsnWrite(DependentFieldSetPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objvalidAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, DEPENDENT_FIELD_SET);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   retval = AsnGenericUserSeqOfAsnWrite(ptr , (AsnWriteFunc) DependentFieldRuleAsnWrite, aip, atp, DEPENDENT_FIELD_SET_E);
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    CommentRuleNew()
*
**************************************************/
NLM_EXTERN 
CommentRulePtr LIBCALL
CommentRuleNew(void)
{
   CommentRulePtr ptr = MemNew((size_t) sizeof(CommentRule));

   ptr -> updated = 0;
   ptr -> require_order = 1;
   ptr -> allow_unlisted = 0;
   return ptr;

}


/**************************************************
*
*    CommentRuleFree()
*
**************************************************/
NLM_EXTERN 
CommentRulePtr LIBCALL
CommentRuleFree(CommentRulePtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> prefix);
   FieldSetFree(ptr -> fields);
   DependentFieldSetFree(ptr -> dependent_rules);
   return MemFree(ptr);
}


/**************************************************
*
*    CommentRuleAsnRead()
*
**************************************************/
NLM_EXTERN 
CommentRulePtr LIBCALL
CommentRuleAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   CommentRulePtr ptr;

   if (! loaded)
   {
      if (! objvalidAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* CommentRule ::= (self contained) */
      atp = AsnReadId(aip, amp, COMMENT_RULE);
   } else {
      atp = AsnLinkType(orig, COMMENT_RULE);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = CommentRuleNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == COMMENT_RULE_prefix) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> prefix = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == COMMENT_RULE_updated) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> updated = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == COMMENT_RULE_fields) {
      ptr -> fields = FieldSetAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == COMMENT_RULE_require_order) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> require_order = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == COMMENT_RULE_allow_unlisted) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> allow_unlisted = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == COMMENT_RULE_dependent_rules) {
      ptr -> dependent_rules = DependentFieldSetAsnRead(aip, atp);
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
   ptr = CommentRuleFree(ptr);
   goto ret;
}



/**************************************************
*
*    CommentRuleAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
CommentRuleAsnWrite(CommentRulePtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objvalidAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, COMMENT_RULE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> prefix != NULL) {
      av.ptrvalue = ptr -> prefix;
      retval = AsnWrite(aip, COMMENT_RULE_prefix,  &av);
   }
   av.boolvalue = ptr -> updated;
   retval = AsnWrite(aip, COMMENT_RULE_updated,  &av);
   if (ptr -> fields != NULL) {
      if ( ! FieldSetAsnWrite(ptr -> fields, aip, COMMENT_RULE_fields)) {
         goto erret;
      }
   }
   av.boolvalue = ptr -> require_order;
   retval = AsnWrite(aip, COMMENT_RULE_require_order,  &av);
   av.boolvalue = ptr -> allow_unlisted;
   retval = AsnWrite(aip, COMMENT_RULE_allow_unlisted,  &av);
   if (ptr -> dependent_rules != NULL) {
      if ( ! DependentFieldSetAsnWrite(ptr -> dependent_rules, aip, COMMENT_RULE_dependent_rules)) {
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
*    CommentSetFree()
*
**************************************************/
NLM_EXTERN 
CommentSetPtr LIBCALL
CommentSetFree(CommentSetPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericUserSeqOfFree(ptr,  (AsnOptFreeFunc) CommentRuleFree);
   return NULL;
}


/**************************************************
*
*    CommentSetAsnRead()
*
**************************************************/
NLM_EXTERN 
CommentSetPtr LIBCALL
CommentSetAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   CommentSetPtr ptr;

   if (! loaded)
   {
      if (! objvalidAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* CommentSet ::= (self contained) */
      atp = AsnReadId(aip, amp, COMMENT_SET);
   } else {
      atp = AsnLinkType(orig, COMMENT_SET);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   func = NULL;

   ptr  = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) CommentRuleAsnRead, (AsnOptFreeFunc) CommentRuleFree);
   if (isError && ptr  == NULL) {
      goto erret;
   }



ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = CommentSetFree(ptr);
   goto ret;
}



/**************************************************
*
*    CommentSetAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
CommentSetAsnWrite(CommentSetPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objvalidAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, COMMENT_SET);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   retval = AsnGenericUserSeqOfAsnWrite(ptr , (AsnWriteFunc) CommentRuleAsnWrite, aip, atp, COMMENT_SET_E);
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}

