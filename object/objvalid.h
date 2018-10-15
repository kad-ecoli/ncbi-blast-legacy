#ifndef _objvalid_ 
#define _objvalid_ 

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif


#ifdef __cplusplus
extern "C" { /* } */
#endif


/**************************************************
*
*    Generated objects for Module NCBI-Structured-comment-validation
*    Generated using ASNCODE Revision: 6.16 at Dec 11, 2009  7:16 AM
*
**************************************************/

NLM_EXTERN Boolean LIBCALL
objvalidAsnLoad PROTO((void));
/* following #defines are for enumerated type, not used by object loaders */
#define Severity_level_none 0
#define Severity_level_info 1
#define Severity_level_warning 2
#define Severity_level_error 3
#define Severity_level_reject 4
#define Severity_level_fatal 5



/**************************************************
*
*    FieldRule
*
**************************************************/
typedef struct struct_Field_rule {
   struct struct_Field_rule PNTR next;
   CharPtr   field_name;
   CharPtr   match_expression;
   Uint1   required;
   Uint2   severity;
} FieldRule, PNTR FieldRulePtr;


NLM_EXTERN FieldRulePtr LIBCALL FieldRuleFree PROTO ((FieldRulePtr ));
NLM_EXTERN FieldRulePtr LIBCALL FieldRuleNew PROTO (( void ));
NLM_EXTERN FieldRulePtr LIBCALL FieldRuleAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL FieldRuleAsnWrite PROTO (( FieldRulePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    FieldSet
*
**************************************************/
typedef struct struct_Field_rule FieldSet;
typedef struct struct_Field_rule PNTR FieldSetPtr;
#define FieldSetNew() Field_ruleNew() 

#ifdef NLM_GENERATED_CODE_PROTO

NLM_EXTERN FieldSetPtr LIBCALL FieldSetFree PROTO ((FieldSetPtr ));
NLM_EXTERN FieldSetPtr LIBCALL FieldSetNew PROTO (( void ));
NLM_EXTERN FieldSetPtr LIBCALL FieldSetAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL FieldSetAsnWrite PROTO (( FieldSetPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */



/**************************************************
*
*    DependentFieldRule
*
**************************************************/
typedef struct struct_Dependent_field_rule {
   struct struct_Dependent_field_rule PNTR next;
   CharPtr   match_name;
   CharPtr   value_constraint;
   struct struct_Field_rule PNTR   other_fields;
   struct struct_Field_rule PNTR   disallowed_fields;
} DependentFieldRule, PNTR DependentFieldRulePtr;


NLM_EXTERN DependentFieldRulePtr LIBCALL DependentFieldRuleFree PROTO ((DependentFieldRulePtr ));
NLM_EXTERN DependentFieldRulePtr LIBCALL DependentFieldRuleNew PROTO (( void ));
NLM_EXTERN DependentFieldRulePtr LIBCALL DependentFieldRuleAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL DependentFieldRuleAsnWrite PROTO (( DependentFieldRulePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    DependentFieldSet
*
**************************************************/
typedef struct struct_Dependent_field_rule DependentFieldSet;
typedef struct struct_Dependent_field_rule PNTR DependentFieldSetPtr;
#define DependentFieldSetNew() Dependent_field_ruleNew() 

#ifdef NLM_GENERATED_CODE_PROTO

NLM_EXTERN DependentFieldSetPtr LIBCALL DependentFieldSetFree PROTO ((DependentFieldSetPtr ));
NLM_EXTERN DependentFieldSetPtr LIBCALL DependentFieldSetNew PROTO (( void ));
NLM_EXTERN DependentFieldSetPtr LIBCALL DependentFieldSetAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL DependentFieldSetAsnWrite PROTO (( DependentFieldSetPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */



/**************************************************
*
*    CommentRule
*
**************************************************/
typedef struct struct_Comment_rule {
   struct struct_Comment_rule PNTR next;
   CharPtr   prefix;
   Uint1   updated;
   struct struct_Field_rule PNTR   fields;
   Uint1   require_order;
   Uint1   allow_unlisted;
   struct struct_Dependent_field_rule PNTR   dependent_rules;
} CommentRule, PNTR CommentRulePtr;


NLM_EXTERN CommentRulePtr LIBCALL CommentRuleFree PROTO ((CommentRulePtr ));
NLM_EXTERN CommentRulePtr LIBCALL CommentRuleNew PROTO (( void ));
NLM_EXTERN CommentRulePtr LIBCALL CommentRuleAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL CommentRuleAsnWrite PROTO (( CommentRulePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    CommentSet
*
**************************************************/
typedef struct struct_Comment_rule CommentSet;
typedef struct struct_Comment_rule PNTR CommentSetPtr;
#define CommentSetNew() Comment_ruleNew() 

#ifdef NLM_GENERATED_CODE_PROTO

NLM_EXTERN CommentSetPtr LIBCALL CommentSetFree PROTO ((CommentSetPtr ));
NLM_EXTERN CommentSetPtr LIBCALL CommentSetNew PROTO (( void ));
NLM_EXTERN CommentSetPtr LIBCALL CommentSetAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL CommentSetAsnWrite PROTO (( CommentSetPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */

#ifdef __cplusplus
/* { */ }
#endif

#endif /* _objvalid_ */

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

