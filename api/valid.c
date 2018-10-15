/*  valid.c
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================
*
* File Name:  valid.c
*
* Author:  James Ostell
*
* Version Creation Date: 1/1/94
*
* $Revision: 6.1569 $
*
* File Description:  Sequence editing utilities
*
* Modifications:
* --------------------------------------------------------------------------
* Date       Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
*
*
* ==========================================================================
*/

static char    *this_module = "valid";

#define THIS_MODULE this_module

static char    *this_file = __FILE__;

#define THIS_FILE this_file

#include <ncbi.h>
#include <objfdef.h>
#include <valid.h>
#include <validerr.h>
#include <sqnutils.h>
#include <gbftdef.h>
#include <gbfeat.h>
#include <objsub.h>
#include <asn2gnbi.h>
#include <explore.h>
#include <gather.h>
#include <subutil.h>
#include <tofasta.h>
#include <findrepl.h>
#define NLM_GENERATED_CODE_PROTO
#include <objmacro.h>
#include <macroapi.h>
#include <objvalid.h>
#include <valapi.h>
#include "ecnum_specific.inc"
#include "ecnum_ambiguous.inc"
#include "ecnum_deleted.inc"
#include "ecnum_replaced.inc"

/*****************************************************************************
*
*   NOTE: look at all the ValidErr calls with severity=0. Some should be
*   bumped up later. Look also for string "PARSER"
*
*****************************************************************************/



#ifdef VAR_ARGS
#include <varargs.h>
#else
#include <stdarg.h>
#endif

static ValidStructPtr globalvsp;        /* for spell checker */

NLM_EXTERN void CDECL ValidErr VPROTO ((ValidStructPtr vsp, int severity, int code1, int code2, const char *fmt, ...));
static void     ValidateBioseqInst (GatherContextPtr gcp);
static void     ValidateBioseqContext (GatherContextPtr gcp);
static void     ValidateBioseqSet (GatherContextPtr gcp);
static void     ValidateGraphsOnBioseq (GatherContextPtr gcp);
static void     ValidateBioseqHist (GatherContextPtr gcp);
static void     SpellCheckSeqDescr (GatherContextPtr gcp);
NLM_EXTERN void CdTransCheck (ValidStructPtr vsp, SeqFeatPtr sfp);
NLM_EXTERN void MrnaTransCheck (ValidStructPtr vsp, SeqFeatPtr sfp);
NLM_EXTERN void ValidateSeqFeat (GatherContextPtr gcp);
NLM_EXTERN void ValidateSeqLoc (ValidStructPtr vsp, SeqLocPtr slp, CharPtr prefix);
NLM_EXTERN Boolean PatchBadSequence (BioseqPtr bsp);
NLM_EXTERN CharPtr FindIDForEntry (SeqEntryPtr sep, CharPtr buf);
NLM_EXTERN void SpellCheckSeqFeat (GatherContextPtr gcp);
NLM_EXTERN void SpellCheckString (ValidStructPtr vsp, CharPtr str);
NLM_EXTERN void SpliceCheck (ValidStructPtr vsp, SeqFeatPtr sfp);
static void     CdConflictCheck (ValidStructPtr vsp, SeqFeatPtr sfp);
static void     SpliceCheckEx (ValidStructPtr vsp, SeqFeatPtr sfp, Boolean checkAll);
static void     CdsProductIdCheck (ValidStructPtr vsp, SeqFeatPtr sfp);
static void     ValidateBioSource (ValidStructPtr vsp, GatherContextPtr gcp, BioSourcePtr biop, SeqFeatPtr sfp, ValNodePtr sdp);
static void     ValidatePubdesc (ValidStructPtr vsp, GatherContextPtr gcp, PubdescPtr pdp);
static void     LookForMultiplePubs (ValidStructPtr vsp, GatherContextPtr gcp, SeqDescrPtr sdp);
static void     ValidateSfpCit (ValidStructPtr vsp, GatherContextPtr gcp, SeqFeatPtr sfp);
static void     ValidateAffil (ValidStructPtr vsp, AffilPtr ap);
static TextFsaPtr GetSpecificECNumberFSA (void);
static TextFsaPtr GetAmbiguousECNumberFSA (void);
static TextFsaPtr GetDeletedECNumberFSA (void);
static TextFsaPtr GetReplacedECNumberFSA (void);
static void ValidateCitSub (ValidStructPtr vsp, CitSubPtr csp);

static Boolean HasFeatId(SeqFeatPtr sfp, Int4 num)
{
  Boolean rval = FALSE;
  ObjectIdPtr oip;

  if (sfp == NULL) {
    return FALSE;
  }
  if (sfp->id.choice == 3) {
    oip = (ObjectIdPtr) sfp->id.value.ptrvalue;
    if (oip->id == num) {
      rval = TRUE;
    }
  }
  return rval;
}


/* alignment validator */
NLM_EXTERN Boolean ValidateSeqAlignWithinValidator (ValidStructPtr vsp, SeqEntryPtr sep, Boolean find_remote_bsp, Boolean do_hist_assembly);

static ValNodePtr genetic_code_name_list = NULL;

/*****************************************************************************
*
*   Perform Validation Checks on a SeqEntry
*
*****************************************************************************/

NLM_EXTERN void ValidStructClear (ValidStructPtr vsp)
{                               /* 0 out a ValidStruct */
  CharPtr         errbuf;
  Int2            cutoff;
  Boolean         patch_seq;
  SpellCheckFunc  spellfunc;
  SpellCallBackFunc spellcallback;
  Boolean         onlyspell;
  Boolean         justwarnonspell;
  Boolean         useSeqMgrIndexes;
  Boolean         suppressContext;
  Boolean         validateAlignments;
  Boolean         farIDsInAlignments;
  Boolean         alignFindRemoteBsp;
  Boolean         doSeqHistAssembly;
  Boolean         alwaysRequireIsoJTA;
  Boolean         farFetchCDSproducts;
  Boolean         farFetchMRNAproducts;
  Boolean         locusTagGeneralMatch;
  Boolean         validateIDSet;
  Boolean         seqSubmitParent;
  Boolean         justShowAccession;
  Boolean         ignoreExceptions;
  Boolean         validateExons;
  Boolean         inferenceAccnCheck;
  Boolean         testLatLonSubregion;
  Boolean         strictLatLonCountry;
  Boolean         rubiscoTest;
  Boolean         indexerVersion;
  Boolean         disableSuppression;
  Int2            validationLimit;
  ValidErrorFunc  errfunc;
  Pointer         userdata;
  Boolean         convertGiToAccn;
  TextFsaPtr      sourceQualTags;
  TextFsaPtr      modifiedBases;
  TextFsaPtr      sgmlStrings;
  Boolean         is_htg_in_sep;
  Boolean         is_barcode_sep;
  Boolean         is_refseq_in_sep;
  Boolean         is_gpipe_in_sep;
  Boolean         is_gps_in_sep;
  Boolean         is_small_genome_set;
  Boolean         is_embl_ddbj_in_sep;
  Boolean         is_old_gb_in_sep;
  Boolean         is_patent_in_sep;
  Boolean         other_sets_in_sep;
  Boolean         is_insd_in_sep;
  Boolean         only_lcl_gnl_in_sep;
  Boolean         has_gnl_prot_sep;
  Boolean         bsp_genomic_in_sep;
  Boolean         is_smupd_in_sep;
  Boolean         feat_loc_has_gi;
  Boolean         feat_prod_has_gi;
  Boolean         has_multi_int_genes;
  Boolean         has_seg_bioseqs;
  Boolean         far_fetch_failure;

  if (vsp == NULL)
    return;

  errbuf = vsp->errbuf;
  cutoff = vsp->cutoff;
  patch_seq = vsp->patch_seq;
  spellfunc = vsp->spellfunc;
  spellcallback = vsp->spellcallback;
  onlyspell = vsp->onlyspell;
  justwarnonspell = vsp->justwarnonspell;
  useSeqMgrIndexes = vsp->useSeqMgrIndexes;
  suppressContext = vsp->suppressContext;
  validateAlignments = vsp->validateAlignments;
  farIDsInAlignments = vsp->farIDsInAlignments;
  alignFindRemoteBsp = vsp->alignFindRemoteBsp;
  doSeqHistAssembly = vsp->doSeqHistAssembly;
  alwaysRequireIsoJTA = vsp->alwaysRequireIsoJTA;
  farFetchCDSproducts = vsp->farFetchCDSproducts;
  farFetchMRNAproducts = vsp->farFetchMRNAproducts;
  locusTagGeneralMatch = vsp->locusTagGeneralMatch;
  validateIDSet = vsp->validateIDSet;
  seqSubmitParent = vsp->seqSubmitParent;
  justShowAccession = vsp->justShowAccession;
  ignoreExceptions = vsp->ignoreExceptions;
  validateExons = vsp->validateExons;
  inferenceAccnCheck = vsp->inferenceAccnCheck;
  testLatLonSubregion = vsp->testLatLonSubregion;
  strictLatLonCountry = vsp->strictLatLonCountry;
  rubiscoTest = vsp->rubiscoTest;
  indexerVersion = vsp->indexerVersion;
  disableSuppression = vsp->disableSuppression;
  validationLimit = vsp->validationLimit;
  errfunc = vsp->errfunc;
  userdata = vsp->userdata;
  convertGiToAccn = vsp->convertGiToAccn;
  sourceQualTags = vsp->sourceQualTags;
  modifiedBases = vsp->modifiedBases;
  sgmlStrings = vsp->sgmlStrings;
  is_htg_in_sep = vsp->is_htg_in_sep;
  is_barcode_sep = vsp->is_barcode_sep;
  is_refseq_in_sep = vsp->is_refseq_in_sep;
  is_gpipe_in_sep = vsp->is_gpipe_in_sep;
  is_gps_in_sep = vsp->is_gps_in_sep;
  is_small_genome_set = vsp->is_small_genome_set;
  other_sets_in_sep = vsp->other_sets_in_sep;
  is_embl_ddbj_in_sep = vsp->is_embl_ddbj_in_sep;
  is_old_gb_in_sep = vsp->is_old_gb_in_sep;
  is_patent_in_sep = vsp->is_patent_in_sep;
  is_insd_in_sep = vsp->is_insd_in_sep;
  only_lcl_gnl_in_sep = vsp->only_lcl_gnl_in_sep;
  has_gnl_prot_sep = vsp->has_gnl_prot_sep;
  bsp_genomic_in_sep = vsp->bsp_genomic_in_sep;
  is_smupd_in_sep = vsp->is_smupd_in_sep;
  feat_loc_has_gi = vsp->feat_loc_has_gi;
  feat_prod_has_gi = vsp->feat_prod_has_gi;
  has_multi_int_genes = vsp->has_multi_int_genes;
  has_seg_bioseqs = vsp->has_seg_bioseqs;
  far_fetch_failure = vsp->far_fetch_failure;
  MemSet ((VoidPtr) vsp, 0, sizeof (ValidStruct));
  vsp->errbuf = errbuf;
  vsp->cutoff = cutoff;
  vsp->patch_seq = patch_seq;
  vsp->spellfunc = spellfunc;
  vsp->spellcallback = spellcallback;
  vsp->onlyspell = onlyspell;
  vsp->justwarnonspell = justwarnonspell;
  vsp->useSeqMgrIndexes = useSeqMgrIndexes;
  vsp->suppressContext = suppressContext;
  vsp->validateAlignments = validateAlignments;
  vsp->farIDsInAlignments = farIDsInAlignments;
  vsp->alignFindRemoteBsp = alignFindRemoteBsp;
  vsp->doSeqHistAssembly = doSeqHistAssembly;
  vsp->alwaysRequireIsoJTA = alwaysRequireIsoJTA;
  vsp->farFetchCDSproducts = farFetchCDSproducts;
  vsp->farFetchMRNAproducts = farFetchMRNAproducts;
  vsp->locusTagGeneralMatch = locusTagGeneralMatch;
  vsp->validateIDSet = validateIDSet;
  vsp->seqSubmitParent = seqSubmitParent;
  vsp->justShowAccession = justShowAccession;
  vsp->ignoreExceptions = ignoreExceptions;
  vsp->validateExons = validateExons;
  vsp->inferenceAccnCheck = inferenceAccnCheck;
  vsp->testLatLonSubregion = testLatLonSubregion;
  vsp->strictLatLonCountry = strictLatLonCountry;
  vsp->rubiscoTest = rubiscoTest;
  vsp->indexerVersion = indexerVersion;
  vsp->disableSuppression = disableSuppression;
  vsp->validationLimit = validationLimit;
  vsp->errfunc = errfunc;
  vsp->userdata = userdata;
  vsp->convertGiToAccn = convertGiToAccn;
  vsp->sourceQualTags = sourceQualTags;
  vsp->modifiedBases = modifiedBases;
  vsp->sgmlStrings = sgmlStrings;
  vsp->is_htg_in_sep = is_htg_in_sep;
  vsp->is_barcode_sep = is_barcode_sep;
  vsp->is_refseq_in_sep = is_refseq_in_sep;
  vsp->is_gpipe_in_sep = is_gpipe_in_sep;
  vsp->is_gps_in_sep = is_gps_in_sep;
  vsp->is_small_genome_set = is_small_genome_set;
  vsp->other_sets_in_sep = other_sets_in_sep;
  vsp->is_embl_ddbj_in_sep = is_embl_ddbj_in_sep;
  vsp->is_old_gb_in_sep = is_old_gb_in_sep;
  vsp->is_patent_in_sep = is_patent_in_sep;
  vsp->is_insd_in_sep = is_insd_in_sep;
  vsp->only_lcl_gnl_in_sep = only_lcl_gnl_in_sep;
  vsp->has_gnl_prot_sep = has_gnl_prot_sep;
  vsp->bsp_genomic_in_sep = bsp_genomic_in_sep;
  vsp->is_smupd_in_sep = is_smupd_in_sep;
  vsp->feat_loc_has_gi = feat_loc_has_gi;
  vsp->feat_prod_has_gi = feat_prod_has_gi;
  vsp->has_multi_int_genes = has_multi_int_genes;
  vsp->has_seg_bioseqs = has_seg_bioseqs;
  vsp->far_fetch_failure = far_fetch_failure;
  return;
}

NLM_EXTERN ValidStructPtr ValidStructNew (void)
{
  ValidStructPtr  vsp;

  vsp = (ValidStructPtr) MemNew (sizeof (ValidStruct));
  return vsp;
}

NLM_EXTERN ValidStructPtr ValidStructFree (ValidStructPtr vsp)
{
  if (vsp == NULL)
    return vsp;

  MemFree (vsp->errbuf);
  TextFsaFree (vsp->sourceQualTags);
  TextFsaFree (vsp->modifiedBases);
  TextFsaFree (vsp->sgmlStrings);
  return (ValidStructPtr) MemFree (vsp);
}

/*****************************************************************************
*
*   ValidErr()
*
*****************************************************************************/

static void ChangeSeqIdToBestID (SeqIdPtr sip)
{
  BioseqPtr       bsp;
  SeqIdPtr        id;
  Pointer         pnt;

  if (sip == NULL)
    return;
  bsp = BioseqFindCore (sip);
  if (bsp == NULL)
    return;
  id = SeqIdDup (SeqIdFindWorst (bsp->id));
  if (id == NULL)
    return;
  /* now remove SeqId contents to reuse SeqId valnode */
  pnt = sip->data.ptrvalue;
  switch (sip->choice) {
  case SEQID_LOCAL:            /* local */
    ObjectIdFree ((ObjectIdPtr) pnt);
    break;
  case SEQID_GIBBSQ:           /* gibbseq */
  case SEQID_GIBBMT:           /* gibbmt */
    break;
  case SEQID_GIIM:             /* giimid */
    GiimFree ((GiimPtr) pnt);
    break;
  case SEQID_GENBANK:          /* genbank */
  case SEQID_EMBL:             /* embl */
  case SEQID_PIR:              /* pir   */
  case SEQID_SWISSPROT:        /* swissprot */
  case SEQID_OTHER:            /* other */
  case SEQID_DDBJ:
  case SEQID_PRF:
  case SEQID_TPG:
  case SEQID_TPE:
  case SEQID_TPD:
  case SEQID_GPIPE:
    TextSeqIdFree ((TextSeqIdPtr) pnt);
    break;
  case SEQID_PATENT:           /* patent seq id */
    PatentSeqIdFree ((PatentSeqIdPtr) pnt);
    break;
  case SEQID_GENERAL:          /* general */
    DbtagFree ((DbtagPtr) pnt);
    break;
  case SEQID_GI:               /* gi */
    break;
  case SEQID_PDB:
    PDBSeqIdFree ((PDBSeqIdPtr) pnt);
    break;
  }
  sip->choice = id->choice;
  sip->data.ptrvalue = id->data.ptrvalue;
  SeqIdStripLocus (sip);
}

static void ChangeSeqLocToBestID (SeqLocPtr slp)
{
  SeqLocPtr       loc;
  PackSeqPntPtr   psp;
  SeqBondPtr      sbp;
  SeqIntPtr       sinp;
  SeqIdPtr        sip;
  SeqPntPtr       spp;

  while (slp != NULL) {
    switch (slp->choice) {
    case SEQLOC_NULL:
      break;
    case SEQLOC_EMPTY:
    case SEQLOC_WHOLE:
      sip = (SeqIdPtr) slp->data.ptrvalue;
      ChangeSeqIdToBestID (sip);
      break;
    case SEQLOC_INT:
      sinp = (SeqIntPtr) slp->data.ptrvalue;
      if (sinp != NULL) {
        sip = sinp->id;
        ChangeSeqIdToBestID (sip);
      }
      break;
    case SEQLOC_PNT:
      spp = (SeqPntPtr) slp->data.ptrvalue;
      if (spp != NULL) {
        sip = spp->id;
        ChangeSeqIdToBestID (sip);
      }
      break;
    case SEQLOC_PACKED_PNT:
      psp = (PackSeqPntPtr) slp->data.ptrvalue;
      if (psp != NULL) {
        sip = psp->id;
        ChangeSeqIdToBestID (sip);
      }
      break;
    case SEQLOC_PACKED_INT:
    case SEQLOC_MIX:
    case SEQLOC_EQUIV:
      loc = (SeqLocPtr) slp->data.ptrvalue;
      while (loc != NULL) {
        ChangeSeqLocToBestID (loc);
        loc = loc->next;
      }
      break;
    case SEQLOC_BOND:
      sbp = (SeqBondPtr) slp->data.ptrvalue;
      if (sbp != NULL) {
        spp = (SeqPntPtr) sbp->a;
        if (spp != NULL) {
          sip = spp->id;
          ChangeSeqIdToBestID (sip);
        }
        spp = (SeqPntPtr) sbp->b;
        if (spp != NULL) {
          sip = spp->id;
          ChangeSeqIdToBestID (sip);
        }
      }
      break;
    case SEQLOC_FEAT:
      break;
    default:
      break;
    }
    slp = slp->next;
  }
}

static Int2 WorstBioseqLabel (BioseqPtr bsp, CharPtr buffer, Int2 buflen, Uint1 content)
{
  CharPtr         tmp;
  Char            label[60];
  Int2            diff, len;
  SeqIdPtr        sip;
  AsnModulePtr    amp;
  AsnTypePtr      ratp, matp;

  if ((bsp == NULL) || (buflen < 1))
    return 0;

  len = buflen;
  label[0] = '\0';

  if (content != OM_LABEL_TYPE) {
    sip = SeqIdStripLocus (SeqIdDup (SeqIdFindWorst (bsp->id)));
    SeqIdWrite (sip, label, PRINTID_FASTA_SHORT, 39);
    SeqIdFree (sip);
    if (content == OM_LABEL_CONTENT)
      return LabelCopy (buffer, label, buflen);

    diff = LabelCopyExtra (buffer, label, buflen, NULL, ": ");
    buflen -= diff;
    buffer += diff;
  }

  amp = AsnAllModPtr ();
  ratp = AsnTypeFind (amp, "Seq-inst.repr");
  matp = AsnTypeFind (amp, "Seq-inst.mol");

  label[0] = '\0';
  tmp = label;
  tmp = StringMove (tmp, AsnEnumTypeStr (ratp, (Int2) (bsp->repr)));
  tmp = StringMove (tmp, ", ");
  tmp = StringMove (tmp, AsnEnumTypeStr (matp, (Int2) (bsp->mol)));
  sprintf (tmp, " len= %ld", (long) (bsp->length));
  diff = LabelCopy (buffer, label, buflen);
  buflen -= diff;
  buffer += diff;

  if (content != OM_LABEL_SUMMARY)
    return (len - buflen);

  return (len - buflen);        /* SUMMARY not done yet */
}

static CharPtr categoryLabel [] = {
  NULL, "SEQ_INST", "SEQ_DESCR", "GENERIC", "SEQ_PKG", "SEQ_FEAT", "SEQ_ALIGN", "SEQ_GRAPH", "SEQ_ANNOT"
};

NLM_EXTERN CharPtr GetValidCategoryName (int errcode)

{
  if (errcode >= 1 && errcode < sizeof (categoryLabel)) return categoryLabel [errcode];
  return NULL;
}

static CharPtr err1Label [] = {
  NULL,
  "ExtNotAllowed",
  "ExtBadOrMissing",
  "SeqDataNotFound",
  "SeqDataNotAllowed",
  "ReprInvalid",
  "CircularProtein",
  "DSProtein",
  "MolNotSet",
  "MolOther",
  "FuzzyLen",
  "InvalidLen",
  "InvalidAlphabet",
  "SeqDataLenWrong",
  "SeqPortFail",
  "InvalidResidue",
  "StopInProtein",
  "PartialInconsistent",
  "ShortSeq",
  "NoIdOnBioseq",
  "BadDeltaSeq",
  "LongHtgsSequence",
  "LongLiteralSequence",
  "SequenceExceeds350kbp",
  "ConflictingIdsOnBioseq",
  "MolNuclAcid",
  "ConflictingBiomolTech",
  "SeqIdNameHasSpace",
  "IdOnMultipleBioseqs",
  "DuplicateSegmentReferences",
  "TrailingX",
  "BadSeqIdFormat",
  "PartsOutOfOrder",
  "BadSecondaryAccn",
  "ZeroGiNumber",
  "RnaDnaConflict",
  "HistoryGiCollision",
  "GiWithoutAccession",
  "MultipleAccessions",
  "HistAssemblyMissing",
  "TerminalNs",
  "UnexpectedIdentifierChange",
  "InternalNsInSeqLit",
  "SeqLitGapLength0",
  "TpaAssmeblyProblem",
  "SeqLocLength",
  "MissingGaps",
  "CompleteTitleProblem",
  "CompleteCircleProblem",
  "BadHTGSeq",
  "GapInProtein",
  "BadProteinStart",
  "TerminalGap",
  "OverlappingDeltaRange",
  "LeadingX",
  "InternalNsInSeqRaw",
  "InternalNsAdjacentToGap",
  "CaseDifferenceInSeqID",
  "DeltaComponentIsGi0",
  "FarFetchFailure",
  "InternalGapsInSeqRaw",
  "SelfReferentialSequence",
  "WholeComponent",
  "TSAHistAssemblyMissing",
  "ProteinsHaveGeneralID",
  "HighNContent",
  "SeqLitDataLength0",
  "DSmRNA",
  "HighNContentStretch",
  "HighNContentPercent",
  "BadSegmentedSeq",
  "SeqLitGapFuzzNot100"
};

static CharPtr err2Label [] = {
  NULL,
  "BioSourceMissing",
  "InvalidForType",
  "FileOpenCollision",
  "Unknown",
  "NoPubFound",
  "NoOrgFound",
  "MultipleBioSources",
  "NoMolInfoFound",
  "BadCountryCode",
  "NoTaxonID",
  "InconsistentBioSources",
  "MissingLineage",
  "SerialInComment",
  "BioSourceNeedsFocus",
  "BadOrganelle",
  "MultipleChromosomes",
  "BadSubSource",
  "BadOrgMod",
  "InconsistentProteinTitle",
  "Inconsistent",
  "ObsoleteSourceLocation",
  "ObsoleteSourceQual",
  "StructuredSourceNote",
  "UnnecessaryBioSourceFocus",
  "RefGeneTrackingWithoutStatus",
  "UnwantedCompleteFlag",
  "CollidingPublications",
  "TransgenicProblem",
  "TaxonomyLookupProblem",
  "MultipleTitles",
  "RefGeneTrackingOnNonRefSeq",
  "BioSourceInconsistency",
  "FastaBracketTitle",
  "MissingText",
  "BadCollectionDate",
  "BadPCRPrimerSequence",
  "BadPunctuation",
  "BadPCRPrimerName",
  "BioSourceOnProtein",
  "BioSourceDbTagConflict",
  "DuplicatePCRPrimerSequence",
  "MultipleNames",
  "MultipleComments",
  "LatLonProblem",
  "LatLonFormat",
  "LatLonRange",
  "LatLonValue",
  "LatLonCountry",
  "LatLonState",
  "BadSpecificHost",
  "RefGeneTrackingIllegalStatus",
  "ReplacedCountryCode",
  "BadInstitutionCode",
  "BadCollectionCode",
  "BadVoucherID",
  "UnstructuredVoucher",
  "ChromosomeLocation",
  "MultipleSourceQualifiers",
  "UnbalancedParentheses",
  "MultipleSourceVouchers",
  "BadCountryCapitalization",
  "WrongVoucherType",
  "UserObjectProblem",
  "TitleHasPMID",
  "BadKeyword",
  "NoOrganismInTitle",
  "MissingChromosome",
  "LatLonAdjacent",
  "BadStrucCommInvalidFieldName",
  "BadStrucCommInvalidFieldValue",
  "BadStrucCommMissingField",
  "BadStrucCommFieldOutOfOrder",
  "BadStrucCommMultipleFields",
  "BioSourceNeedsChromosome",
  "MolInfoConflictsWithBioSource",
  "MissingKeyword",
  "FakeStructuredComment",
  "StructuredCommentPrefixOrSuffixMissing",
  "LatLonWater",
  "LatLonOffshore",
  "MissingPersonalCollectionName",
  "LatLonPrecision",
  "DBLinkProblem",
  "FinishedStatusForWGS"
};

static CharPtr err3Label [] = {
  NULL,
  "NonAsciiAsn",
  "Spell",
  "AuthorListHasEtAl",
  "MissingPubInfo",
  "UnnecessaryPubEquiv",
  "BadPageNumbering",
  "MedlineEntryPub",
  "BadDate",
  "StructuredCitGenCit",
  "CollidingSerialNumbers",
  "EmbeddedScript",
  "PublicationInconsistency",
  "SgmlPresentInText",
  "UnexpectedPubStatusComment",
  "PastReleaseDate"
};

static CharPtr err4Label [] = {
  NULL,
  "NoCdRegionPtr",
  "NucProtProblem",
  "SegSetProblem",
  "EmptySet",
  "NucProtNotSegSet",
  "SegSetNotParts",
  "SegSetMixedBioseqs",
  "PartsSetMixedBioseqs",
  "PartsSetHasSets",
  "FeaturePackagingProblem",
  "GenomicProductPackagingProblem",
  "InconsistentMolInfoBiomols",
  "ArchaicFeatureLocation",
  "ArchaicFeatureProduct",
  "GraphPackagingProblem",
  "InternalGenBankSet",
  "ConSetProblem",
  "NoBioseqFound",
  "INSDRefSeqPackaging",
  "GPSnonGPSPackaging",
  "RefSeqPopSet",
  "BioseqSetClassNotSet",
  "OrphanedProtein",
  "MissingSetTitle",
  "NucProtSetHasTitle",
  "ComponentMissingTitle",
  "SingleItemSet",
  "MisplacedMolInfo",
  "ImproperlyNestedSets"
};

static CharPtr err5Label [] = {
  NULL,
  "InvalidForType",
  "PartialProblem",
  "InvalidType",
  "Range",
  "MixedStrand",
  "SeqLocOrder",
  "CdTransFail",
  "StartCodon",
  "InternalStop",
  "NoProtein",
  "MisMatchAA",
  "TransLen",
  "NoStop",
  "TranslExcept",
  "NoProtRefFound",
  "NotSpliceConsensus",
  "OrfCdsHasProduct",
  "GeneRefHasNoData",
  "ExceptInconsistent",
  "ProtRefHasNoData",
  "GenCodeMismatch",
  "RNAtype0",
  "UnknownImpFeatKey",
  "UnknownImpFeatQual",
  "WrongQualOnImpFeat",
  "MissingQualOnImpFeat",
  "PseudoCdsHasProduct",
  "IllegalDbXref",
  "FarLocation",
  "DuplicateFeat",
  "UnnecessaryGeneXref",
  "TranslExceptPhase",
  "TrnaCodonWrong",
  "BothStrands",
  "CDSgeneRange",
  "CDSmRNArange",
  "OverlappingPeptideFeat",
  "SerialInComment",
  "MultipleCDSproducts",
  "FocusOnBioSourceFeature",
  "PeptideFeatOutOfFrame",
  "InvalidQualifierValue",
  "MultipleMRNAproducts",
  "mRNAgeneRange",
  "TranscriptLen",
  "TranscriptMismatches",
  "CDSproductPackagingProblem",
  "DuplicateInterval",
  "PolyAsiteNotPoint",
  "ImpFeatBadLoc",
  "LocOnSegmentedBioseq",
  "UnnecessaryCitPubEquiv",
  "ImpCDShasTranslation",
  "ImpCDSnotPseudo",
  "MissingMRNAproduct",
  "AbuttingIntervals",
  "CollidingGeneNames",
  "MultiIntervalGene",
  "FeatContentDup",
  "BadProductSeqId",
  "RnaProductMismatch",
  "MissingCDSproduct",
  "BadTrnaCodon",
  "BadTrnaAA",
  "OnlyGeneXrefs",
  "UTRdoesNotAbutCDS",
  "BadConflictFlag",
  "ConflictFlagSet",
  "LocusTagProblem",
  "CollidingLocusTags",
  "AltStartCodon",
  "PartialsInconsistent",
  "GenesInconsistent",
  "DuplicateTranslExcept",
  "TranslExceptAndRnaEditing",
  "NoNameForProtein",
  "TaxonDbxrefOnFeature",
  "UnindexedFeature",
  "CDSmRNAmismatch",
  "UnnecessaryException",
  "LocusTagProductMismatch",
  "MrnaTransFail",
  "PseudoCdsViaGeneHasProduct",
  "MissingGeneXref",
  "FeatureCitationProblem",
  "NestedSeqLocMix",
  "WrongQualOnFeature",
  "MissingQualOnFeature",
  "CodonQualifierUsed",
  "UnknownFeatureQual",
  "BadCharInAuthorName",
  "PolyATail",
  "ProteinNameEndsInBracket",
  "CDSwithMultipleMRNAs",
  "MultipleEquivBioSources",
  "MultipleEquivPublications",
  "BadFullLengthFeature",
  "RedundantFields",
  "CDSwithNoMRNAOverlap",
  "FeatureProductInconsistency",
  "ImproperBondLocation",
  "GeneXrefWithoutGene",
  "SeqFeatXrefProblem",
  "ProductFetchFailure",
  "SuspiciousGeneXref",
  "MissingTrnaAA",
  "CollidingFeatureIDs",
  "ExceptionProblem",
  "PolyAsignalNotRange",
  "OldLocusTagMismtach",
  "DuplicateGeneOntologyTerm",
  "InvalidInferenceValue",
  "HpotheticalProteinMismatch",
  "FeatureRefersToAccession",
  "SelfReferentialProduct",
  "ITSdoesNotAbutRRNA",
  "FeatureSeqIDCaseDifference",
  "FeatureLocationIsGi0",
  "GapFeatureProblem",
  "PseudoCdsHasProtXref",
  "ErroneousException",
  "SegmentedGeneProblem",
  "WholeLocation",
  "BadEcNumberFormat",
  "BadEcNumberValue",
  "EcNumberProblem",
  "VectorContamination",
  "MinusStrandProtein",
  "BadProteinName",
  "GeneXrefWithoutLocus",
  "UTRdoesNotExtendToEnd",
  "CDShasTooManyXs",
  "SuspiciousFrame",
  "TerminalXDiscrepancy",
  "UnnecessaryTranslExcept",
  "SuspiciousQualifierValue",
  "NotSpliceConsensusDonor",
  "NotSpliceConsensusAcceptor",
  "RareSpliceConsensusDonor",
  "SeqFeatXrefNotReciprocal",
  "SeqFeatXrefFeatureMissing",
  "FeatureInsideGap",
  "FeatureCrossesGap",
  "BadAuthorSuffix",
  "BadAnticodonAA",
  "BadAnticodonCodon",
  "BadAnticodonStrand",
  "UndesiredGeneSynonym",
  "UndesiredProteinName",
  "FeatureBeginsOrEndsInGap",
  "GeneOntologyTermMissingGOID",
  "PseudoRnaHasProduct",
  "PseudoRnaViaGeneHasProduct",
  "BadRRNAcomponentOrder",
  "BadRRNAcomponentOverlap",
  "MissingGeneLocusTag",
  "MultipleProtRefs",
  "BadInternalCharacter",
  "BadTrailingCharacter",
  "BadTrailingHyphen",
  "MultipleGeneOverlap",
  "BadCharInAuthorLastName",
  "PseudoCDSmRNArange",
  "ExtendablePartialProblem",
  "GeneXrefNeeded",
  "RubiscoProblem",
  "UnqualifiedException",
  "ProteinNameHasPMID",
  "BadGeneOntologyFormat",
  "InconsistentGeneOntologyTermAndId",
  "MultiplyAnnotatedGenes",
  "ReplicatedGeneSequence",
  "ShortIntron",
  "GeneXrefStrandProblem",
  "CDSmRNAXrefLocationProblem",
  "LocusCollidesWithLocusTag",
  "IdenticalGeneSymbolAndSynonym",
  "NeedsNote",
  "RptUnitRangeProblem",
  "TooManyInferenceAccessions",
  "IntervalBeginsOrEndsInGap"
};

static CharPtr err6Label [] = {
  NULL,
  "SeqIdProblem",
  "StrandRev",
  "DensegLenStart",
  "StartLessthanZero",
  "StartMorethanBiolen",
  "EndLessthanZero",
  "EndMorethanBiolen",
  "LenLessthanZero",
  "LenMorethanBiolen",
  "SumLenStart",
  "AlignDimSeqIdNotMatch",
  "SegsDimSeqIdNotMatch",
  "FastaLike",
  "NullSegs",
  "SegmentGap",
  "SegsDimOne",
  "AlignDimOne",
  "Segtype",
  "BlastAligns",
  "PercentIdentity",
  "ShortAln",
  "UnexpectedAlignmentType"
};

static CharPtr err7Label [] = {
  NULL,
  "GraphMin",
  "GraphMax",
  "GraphBelow",
  "GraphAbove",
  "GraphByteLen",
  "GraphOutOfOrder",
  "GraphBioseqLen",
  "GraphSeqLitLen",
  "GraphSeqLocLen",
  "GraphStartPhase",
  "GraphStopPhase",
  "GraphDiffNumber",
  "GraphACGTScore",
  "GraphNScore",
  "GraphGapScore",
  "GraphOverlap",
  "GraphBioseqId",
  "GraphACGTScoreMany",
  "GraphNScoreMany",
  "GraphLocInvalid"
};

static CharPtr err8Label [] = {
  NULL,
  "AnnotIDs",
  "AnnotLOCs"
};

NLM_EXTERN CharPtr GetValidErrorName (int errcode, int subcode)

{
  if (errcode < 1 || errcode >= sizeof (categoryLabel)) return NULL;
  switch (errcode) {
    case 1 :
      if (subcode >= 1 && subcode < sizeof (err1Label)) return err1Label [subcode];
      break;
    case 2 :
      if (subcode >= 1 && subcode < sizeof (err2Label)) return err2Label [subcode];
      break;
    case 3 :
      if (subcode >= 1 && subcode < sizeof (err3Label)) return err3Label [subcode];
      break;
    case 4 :
      if (subcode >= 1 && subcode < sizeof (err4Label)) return err4Label [subcode];
      break;
    case 5 :
      if (subcode >= 1 && subcode < sizeof (err5Label)) return err5Label [subcode];
      break;
    case 6 :
      if (subcode >= 1 && subcode < sizeof (err6Label)) return err6Label [subcode];
      break;
    case 7 :
      if (subcode >= 1 && subcode < sizeof (err7Label)) return err7Label [subcode];
      break;
    case 8 :
      if (subcode >= 1 && subcode < sizeof (err8Label)) return err8Label [subcode];
      break;
    default :
      break;
  }
  return NULL;
}

NLM_EXTERN CharPtr GetValidExplanation (int errcode, int subcode)

{
  return Nlm_GetErrLongText (THIS_MODULE, errcode, subcode);
}

static void CustValErr (ValidStructPtr vsp, ErrSev severity, int errcode, int subcode)

{
  CharPtr           accession = NULL, context = NULL, label = NULL, location = NULL,
                    message = NULL, objtype = NULL, product = NULL, featureID = NULL;
  BioseqPtr         bsp;
  BioseqSetPtr      bssp;
  Int2              buflen, diff, wrklen;
  CharPtr           ctmp, tmp;
  DbtagPtr          dbt;
  Uint2             entityID = 0, itemtype = 0;
  ValidErrorFunc    errfunc;
  GatherContextPtr  gcp;
  Char              id [64], numbuf [15];
  Uint4             itemID = 0;
  ObjectIdPtr       oip;
  ObjValNodePtr     ovp;
  SeqDescrPtr       sdp;
  SeqEntryPtr       sep;
  SeqFeatPtr        sfp = NULL;
  SeqIdPtr          sip;
  SeqLocPtr         slp;

  if (vsp == NULL) return;
  errfunc = vsp->errfunc;
  if (errfunc == NULL) return;

  gcp = vsp->gcp;
  if (gcp != NULL) {
    entityID = gcp->entityID;
    itemtype = gcp->thistype;
    itemID = gcp->itemID;
  }

  if (severity < SEV_NONE || severity > SEV_MAX) {
    severity = SEV_MAX;
  }

  sip = NULL;
  if (vsp->sfp != NULL) {
    sfp = vsp->sfp;
    bsp = BioseqFindFromSeqLoc (sfp->location);
    if (bsp != NULL) {
      sip = SeqIdFindWorst (bsp->id);
    }
  } else if (vsp->descr != NULL) {
    sdp = vsp->descr;
    if (sdp != NULL && sdp->extended != 0) {
      ovp = (ObjValNodePtr) sdp;
      if (ovp->idx.parenttype == OBJ_BIOSEQ) {
        bsp = (BioseqPtr) ovp->idx.parentptr;
        if (bsp != NULL) {
          sip = SeqIdFindWorst (bsp->id);
        }
      } else if (ovp->idx.parenttype == OBJ_BIOSEQSET) {
        bssp = (BioseqSetPtr) ovp->idx.parentptr;
        if (bssp != NULL) {
          sep = bssp->seqentry;
          if (sep != NULL) {
            sep = FindNthBioseq (sep, 1);
            if (sep != NULL) {
              bsp = (BioseqPtr) sep->data.ptrvalue;
              if (bsp != NULL) {
                sip = SeqIdFindWorst (bsp->id);
              }
            }
          }
        }
      }
    }
  } else if (vsp->bsp != NULL) {
    bsp = vsp->bsp;
    sip = SeqIdFindWorst (bsp->id);
  } else if (vsp->bssp != NULL) {
    bssp = vsp->bssp;
    sep = bssp->seqentry;
    if (sep != NULL) {
      sep = FindNthBioseq (sep, 1);
      if (sep != NULL) {
        bsp = (BioseqPtr) sep->data.ptrvalue;
        if (bsp != NULL) {
          sip = SeqIdFindWorst (bsp->id);
        }
      }
    }
  }
  if (sip != NULL) {
    SeqIdWrite (sip, id, PRINTID_REPORT, sizeof (id) - 1);
    accession = id;
  }

  if (vsp->sfp != NULL) {
    objtype = "FEATURE";
  } else if (vsp->descr != NULL) {
    objtype = "DESCRIPTOR";
  } else if (vsp->bsp != NULL) {
    objtype = "BIOSEQ";
  } else if (vsp->bssp != NULL) {
    objtype = "BIOSEQ-SET";
  }

  message = vsp->errbuf;

  tmp = vsp->errbuf;
  buflen = 4000;
  while (*tmp != '\0') {
    buflen--;
    tmp++;
  }
  tmp++;
  *tmp = '\0';

  wrklen = buflen;
  if (wrklen > 2000) {
     wrklen -= 1000;
  }

  if (vsp->sfp != NULL) {
    sfp = vsp->sfp;
    label = tmp;
    diff = FeatDefLabel (sfp, tmp, wrklen, OM_LABEL_BOTH);
    buflen -= diff;
    tmp += diff;
    *tmp = '\0';
    tmp++;
    *tmp = '\0';

    featureID = tmp;
    dbt = NULL;
    oip = NULL;
    if (sfp->id.choice == 3) {
      oip = (ObjectIdPtr) sfp->id.value.ptrvalue;
    } else if (sfp->id.choice == 4) {
      dbt = (DbtagPtr) sfp->id.value.ptrvalue;
      if (dbt != NULL) {
        oip = dbt->tag;
      }
    }
    if (oip != NULL) {
      if (dbt != NULL && dbt->db != NULL) {
        diff = LabelCopyExtra (tmp, dbt->db, buflen, NULL, ":");
        buflen -= diff;
        tmp += diff;
        *tmp = '\0';
      }
      if (oip->str != NULL) {
        diff = LabelCopyExtra (tmp, oip->str, buflen, NULL, NULL);
      } else {
        sprintf (numbuf, "%ld", (long) oip->id);
        diff = LabelCopyExtra (tmp, numbuf, buflen, NULL, NULL);
      }
      buflen -= diff;
      tmp += diff;
      *tmp = '\0';
    }

    tmp++;
    *tmp = '\0';
  } else if (vsp->descr != NULL) {
    label = tmp;
    diff = SeqDescLabel (vsp->descr, tmp, wrklen, OM_LABEL_BOTH);

    if (diff > 100 && vsp->descr->choice == Seq_descr_comment && errcode == 2 && subcode == 77) {
      diff = 100;
      *(tmp + diff - 3) = '.';
      *(tmp + diff - 2) = '.';
      *(tmp + diff - 1) = '.';
    }
    buflen -= diff;
    tmp += diff;
    *tmp = '\0';
    tmp++;
    *tmp = '\0';
  } else if (vsp->bsp != NULL) {
    label = tmp;
    if (vsp->convertGiToAccn) {
      diff = WorstBioseqLabel (vsp->bsp, tmp, wrklen, OM_LABEL_CONTENT);
    } else {
      diff = BioseqLabel (vsp->bsp, tmp, wrklen, OM_LABEL_BOTH);
    }
    buflen -= diff;
    tmp += diff;
    *tmp = '\0';
    tmp++;
    *tmp = '\0';
  } else if (vsp->bssp != NULL) {
    label = tmp;
    diff = BioseqSetLabel (vsp->bssp, tmp, wrklen, OM_LABEL_BOTH);
    buflen -= diff;
    tmp += diff;
    *tmp = '\0';
    tmp++;
    *tmp = '\0';
  }

  if (vsp->sfp != NULL) {
    sfp = vsp->sfp;

    if (sfp->location != NULL) {
      ctmp = NULL;
      slp = NULL;
      /*
      if (vsp->suppressContext) {
        slp = AsnIoMemCopy (sfp->location, (AsnReadFunc) SeqLocAsnRead, (AsnWriteFunc) SeqLocAsnWrite);
        ChangeSeqLocToBestID (slp);
        ctmp = SeqLocPrint (slp);
        SeqLocFree (slp);
      } else {
        ctmp = SeqLocPrint (sfp->location);
      }
      */
      slp = AsnIoMemCopy (sfp->location, (AsnReadFunc) SeqLocAsnRead, (AsnWriteFunc) SeqLocAsnWrite);
      ChangeSeqLocToBestID (slp);
      ctmp = SeqLocPrint (slp);
      SeqLocFree (slp);
      if (ctmp != NULL) {
        if (StringLen (ctmp) > 800) {
          StringCpy (ctmp + 797, "...");
        }
        location = tmp;
        diff = LabelCopyExtra (tmp, ctmp, buflen, "[", "]");
        buflen -= diff;
        tmp += diff;
        MemFree (ctmp);
        *tmp = '\0';
        tmp++;
        *tmp = '\0';

        sip = SeqLocId (sfp->location);
        if (sip != NULL) {
          bsp = BioseqFind (sip);
          if (bsp != NULL) {
            context = tmp;
            diff = LabelCopy (tmp, "[", buflen);
            buflen -= diff;
            tmp += diff;

            diff = BioseqLabel (bsp, tmp, buflen, OM_LABEL_BOTH);
            buflen -= diff;
            tmp += diff;

            diff = LabelCopy (tmp, "]", buflen);
            buflen -= diff;
            tmp += diff;
          }
        }
        *tmp = '\0';
        tmp++;
        *tmp = '\0';
      }
    }

    if (sfp->product != NULL) {
      ctmp = NULL;
      slp = NULL;
      /*
      if (vsp->suppressContext) {
        slp = AsnIoMemCopy (sfp->product, (AsnReadFunc) SeqLocAsnRead, (AsnWriteFunc) SeqLocAsnWrite);
        ChangeSeqLocToBestID (slp);
        ctmp = SeqLocPrint (slp);
        SeqLocFree (slp);
      } else {
        ctmp = SeqLocPrint (sfp->product);
      }
      */
      slp = AsnIoMemCopy (sfp->product, (AsnReadFunc) SeqLocAsnRead, (AsnWriteFunc) SeqLocAsnWrite);
      ChangeSeqLocToBestID (slp);
      ctmp = SeqLocPrint (slp);
      SeqLocFree (slp);
      if (ctmp != NULL) {
        if (StringLen (ctmp) > 800) {
          StringCpy (ctmp + 797, "...");
        }
        product = tmp;
        diff = LabelCopyExtra (tmp, ctmp, buflen, "[", "]");
        buflen -= diff;
        tmp += diff;
        *tmp = '\0';
        tmp++;
        *tmp = '\0';
        MemFree (ctmp);
      }
    }
  } else if (vsp->descr != NULL) {
    if (vsp->bsp != NULL) {
      context = tmp;
      diff = LabelCopy (tmp, "BIOSEQ: ", buflen);
      buflen -= diff;
      tmp += diff;
      if (vsp->suppressContext || vsp->convertGiToAccn) {
        diff = WorstBioseqLabel (vsp->bsp, tmp, buflen, OM_LABEL_CONTENT);
      } else {
        diff = BioseqLabel (vsp->bsp, tmp, buflen, OM_LABEL_BOTH);
      }
      buflen -= diff;
      tmp += diff;
      *tmp = '\0';
      tmp++;
      *tmp = '\0';
    } else if (vsp->bssp != NULL) {
      context = tmp;
      diff = LabelCopy (tmp, "BIOSEQ-SET: ", buflen);
      buflen -= diff;
      tmp += diff;

      if (vsp->suppressContext || vsp->convertGiToAccn) {
        diff = BioseqSetLabel (vsp->bssp, tmp, buflen, OM_LABEL_CONTENT);
      } else {
        diff = BioseqSetLabel (vsp->bssp, tmp, buflen, OM_LABEL_BOTH);
      }
      buflen -= diff;
      tmp += diff;
      *tmp = '\0';
      tmp++;
      *tmp = '\0';
    }
  }

  (*errfunc) (severity, errcode, subcode, entityID, itemtype, itemID, accession,
              featureID, message, objtype, label, context, location, product, vsp->userdata);
}


/* framework for suppressing validator errors using a list-based strategy */
typedef Boolean (*ValidErrSuppressFunc) PROTO ((ValidStructPtr));

static Boolean IsGenomicPipeline (ValidStructPtr vsp)
{
  if (vsp == NULL) {
    return FALSE;
  } else if (vsp->bsp_genomic_in_sep && vsp->is_gpipe_in_sep) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static Boolean IsUnclassifiedExcept (ValidStructPtr vsp)
{
  Boolean rval = FALSE;
  if (vsp == NULL || vsp->sfp == NULL) {
    return FALSE;
  }
  if (vsp->sfp->excpt && (! vsp->ignoreExceptions)) {
    if (vsp->sfp->data.choice == SEQFEAT_CDREGION) {
      if (StringStr (vsp->sfp->except_text, "unclassified translation discrepancy") != NULL) {
        rval = TRUE;
      }
    } else if (vsp->sfp->idx.subtype == FEATDEF_mRNA) {
      if (StringStr (vsp->sfp->except_text, "unclassified transcription discrepancy") != NULL) {
        rval = TRUE;
      }
    }
  }
  return rval;
}


static Boolean IsNotUnclassifiedExcept (ValidStructPtr vsp)
{
  return !IsUnclassifiedExcept(vsp);
}


static Boolean IsUnclassifedExceptAndGenomicPipeline (ValidStructPtr vsp)
{
  if (IsGenomicPipeline(vsp) && IsUnclassifiedExcept(vsp)) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static Boolean NonconsensusExcept (ValidStructPtr vsp)
{
  Boolean rval = FALSE;
  if (vsp == NULL || vsp->sfp == NULL) {
    return FALSE;
  }
  if (vsp->sfp->excpt && (! vsp->ignoreExceptions)) {
    if (StringISearch (vsp->sfp->except_text, "nonconsensus splice site") != NULL ||
        StringISearch (vsp->sfp->except_text, "heterogeneous population sequenced") != NULL ||
        StringISearch (vsp->sfp->except_text, "low-quality sequence region") != NULL ||
        StringISearch (vsp->sfp->except_text, "artificial location") != NULL) {
      rval = TRUE;
    }
  }
  return rval;
}


typedef struct validerrsuppression {
  int code1;
  int code2;
  CharPtr search_phrase;
  CharPtr exclude_phrase;
  ValidErrSuppressFunc func;
} ValidErrSuppressionData, PNTR ValidErrSuppressionPtr;

static ValidErrSuppressionData valid_suppress[] = {
  {ERR_SEQ_FEAT_PartialProblem, "When SeqFeat.product is a partial Bioseq, SeqFeat.location should also be partial", NULL, IsGenomicPipeline },
  {ERR_SEQ_FEAT_PartialProblem, "End of location should probably be partial", NULL, IsGenomicPipeline},
  {ERR_SEQ_FEAT_PartialProblem, "This SeqFeat should not be partial", NULL, IsGenomicPipeline},
  {ERR_SEQ_FEAT_PartialProblem, "AND is not at consensus splice site", NULL, IsGenomicPipeline},
  {ERR_SEQ_FEAT_PartialProblem, "PartialLocation: Internal partial intervals do not include first/last residue of sequence", NULL, IsGenomicPipeline},
  {ERR_SEQ_FEAT_PartialProblem, "AND is not at consensus splice site", NULL, NonconsensusExcept},
  {ERR_SEQ_FEAT_PartialProblem, "(but is at consensus splice site)", NULL, IsGenomicPipeline},
  {ERR_SEQ_FEAT_PartialProblem, "PartialLocation: Start does not include first/last residue of sequence", NULL, IsGenomicPipeline},
  {ERR_SEQ_FEAT_PartialProblem, "PartialLocation: Stop does not include first/last residue of sequence", NULL, IsGenomicPipeline},
  {ERR_SEQ_FEAT_PartialsInconsistent, NULL, NULL, IsGenomicPipeline },
  {ERR_SEQ_FEAT_PolyATail, NULL, NULL, IsGenomicPipeline },
  {ERR_SEQ_FEAT_InternalStop, NULL, NULL, IsUnclassifedExceptAndGenomicPipeline},
  {ERR_SEQ_FEAT_StartCodon , NULL, NULL, IsUnclassifiedExcept}

};

const Int4 kNumSuppressionRules = sizeof (valid_suppress) / sizeof (ValidErrSuppressionData);

static Boolean ShouldSuppressValidErr (ValidStructPtr vsp, int code1, int code2, const char *fmt)
{
  Int4 i;
  Boolean rval = FALSE;

  if (vsp->disableSuppression) return FALSE;

  for (i = 0; i < kNumSuppressionRules && !rval; i++) {
    if (code1 == valid_suppress[i].code1 && code2 == valid_suppress[i].code2
        && (valid_suppress[i].search_phrase == NULL || StringISearch (fmt, valid_suppress[i].search_phrase) != NULL)
        && (valid_suppress[i].func == NULL || valid_suppress[i].func(vsp))
        && (valid_suppress[i].exclude_phrase == NULL || StringISearch (fmt, valid_suppress[i].exclude_phrase) == NULL)) {
      rval = TRUE;
    }
  }

  return rval;
}


/* framework for changing validator warnings using a list-based strategy */
typedef int (*ValidErrSevChangeFunc) PROTO ((int, ValidStructPtr));

typedef struct validerrsevchange {
  int code1;
  int code2;
  CharPtr search_phrase;
  CharPtr exclude_phrase;
  ValidErrSevChangeFunc func;
} ValidErrSevChangeData, PNTR ValidErrSevChangePtr;


static int LowerToInfoForGenomic (int severity, ValidStructPtr vsp)
{
  if (IsGenomicPipeline(vsp)) {
    return SEV_INFO;
  } else {
    return severity;
  }
}


static int WarnForGPSOrRefSeq (int severity, ValidStructPtr vsp)
{
  Boolean         gpsOrRefSeq = FALSE;
  SeqEntryPtr     sep;
  SeqFeatPtr      sfp;
  BioseqSetPtr    bssp;
  SeqLocPtr       head, slp = NULL, nxt;
  SeqIdPtr        sip, id;
  BioseqPtr       bsp;
  TextSeqIdPtr    tsip;

  sep = vsp->sep;
  if (sep != NULL && IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && bssp->_class == BioseqseqSet_class_gen_prod_set) {
      gpsOrRefSeq = TRUE;
    }
  }

  if (!gpsOrRefSeq) {
    sfp = vsp->sfp;
    head = sfp->location;
    slp = SeqLocFindPart (head, slp, EQUIV_IS_ONE);
    while (slp != NULL && !gpsOrRefSeq) {
      sip = SeqLocId (slp);
      if (sip == NULL)
        break;
      nxt = SeqLocFindPart (head, slp, EQUIV_IS_ONE);

      /* genomic product set or NT_ contig always relaxes to SEV_WARNING */
      bsp = BioseqFind (sip);
      if (bsp != NULL) {
        for (id = bsp->id; id != NULL; id = id->next) {
          if (id->choice == SEQID_OTHER) {
            tsip = (TextSeqIdPtr) id->data.ptrvalue;
            if (tsip != NULL && tsip->accession != NULL) {
              gpsOrRefSeq = TRUE;
            }
          }
        }
      }

      slp = nxt;
    }
  }
  if (gpsOrRefSeq) {
    if (severity > SEV_WARNING) {
      severity = SEV_WARNING;
    }
  }
  return severity;
}


static ValidErrSevChangeData valid_sevchange[] = {
  {ERR_SEQ_FEAT_NotSpliceConsensusDonor, "Splice donor consensus (GT) not found at start of intron, position", NULL, LowerToInfoForGenomic},
  {ERR_SEQ_FEAT_NotSpliceConsensusAcceptor, "Splice acceptor consensus (AG) not found at end of intron, position", NULL, LowerToInfoForGenomic},
  {ERR_SEQ_FEAT_NotSpliceConsensusDonor, "Splice donor consensus (GT) not found after exon", NULL, LowerToInfoForGenomic},
  {ERR_SEQ_FEAT_NotSpliceConsensusDonor, "Splice donor consensus (GT) not found after exon", NULL, WarnForGPSOrRefSeq},
  {ERR_SEQ_FEAT_NotSpliceConsensusAcceptor, "Splice acceptor consensus (AG) not found before exon", NULL, LowerToInfoForGenomic},
  {ERR_SEQ_FEAT_NotSpliceConsensusAcceptor, "Splice acceptor consensus (AG) not found before exon", NULL, WarnForGPSOrRefSeq},
};

const Int4 kNumSevChangeRules = sizeof (valid_sevchange) / sizeof (ValidErrSevChangeData);

static int AdjustSeverity (int severity, ValidStructPtr vsp, int code1, int code2, const char *fmt)
{
  Int4 i;
  int rval = severity;

  for (i = 0; i < kNumSevChangeRules; i++) {
    if (code1 == valid_sevchange[i].code1 && code2 == valid_sevchange[i].code2
        && (valid_sevchange[i].search_phrase == NULL || StringISearch (fmt, valid_sevchange[i].search_phrase) != NULL)
        && (valid_sevchange[i].exclude_phrase == NULL || StringISearch (fmt, valid_sevchange[i].exclude_phrase) == NULL)
        && valid_sevchange[i].func != NULL) {
      rval = (valid_sevchange[i].func)(rval, vsp);
    }
  }

  return rval;
}


#ifdef VAR_ARGS
NLM_EXTERN void CDECL ValidErr (vsp, severity, code1, code2, fmt, va_alist)
     ValidStructPtr vsp;
     int severity;
     int code1;
     int code2;
     const char     *fmt;
     va_dcl
#else
NLM_EXTERN void CDECL ValidErr (ValidStructPtr vsp, int severity, int code1, int code2, const char *fmt, ...)
#endif
{
  va_list           args;
  BioseqPtr         bsp;
  BioseqSetPtr      bssp;
  Int2              buflen, diff;
  CharPtr           ctmp, tmp;
  GatherContextPtr  gcp;
  Char              id [64];
  SeqLocPtr         loc = NULL;
  ObjValNodePtr     ovp;
  SeqDescrPtr       sdp;
  SeqEntryPtr       sep;
  SeqFeatPtr        sfp;
  SeqIdPtr          sip;

  if (vsp == NULL || severity < vsp->cutoff || ShouldSuppressValidErr(vsp, code1, code2, fmt))
    return;

  severity = AdjustSeverity(severity, vsp, code1, code2, fmt);

  if (vsp->errbuf == NULL) {
    vsp->errbuf = MemNew (8192);
    if (vsp->errbuf == NULL)
      AbnormalExit (1);
  }
  tmp = vsp->errbuf;

  vsp->errors[severity]++;

#ifdef VAR_ARGS
  va_start (args);
#else
  va_start (args, fmt);
#endif

  gcp = vsp->gcp;
  buflen = 1023;
  vsprintf (tmp, fmt, args);
  while (*tmp != '\0') {
    buflen--;
    tmp++;
  }

  va_end (args);

  if (vsp->errfunc != NULL) {
    CustValErr (vsp, (ErrSev) (severity), code1, code2);
    vsp->errbuf[0] = '\0';
    return;
  }

  if (vsp->justShowAccession) {
    vsp->errbuf[0] = '\0';
    tmp = vsp->errbuf;
    sip = NULL;

    if (vsp->sfp != NULL) {
      sfp = vsp->sfp;
      bsp = BioseqFindFromSeqLoc (sfp->location);
      if (bsp != NULL) {
        sip = SeqIdFindWorst (bsp->id);
      }
    } else if (vsp->descr != NULL) {
      sdp = vsp->descr;
      if (sdp != NULL && sdp->extended != 0) {
        ovp = (ObjValNodePtr) sdp;
        if (ovp->idx.parenttype == OBJ_BIOSEQ) {
          bsp = (BioseqPtr) ovp->idx.parentptr;
          if (bsp != NULL) {
            sip = SeqIdFindWorst (bsp->id);
          }
        } else if (ovp->idx.parenttype == OBJ_BIOSEQSET) {
          bssp = (BioseqSetPtr) ovp->idx.parentptr;
          if (bssp != NULL) {
            sep = bssp->seqentry;
            if (sep != NULL) {
              sep = FindNthBioseq (sep, 1);
              if (sep != NULL) {
                bsp = (BioseqPtr) sep->data.ptrvalue;
                if (bsp != NULL) {
                  sip = SeqIdFindWorst (bsp->id);
                }
              }
            }
          }
        }
      }
    } else if (vsp->bsp != NULL) {
      bsp = vsp->bsp;
      sip = SeqIdFindWorst (bsp->id);
    } else if (vsp->bssp != NULL) {
      bssp = vsp->bssp;
      sep = bssp->seqentry;
      if (sep != NULL) {
        sep = FindNthBioseq (sep, 1);
        if (sep != NULL) {
          bsp = (BioseqPtr) sep->data.ptrvalue;
          if (bsp != NULL) {
            sip = SeqIdFindWorst (bsp->id);
          }
        }
      }
    }

    if (sip != NULL) {
      SeqIdWrite (sip, id, PRINTID_REPORT, sizeof (id) - 1);
      diff = LabelCopy (tmp, id, buflen);
      buflen -= diff;
      tmp += diff;
    }

    ErrPostItem ((ErrSev) (severity), code1, code2, "%s", vsp->errbuf);
    vsp->errbuf[0] = '\0';
    return;
  }

  if (vsp->sfp != NULL) {
    diff = LabelCopy (tmp, " FEATURE: ", buflen);
    buflen -= diff;
    tmp += diff;

    diff = FeatDefLabel (vsp->sfp, tmp, buflen, OM_LABEL_BOTH);
    buflen -= diff;
    tmp += diff;

    if (vsp->suppressContext) {
      loc = AsnIoMemCopy (vsp->sfp->location, (AsnReadFunc) SeqLocAsnRead, (AsnWriteFunc) SeqLocAsnWrite);
      ChangeSeqLocToBestID (loc);
      ctmp = SeqLocPrint (loc);
      SeqLocFree (loc);
    } else {
      ctmp = SeqLocPrint (vsp->sfp->location);
    }
    if (ctmp != NULL && StringLen (ctmp) > 800) {
      StringCpy (ctmp + 797, "...");
    }
    if (ctmp != NULL) {
      diff = LabelCopyExtra (tmp, ctmp, buflen, " [", "]");
      buflen -= diff;
      tmp += diff;
      MemFree (ctmp);
    }

    if (!vsp->suppressContext) {
      sip = SeqLocId (vsp->sfp->location);
      if (sip != NULL) {
        bsp = BioseqFind (sip);
        if (bsp != NULL) {
          diff = LabelCopy (tmp, " [", buflen);
          buflen -= diff;
          tmp += diff;

          diff = BioseqLabel (bsp, tmp, buflen, OM_LABEL_BOTH);
          buflen -= diff;
          tmp += diff;

          diff = LabelCopy (tmp, "]", buflen);
          buflen -= diff;
          tmp += diff;
        }
      }
    }
    if (vsp->sfp->product != NULL) {
      if (vsp->suppressContext) {
        loc = AsnIoMemCopy (vsp->sfp->product, (AsnReadFunc) SeqLocAsnRead, (AsnWriteFunc) SeqLocAsnWrite);
        ChangeSeqLocToBestID (loc);
        ctmp = SeqLocPrint (loc);
        SeqLocFree (loc);
      } else {
        ctmp = SeqLocPrint (vsp->sfp->product);
      }
      if (ctmp != NULL && StringLen (ctmp) > 800) {
        StringCpy (ctmp + 797, "...");
      }
      if (ctmp != NULL) {
        diff = LabelCopyExtra (tmp, ctmp, buflen, " -> [", "]");
        buflen -= diff;
        tmp += diff;
        MemFree (ctmp);
      }
    }
  } else if (vsp->descr != NULL) {
    diff = LabelCopy (tmp, " DESCRIPTOR: ", buflen);
    buflen -= diff;
    tmp += diff;

    if (vsp->descr->choice == Seq_descr_comment) {
      diff = SeqDescLabel (vsp->descr, tmp, buflen, OM_LABEL_BOTH);
      if (diff > 100) {
        /* truncate long comment in message */
        tmp [94] = ' ';
        tmp [95] = '.';
        tmp [96] = '.';
        tmp [97] = '.';
        tmp [98] = '\0';
        diff = 98;
        buflen -= diff;
        tmp += diff;
      } else {
        buflen -= diff;
        tmp += diff;
      }
    } else {
      diff = SeqDescLabel (vsp->descr, tmp, buflen, OM_LABEL_BOTH);
      buflen -= diff;
      tmp += diff;
    }
  }

  /*
     if (vsp->suppressContext)
     {
     }
     else */
  if (vsp->sfp == NULL) {       /* sfp adds its own context */
    if (vsp->bsp != NULL) {
      diff = LabelCopy (tmp, " BIOSEQ: ", buflen);
      buflen -= diff;
      tmp += diff;

      if (vsp->bsp == NULL) {
        diff = LabelCopy (tmp, "??", buflen);
      } else if (vsp->suppressContext) {
        diff = WorstBioseqLabel (vsp->bsp, tmp, buflen, OM_LABEL_CONTENT);
      } else {
        diff = BioseqLabel (vsp->bsp, tmp, buflen, OM_LABEL_BOTH);
      }
      buflen -= diff;
      tmp += diff;
    } else if (vsp->bssp != NULL) {
      diff = LabelCopy (tmp, " BIOSEQ-SET: ", buflen);
      buflen -= diff;
      tmp += diff;

      if (vsp->suppressContext) {
        diff = BioseqSetLabel (vsp->bssp, tmp, buflen, OM_LABEL_CONTENT);
      } else {
        diff = BioseqSetLabel (vsp->bssp, tmp, buflen, OM_LABEL_BOTH);
      }
      buflen -= diff;
      tmp += diff;
    }
  }

  ErrPostItem ((ErrSev) (severity), code1, code2, "%s", vsp->errbuf);
  vsp->errbuf[0] = '\0';
}


static CharPtr GetUserFieldLabelString (UserFieldPtr ufp)
{
  Char buf[15];

  if (ufp == NULL || ufp->label == NULL) {
    return StringSave ("Unlabeled field");
  } else if (ufp->label->id > 0) {
    sprintf (buf, "%d", ufp->label->id);
    return StringSave (buf);
  } else {
    return StringSave (ufp->label->str);
  }
}


static CharPtr GetUserFieldValueString (UserFieldPtr ufp)
{
  Char buf[15];

  if (ufp == NULL) {
    return StringSave ("Value is missing");
  } else if (ufp->choice == 1) {
    return StringSave (ufp->data.ptrvalue);
  } else if (ufp->choice == 2) {
    sprintf (buf, "%d", ufp->data.intvalue);
    return StringSave (buf);
  } else {
    return StringSave ("Bad format for value");
  }
}


static ErrSev ErrorLevelFromFieldRuleSev (Uint2 severity)
{
  ErrSev sev = SEV_ERROR;
  switch (severity) {
    case Severity_level_none:
      sev = SEV_NONE;
      break;
    case Severity_level_info:
      sev = SEV_INFO;
      break;
    case Severity_level_warning:
      sev = SEV_WARNING;
      break;
    case Severity_level_error:
      sev = SEV_ERROR;
      break;
    case Severity_level_reject:
      sev = SEV_REJECT;
      break;
    case Severity_level_fatal:
      sev = SEV_FATAL;
      break;
  }
  return sev;
}

static void StructuredCommentError (EFieldValid err_code, FieldRulePtr field_rule, UserFieldPtr ufp, UserFieldPtr depend_ufp, Pointer data)
{
  ValidStructPtr     vsp;
  CharPtr            label, val;
  CharPtr            depend_label, depend_val, depend_str = NULL;
  CharPtr            depend_fmt = " when %s has value '%s'";
  ErrSev             sev = SEV_ERROR;

  if ((vsp = (ValidStructPtr) data) == NULL) {
    return;
  }

  if (field_rule != NULL) {
    sev = ErrorLevelFromFieldRuleSev(field_rule->severity);
  }

  if (depend_ufp != NULL) {
    depend_label = GetUserFieldLabelString  (depend_ufp);
    depend_val = GetUserFieldValueString (depend_ufp);
    depend_str = (CharPtr) MemNew (sizeof (Char) * (StringLen (depend_fmt) + StringLen (depend_label) + StringLen (depend_val)));
    sprintf (depend_str, depend_fmt, depend_label, depend_val);
    depend_val = MemFree (depend_val);
    depend_label = MemFree (depend_label);
  }

  switch (err_code) {
    case eFieldValid_Invalid:
      label = GetUserFieldLabelString  (ufp);
      if (field_rule == NULL) {
        ValidErr (vsp, sev, ERR_SEQ_DESCR_BadStrucCommInvalidFieldName, "%s is not a valid field name%s", label, depend_str == NULL ? "" : depend_str);
      } else {
        val = GetUserFieldValueString (ufp);
        ValidErr (vsp, sev, ERR_SEQ_DESCR_BadStrucCommInvalidFieldValue, "%s is not a valid value for %s%s", val, label, depend_str == NULL ? "" : depend_str);
        val = MemFree (val);
      }
      label = MemFree (label);
      break;
    case eFieldValid_MissingRequiredField:
      ValidErr (vsp, sev, ERR_SEQ_DESCR_BadStrucCommMissingField, "Required field %s is missing%s", field_rule->field_name, depend_str == NULL ? "" : depend_str);
      break;
    case eFieldValid_FieldOutOfOrder:
      ValidErr (vsp, sev, ERR_SEQ_DESCR_BadStrucCommFieldOutOfOrder, "%s field is out of order%s", field_rule->field_name, depend_str == NULL ? "" : depend_str);
      break;
    case eFieldValid_DuplicateField:
      ValidErr (vsp, sev, ERR_SEQ_DESCR_BadStrucCommMultipleFields, "Multiple values for %s field%s", field_rule->field_name, depend_str == NULL ? "" : depend_str);
      break;
    case eFieldValid_Disallowed:
      label = GetUserFieldLabelString  (ufp);
      ValidErr (vsp, sev, ERR_SEQ_DESCR_BadStrucCommInvalidFieldName, "%s is not a valid field name%s", label, depend_str == NULL ? "" : depend_str);
      label = MemFree (label);
      break;
    default:
      /* do nothing */
      break;
  }
  depend_str = MemFree (depend_str);
}


static Boolean StringLooksLikeFakeStructuredComment (CharPtr str)
{
  if (StringHasNoText (str)) {
    return FALSE;
  }
  if (StringSearch (str, "::") != NULL) {
    return TRUE;
  }
  return FALSE;
}


/*****************************************************************************
*
*   Valid1GatherProc(gcp)
*     top level gather callback
*     dispatches to other levels
*
*****************************************************************************/
static Boolean Valid1GatherProc (GatherContextPtr gcp)
{
  ValidStructPtr     vsp;
  UserFieldPtr       curr;
  AnnotDescrPtr      desc;
  CharPtr            field;
  SeqAnnotPtr        sap;
  Boolean            is_blast_align;
  Int2               limit;
  SeqFeatPtr         sfp;
  ValNodePtr         sdp;
  SeqGraphPtr        sgp;
  BioSourcePtr       biop;
  ObjectIdPtr        oip;
  PubdescPtr         pdp;
  CharPtr            ptr;
  BioseqPtr          bsp;
  SeqIdPtr           sip;
  CharPtr            str;
  Char               buf [64];
  Char               tmp [64];
  ValNodePtr         vnp2;
  SeqMgrFeatContext  context;
  UserObjectPtr      uop;
  EFieldValid        sc_valid;

  vsp = (ValidStructPtr) (gcp->userdata);
  vsp->gcp = gcp;               /* needed for ValidErr */

  limit = vsp->validationLimit;

  switch (gcp->thistype) {
  case OBJ_BIOSEQ:
    if (!vsp->onlyspell) {
      if (limit == VALIDATE_ALL || limit == VALIDATE_INST) {
        ValidateBioseqInst (gcp);
      }
      if (limit == VALIDATE_ALL || limit == VALIDATE_CONTEXT) {
        ValidateBioseqContext (gcp);
      }
      if (limit == VALIDATE_ALL || limit == VALIDATE_INST) {
        ValidateBioseqHist (gcp);
      }
      if (limit == VALIDATE_ALL || limit == VALIDATE_GRAPH) {
        ValidateGraphsOnBioseq (gcp);
      }
    }
    break;
  case OBJ_BIOSEQSET:
    if (!vsp->onlyspell) {
      if (limit == VALIDATE_ALL || limit == VALIDATE_SET) {
        ValidateBioseqSet (gcp);
      }
    }
    break;
  case OBJ_SEQANNOT:
    if (!vsp->onlyspell) {
      if (limit == VALIDATE_ALL) {
        sap = (SeqAnnotPtr) gcp->thisitem;
        if (sap != NULL) {
          if (sap->type == 2) {
            is_blast_align = FALSE;
            desc = NULL;
            while ((desc = ValNodeFindNext (sap->desc, desc, Annot_descr_user)) != NULL) {
              if (desc->data.ptrvalue != NULL) {
                oip = ((UserObjectPtr) desc->data.ptrvalue)->type;
                if (oip != NULL && StringCmp (oip->str, "Blast Type") == 0) {
                  is_blast_align = TRUE;
                }
              }
            }
            if (is_blast_align) {
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_ALIGN_BlastAligns, "Record contains BLAST alignments");
            }
          }
          if (sap->type == 4) {
            vsp->bssp = NULL;
            vsp->bsp = NULL;
            vsp->descr = NULL;
            vsp->sfp = NULL;
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_ANNOT_AnnotIDs, "Record contains Seq-annot.data.ids");
          }
          if (sap->type == 5) {
            vsp->bssp = NULL;
            vsp->bsp = NULL;
            vsp->descr = NULL;
            vsp->sfp = NULL;
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_ANNOT_AnnotLOCs, "Record contains Seq-annot.data.locs");
          }
        }
      }
    }
    break;
  case OBJ_SEQFEAT:
    if (!vsp->onlyspell) {
      if (limit == VALIDATE_ALL || limit == VALIDATE_FEAT) {
        ValidateSeqFeat (gcp);
        sfp = (SeqFeatPtr) (gcp->thisitem);
        if (sfp != NULL) {
          if (sfp->data.choice == SEQFEAT_BIOSRC) {
            biop = (BioSourcePtr) sfp->data.value.ptrvalue;
            ValidateBioSource (vsp, gcp, biop, sfp, NULL);
          }
          if (sfp->data.choice == SEQFEAT_PUB) {
            pdp = (PubdescPtr) sfp->data.value.ptrvalue;
            ValidatePubdesc (vsp, gcp, pdp);
          }
          if (sfp->cit != NULL) {
            ValidateSfpCit (vsp, gcp, sfp);
          }
          if (vsp->useSeqMgrIndexes) {
            if (SeqMgrGetDesiredFeature (gcp->entityID, NULL, 0, 0, sfp, &context) == NULL) {
              StringCpy (buf, "?");
              bsp = vsp->bsp;
              if (bsp != NULL) {
                SeqIdWrite (bsp->id, buf, PRINTID_FASTA_LONG, sizeof (buf) - 1);
              }
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_UnindexedFeature, "Feature is not indexed on Bioseq %s", buf);
            } else {
              bsp = BioseqFindFromSeqLoc (sfp->location);
              if (bsp != NULL) {
                sip = SeqLocId (sfp->location);
                if (sip != NULL && sip->choice != SEQID_GI && sip->choice != SEQID_GIBBSQ && sip->choice != SEQID_GIBBMT) {
                  SeqIdWrite (sip, buf, PRINTID_FASTA_SHORT, sizeof (buf) - 1);
                  for (sip = bsp->id; sip != NULL; sip = sip->next) {
                    if (sip->choice == SEQID_GI || sip->choice == SEQID_GIBBSQ || sip->choice == SEQID_GIBBMT) continue;
                    SeqIdWrite (sip, tmp, PRINTID_FASTA_SHORT, sizeof (tmp) - 1);
                    if (StringICmp (buf, tmp) != 0) continue;
                    if (StringCmp (buf, tmp) == 0) continue;
                    ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_FeatureSeqIDCaseDifference,
                              "Sequence identifier in feature location differs in capitalization with identifier on Bioseq");
                  }
                }
              }
            }
          }
        }
      }
    }
    if (limit == VALIDATE_ALL || limit == VALIDATE_FEAT) {
      SpellCheckSeqFeat (gcp);
    }
    break;
  case OBJ_SEQGRAPH :
    if (!vsp->onlyspell) {
      if (limit == VALIDATE_ALL || limit == VALIDATE_GRAPH) {
        sgp = (SeqGraphPtr) gcp->thisitem;
        if (sgp != NULL) {
          if (StringICmp (sgp->title, "Phrap Quality") == 0 ||
              StringICmp (sgp->title, "Phred Quality") == 0 ||
              StringICmp (sgp->title, "Gap4") == 0) {
            if (sgp->flags[2] == 3) {
              sip = SeqLocId (sgp->loc);
              if (sip != NULL) {
                if (BioseqFindCore (sip) == NULL) {
                  SeqIdWrite (sip, buf, PRINTID_FASTA_LONG, sizeof (buf) - 1);
                  ValidErr (vsp, SEV_WARNING, ERR_SEQ_GRAPH_GraphBioseqId, "Bioseq not found for Graph location %s", buf);
                }
              }
            }
          }
        }
      }
    }
    break;
  case OBJ_SEQDESC:
    if (limit == VALIDATE_ALL || limit == VALIDATE_DESC) {
      SpellCheckSeqDescr (gcp);
                          /**
              ValidateSeqDescr (gcp);
              **/
      sdp = (ValNodePtr) (gcp->thisitem);
      if (sdp != NULL) {
        if (sdp->choice == Seq_descr_source) {
          biop = (BioSourcePtr) sdp->data.ptrvalue;
          ValidateBioSource (vsp, gcp, biop, NULL, sdp);
        }
        if (sdp->choice == Seq_descr_pub) {
          pdp = (PubdescPtr) sdp->data.ptrvalue;
          ValidatePubdesc (vsp, gcp, pdp);
          LookForMultiplePubs (vsp, gcp, sdp);
        }
        if (sdp->choice == Seq_descr_user) {
          uop = sdp->data.ptrvalue;
          if (uop != NULL && uop->type != NULL && StringICmp (uop->type->str, "StructuredComment") == 0) {
            sc_valid = IsStructuredCommentValid (uop, StructuredCommentError, vsp);
            /* report ? */
            for (curr = uop->data; curr != NULL; curr = curr->next) {
              if (curr->choice != 1) continue;
              oip = curr->label;
              if (oip == NULL) continue;
              field = oip->str;
              if (StringStr (field, "::") != NULL) {
                ValidErr (vsp, SEV_REJECT, ERR_SEQ_DESCR_BadStrucCommInvalidFieldName, "Structured comment field '%s' contains double colons", field);
              }
              str = (CharPtr) curr->data.ptrvalue;
              if (StringStr (str, "::") != NULL) {
                ValidErr (vsp, SEV_REJECT, ERR_SEQ_DESCR_BadStrucCommInvalidFieldValue, "Structured comment value '%s' contains double colons", str);
              }
            }
          }
        }
        if (sdp->choice == Seq_descr_comment) {
          str = (CharPtr) sdp->data.ptrvalue;
          if (StringHasNoText (str)) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_MissingText, "Comment descriptor needs text");
          }
          if (SerialNumberInString (str)) {
            ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_SerialInComment,
                      "Comment may refer to reference by serial number - attach reference specific comments to the reference REMARK instead.");
          }
          if (StringLooksLikeFakeStructuredComment (str)) {
            ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_FakeStructuredComment,
                      "Comment may be formatted to look like a structured comment.");
          }
          for (vnp2 = sdp->next; vnp2 != NULL; vnp2 = vnp2->next) {
            if (vnp2->choice == Seq_descr_comment) {
              ptr = (CharPtr) vnp2->data.ptrvalue;
              if (StringDoesHaveText (ptr) && StringICmp (str, ptr) == 0) {
                ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MultipleComments, "Undesired multiple comment descriptors, identical text");
              }
            }
          }
        }
        if (sdp->choice == Seq_descr_mol_type) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_InvalidForType, "MolType descriptor is obsolete");
        }
        if (sdp->choice == Seq_descr_modif) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_InvalidForType, "Modif descriptor is obsolete");
        }
        if (sdp->choice == Seq_descr_method) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_InvalidForType, "Method descriptor is obsolete");
        }
        if (sdp->choice == Seq_descr_org) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_InvalidForType, "OrgRef descriptor is obsolete");
        }
      }
    }
    break;
  default:
    break;

  }
  return TRUE;
}


static void DiscrepanciesToValidationErrs (ValNodePtr discrepancy_list, Uint4 item_type, ValidStructPtr vsp, int severity, int code1, int code2, char *msg)
{
  ValNodePtr vnp, obj;
  ClickableItemPtr cip;

  if (discrepancy_list == NULL || vsp == NULL) {
    return;
  }
  vsp->bssp = NULL;
  vsp->bsp = NULL;
  vsp->sfp = NULL;
  vsp->descr = NULL;
  for (vnp = discrepancy_list; vnp != NULL; vnp = vnp->next) {
    cip = (ClickableItemPtr) vnp->data.ptrvalue;
    if (cip != NULL) {
      if (cip->clickable_item_type == item_type) {
        if (cip->item_list == NULL) {
          DiscrepanciesToValidationErrs (cip->subcategories, item_type, vsp, severity, code1, code2, msg);
        } else {
          for (obj = cip->item_list; obj != NULL; obj = obj->next) {
            if (obj->choice == OBJ_SEQFEAT) {
              vsp->sfp = obj->data.ptrvalue;
              vsp->gcp->entityID = vsp->sfp->idx.entityID;
              vsp->gcp->thistype = OBJ_SEQFEAT;
              vsp->gcp->itemID = vsp->sfp->idx.itemID;

              ValidErr (vsp, severity, code1, code2, msg);
              vsp->sfp = NULL;
            }
          }
        }
      }
    }
  }
  vsp->sfp = NULL;
}


static void ValidateGeneLocusTags (SeqEntryPtr sep, ValidStructPtr vsp)
{
  ValNode vn;
  ValNodePtr discrepancy_list = NULL;

  if (sep == NULL || vsp == NULL) {
    return;
  }

  vn.choice = 0;
  vn.data.ptrvalue = sep;
  vn.next = NULL;

  AddDiscrepanciesForMissingOrNonUniqueGeneLocusTagsEx (&discrepancy_list, &vn, TRUE);

  DiscrepanciesToValidationErrs (discrepancy_list, DISC_GENE_MISSING_LOCUS_TAG, vsp, SEV_WARNING, ERR_SEQ_FEAT_MissingGeneLocusTag, "Missing gene locus tag");

  discrepancy_list = FreeClickableList (discrepancy_list);
}


static void ValidateShortIntrons (SeqEntryPtr sep, ValidStructPtr vsp)
{
  ValNode vn;
  ValNodePtr discrepancy_list = NULL;

  if (sep == NULL || vsp == NULL) {
    return;
  }

  vn.choice = 0;
  vn.data.ptrvalue = sep;
  vn.next = NULL;

  FindShortIntrons (&discrepancy_list, &vn);

  DiscrepanciesToValidationErrs (discrepancy_list, DISC_SHORT_INTRON, vsp, SEV_WARNING, ERR_SEQ_FEAT_ShortIntron, "Introns should be at least 10 nt long");

  discrepancy_list = FreeClickableList (discrepancy_list);
}


static void LookForAnyPubAndOrg (SeqEntryPtr sep, BoolPtr no_pub, BoolPtr no_biosrc)
{
  BioseqPtr       bsp;
  BioseqSetPtr    bssp;
  SeqAnnotPtr     sap = NULL;
  ValNodePtr      sdp = NULL;
  SeqFeatPtr      sfp;
  SeqEntryPtr     tmp;

  if (sep == NULL || no_pub == NULL || no_biosrc == NULL)
    return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL)
      return;
    sap = bsp->annot;
    sdp = bsp->descr;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL)
      return;
    for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
      LookForAnyPubAndOrg (tmp, no_pub, no_biosrc);
    }
    sap = bssp->annot;
    sdp = bssp->descr;
  } else
    return;
  while (sap != NULL) {
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      while (sfp != NULL) {
        if (sfp->data.choice == SEQFEAT_PUB) {
          *no_pub = FALSE;
        } else if (sfp->data.choice == SEQFEAT_BIOSRC) {
          *no_biosrc = FALSE;
        }
        sfp = sfp->next;
      }
    }
    sap = sap->next;
  }
  while (sdp != NULL) {
    if (sdp->choice == Seq_descr_pub) {
      *no_pub = FALSE;
    } else if (sdp->choice == Seq_descr_source) {
      *no_biosrc = FALSE;
    }
    sdp = sdp->next;
  }
}

typedef struct ftprob {
  Uint4    num_misplaced_features;
  Uint4    num_small_genome_set_misplaced;
  Uint4    num_archaic_locations;
  Uint4    num_archaic_products;
  Uint4    num_misplaced_graphs;
  Uint4    num_gene_feats;
  Uint4    num_gene_xrefs;
  Uint4    num_tpa_with_hist;
  Uint4    num_tpa_without_hist;
  Boolean  has_gi;
  Boolean  loc_has_gi;
  Boolean  loc_has_just_accn;
  Boolean  loc_has_accn_ver;
  Boolean  prod_has_gi;
  Boolean  prod_has_just_accn;
  Boolean  prod_has_accn_ver;
} FeatProb, PNTR FeatProbPtr;

static void CheckFeatPacking (BioseqPtr bsp, SeqFeatPtr sfp, Uint4Ptr num_misplaced_features, Uint4Ptr num_small_genome_set_misplaced)
{
  SeqAnnotPtr     sap;
  BioseqSetPtr    bssp, parent;
  BioseqPtr       par;

  if (sfp->idx.parenttype == OBJ_SEQANNOT) {
    sap = (SeqAnnotPtr) sfp->idx.parentptr;
    if (sap == NULL)
      return;
    if (sap->idx.parenttype == OBJ_BIOSEQ) {
      /* if feature packaged on bioseq, must be target bioseq */
      par = (BioseqPtr) sap->idx.parentptr;
      if (par != bsp && SeqMgrGetParentOfPart (par, NULL) != bsp) {
        /* generated gap feature is an exception */
        if (par == NULL || par->id != NULL) {
          (*num_misplaced_features)++;
        }
      }
      return;
    }
    if (sap->idx.parenttype == OBJ_BIOSEQSET) {
      /* if feature packaged on set, set must contain bioseq */
      bssp = (BioseqSetPtr) sap->idx.parentptr;
      if (bssp == NULL)
        return;
      if (bsp->idx.parenttype == OBJ_BIOSEQSET) {
        parent = (BioseqSetPtr) bsp->idx.parentptr;
        while (parent != NULL) {
          if (parent == bssp) return;
          if (parent->idx.parenttype == OBJ_BIOSEQSET && parent->_class == BioseqseqSet_class_small_genome_set) {
            (*num_small_genome_set_misplaced)++;
            return;
          }
          if (parent->idx.parenttype != OBJ_BIOSEQSET) {
            (*num_misplaced_features)++;
            return;
          }
          parent = (BioseqSetPtr) parent->idx.parentptr;
        }
        (*num_misplaced_features)++;
      }
    }
  }
}

static Boolean IdIsArchaic (SeqIdPtr sip)

{
  BioseqPtr  bsp;
  DbtagPtr   dbt;
  SeqIdPtr   id;

  if (sip == NULL) return FALSE;
  if (sip->choice != SEQID_LOCAL && sip->choice != SEQID_GENERAL) return FALSE;
  bsp = BioseqFind (sip);
  if (bsp == NULL) return FALSE;
  for (id = bsp->id; id != NULL; id = id->next) {
    switch (id->choice) {
      case SEQID_GENERAL :
        if (sip->choice == SEQID_LOCAL) {
          dbt = (DbtagPtr) id->data.ptrvalue;
          if (dbt != NULL && !IsSkippableDbtag(dbt)) {
            return TRUE;
          }
        }
        break;
      case SEQID_GI :
      case SEQID_GENBANK :
      case SEQID_EMBL :
      case SEQID_PATENT :
      case SEQID_OTHER :
      case SEQID_DDBJ :
      case SEQID_TPG :
      case SEQID_TPE :
      case SEQID_TPD :
      case SEQID_GPIPE :
        return TRUE;
      default :
        break;
    }
  }
  return FALSE;
}

static void CheckFeatLocAndProd (SeqFeatPtr sfp, FeatProbPtr fpp)

{
  SeqLocPtr  slp;

  if (sfp == NULL || fpp == NULL) return;
  if (sfp->product != NULL && IdIsArchaic (SeqLocId (sfp->product))) {
    (fpp->num_archaic_products)++;
  }
  slp = SeqLocFindNext (sfp->location, NULL);
  while (slp != NULL) {
    if (IdIsArchaic (SeqLocId (slp))) {
      (fpp->num_archaic_locations)++;
      return;
    }
    slp = SeqLocFindNext (sfp->location, slp);
  }
}

static void CheckGraphPacking (SeqGraphPtr sgp, Pointer userdata)

{
  BioseqPtr    bsp;
  FeatProbPtr  fpp;
  SeqAnnotPtr  sap;
  BioseqPtr    par;

  if (sgp == NULL || userdata == NULL) return;
  fpp = (FeatProbPtr) userdata;
  bsp = BioseqFindFromSeqLoc (sgp->loc);
  if (sgp->idx.parenttype == OBJ_SEQANNOT) {
    sap = (SeqAnnotPtr) sgp->idx.parentptr;
    if (sap == NULL) return;
    if (sap->idx.parenttype == OBJ_BIOSEQ) {
      /* if graph packaged on bioseq, must be target bioseq */
      par = (BioseqPtr) sap->idx.parentptr;
      if (par != bsp && SeqMgrGetParentOfPart (par, NULL) != bsp) {
        (fpp->num_misplaced_graphs)++;
      }
      return;
    }
    (fpp->num_misplaced_graphs)++;
  }
}

static Boolean LIBCALLBACK CountMisplacedFeatures (BioseqPtr bsp, SeqMgrBioseqContextPtr bcontext)

{
  SeqMgrFeatContext  fcontext;
  FeatProbPtr        fpp;
  SeqFeatPtr         sfp;

  fpp = (FeatProbPtr) bcontext->userdata;
  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  while (sfp != NULL) {
    if (! fcontext.ts_image) {
      CheckFeatPacking (bsp, sfp, &(fpp->num_misplaced_features), &(fpp->num_small_genome_set_misplaced));
      CheckFeatLocAndProd (sfp, fpp);
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
  }

  return TRUE;
}

static void CountGeneXrefs (SeqFeatPtr sfp, Pointer userdata)

{
  FeatProbPtr  fpp;
  GeneRefPtr   grp;

  if (sfp == NULL || userdata == NULL) return;
  fpp = (FeatProbPtr) userdata;

  if (sfp->data.choice == SEQFEAT_GENE) {
    (fpp->num_gene_feats)++;
  }

  grp = SeqMgrGetGeneXref (sfp);
  if (grp == NULL || SeqMgrGeneIsSuppressed (grp)) return;

  (fpp->num_gene_xrefs)++;
}

static void CountSfpLocIdTypes (SeqIdPtr sip, Pointer userdata)

{
  FeatProbPtr   fpp;
  TextSeqIdPtr  tsip;

  if (sip == NULL || userdata == NULL) return;
  fpp = (FeatProbPtr) userdata;

  switch (sip->choice) {
    case SEQID_GI :
      fpp->loc_has_gi = TRUE;
      break;
    case SEQID_GENBANK :
    case SEQID_EMBL :
    case SEQID_DDBJ :
    case SEQID_TPG :
    case SEQID_TPE :
    case SEQID_TPD :
    case SEQID_OTHER :
      tsip = (TextSeqIdPtr) sip->data.ptrvalue;
      if (tsip != NULL) {
        if (StringDoesHaveText (tsip->accession)) {
          if (tsip->version < 1) {
            fpp->loc_has_just_accn = TRUE;
          } else {
            fpp->loc_has_accn_ver = TRUE;
          }
        }
      }
      break;
    default :
      break;
  }
}

static void CountSfpProdIdTypes (SeqIdPtr sip, Pointer userdata)

{
  FeatProbPtr   fpp;
  TextSeqIdPtr  tsip;

  if (sip == NULL || userdata == NULL) return;
  fpp = (FeatProbPtr) userdata;

  switch (sip->choice) {
    case SEQID_GI :
      fpp->prod_has_gi = TRUE;
      break;
    case SEQID_GENBANK :
    case SEQID_EMBL :
    case SEQID_DDBJ :
    case SEQID_TPG :
    case SEQID_TPE :
    case SEQID_TPD :
    case SEQID_OTHER :
      tsip = (TextSeqIdPtr) sip->data.ptrvalue;
      if (tsip != NULL) {
        if (StringDoesHaveText (tsip->accession)) {
          if (tsip->version < 1) {
            fpp->prod_has_just_accn = TRUE;
          } else {
            fpp->prod_has_accn_ver = TRUE;
          }
        }
      }
      break;
    default :
      break;
  }
}

static void CountFeatLocIdTypes (SeqFeatPtr sfp, Pointer userdata)

{
  if (sfp == NULL || userdata == NULL) return;

  VisitSeqIdsInSeqLoc (sfp->location, userdata, CountSfpLocIdTypes);
  VisitSeqIdsInSeqLoc (sfp->product, userdata, CountSfpProdIdTypes);
}

NLM_EXTERN Boolean HasTpaUserObject (BioseqPtr bsp)

{
  SeqMgrDescContext  context;
  UserObjectPtr      uop;
  ObjectIdPtr        oip;
  ValNodePtr         vnp;

  if (bsp == NULL) return FALSE;
  vnp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &context);
  while (vnp != NULL) {
    uop = (UserObjectPtr) vnp->data.ptrvalue;
    if (uop != NULL) {
      oip = uop->type;
      if (oip != NULL && StringICmp (oip->str, "TpaAssembly") == 0) return TRUE;
    }
    vnp = SeqMgrGetNextDescriptor (bsp, vnp, Seq_descr_user, &context);
  }
  return FALSE;
}

static void CheckTpaHist (BioseqPtr bsp, Pointer userdata)

{
  FeatProbPtr  fpp;
  SeqHistPtr   shp;
  SeqIdPtr     sip;

  if (bsp == NULL || userdata == NULL) return;
  fpp = (FeatProbPtr) userdata;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GI) {
      fpp->has_gi = TRUE;
    }
  }
  if (! HasTpaUserObject (bsp)) return;
  shp = bsp->hist;
  if (shp != NULL && shp->assembly != NULL) {
    (fpp->num_tpa_with_hist)++;
  } else {
    (fpp->num_tpa_without_hist)++;
  }
}

static Boolean IsNoncuratedRefSeq (BioseqPtr bsp, ErrSev *sev)

{
  SeqIdPtr      sip;
  TextSeqIdPtr  tsip;

  if (bsp == NULL) return FALSE;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_OTHER) {
      tsip = (TextSeqIdPtr) sip->data.ptrvalue;
      if (tsip != NULL && tsip->accession != NULL) {
        if (StringNCmp (tsip->accession, "NM_", 3) == 0 ||
            StringNCmp (tsip->accession, "NP_", 3) == 0 ||
            StringNCmp (tsip->accession, "NG_", 3) == 0 ||
            StringNCmp (tsip->accession, "NR_", 3) == 0) {
          *sev = SEV_WARNING;
          return FALSE;
        }
        return TRUE;
      }
    }
  }
  return FALSE;
}

static Boolean IsGpipe (BioseqPtr bsp)

{
  SeqIdPtr  sip;

  if (bsp == NULL) return FALSE;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GPIPE) return TRUE;
  }
  return FALSE;
}

typedef struct vfcdata {
  ValNodePtr      uids;
  ValNodePtr      unpub;
  ValNodePtr      publshd;
  ValNodePtr      serial;
  ValidStructPtr  vsp;
} VfcData, PNTR VfcPtr;

static Boolean SkipSerialOrUIDPub (ValNodePtr vnp)

{
  CitGenPtr  cgp;

  if (vnp == NULL || vnp->next == NULL) return FALSE;
  if (vnp->choice == PUB_Muid || vnp->choice == PUB_Muid) return TRUE;
  if (vnp->choice != PUB_Gen) return FALSE;
  cgp = (CitGenPtr) vnp->data.ptrvalue;
  if (cgp == NULL) return FALSE;
  if (StringNICmp ("BackBone id_pub", cgp->cit, 15) == 0) return FALSE;
  if (cgp->cit == NULL && cgp->journal == NULL && cgp->date == NULL && cgp->serial_number) return TRUE;
  return FALSE;
}

static void MakePubTags (PubdescPtr pdp, Pointer userdata)

{
  Char        buf [1024];
  CitGenPtr   cgp;
  Int4        muid = 0, pmid = 0;
  VfcPtr      vfp;
  ValNodePtr  vnp, tmp;

  if (pdp == NULL || userdata == NULL) return;
  vfp = (VfcPtr) userdata;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == PUB_Muid) {
      muid = vnp->data.intvalue;
    } else if (vnp->choice == PUB_PMid) {
      pmid = vnp->data.intvalue;
    } else if (vnp->choice == PUB_Gen) {
      cgp = (CitGenPtr) vnp->data.ptrvalue;
      if (cgp != NULL && cgp->serial_number > 0) {
        tmp = ValNodeNew (NULL);
        if (tmp != NULL) {
          tmp->data.intvalue = (Int4) cgp->serial_number;
          tmp->next = vfp->serial;
          vfp->serial = tmp;
        }
      }
    }
  }

  if (pmid != 0) {
    vnp = ValNodeNew (NULL);
    if (vnp != NULL) {
      vnp->choice = 1;
      vnp->data.intvalue = pmid;
      vnp->next = vfp->uids;
      vfp->uids = vnp;
    }
  }
  if (muid != 0) {
    vnp = ValNodeNew (NULL);
    if (vnp != NULL) {
      vnp->choice = 2;
      vnp->data.intvalue = muid;
      vnp->next = vfp->uids;
      vfp->uids = vnp;
    }
  }

  vnp = pdp->pub;
  while (vnp != NULL && SkipSerialOrUIDPub (vnp)) {
    vnp = vnp->next;
  }
  if (vnp != NULL && PubLabelUnique (vnp, buf, sizeof (buf) - 1, OM_LABEL_CONTENT, TRUE) > 0) {
    tmp = ValNodeCopyStr (NULL, 0, buf);
    if (tmp != NULL) {
      if (pmid != 0 || muid != 0) {
        tmp->next = vfp->publshd;
        vfp->publshd = tmp;
      } else {
        tmp->next = vfp->unpub;
        vfp->unpub = tmp;
      }
    }
  }
}

static void CheckOneCit (SeqFeatPtr sfp, ValNodePtr ppr, VfcPtr vfp)

{
  Char              buf [1024];
  GatherContextPtr  gcp;
  size_t            len, lgth;
  CharPtr           str;
  Int4              uid;
  ValNodePtr        vnp;
  ValidStructPtr    vsp;

  if (sfp == NULL || ppr == NULL || vfp == NULL) return;
  vsp = vfp->vsp;
  if (vsp == NULL) return;
  gcp = vsp->gcp;

  if (gcp != NULL) {
    gcp->entityID = sfp->idx.entityID;
    gcp->itemID = sfp->idx.itemID;
    gcp->thistype = OBJ_SEQFEAT;
  }
  vsp->sfp = sfp;

  if (ppr->choice == PUB_PMid || ppr->choice == PUB_Muid) {
    uid = ppr->data.intvalue;
    for (vnp = vfp->uids; vnp != NULL; vnp = vnp->next) {
      if (uid == vnp->data.intvalue) return;
    }
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FeatureCitationProblem,
              "Citation on feature refers to uid [%ld] not on a publication in the record", (long) uid);
    vsp->sfp = NULL;

  } else if (ppr->choice == PUB_Equiv) {
    return;

  } else {
    PubLabelUnique (ppr, buf, sizeof (buf) - 1, OM_LABEL_CONTENT, TRUE);
    lgth = StringLen (buf);
    if (lgth > 0 && buf [lgth - 1] == '>') {
      buf [lgth - 1] = '\0';
     lgth--;
    }
    for (vnp = vfp->unpub; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      len = MIN (lgth, StringLen (str));
      if (StringNICmp (str, buf, len) == 0) return;
    }
    for (vnp = vfp->publshd; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      len = MIN (lgth, StringLen (str));
      if (StringNICmp (str, buf, len) == 0) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FeatureCitationProblem,
                  "Citation on feature needs to be updated to published uid");
        vsp->sfp = NULL;
        return;
      }
    }
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FeatureCitationProblem,
              "Citation on feature refers to a publication not in the record");
    vsp->sfp = NULL;
  }
}

static void CheckFeatCits (SeqFeatPtr sfp, Pointer userdata)

{
  ValNodePtr  ppr, vnp;
  VfcPtr      vfp;

  if (sfp == NULL || sfp->cit == NULL || userdata == NULL) return;
  vfp = (VfcPtr) userdata;

  vnp = sfp->cit;
  for (ppr = vnp->data.ptrvalue; ppr != NULL; ppr = ppr->next) {
    CheckOneCit (sfp, ppr, vfp);
  }
}

static void CheckForCollidingSerials (
  ValidStructPtr vsp,
  GatherContextPtr gcp,
  ValNodePtr list
)

{
  Int4        curr, last;
  Uint2       olditemtype = 0;
  Uint4       olditemid = 0;
  ValNodePtr  vnp, vnp_next;

  if (vsp == NULL || gcp == NULL || list == NULL) return;

  olditemid = gcp->itemID;
  olditemtype = gcp->thistype;
  gcp->itemID = 0;
  gcp->thistype = 0;

  last = (Int4) list->data.intvalue;
  for (vnp = list->next; vnp != NULL; vnp = vnp_next) {
    vnp_next = vnp->next;
    curr = (Int4) vnp->data.intvalue;
    if (last == curr) {
      ValidErr (vsp, SEV_WARNING, ERR_GENERIC_CollidingSerialNumbers,
                "Multiple publications have serial number %ld", (long) curr);
      while (vnp != NULL && vnp->data.intvalue == last) {
        vnp = vnp->next;
      }
      if (vnp == NULL) {
        vnp_next = NULL;
      } else {
        last = vnp->data.intvalue;
        vnp_next = vnp->next;
      }
    } else {
      last = curr;
    }
  }

  gcp->itemID = olditemid;
  gcp->thistype = olditemtype;
}

static void ValidateFeatCits (SeqEntryPtr sep, ValidStructPtr vsp)

{
  SeqEntryPtr    bsep;
  BioseqPtr      bsp = NULL;
  GatherContext  gc;
  VfcData        vfd;

  if (vsp == NULL || sep == NULL) return;

  bsep = FindNthBioseq (sep, 1);
  if (bsep != NULL && IS_Bioseq (bsep)) {
    bsp = (BioseqPtr) bsep->data.ptrvalue;
  }

  vsp->gcp = &gc;
  vsp->bssp = NULL;
  vsp->bsp = bsp;
  vsp->sfp = NULL;
  vsp->descr = NULL;
  MemSet ((Pointer) &gc, 0, sizeof (GatherContext));
  MemSet ((Pointer) &vfd, 0, sizeof (VfcData));
  vfd.vsp = vsp;

  VisitPubdescsInSep (sep, (Pointer) &vfd, MakePubTags);

  VisitFeaturesInSep (sep, (Pointer) &vfd, CheckFeatCits);

  vsp->bssp = NULL;
  vsp->bsp = bsp;
  vsp->sfp = NULL;
  vsp->descr = NULL;
  vfd.serial = ValNodeSort (vfd.serial, SortByIntvalue);
  CheckForCollidingSerials (vsp, vsp->gcp, vfd.serial);

  ValNodeFree (vfd.uids);
  ValNodeFreeData (vfd.unpub);
  ValNodeFreeData (vfd.publshd);
  ValNodeFree (vfd.serial);
}

static void ValidateFeatIDs (SeqEntryPtr sep, Uint2 entityID, ValidStructPtr vsp)

{
  SMFidItemPtr PNTR  array;
 SeqEntryPtr         bsep;
  BioseqPtr           bsp = NULL;
  BioseqExtraPtr     bspextra;
  SMFeatItemPtr      feat;
  GatherContext      gc;
  GatherContextPtr   gcp;
  SMFidItemPtr       item;
  Int4               j;
  CharPtr            last = NULL;
  Int4               num;
  ObjMgrDataPtr      omdp;
  SeqFeatPtr         sfp;

  if (sep == NULL || entityID < 1 || vsp == NULL) return;
  omdp = ObjMgrGetData (entityID);
  if (omdp == NULL) return;
  bspextra = (BioseqExtraPtr) omdp->extradata;
  if (bspextra == NULL) return;
  array = bspextra->featsByFeatID;
  num = bspextra->numfids;
  if (array == NULL || num < 1) return;

  bsep = FindNthBioseq (sep, 1);
  if (bsep != NULL && IS_Bioseq (bsep)) {
    bsp = (BioseqPtr) bsep->data.ptrvalue;
  }

  vsp->gcp = &gc;
  vsp->bssp = NULL;
  vsp->bsp = bsp;
  vsp->sfp = NULL;
  vsp->descr = NULL;
  MemSet ((Pointer) &gc, 0, sizeof (GatherContext));

  for (j = 0; j < num; j++) {
    item = array [j];
    if (item == NULL) continue;
    if (StringDoesHaveText (last)) {
      if (StringICmp (item->fid, last) == 0) {
        feat = item->feat;
        if (feat == NULL) continue;
        sfp = feat->sfp;
        if (sfp == NULL) continue;
        gcp = &gc;
        gcp->entityID = sfp->idx.entityID;
        gcp->itemID = sfp->idx.itemID;
        gcp->thistype = OBJ_SEQFEAT;
        vsp->sfp = sfp;
        ValidErr (vsp, SEV_REJECT, ERR_SEQ_FEAT_CollidingFeatureIDs,
                  "Colliding feature ID %s", last);
      }
    }
    last = item->fid;
  }
}

typedef struct vsicdata {
  ValidStructPtr  vsp;
  ValNodePtr      headid;
  ValNodePtr      tailid;
} VsicData, PNTR VsicDataPtr;

static void CaptureTextSeqIDs (BioseqPtr bsp, Pointer userdata)

{
  Char         buf [200];
  SeqIdPtr     sip;
  VsicDataPtr  vdp;
  ValNodePtr   vnp;

  if (bsp == NULL || userdata == NULL) return;
  vdp = (VsicDataPtr) userdata;

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GI || sip->choice == SEQID_GIBBSQ || sip->choice == SEQID_GIBBMT) continue;
    if (IsNCBIFileID (sip)) continue;
    SeqIdWrite (sip, buf, PRINTID_FASTA_SHORT, sizeof (buf) - 1);
    vnp = ValNodeCopyStr (&(vdp->tailid), 0, buf);
    if (vdp->headid == NULL) {
      vdp->headid = vnp;
    }
    vdp->tailid = vnp;
  }
}

static ValNodePtr UniqueValNodeCaseSensitive (ValNodePtr list)

{
  CharPtr       last;
  ValNodePtr    next;
  Pointer PNTR  prev;
  CharPtr       str;
  ValNodePtr    vnp;

  if (list == NULL) return NULL;
  last = (CharPtr) list->data.ptrvalue;
  vnp = list->next;
  prev = (Pointer PNTR) &(list->next);
  while (vnp != NULL) {
    next = vnp->next;
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringCmp (last, str) == 0) {
      vnp->next = NULL;
      *prev = next;
      ValNodeFreeData (vnp);
    } else {
      last = (CharPtr) vnp->data.ptrvalue;
      prev = (Pointer PNTR) &(vnp->next);
    }
    vnp = next;
  }

  return list;
}

static void ValidateSeqIdCase (SeqEntryPtr sep, ValidStructPtr vsp)

{
  SeqEntryPtr       bsep;
  BioseqPtr         bsp = NULL;
  CharPtr           curr;
  GatherContext     gc;
  GatherContextPtr  gcp;
  CharPtr           prev;
  VsicData          vd;
  ValNodePtr        vnp;

  if (vsp == NULL || sep == NULL) return;

  bsep = FindNthBioseq (sep, 1);
  if (bsep != NULL && IS_Bioseq (bsep)) {
    bsp = (BioseqPtr) bsep->data.ptrvalue;
  }

  MemSet ((Pointer) &gc, 0, sizeof (GatherContext));
  MemSet ((Pointer) &vd, 0, sizeof (VsicData));

  gcp = &gc;
  vsp->gcp = &gc;
  vsp->bssp = NULL;
  vsp->bsp = bsp;
  vsp->sfp = NULL;
  vsp->descr = NULL;
  vd.vsp = vsp;

  VisitBioseqsInSep (sep, (Pointer) &vd, CaptureTextSeqIDs);
  vd.headid = ValNodeSort (vd.headid, SortVnpByString);
  vd.headid = UniqueValNodeCaseSensitive (vd.headid);

  curr = NULL;
  prev = NULL;
  for (vnp = vd.headid; vnp != NULL; vnp = vnp->next, prev = curr) {
    curr = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (curr)) continue;
    if (StringHasNoText (prev)) continue;
    if (StringICmp (curr, prev) != 0) continue;
    if (StringCmp (curr, prev) == 0) continue;
    ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_CaseDifferenceInSeqID,
              "Sequence identifier differs only by case - %s and %s", curr, prev);
  }

  vsp->bssp = NULL;
  vsp->bsp = NULL;
  vsp->sfp = NULL;
  vsp->descr = NULL;

  ValNodeFreeData (vd.headid);
}

static void LookForBioseqFields (BioseqPtr bsp, Pointer userdata)

{
  DbtagPtr        dbt;
  Boolean         has_lcl_gnl = FALSE;
  Boolean         has_others = FALSE;
  SeqIdPtr        sip;
  TextSeqIdPtr    tsip;
  ValidStructPtr  vsp;

  if (bsp == NULL || userdata == NULL) return;
  vsp = (ValidStructPtr) userdata;

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    switch (sip->choice) {
    case SEQID_EMBL:
    case SEQID_DDBJ:
      vsp->is_embl_ddbj_in_sep = TRUE;
      /* and fall through */
    case SEQID_GENBANK:
    case SEQID_TPG:
      vsp->is_insd_in_sep = TRUE;
      tsip = (TextSeqIdPtr) sip->data.ptrvalue;
      if (tsip != NULL) {
        if (StringLen (tsip->accession) == 6) {
          vsp->is_old_gb_in_sep = TRUE;
        }
      }
      break;
    case SEQID_TPE:
    case SEQID_TPD:
      vsp->is_insd_in_sep = TRUE;
      break;
    case SEQID_PATENT:
      vsp->is_patent_in_sep = TRUE;
      break;
    case SEQID_OTHER:
      vsp->is_refseq_in_sep = TRUE;
      break;
    case SEQID_GPIPE:
      vsp->is_gpipe_in_sep = TRUE;
      break;
    case SEQID_GENERAL:
      if (ISA_aa (bsp->mol)) {
        dbt = (DbtagPtr) sip->data.ptrvalue;
        if (dbt == NULL) break;
        if (IsSkippableDbtag (dbt)) break;
        vsp->has_gnl_prot_sep = TRUE;
      }
      break;
    default:
      break;
    }
    if (sip->choice == SEQID_LOCAL || sip->choice == SEQID_GENERAL) {
      has_lcl_gnl = TRUE;
    } else {
      has_others = TRUE;
    }
  }
  if (has_lcl_gnl && ! has_others) {
    vsp->only_lcl_gnl_in_sep = TRUE;
  }
}

static void LookForBioseqSetFields (BioseqSetPtr bssp, Pointer userdata)

{
  ValidStructPtr  vsp;

  if (bssp == NULL || userdata == NULL) return;
  vsp = (ValidStructPtr) userdata;

  /* the switch statement is erroneously reporting gen_prod_set for a pop_set under Xcode */
  /*
  switch (bssp->_class) {
  case BioseqseqSet_class_gen_prod_set:
    vsp->is_gps_in_sep = TRUE;
    break;
  case BioseqseqSet_class_mut_set:
  case BioseqseqSet_class_pop_set:
  case BioseqseqSet_class_phy_set:
  case BioseqseqSet_class_eco_set:
  case BioseqseqSet_class_wgs_set:
  case BioseqseqSet_class_small_genome_set:
    break;
    vsp->other_sets_in_sep = TRUE;
  default:
    break;
  }
  */

  if (bssp->_class == BioseqseqSet_class_gen_prod_set) {
    vsp->is_gps_in_sep = TRUE;
  } else if (bssp->_class == BioseqseqSet_class_mut_set ||
             bssp->_class == BioseqseqSet_class_pop_set ||
             bssp->_class == BioseqseqSet_class_phy_set ||
             bssp->_class == BioseqseqSet_class_eco_set ||
             bssp->_class == BioseqseqSet_class_wgs_set) {
    vsp->other_sets_in_sep = TRUE;
  } else if (bssp->_class == BioseqseqSet_class_small_genome_set) {
    vsp->is_small_genome_set = TRUE;
  }
}

static void LookForSeqDescrFields (SeqDescrPtr sdp, Pointer userdata)

{
  BioSourcePtr    biop;
  MolInfoPtr      mip;
  ObjectIdPtr     oip;
  UserObjectPtr   uop;
  ValidStructPtr  vsp;

  if (sdp == NULL || userdata == NULL) return;
  vsp = (ValidStructPtr) userdata;

  switch (sdp->choice) {
  case Seq_descr_user:
    uop = (UserObjectPtr) sdp->data.ptrvalue;
    if (uop == NULL) break;
    if (StringICmp (uop->_class, "SMART_V1.0") == 0) {
      vsp->is_smupd_in_sep = TRUE;
    }
    oip = uop->type;
    if (oip != NULL) {
      if (StringICmp (oip->str, "GenomeBuild") == 0) {
        vsp->is_gpipe_in_sep = TRUE;
      }
    }
    break;
  case Seq_descr_source:
    biop = (BioSourcePtr) sdp->data.ptrvalue;
    if (biop == NULL) break;
    if (biop->genome == GENOME_genomic) {
      vsp->bsp_genomic_in_sep = TRUE;
    }
    break;
  case Seq_descr_molinfo:
    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip == NULL) break;
    switch (mip->tech) {
    case MI_TECH_htgs_1:
    case MI_TECH_htgs_2:
    case MI_TECH_htgs_3:
    case MI_TECH_htgs_0:
      vsp->is_htg_in_sep = TRUE;
      break;
    case MI_TECH_barcode:
      vsp->is_barcode_sep = TRUE;
      break;
    default:
      break;
    }
    break;
  default:
    break;
  }
}

static void FindMultiIntervalGenes (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  BoolPtr    multiIntervalGenesP;
  SeqLocPtr  slp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_GENE) return;
  multiIntervalGenesP = (BoolPtr) userdata;
  if (multiIntervalGenesP == NULL) return;

  slp = sfp->location;
  if (slp == NULL) return;
  switch (slp->choice) {
    case SEQLOC_PACKED_INT :
    case SEQLOC_PACKED_PNT :
    case SEQLOC_MIX :
    case SEQLOC_EQUIV :
      *multiIntervalGenesP = TRUE;
      break;
    default :
      break;
  }
}

static void FindSegmentedBioseqs (
  BioseqPtr bsp,
  Pointer userdata
)

{
  BoolPtr  segmentedBioseqsP;

  if (bsp == NULL || bsp->repr != Seq_repr_seg) return;
  segmentedBioseqsP = (BoolPtr) userdata;
  if (segmentedBioseqsP == NULL) return;
  *segmentedBioseqsP = TRUE;
}
static void SetPubScratchData (SeqDescrPtr sdp, Pointer userdata)

{
  AuthListPtr    alp;
  Char           buf [2048];
  CitGenPtr      cgp;
  CharPtr        consortium, str, tmp;
  ValNodePtr     vnp;
  ObjValNodePtr  ovp;
  PubdescPtr     pdp;

  if (sdp == NULL || sdp->choice != Seq_descr_pub || sdp->extended == 0) return;
  ovp = (ObjValNodePtr) sdp;
  pdp = (PubdescPtr) sdp->data.ptrvalue;
  if (pdp == NULL) return;

  vnp = pdp->pub;

  /* skip over just serial number */

  if (vnp != NULL && vnp->choice == PUB_Gen && vnp->next != NULL) {
    cgp = (CitGenPtr) vnp->data.ptrvalue;
    if (cgp != NULL) {
      if (StringNICmp ("BackBone id_pub", cgp->cit, 15) != 0) {
        if (cgp->cit == NULL && cgp->journal == NULL && cgp->date == NULL && cgp->serial_number) {
          vnp = vnp->next;
        }
      }
    }
  }

  if (PubLabelUnique (vnp, buf, sizeof (buf) - 1, OM_LABEL_CONTENT, TRUE) > 0) {
    alp = GetAuthListPtr (pdp, NULL);
    if (alp != NULL) {
      consortium = NULL;
      str = GetAuthorsString (GENBANK_FMT, alp, &consortium, NULL, NULL);
      tmp = MemNew (StringLen (buf) + StringLen (str) + StringLen (consortium) + 10);
      if (tmp != NULL) {
        StringCpy (tmp, buf);
        if (StringDoesHaveText (str)) {
          StringCat (tmp, "; ");
          StringCat (tmp, str);
        }
        if (StringDoesHaveText (consortium)) {
          StringCat (tmp, "; ");
          StringCat (tmp, consortium);
        }
        ovp->idx.scratch = tmp;
      }
      MemFree (str);
      MemFree (consortium);
    }
  }
}

static void ClearPubScratchData (SeqDescrPtr sdp, Pointer userdata)

{
  ObjValNodePtr  ovp;

  if (sdp == NULL || sdp->choice != Seq_descr_pub || sdp->extended == 0) return;
  ovp = (ObjValNodePtr) sdp;
  ovp->idx.scratch = MemFree (ovp->idx.scratch);
}

static ValNodePtr SetUpValidateGeneticCodes (void)

{
  Char            ch;
  GeneticCodePtr  codes;
  GeneticCodePtr  gcp;
  ValNodePtr      gencodelist = NULL;
  Int2            i;
  Int4            id;
  Int2            j;
  Char            name [64];
  CharPtr         ptr;
  Char            str [256];
  ValNodePtr      tmp;

  codes = GeneticCodeTableLoad ();
  if (codes != NULL) {
    for (gcp = codes; gcp != NULL; gcp = gcp->next) {
      id = 0;
      str [0] = '\0';
      for (tmp = (ValNodePtr) gcp->data.ptrvalue; tmp != NULL; tmp = tmp->next) {
        switch (tmp->choice) {
          case 1 :
            if (StringLen (str) < 1) {
              StringNCpy_0 (str, (CharPtr) tmp->data.ptrvalue, sizeof (str));
              ptr = str;
              ch = *ptr;
              while (ch != '\0') {
                if (ch == '/') {
                  *ptr = '-';
                }
                ptr++;
                ch = *ptr;
              }
            }
            break;
          case 2 :
            id = tmp->data.intvalue;
            break;
          default :
            break;
        }
      }
      if (id != 7 && id != 8) {
        if (id > 0 /* && id < 30 */ ) {
          i = 0;
          if (StringLen (str + i) > 0) {
            ch = str [i];
            while (ch == ' ' || ch == ';') {
              i++;
              ch = str [i];
            }
            j = 0;
            ch = str [i + j];
            while (ch != '\0' && ch != ';') {
              name [j] = ch;
              j++;
              ch = str [i + j];
            }
            name [j] = '\0';
            i += j;
            if (ch == ';') {
              StringCat (name, ", etc.");
            }
            ValNodeCopyStr (&gencodelist, (Uint1) id, name);
          }
        }
      }
    }
  }
  return gencodelist;
}

typedef struct frd {
  ValidStructPtr    vsp;
  GatherContextPtr  gcp;
  /*
  CharPtr           string;
  */
} FindRepData, PNTR FindRepPtr;

static void FindRepValidate (Uint2 entityID, Uint4 itemID, Uint2 itemtype, Pointer userdata)

{
  FindRepPtr        frp;
  GatherContextPtr  gcp;
  ValidStructPtr    vsp;

  frp = (FindRepPtr) userdata;
  vsp = frp->vsp;
  gcp = frp->gcp;

  gcp->entityID = entityID;
  gcp->itemID = itemID;
  gcp->thistype = itemtype;

  ValidErr (vsp, SEV_ERROR, ERR_GENERIC_EmbeddedScript, "Script tag found in item");
}

static CharPtr findrepstrs [] = {
  "<script", "<object", "<applet", "<embed", "<form", "javascript:", "vbscript:", NULL
};

typedef struct vvmdata {
  Int2        num_mrnas;
  Boolean     accounted_for;
  Boolean     products_unique;
  Boolean     featid_matched;
  SeqFeatPtr  nearbygene;
  SeqFeatPtr  nearbycds;
  SeqFeatPtr  nearbymrna;
} VvmData, PNTR VvmDataPtr;

static void AddScratchToFeatures (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  sfp->idx.scratch = (Pointer) MemNew (sizeof (VvmData));
}

static void ClearScratchOnFeatures (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  sfp->idx.scratch = MemFree (sfp->idx.scratch);
}

static void SetupFeatureScratchData (
  BioseqPtr bsp,
  Pointer userdata
)

{
  SeqFeatPtr         currcds = NULL, currmrna = NULL, currgene = NULL;
  SeqMgrFeatContext  fcontext;
  SeqFeatPtr         sfp;
  VvmDataPtr         vdp;

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  while (sfp != NULL) {
    switch (sfp->idx.subtype) {
      case FEATDEF_GENE :
        currgene = sfp;
        break;
      case FEATDEF_CDS :
        currcds = sfp;
        vdp = (VvmDataPtr) sfp->idx.scratch;
        if (vdp != NULL) {
          if (vdp->nearbygene == NULL) {
            vdp->nearbygene = currgene;
          }
          if (vdp->nearbymrna == NULL) {
            vdp->nearbymrna = currmrna;
          }
        }
        if (currgene != NULL) {
          vdp = (VvmDataPtr) currgene->idx.scratch;
          if (vdp != NULL) {
            if (vdp->nearbycds == NULL) {
              vdp->nearbycds = currcds;
            }
          }
        }
        if (currmrna != NULL) {
          vdp = (VvmDataPtr) currmrna->idx.scratch;
          if (vdp != NULL) {
            if (vdp->nearbycds == NULL) {
              vdp->nearbycds = currcds;
            }
          }
        }
        break;
      case FEATDEF_mRNA :
        currmrna = sfp;
        vdp = (VvmDataPtr) sfp->idx.scratch;
        if (vdp != NULL) {
          if (vdp->nearbygene == NULL) {
            vdp->nearbygene = currgene;
          }
        }
        if (currgene != NULL) {
          vdp = (VvmDataPtr) currgene->idx.scratch;
          if (vdp != NULL) {
            if (vdp->nearbymrna == NULL) {
              vdp->nearbymrna = currmrna;
            }
          }
        }
        break;
      default :
        vdp = (VvmDataPtr) sfp->idx.scratch;
        if (vdp != NULL) {
          if (vdp->nearbygene == NULL) {
            vdp->nearbygene = currgene;
          }
          if (vdp->nearbymrna == NULL) {
            vdp->nearbymrna = currmrna;
          }
          if (vdp->nearbycds == NULL) {
            vdp->nearbycds = currcds;
          }
        }
        break;
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
  }
}

static void TestDeletedOrReplacedECnumbers (ValidStructPtr vsp)

{
  FileCache   fc;
  FILE        *fp = NULL;
  TextFsaPtr  fsa;
  Char        line [512];
  Char        path [PATH_MAX];
  CharPtr     ptr;
  ErrSev      sev;
  CharPtr     str;
  CharPtr     tmp;

  /* only check first time program runs validator */

  fsa = (TextFsaPtr) GetAppProperty ("ReplacedEECNumberFSA");
  if (fsa != NULL) return;

  GetSpecificECNumberFSA ();
  GetAmbiguousECNumberFSA ();
  GetDeletedECNumberFSA ();
  GetReplacedECNumberFSA ();

  if (FindPath ("ncbi", "ncbi", "data", path, sizeof (path))) {
    FileBuildPath (path, NULL, "ecnum_replaced.txt");
    sev = ErrSetMessageLevel (SEV_ERROR);
    fp = FileOpen (path, "r");
    ErrSetMessageLevel (sev);
    if (fp != NULL) {
      FileCacheSetup (&fc, fp);

      str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
      while (str != NULL) {
        if (StringDoesHaveText (str)) {
          ptr = StringChr (str, '\t');
          if (ptr != NULL) {
            *ptr = '\0';
            ptr++;
            if (! ECnumberNotInList (str)) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_EcNumberProblem, "Replaced EC number %s still in live list", str);
            }
            if (ECnumberWasDeleted (str)) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_EcNumberProblem, "Replaced EC number %s in deleted list", str);
            }
            while (StringDoesHaveText (ptr)) {
              tmp = StringChr (ptr, '\t');
              if (tmp != NULL) {
                *tmp = '\0';
                tmp++;
              }
              if (ECnumberNotInList (ptr)) {
                ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_EcNumberProblem, "Replacement EC number %s not in live list", ptr);
              }
              if (ECnumberWasDeleted (ptr)) {
                ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_EcNumberProblem, "Replacement EC number %s in deleted list", ptr);
              }
              ptr = tmp;
            }
          }
        }
        str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
      }

      FileClose (fp);
    }
  }

  if (FindPath ("ncbi", "ncbi", "data", path, sizeof (path))) {
    FileBuildPath (path, NULL, "ecnum_deleted.txt");
    sev = ErrSetMessageLevel (SEV_ERROR);
    fp = FileOpen (path, "r");
    ErrSetMessageLevel (sev);
    if (fp != NULL) {
      FileCacheSetup (&fc, fp);

      str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
      while (str != NULL) {
        if (StringDoesHaveText (str)) {
          ptr = StringChr (str, '\t');
          if (ptr != NULL) {
            *ptr = '\0';
            ptr++;
            if (! ECnumberNotInList (str)) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_EcNumberProblem, "Deleted EC number %s still in live list", str);
            }
          }
        }
        str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
      }

      FileClose (fp);
    }
  }
}


typedef struct collisioninfo {
  CharPtr str;
  SeqIdPtr sip;
  BioseqPtr bsp;
} CollisionInfoData, PNTR CollisionInfoPtr;


static CollisionInfoPtr CollisionInfoNew (SeqIdPtr sip, BioseqPtr bsp)
{
  CollisionInfoPtr cip = (CollisionInfoPtr) MemNew (sizeof (CollisionInfoData));
  cip->sip = sip;
  cip->bsp = bsp;
  cip->str = SeqIdWholeLabel (sip, PRINTID_FASTA_SHORT);
  return cip;
}


static CollisionInfoPtr CollisionInfoFree (CollisionInfoPtr cip)
{
  if (cip != NULL) {
    cip->str = MemFree (cip->str);
    cip = MemFree (cip);
  }
  return cip;
}


static void LongCollisionCallback (BioseqPtr bsp, Pointer data)
{
  SeqIdPtr sip;

  if (bsp == NULL || data == NULL) {
    return;
  }

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (!IsNCBIFileID(sip)) {
      ValNodeAddPointer ((ValNodePtr PNTR) data, 0, CollisionInfoNew (sip, bsp));
    }
  }
}


static int LIBCALLBACK SortVnpByCollisionInfo (VoidPtr ptr1, VoidPtr ptr2)

{
  CollisionInfoPtr  cip1, cip2;
  ValNodePtr    vnp1, vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  cip1 = (CollisionInfoPtr) vnp1->data.ptrvalue;
  cip2 = (CollisionInfoPtr) vnp2->data.ptrvalue;

  if (cip1 == NULL || cip2 == NULL) return 0;
  return StringCmp (cip1->str, cip2->str);
}


static void FindLongIdsThatCollideWhenTruncated (SeqEntryPtr sep, ValidStructPtr vsp, Int4 trunc_len)
{
  ValNodePtr id_list = NULL, vnp, vnp_c;
  CollisionInfoPtr  cip1, cip2;
  BioseqPtr         oldbsp;

  VisitBioseqsInSep (sep, &id_list, LongCollisionCallback);
  id_list = ValNodeSort (id_list, SortVnpByCollisionInfo);
  oldbsp = vsp->bsp;

  for (vnp = id_list; vnp != NULL; vnp = vnp->next) {
    cip1 = (CollisionInfoPtr) vnp->data.ptrvalue;
    vnp_c = vnp->next;
    while (vnp_c != NULL && (cip2 = (CollisionInfoPtr) vnp_c->data.ptrvalue) != NULL
          && StringNCmp (cip1->str, cip2->str, trunc_len) == 0) {
      vsp->bsp = cip2->bsp;
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_BadSeqIdFormat,
                "First %d characters of %s and %s are identical", trunc_len, cip1->str, cip2->str);
      vnp = vnp_c;
      vnp_c = vnp_c->next;
    }
  }

  vsp->bsp = oldbsp;
  vnp = id_list;
  while (vnp != NULL) {
    vnp_c = vnp->next;
    vnp->data.ptrvalue = CollisionInfoFree (vnp->data.ptrvalue);
    vnp->next = NULL;
    vnp = ValNodeFree (vnp);
    vnp = vnp_c;
  }
}


static void ValLookForBigFarSeqs (
  BioseqPtr bsp,
  Pointer userdata
)

{
  Int4         count = 0;
  DeltaSeqPtr  dsp;
  Boolean      is_ddbj = FALSE;
  SeqIdPtr     sip;
  BoolPtr      toomanyfarP;

  if (bsp == NULL || userdata == NULL) return;

  if (bsp->repr != Seq_repr_delta) return;
  if (bsp->seq_ext_type != 4) return;

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_DDBJ) {
      is_ddbj = TRUE;
    }
  }

  if (! is_ddbj) return;

  for (dsp = (DeltaSeqPtr) bsp->seq_ext; dsp != NULL; dsp = dsp->next) {
    if (dsp->choice == 1) {
      count++;
    }
  }

  if (count > 10000) {
    toomanyfarP = (BoolPtr) userdata;
    *toomanyfarP = TRUE;
  }
}

static Boolean ValTooManyFarComponents (
  SeqEntryPtr sep
)

{
  Boolean  toomanyfar = FALSE;

  if (sep == NULL) return FALSE;

  VisitBioseqsInSep (sep, (Pointer) &toomanyfar, ValLookForBigFarSeqs);

  return toomanyfar;
}

static CharPtr inferencePrefix [] = {
  "",
  "similar to sequence",
  "similar to AA sequence",
  "similar to DNA sequence",
  "similar to RNA sequence",
  "similar to RNA sequence, mRNA",
  "similar to RNA sequence, EST",
  "similar to RNA sequence, other RNA",
  "profile",
  "nucleotide motif",
  "protein motif",
  "ab initio prediction",
  "alignment",
  NULL
};


static CharPtr NextColonOrVerticalBarPtr (CharPtr ptr)

{
  Char  ch = '\0';

  if (ptr == NULL) return NULL;

  ch = *ptr;
  while (ch != '\0') {
    if (ch == ':' || ch == '|') return ptr;
    ptr++;
    ch = *ptr;
  }

  return NULL;
}

typedef struct valcountdata {
  Int4  numInferences;
  Int4  numAccessions;
} ValCountData, PNTR ValCountPtr;

static void ValCountInfAccnVer (SeqFeatPtr sfp, Pointer userdata)

{
  Int2         best, j;
  Char         ch;
  GBQualPtr    gbq;
  size_t       len;
  CharPtr      nxt;
  CharPtr      ptr;
  CharPtr      rest;
  CharPtr      str;
  CharPtr      tmp;
  ValCountPtr  vcp;


  if (sfp == NULL || userdata == NULL) return;
  vcp = (ValCountPtr) userdata;

  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    if (StringICmp (gbq->qual, "inference") != 0) continue;
    if (StringHasNoText (gbq->val)) continue;

    (vcp->numInferences)++;

    rest = NULL;
    best = -1;
    for (j = 0; inferencePrefix [j] != NULL; j++) {
      len = StringLen (inferencePrefix [j]);
      if (StringNICmp (gbq->val, inferencePrefix [j], len) != 0) continue;
      rest = gbq->val + len;
      best = j;
    }
    if (best < 0 || inferencePrefix [best] == NULL) continue;
    if (rest == NULL) continue;

    ch = *rest;
    while (IS_WHITESP (ch)) {
      rest++;
      ch = *rest;
    }
    if (StringNICmp (rest, "(same species)", 14) == 0) {
      rest += 14;
    }
    ch = *rest;
    while (IS_WHITESP (ch) || ch == ':') {
      rest++;
      ch = *rest;
    }
    if (StringHasNoText (rest)) continue;

    str = StringSave (rest);

    ptr = str;
    if (best == 12) {
      ptr = StringRChr (str, ':');
      if (ptr != NULL) {
        *ptr = '\0';
        ptr++;
      }
    }
    while (ptr != NULL) {
      nxt = StringChr (ptr, ',');
      if (nxt != NULL) {
        *nxt = '\0';
        nxt++;
      }
      tmp = NextColonOrVerticalBarPtr (ptr);
      if (tmp != NULL) {
        *tmp = '\0';
        tmp++;
        TrimSpacesAroundString (ptr);
        TrimSpacesAroundString (tmp);
        if (StringDoesHaveText (tmp)) {
          if (StringICmp (ptr, "INSD") == 0 || StringICmp (ptr, "RefSeq") == 0) {
            (vcp->numAccessions)++;
          }
        }
      }
      ptr = nxt;
    }

    MemFree (str);
  }
}

NLM_EXTERN Boolean TooManyInferenceAccessions (
  SeqEntryPtr sep,
  Int4Ptr numInferences,
  Int4Ptr numAccessions
)

{
  ValCountData  vcd;

  if (numInferences != NULL) {
    *numInferences = 0;
  }
  if (numAccessions != NULL) {
    *numAccessions = 0;
  }
  if (sep == NULL) return FALSE;

  vcd.numInferences = 0;
  vcd.numAccessions = 0;

  VisitFeaturesInSep (sep, (Pointer) &vcd, ValCountInfAccnVer);

  if (numInferences != NULL) {
    *numInferences = vcd.numInferences;
  }
  if (numAccessions != NULL) {
    *numAccessions = vcd.numAccessions;
  }

  if (vcd.numInferences > 1000 || vcd.numAccessions > 1000) return TRUE;

  return FALSE;
}

NLM_EXTERN Boolean ValidateSeqEntry (SeqEntryPtr sep, ValidStructPtr vsp)

{
  AuthListPtr     alp;
  AuthorPtr       ap;
  DatePtr         cd, dp;
  ContactInfoPtr  cip;
  CitSubPtr       csp;
  Uint2           entityID = 0;
  GatherScope     gs;
  BioseqSetPtr    bssp;
  SeqSubmitPtr    ssp = NULL;
  Boolean         do_many = FALSE;
  Boolean         mult_subs = FALSE;
  Boolean         farFetchProd;
  Boolean         first = TRUE;
  Int4            errors[6];
  Int2            i;
  Boolean         inferenceAccnCheck;
  Boolean         suppress_no_pubs = TRUE;
  Boolean         suppress_no_biosrc = TRUE;
  FeatProb        featprob;
  GatherContextPtr gcp = NULL;
  GatherContext   gc;
  SeqEntryPtr     fsep;
  BioseqPtr       fbsp = NULL;
  Int2            limit;
  SeqEntryPtr     oldsep;
  ErrSev          oldsev;
  ObjMgrDataPtr   omdp;
  SeqEntryPtr     topsep = NULL;
  SeqEntryPtr     tmp;
  ValNodePtr      bsplist;
  SubmitBlockPtr  sbp;
  ErrSev          sev;
  SeqIdPtr        sip;
  Boolean         has_multi_int_genes = FALSE;
  Boolean         has_seg_bioseqs = FALSE;
  Boolean         isGPS = FALSE;
  Boolean         isPatent = FALSE;
  Boolean         isPDB = FALSE;
  FindRepData     frd;
  Int4            numInferences;
  Int4            numAccessions;

  if (sep == NULL || vsp == NULL) return FALSE;

  genetic_code_name_list = SetUpValidateGeneticCodes ();

  vsp->useSeqMgrIndexes = TRUE; /* now always use indexing */

  for (i = 0; i < 6; i++)       /* keep errors between clears */
    errors[i] = 0;

  MemSet ((Pointer) &featprob, 0, sizeof (FeatProb));

  if (vsp->useSeqMgrIndexes) {
    entityID = ObjMgrGetEntityIDForChoice (sep);

    if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
      oldsev = ErrSetMessageLevel (SEV_MAX);
      SeqMgrIndexFeatures (entityID, NULL);
      ErrSetMessageLevel (oldsev);
    }
    SeqMgrExploreBioseqs (entityID, NULL, (Pointer) &featprob, CountMisplacedFeatures, TRUE, TRUE, TRUE);

    topsep = GetTopSeqEntryForEntityID (entityID);
    VisitGraphsInSep (topsep, (Pointer) &featprob, CheckGraphPacking);
    VisitFeaturesInSep (topsep, (Pointer) &featprob, CountGeneXrefs);
    VisitFeaturesInSep (topsep, (Pointer) &featprob, CountFeatLocIdTypes);
    VisitBioseqsInSep (topsep, (Pointer) &featprob, CheckTpaHist);
  } else {

    /* if not using indexing, still need feature->idx.subtype now */

    entityID = ObjMgrGetEntityIDForChoice (sep);
    AssignIDsInEntity (entityID, 0, NULL);
  }

  /* Seq-submit can have multiple entries with no Bioseq-set wrapper */

  omdp = ObjMgrGetData (entityID);
  if (omdp != NULL && omdp->datatype == OBJ_SEQSUB) {
    ssp = (SeqSubmitPtr) omdp->dataptr;
    if (ssp != NULL && ssp->data != NULL) {
      if (sep->next != NULL) {
        do_many = TRUE;
        mult_subs = TRUE;
      }
    }
  }

  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) (sep->data.ptrvalue);
    switch (bssp->_class) {
    /* case BioseqseqSet_class_genbank: */
    case BioseqseqSet_class_pir:
    case BioseqseqSet_class_gibb:
    case BioseqseqSet_class_gi:
    case BioseqseqSet_class_swissprot:
      sep = bssp->seq_set;
      do_many = TRUE;
      break;
    case BioseqseqSet_class_gen_prod_set:
      isGPS = TRUE;
    default:
      break;
    }
  }

  /* if no pubs or biosource, only one message, not one per bioseq */

  if (mult_subs) {
    for (tmp = sep; tmp != NULL; tmp = tmp->next) {
      LookForAnyPubAndOrg (tmp, &suppress_no_pubs, &suppress_no_biosrc);
    }
  } else {
    LookForAnyPubAndOrg (sep, &suppress_no_pubs, &suppress_no_biosrc);
  }

  if (GetAppProperty ("ValidateExons") != NULL) {
    vsp->validateExons = TRUE;
  }

  vsp->is_htg_in_sep = FALSE;
  vsp->is_barcode_sep = FALSE;
  vsp->is_refseq_in_sep = FALSE;
  vsp->is_gpipe_in_sep = FALSE;
  vsp->is_gps_in_sep = FALSE;
  vsp->other_sets_in_sep = FALSE;
  vsp->is_embl_ddbj_in_sep = FALSE;
  vsp->is_old_gb_in_sep = FALSE;
  vsp->is_insd_in_sep = FALSE;
  vsp->only_lcl_gnl_in_sep = FALSE;
  vsp->has_gnl_prot_sep = FALSE;
  vsp->bsp_genomic_in_sep = FALSE;
  vsp->is_smupd_in_sep = FALSE;

  VisitBioseqsInSep (sep, (Pointer) vsp, LookForBioseqFields);
  VisitSetsInSep (sep, (Pointer) vsp, LookForBioseqSetFields);
  VisitDescriptorsInSep (sep, (Pointer) vsp, LookForSeqDescrFields);

  VisitFeaturesInSep (sep, (Pointer) &has_multi_int_genes, FindMultiIntervalGenes);
  vsp->has_multi_int_genes = has_multi_int_genes;
  VisitBioseqsInSep (sep, (Pointer) &has_seg_bioseqs, FindSegmentedBioseqs);
  vsp->has_seg_bioseqs = has_seg_bioseqs;

  /*
  vsp->is_htg_in_sep = FALSE;
  VisitDescriptorsInSep (sep, (Pointer) &(vsp->is_htg_in_sep), LookForHTG);
  vsp->is_barcode_sep = FALSE;
  VisitDescriptorsInSep (sep, (Pointer) &(vsp->is_barcode_sep), LookForBarcode);
  vsp->is_smupd_in_sep = FALSE;
  VisitDescriptorsInSep (sep, (Pointer) &(vsp->is_smupd_in_sep), LookForSMUPD);
  vsp->is_gps_in_sep = FALSE;
  SeqEntryExplore (sep, (Pointer) &(vsp->is_gps_in_sep), LookForGPS);
  vsp->other_sets_in_sep = FALSE;
  SeqEntryExplore (sep, (Pointer) &(vsp->other_sets_in_sep), LookForNonGPS);
  vsp->is_refseq_in_sep = FALSE;
  VisitBioseqsInSep (sep, (Pointer) &(vsp->is_refseq_in_sep), LookForNC);
  vsp->is_embl_ddbj_in_sep = FALSE;
  VisitBioseqsInSep (sep, (Pointer) &(vsp->is_embl_ddbj_in_sep), LookForEmblDdbj);
  vsp->is_insd_in_sep = FALSE;
  VisitBioseqsInSep (sep, (Pointer) &(vsp->is_insd_in_sep), LookForGEDseqID);
  vsp->only_lcl_gnl_in_sep = FALSE;
  VisitBioseqsInSep (sep, (Pointer) &(vsp->only_lcl_gnl_in_sep), LookForLclGnl);
  vsp->has_gnl_prot_sep = FALSE;
  VisitBioseqsInSep (sep, (Pointer) &(vsp->has_gnl_prot_sep), LookForProteinGnl);
  */

  vsp->feat_loc_has_gi = featprob.loc_has_gi;
  vsp->feat_prod_has_gi = featprob.prod_has_gi;

  globalvsp = vsp;              /* for spell checker */

  inferenceAccnCheck = vsp->inferenceAccnCheck;

  while (sep != NULL) {
    vsp->far_fetch_failure = FALSE;

    /* calculate strings for LookForMultipleUnpubPubs test only once for genome product set efficiency */
    VisitDescriptorsInSep (sep, NULL, SetPubScratchData);

    MemSet (&gs, 0, sizeof (GatherScope));
    gs.scope = sep;             /* default is to scope to this set */

    ValidStructClear (vsp);
    vsp->sep = sep;

    MemSet ((Pointer) &gc, 0, sizeof (GatherContext));
    gcp = &gc;
    gc.entityID = ObjMgrGetEntityIDForChoice (sep);
    gc.itemID = 1;
    if (IS_Bioseq (sep)) {
      gc.thistype = OBJ_BIOSEQ;
    } else {
      gc.thistype = OBJ_BIOSEQSET;
    }
    vsp->gcp = gcp;             /* above needed for ValidErr */
    vsp->suppress_no_pubs = suppress_no_pubs;
    vsp->suppress_no_biosrc = suppress_no_biosrc;

    if (vsp->is_refseq_in_sep && vsp->is_insd_in_sep) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_PKG_INSDRefSeqPackaging,
                "INSD and RefSeq records should not be present in the same set");
    }

    if (vsp->is_gps_in_sep && vsp->other_sets_in_sep) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_PKG_GPSnonGPSPackaging,
                "Genomic product set and mut/pop/phy/eco set records should not be present in the same set");
    }

    /* build seqmgr feature indices if not already done */

    bsplist = NULL;
    if (vsp->useSeqMgrIndexes) {
      entityID = ObjMgrGetEntityIDForChoice (sep);

      if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
        oldsev = ErrSetMessageLevel (SEV_MAX);
        SeqMgrIndexFeatures (entityID, NULL);
        ErrSetMessageLevel (oldsev);
      }

      /* lock all remote genome components, locations, and products in advance */

      limit = vsp->validationLimit;
      if (! ValTooManyFarComponents (sep)) {
        if (limit == VALIDATE_ALL || limit == VALIDATE_INST || limit == VALIDATE_HIST) {
          farFetchProd = (Boolean) (vsp->farFetchCDSproducts || vsp->farFetchMRNAproducts);
          oldsev = ErrSetMessageLevel (SEV_WARNING);
          bsplist = LockFarComponentsEx (sep, TRUE, TRUE, farFetchProd, NULL);
          ErrSetMessageLevel (oldsev);
        }
      }
    }

    fsep = FindNthBioseq (sep, 1);
    fbsp = NULL;
    if (fsep != NULL && IS_Bioseq (fsep)) {
      fbsp = (BioseqPtr) fsep->data.ptrvalue;
      /* report context as first bioseq */
      vsp->bsp = fbsp;
    }

    if (fbsp == NULL) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_PKG_NoBioseqFound, "No Bioseqs in this entire record.");
    } else {

      for (sip = fbsp->id; sip != NULL; sip = sip->next) {
        if (sip->choice == SEQID_PATENT) {
          isPatent = TRUE;
        } else if (sip->choice == SEQID_PDB) {
          isPDB = TRUE;
        }
      }

      if (first) {
        TestDeletedOrReplacedECnumbers (vsp);

        if (suppress_no_pubs && (! vsp->seqSubmitParent)) {
          omdp = ObjMgrGetData (gc.entityID);
          if (omdp == NULL || omdp->datatype != OBJ_SEQSUB) {
            sev = SEV_ERROR;
            if ((!isGPS) && (!IsNoncuratedRefSeq (fbsp, &sev)) && (! IsGpipe (fbsp))) {
              ValidErr (vsp, sev, ERR_SEQ_DESCR_NoPubFound, "No publications anywhere on this entire record.");
            }
          }
        }
        if (suppress_no_biosrc) {
          if ((!isPatent) && ((!isPDB))) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_NoOrgFound, "No organism name anywhere on this entire record.");
          }
        }

        if (featprob.num_misplaced_features > 1) {
          ValidErr (vsp, SEV_REJECT, ERR_SEQ_PKG_FeaturePackagingProblem, "There are %d mispackaged features in this record.", (int) featprob.num_misplaced_features);
        } else if (featprob.num_misplaced_features == 1) {
          ValidErr (vsp, SEV_REJECT, ERR_SEQ_PKG_FeaturePackagingProblem, "There is %d mispackaged feature in this record.", (int) featprob.num_misplaced_features);
        }
        if (featprob.num_small_genome_set_misplaced > 1) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_FeaturePackagingProblem, "There are %d mispackaged features in this small genome set record.", (int) featprob.num_small_genome_set_misplaced);
        } else if (featprob.num_small_genome_set_misplaced == 1) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_FeaturePackagingProblem, "There is %d mispackaged feature in this small genome set record.", (int) featprob.num_small_genome_set_misplaced);
        }

        if (featprob.num_misplaced_graphs > 1) {
          ValidErr (vsp, SEV_REJECT, ERR_SEQ_PKG_GraphPackagingProblem, "There are %d mispackaged graphs in this record.", (int) featprob.num_misplaced_graphs);
        } else if (featprob.num_misplaced_graphs == 1) {
          ValidErr (vsp, SEV_REJECT, ERR_SEQ_PKG_GraphPackagingProblem, "There is %d mispackaged graph in this record.", (int) featprob.num_misplaced_graphs);
        }

        /*
        if (featprob.num_archaic_locations > 1) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_ArchaicFeatureLocation, "There are %d archaic feature locations in this record.", (int) featprob.num_archaic_locations);
        } else if (featprob.num_archaic_locations == 1) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_ArchaicFeatureLocation, "There is %d archaic feature location in this record.", (int) featprob.num_archaic_locations);
        }

        if (featprob.num_archaic_products > 1) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_ArchaicFeatureProduct, "There are %d archaic feature products in this record.", (int) featprob.num_archaic_products);
        } else if (featprob.num_archaic_products == 1) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_ArchaicFeatureProduct, "There is %d archaic feature product in this record.", (int) featprob.num_archaic_products);
        }
        */

        if (featprob.num_gene_feats == 0 && featprob.num_gene_xrefs > 0) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_OnlyGeneXrefs, "There are %ld gene xrefs and no gene features in this record.", (long) featprob.num_gene_xrefs);
        }

        if (featprob.num_tpa_with_hist > 0 && featprob.num_tpa_without_hist > 0) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_TpaAssmeblyProblem, "There are %ld TPAs with history and %ld without history in this record.",
                    (long) featprob.num_tpa_with_hist, (long) featprob.num_tpa_without_hist);
        }

        if (featprob.has_gi && featprob.num_tpa_without_hist > 0) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_TpaAssmeblyProblem, "There are %ld TPAs without history in this record, but the record has a gi number assignment.",
                    (long) featprob.num_tpa_without_hist);
        }

        if (vsp->indexerVersion && vsp->has_gnl_prot_sep && (! vsp->is_refseq_in_sep)) {
          if (FindNucBioseq (sep) != NULL) {
            ValidErr (vsp, SEV_INFO, ERR_SEQ_INST_ProteinsHaveGeneralID, "INDEXER_ONLY - Protein bioseqs have general seq-id.");
          }
        }

        first = FALSE;
      }

      vsp->bsp = NULL;

      topsep = GetTopSeqEntryForEntityID (gc.entityID);
      oldsep = SeqEntrySetScope (topsep);

      /* disabled for now
      FindLongIdsThatCollideWhenTruncated (topsep, vsp, 30);
      */

      /* do validator tests using Discrepancy Report */
      ValidateGeneLocusTags (topsep, vsp);
      ValidateShortIntrons (topsep, vsp);

      VisitFeaturesInSep (sep, NULL, AddScratchToFeatures);
      VisitBioseqsInSep (sep, NULL, SetupFeatureScratchData);

      /* AssignIDsInEntity (gc.entityID, 0, NULL); */

      if (inferenceAccnCheck) {
        numInferences = 0;
        numAccessions = 0;
        if (TooManyInferenceAccessions (sep, &numInferences, &numAccessions)) {
          ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_TooManyInferenceAccessions,
                    "Skipping validation of %ld /inference qualifiers with %ld accessions",
                    (long) numInferences, (long) numAccessions);

          /* suppress inference accession.version check for this record */
          vsp->inferenceAccnCheck = FALSE;
        }
      }

      GatherSeqEntry (sep, (Pointer) vsp, Valid1GatherProc, &gs);

      /* restore inferenceAccnCheck flag for next record */
      vsp->inferenceAccnCheck = inferenceAccnCheck;

      if (ssp != NULL) {
        if (ssp->datatype == 1) {
          vsp->bsp = NULL;
          vsp->bssp = NULL;
          vsp->sfp = NULL;
          vsp->descr = NULL;
          vsp->gcp = NULL;
          sbp = ssp->sub;
          if (sbp != NULL) {
            csp = sbp->cit;
            if (csp != NULL) {
              alp = csp->authors;
              if (alp != NULL) {
                ValidateAffil (vsp, alp->affil);
              }
              ValidateCitSub (vsp, csp);
            }
            cip = sbp->contact;
            if (cip != NULL) {
              ap = cip->contact;
              if (ap != NULL) {
                ValidateAffil (vsp, ap->affil);
              }
            }
            if (sbp->hup) {
              dp = sbp->reldate;
              cd = DateCurr ();
              if (dp != NULL && cd != NULL) {
                if (DateMatch (dp, cd, FALSE) == -1) {
                  ValidErr (vsp, SEV_WARNING, ERR_GENERIC_PastReleaseDate,
                            "Record release date has already passed");
                }
              }
              DateFree (cd);
            }
          }
        }
      }

      vsp->gcp = NULL;
      ValidateFeatCits (sep, vsp);
      vsp->gcp = NULL;

      vsp->gcp = NULL;
      ValidateFeatIDs (sep, gc.entityID, vsp);
      vsp->gcp = NULL;

      vsp->gcp = NULL;
      ValidateSeqIdCase (sep, vsp);
      vsp->gcp = NULL;

      if (vsp->validateAlignments) {
        vsp->gcp = NULL;
        ValidateSeqAlignWithinValidator (vsp, sep, vsp->alignFindRemoteBsp, vsp->doSeqHistAssembly);
        vsp->gcp = NULL;
      }

      if (vsp->far_fetch_failure) {
        vsp->gcp = NULL;
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_FarFetchFailure, "Far fetch failures caused some validator tests to be bypassed");
      }

      VisitFeaturesInSep (sep, NULL, ClearScratchOnFeatures);

      SeqEntrySetScope (oldsep);

      VisitDescriptorsInSep (sep, NULL, ClearPubScratchData);
    }

    if (vsp->useSeqMgrIndexes) {

      /* unlock all pre-locked remote genome components */

      bsplist = UnlockFarComponents (bsplist);
    }

    if (do_many) {
      for (i = 0; i < 6; i++)
        errors[i] += vsp->errors[i];
      sep = sep->next;
    } else
      sep = NULL;
  }

  MemSet ((Pointer) &gc, 0, sizeof (GatherContext));
  gcp = &gc;
  gc.entityID = ObjMgrGetEntityIDForChoice (sep);
  vsp->gcp = gcp;
  frd.vsp = vsp;
  frd.gcp = gcp;

  limit = vsp->validationLimit;
  if (limit == VALIDATE_ALL) {
    /*
    frd.string = "?";
    */
    FindStringsInEntity (entityID, findrepstrs, FALSE, FALSE, FALSE, UPDATE_NEVER,
                         NULL, NULL, NULL, TRUE, FindRepValidate, (Pointer) &frd);
  }

  if (do_many) {
    for (i = 0; i < 6; i++)
      vsp->errors[i] = errors[i];
  }

  genetic_code_name_list = ValNodeFreeData (genetic_code_name_list);

  return TRUE;
}


static void ValidateSetContents (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
  BioseqPtr       bsp;
  ValidStructPtr  vsp;

  vsp = (ValidStructPtr) data;

  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) (sep->data.ptrvalue);
    if (ISA_aa (bsp->mol))
      vsp->protcnt++;
    else
      vsp->nuccnt++;
    if (bsp->repr == Seq_repr_seg)
      vsp->segcnt++;

  }
  return;
}


static CharPtr GetBioseqSetClass (Uint1 cl)
{
  if (cl == BioseqseqSet_class_nuc_prot)
    return ("nuc-prot");
  if (cl == BioseqseqSet_class_segset)
    return ("segset");
  if (cl == BioseqseqSet_class_conset)
    return ("conset");
  if (cl == BioseqseqSet_class_parts)
    return ("parts");
  if (cl == BioseqseqSet_class_gibb)
    return ("gibb");
  if (cl == BioseqseqSet_class_gi)
    return ("gi");
  if (cl == BioseqseqSet_class_genbank)
    return ("genbank");
  if (cl == BioseqseqSet_class_pir)
    return ("pir");
  if (cl == BioseqseqSet_class_pub_set)
    return ("pub-set");
  if (cl == BioseqseqSet_class_equiv)
    return ("equiv");
  if (cl == BioseqseqSet_class_swissprot)
    return ("swissprot");
  if (cl == BioseqseqSet_class_pdb_entry)
    return ("pdb-entry");
  if (cl == BioseqseqSet_class_mut_set)
    return ("mut-set");
  if (cl == BioseqseqSet_class_pop_set)
    return ("pop-set");
  if (cl == BioseqseqSet_class_phy_set)
    return ("phy-set");
  if (cl == BioseqseqSet_class_eco_set)
    return ("eco-set");
  if (cl == BioseqseqSet_class_gen_prod_set)
    return ("gen-prod-set");
  if (cl == BioseqseqSet_class_wgs_set)
    return ("wgs-set");
  if (cl == BioseqseqSet_class_small_genome_set)
    return ("small-genome-set");
  if (cl == BioseqseqSet_class_other)
    return ("other");
  return ("not-set");
}


static BioseqSetPtr FindGenProdSetParentOfBioseqSet (BioseqSetPtr bssp)
{
  if (bssp == NULL) {
    return NULL;
  } else if (bssp->idx.parenttype != OBJ_BIOSEQSET) {
    return NULL;
  } else if ((bssp = (BioseqSetPtr)bssp->idx.parentptr) == NULL) {
    return NULL;
  } else if (bssp->_class == BioseqseqSet_class_gen_prod_set) {
    return bssp;
  } else {
    return FindGenProdSetParentOfBioseqSet (bssp);
  }
}


static BioseqSetPtr FindGenProdSetParentOfBioseq (BioseqPtr bsp)
{
  BioseqSetPtr bssp;
  if (bsp == NULL) {
    return NULL;
  } else if (bsp->idx.parenttype != OBJ_BIOSEQSET) {
    return NULL;
  } else if ((bssp = (BioseqSetPtr)bsp->idx.parentptr) == NULL) {
    return NULL;
  } else if (bssp->_class == BioseqseqSet_class_gen_prod_set) {
    return bssp;
  } else {
    return FindGenProdSetParentOfBioseqSet (bssp);
  }
}


static void IfInGPSmustBeMrnaProduct (ValidStructPtr vsp, BioseqPtr bsp)

{
  /* see if in genomic product */
  if (FindGenProdSetParentOfBioseq(bsp) != NULL) {
    if (SeqMgrGetRNAgivenProduct (bsp, NULL) == NULL) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_GenomicProductPackagingProblem, "Nucleotide bioseq should be product of mRNA feature on contig, but is not");
    }
  }
}

static void IfInGPSmustBeCDSProduct (ValidStructPtr vsp, BioseqPtr bsp)

{
  BioseqSetPtr  bssp;
  BioseqPtr     contig;
  ValNodePtr    head, vnp;
  SeqEntryPtr   sep;
  SeqFeatPtr    sfp;

  /* see if in genomic product */
  if ((bssp = FindGenProdSetParentOfBioseq(bsp)) != NULL) {
    sep = bssp->seq_set;
    if (sep == NULL) return;
    if (! IS_Bioseq (sep)) return;
    contig = (BioseqPtr) sep->data.ptrvalue;
    if (contig == NULL) return;
    head = SeqMgrGetSfpProductList (bsp);
    for (vnp = head; vnp != NULL; vnp = vnp->next) {
      sfp = (SeqFeatPtr) vnp->data.ptrvalue;
      if (sfp == NULL) continue;
      if (BioseqFindFromSeqLoc (sfp->location) == contig) return;
    }
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_GenomicProductPackagingProblem, "Protein bioseq should be product of CDS feature on contig, but is not");
  }
}

NLM_EXTERN ValNodePtr BioseqGetSeqDescr(BioseqPtr bsp, Int2 type, ValNodePtr curr);


static void ValidateNucProtSet (BioseqSetPtr bssp, ValidStructPtr vsp)

{
  SeqDescrPtr   sdp;
  SeqEntryPtr   sep;
  BioSourcePtr  biop;
  BioseqPtr     bsp;
  BioseqSetPtr  bssp1;
  OrgRefPtr     orp;
  Int4          prot_biosource = 0;

  if (bssp->_class != BioseqseqSet_class_nuc_prot)
    return;

  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    if (IS_Bioseq (sep)) {
      bsp = (BioseqPtr) sep->data.ptrvalue;
      if (bsp == NULL) continue;
      if (ISA_na (bsp->mol)) {
        IfInGPSmustBeMrnaProduct (vsp, bsp);
      } else if (ISA_aa (bsp->mol)) {
        IfInGPSmustBeCDSProduct (vsp, bsp);
        sdp = BioseqGetSeqDescr (bsp, Seq_descr_source, NULL);
        if (sdp != NULL) {
          prot_biosource++;
        }
      }
    }

    if (!IS_Bioseq_set (sep))
      continue;

    bssp1 = sep->data.ptrvalue;
    if (bssp1 == NULL)
      continue;

    if (bssp1->_class != BioseqseqSet_class_segset) {
      ValidErr (vsp, SEV_REJECT, ERR_SEQ_PKG_NucProtNotSegSet,
                "Nuc-prot Bioseq-set contains wrong Bioseq-set, its class is \"%s\".", GetBioseqSetClass (bssp1->_class));
      break;
    }
  }

  if (prot_biosource > 1) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceOnProtein,
              "Nuc-prot set has %ld proteins with a BioSource descriptor", (long) prot_biosource);
  } else if (prot_biosource > 0) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceOnProtein,
              "Nuc-prot set has %ld protein with a BioSource descriptor", (long) prot_biosource);
  }

  for (sdp = bssp->descr; sdp != NULL; sdp = sdp->next) {
    if (sdp->choice == Seq_descr_title) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_NucProtSetHasTitle,
                "Nuc-prot set should not have title descriptor");
    }
  }

  for (sdp = bssp->descr; sdp != NULL; sdp = sdp->next) {
    if (sdp->choice == Seq_descr_source) {
      biop = (BioSourcePtr) sdp->data.ptrvalue;
      if (biop != NULL) {
        orp = biop->org;
        if (orp != NULL && StringDoesHaveText (orp->taxname)) return;
      }
    }
  }

  sep = vsp->sep;
  if (sep != NULL && IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && bssp->_class == BioseqseqSet_class_genbank) {
      sep = bssp->seq_set;
      if (sep != NULL && IS_Bioseq_set (sep)) {
        bssp = (BioseqSetPtr) sep->data.ptrvalue;
      }
    }
    if (bssp != NULL && bssp->_class == BioseqseqSet_class_gen_prod_set) {
      return;
    }
  }

  ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceMissing,
            "Nuc-prot set does not contain expected BioSource descriptor");
}

typedef struct incons {
  Boolean     diffs;
  MolInfoPtr  mip;
} Incons, PNTR InconsPtr;

static void FindInconsistMolInfos (SeqDescrPtr sdp, Pointer userdata)

{
  InconsPtr   icp;
  MolInfoPtr  mip;

  if (sdp == NULL || sdp->choice != Seq_descr_molinfo) return;
  icp = (InconsPtr) userdata;
  mip = (MolInfoPtr) sdp->data.ptrvalue;
  if (icp == NULL || mip == NULL) return;
  if (icp->mip == NULL) {
    icp->mip = mip;
  } else {
    if (icp->mip->biomol != mip->biomol) {
      icp->diffs = TRUE;
    }
  }
}

static void ValidateSegmentedSet (BioseqSetPtr bssp, ValidStructPtr vsp)

{
  SeqEntryPtr     sep;
  BioseqSetPtr    bssp1;
  BioseqPtr       bsp;
  Incons          inc;
  Uint1           mol = 0;

  if (bssp->_class != BioseqseqSet_class_segset)
    return;

  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    if (IS_Bioseq (sep)) {
      bsp = (BioseqPtr) sep->data.ptrvalue;
      if (bsp != NULL) {
        if (mol == 0 || mol == Seq_mol_other) {
          mol = bsp->mol;
        } else if (bsp->mol != Seq_mol_other) {
          if (ISA_na (bsp->mol) != ISA_na (mol)) {
            ValidErr (vsp, SEV_REJECT, ERR_SEQ_PKG_SegSetMixedBioseqs, "Segmented set contains mixture of nucleotides and proteins");
          }
        }
      }
    }

    if (!IS_Bioseq_set (sep))
      continue;

    bssp1 = sep->data.ptrvalue;
    if (bssp1 == NULL)
      continue;

    if (bssp1->_class != BioseqseqSet_class_parts) {
      ValidErr (vsp, SEV_REJECT, ERR_SEQ_PKG_SegSetNotParts,
                "Segmented set contains wrong Bioseq-set, its class is \"%s\".", GetBioseqSetClass (bssp1->_class));
      break;
    }
  }

  inc.diffs = FALSE;
  inc.mip = NULL;
  VisitDescriptorsInSet (bssp, (Pointer) &inc, FindInconsistMolInfos);
  if (inc.diffs) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_PKG_InconsistentMolInfoBiomols, "Segmented set contains inconsistent MolInfo biomols");
  }
}

static void ValidatePartsSet (BioseqSetPtr bssp, ValidStructPtr vsp)

{
  SeqEntryPtr     sep;
  BioseqSetPtr    bssp1;
  BioseqPtr       bsp;
  Uint1           mol = 0;

  if (bssp->_class != BioseqseqSet_class_parts)
    return;

  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    if (IS_Bioseq (sep)) {
      bsp = (BioseqPtr) sep->data.ptrvalue;
      if (bsp != NULL) {
        if (mol == 0 || mol == Seq_mol_other) {
          mol = bsp->mol;
        } else if (bsp->mol != Seq_mol_other) {
          if (ISA_na (bsp->mol) != ISA_na (mol)) {
            ValidErr (vsp, SEV_REJECT, ERR_SEQ_PKG_PartsSetMixedBioseqs, "Parts set contains mixture of nucleotides and proteins");
            break;
          }
        }
      }
    }
  }

  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    if (IS_Bioseq_set (sep)) {
      bssp1 = sep->data.ptrvalue;
      if (bssp1 == NULL)
        continue;

      ValidErr (vsp, SEV_REJECT, ERR_SEQ_PKG_PartsSetHasSets,
                "Parts set contains unwanted Bioseq-set, its class is \"%s\".", GetBioseqSetClass (bssp1->_class));
      break;
    }
  }
}

static Boolean CheckForInconsistentBiosources (SeqEntryPtr sep, ValidStructPtr vsp, OrgRefPtr PNTR orpp, BioseqSetPtr top)

{
  BioseqPtr       bsp;
  BioseqSetPtr    bssp;
  SeqEntryPtr     tmp;
  ValNodePtr      sdp;
  SeqFeatPtr      sfp;
  SeqMgrDescContext dcontext;
  SeqMgrFeatContext fcontext;
  BioSourcePtr    biop;
  OrgRefPtr       orp;
  OrgRefPtr       firstorp;
  GatherContextPtr gcp;
  Uint2           entityID = 0, oldEntityID;
  Uint4           itemID = 0, oldItemID;
  Uint2           itemtype = 0, oldItemtype;
  size_t          len, len1, len2;
  ErrSev          sev;
  CharPtr         sp;

  if (sep == NULL || vsp == NULL || orpp == NULL)
    return FALSE;
  gcp = vsp->gcp;

  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL)
      return FALSE;
    for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
      if (CheckForInconsistentBiosources (tmp, vsp, orpp, top))
        return TRUE;
    }
    return FALSE;
  }

  if (!IS_Bioseq (sep))
    return FALSE;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL)
    return FALSE;

  biop = NULL;
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp != NULL) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
    entityID = dcontext.entityID;
    itemID = dcontext.itemID;
    itemtype = OBJ_SEQDESC;
  } else {
    sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_BIOSRC, 0, &fcontext);
    if (sfp != NULL) {
      biop = (BioSourcePtr) sfp->data.value.ptrvalue;
      entityID = fcontext.entityID;
      itemID = fcontext.itemID;
      itemtype = OBJ_SEQFEAT;
    }
  }
  if (biop == NULL)
    return FALSE;
  orp = biop->org;
  if (orp == NULL)
    return FALSE;

  firstorp = *orpp;
  if (firstorp == NULL) {
    *orpp = orp;
    return FALSE;
  }

  if (StringICmp (orp->taxname, firstorp->taxname) == 0)
    return FALSE;

  sev = SEV_ERROR;
  sp = StringStr (orp->taxname, " sp. ");
  if (sp != NULL) {
    len = sp - orp->taxname + 5;
    if (StringNCmp (orp->taxname, firstorp->taxname, len) == 0) {
      sev = SEV_WARNING;
    }
  }

  if (sev == SEV_ERROR) {
    len1 = StringLen (orp->taxname);
    len2 = StringLen (firstorp->taxname);
    len = MIN (len1, len2);
    if (len > 0 && StringNCmp (orp->taxname, firstorp->taxname, len) == 0) {
      sev = SEV_WARNING;
    }
  }

  oldEntityID = gcp->entityID;
  oldItemID = gcp->itemID;
  oldItemtype = gcp->thistype;

  gcp->entityID = entityID;
  gcp->itemID = itemID;
  gcp->thistype = itemtype;

  if (top != NULL) {
    gcp->entityID = top->idx.entityID;
    gcp->itemID = top->idx.itemID;
    gcp->thistype = OBJ_BIOSEQSET;
  }

  /* only report the first one that doesn't match - but might be lower severity if not all are sp. */

  ValidErr (vsp, sev, ERR_SEQ_DESCR_InconsistentBioSources, "Population set contains inconsistent organisms.");

  gcp->entityID = oldEntityID;
  gcp->itemID = oldItemID;
  gcp->thistype = oldItemtype;

  return TRUE;
}

static Boolean CheckForInconsistentMolInfos (SeqEntryPtr sep, ValidStructPtr vsp, MolInfoPtr PNTR mipp, BioseqSetPtr top)

{
  BioseqPtr          bsp;
  BioseqSetPtr       bssp;
  SeqMgrDescContext  dcontext;
  Uint2              entityID = 0, oldEntityID;
  MolInfoPtr         firstmip;
  GatherContextPtr   gcp;
  Uint4              itemID = 0, oldItemID;
  Uint2              itemtype = 0, oldItemtype;
  MolInfoPtr         mip;
  ValNodePtr         sdp;
  SeqEntryPtr        tmp;

  if (sep == NULL || vsp == NULL || mipp == NULL)
    return FALSE;
  gcp = vsp->gcp;

  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL)
      return FALSE;
    for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
      if (CheckForInconsistentMolInfos (tmp, vsp, mipp, top))
        return TRUE;
    }
    return FALSE;
  }

  if (!IS_Bioseq (sep))
    return FALSE;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL)
    return FALSE;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
  if (sdp == NULL) return FALSE;
  mip = (MolInfoPtr) sdp->data.ptrvalue;
  if (mip == NULL || mip->biomol == MOLECULE_TYPE_PEPTIDE) return FALSE;

  firstmip = *mipp;
  if (firstmip == NULL) {
    *mipp = mip;
    return FALSE;
  }

  if (mip->biomol == firstmip->biomol) return FALSE;

  oldEntityID = gcp->entityID;
  oldItemID = gcp->itemID;
  oldItemtype = gcp->thistype;

  gcp->entityID = entityID;
  gcp->itemID = itemID;
  gcp->thistype = itemtype;

  if (top != NULL) {
    gcp->entityID = top->idx.entityID;
    gcp->itemID = top->idx.itemID;
    gcp->thistype = OBJ_BIOSEQSET;
  }

  /* only report the first one that doesn't match */

  ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_InconsistentMolInfoBiomols, "Pop/phy/mut/eco set contains inconsistent MolInfo biomols");

  gcp->entityID = oldEntityID;
  gcp->itemID = oldItemID;
  gcp->thistype = oldItemtype;

  return TRUE;
}

static void LookForMolInfoInconsistency (BioseqSetPtr bssp, ValidStructPtr vsp)

{
  MolInfoPtr    mip = NULL;
  SeqEntryPtr   sep;

  if (bssp == NULL) return;

  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    if (CheckForInconsistentMolInfos (sep, vsp, &mip, bssp))
      return;
  }
}

static Boolean SetHasMolInfo (BioseqSetPtr bssp)

{
  SeqDescrPtr  sdp;

  if (bssp == NULL) return FALSE;

  for (sdp = bssp->descr; sdp != NULL; sdp = sdp->next) {
    if (sdp->choice == Seq_descr_molinfo) return TRUE;
  }

  return FALSE;
}

static void ValidatePopSet (BioseqSetPtr bssp, ValidStructPtr vsp)

{
  BioseqSetPtr  bssp1;
  OrgRefPtr     orp = NULL;
  SeqEntryPtr   sep;

  if (bssp->_class != BioseqseqSet_class_pop_set)
    return;

  if (vsp->is_refseq_in_sep) {
    ValidErr (vsp, SEV_REJECT, ERR_SEQ_PKG_RefSeqPopSet,
              "RefSeq record should not be a Pop-set");
  }

  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    if (!IS_Bioseq_set (sep)) continue;
    bssp1 = sep->data.ptrvalue;
    if (bssp1 == NULL) continue;

    if (bssp1->_class == BioseqseqSet_class_genbank) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_InternalGenBankSet,
                "Bioseq-set contains internal GenBank Bioseq-set");
    }
  }

  if (SetHasMolInfo (bssp)) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_MisplacedMolInfo, "Pop set has MolInfo on set");
  }

  LookForMolInfoInconsistency (bssp, vsp);

  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    if (CheckForInconsistentBiosources (sep, vsp, &orp, bssp))
      return;
  }
}

typedef struct mutsetsrcdata {
  CharPtr  taxname;
  Int2     num_not_mut_origin;
  Boolean  failed;
} MutSetSrcData, PNTR MutSetSrcPtr;

static void CheckMutSetSources (BioSourcePtr biop, Pointer userdata)

{
  MutSetSrcPtr  mssp;
  OrgRefPtr     orp;

  if (biop == NULL || userdata == NULL) return;
  mssp = (MutSetSrcPtr) userdata;

  orp = biop->org;
  if (orp == NULL || StringHasNoText (orp->taxname)) return;
  if (mssp->taxname == NULL) {
    mssp->taxname = orp->taxname;
  } else if (StringCmp (mssp->taxname, orp->taxname) != 0) {
    mssp->failed = TRUE;
  }
  if (biop->origin != ORG_MUT) {
    (mssp->num_not_mut_origin)++;
  }
}

static void ValidateMutSet (BioseqSetPtr bssp, ValidStructPtr vsp)

{
  BioseqSetPtr   bssp1;
/*  MutSetSrcData  mssd; */
  SeqEntryPtr    sep;

  if (bssp->_class != BioseqseqSet_class_mut_set)
    return;

  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    if (!IS_Bioseq_set (sep)) continue;
    bssp1 = sep->data.ptrvalue;
    if (bssp1 == NULL) continue;

    if (bssp1->_class == BioseqseqSet_class_genbank) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_InternalGenBankSet,
                "Bioseq-set contains internal GenBank Bioseq-set");
    }
  }

  if (SetHasMolInfo (bssp)) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_MisplacedMolInfo, "Mut set has MolInfo on set");
  }

  LookForMolInfoInconsistency (bssp, vsp);

  /* error is currently suppressed
  MemSet ((Pointer) &mssd, 0, sizeof (MutSetSrcData));
  VisitBioSourcesInSet (bssp, (Pointer) &mssd, CheckMutSetSources);
  if (mssd.failed) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_InconsistentBioSources, "Mutation set contains inconsistent organisms.");
  }
  if (mssd.num_not_mut_origin > 1) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_InconsistentBioSources, "Mutation set contains more than one non-mutant organism.");
  }
  */
}

static void ValidateGenbankSet (BioseqSetPtr bssp, ValidStructPtr vsp)

{
  BioseqSetPtr    bssp1;
  SeqEntryPtr     sep;

  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    if (!IS_Bioseq_set (sep)) continue;
    bssp1 = sep->data.ptrvalue;
    if (bssp1 == NULL) continue;

    if (bssp1->_class == BioseqseqSet_class_genbank) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_InternalGenBankSet,
                "Bioseq-set contains internal GenBank Bioseq-set");
    }
  }

  if (SetHasMolInfo (bssp)) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_MisplacedMolInfo, "Genbank set has MolInfo on set");
  }
}

static void ValidatePhyEcoWgsSet (BioseqSetPtr bssp, ValidStructPtr vsp)

{
  BioseqSetPtr  bssp1;
  SeqEntryPtr   sep;

  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    if (!IS_Bioseq_set (sep)) continue;
    bssp1 = sep->data.ptrvalue;
    if (bssp1 == NULL) continue;

    if (bssp1->_class == BioseqseqSet_class_genbank) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_InternalGenBankSet,
                "Bioseq-set contains internal GenBank Bioseq-set");
    }
  }

  if (SetHasMolInfo (bssp)) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_MisplacedMolInfo, "Phy/eco/wgs set has MolInfo on set");
  }

  LookForMolInfoInconsistency (bssp, vsp);
}

static void ValidateGenProdSet (BioseqSetPtr bssp, ValidStructPtr vsp)

{
  BioseqPtr       bsp;
  BioseqPtr       cdna;
  SeqMgrFeatContext fcontext;
  GatherContextPtr gcp = NULL;
  CharPtr         loc = NULL;
  SeqFeatPtr      mrna;
  Uint2           olditemtype = 0;
  Uint4           olditemid = 0;
  SeqEntryPtr     sep;
  SeqIdPtr        sip;

  if (bssp->_class != BioseqseqSet_class_gen_prod_set)
    return;

  if (bssp->annot != NULL) {
    ValidErr (vsp, SEV_REJECT, ERR_SEQ_PKG_GenomicProductPackagingProblem, "Seq-annot packaged directly on genomic product set");
  }

  sep = bssp->seq_set;
  if (!IS_Bioseq (sep))
    return;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL)
    return;

  gcp = vsp->gcp;
  if (gcp == NULL)
    return;
  olditemid = gcp->itemID;
  olditemtype = gcp->thistype;

  if (vsp->useSeqMgrIndexes) {
    mrna = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_mRNA, &fcontext);
    while (mrna != NULL) {
      cdna = BioseqFindFromSeqLoc (mrna->product);
      if (cdna == NULL) {
        gcp->itemID = mrna->idx.itemID;
        gcp->thistype = OBJ_SEQFEAT;
        loc = SeqLocPrint (mrna->product);
        if (loc == NULL) {
          loc = StringSave ("?");
        }
        sip = SeqLocId (mrna->product);
        /* okay to have far RefSeq product */
        if (sip == NULL || sip->choice != SEQID_OTHER) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_GenomicProductPackagingProblem, "Product of mRNA feature (%s) not packaged in genomic product set", loc);
        }
        MemFree (loc);
      }
      mrna = SeqMgrGetNextFeature (bsp, mrna, 0, FEATDEF_mRNA, &fcontext);
    }
  }

  if (SetHasMolInfo (bssp)) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_MisplacedMolInfo, "GenProd set has MolInfo on set");
  }

  gcp->itemID = olditemid;
  gcp->thistype = olditemtype;
}

static void NestedSetProc (BioseqSetPtr bssp, Pointer userdata)

{
  ValidStructPtr   vsp;
  GatherContextPtr gcp = NULL;

  if (bssp == NULL) return;

  /* pop/phy/mut/eco set can contain up to nuc-prot sets */
  switch (bssp->_class) {
  case BioseqseqSet_class_nuc_prot:
  case BioseqseqSet_class_segset:
  case BioseqseqSet_class_parts:
    return;
  default:
    break;
  }

  vsp = (ValidStructPtr) userdata;
  if (vsp == NULL) return;
  gcp = vsp->gcp;
  if (gcp == NULL) return;

  ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_ImproperlyNestedSets, "Nested sets within Pop/Phy/Mut/Eco/Wgs set");
}

static void CheckForNestedSets (BioseqSetPtr bssp, Pointer userdata)

{
  SeqEntryPtr  sep;

  if (bssp == NULL) return;

  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    if (!IS_Bioseq_set (sep)) continue;
    VisitSetsInSep (sep, userdata, NestedSetProc);
  }
}

static void FindDBlinkUserObject (SeqDescrPtr sdp, Pointer userdata)

{
  GatherContextPtr  gcp;
  ObjectIdPtr       oip;
  UserObjectPtr     uop;
  ValidStructPtr    vsp;

  if (sdp == NULL || sdp->choice != Seq_descr_user) return;
  uop = (UserObjectPtr) sdp->data.ptrvalue;
  if (uop == NULL) return;
  oip = uop->type;
  if (oip == NULL) return;

  if (StringICmp (oip->str, "DBLink") != 0) return;

  vsp = (ValidStructPtr) userdata;
  if (vsp == NULL) return;
  gcp = vsp->gcp;
  if (gcp == NULL) return;

  ValidErr (vsp, SEV_REJECT, ERR_SEQ_DESCR_DBLinkProblem, "DBLink user object should not be on this set");
}

static void ShouldHaveNoDblink (BioseqSetPtr bssp, Pointer userdata)

{
  ValidStructPtr  vsp;

  if (bssp == NULL) return;
  vsp = (ValidStructPtr) userdata;
  if (vsp == NULL) return;

  VisitDescriptorsOnSet (bssp, vsp, FindDBlinkUserObject);
}

static void ValidateBioseqSet (GatherContextPtr gcp)

{
  BioseqSetPtr    bssp;
  ValidStructPtr  vsp;
  SeqEntryPtr     sep;
  SeqDescrPtr     sdp;
  CharPtr         str;
  Boolean         has_title = FALSE;

  vsp = (ValidStructPtr) (gcp->userdata);
  bssp = (BioseqSetPtr) (gcp->thisitem);
  vsp->bssp = bssp;
  vsp->bsp = NULL;
  vsp->descr = NULL;
  vsp->sfp = NULL;

  if (vsp->non_ascii_chars) {   /* non_ascii chars in AsnRead step */
    ValidErr (vsp, SEV_ERROR, ERR_GENERIC_NonAsciiAsn, "Non-ascii chars in input ASN.1 strings");
    vsp->non_ascii_chars = FALSE;       /* only do once */
  }

  vsp->nuccnt = 0;
  vsp->segcnt = 0;
  vsp->protcnt = 0;

  sep = gcp->sep;

  SeqEntryExplore (sep, (Pointer) vsp, ValidateSetContents);

  switch (bssp->_class) {
  case BioseqseqSet_class_not_set:
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_BioseqSetClassNotSet, "Bioseq_set class not set");
    break;
  case BioseqseqSet_class_nuc_prot:
    if (vsp->nuccnt == 0) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_PKG_NucProtProblem, "No nucleotides in nuc-prot set");
    }
    if (vsp->protcnt == 0) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_PKG_NucProtProblem, "No proteins in nuc-prot set");
    }
    if (vsp->nuccnt > 1 && vsp->segcnt == 0) {
      ValidErr (vsp, SEV_REJECT, ERR_SEQ_PKG_NucProtProblem, "Multiple unsegmented nucleotides in nuc-prot set");
    }
    ValidateNucProtSet (bssp, vsp);
    break;
  case BioseqseqSet_class_segset:
    if (vsp->segcnt == 0) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_PKG_SegSetProblem, "No segmented Bioseq in segset");
    }
    ValidateSegmentedSet (bssp, vsp);
    break;
  case BioseqseqSet_class_conset:
    if (vsp->indexerVersion && (! vsp->is_refseq_in_sep)) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_PKG_ConSetProblem, "Set class should not be conset");
    }
    break;
  case BioseqseqSet_class_parts:
    ValidatePartsSet (bssp, vsp);
    break;
  case BioseqseqSet_class_genbank:
    ValidateGenbankSet (bssp, vsp);
    ShouldHaveNoDblink (bssp, vsp);
    break;
  case BioseqseqSet_class_pop_set:
    ValidatePopSet (bssp, vsp);
    CheckForNestedSets (bssp, vsp);
    ShouldHaveNoDblink (bssp, vsp);
    break;
  case BioseqseqSet_class_mut_set:
    ValidateMutSet (bssp, vsp);
    CheckForNestedSets (bssp, vsp);
    ShouldHaveNoDblink (bssp, vsp);
    break;
  case BioseqseqSet_class_phy_set:
  case BioseqseqSet_class_eco_set:
  case BioseqseqSet_class_wgs_set:
    ValidatePhyEcoWgsSet (bssp, vsp);
    CheckForNestedSets (bssp, vsp);
    ShouldHaveNoDblink (bssp, vsp);
    break;
  case BioseqseqSet_class_small_genome_set:
    ValidatePhyEcoWgsSet (bssp, vsp);
    CheckForNestedSets (bssp, vsp);
    break;
  case BioseqseqSet_class_gen_prod_set:
    ValidateGenProdSet (bssp, vsp);
    break;
  /*
  case BioseqseqSet_class_other:
    ValidErr (vsp, SEV_REJECT, ERR_SEQ_PKG_GenomicProductPackagingProblem, "Genomic product set class incorrectly set to other");
    break;
  */
  default:
    if (!((vsp->nuccnt) || (vsp->protcnt))) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_EmptySet, "No Bioseqs in this set");
    }
    break;
  }

  switch (bssp->_class) {
  case BioseqseqSet_class_pop_set:
  case BioseqseqSet_class_mut_set:
  case BioseqseqSet_class_phy_set:
  case BioseqseqSet_class_eco_set:
    for (sdp = bssp->descr; sdp != NULL; sdp = sdp->next) {
      if (sdp->choice != Seq_descr_title) continue;
      str = (CharPtr) sdp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      has_title = TRUE;
    }
    if (! has_title && (vsp->is_insd_in_sep || vsp->is_refseq_in_sep)) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_MissingSetTitle, "Pop/Phy/Mut/Eco set does not have title");
    }
    sep = bssp->seq_set;
    if (sep == NULL) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_EmptySet, "Pop/Phy/Mut/Eco set has no components");
    } else if (sep->next == NULL) {
      if (VisitAlignmentsInSep (gcp->sep, NULL, NULL) == 0) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_SingleItemSet, "Pop/Phy/Mut/Eco set has only one component and no alignments");
      }
    }
    break;
  default:
    break;
  }
}

static Boolean SuppressTrailingXMessage (BioseqPtr bsp)
{
  ByteStorePtr    bs;
  SeqFeatPtr      cds;
  Boolean         hasstar;
  Int4            len;
  MolInfoPtr      mip;
  SeqDescrPtr     sdp;
  CharPtr         str;

  cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
  if (cds != NULL) {
    bs = ProteinFromCdRegionEx (cds, TRUE, FALSE);
    if (bs != NULL) {
      str = BSMerge (bs, NULL);
      BSFree (bs);
      hasstar = FALSE;
      if (str != NULL) {
        len = StringLen (str);
        if (len > 1 && str[len - 1] == '*') {
          hasstar = TRUE;
        }
      }
      MemFree (str);
      return hasstar;
    }
  }
  sdp = BioseqGetSeqDescr (bsp, Seq_descr_molinfo, NULL);
  if (sdp != NULL) {
    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip != NULL) {
      if (mip->completeness == 4 || mip->completeness == 5)
        return TRUE;
    }
  }
  return FALSE;
}

static void LookForSecondaryConflict (ValidStructPtr vsp, GatherContextPtr gcp, CharPtr accn, ValNodePtr extra_acc)
{
  CharPtr         str;
  ValNodePtr      vnp;

  if (vsp == NULL || gcp == NULL)
    return;
  if (StringHasNoText (accn))
    return;
  for (vnp = extra_acc; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str))
      continue;
    if (StringICmp (accn, str) == 0) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_BadSecondaryAccn, "%s used for both primary and secondary accession", accn);
    }
  }
}

static void CheckSegBspAgainstParts (ValidStructPtr vsp, GatherContextPtr gcp, BioseqPtr bsp)
{
  BioseqSetPtr    bssp;
  Boolean         is_odd;
  BioseqPtr       part;
  SeqEntryPtr     sep;
  SeqIdPtr        sip;
  SeqLocPtr       slp;
  BioseqPtr       vbsp;

  if (vsp == NULL || gcp == NULL || bsp == NULL)
    return;
  if (!vsp->useSeqMgrIndexes)
    return;

  if (bsp->repr != Seq_repr_seg || bsp->seq_ext_type != 1 || bsp->seq_ext == NULL)
    return;

  sep = bsp->seqentry;
  if (sep == NULL)
    return;
  sep = sep->next;
  if (sep == NULL)
    return;
  if (!IS_Bioseq_set (sep))
    return;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp == NULL)
    return;
  if (bssp->_class != BioseqseqSet_class_parts)
    return;

  is_odd = FALSE;
  for (slp = (ValNodePtr) bsp->seq_ext; slp != NULL; slp = slp->next) {
    is_odd = (! is_odd);
    if (is_odd) {
      if (slp->choice == SEQLOC_NULL) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_BadSegmentedSeq, "Odd segmented component is not expected to be NULL");
      }
    } else {
      if (slp->choice != SEQLOC_NULL) {
        vbsp = BioseqFindFromSeqLoc (slp);
        if (vbsp != NULL) {
          if (vbsp->repr != Seq_repr_virtual) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_BadSegmentedSeq, "Even segmented component is expected to be NULL or VIRTUAL");
          }
        }
      }
    }
  }

  sep = bssp->seq_set;
  for (slp = (ValNodePtr) bsp->seq_ext; slp != NULL; slp = slp->next) {
    if (slp->choice == SEQLOC_NULL)
      continue;
    if (sep == NULL) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_PartsOutOfOrder, "Parts set does not contain enough Bioseqs");
      return;
    }
    if (IS_Bioseq (sep)) {
      part = (BioseqPtr) sep->data.ptrvalue;
      sip = SeqLocId (slp);
      if (sip != NULL && part != NULL) {
        if (!SeqIdIn (sip, part->id)) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_PartsOutOfOrder, "Segmented bioseq seq_ext does not correspond to parts packaging order");
          return;
        }
      }
    } else {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_PartsOutOfOrder, "Parts set component is not Bioseq");
      return;
    }
    sep = sep->next;
  }
  if (sep != NULL) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_PartsOutOfOrder, "Parts set contains too many Bioseqs");
  }
}

/*****************************************************************************
*
*   ValidateBioseqHist(gcp)
*      Validate one Bioseq Seq-hist
*
*****************************************************************************/
static void ValidateBioseqHist (GatherContextPtr gcp)

{
  BioseqPtr       bsp;
  Int4            gi = 0;
  SeqHistPtr      hist;
  SeqIdPtr        sip;
  ValidStructPtr  vsp;

  if (gcp == NULL) return;
  vsp = (ValidStructPtr) (gcp->userdata);
  bsp = (BioseqPtr) (gcp->thisitem);
  vsp->bsp = bsp;
  vsp->descr = NULL;
  vsp->sfp = NULL;
  vsp->bssp = (BioseqSetPtr) (gcp->parentitem);
  vsp->bsp_partial_val = 0;

  if (bsp == NULL) return;
  hist = bsp->hist;
  if (hist == NULL) return;

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GI) {
      gi = (Int4) sip->data.intvalue;
    }
  }
  if (gi == 0) return;

  if (hist->replaced_by_ids != NULL && hist->replaced_by_date != NULL) {

    for (sip = hist->replaced_by_ids; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_GI) {
        if (gi == (Int4) sip->data.intvalue) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_HistoryGiCollision, "Replaced by gi (%ld) is same as current Bioseq", (long) gi);
        }
      }
    }
  }

  if (hist->replace_ids != NULL && hist->replace_date != NULL) {

    for (sip = hist->replace_ids; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_GI) {
        if (gi == (Int4) sip->data.intvalue) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_HistoryGiCollision, "Replaces gi (%ld) is same as current Bioseq", (long) gi);
        }
      }
    }
  }
}

/*****************************************************************************
*
*   ValidateBioseqInst(gcp)
*      Validate one Bioseq Seq-inst
*
*****************************************************************************/
static Boolean IsTpa (
  BioseqPtr bsp,
  Boolean has_tpa_assembly,
  BoolPtr isRefSeqP
)

{
  DbtagPtr  dbt;
  Boolean   has_bankit = FALSE;
  Boolean   has_genbank = FALSE;
  Boolean   has_gi = FALSE;
  Boolean   has_local = FALSE;
  Boolean   has_refseq = FALSE;
  Boolean   has_smart = FALSE;
  Boolean   has_tpa = FALSE;
  SeqIdPtr  sip;

  if (bsp == NULL || bsp->id == NULL) return FALSE;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    switch (sip->choice) {
      case SEQID_LOCAL :
        has_local = TRUE;
        break;
      case SEQID_GENBANK :
      case SEQID_EMBL :
      case SEQID_DDBJ :
        has_genbank = TRUE;
        break;
      case SEQID_OTHER :
        has_refseq = TRUE;
        if (isRefSeqP != NULL) {
          *isRefSeqP = TRUE;
        }
        break;
      case SEQID_GI :
        has_gi = TRUE;
        break;
      case SEQID_TPG :
      case SEQID_TPE :
      case SEQID_TPD :
        has_tpa = TRUE;
        break;
      case SEQID_GENERAL :
        dbt = (DbtagPtr) sip->data.ptrvalue;
        if (dbt != NULL) {
          if (StringICmp (dbt->db, "BankIt") == 0) {
            has_bankit = TRUE;
          }
          if (StringICmp (dbt->db, "TMSMART") == 0) {
            has_smart = TRUE;
          }
        }
        break;
      case SEQID_GPIPE :
        break;
      default :
        break;
    }
  }

  if (has_genbank) return FALSE;
  if (has_tpa) return TRUE;
  if (has_refseq) return FALSE;
  if (has_bankit && has_tpa_assembly) return TRUE;
  if (has_smart && has_tpa_assembly) return TRUE;
  if (has_gi) return FALSE;
  if (has_local && has_tpa_assembly) return TRUE;

  return FALSE;
}

static void ValidateIDSetAgainstDb (GatherContextPtr gcp, ValidStructPtr vsp, BioseqPtr bsp)

{
  SeqIdPtr        sip, sipset;
  SeqIdPtr        gbId = NULL;
  SeqIdPtr        dbGbId;
  DbtagPtr        generalID = NULL;
  DbtagPtr        dbGeneralID;
  Int4            gi = 0;
  Int4            dbGI;
  Char            oldGenID [128], newGenID [128];

  if (gcp != NULL && vsp != NULL && bsp != NULL && vsp->validateIDSet) {
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      switch (sip->choice) {
        case SEQID_GENBANK:
          gbId = sip;
          break;
        case SEQID_GI :
          gi = (Int4) sip->data.intvalue;
          break;
        case SEQID_GENERAL :
          generalID = (DbtagPtr) sip->data.ptrvalue;
          break;
        default :
          break;
      }
    }
    if (gi == 0 && gbId != NULL) {
      gi = GetGIForSeqId (gbId);
    }
    if (gi > 0) {
      sipset = GetSeqIdSetForGI (gi);
      if (sipset != NULL) {
        dbGI = 0;
        dbGbId = NULL;
        dbGeneralID = NULL;
        oldGenID [0] = '\0';
        newGenID [0] = '\0';
        for (sip = sipset; sip != NULL; sip = sip->next) {
          switch (sip->choice) {
            case SEQID_GI :
              dbGI = (Int4) sip->data.intvalue;
              break;
            case SEQID_GENBANK:
              dbGbId = sip;
              break;
            case SEQID_GENERAL :
              dbGeneralID = (DbtagPtr) sip->data.ptrvalue;
              break;
            default :
              break;
          }
        }
        if (dbGI != gi) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_UnexpectedIdentifierChange, "New gi number (%ld) does not match one in NCBI sequence repository (%ld)", (long) gi, (long) dbGI);
        }
        if (gbId != NULL && dbGbId != NULL) {
          if (! SeqIdMatch (gbId, dbGbId)) {
            SeqIdWrite (dbGbId, oldGenID, PRINTID_FASTA_SHORT, sizeof (oldGenID));
            SeqIdWrite (gbId, newGenID, PRINTID_FASTA_SHORT, sizeof (newGenID));
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_UnexpectedIdentifierChange, "New accession (%s) does not match one in NCBI sequence repository (%s) on gi (%ld)", newGenID, oldGenID, (long) gi);
          }
        } else if (gbId != NULL) {
          SeqIdWrite (gbId, newGenID, PRINTID_FASTA_SHORT, sizeof (newGenID));
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_UnexpectedIdentifierChange, "Gain of accession (%s) on gi (%ld) compared to the NCBI sequence repository", newGenID, (long) gi);
        } else if (dbGbId != NULL) {
          SeqIdWrite (dbGbId, oldGenID, PRINTID_FASTA_SHORT, sizeof (oldGenID));
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_UnexpectedIdentifierChange, "Loss of accession (%s) on gi (%ld) compared to the NCBI sequence repository", oldGenID, (long) gi);
        }
        if (generalID != NULL && dbGeneralID != NULL) {
          if (! DbtagMatch (generalID, dbGeneralID)) {
            DbtagLabel (dbGeneralID, oldGenID, sizeof (oldGenID));
            DbtagLabel (generalID, newGenID, sizeof (newGenID));
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_UnexpectedIdentifierChange, "New general ID (%s) does not match one in NCBI sequence repository (%s) on gi (%ld)", newGenID, oldGenID, (long) gi);
          }
        } else if (generalID != NULL) {
          DbtagLabel (generalID, newGenID, sizeof (newGenID));
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_UnexpectedIdentifierChange, "Gain of general ID (%s) on gi (%ld) compared to the NCBI sequence repository", newGenID, (long) gi);
        } else if (dbGeneralID != NULL) {
          DbtagLabel (dbGeneralID, oldGenID, sizeof (oldGenID));
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_UnexpectedIdentifierChange, "Loss of general ID (%s) on gi (%ld) compared to the NCBI sequence repository", oldGenID, (long) gi);
        }
      }
      SeqIdSetFree (sipset);
    }
  }
}

typedef struct enrun {
  GatherContextPtr  gcp;
  ValidStructPtr    vsp;
  Int4              ncount;
  Int4              maxrun;
  Int4              seqpos;
  Int4              gapcount;
  Boolean           showAll;
  Boolean           inNrun;
  Boolean           isWGS;
} RunOfNs, PNTR RunOfNsPtr;

static void LIBCALLBACK CountAdjacentProc (CharPtr sequence, Pointer userdata)

{
  Char              ch;
  GatherContextPtr  gcp;
  RunOfNsPtr        ronp;
  CharPtr           str;
  ValidStructPtr    vsp;

  ronp = (RunOfNsPtr) userdata;
  if (sequence == NULL || ronp == NULL) return;

  str = sequence;
  ch = *str;
  while (ch != '\0') {
    (ronp->seqpos)++;
    if (ch == 'N') {
      (ronp->ncount)++;
      if (ronp->ncount > ronp->maxrun) {
        ronp->maxrun = ronp->ncount;
      }
      ronp->inNrun = TRUE;
    } else {
      if (ch == '-') {
        (ronp->gapcount)++;
      }
      if (ronp->inNrun && ronp->showAll && ronp->isWGS && ronp->ncount >= 20 && ronp->seqpos > ronp->ncount + 1) {
        vsp = ronp->vsp;
        gcp = ronp->gcp;
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_InternalNsInSeqRaw, "Run of %ld Ns in raw sequence starting at base %ld",
                  (long) ronp->ncount, (long) (ronp->seqpos - ronp->ncount));
      } else if (ronp->inNrun && ronp->showAll && ronp->ncount >= 100 && ronp->seqpos > ronp->ncount + 1) {
        vsp = ronp->vsp;
        gcp = ronp->gcp;
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_InternalNsInSeqRaw, "Run of %ld Ns in raw sequence starting at base %ld",
                  (long) ronp->ncount, (long) (ronp->seqpos - ronp->ncount));
      }
      ronp->ncount = 0;
      ronp->inNrun = FALSE;
    }
    str++;
    ch = *str;
  }
}


static Int4 CountAdjacentNsInInterval (GatherContextPtr gcp, BioseqPtr bsp, Int4 from, Int4 to)
{
  SeqLocPtr slp;
  RunOfNs   ron;

  if (bsp == NULL || from < 0 || to < from || ISA_aa (bsp->mol)) {
    return 0;
  }

  slp = SeqLocIntNew (from, to, Seq_strand_plus, bsp->id);
  ron.gcp = gcp;
  ron.vsp = (ValidStructPtr) (gcp->userdata);
  ron.ncount = 0;
  ron.maxrun = 0;
  ron.seqpos = 0;
  ron.gapcount = 0;
  ron.showAll = FALSE;
  ron.inNrun = FALSE;
  ron.isWGS = FALSE;
  SeqPortStreamLoc (slp, STREAM_EXPAND_GAPS, (Pointer) &ron, CountAdjacentProc);
  slp = SeqLocFree (slp);
  return ron.maxrun;
}


static Boolean HasUnparsedBrackets (CharPtr title)

{
  CharPtr  str;

  if (StringHasNoText (title)) return FALSE;

  str = StringChr (title, '[');
  if (str == NULL) return FALSE;
  str = StringChr (str, '=');
  if (str == NULL) return FALSE;
  str = StringChr (str, ']');
  if (str == NULL) return FALSE;
  return TRUE;
}

static CharPtr GetSequencePlusGapByFeature (SeqFeatPtr sfp)

{
  Int4     len;
  CharPtr  str = NULL;

  if (sfp == NULL) return NULL;
  len = SeqLocLen (sfp->location);
  if (len > 0 && len < MAXALLOC) {
    str = MemNew (sizeof (Char) * (len + 2));
    if (str != NULL) {
      SeqPortStreamLoc (sfp->location, EXPAND_GAPS_TO_DASHES, (Pointer) str, NULL);
    }
  }

  return str;
}

static Boolean IsWgsIntermediate (SeqEntryPtr sep)

{
  BioseqPtr    bsp;
  Boolean      has_gi = FALSE, is_other = FALSE, is_wgs = FALSE;
  MolInfoPtr   mip;
  SeqDescrPtr  sdp;
  SeqIdPtr     sip;

  bsp = FindNucBioseq (sep);
  if (bsp == NULL) return FALSE;

  for (sdp = bsp->descr; sdp != NULL; sdp = sdp->next) {
    if (sdp->choice != Seq_descr_molinfo) continue;
    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip == NULL) continue;
    if (mip->tech == MI_TECH_wgs) {
      is_wgs = TRUE;
    }
  }
  if (! is_wgs) return FALSE;

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_OTHER) {
      is_other = TRUE;
    } else if (sip->choice == SEQID_GI) {
      has_gi = TRUE;
    }
  }
  if (! is_other) return FALSE;
  if (has_gi) return FALSE;

  return TRUE;
}

typedef struct reusedata {
  CharPtr  seqidstr;
  Int4     from;
  Int4     to;
} ReuseData, PNTR ReuseDataPtr;

static int LIBCALLBACK SortVnpByDeltaLoc (VoidPtr ptr1, VoidPtr ptr2)

{
  int           compare;
  ReuseDataPtr  rdp1, rdp2;
  ValNodePtr    vnp1, vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  rdp1 = (ReuseDataPtr) vnp1->data.ptrvalue;
  rdp2 = (ReuseDataPtr) vnp2->data.ptrvalue;
  if (rdp1 == NULL || rdp2 == NULL) return 0;

  compare = StringICmp (rdp1->seqidstr, rdp2->seqidstr);
  if (compare > 0) {
    return 1;
  } else if (compare < 0) {
    return -1;
  }

  if (rdp1->from > rdp2->from) {
    return 1;
  } else if (rdp1->from < rdp2->from) {
    return -1;
  }

  if (rdp1->to > rdp2->to) {
    return 1;
  } else if (rdp1->to < rdp2->to) {
    return -1;
  }

  return 0;
}

static void CheckDeltaForReuse (ValidStructPtr vsp, GatherContextPtr gcp, BioseqPtr bsp)

{
  Char          buf [80];
  ValNodePtr    head = NULL;
  ValNodePtr    last = NULL;
  ReuseDataPtr  lastrdp = NULL;
  ReuseDataPtr  rdp;
  SeqIntPtr     sintp;
  SeqIdPtr      sip;
  SeqLocPtr     slp;
  ValNodePtr    vnp_dsp, vnp_r;

  if (vsp == NULL || gcp == NULL || bsp == NULL) return;

  for (vnp_dsp = (ValNodePtr) bsp->seq_ext; vnp_dsp != NULL; vnp_dsp = vnp_dsp->next) {
    if (vnp_dsp->choice != 1) continue;
    slp = (SeqLocPtr) vnp_dsp->data.ptrvalue;
    if (slp == NULL) continue;
    if (slp->choice != SEQLOC_INT) continue;
    sintp = (SeqIntPtr) slp->data.ptrvalue;
    if (sintp == NULL) continue;
    sip = sintp->id;
    if (sip == NULL) continue;
    if (! SeqIdWrite (sip, buf, PRINTID_FASTA_SHORT, sizeof (buf) - 1)) continue;
    rdp = (ReuseDataPtr) MemNew (sizeof (ReuseData));
    if (rdp == NULL) continue;
    rdp->seqidstr = StringSave (buf);
    rdp->from = sintp->from;
    rdp->to = sintp->to;
    vnp_r = ValNodeAddPointer (&last, 0, (Pointer) rdp);
    if (head == NULL) {
      head = vnp_r;
    }
    last = vnp_r;
  }

  if (head == NULL) return;

  head = ValNodeSort (head, SortVnpByDeltaLoc);

  for (vnp_r = head; vnp_r != NULL; vnp_r = vnp_r->next) {
    rdp = (ReuseDataPtr) vnp_r->data.ptrvalue;
    if (rdp == NULL) continue;
    if (lastrdp != NULL) {
      if (StringICmp (lastrdp->seqidstr, rdp->seqidstr) == 0) {
        if (lastrdp->to >= rdp->from && lastrdp->from <= rdp->to) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_OverlappingDeltaRange,
                    "Overlapping delta range %ld-%ld and %ld-%ld on a Bioseq %s",
                    (long) rdp->from + 1, (long) rdp->to + 1, (long) lastrdp->from + 1,
                    (long) lastrdp->to + 1, rdp->seqidstr);
        }
      }
    }
    lastrdp = rdp;
  }

  for (vnp_r = head; vnp_r != NULL; vnp_r = vnp_r->next) {
    rdp = (ReuseDataPtr) vnp_r->data.ptrvalue;
    if (rdp == NULL) continue;
    rdp->seqidstr = MemFree (rdp->seqidstr);
  }
  ValNodeFreeData (head);
}

static CharPtr legal_refgene_status_strings [] = {
  "Inferred",
  "Provisional",
  "Predicted",
  "Validated",
  "Reviewed",
  "Model",
  "WGS",
  "Pipeline",
  NULL
};


static void ReportLongSeqId (SeqIdPtr sip, ValidStructPtr vsp, Int4 max_len)
{
  Int4 id_len = 0;
  CharPtr id_txt;

  if (sip == NULL || vsp == NULL || IsNCBIFileID(sip)) {
    return;
  }

  id_len = SeqIdLabelLen(sip, PRINTID_FASTA_SHORT);
  if (id_len > max_len) {
    id_txt = SeqIdWholeLabel (sip, PRINTID_FASTA_SHORT);
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_BadSeqIdFormat, "Sequence ID is unusually long (%d): %s", id_len, id_txt);
    id_txt = MemFree (id_txt);
  }

}


static Boolean SequenceHasGaps (BioseqPtr bsp)
{
  SeqMgrFeatContext context;
  SeqFeatPtr sfp;

  if (bsp == NULL) {
    return FALSE;
  }
  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_gap, &context);
  if (sfp == NULL) {
    return FALSE;
  } else {
    return TRUE;
  }
}


static void ValidateBioseqInst (GatherContextPtr gcp)
{
  Boolean         retval = TRUE;
  Int2            i, start_at, num;
  Boolean         errors[4], check_alphabet;
  static char    *repr[8] = {
    "virtual", "raw", "segmented", "constructed",
    "reference", "consensus", "map", "delta"
  };
  /*
  SeqPortPtr      spp;
  */
  Int2            residue, x, termination, gapchar;
  Boolean         gapatstart;
  Int4            len, divisor = 1, len2, len3;
  ValNode         head, vn;
  ValNodePtr      vnp, idlist;
  BioseqContextPtr bcp;
  Boolean         got_partial, is_invalid;
  int             seqtype, terminations, dashes;
  ValidStructPtr  vsp;
  BioseqPtr       bsp, bsp2;
  SeqIdPtr        sip1, sip2, sip3;
  SeqLocPtr       slp;
  SeqIntPtr       sintp;
  Char            buf1[41], buf2[41];
  SeqLitPtr       slitp;
  SeqCodeTablePtr sctp;
  MolInfoPtr      mip = NULL;
  BioSourcePtr    biop = NULL;
  OrgRefPtr       orp;
  SeqMgrDescContext context;
  SeqFeatPtr      cds;
  CdRegionPtr     crp;
  GBBlockPtr      gbp;
  GeneRefPtr      grp;
  SeqFeatPtr      gene;
  SeqMgrFeatContext genectxt;
  CharPtr         genelbl = NULL;
  SeqFeatPtr      prot;
  SeqMgrFeatContext protctxt;
  CharPtr         protlbl = NULL;
  TextSeqIdPtr    tsip;
  CharPtr         ptr, last, str, title, buf, bufplus;
  Uint1           lastchoice;
  Char            ch;
  Boolean         multitoken;
  Boolean         hasGi = FALSE;
  SeqHistPtr      hist;
  Boolean         hist_asm_missing = FALSE;
  IntFuzzPtr      ifp;
  Boolean         in_gap;
  Boolean         in_N;
  Boolean         in_nps;
  Boolean         isActiveFin = FALSE;
  Boolean         isDDBJ = FALSE;
  Boolean         isDraft = FALSE;
  Boolean         isEMBL = FALSE;
  Boolean         isFullTop = FALSE;
  Boolean         isGB = FALSE;
  Boolean         isGIBBMT = FALSE;
  Boolean         isGIBBSQ = FALSE;
  Boolean         isPatent = FALSE;
  Boolean         isPDB = FALSE;
  Boolean         isPreFin = FALSE;
  Boolean         isNC = FALSE;
  Boolean         isNG = FALSE;
  Boolean         isNTorNC = FALSE;
  Boolean         isNZ;
  Boolean         is_gps = FALSE;
  Boolean         isRefSeq = FALSE;
  Boolean         isSwissProt = FALSE;
  Boolean         is_genome_assembly;
  Boolean         is_finished_status;
  Boolean         only_local = TRUE;
  Boolean         isLRG = FALSE;
  ValNodePtr      keywords;
  Boolean         last_is_gap;
  Boolean         non_interspersed_gaps;
  Int2            num_adjacent_gaps;
  Int2            num_gaps;
  Boolean         reportFastaBracket;
  SeqFeatPtr      sfp;
  SeqEntryPtr     sep;
  ErrSev          sev;
  DbtagPtr        dbt;
  SeqIdPtr        sip;
  Int2            trailingX = 0;
  Int2            numletters, numdigits, numunderscores;
  Boolean         letterAfterDigit, badIDchars;
  EMBLBlockPtr    ebp;
  SeqDescrPtr     sdp;
  SeqMgrDescContext dcontext;
  Uint2           oldEntityID, oldItemtype;
  Uint4           oldItemID;
  size_t          buflen = 1001;
  ItemInfo        ii;
  Uint1           tech;
  Uint2           olditemtype = 0;
  Uint4           olditemid = 0;
  ObjValNodePtr   ovp;
  BioseqSetPtr    bssp;
  UserObjectPtr   uop;
  UserFieldPtr    ufp;
  ObjectIdPtr     oip;
  Boolean         hasRefGeneTracking = FALSE;
  Boolean         hasRefTrackStatus;
  Boolean         hasLegalStatus;
  Int2            accn_count = 0;
  Int2            gi_count = 0;
  Int4            runsofn;
  Int4            segnum;
  StreamCache     sc;
  RunOfNs         ron;
  Boolean         leadingX;
  Boolean         isLower;
  Boolean         isFirst;
  CharPtr         bases;
  Int4            dnalen;
  Int4            total;
  Boolean         has_barcode_keyword = FALSE;
  CharPtr         keyword;
  Int4            count;
  Boolean         doNotSkip;
  SeqMgrPtr       smp;
  Int4            dblink_count = 0;
  Int4            taa_count = 0;
  Int4            bs_count = 0;
  Int4            pdb_count = 0;
  Int4            sra_count = 0;
  Int4            bp_count = 0;
  Int4            unknown_count = 0;

  /* set up data structures */

  vsp = (ValidStructPtr) (gcp->userdata);
  bsp = (BioseqPtr) (gcp->thisitem);
  vsp->bsp = bsp;
  vsp->descr = NULL;
  vsp->sfp = NULL;
  vsp->bssp = (BioseqSetPtr) (gcp->parentitem);
  vsp->bsp_partial_val = 0;

  sep = vsp->sep;

  if (vsp->non_ascii_chars) {   /* non_ascii chars in AsnRead step */
    ValidErr (vsp, SEV_REJECT, ERR_GENERIC_NonAsciiAsn, "Non-ascii chars in input ASN.1 strings");
    vsp->non_ascii_chars = FALSE;       /* only do once */
  }

  if (bsp->id == NULL) {
    ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_NoIdOnBioseq, "No ids on a Bioseq");
    return;
  }

  for (sip1 = bsp->id; sip1 != NULL; sip1 = sip1->next) {
    if (sip1->choice != SEQID_LOCAL) {
      only_local = FALSE;
    }
    if (sip1->choice == SEQID_OTHER) {
      isRefSeq = TRUE;
      tsip = (TextSeqIdPtr) sip1->data.ptrvalue;
      if (tsip != NULL && tsip->accession != NULL) {
        if (StringNICmp (tsip->accession, "NT_", 3) == 0) {
          isNTorNC = TRUE;
        } else if (StringNICmp (tsip->accession, "NC_", 3) == 0) {
          isNTorNC = TRUE;
          isNC = TRUE;
        } else if (StringNICmp (tsip->accession, "NG_", 3) == 0) {
          isNG = TRUE;
        }
      }
    } else if (sip1->choice == SEQID_GI) {
      hasGi = TRUE;
    } else if (sip1->choice == SEQID_GENBANK) {
      isGB = TRUE;
    } else if (sip1->choice == SEQID_EMBL) {
      isEMBL = TRUE;
    } else if (sip1->choice == SEQID_DDBJ) {
      isDDBJ = TRUE;
    } else if (sip1->choice == SEQID_SWISSPROT) {
      isSwissProt = TRUE;
    } else if (sip1->choice == SEQID_GIBBSQ) {
      isGIBBSQ = TRUE;
    } else if (sip1->choice == SEQID_GIBBMT) {
      isGIBBMT = TRUE;
    }

    for (sip2 = sip1->next; sip2 != NULL; sip2 = sip2->next) {
      if (SeqIdComp (sip1, sip2) != SIC_DIFF) {
        SeqIdWrite (sip1, buf1, PRINTID_FASTA_SHORT, 40);
        SeqIdWrite (sip2, buf2, PRINTID_FASTA_SHORT, 40);
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_ConflictingIdsOnBioseq, "Conflicting ids on a Bioseq: (%s - %s)", buf1, buf2);
      }
    }
  }

  for (sip1 = bsp->id; sip1 != NULL; sip1 = sip1->next) {
    /* disabled for now
    ReportLongSeqId (sip1, vsp, 40);
    */
    switch (sip1->choice) {
    case SEQID_TPG:
    case SEQID_TPE:
    case SEQID_TPD:
      hist = bsp->hist;
      if (hist == NULL || hist->assembly == NULL) {
        if (ISA_na (bsp->mol) && bsp->repr != Seq_repr_seg) {
          hist_asm_missing = TRUE;
          keywords = NULL;
          vnp = GetNextDescriptorUnindexed (bsp, Seq_descr_genbank, NULL);
          if (vnp != NULL && vnp->choice == Seq_descr_genbank) {
            gbp = (GBBlockPtr) vnp->data.ptrvalue;
            if (gbp != NULL) {
              keywords = gbp->keywords;
            }
          }
          if (keywords == NULL) {
            vnp = GetNextDescriptorUnindexed (bsp, Seq_descr_embl, NULL);
            if (vnp != NULL && vnp->choice == Seq_descr_embl) {
              ebp = (EMBLBlockPtr) vnp->data.ptrvalue;
              if (ebp != NULL) {
                keywords = ebp->keywords;
              }
            }
          }
          if (keywords != NULL) {
            for (vnp = keywords; vnp != NULL; vnp = vnp->next) {
              str = (CharPtr) vnp->data.ptrvalue;
              if (StringHasNoText (str)) continue;
              if (StringICmp (str, "TPA:reassembly") == 0) {
                hist_asm_missing = FALSE;
              }
            }
          }
          if (hist_asm_missing) {
            SeqIdWrite (bsp->id, buf1, PRINTID_FASTA_SHORT, 40);
            ValidErr (vsp, SEV_INFO, ERR_SEQ_INST_HistAssemblyMissing, "TPA record %s should have Seq-hist.assembly for PRIMARY block", buf1);
          }
        }
      }
      /* continue falling through */
    case SEQID_GENBANK:
    case SEQID_EMBL:
    case SEQID_DDBJ:
      tsip = (TextSeqIdPtr) sip1->data.ptrvalue;
      if (tsip != NULL && tsip->accession != NULL) {
        numletters = 0;
        numdigits = 0;
        letterAfterDigit = FALSE;
        badIDchars = FALSE;
        for (ptr = tsip->accession, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
          if (IS_UPPER (ch)) {
            numletters++;
            if (numdigits > 0) {
              letterAfterDigit = TRUE;
            }
          } else if (IS_DIGIT (ch)) {
            numdigits++;
          } else {
            badIDchars = TRUE;
          }
        }
        if (letterAfterDigit || badIDchars) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_BadSeqIdFormat, "Bad accession %s", tsip->accession);
        } else if (numletters == 1 && numdigits == 5 && ISA_na (bsp->mol)) {
        } else if (numletters == 2 && numdigits == 6 && ISA_na (bsp->mol)) {
        } else if (numletters == 3 && numdigits == 5 && ISA_aa (bsp->mol)) {
        } else if (numletters == 2 && numdigits == 6 && ISA_aa (bsp->mol) && bsp->repr == Seq_repr_seg) {
        } else if (numletters == 4 && numdigits == 8 && ISA_na (bsp->mol) &&
                   (sip1->choice == SEQID_GENBANK || sip1->choice == SEQID_EMBL || sip1->choice == SEQID_DDBJ ||
                    sip1->choice == SEQID_TPG || sip1->choice == SEQID_TPE || sip1->choice == SEQID_TPD)) {
        } else if (numletters == 4 && numdigits == 9 && ISA_na (bsp->mol) &&
                   (sip1->choice == SEQID_GENBANK || sip1->choice == SEQID_EMBL || sip1->choice == SEQID_DDBJ ||
                    sip1->choice == SEQID_TPG || sip1->choice == SEQID_TPE || sip1->choice == SEQID_TPD)) {
        } else if (numletters == 5 && numdigits == 7 && ISA_na (bsp->mol) &&
                   (sip1->choice == SEQID_GENBANK || sip1->choice == SEQID_EMBL || sip1->choice == SEQID_DDBJ)) {
        } else {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_BadSeqIdFormat, "Bad accession %s", tsip->accession);
        }
        if (vsp->useSeqMgrIndexes) {
          vnp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_genbank, &context);
          if (vnp != NULL) {
            gbp = (GBBlockPtr) vnp->data.ptrvalue;
            if (gbp != NULL) {
              LookForSecondaryConflict (vsp, gcp, tsip->accession, gbp->extra_accessions);
            }
          }
          vnp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_embl, &context);
          if (vnp != NULL) {
            ebp = (EMBLBlockPtr) vnp->data.ptrvalue;
            if (ebp != NULL) {
              LookForSecondaryConflict (vsp, gcp, tsip->accession, ebp->extra_acc);
            }
          }
        }
        if (hasGi) {
          if (tsip->version == 0) {
            ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_BadSeqIdFormat, "Accession %s has 0 version", tsip->accession);
          }
        }
      }
      /* and keep going with further test */
    case SEQID_OTHER:
      tsip = (TextSeqIdPtr) sip1->data.ptrvalue;
      if (tsip != NULL && tsip->name != NULL) {
        multitoken = FALSE;
        for (ptr = tsip->name, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
          if (IS_WHITESP (ch)) {
            multitoken = TRUE;
          }
        }
        if (multitoken) {
          ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_SeqIdNameHasSpace, "Seq-id.name '%s' should be a single word without any spaces", tsip->name);
        }
      }
      if (tsip != NULL && tsip->accession != NULL && sip1->choice == SEQID_OTHER) {
        numletters = 0;
        numdigits = 0;
        numunderscores = 0;
        letterAfterDigit = FALSE;
        badIDchars = FALSE;
        ptr = tsip->accession;
        isNZ = (Boolean) (StringNCmp (ptr, "NZ_", 3) == 0);
        if (isNZ) {
          ptr += 3;
        }
        for (ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
          if (IS_UPPER (ch)) {
            numletters++;
            if (numdigits > 0 || numunderscores > 0) {
              letterAfterDigit = TRUE;
            }
          } else if (IS_DIGIT (ch)) {
            numdigits++;
          } else if (ch == '_') {
            numunderscores++;
            if (numdigits > 0 || numunderscores > 1) {
              letterAfterDigit = TRUE;
            }
          } else {
            badIDchars = TRUE;
          }
        }
        if (letterAfterDigit || badIDchars) {
          ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_BadSeqIdFormat, "Bad accession %s", tsip->accession);
        } else if (isNZ && numletters == 4 && (numdigits == 8 || numdigits == 9) && numunderscores == 0) {
        } else if (isNZ && ValidateAccn (tsip->accession) == 0) {
        } else if (numletters == 2 && numdigits == 6 && numunderscores == 1) {
        } else if (numletters == 2 && numdigits == 8 && numunderscores == 1) {
        } else if (numletters == 2 && numdigits == 9 && numunderscores == 1) {
        } else {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_BadSeqIdFormat, "Bad accession %s", tsip->accession);
        }
      }
      if (hasGi && tsip != NULL && tsip->accession == NULL && (! StringHasNoText (tsip->name))) {
        if (sip1->choice == SEQID_DDBJ && bsp->repr == Seq_repr_seg) {
          sev = SEV_WARNING;
          /*
          ValidErr (vsp, sev, ERR_SEQ_INST_BadSeqIdFormat, "Missing accession for %s", tsip->name);
          */
        } else {
          sev = SEV_REJECT;
          ValidErr (vsp, sev, ERR_SEQ_INST_BadSeqIdFormat, "Missing accession for %s", tsip->name);
        }
      }
      /* and keep going with additional test */
    case SEQID_PIR:
    case SEQID_SWISSPROT:
    case SEQID_PRF:
      tsip = (TextSeqIdPtr) sip1->data.ptrvalue;
      if (tsip != NULL && StringHasNoText (tsip->accession) && ISA_na (bsp->mol)) {
        if (bsp->repr != Seq_repr_seg || hasGi) {
          if (sip1->choice != SEQID_DDBJ || bsp->repr != Seq_repr_seg) {
            SeqIdWrite (bsp->id, buf1, PRINTID_FASTA_LONG, sizeof (buf1) - 1);
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_BadSeqIdFormat, "Missing accession for %s", buf1);
          }
        }
      }
      if (tsip != NULL && StringHasNoText (tsip->accession) &&
          StringHasNoText (tsip->name) && ISA_aa (bsp->mol)) {
        if (sip1->choice == SEQID_PIR || sip1->choice == SEQID_SWISSPROT || sip1->choice == SEQID_PRF) {
          SeqIdWrite (bsp->id, buf1, PRINTID_FASTA_LONG, sizeof (buf1) - 1);
          ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_BadSeqIdFormat, "Missing identifier for %s", buf1);
        }
      }
      accn_count++;
      break;
    case SEQID_GPIPE:
      break;
    case SEQID_PATENT:
      isPatent = TRUE;
      break;
    case SEQID_PDB:
      isPDB = TRUE;
      break;
    case SEQID_GI:
      if (sip1->data.intvalue <= 0) {
        ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_ZeroGiNumber, "Invalid GI number");
      }
      gi_count++;
      break;
    case SEQID_GENERAL:
      dbt = (DbtagPtr) sip1->data.ptrvalue;
      if (dbt != NULL) {
        if (StringICmp (dbt->db, "LRG") == 0) {
          isLRG = TRUE;
        }
        if (StringLen (dbt->db) > 20) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_BadSeqIdFormat, "Database name longer than 20 characters");
        }
      }
      break;
    default:
      break;
    }
  }

  if (isLRG) {
    if (! isNG) {
      ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_ConflictingIdsOnBioseq, "LRG sequence needs NG_ accession");
    }
  }

  if (gi_count > 0 && accn_count == 0 && (! isPDB) && bsp->repr != Seq_repr_virtual) {
    if (vsp->seqSubmitParent) {
      sev = SEV_WARNING;
    } else {
      sev = SEV_ERROR;
    }
    ValidErr (vsp, sev, ERR_SEQ_INST_GiWithoutAccession, "No accession on sequence with gi number");
  }
  if (gi_count > 0 && accn_count > 1) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_MultipleAccessions, "Multiple accessions on sequence with gi number");
  }

  /* optionally check IDs against older version in database */

  if (vsp->validateIDSet) {
    ValidateIDSetAgainstDb (gcp, vsp, bsp);
  }

  vnp = NULL;
  if (vsp->useSeqMgrIndexes) {
    vnp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &context);
  } else {
    bcp = BioseqContextNew (bsp);
    vnp = BioseqContextGetSeqDescr (bcp, Seq_descr_molinfo, NULL, NULL);
    BioseqContextFree (bcp);
  }
  if (vnp != NULL) {
    mip = (MolInfoPtr) vnp->data.ptrvalue;
  }

  if (vsp->useSeqMgrIndexes) {
    vnp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &context);
    while (vnp != NULL) {
      uop = (UserObjectPtr) vnp->data.ptrvalue;
      if (uop != NULL) {
        oip = uop->type;
        if (oip != NULL && StringICmp (oip->str, "TpaAssembly") == 0) {
          if (! IsTpa (bsp, TRUE, &isRefSeq)) {
            olditemid = gcp->itemID;
            olditemtype = gcp->thistype;
            gcp->itemID = context.itemID;
            gcp->thistype = OBJ_SEQDESC;
            SeqIdWrite (bsp->id, buf1, PRINTID_FASTA_SHORT, 40);
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_InvalidForType, "Non-TPA record %s should not have TpaAssembly object", buf1);
            gcp->itemID = olditemid;
            gcp->thistype = olditemtype;
          }
        } else if (oip != NULL && StringICmp (oip->str, "RefGeneTracking") == 0) {
          hasRefGeneTracking = TRUE;
          hasRefTrackStatus = FALSE;
          hasLegalStatus = FALSE;
          for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
            oip = ufp->label;
            if (oip != NULL && StringCmp (oip->str, "Status") == 0) {
              hasRefTrackStatus = TRUE;
              str = (CharPtr) ufp->data.ptrvalue;
              if (StringHasNoText (str)) {
                str = "?";
              }
              for (i = 0; legal_refgene_status_strings [i] != NULL; i++) {
                if (StringICmp (str, legal_refgene_status_strings [i]) == 0) {
                  hasLegalStatus = TRUE;
                  break;
                }
              }
              if (! hasLegalStatus) {
                olditemid = gcp->itemID;
                olditemtype = gcp->thistype;
                gcp->itemID = context.itemID;
                gcp->thistype = OBJ_SEQDESC;
                vsp->descr = vnp;
                ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_RefGeneTrackingIllegalStatus, "RefGeneTracking object has illegal Status '%s'", str);
                vsp->descr = NULL;
                gcp->itemID = olditemid;
                gcp->thistype = olditemtype;
              }
            }
          }
          if (! hasRefTrackStatus) {
            olditemid = gcp->itemID;
            olditemtype = gcp->thistype;
            gcp->itemID = context.itemID;
            gcp->thistype = OBJ_SEQDESC;
            vsp->descr = vnp;
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_RefGeneTrackingWithoutStatus, "RefGeneTracking object needs to have Status set");
            vsp->descr = NULL;
            gcp->itemID = olditemid;
            gcp->thistype = olditemtype;
          }
          if (! isRefSeq && ! vsp->is_refseq_in_sep) {
            olditemid = gcp->itemID;
            olditemtype = gcp->thistype;
            gcp->itemID = context.itemID;
            gcp->thistype = OBJ_SEQDESC;
            vsp->descr = vnp;
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_RefGeneTrackingOnNonRefSeq, "RefGeneTracking object should only be in RefSeq record");
            vsp->descr = NULL;
            gcp->itemID = olditemid;
            gcp->thistype = olditemtype;
          }
        } else if (oip != NULL && StringICmp (oip->str, "DBLink") == 0) {
          dblink_count++;
          for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
            oip = ufp->label;
            if (oip == NULL || oip->str == NULL) continue;
            if (StringICmp (oip->str, "Trace Assembly Archive") == 0 && ufp->choice == 8) {
              taa_count++;
            } else if (StringICmp (oip->str, "BioSample") == 0 && ufp->choice == 7) {
              bs_count++;
            } else if (StringICmp (oip->str, "ProbeDB") == 0 && ufp->choice == 7) {
              pdb_count++;
            } else if (StringICmp (oip->str, "Sequence Read Archive") == 0 && ufp->choice == 7) {
              sra_count++;
            } else if (StringICmp (oip->str, "BioProject") == 0 && ufp->choice == 7) {
              bp_count++;
            } else {
              unknown_count++;
            }
          }
        } else if (oip != NULL && StringICmp (oip->str, "StructuredComment") == 0) {
          is_genome_assembly = FALSE;
          is_finished_status = FALSE;
          for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
            oip = ufp->label;
            if (oip == NULL || oip->str == NULL) continue;
            if (StringICmp (oip->str, "StructuredCommentPrefix") == 0) {
              if (StringCmp ((CharPtr) ufp->data.ptrvalue, "##Genome-Assembly-Data-START##") == 0) {
                is_genome_assembly = TRUE;
              }
             } else if (StringICmp (oip->str, "Current Finishing Status") == 0) {
               if (StringCmp ((CharPtr) ufp->data.ptrvalue, "Finished") == 0) {
                 is_finished_status = TRUE;
               }
             }
          }
          if (is_genome_assembly && is_finished_status && mip != NULL && mip->tech == MI_TECH_wgs) {
            olditemid = gcp->itemID;
            olditemtype = gcp->thistype;
            gcp->itemID = context.itemID;
            gcp->thistype = OBJ_SEQDESC;
            SeqIdWrite (bsp->id, buf1, PRINTID_FASTA_SHORT, 40);
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_FinishedStatusForWGS, "WGS record %s should not have Finished status", buf1);
            gcp->itemID = olditemid;
            gcp->thistype = olditemtype;
          }
        }
        if ((keyword = KeywordForStructuredCommentName (uop)) != NULL) {
          if (IsStructuredCommentValid(uop, NULL, NULL) == eFieldValid_Valid) {
            if (HasKeywordForStructuredCommentName(bsp, uop)) {
              /* as it should be */
            } else {
              ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_MissingKeyword, "Structured Comment compliant, keyword should be added");
            }
          } else {
            if (HasKeywordForStructuredCommentName(bsp, uop)) {
              ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_BadKeyword, "Structured Comment is non-compliant, keyword should be removed");
            }
          }
        }
      }
      vnp = SeqMgrGetNextDescriptor (bsp, vnp, Seq_descr_user, &context);
    }
  }

  if (dblink_count > 1) {
    ValidErr (vsp, SEV_REJECT, ERR_SEQ_DESCR_DBLinkProblem, "%ld DBLink user objects apply to a Bioseq", (long) dblink_count);
  }
  if (taa_count > 1) {
    ValidErr (vsp, SEV_REJECT, ERR_SEQ_DESCR_DBLinkProblem, "DBLink user object has %ld Trace Assembly Archive entries", (long) taa_count);
  }
  if (bs_count > 1) {
    ValidErr (vsp, SEV_REJECT, ERR_SEQ_DESCR_DBLinkProblem, "DBLink user object has %ld BioSample entries", (long) bs_count);
  }
  if (pdb_count > 1) {
    ValidErr (vsp, SEV_REJECT, ERR_SEQ_DESCR_DBLinkProblem, "DBLink user object has %ld ProbeDB entries", (long) pdb_count);
  }
  if (sra_count > 1) {
    ValidErr (vsp, SEV_REJECT, ERR_SEQ_DESCR_DBLinkProblem, "DBLink user object has %ld Sequence Read Archive entries", (long) sra_count);
  }
  if (bp_count > 1) {
    ValidErr (vsp, SEV_REJECT, ERR_SEQ_DESCR_DBLinkProblem, "DBLink user object has %ld BioProject entries", (long) bp_count);
  }
  if (unknown_count > 0) {
    ValidErr (vsp, SEV_REJECT, ERR_SEQ_DESCR_DBLinkProblem, "DBLink user object has %ld unrecognized entries", (long) unknown_count);
  }

  for (sip1 = bsp->id; sip1 != NULL; sip1 = sip1->next) {
    bsp2 = BioseqFindSpecial (sip1);
    if (bsp2 == NULL) {
      if (!isPatent) {
        SeqIdWrite (sip1, buf1, PRINTID_FASTA_SHORT, 40);
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_IdOnMultipleBioseqs, "BioseqFind (%s) unable to find itself - possible internal error", buf1);
      }
    } else if (bsp2 != bsp) {
      if (sip1->choice == SEQID_GENERAL) {
        dbt = (DbtagPtr) sip1->data.ptrvalue;
        if (dbt != NULL && StringICmp (dbt->db, "NCBIFILE") == 0) continue;
      }
      SeqIdWrite (sip1, buf1, PRINTID_FASTA_SHORT, 40);
      ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_IdOnMultipleBioseqs, "SeqID %s is present on multiple Bioseqs in record", buf1);
    }
  }

  for (i = 0; i < 4; i++)
    errors[i] = FALSE;

  switch (bsp->repr) {
  case Seq_repr_virtual:
    if ((bsp->seq_ext_type) || (bsp->seq_ext != NULL))
      errors[0] = TRUE;
    if ((bsp->seq_data_type) || (bsp->seq_data != NULL))
      errors[3] = TRUE;
    break;
  case Seq_repr_map:
    if ((bsp->seq_ext_type != 3) || (bsp->seq_ext == NULL))
      errors[1] = TRUE;
    if ((bsp->seq_data_type) || (bsp->seq_data != NULL))
      errors[3] = TRUE;
    break;
  case Seq_repr_ref:
    if ((bsp->seq_ext_type != 2) || (bsp->seq_ext == NULL))
      errors[1] = TRUE;
    if ((bsp->seq_data_type) || (bsp->seq_data != NULL))
      errors[3] = TRUE;
    break;
  case Seq_repr_seg:
    if ((bsp->seq_ext_type != 1) || (bsp->seq_ext == NULL))
      errors[1] = TRUE;
    if ((bsp->seq_data_type) || (bsp->seq_data != NULL))
      errors[3] = TRUE;
    break;
  case Seq_repr_raw:
  case Seq_repr_const:
    if ((bsp->seq_ext_type) || (bsp->seq_ext != NULL))
      errors[0] = TRUE;
    if ((bsp->seq_data_type < 1) || (bsp->seq_data_type > 11)
        || (bsp->seq_data == NULL))
      errors[2] = TRUE;
    break;
  case Seq_repr_delta:
    if ((bsp->seq_ext_type != 4) || (bsp->seq_ext == NULL))
      errors[1] = TRUE;
    if ((bsp->seq_data_type) || (bsp->seq_data != NULL))
      errors[3] = TRUE;
    break;
  default:
    ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_ReprInvalid, "Invalid Bioseq->repr = %d", (int) (bsp->repr));
    return;
  }

  if (errors[0] == TRUE) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_ExtNotAllowed, "Bioseq-ext not allowed on %s Bioseq", repr[bsp->repr - 1]);
    retval = FALSE;
  }

  if (errors[1] == TRUE) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_ExtBadOrMissing, "Missing or incorrect Bioseq-ext on %s Bioseq", repr[bsp->repr - 1]);
    retval = FALSE;
  }

  if (errors[2] == TRUE) {
    ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_SeqDataNotFound, "Missing Seq-data on %s Bioseq", repr[bsp->repr - 1]);
    retval = FALSE;
  }

  if (errors[3] == TRUE) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_SeqDataNotAllowed, "Seq-data not allowed on %s Bioseq", repr[bsp->repr - 1]);
    retval = FALSE;
  }

  if (!retval)
    return;

  oldEntityID = gcp->entityID;
  oldItemID = gcp->itemID;
  oldItemtype = gcp->thistype;

  if (ISA_aa (bsp->mol)) {
    if (bsp->topology > 1) {    /* not linear */
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_CircularProtein, "Non-linear topology set on protein");
    }
    if (bsp->strand > 1) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_DSProtein, "Protein not single stranded");
    }

  } else {
    if (!bsp->mol)
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_MolNotSet, "Bioseq.mol is 0");
    else if (bsp->mol == Seq_mol_other)
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_MolOther, "Bioseq.mol is type other");
    else if (bsp->mol == Seq_mol_na)
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_MolNuclAcid, "Bioseq.mol is type na");
  }

  if (ISA_na (bsp->mol)) {
    if (bsp->strand > 1 && mip != NULL) {
      if (mip->biomol == MOLECULE_TYPE_MRNA) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_DSmRNA, "mRNA not single stranded");
      }
    }
  }

  gcp->entityID = oldEntityID;
  gcp->itemID = oldItemID;
  gcp->thistype = oldItemtype;

  /* check sequence alphabet */
  if (((bsp->repr == Seq_repr_raw) || (bsp->repr == Seq_repr_const)) && bsp->seq_data_type != Seq_code_gap) {
    if (bsp->fuzz != NULL) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_FuzzyLen, "Fuzzy length on %s Bioseq", repr[bsp->repr - 1]);
    }

    if (bsp->length < 1) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_InvalidLen, "Invalid Bioseq length [%ld]", (long) bsp->length);
    }

    seqtype = (int) (bsp->seq_data_type);
    switch (seqtype) {
    case Seq_code_iupacna:
    case Seq_code_ncbi2na:
    case Seq_code_ncbi4na:
    case Seq_code_ncbi8na:
    case Seq_code_ncbipna:
      if (ISA_aa (bsp->mol)) {
        ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_InvalidAlphabet, "Using a nucleic acid alphabet on a protein sequence");
        return;
      }
      break;
    case Seq_code_iupacaa:
    case Seq_code_ncbi8aa:
    case Seq_code_ncbieaa:
    case Seq_code_ncbipaa:
    case Seq_code_ncbistdaa:
      if (ISA_na (bsp->mol)) {
        ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_InvalidAlphabet, "Using a protein alphabet on a nucleic acid");
        return;
      }
      break;
    case Seq_code_gap:
      break;
    default:
      ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_InvalidAlphabet, "Using illegal sequence alphabet [%d]", (int) bsp->seq_data_type);
      return;
    }

    check_alphabet = FALSE;
    switch (seqtype) {
    case Seq_code_iupacaa:
    case Seq_code_iupacna:
    case Seq_code_ncbieaa:
    case Seq_code_ncbistdaa:
      check_alphabet = TRUE;

    case Seq_code_ncbi8na:
    case Seq_code_ncbi8aa:
      divisor = 1;
      break;

    case Seq_code_ncbi4na:
      divisor = 2;
      break;

    case Seq_code_ncbi2na:
      divisor = 4;
      break;

    case Seq_code_ncbipna:
      divisor = 5;
      break;

    case Seq_code_ncbipaa:
      divisor = 21;
      break;
    }

    len = bsp->length;
    if (len % divisor)
      len += divisor;
    len /= divisor;
    len2 = BSLen ((ByteStorePtr) bsp->seq_data);
    if (len > len2) {
      ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_SeqDataLenWrong, "Bioseq.seq_data too short [%ld] for given length [%ld]", (long) (len2 * divisor),
                (long) bsp->length);
      return;
    } else if (len < len2) {
      ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_SeqDataLenWrong, "Bioseq.seq_data is larger [%ld] than given length [%ld]", (long) (len2 * divisor),
                (long) bsp->length);
    }

    if (check_alphabet) {       /* check 1 letter alphabets */
      switch (seqtype) {
      case Seq_code_iupacaa:
      case Seq_code_ncbieaa:
        termination = '*';
        gapchar = '-';
        break;
      case Seq_code_ncbistdaa:
        termination = 25;
        gapchar = 0;
        break;
      default:
        termination = '\0';
        gapchar = '\0';
        break;
      }
      if (! StreamCacheSetup (bsp, NULL, STREAM_EXPAND_GAPS, &sc)) {
        ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_SeqPortFail, "Can't open StreamCache");
        return;
      }
      /*
      spp = SeqPortNew (bsp, 0, -1, 0, 0);
      if (spp == NULL) {
        ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_SeqPortFail, "Can't open SeqPort");
        return;
      }
      */
      i = 0;
      terminations = 0;
      dashes = 0;
      gapatstart = FALSE;
      trailingX = 0;
      leadingX = FALSE;
      isLower = FALSE;
      isFirst = TRUE;
      for (len = 0; len < bsp->length; len++) {
        residue = StreamCacheGetResidue (&sc);
        /*
        residue = SeqPortGetResidue (spp);
        */
        if (!IS_residue (residue)) {
          i++;
          if (i > 10) {
            ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_InvalidResidue, "More than 10 invalid residues. Checking stopped");
            /*
            SeqPortFree (spp);
            */
            if (vsp->patch_seq)
              PatchBadSequence (bsp);
            return;
          } else {
            BSSeek ((ByteStorePtr) bsp->seq_data, len, SEEK_SET);
            x = BSGetByte ((ByteStorePtr) bsp->seq_data);
            if (bsp->seq_data_type == Seq_code_ncbistdaa) {
              ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_InvalidResidue, "Invalid residue [%d] at position [%ld]", (int) x, (long) (len + 1));
            } else if (IS_ALPHA (x)) {
              ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_InvalidResidue, "Invalid residue '%c' at position [%ld]", (char) x, (long) (len + 1));
            } else {
              ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_InvalidResidue, "Invalid residue [%d] at position [%ld]", (int) x, (long) (len + 1));
            }
          }
        } else if (residue == termination) {
          terminations++;
          trailingX = 0;        /* suppress if followed by terminator */
        } else if (residue == gapchar) {
          dashes++;
          if (len == 0) {
            gapatstart = TRUE;
          }
        } else if (residue == 'X') {
          if (ISA_na (bsp->mol)) {
            ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_InvalidResidue, "Invalid nucleotide residue '%c' at position [%ld]", (char) residue, (long) (len + 1));
          } else {
            trailingX++;
            if (isFirst) {
              leadingX = TRUE;
            }
          }
        } else if (ISA_na (bsp->mol) && StringChr ("EFIJLOPQZ", (Char) residue) != NULL) {
          ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_InvalidResidue, "Invalid nucleotide residue '%c' at position [%ld]", (char) residue, (long) (len + 1));
        } else if (! IS_ALPHA ((Char) residue)) {
          ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_InvalidResidue, "Invalid residue '%c' at position [%ld]", (char) residue, (long) (len + 1));
        } else {
          trailingX = 0;
          if (IS_LOWER ((Char) residue)) {
            isLower = TRUE;
          }
        }
        isFirst = FALSE;
      }
      /*
      SeqPortFree (spp);
      */
      if (ISA_aa (bsp->mol) && (leadingX || trailingX > 0)) {
        /* only show leading or trailing X if product of NNN in nucleotide */
        cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
        if (cds != NULL) {
          crp = (CdRegionPtr) cds->data.value.ptrvalue;
          if (crp != NULL) {
            dnalen = SeqLocLen (cds->location);
            if (dnalen > 5) {
              bases = ReadCodingRegionBases (cds->location, dnalen, crp->frame, &total);
              len = StringLen (bases);
              if (len > 5) {
                if (StringNICmp (bases, "NNN", 3) != 0) {
                  leadingX = FALSE;
                }
                if (StringNICmp (bases + len - 3, "NNN", 3) != 0) {
                  trailingX = 0;
                }
              }
              MemFree (bases);
            }
          }
        }
      }
      if (leadingX) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_LeadingX, "Sequence starts with leading X", (int) leadingX);
      }
      if (trailingX > 0 && SuppressTrailingXMessage (bsp)) {
        /* suppress if cds translation ends in '*' or 3' partial */
      } else if (trailingX > 1) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_TrailingX, "Sequence ends in %d trailing Xs", (int) trailingX);
      } else if (trailingX > 0) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_TrailingX, "Sequence ends in %d trailing X", (int) trailingX);
      }
      if (isLower) {
        ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_InvalidResidue, "Sequence contains lower-case characters");
      }
      if (terminations > 0 || dashes > 0) {
        cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
        grp = SeqMgrGetGeneXref (cds);
        genelbl = NULL;
        if (grp == NULL && cds != NULL) {
          gene = SeqMgrGetOverlappingGene (cds->location, &genectxt);
          if (gene != NULL) {
            grp = (GeneRefPtr) gene->data.value.ptrvalue;
          }
        }
        if (grp != NULL && (!SeqMgrGeneIsSuppressed (grp))) {
          if (grp->locus != NULL)
            genelbl = (grp->locus);
          else if (grp->locus_tag != NULL)
            genelbl = (grp->locus_tag);
          else if (grp->desc != NULL)
            genelbl = (grp->desc);
          else if (grp->syn != NULL)
            genelbl = (CharPtr) (grp->syn->data.ptrvalue);
        }
        prot = SeqMgrGetBestProteinFeature (bsp, &protctxt);
        protlbl = protctxt.label;
      }
      if (StringHasNoText (genelbl)) {
        genelbl = "gene?";
      }
      if (StringHasNoText (protlbl)) {
        protlbl = "prot?";
      }
      if (dashes > 0) {
        if (gapatstart && dashes == 1) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_BadProteinStart, "gap symbol at start of protein sequence (%s - %s)", genelbl, protlbl);
        } else if (gapatstart) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_BadProteinStart, "gap symbol at start of protein sequence (%s - %s)", genelbl, protlbl);
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_GapInProtein, "[%d] internal gap symbols in protein sequence (%s - %s)", (dashes - 1), genelbl, protlbl);
        } else {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_GapInProtein, "[%d] internal gap symbols in protein sequence (%s - %s)", dashes, genelbl, protlbl);
        }
      }
      if (terminations) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_StopInProtein, "[%d] termination symbols in protein sequence (%s - %s)", terminations, genelbl, protlbl);
        if (!i)
          return;
      }
      if (i) {
        if (vsp->patch_seq)
          PatchBadSequence (bsp);
        return;
      }

    }
  }

  if (ISA_na (bsp->mol) && bsp->repr == Seq_repr_delta && DeltaLitOnly (bsp)) {
    if (! StreamCacheSetup (bsp, NULL, EXPAND_GAPS_TO_DASHES, &sc)) {
      ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_SeqPortFail, "Can't open StreamCache");
      return;
    }
    in_gap = FALSE;
    in_N = FALSE;
    for (len = 0; len < bsp->length; len++) {
      residue = StreamCacheGetResidue (&sc);
      if (residue == '-') {
        if (in_N) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_InternalNsAdjacentToGap,
                    "Ambiguous residue N is adjacent to a gap around position %ld",
                    (long) len + 1);
        }
        in_N = FALSE;
        in_gap = TRUE;
      } else if (residue == 'N') {
        if (in_gap) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_InternalNsAdjacentToGap,
                    "Ambiguous residue N is adjacent to a gap around position %ld",
                    (long) len + 1);
        }
        in_gap = FALSE;
        in_N = TRUE;
      } else {
        in_gap = FALSE;
        in_N = FALSE;
      }
    }
  }

  if ((bsp->repr == Seq_repr_seg) || (bsp->repr == Seq_repr_ref)) {     /* check segmented sequence */
    head.choice = SEQLOC_MIX;
    head.data.ptrvalue = bsp->seq_ext;
    head.next = NULL;
    ValidateSeqLoc (vsp, (SeqLocPtr) & head, "Segmented Bioseq");
    /* check the length */
    len = 0;
    vnp = NULL;
    while ((vnp = SeqLocFindNext (&head, vnp)) != NULL) {
      len2 = SeqLocLen (vnp);
      if (len2 > 0)
        len += len2;
    }
    if (bsp->length > len) {
      ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_SeqDataLenWrong, "Bioseq.seq_data too short [%ld] for given length [%ld]", (long) (len), (long) bsp->length);
    } else if (bsp->length < len) {
      ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_SeqDataLenWrong, "Bioseq.seq_data is larger [%ld] than given length [%ld]", (long) (len), (long) bsp->length);
    }

    vnp = NULL;
    idlist = NULL;
    while ((vnp = SeqLocFindNext (&head, vnp)) != NULL) {
      sip1 = SeqLocId (vnp);
      if (sip1 != NULL) {
        SeqIdWrite (sip1, buf1, PRINTID_FASTA_SHORT, 40);
        ValNodeCopyStr (&idlist, vnp->choice, buf1);
      }
    }
    if (idlist != NULL) {
      idlist = ValNodeSort (idlist, SortVnpByString);
      last = (CharPtr) idlist->data.ptrvalue;
      lastchoice = (Uint1) idlist->choice;
      vnp = idlist->next;
      while (vnp != NULL) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringICmp (last, str) == 0) {
          if (vnp->choice == lastchoice && lastchoice == SEQLOC_WHOLE) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_DuplicateSegmentReferences, "Segmented sequence has multiple references to %s", str);
          } else {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_DuplicateSegmentReferences,
                      "Segmented sequence has multiple references to %s that are not SEQLOC_WHOLE", str);
          }
        } else {
          last = (CharPtr) vnp->data.ptrvalue;
          lastchoice = (Uint1) vnp->choice;
        }
        vnp = vnp->next;
      }
      ValNodeFreeData (idlist);
    }

    vsp->bsp_partial_val = SeqLocPartialCheck ((SeqLocPtr) (&head));
    if (ISA_aa (bsp->mol)) {
      got_partial = FALSE;
      if (mip != NULL) {
          switch (mip->completeness) {
          case 2:             /* partial */
            got_partial = TRUE;
            if (!vsp->bsp_partial_val) {
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_PartialInconsistent, "Complete segmented sequence with MolInfo partial");
            }
            break;
          case 3:             /* no-left */
            if (!(vsp->bsp_partial_val & SLP_START) || (vsp->bsp_partial_val && SLP_STOP))
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_PartialInconsistent, "No-left inconsistent with segmented SeqLoc");
            got_partial = TRUE;
            break;
          case 4:             /* no-right */
            if (!(vsp->bsp_partial_val & SLP_STOP) || (vsp->bsp_partial_val && SLP_START))
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_PartialInconsistent, "No-right inconsistent with segmented SeqLoc");
            got_partial = TRUE;
            break;
          case 5:             /* no-ends */
            if ((!(vsp->bsp_partial_val & SLP_STOP)) || (!(vsp->bsp_partial_val & SLP_START)))
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_PartialInconsistent, "No-ends inconsistent with segmented SeqLoc");
            got_partial = TRUE;
            break;
          default:
            break;
          }
      }
      if (!got_partial && vsp->bsp_partial_val) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_PartialInconsistent, "Partial segmented sequence without MolInfo partial");
      }
    }
  }

  if (bsp->repr == Seq_repr_delta || bsp->repr == Seq_repr_raw) {

    vnp = NULL;
    if (vsp->useSeqMgrIndexes) {
      sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_genbank, &context);
    } else {
      bcp = BioseqContextNew (bsp);
      sdp = BioseqContextGetSeqDescr (bcp, Seq_descr_genbank, NULL, NULL);
      BioseqContextFree (bcp);
    }
    if (sdp != NULL) {
      gbp = (GBBlockPtr) sdp->data.ptrvalue;
      if (gbp != NULL) {
        for (vnp = gbp->keywords; vnp != NULL; vnp = vnp->next) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (StringICmp (str, "HTGS_ACTIVEFIN") == 0) {
            isActiveFin = TRUE;
          } else if (StringICmp (str, "HTGS_DRAFT") == 0) {
            isDraft = TRUE;
          } else if (StringICmp (str, "HTGS_FULLTOP") == 0) {
            isFullTop = TRUE;
          } else if (StringICmp (str, "HTGS_PREFIN") == 0) {
            isPreFin = TRUE;
          } else if (StringICmp (str, "BARCODE") == 0) {
            has_barcode_keyword = TRUE;
            if (/* ! vsp->is_barcode_sep */ mip == NULL || mip->tech != MI_TECH_barcode) {
              olditemid = gcp->itemID;
              olditemtype = gcp->thistype;
              if (sdp->extended != 0) {
                ovp = (ObjValNodePtr) sdp;
                gcp->itemID = ovp->idx.itemID;
                gcp->thistype = OBJ_SEQDESC;
              }
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadKeyword, "BARCODE keyword without Molinfo.tech barcode");
              gcp->itemID = olditemid;
              gcp->thistype = olditemtype;
            }
          }
        }
      }
    }
    if (mip != NULL && mip->tech == MI_TECH_barcode && !has_barcode_keyword) {
      ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_BadKeyword, "Molinfo.tech barcode without BARCODE keyword");
    }
  }


  if (bsp->repr == Seq_repr_delta) {
    len = 0;
    count = 0;
    doNotSkip = TRUE;
    for (vnp = (ValNodePtr) (bsp->seq_ext); vnp != NULL; vnp = vnp->next) {
      if (vnp->choice == 1) {
        count++;
      }
    }
    if (count > 10000) {
      smp = SeqMgrGet ();
      if (smp != NULL) {
        if (smp->seq_len_lookup_func == NULL) {
          doNotSkip = FALSE;
        }
      }
    }
    if (doNotSkip) {
      for (vnp = (ValNodePtr) (bsp->seq_ext), segnum = 1; vnp != NULL; vnp = vnp->next, segnum++) {
        if (vnp->data.ptrvalue == NULL)
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_SeqDataLenWrong, "NULL pointer in delta seq_ext valnode");
        else {
          switch (vnp->choice) {
          case 1:                /* SeqLocPtr */
            slp = (SeqLocPtr) (vnp->data.ptrvalue);
            if (slp != NULL && slp->choice == SEQLOC_WHOLE) {
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_WholeComponent, "Delta seq component should not be of type whole");
            }
            sip3 = SeqLocId (slp);
            if (sip3 != NULL) {
              if (sip3->choice == SEQID_GI && sip3->data.intvalue <= 0) {
                ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_DeltaComponentIsGi0, "Delta component is gi|0");
              }
              for (sip1 = bsp->id; sip1 != NULL; sip1 = sip1->next) {
                if (SeqIdComp (sip1, sip3) == SIC_YES) {
                  ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_SelfReferentialSequence,
                            "Self-referential delta sequence");
                }
              }
            }
            len2 = SeqLocLen (slp);
            if (len2 < 0)
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_SeqDataLenWrong, "-1 length on seq-loc of delta seq_ext");
            else
              len += len2;
            sip3 = SeqLocId (slp);
            if (sip3 != NULL && slp != NULL && slp->choice == SEQLOC_INT) {
              sintp = (SeqIntPtr) slp->data.ptrvalue;
              if (sintp != NULL && (sip3->choice == SEQID_GI ||
                                   sip3->choice == SEQID_GENBANK ||
                                   sip3->choice == SEQID_EMBL ||
                                   sip3->choice == SEQID_DDBJ ||
                                   sip3->choice == SEQID_TPG ||
                                   sip3->choice == SEQID_TPE ||
                                   sip3->choice == SEQID_TPD ||
                                   sip3->choice == SEQID_OTHER)) {
                vn.choice = SEQLOC_WHOLE;
                vn.data.ptrvalue = sip3;
                vn.next = NULL;
                len3 = SeqLocLen (&vn);
                /* -1 signifies failure to lookup or not connected to lookup function */
                if (len3 != -1) {
                  if (sintp->to >= len3) {
                    SeqIdWrite (sip3, buf1, PRINTID_FASTA_SHORT, 40);
                    ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_SeqDataLenWrong,
                              "Seq-loc extent (%ld) greater than length of %s (%ld)",
                              (long) (sintp->to + 1), buf1, (long) len3);
                  }
                }
              }
            }
            if (len2 <= 10) {
              str = SeqLocPrint ((SeqLocPtr) (vnp->data.ptrvalue));
              if (str == NULL) {
                str = StringSave ("?");
              }
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_SeqLocLength, "Short length (%ld) on seq-loc (%s) of delta seq_ext", (long) len2, str);
              MemFree (str);
            }
            break;
          case 2:                /* SeqLitPtr */
            slitp = (SeqLitPtr) (vnp->data.ptrvalue);
            if (slitp->seq_data != NULL && slitp->seq_data_type != Seq_code_gap) {
              if (slitp->length == 0) {
                ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_SeqLitDataLength0, "Seq-lit of length 0 in delta chain");
              }
              sctp = SeqCodeTableFind (slitp->seq_data_type);
              if (sctp == NULL) {
                ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_InvalidAlphabet, "Using illegal sequence alphabet [%d] in SeqLitPtr", (int) slitp->seq_data_type);
                len += slitp->length;
                break;
              }

              start_at = (Int2) (sctp->start_at);
              num = (Int2) (sctp->num);

              switch (slitp->seq_data_type) {
              case Seq_code_iupacaa:
              case Seq_code_iupacna:
              case Seq_code_ncbieaa:
              case Seq_code_ncbistdaa:
                BSSeek ((ByteStorePtr) slitp->seq_data, 0, SEEK_SET);
                for (len2 = 1; len2 <= (slitp->length); len2++) {
                  is_invalid = FALSE;
                  residue = BSGetByte ((ByteStorePtr) slitp->seq_data);
                  i = residue - start_at;
                  if ((i < 0) || (i >= num))
                    is_invalid = TRUE;
                  else if (*(sctp->names[i]) == '\0')
                    is_invalid = TRUE;
                  if (is_invalid) {
                    if (slitp->seq_data_type == Seq_code_ncbistdaa)
                      ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_InvalidResidue, "Invalid residue [%d] at position [%ld]", (int) residue, (long) (len + len2));
                    else
                      ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_InvalidResidue, "Invalid residue [%c] at position [%ld]", (char) residue, (long) (len + len2));
                  }
                }
                break;
              default:
                break;
              }
              if (mip != NULL) {
                if (mip->tech == MI_TECH_htgs_1 || mip->tech == MI_TECH_htgs_2) {
                  runsofn = CountAdjacentNsInInterval (gcp, bsp, len, len + slitp->length - 1);
                  if (runsofn > 80) {
                    ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_InternalNsInSeqLit, "Run of %ld Ns in delta component %d that starts at base %ld", (long) runsofn, segnum, (long) (len + 1));
                  }
                } else if (mip->tech == MI_TECH_wgs) {
                  runsofn = CountAdjacentNsInInterval (gcp, bsp, len, len + slitp->length - 1);
                  if (runsofn >= 20) {
                    ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_InternalNsInSeqLit, "Run of %ld Ns in delta component %d that starts at base %ld", (long) runsofn, segnum, (long) (len + 1));
                  }
                } else if (mip->tech == MI_TECH_composite_wgs_htgs) {
                  runsofn = CountAdjacentNsInInterval (gcp, bsp, len, len + slitp->length - 1);
                  if (runsofn > 80) {
                    ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_InternalNsInSeqLit, "Run of %ld Ns in delta component %d that starts at base %ld", (long) runsofn, segnum, (long) (len + 1));
                  }
                } else {
                  runsofn = CountAdjacentNsInInterval (gcp, bsp, len, len + slitp->length - 1);
                  if (runsofn > 100) {
                    ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_InternalNsInSeqLit, "Run of %ld Ns in delta component %d that starts at base %ld", (long) runsofn, segnum, (long) (len + 1));
                  }
                }
              }
            } else if (slitp->length == 0) {
              if (isSwissProt) {
                sev = SEV_WARNING;
              } else {
                sev = SEV_ERROR;
              }
              ifp = slitp->fuzz;
              if (ifp == NULL || ifp->choice != 4 || ifp->a != 0) {
                ValidErr (vsp, sev, ERR_SEQ_INST_SeqLitGapLength0, "Gap of length 0 in delta chain");
              } else {
                ValidErr (vsp, sev, ERR_SEQ_INST_SeqLitGapLength0, "Gap of length 0 with unknown fuzz in delta chain");
              }
            } else if (slitp->seq_data == NULL || slitp->seq_data_type == Seq_code_gap) {
              if (slitp->length != 100 && slitp->fuzz != NULL) {
                ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_SeqLitGapFuzzNot100, "Gap of unknown length should have length 100");
              }
            }
            len += slitp->length;
            break;
          default:
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_ExtNotAllowed, "Illegal choice [%d] in delta chain", (int) (vnp->choice));
            break;
          }
        }
      }
      if (bsp->length > len) {
        ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_SeqDataLenWrong, "Bioseq.seq_data too short [%ld] for given length [%ld]", (long) (len), (long) bsp->length);
      } else if (bsp->length < len) {
        ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_SeqDataLenWrong, "Bioseq.seq_data is larger [%ld] than given length [%ld]", (long) (len), (long) bsp->length);
      }

    } else {

      for (vnp = (ValNodePtr) (bsp->seq_ext); vnp != NULL; vnp = vnp->next) {
        if (vnp->data.ptrvalue == NULL) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_SeqDataLenWrong, "NULL pointer in delta seq_ext valnode");
          continue;
        }
        switch (vnp->choice) {
          case 1 :
            slp = (SeqLocPtr) vnp->data.ptrvalue;
            if (slp->choice == SEQLOC_WHOLE) {
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_WholeComponent, "Delta seq component should not be of type whole");
              break;
            }
            if (slp->choice == SEQLOC_INT) {
              sintp = (SeqIntPtr) slp->data.ptrvalue;
              if (sintp != NULL) {
                len2 = sintp->to - sintp->from + 1;
                if (len2 < 0) {
                  ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_SeqDataLenWrong, "-1 length on seq-loc of delta seq_ext");
                } else {
                  len += len2;
                }
              }
            }
            break;
          case 2 :
            slitp = (SeqLitPtr) vnp->data.ptrvalue;
            len += slitp->length;
            break;
          default :
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_ExtNotAllowed, "Illegal choice [%d] in delta chain", (int) (vnp->choice));
            break;
        }
      }
      if (bsp->length > len) {
        ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_SeqDataLenWrong, "Bioseq.seq_data too short [%ld] for given length [%ld]", (long) (len), (long) bsp->length);
      } else if (bsp->length < len) {
        ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_SeqDataLenWrong, "Bioseq.seq_data is larger [%ld] than given length [%ld]", (long) (len), (long) bsp->length);
      }
    }
    if (mip != NULL) {
      is_gps = FALSE;
      sep = vsp->sep;
      if (sep != NULL && IS_Bioseq_set (sep)) {
        bssp = (BioseqSetPtr) sep->data.ptrvalue;
        if (bssp != NULL && bssp->_class == BioseqseqSet_class_gen_prod_set) {
          is_gps = TRUE;
        }
      }
      if ((!isNTorNC) && (! is_gps) && mip->tech != MI_TECH_htgs_0 && mip->tech != MI_TECH_htgs_1 &&
          mip->tech != MI_TECH_htgs_2 && mip->tech != MI_TECH_htgs_3 && mip->tech != MI_TECH_wgs &&
          mip->tech != MI_TECH_composite_wgs_htgs && mip->tech != MI_TECH_unknown && mip->tech != MI_TECH_standard
          && mip->tech != MI_TECH_htc && mip->tech != MI_TECH_barcode && mip->tech != MI_TECH_tsa) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_BadDeltaSeq, "Delta seq technique should not be [%d]", (int) (mip->tech));
      }
    }
  } else if (bsp->repr == Seq_repr_raw) {
    ron.gcp = gcp;
    ron.vsp = vsp;
    ron.ncount = 0;
    ron.maxrun = 0;
    ron.seqpos = 0;
    ron.gapcount = 0;
    ron.showAll = TRUE;
    ron.inNrun = FALSE;
    ron.isWGS = FALSE;
    if (mip == NULL) {
      vnp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &context);
      if (vnp != NULL) {
        mip = (MolInfoPtr) vnp->data.ptrvalue;
      }
    }
    if (mip != NULL && mip->tech == MI_TECH_wgs) {
      ron.isWGS = TRUE;
    }

    SeqPortStream (bsp, EXPAND_GAPS_TO_DASHES, (Pointer) &ron, CountAdjacentProc);

    /*
    if (ron.inNrun && ron.showAll && ron.ncount >= 100) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_InternalNsInSeqRaw, "Run of %ld Ns in raw sequence starting at base %ld",
                (long) ron.ncount, (long) (ron.seqpos - ron.ncount + 1));
    }
    */

    if (ron.gapcount > 0 && ISA_na (bsp->mol)) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_InternalGapsInSeqRaw, "Raw nucleotide should not contain gap characters");
    }

    /*
    if (ron.maxrun >= 100) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_InternalNsInSeqRaw, "Run of %ld Ns in raw sequence", (long) ron.maxrun);
    }
    */
  }

  if (bsp->repr == Seq_repr_delta) {
    CheckDeltaForReuse (vsp, gcp, bsp);
  }

  sev = SEV_ERROR;
  if (mip != NULL) {
    if (mip->tech != MI_TECH_htgs_0 && mip->tech != MI_TECH_htgs_1 &&
        mip->tech != MI_TECH_htgs_2 && mip->tech != MI_TECH_htgs_3) {
      sev = SEV_WARNING;
    }
  }

  if (bsp->repr == Seq_repr_delta && bsp->seq_ext_type == 4 && bsp->seq_ext != NULL) {
    vnp = (DeltaSeqPtr) bsp->seq_ext;
    if (vnp != NULL && vnp->choice == 2) {
      slitp = (SeqLitPtr) vnp->data.ptrvalue;
      if (slitp != NULL && (slitp->seq_data == NULL || slitp->seq_data_type == Seq_code_gap)) {
        ValidErr (vsp, sev, ERR_SEQ_INST_BadDeltaSeq, "First delta seq component is a gap");
      }
    }
    last_is_gap = FALSE;
    num_adjacent_gaps = 0;
    num_gaps = 0;
    non_interspersed_gaps = FALSE;
    while (vnp->next != NULL) {
      vnp = vnp->next;
      if (vnp != NULL && vnp->choice == 2) {
        slitp = (SeqLitPtr) vnp->data.ptrvalue;
        if (slitp != NULL && (slitp->seq_data == NULL || slitp->seq_data_type == Seq_code_gap)) {
          if (last_is_gap) {
            num_adjacent_gaps++;
          }
          last_is_gap = TRUE;
          num_gaps++;
        } else {
          if (! last_is_gap) {
            non_interspersed_gaps = TRUE;
          }
          last_is_gap = FALSE;
        }
      } else {
        if (! last_is_gap) {
          non_interspersed_gaps = TRUE;
        }
        last_is_gap = FALSE;
      }
    }
    if (non_interspersed_gaps && (! hasGi) && mip != NULL &&
        (mip->tech == MI_TECH_htgs_0 || mip->tech == MI_TECH_htgs_1 || mip->tech == MI_TECH_htgs_2 /* || mip->tech == MI_TECH_htgs_3 */)) {
      if (hasRefGeneTracking) {
        ValidErr (vsp, SEV_INFO, ERR_SEQ_INST_MissingGaps, "HTGS delta seq should have gaps between all sequence runs");
      } else {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_MissingGaps, "HTGS delta seq should have gaps between all sequence runs");
      }
    }
    if (vnp != NULL && vnp->choice == 2) {
      slitp = (SeqLitPtr) vnp->data.ptrvalue;
      if (slitp != NULL && (slitp->seq_data == NULL || slitp->seq_data_type == Seq_code_gap)) {
        ValidErr (vsp, sev, ERR_SEQ_INST_BadDeltaSeq, "Last delta seq component is a gap");
      }
    }
    if (num_gaps == 0 && mip != NULL) {
      if (/* mip->tech == MI_TECH_htgs_1 || */ mip->tech == MI_TECH_htgs_2) {
        if (VisitGraphsInSep (sep, NULL, NULL) == 0) {
          if (! isActiveFin) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_BadHTGSeq, "HTGS 2 delta seq has no gaps and no graphs");
          }
        }
      }
    }
    sev = SEV_ERROR;
    if (isRefSeq) {
      sev = SEV_WARNING;
    }
    if (num_adjacent_gaps > 1) {
      ValidErr (vsp, sev, ERR_SEQ_INST_BadDeltaSeq, "There are %d adjacent gaps in delta seq", (int) num_adjacent_gaps);
    } else if (num_adjacent_gaps > 0) {
      ValidErr (vsp, sev, ERR_SEQ_INST_BadDeltaSeq, "There is %d adjacent gap in delta seq", (int) num_adjacent_gaps);
    }
  }

  if (bsp->repr == Seq_repr_raw) {
    if (mip != NULL) {
      if (/* mip->tech == MI_TECH_htgs_1 || */ mip->tech == MI_TECH_htgs_2) {
        if (VisitGraphsInSep (sep, NULL, NULL) == 0) {
          if (! isActiveFin) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_BadHTGSeq, "HTGS 2 raw seq has no gaps and no graphs");
          }
        }
      }
    }
  }

  if (mip != NULL && mip->tech == MI_TECH_htgs_3) {
    if (isDraft) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_BadHTGSeq, "HTGS 3 sequence should not have HTGS_DRAFT keyword");
    }
    if (isPreFin) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_BadHTGSeq, "HTGS 3 sequence should not have HTGS_PREFIN keyword");
    }
    if (isActiveFin) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_BadHTGSeq, "HTGS 3 sequence should not have HTGS_ACTIVEFIN keyword");
    }
    if (isFullTop) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_BadHTGSeq, "HTGS 3 sequence should not have HTGS_FULLTOP keyword");
    }
  }

  if (ISA_aa (bsp->mol)) {
    if ((bsp->length <= 3) && (bsp->length >= 0) && (!isPDB)) {
      if (mip == NULL || mip->completeness < 2 || mip->completeness > 5) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_ShortSeq, "Sequence only %ld residues", (long) (bsp->length));
      }
    }

  } else {
    if ((bsp->length <= 10) && (bsp->length >= 0) && (!isPDB)) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_ShortSeq, "Sequence only %ld residues", (long) (bsp->length));
    }
  }

#if 0
  if (bsp->length > 350000 && (! isNTorNC)) {
    Boolean         isGenBankEMBLorDDBJ;
    Boolean         litHasData;
    if (bsp->repr == Seq_repr_delta) {
      isGenBankEMBLorDDBJ = FALSE;
      /* suppress this for data from genome annotation project */
      VisitBioseqsInSep (vsp->sep, (Pointer) &isGenBankEMBLorDDBJ, LookForGEDseqID);
      if (mip != NULL && isGenBankEMBLorDDBJ) {
        if (mip->tech == MI_TECH_htgs_0 || mip->tech == MI_TECH_htgs_1 || mip->tech == MI_TECH_htgs_2) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_LongHtgsSequence, "Phase 0, 1 or 2 HTGS sequence exceeds 350kbp limit");
        } else if (mip->tech == MI_TECH_htgs_3) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_SequenceExceeds350kbp, "Phase 3 HTGS sequence exceeds 350kbp limit");
        } else if (mip->tech == MI_TECH_wgs) {
          /*
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_SequenceExceeds350kbp, "WGS sequence exceeds 350kbp limit");
          */
        } else {
          len = 0;
          litHasData = FALSE;
          for (vnp = (ValNodePtr) (bsp->seq_ext); vnp != NULL; vnp = vnp->next) {
            if (vnp->choice == 2) {
              slitp = (SeqLitPtr) (vnp->data.ptrvalue);
              if (slitp != NULL) {
                if (slitp->seq_data != NULL) {
                  litHasData = TRUE;
                }
                len += slitp->length;
              }
            }
          }
          if (len > 500000 && litHasData) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_LongLiteralSequence, "Length of sequence literals exceeds 500kbp limit");
          }
        }
      }
    } else if (bsp->repr == Seq_repr_raw) {
      vnp = NULL;
      if (vsp->useSeqMgrIndexes) {
        vnp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &context);
      } else {
        bcp = BioseqContextNew (bsp);
        vnp = BioseqContextGetSeqDescr (bcp, Seq_descr_molinfo, NULL, NULL);
        BioseqContextFree (bcp);
      }
      if (vnp != NULL) {
        mip = (MolInfoPtr) vnp->data.ptrvalue;
      }
      if (mip != NULL) {
        if (mip->tech == MI_TECH_htgs_0 || mip->tech == MI_TECH_htgs_1 || mip->tech == MI_TECH_htgs_2) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_LongHtgsSequence, "Phase 0, 1 or 2 HTGS sequence exceeds 350kbp limit");
        } else if (mip->tech == MI_TECH_htgs_3) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_SequenceExceeds350kbp, "Phase 3 HTGS sequence exceeds 350kbp limit");
        } else if (mip->tech == MI_TECH_wgs) {
          /*
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_SequenceExceeds350kbp, "WGS sequence exceeds 350kbp limit");
          */
        } else {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_SequenceExceeds350kbp, "Length of sequence exceeds 350kbp limit");
        }
      } else {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_SequenceExceeds350kbp, "Length of sequence exceeds 350kbp limit");
      }
    } else {
      /* Could be a segset header bioseq that is > 350kbp */
      /* No-op for now? Or generate a warning? */
    }
  }
#endif

  if (bsp->repr == Seq_repr_seg) {
    CheckSegBspAgainstParts (vsp, gcp, bsp);
  }

  if (ISA_na (bsp->mol) || ISA_aa (bsp->mol)) {
    vnp = BioseqGetSeqDescr (bsp, Seq_descr_title, NULL);
    if (vnp != NULL) {
      title = (CharPtr) vnp->data.ptrvalue;
      if (StringDoesHaveText (title)) {
        if (HasUnparsedBrackets (title)) {
          reportFastaBracket = TRUE;
          for (sip = bsp->id; sip != NULL; sip = sip->next) {
            if (sip->choice != SEQID_GENERAL) continue;
            dbt = (DbtagPtr) sip->data.ptrvalue;
            if (dbt == NULL) continue;
            if (StringICmp (dbt->db, "TMSMART") == 0) {
              reportFastaBracket = FALSE;
            }
            if (StringICmp (dbt->db, "BankIt") == 0) {
              reportFastaBracket = FALSE;
            }
          }
          if (reportFastaBracket) {
            sdp = NULL;
            if (vsp->useSeqMgrIndexes) {
              sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
            } else {
              bcp = BioseqContextNew (bsp);
              sdp = BioseqContextGetSeqDescr (bcp, Seq_descr_source, NULL, NULL);
              BioseqContextFree (bcp);
            }
            if (sdp != NULL) {
              biop = (BioSourcePtr) sdp->data.ptrvalue;
              if (biop != NULL) {
                orp = biop->org;
                if (orp != NULL) {
                  if (StringDoesHaveText (orp->taxname)) {
                    if (StringChr (orp->taxname, '=') != NULL) {
                      if (StringISearch (title, orp->taxname) != NULL) {
                        reportFastaBracket = FALSE;
                      }
                    }
                  }
                }
              }
            }
          }
          if (reportFastaBracket) {
            olditemid = gcp->itemID;
            olditemtype = gcp->thistype;
            if (vnp->extended != 0) {
              ovp = (ObjValNodePtr) vnp;
              gcp->itemID = ovp->idx.itemID;
              gcp->thistype = OBJ_SEQDESC;
            }
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_FastaBracketTitle, "Title may have unparsed [...=...] construct");
            gcp->itemID = olditemid;
            gcp->thistype = olditemtype;
          }
        }

        if (StringISearch (title, "complete genome") != NULL && SequenceHasGaps (bsp)) {
          /* warning if title contains complete genome but sequence contains gap features */
          olditemid = gcp->itemID;
          olditemtype = gcp->thistype;
          gcp->itemID = bsp->idx.itemID;
          gcp->thistype = OBJ_BIOSEQ;
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_CompleteTitleProblem, "Title contains 'complete genome' but sequence has gaps");
          gcp->itemID = olditemid;
          gcp->thistype = olditemtype;
        }
      }
    } else {
      if (ISA_na (bsp->mol) && vsp->other_sets_in_sep && (vsp->is_insd_in_sep || vsp->is_refseq_in_sep) && vsp->indexerVersion) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_ComponentMissingTitle,
                  "Nucleotide component of pop/phy/mut/eco/wgs set is missing its title");
      }
    }
  }

  if (ISA_aa (bsp->mol) && vsp->useSeqMgrIndexes) {
    vnp = BioseqGetSeqDescr (bsp, Seq_descr_title, NULL);
    if (vnp != NULL) {
      if (bsp->idx.parenttype == OBJ_BIOSEQSET) {
        bssp = (BioseqSetPtr) bsp->idx.parentptr;
        while (bssp != NULL && bssp->_class != BioseqseqSet_class_nuc_prot) {
          if (bssp->idx.parenttype == OBJ_BIOSEQSET) {
            bssp = (BioseqSetPtr) bssp->idx.parentptr;
          } else {
            bssp = NULL;
          }
        }
        if (bssp != NULL && bssp->_class == BioseqseqSet_class_nuc_prot) {
          title = (CharPtr) vnp->data.ptrvalue;
          tech = 0;
          sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
          if (sdp != NULL) {
            mip = (MolInfoPtr) sdp->data.ptrvalue;
            if (mip != NULL) {
              tech = mip->tech;
            }
          }
          buf = MemNew (sizeof (Char) * (buflen + 1));
          MemSet ((Pointer) (&ii), 0, sizeof (ItemInfo));
          /* check generated protein defline with first prp->name - new convention */
          if (buf != NULL && NewCreateDefLineBuf (&ii, bsp, buf, buflen, TRUE, FALSE)) {
            if (StringICmp (buf, title) != 0) {
              /* okay if instantiated title has single trailing period */
              len2 = StringLen (buf);
              len3 = StringLen (title);
              if (len3 == len2 + 1 && title [len3 - 1] == '.' && len3 > 3 && title [len3 - 2] != '.') {
                StringCat (buf, ".");
              }
            }
            if (StringICmp (buf, title) != 0) {
              /* also check generated protein defline with all prp->names - old convention */
              if (NewCreateDefLineBuf (&ii, bsp, buf, buflen, TRUE, TRUE)) {
                bufplus = buf;
                if (StringNCmp (bufplus, "PREDICTED: ", 11) == 0) {
                  bufplus += 11;
                } else if (StringNCmp (bufplus, "UNVERIFIED: ", 12) == 0) {
                  bufplus += 12;
                }
                if (StringNCmp (title, "PREDICTED: ", 11) == 0) {
                  title += 11;
                } else if (StringNCmp (title, "UNVERIFIED: ", 12) == 0) {
                  title += 12;
                }
                if (StringICmp (bufplus, title) != 0) {
                  olditemid = gcp->itemID;
                  olditemtype = gcp->thistype;
                  if (vnp->extended != 0) {
                    ovp = (ObjValNodePtr) vnp;
                    gcp->itemID = ovp->idx.itemID;
                    gcp->thistype = OBJ_SEQDESC;
                  }
                  ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_InconsistentProteinTitle,
                            "Instantiated protein title does not match automatically generated title");
                  gcp->itemID = olditemid;
                  gcp->thistype = olditemtype;
                }
              }
            }
          }
          MemFree (buf);
        }
      }
    }
  }

  if (ISA_aa (bsp->mol) && vsp->useSeqMgrIndexes) {
    in_nps = FALSE;
    if (bsp->idx.parenttype == OBJ_BIOSEQSET) {
      bssp = (BioseqSetPtr) bsp->idx.parentptr;
      while (bssp != NULL && bssp->_class != BioseqseqSet_class_nuc_prot) {
        if (bssp->idx.parenttype == OBJ_BIOSEQSET) {
          bssp = (BioseqSetPtr) bssp->idx.parentptr;
        } else {
          bssp = NULL;
        }
      }
      if (bssp != NULL && bssp->_class == BioseqseqSet_class_nuc_prot) {
        in_nps = TRUE;
      }
    }
    if (! in_nps) {
      if (isGB || isEMBL || isDDBJ || isRefSeq) {
        if (! isGIBBMT && ! isGIBBSQ) {
          olditemid = gcp->itemID;
          olditemtype = gcp->thistype;
          gcp->itemID = bsp->idx.itemID;
          gcp->thistype = OBJ_BIOSEQ;
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_PKG_OrphanedProtein, "Orphaned stand-alone protein");
          gcp->itemID = olditemid;
          gcp->thistype = olditemtype;
        }
      }
    }
  }

  vnp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &context);
  if (vnp != NULL) {
    mip = (MolInfoPtr) vnp->data.ptrvalue;
    if (mip != NULL) {
      if (mip->completeness != 1 && isGB) {
        buf = MemNew (sizeof (Char) * (4097));
        if (buf != NULL && NewCreateDefLineBuf (NULL, bsp, buf, 4096, FALSE, FALSE)) {
          if (StringStr (buf, "complete genome") != NULL) {
            olditemid = gcp->itemID;
            olditemtype = gcp->thistype;
            if (vnp->extended != 0) {
              ovp = (ObjValNodePtr) vnp;
              gcp->itemID = ovp->idx.itemID;
              gcp->thistype = OBJ_SEQDESC;
            }
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_CompleteTitleProblem, "Complete genome in title without complete flag set");
            gcp->itemID = olditemid;
            gcp->thistype = olditemtype;
          }
        }
        MemFree (buf);
      }
      if (mip->completeness != 1 && bsp->topology == 2) {
        olditemid = gcp->itemID;
        olditemtype = gcp->thistype;
        if (vnp->extended != 0) {
          ovp = (ObjValNodePtr) vnp;
          gcp->itemID = ovp->idx.itemID;
          gcp->thistype = OBJ_SEQDESC;
        }
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_CompleteCircleProblem, "Circular topology without complete flag set");
        gcp->itemID = olditemid;
        gcp->thistype = olditemtype;
      }
    }
  }

  if (ISA_na (bsp->mol) && (bsp->repr == Seq_repr_raw || (bsp->repr == Seq_repr_delta && DeltaLitOnly (bsp))) && bsp->length > 10 && bsp->topology != 2) {
    /* check for N bases at start or stop of sequence */
    sfp = (SeqFeatPtr) MemNew (sizeof (SeqFeat));
    if (sfp == NULL) return;
    sfp->data.choice = SEQFEAT_COMMENT;

    sfp->location = AddIntervalToLocation (NULL, bsp->id, 0, 9, FALSE, FALSE);
    str = GetSequencePlusGapByFeature (sfp);
    if (str != NULL) {
      if (str [0] == 'n' || str [0] == 'N') {
        if (isNC || isPatent) {
          sev = SEV_WARNING;
        } else if (bsp->topology == TOPOLOGY_CIRCULAR) {
          sev = SEV_WARNING;
        } else if (only_local) {
          sev = SEV_WARNING;
        } else if (StringICmp (str, "NNNNNNNNNN") == 0) {
          sev = SEV_ERROR;
        } else {
          sev = SEV_WARNING;
        }
        ValidErr (vsp, sev, ERR_SEQ_INST_TerminalNs, "N at beginning of sequence");
      }
      if (str [0] == '-' || str [0] == '-') {
        if (isNC || isPatent) {
          sev = SEV_WARNING;
        } else if (bsp->topology == TOPOLOGY_CIRCULAR) {
          sev = SEV_WARNING;
        } else if (StringICmp (str, "----------") == 0) {
          sev = SEV_ERROR;
        } else {
          sev = SEV_WARNING;
        }
        ValidErr (vsp, sev, ERR_SEQ_INST_TerminalGap, "Gap at beginning of sequence");
      }
    }
    MemFree (str);
    sfp->location = SeqLocFree (sfp->location);

    sfp->location = AddIntervalToLocation (NULL, bsp->id, bsp->length - 10, bsp->length - 1, FALSE, FALSE);
    str = GetSequencePlusGapByFeature (sfp);
    len = StringLen (str);
    if (str != NULL && len > 0) {
      if (str [len - 1] == 'n' || str [len - 1] == 'N') {
        if (isNC || isPatent) {
          sev = SEV_WARNING;
        } else if (bsp->topology == TOPOLOGY_CIRCULAR) {
          sev = SEV_WARNING;
        } else if (only_local) {
          sev = SEV_WARNING;
        } else if (StringICmp (str, "NNNNNNNNNN") == 0) {
          sev = SEV_ERROR;
        } else {
          sev = SEV_WARNING;
        }
        ValidErr (vsp, sev, ERR_SEQ_INST_TerminalNs, "N at end of sequence");
      }
      if (str [len - 1] == '-' || str [len - 1] == '-') {
        if (isNC || isPatent) {
          sev = SEV_WARNING;
        } else if (bsp->topology == TOPOLOGY_CIRCULAR) {
          sev = SEV_WARNING;
        } else if (StringICmp (str, "----------") == 0) {
          sev = SEV_ERROR;
        } else {
          sev = SEV_WARNING;
        }
        ValidErr (vsp, sev, ERR_SEQ_INST_TerminalGap, "Gap at end of sequence");
      }
    }
    MemFree (str);
    sfp->location = SeqLocFree (sfp->location);

    MemFree (sfp);
  }
}

/*****************************************************************************
*
*   ValidatePubdesc(gcp)
*      Check pubdesc for missing information
*
*****************************************************************************/
static Boolean HasNoText (CharPtr str)
{
  Char            ch;

  if (str != NULL) {
    ch = *str;
    while (ch != '\0') {
      if (ch > ' ') {
        return FALSE;
      }
      str++;
      ch = *str;
    }
  }
  return TRUE;
}

static Boolean HasNoName (ValNodePtr name)
{
  AuthorPtr       ap;
  NameStdPtr      nsp;
  PersonIdPtr     pid;

  if (name != NULL) {
    ap = name->data.ptrvalue;
    if (ap != NULL) {
      pid = ap->name;
      if (pid != NULL) {
        if (pid->choice == 2) {
          nsp = pid->data;
          if (nsp != NULL) {
            if (!HasNoText (nsp->names[0])) {
              return FALSE;
            }
          }
        } else if (pid->choice == 5) {
          /* consortium */
          if (!HasNoText ((CharPtr) pid->data)) {
            return FALSE;
          }
        }
      }
    }
  }
  return TRUE;
}

static void ValidateAffil (ValidStructPtr vsp, AffilPtr ap)

{
  if (ap != NULL) {
    if (ap->affil == NULL && ap->div == NULL && ap->street == NULL && ap->city == NULL &&
        ap->sub == NULL && ap->postal_code == NULL && ap->country == NULL &&
        ap->phone == NULL && ap->fax == NULL && ap->email == NULL) {
      /* no affiliation */
    } else {
      if (ap->choice == 2) {
        /*
        if (StringHasNoText (ap->city)) {
          ValidErr (vsp, SEV_WARNING, ERR_GENERIC_MissingPubInfo, "Submission citation affiliation has no city");
        }
        */
        if (StringHasNoText (ap->country)) {
          ValidErr (vsp, SEV_WARNING, ERR_GENERIC_MissingPubInfo, "Submission citation affiliation has no country");
        }
        if (StringCmp (ap->country, "USA") == 0) {
          if (StringHasNoText (ap->sub)) {
            ValidErr (vsp, SEV_WARNING, ERR_GENERIC_MissingPubInfo, "Submission citation affiliation has no state");
          }
        }
      }
    }
  }
}

#define DATE_OKAY        0
#define EMPTY_DATE       1
#define BAD_DATE_STR     2
#define BAD_DATE_YEAR    4
#define BAD_DATE_MONTH   8
#define BAD_DATE_DAY    16
#define BAD_DATE_SEASON 32
#define BAD_DATE_OTHER  64

static Boolean DateIsBad (DatePtr dp, Boolean needFullDate, Int2Ptr baddatep)

{
  Char     ch;
  CharPtr  ptr;
  Int2     rsult = DATE_OKAY;

  if (dp == NULL) {
    rsult = EMPTY_DATE;
  } else if (dp->data [0] == 0) {
    if (dp->str == NULL) {
      rsult = BAD_DATE_STR;
    } else if (StringCmp (dp->str, "?") == 0) {
      rsult = BAD_DATE_STR;
    }
  } else if (dp->data [0] == 1) {
    if (dp->data [1] == 0) {
      rsult |= BAD_DATE_YEAR;
    }
    if (dp->data [2] > 12) {
      rsult |= BAD_DATE_MONTH;
    }
    if (dp->data [3] > 31) {
      rsult |= BAD_DATE_DAY;
    }
    if (needFullDate) {
      if (dp->data [2] == 0) {
        rsult |= BAD_DATE_MONTH;
      }
      if (dp->data [3] == 0) {
        rsult |= BAD_DATE_DAY;
      }
    }
    if (StringDoesHaveText (dp->str)) {
      ptr = dp->str;
      ch = *ptr;
      while (ch != '\0') {
        if (IS_ALPHA (ch) || ch == '-') {
        } else {
          rsult |= BAD_DATE_SEASON;
        }
        ptr++;
        ch = *ptr;
      }
    }
  } else {
    rsult = BAD_DATE_OTHER;
  }

  if (baddatep != NULL) {
    *baddatep = rsult;
  }

  return (Boolean) (rsult > 0);
}

static void PrintBadDateError (ValidStructPtr vsp, Int2 baddate, int severity, int code1, int code2, CharPtr mssg)

{
  Char  buf [256];

  buf [0] = '\0';

  if ((baddate & EMPTY_DATE) != 0) {
    StringCat (buf, "EMPTY_DATE ");
  }
  if ((baddate & BAD_DATE_STR) != 0) {
    StringCat (buf, "BAD_STR ");
  }
  if ((baddate & BAD_DATE_YEAR) != 0) {
    StringCat (buf, "BAD_YEAR ");
  }
  if ((baddate & BAD_DATE_MONTH) != 0) {
    StringCat (buf, "BAD_MONTH ");
  }
  if ((baddate & BAD_DATE_DAY) != 0) {
    StringCat (buf, "BAD_DAY ");
  }
  if ((baddate & BAD_DATE_SEASON) != 0) {
    StringCat (buf, "BAD_SEASON ");
  }
  if ((baddate & BAD_DATE_OTHER) != 0) {
    StringCat (buf, "BAD_OTHER ");
  }

  TrimSpacesAroundString (buf);

  ValidErr (vsp, severity, code1, code2, "%s - %s", mssg, buf);
}

static void ValidateCitSub (ValidStructPtr vsp, CitSubPtr csp)
{
  AffilPtr        ap;
  AuthListPtr     alp;
  Int2            baddate;
  DatePtr         dp;
  ValNodePtr      name;
  Boolean         hasAffil = FALSE;
  Boolean         hasName = FALSE;
  ErrSev          sev;

  if (vsp == NULL || csp == NULL) return;

  sev = SEV_ERROR;
  if (vsp->is_refseq_in_sep) {
    sev = SEV_WARNING;
  }
  if (vsp->is_htg_in_sep) {
    sev = SEV_WARNING;
  }

  alp = csp->authors;
  if (alp != NULL) {
    if (alp->choice == 1) {
      for (name = alp->names; name != NULL; name = name->next) {
        if (!HasNoName (name)) {
          hasName = TRUE;
        }
      }
    } else if (alp->choice == 2 || alp->choice == 3) {
      for (name = alp->names; name != NULL; name = name->next) {
        if (!HasNoText ((CharPtr) name->data.ptrvalue)) {
          hasName = TRUE;
        }
      }
    }
    ap = alp->affil;
    if (ap != NULL) {
      if (ap->affil == NULL && ap->div == NULL && ap->street == NULL && ap->city == NULL &&
           ap->sub == NULL && ap->postal_code == NULL && ap->country == NULL &&
           ap->phone == NULL && ap->fax == NULL && ap->email == NULL) {
        /* no affiliation */
      } else {
        hasAffil = TRUE;
        if (ap->choice == 2) {
          /*
          if (StringHasNoText (ap->city)) {
            ValidErr (vsp, SEV_WARNING, ERR_GENERIC_MissingPubInfo, "Submission citation affiliation has no city");
          }
          */
          if (StringHasNoText (ap->country)) {
            ValidErr (vsp, SEV_WARNING, ERR_GENERIC_MissingPubInfo, "Submission citation affiliation has no country");
          }
          if (StringCmp (ap->country, "USA") == 0) {
            if (StringHasNoText (ap->sub)) {
              ValidErr (vsp, SEV_WARNING, ERR_GENERIC_MissingPubInfo, "Submission citation affiliation has no state");
            }
          }
        }
      }
    }
  }
  if (!hasName) {
    ValidErr (vsp, SEV_ERROR, ERR_GENERIC_MissingPubInfo, "Submission citation has no author names");
  }
  if (!hasAffil) {
    if (! vsp->is_patent_in_sep) {
      ValidErr (vsp, sev, ERR_GENERIC_MissingPubInfo, "Submission citation has no affiliation");
    }
  }
  dp = csp->date;
  if (dp != NULL) {
    if (DateIsBad (dp, FALSE, &baddate)) {
      PrintBadDateError (vsp, baddate, SEV_ERROR, ERR_GENERIC_BadDate, "Submission citation date has error");
    }
  }
}

static void LookForMultiplePubs (ValidStructPtr vsp, GatherContextPtr gcp, SeqDescrPtr sdp)

{
  Bioseq       bs;
  Boolean      collision, otherpub;
  Int4         muid, pmid;
  SeqDescrPtr  nextpub;
  PubdescPtr   pdp;
  ValNodePtr   vnp;


  if (sdp != NULL && sdp->choice == Seq_descr_pub && sdp->extended != 0 && vsp != NULL && gcp != NULL) {
    MemSet ((Pointer) &bs, 0, sizeof (Bioseq));
    pdp = (PubdescPtr) sdp->data.ptrvalue;
    if (pdp != NULL) {
      otherpub = FALSE;
      muid = 0;
      pmid = 0;
      for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
        if (vnp->choice == PUB_Muid) {
          muid = vnp->data.intvalue;
        } else if (vnp->choice == PUB_PMid) {
          pmid = vnp->data.intvalue;
        } else {
          otherpub = TRUE;
        }
      }
      if (otherpub) {
        if (muid > 0 || pmid > 0) {
          collision = FALSE;
          nextpub = GetNextDescriptorUnindexed (&bs, Seq_descr_pub, sdp);
          while (nextpub != NULL) {
            pdp = (PubdescPtr) nextpub->data.ptrvalue;
            if (pdp != NULL) {
              for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
                if (vnp->choice == PUB_Muid) {
                  if (muid > 0 && muid == vnp->data.intvalue) {
                    collision = TRUE;
                  }
                } else if (vnp->choice == PUB_PMid) {
                  if (pmid > 0 && pmid == vnp->data.intvalue) {
                    collision = TRUE;
                  }
                }
              }
            }
            nextpub = GetNextDescriptorUnindexed (&bs, Seq_descr_pub, nextpub);
          }
          if (collision) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_CollidingPublications, "Multiple publications with same identifier");
          }
        }
      }
    }
  }
}

static void LookForMultipleUnpubPubs (ValidStructPtr vsp, GatherContextPtr gcp, BioseqPtr bsp)

{
  Char               buf [2048];
  CharPtr            last, str;
  SeqMgrDescContext  dcontext;
  ValNodePtr         list = NULL, next, vnp;
  ObjValNodePtr      ovp;
  PubdescPtr         pdp;
  SeqDescrPtr        sdp;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_pub, &dcontext);
  while (sdp) {
    pdp = (PubdescPtr) sdp->data.ptrvalue;
    if (pdp != NULL) {
      ovp = (ObjValNodePtr) sdp;
      if (ovp->idx.scratch != NULL) {
        ValNodeCopyStr (&list, 0, ovp->idx.scratch);
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_pub, &dcontext);
  }

  if (list == NULL) return;

  list = ValNodeSort (list, SortVnpByString);
  last = (CharPtr) list->data.ptrvalue;
  vnp = list->next;
  while (vnp != NULL) {
    next = vnp->next;
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringICmp (last, str) == 0) {
      StringNCpy_0 (buf, str, sizeof (buf));
      StringCpy (buf + 100, "...");
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_CollidingPublications,
                "Multiple equivalent publications annotated on this sequence [%s]", buf);
    } else {
      last = (CharPtr) vnp->data.ptrvalue;
    }
    vnp = next;
  }

  ValNodeFreeData (list);
}

static Boolean BadCharsInAuth (CharPtr str, CharPtr PNTR badauthor, Boolean allowcomma, Boolean allowperiod, Boolean last)
{
  Char     ch;
  CharPtr  ptr, stp = NULL;

  if (StringHasNoText (str)) return FALSE;
  if (last) {
      stp = StringISearch(str, "St.");
      if (stp == str) {
          stp += 2;  /* point to the period */
      }
  }

  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    /* success on any of these tests are allowed values */
    if (IS_ALPHA (ch)) {
    } else if (ch == '-' || ch == '\'' || ch == ' ') {
    } else if (ch == ',' && allowcomma) {
    } else if (ch == '.' && (allowperiod || stp == ptr)) {
    } else {
      /* bad character found */
      *badauthor = str;
      return TRUE;
    }
    ptr++;
    ch = *ptr;
  }

  return FALSE;
}

static Boolean BadCharsInName (ValNodePtr name, CharPtr PNTR badauthor, BoolPtr last_name_badP)

{
  AuthorPtr    ap;
  NameStdPtr   nsp;
  PersonIdPtr  pid;

  if (name == NULL) return FALSE;
  ap = name->data.ptrvalue;
  if (ap == NULL) return FALSE;
  pid = ap->name;
  if (pid == NULL) return FALSE;

  if (pid->choice == 2) {
    nsp = pid->data;
    if (nsp == NULL) return FALSE;
    if (StringICmp (nsp->names [0], "et al.") == 0) return FALSE;
    if (BadCharsInAuth (nsp->names [0], badauthor, FALSE, FALSE, TRUE)) {
      if (last_name_badP != NULL) {
        *last_name_badP = TRUE;
      }
      return TRUE; /* last    */
    }
    if (BadCharsInAuth (nsp->names [1], badauthor, FALSE, FALSE, FALSE)) return TRUE; /* first    */
    if (BadCharsInAuth (nsp->names [4], badauthor, FALSE, TRUE, FALSE)) return TRUE;  /* initials */
    if (BadCharsInAuth (nsp->names [5], badauthor, FALSE, TRUE, FALSE)) return TRUE;  /* suffix */
  }

  return FALSE;
}

static CharPtr suffixList [] = {
  "Jr.", "Sr.", "II", "III", "IV", "V", "VI", NULL
};

static void ValidateSuffix (ValidStructPtr vsp, GatherContextPtr gcp, PubdescPtr pdp, ValNodePtr name)

{
  AuthorPtr    ap;
  Int2         i;
  NameStdPtr   nsp;
  PersonIdPtr  pid;
  CharPtr      suffix;

  if (name == NULL) return;
  ap = name->data.ptrvalue;
  if (ap == NULL) return;
  pid = ap->name;
  if (pid == NULL) return;

  if (pid->choice == 2) {
    nsp = pid->data;
    if (nsp == NULL) return;
    suffix = nsp->names [5];
    if (StringHasNoText (suffix)) return;
    for (i = 0; suffixList [i] != NULL; i++) {
      if (StringICmp (suffix, suffixList [i]) == 0) return;
    }
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadAuthorSuffix, "Bad author suffix %s", suffix);
  }
}

static Boolean StringAlreadyInList (
  ValNodePtr head,
  CharPtr str
)

{
  ValNodePtr  vnp;

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    if (StringICmp ((CharPtr) vnp->data.ptrvalue, str) == 0) return TRUE;
  }

  return FALSE;
}

static void ValidatePubdesc (ValidStructPtr vsp, GatherContextPtr gcp, PubdescPtr pdp)
{
  AuthListPtr     alp;
  AuthorPtr       ap;
  CharPtr         badauthor;
  Int2            baddate;
  CitArtPtr       cap = NULL;
  CitGenPtr       cgp;
  CitJourPtr      cjp = NULL;
  ValNodePtr      conslist = NULL;
  CitSubPtr       csp;
  DatePtr         dp;
  Boolean         hasName, hasTitle, hasIsoJTA = FALSE,
                  inPress = FALSE, electronic_journal = FALSE,
                  conflicting_pmids = FALSE, redundant_pmids = FALSE,
                  conflicting_muids = FALSE, redundant_muids = FALSE,
                  unpub = FALSE;
  ImprintPtr      imp;
  Boolean         last_name_bad;
  Int4            muid = 0;
  Boolean         noVol = FALSE, noPages = FALSE;
  ValNodePtr      name;
  PersonIdPtr     pid;
  Int4            pmid = 0;
  CharPtr         ptr;
  ErrSev          sev;
  Int4            start;
  Int4            stop;
  CharPtr         str;
  Char            temp [64];
  Int4            thepmid = 0;
  ValNodePtr      title;
  Int4            uid = 0;
  long int        val;
  ValNodePtr      vnp;

  if (vsp == NULL || pdp == NULL) return;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == PUB_PMid) {
      thepmid = vnp->data.intvalue;
    }
  }
  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
    case PUB_Gen:
      cgp = (CitGenPtr) vnp->data.ptrvalue;
      hasName = FALSE;
      if (cgp != NULL) {
        if (StringDoesHaveText (cgp->cit)) {
          /* skip if just BackBone id number */
          if (StringNICmp (cgp->cit, "BackBone id_pub = ", 18) == 0 && cgp->journal == NULL && cgp->date == NULL && cgp->serial_number < 0) break;
          if (StringNICmp (cgp->cit, "submitted", 8) == 0 ||
              StringNICmp (cgp->cit, "unpublished", 11) == 0 ||
              StringNICmp (cgp->cit, "Online Publication", 18) == 0 ||
              StringNICmp (cgp->cit, "Published Only in DataBase", 26) == 0) {
            unpub = TRUE;
          } else if (StringNICmp (cgp->cit, "(er) ", 5) == 0) {
            unpub = TRUE;
          } else {
            ValidErr (vsp, SEV_ERROR, ERR_GENERIC_MissingPubInfo, "Unpublished citation text invalid");
          }
          if (StringStr (cgp->cit, "Title=") != NULL) {
            ValidErr (vsp, SEV_ERROR, ERR_GENERIC_StructuredCitGenCit, "Unpublished citation has embedded Title");
          }
          if (StringStr (cgp->cit, "Journal=") != NULL) {
            ValidErr (vsp, SEV_ERROR, ERR_GENERIC_StructuredCitGenCit, "Unpublished citation has embedded Journal");
          }
        }
        /* skip if just serial number */
        if (cgp->cit == NULL && cgp->journal == NULL && cgp->date == NULL && cgp->serial_number > -1) break;
        dp = cgp->date;
        if (dp == NULL) {
          if (! unpub) {
            ValidErr (vsp, SEV_WARNING, ERR_GENERIC_MissingPubInfo, "Publication date missing");
          }
        } else if (dp->str != NULL) {
          if (StringCmp (dp->str, "?") == 0) {
            ValidErr (vsp, SEV_WARNING, ERR_GENERIC_MissingPubInfo, "Publication date marked as '?'");
          }
        } else if (dp->data [1] == 0) {
          ValidErr (vsp, SEV_WARNING, ERR_GENERIC_MissingPubInfo, "Publication date not set");
        } else if (DateIsBad (dp, FALSE, &baddate)) {
          PrintBadDateError (vsp, baddate, SEV_ERROR, ERR_GENERIC_BadDate, "Publication date has error");
        }
        alp = cgp->authors;
        if (alp != NULL) {
          if (alp->choice == 1) {
            for (name = alp->names; name != NULL; name = name->next) {
              if (!HasNoName (name)) {
                hasName = TRUE;
              }
            }
          } else if (alp->choice == 2 || alp->choice == 3) {
            for (name = alp->names; name != NULL; name = name->next) {
              if (!HasNoText ((CharPtr) name->data.ptrvalue)) {
                hasName = TRUE;
              }
            }
          }
        }
        if (!hasName) {
          sev = SEV_ERROR;
          if (vsp->is_refseq_in_sep) {
            sev = SEV_WARNING;
          }
          ValidErr (vsp, sev, ERR_GENERIC_MissingPubInfo, "Publication has no author names");
        }
      }
      break;
    case PUB_PMid:
      if (pmid == 0) {
        pmid = vnp->data.intvalue;
      } else if (pmid != vnp->data.intvalue) {
        conflicting_pmids = TRUE;
      } else {
        redundant_pmids = TRUE;
      }
      if (uid == 0) {
        uid = vnp->data.intvalue;
      }
      break;
    case PUB_Muid:
      if (muid == 0) {
        muid = vnp->data.intvalue;
      } else if (muid != vnp->data.intvalue) {
        conflicting_muids = TRUE;
      } else {
        redundant_muids = TRUE;
      }
      if (uid == 0) {
        uid = vnp->data.intvalue;
      }
      break;
    case PUB_Sub :
      csp = (CitSubPtr) vnp->data.ptrvalue;
      if (csp != NULL) {
        ValidateCitSub (vsp, csp);
      }
      break;
    case PUB_Medline:
      ValidErr (vsp, SEV_ERROR, ERR_GENERIC_MedlineEntryPub, "Publication is medline entry");
      break;
    case PUB_Article:
      cap = (CitArtPtr) vnp->data.ptrvalue;
      hasName = FALSE;
      hasTitle = FALSE;
      if (cap != NULL) {
        for (title = cap->title; title != NULL; title = title->next) {
          if (!HasNoText ((CharPtr) title->data.ptrvalue)) {
            hasTitle = TRUE;
          }
        }
        if (!hasTitle) {
          ValidErr (vsp, SEV_ERROR, ERR_GENERIC_MissingPubInfo, "Publication has no title");
        }
        alp = cap->authors;
        if (alp != NULL) {
          if (alp->choice == 1) {
            for (name = alp->names; name != NULL; name = name->next) {
              if (!HasNoName (name)) {
                hasName = TRUE;
              }
            }
          } else if (alp->choice == 2 || alp->choice == 3) {
            for (name = alp->names; name != NULL; name = name->next) {
              if (!HasNoText ((CharPtr) name->data.ptrvalue)) {
                hasName = TRUE;
              }
            }
          }
        }
        if (!hasName) {
          ValidErr (vsp, SEV_ERROR, ERR_GENERIC_MissingPubInfo, "Publication has no author names");
        }
      }

      if (cap != NULL) {
        switch (cap->from) {
        case 1:
          cjp = (CitJourPtr) cap->fromptr;
          if (cjp != NULL) {
            hasTitle = FALSE;
            for (title = cjp->title; title != NULL; title = title->next) {
              if (title->choice == Cit_title_iso_jta) {
                hasIsoJTA = TRUE;
              }
              if (!HasNoText ((CharPtr) title->data.ptrvalue)) {
                hasTitle = TRUE;
                if (title->choice == Cit_title_name) {
                  if (StringNCmp ((CharPtr) title->data.ptrvalue, "(er)", 4) == 0) {
                    electronic_journal = TRUE;
                  }
                }
              }
            }
            if (!hasTitle) {
              ValidErr (vsp, SEV_ERROR, ERR_GENERIC_MissingPubInfo, "Journal title missing");
            }
            imp = cjp->imp;
            if (imp != NULL) {
              if (imp->pubstatus == PUBSTATUS_epublish || imp->pubstatus == PUBSTATUS_aheadofprint) {
                electronic_journal = TRUE;
              }
              if (imp->prepub == 2) {
                inPress = TRUE;
                if (StringDoesHaveText (imp->pages)) {
                  ValidErr (vsp, SEV_WARNING, ERR_GENERIC_PublicationInconsistency, "In-press is not expected to have page numbers");
                }
              }
              if (imp->prepub == 0 && imp->pubstatus != PUBSTATUS_aheadofprint) {
                noVol = StringHasNoText (imp->volume);
                noPages = StringHasNoText (imp->pages);
                sev = SEV_ERROR;
                if (vsp->is_refseq_in_sep) {
                  sev = SEV_WARNING;
                }
                if (imp->pubstatus == PUBSTATUS_epublish) {
                  sev = SEV_WARNING;
                }
                if (noVol && noPages) {
                  ValidErr (vsp, sev, ERR_GENERIC_MissingPubInfo, "Journal volume and pages missing");
                } else if (noVol) {
                  if (! electronic_journal) {
                    ValidErr (vsp, sev, ERR_GENERIC_MissingPubInfo, "Journal volume missing");
                  }
                } else if (noPages) {
                  ValidErr (vsp, sev, ERR_GENERIC_MissingPubInfo, "Journal pages missing");
                }
                if ((! noPages) && (! electronic_journal)) {
                  sev = SEV_WARNING;
                  StringNCpy_0 (temp, imp->pages, sizeof (temp));
                  ptr = StringChr (temp, '-');
                  if (ptr != NULL) {
                    *ptr = '\0';
                    ptr++;
                    if (sscanf (temp, "%ld", &val) == 1) {
                      start = (Int4) val;
                      if (sscanf (ptr, "%ld", &val) == 1) {
                        stop = (Int4) val;
                        if (start == 0 || stop == 0) {
                          ValidErr (vsp, sev, ERR_GENERIC_BadPageNumbering, "Page numbering has zero value");
                        } else if (start < 0 || stop < 0) {
                          ValidErr (vsp, sev, ERR_GENERIC_BadPageNumbering, "Page numbering has negative value");
                        } else if (start > stop) {
                          ValidErr (vsp, sev, ERR_GENERIC_BadPageNumbering, "Page numbering out of order");
                        } else if (stop > start + 50) {
                          ValidErr (vsp, sev, ERR_GENERIC_BadPageNumbering, "Page numbering greater than 50");
                        }
                      } else {
                        ValidErr (vsp, sev, ERR_GENERIC_BadPageNumbering, "Page numbering stop looks strange");
                      }
                    } else if (! IS_ALPHA (temp [0])) {
                      ValidErr (vsp, sev, ERR_GENERIC_BadPageNumbering, "Page numbering start looks strange");
                    }
                  }
                }
                dp = imp->date;
                if (dp == NULL) {
                  ValidErr (vsp, SEV_WARNING, ERR_GENERIC_MissingPubInfo, "Publication date missing");
                } else if (dp->str != NULL) {
                  if (StringCmp (dp->str, "?") == 0) {
                    ValidErr (vsp, SEV_WARNING, ERR_GENERIC_MissingPubInfo, "Publication date marked as '?'");
                  }
                } else if (dp->data [1] == 0) {
                  ValidErr (vsp, SEV_WARNING, ERR_GENERIC_MissingPubInfo, "Publication date not set");
                } else if (DateIsBad (dp, FALSE, &baddate)) {
                  PrintBadDateError (vsp, baddate, SEV_ERROR, ERR_GENERIC_BadDate, "Publication date has error");
                }
              }
              if (imp->pubstatus == PUBSTATUS_aheadofprint && imp->prepub != 2) {
                if (noVol || noPages) {
                } else {
                  ValidErr (vsp, SEV_WARNING, ERR_GENERIC_PublicationInconsistency, "Ahead-of-print without in-press");
                }
              }
              if (imp->pubstatus == PUBSTATUS_epublish && imp->prepub == 2) {
                ValidErr (vsp, SEV_WARNING, ERR_GENERIC_PublicationInconsistency, "Electronic-only publication should not also be in-press");
              }
              if (imp->pubstatus == PUBSTATUS_epublish || imp->pubstatus == PUBSTATUS_ppublish || imp->pubstatus == PUBSTATUS_aheadofprint) {
                if (StringDoesHaveText (pdp->comment)) {
                  if (StringStr (pdp->comment, "Publication Status") != NULL ||
                      StringStr (pdp->comment, "Publication-Status") != NULL ||
                      StringStr (pdp->comment, "Publication_Status") != NULL) {
                    ValidErr (vsp, SEV_ERROR, ERR_GENERIC_UnexpectedPubStatusComment, "Publication status is in comment for pmid %ld", (long) thepmid);
                  }
                }
              }
            }
          }
          break;
        default:
          break;
        }
      }
      break;
    case PUB_Equiv:
      ValidErr (vsp, SEV_WARNING, ERR_GENERIC_UnnecessaryPubEquiv, "Publication has unexpected internal Pub-equiv");
      break;
    default:
      break;
    }
  }

  if (conflicting_pmids) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_CollidingPublications, "Multiple conflicting pmids in a single publication");
  } else if (redundant_pmids) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_CollidingPublications, "Multiple redundant pmids in a single publication");
  }
  if (conflicting_muids) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_CollidingPublications, "Multiple conflicting muids in a single publication");
  } else if (redundant_muids) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_CollidingPublications, "Multiple redundant muids in a single publication");
  }

  if (cap != NULL && cjp != NULL && (uid > 0 || inPress || vsp->alwaysRequireIsoJTA)) {
    if (! hasIsoJTA) {
      if (! electronic_journal) {
        ValidErr (vsp, SEV_WARNING, ERR_GENERIC_MissingPubInfo, "ISO journal title abbreviation missing");
      }
    }
  }

  alp = GetAuthListPtr (pdp, NULL);
  if (alp != NULL) {
    sev = SEV_ERROR;
    if (vsp->is_refseq_in_sep) {
      sev = SEV_WARNING;
    }
    if (alp->choice == 1) {
      for (name = alp->names; name != NULL; name = name->next) {
        badauthor = NULL;
        last_name_bad = FALSE;
        if (BadCharsInName (name, &badauthor, &last_name_bad)) {
          if (StringHasNoText (badauthor)) {
            badauthor = "?";
          }
          if (last_name_bad) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadCharInAuthorLastName, "Bad characters in author %s", badauthor);
          } else {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadCharInAuthorName, "Bad characters in author %s", badauthor);
          }
        }
        ValidateSuffix (vsp, gcp, pdp, name);
        ap = (AuthorPtr) name->data.ptrvalue;
        if (ap == NULL) continue;
        pid = ap->name;
        if (pid == NULL) continue;
        if (pid->choice == 5) {
          str = (CharPtr) pid->data;
          if (StringHasNoText (str)) {
            ValidErr (vsp, SEV_WARNING, ERR_GENERIC_PublicationInconsistency, "Empty consortium");
            continue;
          }
          if (StringAlreadyInList (conslist, str)) {
            ValidErr (vsp, SEV_WARNING, ERR_GENERIC_PublicationInconsistency, "Duplicate consortium '%s'", str);
            continue;
          }
          ValNodeAddPointer (&conslist, 0, (Pointer) str);
        }
      }
    } else if (alp->choice == 2 || alp->choice == 3) {
      for (name = alp->names; name != NULL; name = name->next) {
        badauthor = NULL;
        if (BadCharsInAuth ((CharPtr) name->data.ptrvalue, &badauthor, TRUE, TRUE, FALSE)) {
          if (StringHasNoText (badauthor)) {
            badauthor = "?";
          }
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadCharInAuthorName, "Bad characters in author %s", badauthor);
        }
      }
    }
  }
  ValNodeFree (conslist);
}

static void ValidateSfpCit (ValidStructPtr vsp, GatherContextPtr gcp, SeqFeatPtr sfp)
{
  ValNodePtr      ppr;
  ValNodePtr      psp;

  if (vsp == NULL || sfp == NULL || sfp->cit == NULL)
    return;
  psp = sfp->cit;
  if (psp == NULL)
    return;
  for (ppr = (ValNodePtr) psp->data.ptrvalue; ppr != NULL; ppr = ppr->next) {
    if (ppr->choice == PUB_Equiv) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UnnecessaryCitPubEquiv, "Citation on feature has unexpected internal Pub-equiv");
      return;
    }
  }
}

typedef struct bioseqvalid
{
  ValidStructPtr  vsp;
  Boolean         is_aa;         /* bioseq is protein? */
  Boolean         is_mrna;       /* molinfo is mrna? */
  Boolean         is_prerna;     /* molinfo is precursor rna? */
  Boolean         is_artificial; /* biosource origin is artificial */
  Boolean         is_synthetic;  /* biosource origin synthetic */
  Boolean         is_syn_constr; /* is organism name synthetic construct, plasmid, vector, or SYN division? */
  Boolean         got_a_pub;
  int             last_na_mol, last_na_mod, last_organelle, last_partialness, last_left_right,
                  last_biomol, last_tech, last_completeness,
                  num_full_length_src_feat, num_full_length_prot_ref,
                  num_justprot, num_preprot, num_matpep, num_sigpep, num_transpep;
  ValNodePtr      last_gb, last_embl, last_prf, last_pir, last_sp, last_pdb,
                  last_create, last_update, last_biosrc, last_orgref;
  OrgRefPtr       last_org;
  GatherContextPtr gcp;
  BioseqPtr        bsp;
}
BioseqValidStr , PNTR BioseqValidStrPtr;

static void CheckForNucProt (BioseqSetPtr bssp, Pointer userdata)
{
  BoolPtr         hasPartsP;

  if (bssp->_class == BioseqseqSet_class_nuc_prot) {
    hasPartsP = (BoolPtr) userdata;
    *hasPartsP = TRUE;
  }
}

static void CheckForParts (BioseqSetPtr bssp, Pointer userdata)
{
  BoolPtr         hasPartsP;

  if (bssp->_class == BioseqseqSet_class_parts) {
    hasPartsP = (BoolPtr) userdata;
    *hasPartsP = TRUE;
  }
}

static Boolean DeltaOrFarSeg (SeqEntryPtr sep, SeqLocPtr location)
{
  BioseqPtr       bsp;
  Boolean         hasParts = FALSE;

  bsp = BioseqFindFromSeqLoc (location);
  if (bsp != NULL) {
    if (bsp->repr == Seq_repr_delta) {
      VisitSetsInSep (sep, (Pointer) &hasParts, CheckForNucProt);
      if (!hasParts)
        return TRUE;
    }
    if (bsp->repr == Seq_repr_seg) {
      VisitSetsInSep (sep, (Pointer) &hasParts, CheckForParts);
      if (!hasParts)
        return TRUE;
    }
  }
  return FALSE;
}


static Boolean IsOrganelleBioseq (BioseqPtr bsp)
{
  SeqDescrPtr sdp;
  SeqMgrDescContext dcontext;
  BioSourcePtr biop;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp == NULL || (biop = sdp->data.ptrvalue) == NULL || !IsLocationOrganelle (biop->genome)) {
    return FALSE;
  } else {
    return TRUE;
  }
}


static void
ValidateIntronEndsAtSpliceSiteOrGap
(ValidStructPtr vsp,
 SeqLocPtr slp)
{
  BioseqPtr          bsp;
  SeqIdPtr           sip;
  Uint1              strand;
  Int4               strt, stop, pos;
  Boolean            partial5, partial3;
  Char               buf[3];
  Char               id_buf[150];
  SeqFeatPtr         rna;
  SeqMgrFeatContext  rcontext;

  if (vsp == NULL || slp == NULL) return;
  CheckSeqLocForPartial (slp, &partial5, &partial3);
  if (partial5 && partial3) return;

  /* suppress if contained by rRNA - different consensus splice site */

  rna = SeqMgrGetOverlappingFeature (slp, 0, vsp->rrna_array, vsp->numrrna,
                                     NULL, CONTAINED_WITHIN, &rcontext);
  if (rna != NULL) return;

  /* suppress if contained by tRNA - different consensus splice site */

  rna = SeqMgrGetOverlappingFeature (slp, 0, vsp->trna_array, vsp->numtrna,
                                     NULL, CONTAINED_WITHIN, &rcontext);
  if (rna != NULL) return;


  sip = SeqLocId (slp);
  if (sip == NULL)
    return;

  bsp = NULL;
  if (sip != NULL && (sip->choice != SEQID_GI || sip->data.intvalue > 0)) {
    bsp = BioseqLockById (sip);
  }
  if (bsp == NULL)
    return;

  if (IsOrganelleBioseq(bsp)) {
    BioseqUnlock (bsp);
    return;
  }

  BioseqLabel (bsp, id_buf, sizeof (id_buf) - 1, OM_LABEL_CONTENT);

  strt = SeqLocStart (slp);
  stop = SeqLocStop (slp);

  strand = SeqLocStrand (slp);

  if (!partial5) {
    if (strand == Seq_strand_minus) {
      SeqPortStreamInt (bsp, stop - 1, stop, Seq_strand_minus, EXPAND_GAPS_TO_DASHES, (Pointer) buf, NULL);
      pos = stop;
    } else {
      SeqPortStreamInt (bsp, strt, strt + 1, Seq_strand_plus, EXPAND_GAPS_TO_DASHES, (Pointer) buf, NULL);
      pos = strt;
    }
    if ((buf[0] == '-' && buf[1] == '-')
        || (buf[0] == 'G' && buf[1] == 'T')
        || (buf[0] == 'G' && buf[1] == 'C')) {
      /* location is ok */
    } else if (pos == 0 || pos == bsp->length - 1) {
      ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_NotSpliceConsensusDonor,
                "Splice donor consensus (GT) not found at start of terminal intron, position %ld of %s", (long) (pos + 1), id_buf);
    } else {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_NotSpliceConsensusDonor,
                "Splice donor consensus (GT) not found at start of intron, position %ld of %s", (long) (pos + 1), id_buf);
    }
  }
  if (!partial3) {
    if (strand == Seq_strand_minus) {
      SeqPortStreamInt (bsp, strt, strt + 1, Seq_strand_minus, EXPAND_GAPS_TO_DASHES, (Pointer) buf, NULL);
      pos = strt;
    } else {
      SeqPortStreamInt (bsp, stop - 1, stop, Seq_strand_plus, EXPAND_GAPS_TO_DASHES, (Pointer) buf, NULL);
      pos = stop;
    }
    if ((buf[0] == '-' && buf[1] == '-')
        || (buf[0] == 'A' && buf[1] == 'G')) {
      /* location is ok */
    } else if (pos == 0 || pos == bsp->length - 1) {
      ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_NotSpliceConsensusAcceptor,
                "Splice acceptor consensus (AG) not found at end of terminal intron, position %ld of %s, but at end of sequence", (long) (pos + 1), id_buf);
    } else {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_NotSpliceConsensusAcceptor,
                "Splice acceptor consensus (AG) not found at end of intron, position %ld of %s", (long) (pos + 1), id_buf);
    }
  }
  BioseqUnlock (bsp);
}

static Boolean IsLocInSmallGenomeSet (
  SeqLocPtr loc
)

{
  BioseqPtr  bsp;
  SeqIdPtr   sip;
  SeqLocPtr  slp;

  if (loc == NULL) return FALSE;

  slp = SeqLocFindNext (loc, NULL);
  while (slp != NULL) {
    sip = SeqLocId (slp);
    if (sip == NULL) return FALSE;
    bsp = BioseqFind (sip);
    if (bsp == NULL) return FALSE;
    slp = SeqLocFindNext (loc, slp);
  }

  return TRUE;
}

static Boolean AllPartsInSmallGenomeSet (
  SeqLocPtr loc,
  ValidStructPtr vsp,
  BioseqPtr bsp
)

{
  BioseqSetPtr  bssp;
  SeqEntryPtr   oldscope;
  Boolean       rsult = FALSE;
  SeqEntryPtr   sep;

  if (loc == NULL || vsp == NULL || bsp == NULL) return FALSE;

  sep = vsp->sep;
  if (sep == NULL) return FALSE;
  if (! IS_Bioseq_set (sep)) return FALSE;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp == NULL) return FALSE;

  /* if genbank set wraps everything, go down one set level */
  if (bssp->_class == BioseqseqSet_class_genbank) {
    sep = bssp->seq_set;
    if (sep == NULL) return FALSE;
    if (! IS_Bioseq_set (sep)) return FALSE;
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
  }

  /* check for small genome set */
  if (bssp == NULL || bssp->_class != BioseqseqSet_class_small_genome_set) return FALSE;

  /* scope within small genome set for subsequent BioseqFind calls */
  oldscope = SeqEntrySetScope (sep);

  rsult = IsLocInSmallGenomeSet (loc);

  SeqEntrySetScope (oldscope);

  return rsult;
}


/*****************************************************************************
*
*   ValidateSeqFeatContext(gcp)
*      Gather callback helper function for validating context on a Bioseq
*
*****************************************************************************/
static Boolean ValidateSeqFeatCommon (SeqFeatPtr sfp, BioseqValidStrPtr bvsp, ValidStructPtr vsp,
                                      Int4 left, Int4 right, Int2 numivals, Uint4 featitemid, Boolean farloc, BioseqPtr bsp)
{
  BioseqSetPtr    bssp;
  Boolean         do_error;
  GatherContextPtr gcp = NULL;
  ImpFeatPtr      ifp;
  Uint2           olditemtype = 0;
  Uint4           olditemid = 0;
  ProtRefPtr      prp;
  RnaRefPtr       rrp;
  CharPtr         str;
  SeqLocPtr       slp;
  SeqIdPtr        sip;
  TextSeqIdPtr    tsip;
  Boolean         on_seg = FALSE;
  Boolean         is_emb = FALSE;
  Boolean         is_nc = FALSE;
  Boolean         is_refseq = FALSE;
  ErrSev          sev;
  Boolean         no_nonconsensus_except;


  vsp->descr = NULL;
  vsp->sfp = sfp;

  if (featitemid > 0) {
    gcp = vsp->gcp;
    if (gcp != NULL) {
      olditemid = gcp->itemID;
      olditemtype = gcp->thistype;
      gcp->itemID = featitemid;
      gcp->thistype = OBJ_SEQFEAT;
    }
  }

  if (bsp != NULL) {
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_OTHER) {
        is_refseq = TRUE;
        tsip = (TextSeqIdPtr) sip->data.ptrvalue;
        if (tsip != NULL && tsip->accession != NULL) {
          if (StringNICmp (tsip->accession, "NT_", 3) == 0) {
            is_nc = TRUE;
          }
        }
      } else if (sip->choice == SEQID_EMBL) {
        is_emb = TRUE;
      }
    }
  }

  if (bvsp->is_aa) {
    if (sfp->data.choice == SEQFEAT_PROT) {
      if ((left == 0) && (right == ((vsp->bsp->length) - 1))) {
        bvsp->num_full_length_prot_ref++;
        prp = (ProtRefPtr) sfp->data.value.ptrvalue;
        if (prp != NULL) {
          switch (prp->processed) {
            case 0:
              bvsp->num_justprot++;
              break;
            case 1:
              bvsp->num_preprot++;
              break;
            case 2:
              bvsp->num_matpep++;
              break;
            case 3:
              bvsp->num_sigpep++;
              break;
            case 4:
              bvsp->num_transpep++;
              break;
            default:
              break;
          }
        }
      }
    }

    switch (sfp->data.choice) {
    case SEQFEAT_CDREGION:
    case SEQFEAT_RNA:
    case SEQFEAT_RSITE:
    case SEQFEAT_TXINIT:
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidForType, "Invalid feature for a protein Bioseq.");
      break;
    case SEQFEAT_GENE:
        if (bsp != NULL) {
          do_error = FALSE;
          if (bsp->idx.parenttype == OBJ_BIOSEQSET) {
            bssp = (BioseqSetPtr) bsp->idx.parentptr;
            while (bssp != NULL) {
              switch (bssp->_class) {
              case BioseqseqSet_class_nuc_prot :
              case BioseqseqSet_class_mut_set :
              case BioseqseqSet_class_pop_set :
              case BioseqseqSet_class_phy_set :
              case BioseqseqSet_class_eco_set :
              case BioseqseqSet_class_gen_prod_set :
                do_error = TRUE;
                break;
              default :
                break;
              }
              if (bssp->idx.parenttype == OBJ_BIOSEQSET) {
                bssp = (BioseqSetPtr) bssp->idx.parentptr;
              } else {
                bssp = NULL;
              }
            }
          }
          if (do_error) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidForType, "Invalid feature for a protein Bioseq.");
          }
        }
      break;
    default:
      break;
    }

  } else {
    switch (sfp->data.choice) {
    case SEQFEAT_PROT:
    case SEQFEAT_PSEC_STR:
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidForType, "Invalid feature for a nucleotide Bioseq.");
      break;
    default:
      break;
    }

  }

  if (bvsp->is_mrna) {
    switch (sfp->data.choice) {
    case SEQFEAT_CDREGION:
      if (numivals > 1) {
        if ((! sfp->excpt) ||
            (StringISearch (sfp->except_text, "ribosomal slippage") == NULL)) {
          sev = SEV_ERROR;
          if (is_refseq) {
            sev = SEV_WARNING;
          }
          ValidErr (vsp, sev, ERR_SEQ_FEAT_InvalidForType, "Multi-interval CDS feature is invalid on an mRNA (cDNA) Bioseq.");
        }
      }
      break;
    case SEQFEAT_RNA:
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp != NULL && rrp->type == 2) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidForType, "mRNA feature is invalid on an mRNA (cDNA) Bioseq.");
      }
      break;
    case SEQFEAT_IMP:
      ifp = (ImpFeatPtr) sfp->data.value.ptrvalue;
      if (ifp != NULL && ifp->key != NULL && (!HasNoText (ifp->key))) {
        if (StringCmp (ifp->key, "intron") == 0 || StringCmp (ifp->key, "CAAT_signal") == 0) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidForType, "Invalid feature for an mRNA Bioseq.");
        }
      }
      break;
    default:
      break;
    }
  } else if (bvsp->is_prerna) {
    switch (sfp->data.choice) {
    case SEQFEAT_IMP:
      ifp = (ImpFeatPtr) sfp->data.value.ptrvalue;
      if (ifp != NULL && ifp->key != NULL && (!HasNoText (ifp->key))) {
        if (StringCmp (ifp->key, "CAAT_signal") == 0) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidForType, "Invalid feature for an pre-RNA Bioseq.");
        }
      }
      break;
    default:
      break;
    }
  }

  if (farloc && (! is_nc) && (! is_emb) && (! AllPartsInSmallGenomeSet (sfp->location, vsp, bsp))) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FarLocation, "Feature has 'far' location - accession not packaged in record");
  }

  if ((sfp->data.choice == SEQFEAT_PUB) || (sfp->cit != NULL))
    bvsp->got_a_pub = TRUE;

  str = (CharPtr) sfp->comment;
  if (SerialNumberInString (str)) {
    ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_SerialInComment,
              "Feature comment may refer to reference by serial number - attach reference specific comments to the reference REMARK instead.");
  }

  if (bsp != NULL && bsp->repr == Seq_repr_seg) {
    slp = SeqLocFindNext (sfp->location, NULL);
    while (slp != NULL) {
      sip = SeqLocId (slp);
      if (sip != NULL) {
        if (SeqIdIn (sip, bsp->id)) {
          on_seg = TRUE;
        }
      }
      slp = SeqLocFindNext (sfp->location, slp);
    }
    if (on_seg) {
      sev = SEV_ERROR;
      if (is_nc) {
        sev = SEV_WARNING;
      }
      if (! DeltaOrFarSeg (vsp->sep, sfp->location)) {
        ValidErr (vsp, sev, ERR_SEQ_FEAT_LocOnSegmentedBioseq, "Feature location on segmented bioseq, not on parts");
      }
    }
  }

  if (sfp->idx.subtype == FEATDEF_intron) {
    no_nonconsensus_except = TRUE;
    if (sfp->excpt) {
      if (StringISearch (sfp->except_text, "nonconsensus splice site") != NULL) {
        no_nonconsensus_except = FALSE;
      }
    }
    if (no_nonconsensus_except) {
      ValidateIntronEndsAtSpliceSiteOrGap (vsp, sfp->location);
    }
  }

  if (gcp != NULL) {
    gcp->itemID = olditemid;
    gcp->thistype = olditemtype;
  }

  return TRUE;
}

static Boolean GeneSpansOrigin (SeqMgrFeatContextPtr context, Int4 bsplength)

{
  Int4Ptr  ivals;

  if (context == NULL || bsplength < 1) return FALSE;
  ivals = context->ivals;
  if (ivals == NULL || context->numivals != 2) return FALSE;
  if (context->strand == Seq_strand_minus) {
    if (ivals [1] == 0 && ivals [2] == bsplength - 1) return TRUE;
  } else {
    if (ivals [2] == 0 && ivals [1] == bsplength - 1) return TRUE;
  }
  return FALSE;
}

static void CheckMultiIntervalGene (SeqFeatPtr sfp, SeqMgrFeatContextPtr context, ValidStructPtr vsp, GatherContextPtr gcp)

{
  BioseqPtr     bsp;
  Int4          count;
  SeqLocPtr     mappedloc = NULL;
  Uint2         olditemtype = 0;
  Uint4         olditemid = 0;
  Boolean       segmented = FALSE;
  ErrSev        sev = /* SEV_ERROR */ SEV_WARNING;
  SeqIdPtr      sip;
  SeqLocPtr     slp;
  TextSeqIdPtr  tsip;

  if (sfp == NULL || context == NULL || vsp == NULL) return;
  if (context->numivals < 2) return;

  if (sfp->excpt) {
    if (StringISearch (sfp->except_text, "trans-splicing") != NULL) return;
  }

  if (SeqLocId (sfp->location) == NULL) {
    bsp = context->bsp;
    if (bsp == NULL || bsp->repr != Seq_repr_seg) return;
    mappedloc = SeqLocMerge (bsp, sfp->location, NULL, FALSE, TRUE, FALSE);
    if (mappedloc == NULL) return;
    count = 0;
    slp = SeqLocFindNext (mappedloc, NULL);
    while (slp != NULL) {
      count++;
      slp = SeqLocFindNext (mappedloc, slp);
    }
    SeqLocFree (mappedloc);
    if (count < 2) return;
    segmented = TRUE;
  }

  bsp = context->bsp;
  if (bsp != NULL) {
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_OTHER) {
        tsip = (TextSeqIdPtr) sip->data.ptrvalue;
        if (tsip != NULL && tsip->accession != NULL) {
          if (StringNICmp (tsip->accession, "NC_", 3) == 0) {
            sev = SEV_WARNING;
          }
        }
      } else if (sip->choice == SEQID_EMBL || sip->choice == SEQID_DDBJ) {
        sev = SEV_WARNING;
      }
    }
    if (bsp->topology == 2) {
      if (context->numivals == 2 && GeneSpansOrigin (context, bsp->length)) return;
      sev = SEV_WARNING;
    }
  }

  if (gcp != NULL) {
    olditemid = gcp->itemID;
    olditemtype = gcp->thistype;
    gcp->itemID = context->itemID;
    gcp->thistype = OBJ_SEQFEAT;
  }

  vsp->sfp = sfp;
  if (segmented) {
    ValidErr (vsp, sev, ERR_SEQ_FEAT_SegmentedGeneProblem,
              "Gene feature on segmented sequence should cover all bases within its extremes");
  } else if (vsp->is_small_genome_set) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_MultiIntervalGene,
              "Multiple interval gene feature in small genome set - set trans-splicing exception if appropriate");
  } else {
    ValidErr (vsp, sev, ERR_SEQ_FEAT_MultiIntervalGene,
              "Gene feature on non-segmented sequence should not have multiple intervals");
  }
  vsp->sfp = NULL;

  if (gcp != NULL) {
    gcp->itemID = olditemid;
    gcp->thistype = olditemtype;
  }
}

static Boolean LIBCALLBACK ValidateSeqFeatIndexed (SeqFeatPtr sfp, SeqMgrFeatContextPtr context)
{
  ValidStructPtr  vsp;
  BioseqValidStrPtr bvsp;

  bvsp = (BioseqValidStrPtr) context->userdata;
  vsp = bvsp->vsp;

  if (sfp->data.choice == SEQFEAT_GENE) {
    CheckMultiIntervalGene (sfp, context, vsp, vsp->gcp);
  }

  return ValidateSeqFeatCommon (sfp, bvsp, vsp, context->left, context->right, context->numivals, context->itemID, context->farloc, context->bsp);
}

static void ValidateSeqFeatContext (GatherContextPtr gcp)
{
  ValidStructPtr  vsp;
  BioseqValidStrPtr bvsp;
  SeqFeatPtr      sfp;

  bvsp = (BioseqValidStrPtr) (gcp->userdata);
  vsp = bvsp->vsp;
  sfp = (SeqFeatPtr) (gcp->thisitem);

  ValidateSeqFeatCommon (sfp, bvsp, vsp, gcp->extremes.left, gcp->extremes.right, 0, 0, FALSE, NULL);
}

/*****************************************************************************
*
*   CountryIsValid(name)
*      Validates subsource country against official country names
*
*****************************************************************************/

static CharPtr  Nlm_valid_country_codes [] = {
  "Afghanistan",
  "Albania",
  "Algeria",
  "American Samoa",
  "Andorra",
  "Angola",
  "Anguilla",
  "Antarctica",
  "Antigua and Barbuda",
  "Arctic Ocean",
  "Argentina",
  "Armenia",
  "Aruba",
  "Ashmore and Cartier Islands",
  "Atlantic Ocean",
  "Australia",
  "Austria",
  "Azerbaijan",
  "Bahamas",
  "Bahrain",
  "Baker Island",
  "Baltic Sea",
  "Bangladesh",
  "Barbados",
  "Bassas da India",
  "Belarus",
  "Belgium",
  "Belize",
  "Benin",
  "Bermuda",
  "Bhutan",
  "Bolivia",
  "Borneo",
  "Bosnia and Herzegovina",
  "Botswana",
  "Bouvet Island",
  "Brazil",
  "British Virgin Islands",
  "Brunei",
  "Bulgaria",
  "Burkina Faso",
  "Burundi",
  "Cambodia",
  "Cameroon",
  "Canada",
  "Cape Verde",
  "Cayman Islands",
  "Central African Republic",
  "Chad",
  "Chile",
  "China",
  "Christmas Island",
  "Clipperton Island",
  "Cocos Islands",
  "Colombia",
  "Comoros",
  "Cook Islands",
  "Coral Sea Islands",
  "Costa Rica",
  "Cote d'Ivoire",
  "Croatia",
  "Cuba",
  "Cyprus",
  "Czech Republic",
  "Democratic Republic of the Congo",
  "Denmark",
  "Djibouti",
  "Dominica",
  "Dominican Republic",
  "East Timor",
  "Ecuador",
  "Egypt",
  "El Salvador",
  "Equatorial Guinea",
  "Eritrea",
  "Estonia",
  "Ethiopia",
  "Europa Island",
  "Falkland Islands (Islas Malvinas)",
  "Faroe Islands",
  "Fiji",
  "Finland",
  "France",
  "French Guiana",
  "French Polynesia",
  "French Southern and Antarctic Lands",
  "Gabon",
  "Gambia",
  "Gaza Strip",
  "Georgia",
  "Germany",
  "Ghana",
  "Gibraltar",
  "Glorioso Islands",
  "Greece",
  "Greenland",
  "Grenada",
  "Guadeloupe",
  "Guam",
  "Guatemala",
  "Guernsey",
  "Guinea",
  "Guinea-Bissau",
  "Guyana",
  "Haiti",
  "Heard Island and McDonald Islands",
  "Honduras",
  "Hong Kong",
  "Howland Island",
  "Hungary",
  "Iceland",
  "India",
  "Indian Ocean",
  "Indonesia",
  "Iran",
  "Iraq",
  "Ireland",
  "Isle of Man",
  "Israel",
  "Italy",
  "Jamaica",
  "Jan Mayen",
  "Japan",
  "Jarvis Island",
  "Jersey",
  "Johnston Atoll",
  "Jordan",
  "Juan de Nova Island",
  "Kazakhstan",
  "Kenya",
  "Kerguelen Archipelago",
  "Kingman Reef",
  "Kiribati",
  "Kosovo",
  "Kuwait",
  "Kyrgyzstan",
  "Laos",
  "Latvia",
  "Lebanon",
  "Lesotho",
  "Liberia",
  "Libya",
  "Liechtenstein",
  "Lithuania",
  "Luxembourg",
  "Macau",
  "Macedonia",
  "Madagascar",
  "Malawi",
  "Malaysia",
  "Maldives",
  "Mali",
  "Malta",
  "Marshall Islands",
  "Martinique",
  "Mauritania",
  "Mauritius",
  "Mayotte",
  "Mediterranean Sea",
  "Mexico",
  "Micronesia",
  "Midway Islands",
  "Moldova",
  "Monaco",
  "Mongolia",
  "Montenegro",
  "Montserrat",
  "Morocco",
  "Mozambique",
  "Myanmar",
  "Namibia",
  "Nauru",
  "Navassa Island",
  "Nepal",
  "Netherlands",
  "Netherlands Antilles",
  "New Caledonia",
  "New Zealand",
  "Nicaragua",
  "Niger",
  "Nigeria",
  "Niue",
  "Norfolk Island",
  "North Korea",
  "North Sea",
  "Northern Mariana Islands",
  "Norway",
  "Oman",
  "Pacific Ocean",
  "Pakistan",
  "Palau",
  "Palmyra Atoll",
  "Panama",
  "Papua New Guinea",
  "Paracel Islands",
  "Paraguay",
  "Peru",
  "Philippines",
  "Pitcairn Islands",
  "Poland",
  "Portugal",
  "Puerto Rico",
  "Qatar",
  "Republic of the Congo",
  "Reunion",
  "Romania",
  "Ross Sea",
  "Russia",
  "Rwanda",
  "Saint Helena",
  "Saint Kitts and Nevis",
  "Saint Lucia",
  "Saint Pierre and Miquelon",
  "Saint Vincent and the Grenadines",
  "Samoa",
  "San Marino",
  "Sao Tome and Principe",
  "Saudi Arabia",
  "Senegal",
  "Serbia",
  "Seychelles",
  "Sierra Leone",
  "Singapore",
  "Slovakia",
  "Slovenia",
  "Solomon Islands",
  "Somalia",
  "South Africa",
  "South Georgia and the South Sandwich Islands",
  "South Korea",
  "South Sudan",
  "Southern Ocean",
  "Spain",
  "Spratly Islands",
  "Sri Lanka",
  "Sudan",
  "Suriname",
  "Svalbard",
  "Swaziland",
  "Sweden",
  "Switzerland",
  "Syria",
  "Taiwan",
  "Tajikistan",
  "Tanzania",
  "Tasman Sea",
  "Thailand",
  "Togo",
  "Tokelau",
  "Tonga",
  "Trinidad and Tobago",
  "Tromelin Island",
  "Tunisia",
  "Turkey",
  "Turkmenistan",
  "Turks and Caicos Islands",
  "Tuvalu",
  "Uganda",
  "Ukraine",
  "United Arab Emirates",
  "United Kingdom",
  "Uruguay",
  "USA",
  "Uzbekistan",
  "Vanuatu",
  "Venezuela",
  "Viet Nam",
  "Virgin Islands",
  "Wake Island",
  "Wallis and Futuna",
  "West Bank",
  "Western Sahara",
  "Yemen",
  "Zambia",
  "Zimbabwe",
  NULL
};

static CharPtr  Nlm_formerly_valid_country_codes [] = {
  "Belgian Congo",
  "British Guiana",
  "Burma",
  "Czechoslovakia",
  "Korea",
  "Serbia and Montenegro",
  "Siam",
  "USSR",
  "Yugoslavia",
  "Zaire",
  NULL
};

NLM_EXTERN CharPtr PNTR GetValidCountryList (void)

{
  return (CharPtr PNTR) Nlm_valid_country_codes;
}

NLM_EXTERN Boolean CountryIsValid (CharPtr name, BoolPtr old_countryP, BoolPtr bad_capP)
{
  Int2     L, R, mid;
  CharPtr  ptr;
  Char     str [256];

  if (StringHasNoText (name)) return FALSE;

  StringNCpy_0 (str, name, sizeof (str));
  ptr = StringChr (str, ':');
  if (ptr != NULL) {
    *ptr = '\0';
  }

  L = 0;
  R = sizeof (Nlm_valid_country_codes) / sizeof (Nlm_valid_country_codes [0]) - 1; /* -1 because now NULL terminated */

  while (L < R) {
    mid = (L + R) / 2;
    if (StringICmp (Nlm_valid_country_codes [mid], str) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  if (StringICmp (Nlm_valid_country_codes [R], str) == 0) {
    if (bad_capP != NULL) {
      if (StringCmp (Nlm_valid_country_codes [R], str) != 0) {
        *bad_capP = TRUE;
      }
    }
    return TRUE;
  }

  L = 0;
  R = sizeof (Nlm_formerly_valid_country_codes) / sizeof (Nlm_formerly_valid_country_codes [0]) - 1; /* -1 because now NULL terminated */

  while (L < R) {
    mid = (L + R) / 2;
    if (StringICmp (Nlm_formerly_valid_country_codes [mid], str) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  if (StringICmp (Nlm_formerly_valid_country_codes [R], str) == 0) {
    if (old_countryP != NULL) {
      *old_countryP = TRUE;
    }
    if (bad_capP != NULL) {
      if (StringCmp (Nlm_formerly_valid_country_codes [R], str) != 0) {
        *bad_capP = TRUE;
      }
    }
    return FALSE;
  }

  return FALSE;
}


NLM_EXTERN CharPtr GetCorrectedCountryCapitalization (CharPtr name)
{
  Int2     L, R, mid;
  CharPtr  ptr;
  Char     str [256];

  if (StringHasNoText (name)) return NULL;

  StringNCpy_0 (str, name, sizeof (str));
  ptr = StringChr (str, ':');
  if (ptr != NULL) {
    *ptr = '\0';
  }

  L = 0;
  R = sizeof (Nlm_valid_country_codes) / sizeof (Nlm_valid_country_codes [0]) - 1; /* -1 because now NULL terminated */

  while (L < R) {
    mid = (L + R) / 2;
    if (StringICmp (Nlm_valid_country_codes [mid], str) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  if (StringICmp (Nlm_valid_country_codes [R], str) == 0) {
    return Nlm_valid_country_codes[R];
  }

  return NULL;
}

static CharPtr bodiesOfWater [] = {
  "Basin",
  "Bay",
  "Bight",
  "Canal",
  "Channel",
  "Coastal",
  "Cove",
  "Estuary",
  "Fjord",
  "Freshwater",
  "Gulf",
  "Harbor",
  "Inlet",
  "Lagoon",
  "Lake",
  "Narrows",
  "Ocean",
  "Offshore",
  "Passage",
  "Passages",
  "Reef",
  "River",
  "Sea",
  "Seawater",
  "Sound",
  "Strait",
  "Trench",
  "Trough",
  "Water",
  "Waters",
  NULL
};

static TextFsaPtr GetBodiesOfWaterFSA (void)


{
  TextFsaPtr  fsa;
  Int2        i;
  CharPtr     prop = "BodiesOfWaterFSA";

  fsa = (TextFsaPtr) GetAppProperty (prop);
  if (fsa != NULL) return fsa;

  fsa = TextFsaNew ();
  if (fsa != NULL) {
    for (i = 0; bodiesOfWater [i] != NULL; i++) {
      TextFsaAdd (fsa, bodiesOfWater [i]);
    }
  }

  SetAppProperty (prop, (Pointer) fsa);

  return fsa;
}

NLM_EXTERN Boolean StringContainsBodyOfWater (CharPtr str)

{
  Char        ch;
  TextFsaPtr  fsa;
  CharPtr     ptr;
  Int4        state;
  ValNodePtr  matches;

  if (StringHasNoText (str)) return FALSE;

  fsa = GetBodiesOfWaterFSA ();
  if (fsa == NULL) return FALSE;

  state = 0;
  ptr = str;
  ch = *ptr;

  while (ch != '\0') {
    matches = NULL;
    state = TextFsaNext (fsa, state, ch, &matches);
    ptr++;
    ch = *ptr;
    if (ch == '\0' || ch == ',' || ch == ':' || ch == ';' || ch == ' ') {
      if (matches != NULL) return TRUE;
      state = 0;
    }
  }

  return FALSE;
}

/* BEGINNING OF NEW LATITUDE-LONGITUDE COUNTRY VALIDATION CODE */

/* latitude-longitude to country conversion */

typedef struct ctyblock {
  CharPtr  name;    /* name of country or country: subregion */
  CharPtr  level0;  /* just the country */
  CharPtr  level1;  /* just the subregion */
  Int4     area;    /* pixel area for choosing smallest overlapping subregion */
  Int4     minlat;  /* minimum latitude */
  Int4     maxlat;  /* maximum latitude */
  Int4     minlon;  /* minimum longitude */
  Int4     maxlon;  /* maximum longitude */
} CtyBlock, PNTR CtyBlockPtr;

typedef struct latblock {
  CtyBlockPtr  landmass;   /* points to instance in countries list */
  Int4         lat;        /* latitude (integer in 10ths of a degree) */
  Int4         minlon;     /* minimum longitude */
  Int4         maxlon;     /* maximum longitude */
} LatBlock, PNTR LatBlockPtr;

typedef struct ctryset {
  ValNodePtr        ctyblocks;      /* linked list of country blocks */
  CtyBlockPtr PNTR  ctyarray;       /* country blocks sorted by name */
  Int4              numCtyBlocks;
  ValNodePtr        latblocks;      /* linked list of latitude blocks */
  LatBlockPtr PNTR  latarray;       /* latitude blocks sorted by latitude then longitude */
  Int4              numLatBlocks;
  FloatHi           scale;
} CtrySet, PNTR CtrySetPtr;

static int LIBCALLBACK SortByCountry (
  VoidPtr ptr1,
  VoidPtr ptr2
)

{
  CtyBlockPtr  cbp1;
  CtyBlockPtr  cbp2;
  int          cmp;
  ValNodePtr   vnp1;
  ValNodePtr   vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  cbp1 = (CtyBlockPtr) vnp1->data.ptrvalue;
  cbp2 = (CtyBlockPtr) vnp2->data.ptrvalue;
  if (cbp1 == NULL || cbp2 == NULL) return 0;

  cmp = StringICmp (cbp1->name, cbp2->name);
  if (cmp > 0) {
    return 1;
  } else if (cmp < 0) {
    return -1;
  }

  return 0;
}

static int LIBCALLBACK SortByLatLon (
  VoidPtr ptr1,
  VoidPtr ptr2
)

{
  CtyBlockPtr  cbp1;
  CtyBlockPtr  cbp2;
  int          cmp;
  LatBlockPtr  lbp1;
  LatBlockPtr  lbp2;
  ValNodePtr   vnp1;
  ValNodePtr   vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  lbp1 = (LatBlockPtr) vnp1->data.ptrvalue;
  lbp2 = (LatBlockPtr) vnp2->data.ptrvalue;
  if (lbp1 == NULL || lbp2 == NULL) return 0;

  if (lbp1->lat < lbp2->lat) {
    return -1;
  } else if (lbp1->lat > lbp2->lat) {
    return 1;
  }

  if (lbp1->minlon < lbp2->minlon) {
    return -1;
  } else if (lbp1->minlon > lbp2->minlon) {
    return 1;
  }

  if (lbp1->maxlon < lbp2->maxlon) {
    return 1;
  } else if (lbp1->maxlon > lbp2->maxlon) {
    return -1;
  }

  cbp1 = lbp1->landmass;
  cbp2 = lbp2->landmass;
  if (cbp1 == NULL || cbp2 == NULL) return 0;

  if (cbp1->area < cbp2->area) {
    return -1;
  } else if (cbp1->area > cbp2->area) {
    return 1;
  }

  cmp = StringICmp (cbp1->name, cbp2->name);
  if (cmp > 0) {
    return 1;
  } else if (cmp < 0) {
    return -1;
  }

  return 0;
}

#define EPSILON 0.001

static Int4 ConvertLat (FloatHi lat, FloatHi scale) {

  Int4  val = 0;

  if (lat < -90.0) {
    lat = -90.0;
  }
  if (lat > 90.0) {
    lat = 90.0;
  }

  if (lat > 0) {
    val = (Int4) (lat * scale + EPSILON);
  } else {
    val = (Int4) (-(-lat * scale + EPSILON));
  }

  return val;
}

static Int4 ConvertLon (FloatHi lon, FloatHi scale) {

  Int4  val = 0;

  if (lon < -180.0) {
    lon = -180.0;
  }
  if (lon > 180.0) {
    lon = 180.0;
  }

  if (lon > 0) {
    val = (Int4) (lon * scale + EPSILON);
  } else {
    val = (Int4) (-(-lon * scale + EPSILON));
  }

  return val;
}

static CtrySetPtr FreeLatLonCountryData (
  CtrySetPtr csp
)

{
  CtyBlockPtr  cbp;
  ValNodePtr   vnp;

  if (csp == NULL) return NULL;

  for (vnp = csp->ctyblocks; vnp != NULL; vnp = vnp->next) {
    cbp = (CtyBlockPtr) vnp->data.ptrvalue;
    if (cbp == NULL) continue;
    MemFree (cbp->name);
    MemFree (cbp->level0);
    MemFree (cbp->level1);
  }

  ValNodeFreeData (csp->ctyblocks);
  ValNodeFreeData (csp->latblocks);

  MemFree (csp->ctyarray);
  MemFree (csp->latarray);

  MemFree (csp);

  return NULL;
}

/* Original data source is Natural Earth.  Free vector and raster map data @ http://naturalearthdata.com */

static CharPtr LatLonCountryReadNextLine (
  FileCache PNTR fcp,
  CharPtr buf,
  size_t bufsize,
  CharPtr PNTR local,
  Int4Ptr idxP
)

{
  Int4     idx;
  CharPtr  str = NULL;

  if (fcp != NULL) {
    str = FileCacheReadLine (fcp, buf, bufsize, NULL);
  }

  if (local != NULL && idxP != NULL) {
    idx = *idxP;
    str = local [idx];
    if (str != NULL) {
      StringNCpy_0 (buf, local [idx], bufsize);
      str = buf;
    }
    idx++;
    *idxP = idx;
  }

  return str;
}

static CtrySetPtr ReadLatLonCountryData (
  CharPtr prop,
  CharPtr file,
  CharPtr PNTR local
)

{
  Char              buf [128];
  Char              ch;
  CtyBlockPtr       cbp = NULL;
  CtrySetPtr        csp = NULL;
  CtyBlockPtr PNTR  ctyarray;
  ValNodePtr        ctyblocks = NULL;
  FileCache         fc;
  FileCache PNTR    fcp = NULL;
  FILE              *fp = NULL;
  Int4              i;
  Int4              idx = 0;
  ValNodePtr        lastlatblock = NULL;
  ValNodePtr        lastctyblock = NULL;
  FloatHi           latitude;
  LatBlockPtr PNTR  latarray;
  ValNodePtr        latblocks = NULL;
  LatBlockPtr       lbp;
  Char              line [1024];
  FloatHi           maxlongitude;
  FloatHi           minlongitude;
  Char              path [PATH_MAX];
  CharPtr           ptr;
  CharPtr           recentCountry = NULL;
  FloatHi           scale = 0.0;
  Boolean           scale_not_set = TRUE;
  ErrSev            sev;
  CharPtr           str;
  Char              tmp [128];
  double            val;
  ValNodePtr        vnp;
  CharPtr           wrk;

  if (FindPath ("ncbi", "ncbi", "data", path, sizeof (path))) {
    FileBuildPath (path, NULL, file);
    sev = ErrSetMessageLevel (SEV_ERROR);
    fp = FileOpen (path, "r");
    ErrSetMessageLevel (sev);
  }

  if (fp != NULL) {
    FileCacheSetup (&fc, fp);
    fcp = &fc;
    local = NULL;
  } else if (local == NULL) {
    return NULL;
  }

  for (str = LatLonCountryReadNextLine (fcp, line, sizeof (line), local, &idx);
       str != NULL;
       str = LatLonCountryReadNextLine (fcp, line, sizeof (line), local, &idx)) {
    if (StringHasNoText (str)) continue;

    /* if reading from local copy, str cannot be modified, so copy to local buf and reset pointer */

    StringNCpy_0 (buf, str, sizeof (buf));
    str = buf;

    ch = str [0];

    /* ignore comment lines starting with hyphen */

    if (ch == '-') continue;

    /* Scale should be at top of file, after comments */

    if (IS_DIGIT (ch)) {
      if (scale_not_set && sscanf (str, "%lf", &val) == 1) {
        scale = (FloatHi) val;
        scale_not_set = FALSE;
      }

      continue;
    }

    /* Country starts on first column */

    if (IS_ALPHA (ch)) {

      if (scale_not_set) {
        scale = 20.0;
        scale_not_set = FALSE;
      }

      ptr = StringChr (str, '\t');
      if (ptr != NULL) {
        *ptr = '\0';
      }

      if (StringCmp (str, recentCountry) == 0) continue;

      cbp = (CtyBlockPtr) MemNew (sizeof (CtyBlock));
      if (cbp == NULL) continue;

      TrimSpacesAroundString (str);
      cbp->name = StringSave (str);
      StringNCpy_0 (tmp, str, sizeof (tmp));
      ptr = StringChr (tmp, ':');
      if (ptr != NULL) {
        *ptr = '\0';
        ptr++;
        TrimSpacesAroundString (ptr);
        if (StringDoesHaveText (ptr)) {
          cbp->level1 = StringSave (ptr);
        }
        TrimSpacesAroundString (tmp);
        cbp->level0 = StringSave (tmp);
      } else {
        TrimSpacesAroundString (str);
        cbp->level0 = StringSave (str);
      }
      cbp->area = 0;
      cbp->minlat = INT4_MAX;
      cbp->maxlat = INT4_MIN;
      cbp->minlon = INT4_MAX;
      cbp->maxlon = INT4_MIN;
      vnp = ValNodeAddPointer (&lastctyblock, 0, (Pointer) cbp);
      if (ctyblocks == NULL) {
        ctyblocks = vnp;
      }
      lastctyblock = vnp;

      recentCountry = cbp->name;

      continue;
    }

    /* Latitude with longitude min/max pairs on line starting with tab */

    if (ch != '\t') continue;

    wrk = StringSave (str + 1);
    if (wrk == NULL) continue;

    ptr = StringChr (wrk, '\t');
    if (ptr != NULL) {
      *ptr = '\0';
      ptr++;
      if (sscanf (wrk, "%lf", &val) == 1) {
        latitude = (FloatHi) val;

        str = ptr;
        while (StringDoesHaveText (str)) {
          ptr = StringChr (str, '\t');
          if (ptr != NULL) {
            *ptr = '\0';
            ptr++;
          }
          if (sscanf (str, "%lf", &val) != 1) {
            /* prevent infinite loop if it fails */
            str = NULL;
          } else {
            minlongitude = (FloatHi) val;
            str = ptr;
            if (StringDoesHaveText (str)) {
              ptr = StringChr (str, '\t');
              if (ptr != NULL) {
                *ptr = '\0';
                ptr++;
              }
              if (sscanf (str, "%lf", &val) == 1) {
                maxlongitude = (FloatHi) val;

                lbp = (LatBlockPtr) MemNew (sizeof (LatBlock));
                if (lbp != NULL) {
                  lbp->landmass = cbp;
                  lbp->lat = ConvertLat (latitude, scale);
                  lbp->minlon = ConvertLon (minlongitude, scale);
                  lbp->maxlon = ConvertLon (maxlongitude, scale);

                  vnp = ValNodeAddPointer (&lastlatblock, 0, (Pointer) lbp);
                  if (latblocks == NULL) {
                    latblocks = vnp;
                  }
                  lastlatblock = vnp;
                }
              }
            }
            str = ptr;
          }
        }
      }
    }

    MemFree (wrk);
  }

  if (fp != NULL) {
    FileClose (fp);
  }

  if (ctyblocks == NULL || latblocks == NULL) {
    return NULL;
  }

  csp = (CtrySetPtr) MemNew (sizeof (CtrySet));
  if (csp == NULL) return NULL;

  for (vnp = latblocks; vnp != NULL; vnp = vnp->next) {
    lbp = (LatBlockPtr) vnp->data.ptrvalue;
    if (lbp == NULL) continue;
    cbp = lbp->landmass;
    if (cbp == NULL) continue;
    cbp->area += lbp->maxlon - lbp->minlon + 1;
    if (cbp->minlat > lbp->lat) {
      cbp->minlat = lbp->lat;
    }
    if (cbp->maxlat < lbp->lat) {
      cbp->maxlat = lbp->lat;
    }
    if (cbp->minlon > lbp->minlon) {
      cbp->minlon = lbp->minlon;
    }
    if (cbp->maxlon < lbp->maxlon) {
      cbp->maxlon = lbp->maxlon;
    }
  }

  ctyblocks = ValNodeSort (ctyblocks, SortByCountry);
  csp->ctyblocks = ctyblocks;
  csp->numCtyBlocks = ValNodeLen (ctyblocks);

  latblocks = ValNodeSort (latblocks, SortByLatLon);
  csp->latblocks = latblocks;
  csp->numLatBlocks = ValNodeLen (latblocks);

  if (scale_not_set) {
    scale = 20.0;
  }
  csp->scale = scale;

  ctyarray = (CtyBlockPtr PNTR) MemNew (sizeof (CtyBlockPtr) * (csp->numCtyBlocks + 1));
  if (ctyarray != NULL) {
    for (vnp = ctyblocks, i = 0; vnp != NULL; vnp = vnp->next, i++) {
      cbp = (CtyBlockPtr) vnp->data.ptrvalue;
      ctyarray [i] = cbp;
    }

    csp->ctyarray = ctyarray;
  }

  latarray = (LatBlockPtr PNTR) MemNew (sizeof (LatBlockPtr) * (csp->numLatBlocks + 1));
  if (latarray != NULL) {
    for (vnp = latblocks, i = 0; vnp != NULL; vnp = vnp->next, i++) {
      lbp = (LatBlockPtr) vnp->data.ptrvalue;
      latarray [i] = lbp;
    }

    csp->latarray = latarray;
  }

/*
{
  FILE *fp;
  fp = FileOpen ("ctrymap.txt", "w");
  if (fp != NULL) {
    for (vnp = latblocks; vnp != NULL; vnp = vnp->next) {
      lbp = (LatBlockPtr) vnp->data.ptrvalue;
      if (lbp == NULL) continue;
      cbp = lbp->landmass;
      if (cbp == NULL) continue;
      fprintf (fp, "%s\t[%d]\t%d\t%d\t%d\n", cbp->name, (int) cbp->area,
               (int) lbp->lat, (int) lbp->minlon, (int) lbp->maxlon);
    }
    FileClose (fp);
  }
}
*/

  return csp;
}

static Boolean ctryset_not_found = FALSE;
static Boolean watrset_not_found = FALSE;

extern CharPtr latlon_onedegree [];
extern CharPtr water_onedegree [];

static CtrySetPtr GetLatLonCountryData (void)

{
  CtrySetPtr  csp = NULL;
  CharPtr     prop = "CountryLatLonData";

  csp = (CtrySetPtr) GetAppProperty (prop);
  if (csp != NULL) return csp;

  if (ctryset_not_found) return NULL;

  csp = ReadLatLonCountryData (prop, "lat_lon_country.txt", latlon_onedegree);

  if (csp == NULL) {
    ctryset_not_found = TRUE;
    return NULL;
  }

  SetAppProperty (prop, (Pointer) csp);

  return csp;
}

static CtrySetPtr GetLatLonWaterData (void)

{
  CtrySetPtr  csp = NULL;
  CharPtr     prop = "WaterLatLonData";

  csp = (CtrySetPtr) GetAppProperty (prop);
  if (csp != NULL) return csp;

  if (watrset_not_found) return NULL;

  csp = ReadLatLonCountryData (prop, "lat_lon_water.txt", water_onedegree);

  if (csp == NULL) {
    watrset_not_found = TRUE;
    return NULL;
  }

  SetAppProperty (prop, (Pointer) csp);

  return csp;
}

static CtyBlockPtr GetEntryInLatLonListIndex (
  CharPtr country,
  CtrySetPtr csp
)

{
  CtyBlockPtr PNTR  array;
  CtyBlockPtr       cbp;
  Int2              L, R, mid;

  if (StringHasNoText (country)) return NULL;
  if (csp == NULL) return NULL;

  array = csp->ctyarray;
  if (array == NULL) return NULL;

  L = 0;
  R = csp->numCtyBlocks - 1;

  while (L < R) {
    mid = (L + R) / 2;
    cbp = array [mid];
    if (cbp != NULL && cbp->name != NULL && StringICmp (cbp->name, country) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  cbp = array [R];
  if (cbp != NULL && cbp->name != NULL && StringICmp (cbp->name, country) == 0) return cbp;

  return NULL;
}

NLM_EXTERN Boolean CountryIsInLatLonList (
  CharPtr country
)

{
  CtyBlockPtr  cbp;
  CtrySetPtr   csp;

  if (StringHasNoText (country)) return FALSE;
  csp = GetLatLonCountryData ();
  if (csp == NULL) return FALSE;

  cbp = GetEntryInLatLonListIndex (country, csp);
  if (cbp != NULL && cbp->name != NULL && StringICmp (cbp->name, country) == 0) return TRUE;

  return FALSE;
}

NLM_EXTERN Boolean IsCountryInLatLonList (
  CharPtr country
)

{
  return CountryIsInLatLonList (country);
}

NLM_EXTERN Boolean WaterIsInLatLonList (
  CharPtr country
)

{
  CtyBlockPtr  cbp;
  CtrySetPtr   csp;

  if (StringHasNoText (country)) return FALSE;
  csp = GetLatLonWaterData ();
  if (csp == NULL) return FALSE;

  cbp = GetEntryInLatLonListIndex (country, csp);
  if (cbp != NULL && cbp->name != NULL && StringICmp (cbp->name, country) == 0) return TRUE;

  return FALSE;
}

static int LatLonCmp (
  LatBlockPtr lbp,
  Int2 latitude
)

{
  if (lbp == NULL) return 0;

  if (lbp->lat < latitude) {
    return -1;
  } else if (lbp->lat > latitude) {
    return 1;
  }

  return 0;
}

static Int4 GetLatLonIndex (
  CtrySetPtr csp,
  LatBlockPtr PNTR array,
  Int2 latitude
)

{
  LatBlockPtr  lbp;
  Int4         L, R, mid;

  if (csp == NULL || array == NULL) return 0;

  L = 0;
  R = csp->numLatBlocks - 1;

  while (L < R) {
    mid = (L + R) / 2;
    lbp = array [mid];
    if (lbp != NULL && LatLonCmp (lbp, latitude) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  return R;
}

static Boolean SubregionStringICmp (
  CharPtr region,
  CharPtr country
)

{
  Char     possible [256];
  CharPtr  ptr;

  if (StringHasNoText (region) || StringHasNoText (country)) return FALSE;
  StringNCpy_0 (possible, region, sizeof (possible));
  ptr = StringChr (possible, ':');
  if (ptr == NULL) return FALSE;
  *ptr = '\0';
  if (StringICmp (possible, country) == 0) return TRUE;
  return FALSE;
}

static Boolean RegionContainsLatLon (
  CharPtr country,
  FloatHi lat,
  FloatHi lon,
  CtrySetPtr csp
)

{
  LatBlockPtr PNTR  array;
  CtyBlockPtr       cbp;
  Int4              latitude;
  Int4              longitude;
  LatBlockPtr       lbp;
  Int4              R;

  if (StringHasNoText (country)) return FALSE;
  if (csp == NULL) return FALSE;

  array = csp->latarray;
  if (array == NULL) return FALSE;

  latitude = ConvertLat (lat, csp->scale);
  longitude = ConvertLon (lon, csp->scale);

  for (R = GetLatLonIndex (csp, array, latitude); R < csp->numLatBlocks; R++) {
    lbp = array [R];
    if (lbp == NULL) break;
    if (latitude != lbp->lat) break;

    if (longitude < lbp->minlon) continue;
    if (longitude > lbp->maxlon) continue;

    cbp = lbp->landmass;
    if (cbp == NULL) continue;
    if (StringICmp (cbp->name, country) == 0) return TRUE;
    if (SubregionStringICmp (cbp->name, country)) return TRUE;
  }

  return FALSE;
}

NLM_EXTERN Boolean CountryContainsLatLon (
  CharPtr country,
  FloatHi lat,
  FloatHi lon
)

{
  CtrySetPtr  csp;

  if (StringHasNoText (country)) return FALSE;

  csp = GetLatLonCountryData ();
  if (csp == NULL) return FALSE;

  return RegionContainsLatLon (country, lat, lon, csp);
}

NLM_EXTERN Boolean TestLatLonForCountry (
  CharPtr country,
  FloatHi lat,
  FloatHi lon
)

{
  return CountryContainsLatLon (country, lat, lon);
}

NLM_EXTERN Boolean WaterContainsLatLon (
  CharPtr country,
  FloatHi lat,
  FloatHi lon
)

{
  CtrySetPtr  csp;

  if (StringHasNoText (country)) return FALSE;

  csp = GetLatLonWaterData ();
  if (csp == NULL) return FALSE;

  return RegionContainsLatLon (country, lat, lon, csp);
}

static Boolean NewLatLonCandidateIsBetter (
  CharPtr country,
  CharPtr province,
  CtyBlockPtr best,
  CtyBlockPtr cbp,
  Boolean newer_is_smaller
)

{
  if (cbp == NULL) return FALSE;
  if (best == NULL) return TRUE;

  /* if no preferred country, just look for smallest area */
  if (country == NULL) {
    return newer_is_smaller;
  }

  /* if match to preferred country */
  if (StringICmp (country, cbp->level0) == 0) {

    /* if best was not preferred country, take new match */
    if (StringICmp (country, best->level0) != 0) return TRUE;

    /* if match to preferred province */
    if (province != NULL && StringICmp (province, cbp->level1) == 0) {

      /* if best was not preferred province, take new match */
      if (StringICmp (province, best->level1) != 0) return TRUE;
    }

    /* if both match province, or neither does, or no preferred province, take smallest */
    return newer_is_smaller;
  }

  /* if best matches preferred country, keep */
  if (StringICmp (country, best->level0) == 0) return FALSE;

  /* otherwise take smallest */
  return newer_is_smaller;
}

static CtyBlockPtr LookupRegionByLatLon (
  FloatHi lat,
  FloatHi lon,
  CharPtr country,
  CharPtr province,
  CtrySetPtr csp
)

{
  LatBlockPtr PNTR  array;
  CtyBlockPtr       cbp, best = NULL;
  Int4              latitude;
  Int4              longitude;
  LatBlockPtr       lbp;
  Int4              R;

  if (csp == NULL) return NULL;

  array = csp->latarray;
  if (array == NULL) return NULL;

  latitude = ConvertLat (lat, csp->scale);
  longitude = ConvertLon (lon, csp->scale);

  for (R = GetLatLonIndex (csp, array, latitude); R < csp->numLatBlocks; R++) {
    lbp = array [R];
    if (lbp == NULL) break;
    if (latitude != lbp->lat) break;

    if (longitude < lbp->minlon) continue;
    if (longitude > lbp->maxlon) continue;

    cbp = lbp->landmass;
    if (cbp == NULL) continue;

    if (best == NULL || NewLatLonCandidateIsBetter (country, province, best, cbp, (Boolean) (cbp->area < best->area))) {
      best = cbp;
    }
  }

  return best;
}

static CtyBlockPtr GuessCountryByLatLon (
  FloatHi lat,
  FloatHi lon,
  CharPtr country,
  CharPtr province
)

{
  CtrySetPtr  csp;

  csp = GetLatLonCountryData ();
  if (csp == NULL) return NULL;

  return LookupRegionByLatLon (lat, lon, country, province, csp);
}

static CtyBlockPtr GuessWaterByLatLon (
  FloatHi lat,
  FloatHi lon,
  CharPtr country
)

{
  CtrySetPtr  csp;

  csp = GetLatLonWaterData ();
  if (csp == NULL) return NULL;

  return LookupRegionByLatLon (lat, lon, country, NULL, csp);
}

NLM_EXTERN CharPtr LookupCountryByLatLon (
  FloatHi lat,
  FloatHi lon
)

{
  CtyBlockPtr  cbp;

  cbp = GuessCountryByLatLon (lat, lon, NULL, NULL);
  if (cbp == NULL) return NULL;

  return cbp->name;
}

NLM_EXTERN CharPtr GuessCountryForLatLon (
  FloatHi lat,
  FloatHi lon
)

{
  return LookupCountryByLatLon (lat, lon);
}

NLM_EXTERN CharPtr LookupWaterByLatLon (
  FloatHi lat,
  FloatHi lon
)

{
  CtyBlockPtr  cbp;

  cbp = GuessWaterByLatLon (lat, lon, NULL);
  if (cbp == NULL) return NULL;

  return cbp->name;
}

NLM_EXTERN FloatHi CountryDataScaleIs (void)

{
  CtrySetPtr  csp;

  csp = GetLatLonCountryData ();
  if (csp == NULL) return 0.0;

  return csp->scale;
}

NLM_EXTERN FloatHi WaterDataScaleIs (void)

{
  CtrySetPtr  csp;

  csp = GetLatLonWaterData ();
  if (csp == NULL) return 0.0;

  return csp->scale;
}


static Boolean RegionExtremesOverlap (
  CharPtr first,
  CharPtr second,
  CtrySetPtr csp
)

{
  CtyBlockPtr  cbp1, cbp2;

  if (StringHasNoText (first) || StringHasNoText (second)) return FALSE;
  if (csp == NULL) return FALSE;

  cbp1 = GetEntryInLatLonListIndex (first, csp);
  if (cbp1 == NULL || cbp1->name == NULL || StringICmp (cbp1->name, first) != 0) return FALSE;

  cbp2 = GetEntryInLatLonListIndex (second, csp);
  if (cbp2 == NULL || cbp2->name == NULL || StringICmp (cbp2->name, second) != 0) return FALSE;

  if (cbp1->minlat > cbp2->maxlat) return FALSE;
  if (cbp2->minlat > cbp1->maxlat) return FALSE;
  if (cbp1->minlon > cbp2->maxlon) return FALSE;
  if (cbp2->minlon > cbp1->maxlon) return FALSE;

  return TRUE;
}

NLM_EXTERN Boolean CountryExtremesOverlap (
  CharPtr first,
  CharPtr second
)

{
  CtrySetPtr  csp;

  if (StringHasNoText (first) || StringHasNoText (second)) return FALSE;
  csp = GetLatLonCountryData ();
  if (csp == NULL) return FALSE;

  return RegionExtremesOverlap (first, second, csp);
}

NLM_EXTERN Boolean CountryBoxesOverlap (
  CharPtr country1,
  CharPtr country2
)

{
  return CountryExtremesOverlap (country1, country2);
}

NLM_EXTERN Boolean WaterExtremesOverlap (
  CharPtr first,
  CharPtr second
)

{
  CtrySetPtr  csp;

  if (StringHasNoText (first) || StringHasNoText (second)) return FALSE;
  csp = GetLatLonWaterData ();
  if (csp == NULL) return FALSE;

  return RegionExtremesOverlap (first, second, csp);
}

/*
Distance on a spherical surface calculation adapted from
http://www.linuxjournal.com/magazine/
work-shell-calculating-distance-between-two-latitudelongitude-points
*/

#define EARTH_RADIUS 6371.0 /* average radius of non-spherical earth in kilometers */
#define CONST_PI 3.14159265359

static double DegreesToRadians (
  FloatHi degrees
)

{
  return (degrees * (CONST_PI / 180.0));
}

static FloatHi DistanceOnGlobe (
  FloatHi latA,
  FloatHi lonA,
  FloatHi latB,
  FloatHi lonB
)

{
  double lat1, lon1, lat2, lon2;
  double dLat, dLon, a, c;

  lat1 = DegreesToRadians (latA);
  lon1 = DegreesToRadians (lonA);
  lat2 = DegreesToRadians (latB);
  lon2 = DegreesToRadians (lonB);

  dLat = lat2 - lat1;
  dLon = lon2 - lon1;

   a = sin (dLat / 2) * sin (dLat / 2) +
       cos (lat1) * cos (lat2) * sin (dLon / 2) * sin (dLon / 2);
   c = 2 * atan2 (sqrt (a), sqrt (1 - a));

  return (FloatHi) (EARTH_RADIUS * c);
}

static FloatHi ErrorDistance (
  FloatHi latA,
  FloatHi lonA,
  FloatHi scale)
{
  double lat1, lon1, lat2, lon2;
  double dLat, dLon, a, c;

  lat1 = DegreesToRadians (latA);
  lon1 = DegreesToRadians (lonA);
  lat2 = DegreesToRadians (latA + (1.0 / scale));
  lon2 = DegreesToRadians (lonA + (1.0 / scale));

  dLat = lat2 - lat1;
  dLon = lon2 - lon1;

   a = sin (dLat / 2) * sin (dLat / 2) +
       cos (lat1) * cos (lat2) * sin (dLon / 2) * sin (dLon / 2);
   c = 2 * atan2 (sqrt (a), sqrt (1 - a));

  return (FloatHi) (EARTH_RADIUS * c);
  
}


static CtyBlockPtr RegionClosestToLatLon (
  FloatHi lat,
  FloatHi lon,
  FloatHi range,
  FloatHi PNTR distanceP,
  CtrySetPtr csp
)

{
  LatBlockPtr PNTR  array;
  CtyBlockPtr       cbp, best = NULL;
  FloatHi           closest = EARTH_RADIUS * CONST_PI * 2;
  FloatHi           delta;
  Int4              latitude;
  Int4              longitude;
  Int4              maxDelta;
  LatBlockPtr       lbp;
  Int4              R;
  Int4              x;
  Int4              y;
  Boolean           is_geographically_better;

  if (distanceP != NULL) {
    *distanceP = 0.0;
  }

  if (csp == NULL) return NULL;

  array = csp->latarray;
  if (array == NULL) return NULL;

  latitude = ConvertLat (lat, csp->scale);
  longitude = ConvertLon (lon, csp->scale);

  maxDelta = (Int4) (range * csp->scale + EPSILON);

  for (R = GetLatLonIndex (csp, array, latitude - maxDelta); R < csp->numLatBlocks; R++) {
    lbp = array [R];
    if (lbp == NULL) break;
    if (latitude + maxDelta < lbp->lat) break;

    if (longitude < lbp->minlon - maxDelta) continue;
    if (longitude > lbp->maxlon + maxDelta) continue;

    cbp = lbp->landmass;
    if (cbp == NULL) continue;

    if (longitude < lbp->minlon) {
      x = lbp->minlon;
    } else if (longitude > lbp->maxlon) {
      x = lbp->maxlon;
    } else {
      x = longitude;
    }

    y = lbp->lat;

    delta = DistanceOnGlobe (lat, lon, (FloatHi) (y / csp->scale), (FloatHi) (x / csp->scale));

    is_geographically_better = FALSE;
    if (delta < closest) {
      is_geographically_better = TRUE;
    } else if (delta - closest < 0.000001) {
      if (best == NULL || cbp->area < best->area) {
        is_geographically_better = TRUE;
      }
    }

    if (best == NULL || NewLatLonCandidateIsBetter (NULL, NULL, best, cbp, is_geographically_better)) {
      best = cbp;
      closest = delta;
    }
  }

  if (best != NULL) {
    if (distanceP != NULL) {
      *distanceP = closest;
    }
  }

  return best;
}

static CtyBlockPtr NearestCountryByLatLon (
  FloatHi lat,
  FloatHi lon,
  FloatHi range,
  FloatHi PNTR distanceP
)

{
  CtrySetPtr  csp;

  csp = GetLatLonCountryData ();
  if (csp == NULL) return NULL;

  return RegionClosestToLatLon (lat, lon, range, distanceP, csp);
}

static CtyBlockPtr NearestWaterByLatLon (
  FloatHi lat,
  FloatHi lon,
  FloatHi range,
  FloatHi PNTR distanceP
)

{
  CtrySetPtr  csp;

  csp = GetLatLonWaterData ();
  if (csp == NULL) return NULL;

  return RegionClosestToLatLon (lat, lon, range, distanceP, csp);
}

NLM_EXTERN CharPtr CountryClosestToLatLon (
  FloatHi lat,
  FloatHi lon,
  FloatHi range,
  FloatHi PNTR distanceP
)

{
  CtyBlockPtr  cbp;

  cbp = NearestCountryByLatLon (lat, lon, range, distanceP);
  if (cbp == NULL) return NULL;

  return cbp->name;
}

NLM_EXTERN CharPtr WaterClosestToLatLon (
  FloatHi lat,
  FloatHi lon,
  FloatHi range,
  FloatHi PNTR distanceP
)

{
  CtyBlockPtr  cbp;

  cbp = NearestWaterByLatLon (lat, lon, range, distanceP);
  if (cbp == NULL) return NULL;

  return cbp->name;
}

static CtyBlockPtr RegionIsNearLatLon (
  CharPtr country,
  CharPtr province,
  FloatHi lat,
  FloatHi lon,
  FloatHi range,
  FloatHi PNTR distanceP,
  CtrySetPtr csp
)

{
  LatBlockPtr PNTR  array;
  CtyBlockPtr       cbp, best = NULL;
  FloatHi           closest = EARTH_RADIUS * CONST_PI * 2;
  FloatHi           delta;
  Int4              latitude;
  Int4              longitude;
  Int4              maxDelta;
  LatBlockPtr       lbp;
  Int4              R;
  Int4              x;
  Int4              y;

  if (distanceP != NULL) {
    *distanceP = 0.0;
  }

  if (StringHasNoText (country)) return NULL;
  if (csp == NULL) return NULL;

  array = csp->latarray;
  if (array == NULL) return NULL;

  latitude = ConvertLat (lat, csp->scale);
  longitude = ConvertLon (lon, csp->scale);

  maxDelta = (Int4) (range * csp->scale + EPSILON);

  for (R = GetLatLonIndex (csp, array, latitude - maxDelta); R < csp->numLatBlocks; R++) {
    lbp = array [R];
    if (lbp == NULL) break;
    if (latitude + maxDelta < lbp->lat) break;

    if (longitude < lbp->minlon - maxDelta) continue;
    if (longitude > lbp->maxlon + maxDelta) continue;

    cbp = lbp->landmass;
    if (cbp == NULL) continue;

    if (StringICmp (country, cbp->level0) != 0) continue;
    if (/* province != NULL && */ StringICmp (province, cbp->level1) != 0) continue;

    if (longitude < lbp->minlon) {
      x = lbp->minlon;
    } else if (longitude > lbp->maxlon) {
      x = lbp->maxlon;
    } else {
      x = longitude;
    }

    y = lbp->lat;

    delta = DistanceOnGlobe (lat, lon, (FloatHi) (y / csp->scale), (FloatHi) (x / csp->scale));

    if (best == NULL || delta < closest) {
      best = cbp;
      closest = delta;
    }
  }

  if (best != NULL) {
    if (distanceP != NULL) {
      *distanceP = closest;
    }
  }

  return best;
}

static CtyBlockPtr CountryToLatLonDistance (
  CharPtr country,
  CharPtr province,
  FloatHi lat,
  FloatHi lon,
  FloatHi range,
  FloatHi PNTR distanceP
)

{
  CtrySetPtr  csp;

  csp = GetLatLonCountryData ();
  if (csp == NULL) return NULL;

  return RegionIsNearLatLon (country, province, lat, lon, range, distanceP, csp);
}

static CtyBlockPtr WaterToLatLonDistance (
  CharPtr country,
  FloatHi lat,
  FloatHi lon,
  FloatHi range,
  FloatHi PNTR distanceP
)

{
  CtrySetPtr  csp;

  csp = GetLatLonWaterData ();
  if (csp == NULL) return NULL;

  return RegionIsNearLatLon (country, NULL, lat, lon, range, distanceP, csp);
}

NLM_EXTERN Boolean CountryIsNearLatLon (
  CharPtr country,
  FloatHi lat,
  FloatHi lon,
  FloatHi range,
  FloatHi PNTR distanceP
)

{
  CtyBlockPtr  cbp;

  cbp = CountryToLatLonDistance (country, NULL, lat, lon, range, distanceP);
  if (cbp == NULL) return FALSE;

  return TRUE;
}

NLM_EXTERN Boolean WaterIsNearLatLon (
  CharPtr country,
  FloatHi lat,
  FloatHi lon,
  FloatHi range,
  FloatHi PNTR distanceP
)

{
  CtyBlockPtr  cbp;

  cbp = WaterToLatLonDistance (country, lat, lon, range, distanceP);
  if (cbp == NULL) return FALSE;

  return TRUE;
}

/*
static void WriteLatLonRegionData (
  CtrySetPtr csp,
  FILE* fp
)

{
  Char         buf [150];
  CtyBlockPtr  cbp;
  LatBlockPtr  lbp;
  ValNodePtr   vnp;

  if (csp == NULL || fp == NULL) return;

  for (vnp = csp->latblocks; vnp != NULL; vnp = vnp->next) {
    lbp = (LatBlockPtr) vnp->data.ptrvalue;
    if (lbp == NULL) {
      fprintf (fp, "NULL LatBlockPtr\n");
      continue;
    }
    cbp = lbp->landmass;
    if (cbp == NULL) {
      fprintf (fp, "NULL CtyBlockPtr\n");
      continue;
    }

    if (StringHasNoText (cbp->name)) {
      fprintf (fp, "NULL cbp->name\n");
      continue;
    }

    StringNCpy_0 (buf, cbp->name, 50);
    StringCat (buf, "                                                  ");
    buf [50] = '\0';

    fprintf (fp, "%s %4d : %4d  .. %4d\n", buf, (int) lbp->lat, (int) lbp->minlon, (int) lbp->maxlon);
  }

  fprintf (fp, "\n\n");
}

static void TestLatLonCountryData (void)

{
  CtrySetPtr  csp;
  FILE        *fp;

  fp = FileOpen ("stdout", "w");
  if (fp == NULL) {
    Message (MSG_OK, "Unable to open output file");
    return;
  }

  csp = GetLatLonCountryData ();
  if (csp == NULL) {
    fprintf (fp, "GetLatLonCountryData failed\n");
    FileClose (fp);
    return;
  }

  WriteLatLonRegionData (csp, fp);

  csp = GetLatLonWaterData ();
  if (csp == NULL) {
    fprintf (fp, "GetLatLonWaterData failed\n");
    FileClose (fp);
    return;
  }

  WriteLatLonRegionData (csp, fp);

  FileClose (fp);
}
*/

/* END OF NEW LATITUDE-LONGITUDE COUNTRY VALIDATION CODE */

static Boolean StringListIsUnique (ValNodePtr list)

{
  CharPtr     last;
  ValNodePtr  next;
  CharPtr     str;
  ValNodePtr  vnp;

  if (list == NULL) return TRUE;
  last = (CharPtr) list->data.ptrvalue;
  vnp = list->next;
  while (vnp != NULL) {
    next = vnp->next;
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringICmp (last, str) == 0) {
      return FALSE;
    } else {
      last = (CharPtr) vnp->data.ptrvalue;
    }
    vnp = next;
  }

  return TRUE;
}

static CharPtr modified_base_abbrevs [] = {
  "<ac4c>",
  "<chm5u>",
  "<cm>",
  "<cmnm5s2u>",
  "<cmnm5u>",
  "<d>",
  "<fm>",
  "<gal q>",
  "<gm>",
  "<i>",
  "<i6a>",
  "<m1a>",
  "<m1f>",
  "<m1g>",
  "<m1i>",
  "<m22g>",
  "<m2a>",
  "<m2g>",
  "<m3c>",
  "<m5c>",
  "<m6a>",
  "<m7g>",
  "<mam5u>",
  "<mam5s2u>",
  "<man q>",
  "<mcm5s2u>",
  "<mcm5u>",
  "<mo5u>",
  "<ms2i6a>",
  "<ms2t6a>",
  "<mt6a>",
  "<mv>",
  "<o5u>",
  "<osyw>",
  "<p>",
  "<q>",
  "<s2c>",
  "<s2t>",
  "<s2u>",
  "<s4u>",
  "<t>",
  "<t6a>",
  "<tm>",
  "<um>",
  "<yw>",
  "<x>",
  "<OTHER>",
  NULL
};

static void InitializeModBaseFSA (ValidStructPtr vsp)

{
  Int2  i;

  vsp->modifiedBases = TextFsaNew ();
  for (i = 0; modified_base_abbrevs [i] != NULL; i++) {
    TextFsaAdd (vsp->modifiedBases, modified_base_abbrevs [i]);
  }
}

static Boolean PrimerSeqIsValid (ValidStructPtr vsp, CharPtr name, Char PNTR badch)

{
  Char        ch;
  TextFsaPtr  fsa;
  size_t      len;
  ValNodePtr  matches;
  CharPtr     ptr;
  Int4        state;
  Boolean     first;

  if (badch != NULL) {
    *badch = '\0';
  }

  if (vsp == NULL) return FALSE;
  if (vsp->modifiedBases == NULL) {
    InitializeModBaseFSA (vsp);
  }
  fsa = vsp->modifiedBases;
  if (fsa == NULL) return FALSE;

  if (StringHasNoText (name)) return FALSE;
  len = StringLen (name);
  if (len < 1) return FALSE;

  if (StringChr (name, ',') != NULL) {
    if (name [0] != '(' || name [len - 1] != ')') return FALSE;
  } else {
    if (StringChr (name, '(') != NULL) return FALSE;
    if (StringChr (name, ')') != NULL) return FALSE;
  }

  if (StringChr (name, ';') != NULL) return FALSE;
  /* if (StringChr (name, ' ') != NULL) return FALSE; */

  ptr = name;
  ch = *ptr;
  while (ch != '\0') {
    if (ch == '<') {
      state = 0;
      matches = NULL;
      first = TRUE;
      while (ch != '\0' && ch != '>' && (first || ch != '<')) {
        state = TextFsaNext (fsa, state, ch, &matches);
        ptr++;
        ch = *ptr;
        first = FALSE;
      }
      if (ch != '>' || ch == '<') {
        if (badch != NULL) {
          *badch = ch;
        }
        return FALSE;
      }
      state = TextFsaNext (fsa, state, ch, &matches);
      if (matches == NULL) {
        if (badch != NULL) {
          *badch = ch;
        }
        return FALSE;
      }
    } else {
      if (ch != '(' && ch != ')' && ch != ',' && ch != ':') {
        if (! (IS_ALPHA (ch))) {
          if (badch != NULL) {
            *badch = ch;
          }
          return FALSE;
        }
        ch = TO_UPPER (ch);
        if (StringChr ("ABCDGHKMNRSTVWY", ch) == NULL) {
          if (badch != NULL) {
            ch = TO_LOWER (ch);
            *badch = ch;
          }
          return FALSE;
        }
      }
    }
    ptr++;
    ch = *ptr;
  }

  return TRUE;
}

/*
static ValNodePtr ParsePrimerSeqIntoComponents (
  CharPtr strs
)

{
  Char        ch;
  ValNodePtr  head = NULL;
  CharPtr     ptr, str, tmp;

  if (StringHasNoText (strs)) return NULL;

  tmp = StringSave (strs);
  if (tmp == NULL) return NULL;

  str = tmp;
  while (StringDoesHaveText (str)) {
    ptr = str;
    ch = *ptr;

    while (ch != '\0' && ch != '(' && ch != ')' && ch != ',' && ch != ';' && ch != ':') {
      ptr++;
      ch = *ptr;
    }
    if (ch != '\0' && ptr != NULL) {
      *ptr = '\0';
      ptr++;
    }

    TrimSpacesAroundString (str);
    if (StringDoesHaveText (str)) {
      ValNodeCopyStr (&head, 0, str);
    }

    str = ptr;
  }

  MemFree (tmp);
  return head;
}

static Boolean PrimerSeqHasDuplicates (CharPtr name)

{
  ValNodePtr  head;
  Boolean     rsult = FALSE;

  if (StringHasNoText (name)) return FALSE;

  head = ParsePrimerSeqIntoComponents (name);
  if (head == NULL) return FALSE;
  head = ValNodeSort (head, SortVnpByString);
  if (! StringListIsUnique (head)) {
    rsult = TRUE;
  }
  ValNodeFreeData (head);

  return rsult;
}
*/

static Int2 CountDigits (CharPtr str)

{
  Char  ch;
  Int2  count = 0;

  if (str == NULL) return count;
  ch = *str;
  while (IS_DIGIT (ch)) {
    count++;
    str++;
    ch = *str;
  }
  return count;
}

static Boolean LatLonIsValid (CharPtr name)

{
  Char     ch;
  Int2     count;
  CharPtr  str;

  if (StringHasNoText (name)) return FALSE;
  str = name;

  count = CountDigits (str);
  if (count < 1 || count > 2) return FALSE;
  str += count;

  ch = *str;
  if (ch == '.') {
    str++;
    count = CountDigits (str);
    if (count != 2) return FALSE;
    str += count;
  }

  ch = *str;
  if (ch != ' ') return FALSE;
  str++;
  ch = *str;
  if (ch != 'N' && ch != 'S') return FALSE;
  str++;
  ch = *str;
  if (ch != ' ') return FALSE;
  str++;

  count = CountDigits (str);
  if (count < 1 || count > 3) return FALSE;
  str += count;

  ch = *str;
  if (ch == '.') {
    str++;
    count = CountDigits (str);
    if (count != 2) return FALSE;
    str += count;
  }

  ch = *str;
  if (ch != ' ') return FALSE;
  str++;
  ch = *str;
  if (ch != 'E' && ch != 'W') return FALSE;
  str++;

  ch = *str;
  if (ch != '\0') return FALSE;

  return TRUE;
}

static CharPtr source_qual_prefixes [] = {
  "acronym:",
  "anamorph:",
  "authority:",
  "biotype:",
  "biovar:",
  "bio_material:",
  "breed:",
  "cell_line:",
  "cell_type:",
  "chemovar:",
  "chromosome:",
  "clone:",
  "clone_lib:",
  "collected_by:",
  "collection_date:",
  "common:",
  "country:",
  "cultivar:",
  "culture_collection:",
  "dev_stage:",
  "dosage:",
  "ecotype:",
  "endogenous_virus_name:",
  "environmental_sample:",
  "forma:",
  "forma_specialis:",
  "frequency:",
  "fwd_pcr_primer_name",
  "fwd_pcr_primer_seq",
  "fwd_primer_name",
  "fwd_primer_seq",
  "genotype:",
  "germline:",
  "group:",
  "haplogroup:",
  "haplotype:",
  "identified_by:",
  "insertion_seq_name:",
  "isolate:",
  "isolation_source:",
  "lab_host:",
  "lat_lon:",
  "left_primer:",
  "linkage_group:",
  "map:",
  "mating_type:",
  "metagenome_source:",
  "metagenomic:",
  "nat_host:",
  "pathovar:",
  "placement:",
  "plasmid_name:",
  "plastid_name:",
  "pop_variant:",
  "rearranged:",
  "rev_pcr_primer_name",
  "rev_pcr_primer_seq",
  "rev_primer_name",
  "rev_primer_seq",
  "right_primer:",
  "segment:",
  "serogroup:",
  "serotype:",
  "serovar:",
  "sex:",
  "specimen_voucher:",
  "strain:",
  "subclone:",
  "subgroup:",
  "substrain:",
  "subtype:",
  "sub_species:",
  "synonym:",
  "taxon:",
  "teleomorph:",
  "tissue_lib:",
  "tissue_type:",
  "transgenic:",
  "transposon_name:",
  "type:",
  "variety:",
  NULL
};

static void InitializeSourceQualTags (ValidStructPtr vsp)

{
  Int2  i;

  vsp->sourceQualTags = TextFsaNew ();
  for (i = 0; source_qual_prefixes [i] != NULL; i++) {
    TextFsaAdd (vsp->sourceQualTags, source_qual_prefixes [i]);
  }
}

static void ValidateSourceQualTags (ValidStructPtr vsp, GatherContextPtr gcp, BioSourcePtr biop, CharPtr str)

{
  Char        ch;
  CharPtr     hit;
  Boolean     okay;
  CharPtr     ptr;
  CharPtr     tmp;
  Int4        state;
  ValNodePtr  matches;

  if (vsp->sourceQualTags == NULL || StringHasNoText (str)) return;
  state = 0;
  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    matches = NULL;
    state = TextFsaNext (vsp->sourceQualTags, state, ch, &matches);
    if (matches != NULL) {
      hit = (CharPtr) matches->data.ptrvalue;
      if (StringHasNoText (hit)) {
        hit = "?";
      }
      okay = TRUE;
      tmp = ptr - StringLen (hit);
      if (tmp > str) {
        ch = *tmp;
        if ((! IS_WHITESP (ch)) && ch != ';') {
          okay = FALSE;
        }
      }
      if (okay) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_StructuredSourceNote,
                  "Source note has structured tag '%s'", hit);
      }
    }
    ptr++;
    ch = *ptr;
  }
}


static CharPtr GetOrgModWarning (Uint2 subtype)
{
  CharPtr warning = NULL;

  switch (subtype) {
    /*
    case ORGMOD_biovar:
      warning = "Biovar value specified is not found in taxname";
      break;
    */
    case ORGMOD_forma:
      warning = "Forma value specified is not found in taxname";
      break;
    case ORGMOD_forma_specialis:
      warning = "Forma specialis value specified is not found in taxname";
      break;
    /*
    case ORGMOD_pathovar:
      warning = "Pathovar value specified is not found in taxname";
      break;
    */
    case ORGMOD_sub_species:
      warning = "Subspecies value specified is not found in taxname";
      break;
    case ORGMOD_variety:
      warning = "Variety value specified is not found in taxname";
      break;
  }
  return warning;
}


static Boolean ValidateOrgModInTaxName (ValidStructPtr vsp, OrgModPtr mod, CharPtr taxname, Boolean varietyOK)
{
  CharPtr cp, f, warn;
  Int4    word_len, name_len;

  if (vsp == NULL || mod == NULL) return FALSE;

  name_len = StringLen (mod->subname);

  /* skip first word */
  word_len = StringCSpn (taxname, " ");
  cp = taxname + word_len;
  cp += StringSpn (cp, " ");
  /* skip second word */
  word_len = StringCSpn (cp, " ");
  cp += word_len;
  cp += StringSpn (cp, " ");

  f = StringSearch (cp, mod->subname);
  while (f != NULL && ((f != cp && isalpha (*(f - 1))) || isalpha (*(f + name_len)))) {
    f = StringSearch (f + 1, mod->subname);
  }
  if (f == NULL) {
    warn = GetOrgModWarning (mod->subtype);
    if (warn != NULL) {
      /* variety is sorted before sub_species, so if variety was okay in taxname, can ignore missing sub_species */
      if (mod->subtype == ORGMOD_sub_species && varietyOK) return FALSE;
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, warn);
      return FALSE;
    }
  }

  return TRUE;
}

/* institution:collection is now stored as a ValNode list of strings, sorted and indexed */

static Boolean inst_code_not_found = FALSE;

static ValNodePtr    ic_code_list = NULL;
static CharPtr PNTR  ic_code_data = NULL;
static Uint1 PNTR    ic_code_type = NULL;
static Int4          ic_code_len = 0;

#define BIO_MATERIAL_TYPE       1
#define CULTURE_COLLECTION_TYPE 2
#define SPECIMEN_VOUCHER_TYPE   4

static void SetupInstCollTable (void)

{
  FileCache   fc;
  CharPtr     file = "institution_codes.txt";
  FILE        *fp = NULL;
  Int4        i;
  ValNodePtr  last = NULL;
  Char        line [512];
  Char        path [PATH_MAX];
  CharPtr     ptr;
  ErrSev      sev;
  CharPtr     str;
  CharPtr     tmp;
  Uint1       type;
  ValNodePtr  vnp;

  if (ic_code_data != NULL) return;
  if (inst_code_not_found) return;

  if (FindPath ("ncbi", "ncbi", "data", path, sizeof (path))) {
    FileBuildPath (path, NULL, file);
    sev = ErrSetMessageLevel (SEV_ERROR);
    fp = FileOpen (path, "r");
    ErrSetMessageLevel (sev);
  }

  if (fp == NULL) {
    inst_code_not_found = TRUE;
    return;
  }

  FileCacheSetup (&fc, fp);

  str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
  while (str != NULL) {
    if (StringDoesHaveText (str)) {
      type = 0;
      ptr = StringChr (str, '\t');
      if (ptr != NULL) {
        *ptr = '\0';
        ptr++;
        tmp = StringChr (ptr, '\t');
        if (tmp != NULL) {
          *tmp = '\0';
          if (StringChr (ptr, 'b') != NULL) {
            type |= BIO_MATERIAL_TYPE;
          }
          if (StringChr (ptr, 'c') != NULL) {
            type |= CULTURE_COLLECTION_TYPE;
          }
          if (StringChr (ptr, 's') != NULL) {
            type |= SPECIMEN_VOUCHER_TYPE;
          }
        }
      }
      TrimSpacesAroundString (str);
      vnp = ValNodeCopyStr (&last, type, str);
      if (ic_code_list == NULL) {
        ic_code_list = vnp;
      }
      last = vnp;
    }
    str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
  }

  FileClose (fp);

  ic_code_len = ValNodeLen (ic_code_list);
  if (ic_code_len > 0) {
    ic_code_list = ValNodeSort (ic_code_list, SortVnpByString);
    ic_code_data = (CharPtr PNTR) MemNew (sizeof (CharPtr) * (ic_code_len + 1));
    if (ic_code_data != NULL) {
      for (vnp = ic_code_list, i = 0; vnp != NULL; vnp = vnp->next, i++) {
        str = (CharPtr) vnp->data.ptrvalue;
        ic_code_data [i] = str;
      }
    }

    ic_code_type = (Uint1 PNTR) MemNew (sizeof (Uint1) * (ic_code_len + 1));
    if (ic_code_type != NULL) {
      for (vnp = ic_code_list, i = 0; vnp != NULL; vnp = vnp->next, i++) {
        ic_code_type [i] = vnp->choice;
      }
    }
  }
}

static CharPtr CheckInstCollName (CharPtr name, Uint1Ptr typeP)

{
  Int4     L, R, mid;
  CharPtr  str;

  SetupInstCollTable ();

  if (typeP != NULL) {
    *typeP = 0;
  }

  L = 0;
  R = ic_code_len - 1;
  while (L < R) {
    mid = (L + R) / 2;
    str = ic_code_data [(int) mid];
    if (StringICmp (str, name) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }
  if (R < 0) return NULL;

  if (typeP != NULL) {
    *typeP = ic_code_type [(int) R];
  }

  return ic_code_data [(int) R];
}

static void ValidateOrgModVoucher (ValidStructPtr vsp, OrgModPtr mod)

{
  Char     buf [512];
  CharPtr  inst = NULL, id = NULL, coll = NULL, str;
  size_t   len1, len2;
  Uint1    type;

  if (vsp == NULL || mod == NULL) return;

  StringNCpy_0 (buf, mod->subname, sizeof (buf));
  if (StringChr (buf, ':') == NULL) {
    if (mod->subtype == ORGMOD_culture_collection) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_UnstructuredVoucher, "Culture_collection should be structured, but is not");
    }
    return;
  }
  if (! ParseStructuredVoucher (buf, &inst, &id) || inst == NULL || inst[0] == ':') {
    if (StringHasNoText (inst) || inst [0] == ':') {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadInstitutionCode, "Voucher is missing institution code");
    }
    if (StringHasNoText (id)) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadVoucherID, "Voucher is missing specific identifier");
    }
    return;
  }
  if (inst == NULL) return;

  str = CheckInstCollName (inst, &type);
  if (StringCmp (str, inst) == 0) {
    if ((mod->subtype == ORGMOD_bio_material && (type & BIO_MATERIAL_TYPE) == 0) ||
        (mod->subtype == ORGMOD_culture_collection && (type & CULTURE_COLLECTION_TYPE) == 0) ||
        (mod->subtype == ORGMOD_specimen_voucher && (type & SPECIMEN_VOUCHER_TYPE) == 0)) {
      if ((type & BIO_MATERIAL_TYPE) != 0) {
        ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_WrongVoucherType, "Institution code %s should be bio_material", inst);
      } else if ((type & CULTURE_COLLECTION_TYPE) != 0) {
        ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_WrongVoucherType, "Institution code %s should be culture_collection", inst);
      } else if ((type & SPECIMEN_VOUCHER_TYPE) != 0) {
        ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_WrongVoucherType, "Institution code %s should be specimen_voucher", inst);
      }
    }
    return;
  }

  if (StringICmp (str, inst) == 0) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadInstitutionCode,
              "Institution code %s exists, but correct capitalization is %s", inst, str);
    return;
  }

  /* previously ignored personal collections, now complain if name missing */
  if (StringNICmp (inst, "personal", 8) == 0) {
    if (StringICmp (inst, "personal") == 0 && StringLen (str) > 0) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MissingPersonalCollectionName,
                "Personal collection does not have name of collector");
    }
    return;
  }

  len1 = StringLen (inst);
  len2 = StringLen (str);

  if (len1 < len2) {
    if (StringNICmp (str, inst, len1) == 0 && str [len1] == '<') {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadInstitutionCode, "Institution code %s needs to be qualified with a <COUNTRY> designation", inst);
      return;
    }
  }

  coll = StringChr (inst, ':');
  if (coll == NULL) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadInstitutionCode, "Institution code %s is not in list", inst);
    return;
  }

  *coll = '\0';
  coll++;
  str = CheckInstCollName (inst, &type);
  if (StringCmp (str, inst) == 0) {
    if (StringCmp (coll, "DNA") == 0) {
      /* DNA is a valid collection for any institution (using bio_material) */
      if (mod->subtype != ORGMOD_bio_material) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_WrongVoucherType, "DNA should be bio_material");
      }
    } else {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadCollectionCode,
                "Institution code %s exists, but collection %s:%s is not in list", inst, inst, coll);
    }
    return;
  }

  len1 = StringLen (inst);
  len2 = StringLen (str);

  if (len1 < len2) {
    if (StringNICmp (str, inst, len1) == 0 && str [len1] == '<') {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadInstitutionCode, "Institution code in %s:%s needs to be qualified with a <COUNTRY> designation", inst, coll);
      return;
    }
  }

  ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadInstitutionCode, "Institution code %s:%s is not in list", inst, coll);
}

NLM_EXTERN Boolean VoucherInstitutionIsValid (CharPtr inst)

{
  CharPtr  str;
  Uint1    type;

  if (StringHasNoText (inst)) return FALSE;

  str = CheckInstCollName (inst, &type);
  if (StringCmp (str, inst) == 0) return TRUE;

  return FALSE;
}

/* works on subname copy that it can change */

NLM_EXTERN Boolean ParseStructuredVoucher (
  CharPtr subname,
  CharPtr PNTR inst,
  CharPtr PNTR id
)

{
  CharPtr  ptr;
  CharPtr  tmp;

  if (StringHasNoText (subname)) return FALSE;
  if (StringLen (subname) < 3) return FALSE;
  TrimSpacesAroundString (subname);

  ptr = StringChr (subname, ':');
  if (ptr == NULL) return FALSE;

  *inst = subname;

  tmp = StringChr (ptr + 1, ':');
  if (tmp != NULL) {
    *tmp = '\0';
    tmp++;
    TrimSpacesAroundString (tmp);
    *id = tmp;
  } else {
    *ptr = '\0';
    ptr++;
    TrimSpacesAroundString (ptr);
    *id = ptr;
  }

  if (StringHasNoText (*inst) || StringHasNoText (*id)) return FALSE;

  return TRUE;
}

static void ValidateLatLon (ValidStructPtr vsp, CharPtr lat_lon)

{
  Boolean format_ok = FALSE, lat_in_range = FALSE, lon_in_range = FALSE, precision_ok = FALSE;
  CharPtr ptr;
  Char    tmp [128];

  IsCorrectLatLonFormat (lat_lon, &format_ok, &precision_ok, &lat_in_range, &lon_in_range);

  if (! format_ok) {
    /* may have comma and then altitude, so just get lat_lon component */
    StringNCpy_0 (tmp, lat_lon, sizeof (tmp));
    ptr = StringChr (tmp, ',');
    if (ptr != NULL) {
      *ptr = '\0';
      lat_lon = tmp;
      IsCorrectLatLonFormat (tmp, &format_ok, &precision_ok, &lat_in_range, &lon_in_range);
      if (format_ok) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_LatLonFormat, "lat_lon format has extra text after correct dd.dd N|S ddd.dd E|W format");
      }
    }
  }

  if (!format_ok) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_LatLonFormat, "lat_lon format is incorrect - should be dd.dd N|S ddd.dd E|W");
  } else {
    if (!lat_in_range) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_LatLonRange, "latitude value is out of range - should be between 90.00 N and 90.00 S");
    }
    if (!lon_in_range) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_LatLonRange, "longitude value is out of range - should be between 180.00 E and 180.00 W");
    }
    if (! precision_ok) {
      /*
      ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_LatLonPrecision, "lat_lon precision is incorrect - should only have two digits to the right of the decimal point");
      */
    }
  }
}


static void ValidateLocationForHIV (ValidStructPtr vsp, BioSourcePtr biop, BioseqPtr bsp)
{
  SeqDescrPtr sdp;
  SeqMgrDescContext context;
  MolInfoPtr mip;

  if (vsp == NULL || biop == NULL) {
    return;
  }

  if (bsp != NULL) {
    if (bsp->mol == Seq_mol_dna) {
      if (biop->genome != GENOME_proviral) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "HIV with moltype DNA should be proviral");
      }
    } else if (bsp->mol == Seq_mol_rna) {
      sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &context);
      if (sdp != NULL
          && (mip = (MolInfoPtr) sdp->data.ptrvalue) != NULL
          && mip->biomol == MOLECULE_TYPE_MRNA) {
        ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_BioSourceInconsistency, "HIV with mRNA molecule type is rare");
      }
    }
  }
}

/*****************************************************************************
*
*   ValidateSeqDescrContext(gcp)
*      Gather callback helper function for validating context on a Bioseq
*
*****************************************************************************/
static Boolean UnbalancedParentheses (CharPtr str)

{
  Char  ch;
  Int4  fwd_par = 0, rev_par = 0, fwd_bkt = 0, rev_bkt = 0;

  if (StringHasNoText (str)) return FALSE;

  ch = *str;
  while (ch != '\0') {
    if (ch == '(') {
      fwd_par++;
    } else if (ch == ')') {
      rev_par++;
    } else if (ch == '[') {
      fwd_bkt++;
    } else if (ch == ']') {
      rev_bkt++;
    }
    if (fwd_par < rev_par) return TRUE;
    if (fwd_bkt < rev_bkt) return TRUE;
    str++;
    ch = *str;
  }

  if (fwd_par != rev_par) return TRUE;
  if (fwd_bkt != rev_bkt) return TRUE;

  return FALSE;
}

static CharPtr sgml_strings [] = {
  "&gt;",
  "&lt;",
  "&amp;",
  "&agr;",
  "&Agr;",
  "&bgr;",
  "&Bgr;",
  "&ggr;",
  "&Ggr;",
  "&dgr;",
  "&Dgr;",
  "&egr;",
  "&Egr;",
  "&zgr;",
  "&Zgr;",
  "&eegr;",
  "&EEgr;",
  "&thgr;",
  "&THgr;",
  "&igr;",
  "&Igr;",
  "&kgr;",
  "&Kgr;",
  "&lgr;",
  "&Lgr;",
  "&mgr;",
  "&Mgr;",
  "&ngr;",
  "&Ngr;",
  "&xgr;",
  "&Xgr;",
  "&ogr;",
  "&Ogr;",
  "&pgr;",
  "&Pgr;",
  "&rgr;",
  "&Rgr;",
  "&sgr;",
  "&Sgr;",
  "&sfgr;",
  "&tgr;",
  "&Tgr;",
  "&ugr;",
  "&Ugr;",
  "&phgr;",
  "&PHgr;",
  "&khgr;",
  "&KHgr;",
  "&psgr;",
  "&PSgr;",
  "&ohgr;",
  "&OHgr;",
  NULL
};

static void InitializeSgmlStringsFSA (ValidStructPtr vsp)

{
  Int2  i;

  vsp->sgmlStrings = TextFsaNew ();
  for (i = 0; sgml_strings [i] != NULL; i++) {
    TextFsaAdd (vsp->sgmlStrings, sgml_strings [i]);
  }
}

static Boolean StringHasSgml (ValidStructPtr vsp, CharPtr str)

{
  Int2        ascii_len;
  Char        buf [256];
  Char        ch;
  TextFsaPtr  fsa;
  ValNodePtr  matches;
  Boolean     not_sgml;
  CharPtr     ptr;
  ErrSev      sev;
  Int4        state;

  if (StringHasNoText (str)) return FALSE;
  if (StringChr (str, '&') == NULL) return FALSE;

  if (vsp == NULL) return FALSE;
  if (vsp->sgmlStrings == NULL) {
    InitializeSgmlStringsFSA (vsp);
  }
  fsa = vsp->sgmlStrings;
  if (fsa == NULL) return FALSE;

  not_sgml = TRUE;
  state = 0;
  matches = NULL;
  for (ptr = str, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
    state = TextFsaNext (fsa, state, ch, &matches);
    if (matches != NULL) {
      not_sgml = FALSE;
    }
  }
  if (not_sgml) return FALSE;

  sev = ErrSetMessageLevel (SEV_REJECT);
  ascii_len = Sgml2AsciiLen (str);
  if (ascii_len + 2 >= sizeof (buf)) {
    ErrSetMessageLevel (sev);
    return FALSE;
  }

  buf [0] = '\0';
  Sgml2Ascii (str, buf, ascii_len + 1);
  ErrSetMessageLevel (sev);

  if (StringHasNoText (buf)) return FALSE;
  if (StringCmp (str, buf) == 0) return FALSE;

  return TRUE;
}

static CharPtr valid_sex_values [] = {
  "female",
  "male",
  "hermaphrodite",
  "unisexual",
  "bisexual",
  "asexual",
  "monoecious",
  "monecious",
  "dioecious",
  "diecious",
  NULL
};

static Boolean IsValidSexValue (CharPtr str)

{
  int  i;

  if (StringHasNoText (str)) return FALSE;

  for (i = 0; valid_sex_values [i] != NULL; i++) {
    if (StringICmp (str, valid_sex_values [i]) == 0) return TRUE;
  }

  return FALSE;
}

static Boolean LatLonInRange (
  FloatHi lat,
  FloatHi lon
)

{
  if (lat < -90.0001 || lat > 90.0001) return FALSE;
  if (lon < -180.0001 || lon > 180.0001) return FALSE;

  return TRUE;
}

static Boolean RegionIsClosestToLatLon (
  CharPtr country,
  FloatHi lat,
  FloatHi lon,
  FloatHi range,
  FloatHi PNTR distanceP,
  CtrySetPtr csp
)

{
  LatBlockPtr PNTR  array;
  CtyBlockPtr       cbp;
  FloatHi           closest = EARTH_RADIUS * CONST_PI * 2;
  CharPtr           guess = NULL;
  FloatHi           delta;
  Int4              latitude;
  Int4              longitude;
  Int4              maxDelta;
  LatBlockPtr       lbp;
  Int4              R;
  Int4              x;
  Int4              y;


  if (StringHasNoText (country)) return FALSE;

  if (distanceP != NULL) {
    *distanceP = 0.0;
  }

  if (csp == NULL) return FALSE;

  array = csp->latarray;
  if (array == NULL) return FALSE;

  latitude = ConvertLat (lat, csp->scale);
  longitude = ConvertLon (lon, csp->scale);

  maxDelta = (Int4) (range * csp->scale + EPSILON);

  for (R = GetLatLonIndex (csp, array, latitude - maxDelta); R < csp->numLatBlocks; R++) {
    lbp = array [R];
    if (lbp == NULL) break;
    if (latitude + maxDelta < lbp->lat) break;

    if (longitude < lbp->minlon - maxDelta) continue;
    if (longitude > lbp->maxlon + maxDelta) continue;

    cbp = lbp->landmass;
    if (cbp == NULL) continue;

    if (longitude < lbp->minlon) {
      x = lbp->minlon;
    } else if (longitude > lbp->maxlon) {
      x = lbp->maxlon;
    } else {
      x = longitude;
    }

    y = lbp->lat;

    delta = DistanceOnGlobe (lat, lon, (FloatHi) (y / csp->scale), (FloatHi) (x / csp->scale));

    if (delta < closest) {
      guess = cbp->name;
      closest = delta;
    } else if (delta == closest) {
      if (StringCmp (country, cbp->name) == 0) {
        guess = cbp->name;
      }
    }
  }

  if (guess != NULL) {
    if (distanceP != NULL) {
      *distanceP = closest;
    }
  }

  if (StringCmp (guess, country) == 0) return TRUE;

  return FALSE;
}


static Boolean CountryIsClosestToLatLon (
  CharPtr country,
  FloatHi lat,
  FloatHi lon,
  FloatHi range,
  FloatHi PNTR distanceP
)

{
  CtrySetPtr  csp;

  csp = GetLatLonCountryData ();
  if (csp == NULL) return FALSE;

  return RegionIsClosestToLatLon (country, lat, lon, range, distanceP, csp);
}


static int AdjustAndRoundDistance (
  FloatHi distance,
  FloatHi scale
)

{
  if (scale < 1.1) {
    distance += 111.19;
  } else if (scale > 19.5 && scale < 20.5) {
    distance += 5.56;
  } else if (scale > 99.5 && scale < 100.5) {
    distance += 1.11;
  }

  return (int) (distance + 0.5);
}

typedef struct latlonmap {
  FloatHi  lat;
  FloatHi  lon;
  CharPtr  fullguess;
  CharPtr  guesscountry;
  CharPtr  guessprovince;
  CharPtr  guesswater;
  CharPtr  closestfull;
  CharPtr  closestcountry;
  CharPtr  closestprovince;
  CharPtr  closestwater;
  CharPtr  claimedfull;
  int      landdistance;
  int      waterdistance;
  int      claimeddistance;
} LatLonMap, PNTR LatLonMapPtr;

static void CalculateLatLonMap (
  FloatHi lat,
  FloatHi lon,
  CharPtr country,
  CharPtr province,
  FloatHi scale,
  LatLonMapPtr lmp
)

{
  CtyBlockPtr  cbp;
  FloatHi      landdistance = 0.0, waterdistance = 0.0, claimeddistance = 0.0;
  Boolean      goodmatch = FALSE;

  if (lmp == NULL) return;

  /* initialize result values */
  MemSet ((Pointer) lmp, 0, sizeof (LatLonMap));

  lmp->lat = lat;
  lmp->lon = lon;

  /* lookup region by coordinates, or find nearest region and calculate distance */
  cbp = GuessCountryByLatLon (lat, lon, country, province);
  if (cbp != NULL) {
    /* successfully found inside some country */
    lmp->fullguess = cbp->name;
    lmp->guesscountry = cbp->level0;
    lmp->guessprovince = cbp->level1;
    if (StringICmp (country, lmp->guesscountry) == 0 && (province == NULL || StringICmp (province, lmp->guessprovince) == 0)) {
      goodmatch = TRUE;
    }
  } else {
    /* not inside a country, check water */
    cbp = GuessWaterByLatLon (lat, lon, country);
    if (cbp != NULL) {
      /* found inside water */
      lmp->guesswater = cbp->name;
      if (StringICmp (country, lmp->guesswater) == 0) {
        goodmatch = TRUE;
      }
      /*
      also see if close to land for coastal warning (if country is land)
      or proximity message (if country is water)
      */
      cbp = NearestCountryByLatLon (lat, lon, 5.0, &landdistance);
      if (cbp != NULL) {
        lmp->closestfull = cbp->name;
        lmp->closestcountry = cbp->level0;
        lmp->closestprovince = cbp->level1;
        lmp->landdistance = AdjustAndRoundDistance (landdistance, scale);
        if (StringICmp (country, lmp->closestcountry) == 0 && (province == NULL || StringICmp (province, lmp->closestprovince) == 0)) {
          goodmatch = TRUE;
        }
      }
    } else {
      /* may be coastal inlet, area of data insufficiency */
      cbp = NearestCountryByLatLon (lat, lon, 5.0, &landdistance);
      if (cbp != NULL) {
        lmp->closestfull = cbp->name;
        lmp->closestcountry = cbp->level0;
        lmp->closestprovince = cbp->level1;
        lmp->landdistance = AdjustAndRoundDistance (landdistance, scale);
        if (StringICmp (country, lmp->closestcountry) == 0 && (province == NULL || StringICmp (province, lmp->closestprovince) == 0)) {
          goodmatch = TRUE;
        }
      }
      cbp = NearestWaterByLatLon (lat, lon, 5.0, &waterdistance);
      if (cbp != NULL) {
        lmp->closestwater = cbp->level0;
        lmp->waterdistance = AdjustAndRoundDistance (waterdistance, scale);
        if (StringICmp (country, lmp->closestwater) == 0) {
          goodmatch = TRUE;
        }
      }
    }
  }
  /* if guess is not the provided country or province, calculate distance to claimed country */
  if (! goodmatch) {
    cbp = CountryToLatLonDistance (country, province, lat, lon, 5.0, &claimeddistance);
    if (cbp != NULL) {
      if (claimeddistance < ErrorDistance(lmp->lat, lmp->lon, scale)) {
        lmp->guesscountry = country;
        lmp->guessprovince = province;
        lmp->fullguess = cbp->name;
      } else {
        lmp->claimedfull = cbp->name;
        lmp->claimeddistance = AdjustAndRoundDistance (claimeddistance, scale);
      }
    } else if (province == NULL) {
      cbp = WaterToLatLonDistance (country, lat, lon, 5.0, &claimeddistance);
      if (cbp != NULL) {
        lmp->claimedfull = cbp->name;
        lmp->claimeddistance = AdjustAndRoundDistance (claimeddistance, scale);
      }
    }
  }
}


enum {
  eLatLonClassify_CountryMatch = 1 ,
  eLatLonClassify_ProvinceMatch = 2 ,
  eLatLonClassify_WaterMatch = 4 ,
  eLatLonClassify_CountryClosest = 8 ,
  eLatLonClassify_ProvinceClosest = 16 ,
  eLatLonClassify_WaterClosest = 32 ,
  eLatLonClassify_Error = 256
} ELatLonClassify;


static Uint4 ClassifyLatLonMap (
  CharPtr fullname,
  CharPtr country,
  CharPtr province,
  LatLonMapPtr lmp
)

{
  Uint4 rval = 0;

  if (lmp == NULL) return eLatLonClassify_Error;

  /* compare guesses or closest regions to indicated country and province */
  if (lmp->guesscountry != NULL) {

    /* if top level countries match */
    if (StringICmp (country, lmp->guesscountry) == 0) {
      rval |= eLatLonClassify_CountryMatch;
      /* if both are null, call it a match */
      if (StringICmp (province, lmp->guessprovince) == 0) {
        rval |= eLatLonClassify_ProvinceMatch;
      }
    }
    /* if they don't match, do they overlap or are closest? */
    if (!(rval & eLatLonClassify_CountryMatch)) {
      if (StringICmp (country, lmp->closestcountry) == 0) {
        rval |= eLatLonClassify_CountryClosest;
        if (StringICmp (province, lmp->closestprovince) == 0) {
          rval |= eLatLonClassify_ProvinceClosest;
        }
      }
    } else if (!(rval & eLatLonClassify_ProvinceMatch) && province != NULL) {
      if (StringICmp (province, lmp->closestprovince) == 0) {
        rval |= eLatLonClassify_ProvinceClosest;
      }
    }
  }
  if (lmp->guesswater != NULL) {
    /* was the non-approved body of water correctly indicated? */
    if (StringICmp (country, lmp->guesswater) == 0) {
      rval |= eLatLonClassify_WaterMatch;
    } else if (StringICmp (country, lmp->closestwater) == 0) {
      rval |= eLatLonClassify_WaterClosest;
    } 
  }
  if (lmp->closestcountry != NULL && StringICmp (country, lmp->closestcountry) == 0) {
    if (lmp->guesscountry == NULL && lmp->guesswater == NULL) {
      /* coastal area */
      rval |= eLatLonClassify_CountryMatch;
      lmp->guesscountry = lmp->closestcountry;
      lmp->fullguess = lmp->closestcountry;
      if (lmp->closestprovince != NULL && StringICmp (province, lmp->closestprovince) == 0) {
        rval |= eLatLonClassify_ProvinceMatch;
        lmp->guessprovince = lmp->closestprovince;
        lmp->fullguess = lmp->closestfull;
      }
    } else {      
      rval |= eLatLonClassify_CountryClosest;
      if (lmp->closestprovince != NULL && StringICmp (province, lmp->closestprovince) == 0) {
        rval |= eLatLonClassify_ProvinceClosest;
      }
    }
  }
  return rval;
}


static void LatLonWaterErrors (
  ValidStructPtr vsp,
  LatLonMapPtr lmp,
  Uint4 test,
  FloatHi neardist,
  CharPtr country,
  CharPtr province,
  CharPtr lat_lon,
  CharPtr fullname,
  FloatHi scale
  )
{
  CharPtr fmt = "Lat_lon '%s' is closest to %s'%s' at distance %d km, but in water '%s'";
  CharPtr claimed_fmt = "Lat_lon '%s' is closest to %s'%s' at distance %d km, but in water '%s' - claimed region '%s' is at distance %d km";

  Boolean suppress = FALSE;
  CharPtr reportregion;
  CharPtr nosubphrase = "";
  CharPtr desphrase = "designated subregion ";
  CharPtr subphrase = "another subregion ";
  CharPtr phrase = nosubphrase;
  Boolean show_claimed = FALSE;

  if (test & (eLatLonClassify_CountryClosest | eLatLonClassify_ProvinceClosest)) {

    if (lmp->landdistance < 22) {
      /* for now, will not report */
      /* this is a policy decision */
      suppress = TRUE;
    } else if (StringStr (fullname, "Island") != NULL) {
      suppress = TRUE;
    }

    if (test & eLatLonClassify_ProvinceClosest) {
      reportregion = fullname;
      phrase = desphrase;
    } else {
      /* wasn't closest province, so must be closest country */
      if (province != NULL && vsp->testLatLonSubregion) {
        phrase = subphrase;
        reportregion = lmp->closestfull;
      } else {
        reportregion = lmp->closestcountry;
      }
      if (lmp->claimedfull != NULL) {
        show_claimed = TRUE;
      }
    }

    if (!suppress) {
      if (show_claimed) {
        ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_LatLonWater, claimed_fmt, lat_lon, 
                  phrase, reportregion,
                  lmp->landdistance, lmp->guesswater, 
                  lmp->claimedfull, lmp->claimeddistance);
      } else {
        ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_LatLonWater, 
                  fmt, lat_lon, 
                  phrase, reportregion, 
                  lmp->landdistance, lmp->guesswater);
      }
    }

  } else if (neardist > 0) {
    fmt = "Lat_lon '%s' is in water '%s', '%s' is %d km away";
    ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_LatLonWater, fmt, lat_lon, lmp->guesswater, fullname, AdjustAndRoundDistance (neardist, scale));
  } else {
    fmt = "Lat_lon '%s' is in water '%s'";
    ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_LatLonWater, fmt, lat_lon, lmp->guesswater);
  }
}


static void LatLonLandErrors (
  ValidStructPtr vsp,
  LatLonMapPtr lmp,
  CharPtr country,
  CharPtr province,
  CharPtr lat_lon,
  CharPtr fullname
  )
{
  CharPtr fmt;

  if (lmp->claimedfull != NULL) {
    fmt = "Lat_lon '%s' maps to '%s' instead of '%s' - claimed region '%s' is at distance %d km";
    if (province != NULL) {
      if (StringICmp (lmp->guesscountry, country) == 0) {
        if (vsp->testLatLonSubregion) {
          ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_LatLonState, fmt, lat_lon, lmp->fullguess, fullname, lmp->claimedfull, lmp->claimeddistance);
        }
      } else {
        ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_LatLonCountry, fmt, lat_lon, lmp->fullguess, fullname, lmp->claimedfull, lmp->claimeddistance);
      }
    } else {
      ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_LatLonCountry, fmt, lat_lon, lmp->fullguess, country, lmp->claimedfull, lmp->claimeddistance);
    }
  } else {
    fmt = "Lat_lon '%s' maps to '%s' instead of '%s'";
    if (StringICmp (lmp->guesscountry, country) == 0) {
      if (vsp->testLatLonSubregion) {
        ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_LatLonState, fmt, lat_lon, lmp->fullguess, fullname);
      }
    } else {
      ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_LatLonCountry, fmt, lat_lon, lmp->fullguess, fullname);
    }
  }
}


typedef enum {
  eLatLonAdjust_none = 0 ,
  eLatLonAdjust_flip = 1 ,
  eLatLonAdjust_negate_lat = 2 ,
  eLatLonAdjust_negate_lon = 4
} ELatLonAdjust;

static void NewerValidateCountryLatLon (
  ValidStructPtr vsp,
  GatherContextPtr gcp,
  CharPtr countryname,
  CharPtr lat_lon
)

{
  Char        buf0 [256], buf1 [256], buf2 [256];
  CharPtr     country = NULL, province = NULL, fullname = NULL;
  CtrySetPtr  csp;
  Boolean     format_ok = FALSE, lat_in_range = FALSE, lon_in_range = FALSE, precision_ok = FALSE;
  FloatHi     lat = 0.0;
  FloatHi     lon = 0.0;
  LatLonMap   llm, adjusted;
  CharPtr     ptr;
  FloatHi     scale = 1.0;
  FloatHi     neardist = 0.0;
  ELatLonAdjust adjust = eLatLonAdjust_none;
  Uint4         test, adjust_test = 0;
  CharPtr       fmt;

  if (vsp == NULL || gcp == NULL) return;
  if (StringHasNoText (countryname)) return;
  if (StringHasNoText (lat_lon)) return;

  IsCorrectLatLonFormat (lat_lon, &format_ok, &precision_ok, &lat_in_range, &lon_in_range);
  if (! format_ok) {
    /* may have comma and then altitude, so just get lat_lon component */
    StringNCpy_0 (buf0, lat_lon, sizeof (buf0));
    ptr = StringChr (buf0, ',');
    if (ptr != NULL) {
      *ptr = '\0';
      lat_lon = buf0;
      IsCorrectLatLonFormat (lat_lon, &format_ok, &precision_ok, &lat_in_range, &lon_in_range);
    }
  }

  /* reality checks - do not bail if only precision issue */
  if (! format_ok) {
    /* incorrect lat_lon format should be reported elsewhere */
    return;
  }
  if (! lat_in_range) {
    /* incorrect latitude range should be reported elsewhere */
    return;
  }
  if (! lon_in_range) {
    /* incorrect longitude range should be reported elsewhere */
    return;
  }

  if (! ParseLatLon (lat_lon, &lat, &lon)) {
    /* report unable to parse lat_lon */
    return;
  }

  StringNCpy_0 (buf1, countryname, sizeof (buf1));
  /* trim at comma or semicolon, leaving only country/ocean and possibly state/province */
  ptr = StringChr (buf1, ',');
  if (ptr != NULL) {
    *ptr = '\0';
  }
  ptr = StringChr (buf1, ';');
  if (ptr != NULL) {
    *ptr = '\0';
  }
  TrimSpacesAroundString (buf1);
  if (StringDoesHaveText (buf1)) {
    fullname = buf1;
  }

  StringNCpy_0 (buf2, buf1, sizeof (buf2));
  /* separate country from state/province */
  ptr = StringChr (buf2, ':');
  if (ptr != NULL) {
    if (CountryIsInLatLonList (buf2)) {
      /* store province if in data list as subregion of designated country */
      *ptr = '\0';
      ptr++;
      TrimSpacesAroundString (ptr);
      if (StringDoesHaveText (ptr)) {
        province = ptr;
      }
    } else {
      /* otherwise just truncate country at colon, trimming further descriptive information */
      *ptr = '\0';
      ptr++;
    }
  }
  TrimSpacesAroundString (buf2);
  if (StringDoesHaveText (buf2)) {
    country = buf2;
  }

  if (StringHasNoText (country)) {
    /* report leading colon without country */
    return;
  }

  /* known exceptions - don't even bother calculating any further */
  if (StringCmp (country, "Antarctica") == 0 && lat < -60.0) {
    return;
  }

  if (! CountryIsInLatLonList (country)) {
    if (! WaterIsInLatLonList (country)) {
      /* report unrecognized country */
      return;
    } else {
      /* report that it may refer to specific small body of water */
      /* continue to look for nearby country for proximity report */
      /* (do not return) */
    }
  }

  csp = GetLatLonCountryData ();
  if (csp == NULL) {
    /* report unable to find data */
    return;
  }

  /* scale (reciprocal of degree resolution) needed for adjusting offshore distance calculation */
  scale = csp->scale;

  /* calculate assignment or proximity by coordinates */
  CalculateLatLonMap (lat, lon, country, province, scale, &llm);

  /* compare indicated country/province to guess/proximate country/water */
  test = ClassifyLatLonMap (fullname, country, province, &llm);

  if (!test && CountryIsNearLatLon(country, lat, lon, 2.0, &neardist) && neardist < 5.0) {
    llm.guesscountry = country;
    llm.guessprovince = NULL;
    test = ClassifyLatLonMap (fullname, country, province, &llm);
  }

  if (!test && !CountryIsNearLatLon(country, lat, lon, 20.0, &neardist) && !WaterIsNearLatLon(country, lat, lon, 20.0, &neardist)) {
    CalculateLatLonMap (lon, lat, country, province, scale, &adjusted);
    adjust_test = ClassifyLatLonMap (fullname, country, province, &adjusted);
    if (adjust_test) {
      adjust = eLatLonAdjust_flip;
    } else {
      CalculateLatLonMap (-lat, lon, country, province, scale, &adjusted);
      adjust_test = ClassifyLatLonMap (fullname, country, province, &adjusted);
      if (adjust_test) {
        adjust = eLatLonAdjust_negate_lat;
      } else {
        CalculateLatLonMap (lat, -lon, country, province, scale, &adjusted);
        adjust_test = ClassifyLatLonMap (fullname, country, province, &adjusted);
        if (adjust_test) {
          adjust = eLatLonAdjust_negate_lon;
        }
      }
    }

    if (adjust_test) {
      test = adjust_test;
      MemCopy (&llm, &adjusted, sizeof (LatLonMap));
    }
  }

  if (adjust) {
    if (adjust == eLatLonAdjust_flip) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_LatLonValue, "Latitude and longitude values appear to be exchanged");
    } else if (adjust == eLatLonAdjust_negate_lat) {
      if (lat < 0.0) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_LatLonValue, "Latitude should be set to N (northern hemisphere)");
      } else {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_LatLonValue, "Latitude should be set to S (southern hemisphere)");
      }
    } else if (adjust == eLatLonAdjust_negate_lon) {
      if (lon < 0.0) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_LatLonValue, "Longitude should be set to E (eastern hemisphere)");
      } else {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_LatLonValue, "Longitude should be set to W (western hemisphere)");
      }
    }
  } else {
    if ((test & eLatLonClassify_CountryMatch) && (test & eLatLonClassify_ProvinceMatch)) {
      /* success!  nothing to report */
    } else if (test & eLatLonClassify_WaterMatch) {
      /* success!  nothing to report */
    } else if (test & eLatLonClassify_CountryMatch && province == NULL) {
      if (vsp->testLatLonSubregion) {
        fmt = "Lat_lon %s is in %s (more specific than %s)";
        ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_LatLonState, fmt, lat_lon, llm.fullguess, country);
      }
    } else if (llm.guesswater != NULL) {
      LatLonWaterErrors(vsp, &llm, test, neardist, country, province, lat_lon, fullname, scale);
    } else if (llm.guesscountry != NULL) {
      LatLonLandErrors (vsp, &llm, country, province, lat_lon, fullname);
    } else if (llm.closestcountry != NULL) {
      fmt = "Lat_lon '%s' is closest to '%s' instead of '%s'";
      ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_LatLonCountry, fmt, lat_lon, llm.closestcountry, fullname);
    } else if (llm.closestwater != NULL) {
      fmt = "Lat_lon '%s' is closest to '%s' instead of '%s'";
      ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_LatLonWater, fmt, lat_lon, llm.closestwater, fullname);
    } else {
      fmt = "Unable to determine mapping for lat_lon '%s' and country '%s'";
      ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_LatLonCountry, fmt, lat_lon, fullname);
    }
  }
}


/* note - special case for sex because it prevents a different message from being displayed, do not list here */
static const Uint1 sUnexpectedViralSubSourceQualifiers[] = {
  SUBSRC_cell_line, 
  SUBSRC_cell_type, 
  SUBSRC_tissue_type,
  SUBSRC_dev_stage
};

static const Int4 sNumUnexpectedViralSubSourceQualifiers = sizeof (sUnexpectedViralSubSourceQualifiers) / sizeof (Uint1);


static Boolean IsUnexpectedViralSubSourceQualifier (Uint1 subtype)
{
  Int4 i;
  Boolean rval = FALSE;

  for (i = 0; i < sNumUnexpectedViralSubSourceQualifiers && !rval; i++) {
    if (subtype == sUnexpectedViralSubSourceQualifiers[i]) {
      rval = TRUE;
    }
  }
  return rval;
}

static const Uint1 sUnexpectedViralOrgModQualifiers[] = {
  ORGMOD_breed,
  ORGMOD_cultivar,
  ORGMOD_specimen_voucher
};

static const Int4 sNumUnexpectedViralOrgModQualifiers = sizeof (sUnexpectedViralOrgModQualifiers) / sizeof (Uint1);


static Boolean IsUnexpectedViralOrgModQualifier (Uint1 subtype)
{
  Int4 i;
  Boolean rval = FALSE;

  for (i = 0; i < sNumUnexpectedViralOrgModQualifiers && !rval; i++) {
    if (subtype == sUnexpectedViralOrgModQualifiers[i]) {
      rval = TRUE;
    }
  }
  return rval;
}


static void ValidateBioSource (ValidStructPtr vsp, GatherContextPtr gcp, BioSourcePtr biop, SeqFeatPtr sfp, ValNodePtr sdp)
{
  Char            badch;
  Boolean         bad_cap = FALSE;
  Boolean         bad_frequency;
  BioseqPtr       bsp;
  BioseqSetPtr    bssp;
  Char            ch;
  Boolean         chromconf = FALSE;
  Int2            chromcount = 0;
  SubSourcePtr    chromosome = NULL;
  CharPtr         countryname = NULL;
  CtrySetPtr      csp;
  ValNodePtr      db;
  DbtagPtr        dbt;
  CharPtr         gb_synonym = NULL;
  Boolean         germline = FALSE;
  CharPtr         good;
  Boolean         has_strain = FALSE;
  Boolean         has_fwd_pcr_seq = FALSE;
  Boolean         has_rev_pcr_seq = FALSE;
  Boolean         has_pcr_name = FALSE;
  Boolean         has_metagenome_source = FALSE;
  Boolean         has_plasmid = FALSE;
  Int4            id;
  Boolean         is_env_sample = FALSE;
  Boolean         is_iso_source = FALSE;
  Boolean         is_mating_type = FALSE;
  Boolean         is_metagenomic = FALSE;
  Boolean         is_sex = FALSE;
  Boolean         is_specific_host = FALSE;
  Boolean         is_transgenic = FALSE;
  Boolean         isAnimal = FALSE;
  Boolean         isArchaea = FALSE;
  Boolean         isBacteria = FALSE;
  Boolean         isFungal = FALSE;
  Boolean         isPlant = FALSE;
  Boolean         isViral = FALSE;
  Boolean         is_bc;
  Boolean         is_rf;
  Boolean         is_sc;
  CharPtr         last_db = NULL;
  CharPtr         lat_lon = NULL;
  Int2            num_bio_material = 0;
  Int2            num_culture_collection = 0;
  Int2            num_specimen_voucher = 0;
  Int2            num_country = 0;
  Int2            num_lat_lon = 0;
  Int2            num_fwd_primer_seq = 0;
  Int2            num_rev_primer_seq = 0;
  Int2            num_fwd_primer_name = 0;
  Int2            num_rev_primer_name = 0;
  ObjectIdPtr     oip;
  Boolean         old_country = FALSE;
  OrgNamePtr      onp;
  OrgModPtr       omp, nxtomp;
  OrgRefPtr       orp;
  ObjValNodePtr   ovp;
  Int4            primer_len_before;
  Int4            primer_len_after;
  ValNodePtr      pset;
  Boolean         rearranged = FALSE;
  SeqEntryPtr     sep;
  ErrSev          sev;
  SubSourcePtr    ssp;
  CharPtr         str;
  CharPtr         synonym = NULL;
  Boolean         varietyOK;
  CharPtr         inst1, inst2, id1, id2, coll1, coll2;
  Char            buf1 [512], buf2 [512];
  PCRPrimerPtr      ppp;
  PCRReactionSetPtr prp;

  if (vsp->sourceQualTags == NULL) {
    InitializeSourceQualTags (vsp);
  }
  if (biop == NULL)
    return;
  if (biop->genome == GENOME_transposon || biop->genome == GENOME_insertion_seq) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_ObsoleteSourceLocation,
              "Transposon and insertion sequence are no longer legal locations");
  }

  if (vsp->indexerVersion && biop->genome == GENOME_chromosome) {
    ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_ChromosomeLocation, "INDEXER_ONLY - BioSource location is chromosome");
  }

  orp = biop->org;
  if (orp != NULL) {
    onp = orp->orgname;
    if (onp != NULL) {
      if (StringNICmp (onp->lineage, "Viruses; ", 9) == 0) {
        isViral = TRUE;
      } else if (StringNICmp (onp->lineage, "Eukaryota; Metazoa; ", 20) == 0) {
        isAnimal = TRUE;
      } else if (StringNICmp (onp->lineage, "Eukaryota; Viridiplantae; Streptophyta; Embryophyta; ", 53) == 0 ||
                 StringNICmp (onp->lineage, "Eukaryota; Rhodophyta; ", 23) == 0 ||
                 StringNICmp (onp->lineage, "Eukaryota; stramenopiles; Phaeophyceae; ", 40) == 0) {
        isPlant = TRUE;
      } else if (StringNICmp (onp->lineage, "Bacteria; ", 10) == 0) {
        isBacteria = TRUE;
      } else if (StringNICmp (onp->lineage, "Archaea; ", 9) == 0) {
        isArchaea = TRUE;
      } else if (StringNICmp (onp->lineage, "Eukaryota; Fungi; ", 18) == 0) {
        isFungal = TRUE;
      }
    }
  }

  ssp = biop->subtype;
  while (ssp != NULL) {
    if (ssp->subtype == SUBSRC_country) {
      num_country++;
      if (countryname != NULL) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadCountryCode, "Multiple country names on BioSource");
      }
      countryname = ssp->name;
      if (CountryIsValid (countryname, &old_country, &bad_cap)) {
        if (bad_cap) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadCountryCapitalization, "Bad country capitalization [%s]", countryname);
        }
      } else {
        if (StringHasNoText (countryname)) {
          countryname = "?";
        }
        if (old_country) {
          ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_ReplacedCountryCode, "Replaced country name [%s]", countryname);
        } else {
          /*
          sev = SEV_ERROR;
          if (vsp->is_barcode_sep && vsp->seqSubmitParent) {
            sev = SEV_WARNING;
          }
          */
          sev = SEV_WARNING;
          ValidErr (vsp, sev, ERR_SEQ_DESCR_BadCountryCode, "Bad country name [%s]", countryname);
        }
      }
    } else if (ssp->subtype == SUBSRC_chromosome) {
      chromcount++;
      if (chromosome != NULL) {
        if (StringICmp (ssp->name, chromosome->name) != 0) {
          chromconf = TRUE;
        }
      } else {
        chromosome = ssp;
      }
    } else if (ssp->subtype == SUBSRC_transposon_name || ssp->subtype == SUBSRC_insertion_seq_name) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_ObsoleteSourceQual,
                "Transposon name and insertion sequence name are no longer legal qualifiers");
    } else if (ssp->subtype == 0) {
      ValidErr (vsp, SEV_REJECT, ERR_SEQ_DESCR_BadSubSource, "Unknown subsource subtype %d", (int) (ssp->subtype));
    } else if (ssp->subtype == SUBSRC_other) {
      ValidateSourceQualTags (vsp, gcp, biop, ssp->name);
    } else if (ssp->subtype == SUBSRC_germline) {
      germline = TRUE;
      str = ssp->name;
      if (str == NULL || str [0] != '\0') {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Germline qualifier should not have descriptive text");
      }
    } else if (ssp->subtype == SUBSRC_rearranged) {
      rearranged = TRUE;
      str = ssp->name;
      if (str == NULL || str [0] != '\0') {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Rearranged qualifier should not have descriptive text");
      }
    } else if (ssp->subtype == SUBSRC_transgenic) {
      is_transgenic = TRUE;
      str = ssp->name;
      if (str == NULL || str [0] != '\0') {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Transgenic qualifier should not have descriptive text");
      }
    } else if (ssp->subtype == SUBSRC_environmental_sample) {
      is_env_sample = TRUE;
      str = ssp->name;
      if (str == NULL || str [0] != '\0') {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Environmental_sample qualifier should not have descriptive text");
      }
    } else if (ssp->subtype == SUBSRC_metagenomic) {
      is_metagenomic = TRUE;
      str = ssp->name;
      if (str == NULL || str [0] != '\0') {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Metagenomic qualifier should not have descriptive text");
      }
    } else if (ssp->subtype == SUBSRC_isolation_source) {
      is_iso_source = TRUE;
    } else if (ssp->subtype == SUBSRC_sex) {
      is_sex = TRUE;
      str = ssp->name;
      if (isAnimal || isPlant) {
        /* always use /sex, do not check values at this time */
      } else if (isViral) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Virus has unexpected sex qualifier");
      } else if (isBacteria || isArchaea || isFungal) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Unexpected use of /sex qualifier");
      } else if (IsValidSexValue (str)) {
        /* otherwise expect male or female, or a few others */
      } else {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Unexpected use of /sex qualifier");
      }
    } else if (ssp->subtype == SUBSRC_mating_type) {
      is_mating_type = TRUE;
      str = ssp->name;
      if (isAnimal || isPlant || isViral) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Unexpected use of /mating_type qualifier");
      } else if (IsValidSexValue (str)) {
        /* complain if one of the values that should go in /sex */
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Unexpected use of /mating_type qualifier");
      }
    } else if (ssp->subtype == SUBSRC_plasmid_name) {
      has_plasmid = TRUE;
      if (biop->genome != GENOME_plasmid) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Plasmid subsource but not plasmid location");
      }
    } else if (ssp->subtype == SUBSRC_plastid_name) {
      if (StringCmp (ssp->name, "chloroplast") == 0) {
        if (biop->genome != GENOME_chloroplast) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Plastid name subsource chloroplast but not chloroplast location");
        }
      } else if (StringCmp (ssp->name, "chromoplast") == 0) {
        if (biop->genome != GENOME_chromoplast) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Plastid name subsource chromoplast but not chromoplast location");
        }
      } else if (StringCmp (ssp->name, "kinetoplast") == 0) {
        if (biop->genome != GENOME_kinetoplast) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Plastid name subsource kinetoplast but not kinetoplast location");
        }
      } else if (StringCmp (ssp->name, "plastid") == 0) {
        if (biop->genome != GENOME_plastid) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Plastid name subsource plastid but not plastid location");
        }
      } else if (StringCmp (ssp->name, "apicoplast") == 0) {
        if (biop->genome != GENOME_apicoplast) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Plastid name subsource apicoplast but not apicoplast location");
        }
      } else if (StringCmp (ssp->name, "leucoplast") == 0) {
        if (biop->genome != GENOME_leucoplast) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Plastid name subsource leucoplast but not leucoplast location");
        }
      } else if (StringCmp (ssp->name, "proplastid") == 0) {
        if (biop->genome != GENOME_proplastid) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Plastid name subsource proplastid but not proplastid location");
        }
      } else if (StringCmp (ssp->name, "chromatophore") == 0) {
        if (biop->genome != GENOME_chromatophore) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Plastid name subsource chromatophore but not chromatophore location");
        }
      } else {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Plastid name subsource contains unrecognized value");
      }
    } else if (ssp->subtype == SUBSRC_collection_date) {
      if (! CollectionDateIsValid (ssp->name)) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadCollectionDate, "Collection_date format is not in DD-Mmm-YYYY format");
      }
      else if (CollectionDateIsInTheFuture (ssp->name)) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadCollectionDate, "Collection_date is in the future");
      }
    } else if (ssp->subtype == SUBSRC_fwd_primer_seq) {
      num_fwd_primer_seq++;
      has_fwd_pcr_seq = TRUE;
      if (! PrimerSeqIsValid (vsp, ssp->name, &badch)) {
        if (badch < ' ' || badch > '~') {
          badch = '?';
        }
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadPCRPrimerSequence,
                  "PCR forward primer sequence format is incorrect, first bad character is '%c'", (char) badch);
      }
      /*
      if (PrimerSeqHasDuplicates (ssp->name)) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_DuplicatePCRPrimerSequence,
                  "PCR forward primer sequence has duplicates");
      }
      */
    } else if (ssp->subtype == SUBSRC_rev_primer_seq) {
      num_rev_primer_seq++;
      has_rev_pcr_seq = TRUE;
      if (! PrimerSeqIsValid (vsp, ssp->name, &badch)) {
        if (badch < ' ' || badch > '~') {
          badch = '?';
        }
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadPCRPrimerSequence,
                  "PCR reverse primer sequence format is incorrect, first bad character is '%c'", (char) badch);
      }
      /*
      if (PrimerSeqHasDuplicates (ssp->name)) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_DuplicatePCRPrimerSequence,
                  "PCR reverse primer sequence has duplicates");
      }
      */
    } else if (ssp->subtype == SUBSRC_fwd_primer_name) {
      num_fwd_primer_name++;
      if (StringLen (ssp->name) > 10 && PrimerSeqIsValid (vsp, ssp->name, &badch)) {
        if (badch < ' ' || badch > '~') {
          badch = '?';
        }
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadPCRPrimerName, "PCR primer name appears to be a sequence");
      }
      has_pcr_name = TRUE;
    } else if (ssp->subtype == SUBSRC_rev_primer_name) {
      num_rev_primer_name++;
      if (StringLen (ssp->name) > 10 && PrimerSeqIsValid (vsp, ssp->name, &badch)) {
        if (badch < ' ' || badch > '~') {
          badch = '?';
        }
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadPCRPrimerName, "PCR primer name appears to be a sequence");
      }
      has_pcr_name = TRUE;
    } else if (ssp->subtype == SUBSRC_lat_lon) {
      num_lat_lon++;
      if (lat_lon != NULL) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_LatLonProblem, "Multiple lat_lon on BioSource");
      }
      lat_lon = ssp->name;
      ValidateLatLon (vsp, lat_lon);
    } else if (ssp->subtype == SUBSRC_frequency) {
      str = ssp->name;
      if (StringDoesHaveText (str)) {
        bad_frequency = FALSE;
        if (StringCmp (str, "0") == 0) {
          /* ignore */
        } else if (StringCmp (str, "1") == 0) {
          ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_BioSourceInconsistency, "bad frequency qualifier value %s", ssp->name);
        } else {
          ch = *str;
          if (ch == '0') {
            str++;
            ch = *str;
          }
          if (ch == '.') {
            str++;
            ch = *str;
            if (! IS_DIGIT (ch)) {
              bad_frequency = TRUE;
            } else {
              while (ch != '\0') {
                if (! IS_DIGIT (ch)) {
                  bad_frequency = TRUE;
                }
                str++;
                ch = *str;
              }
            }
          } else {
            bad_frequency = TRUE;
          }
        }
        if (bad_frequency) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "bad frequency qualifier value %s", ssp->name);
        }
      }
    }

    if (isViral && IsUnexpectedViralSubSourceQualifier(ssp->subtype)) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Virus has unexpected %s qualifier", GetSubsourceQualName (ssp->subtype));
    }
    ssp = ssp->next;
  }

  if (biop->genome == GENOME_plasmid) {
    if (! has_plasmid) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Plasmid location but not plasmid subsource");
    }
  }

  if (num_country > 1) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MultipleSourceQualifiers, "Multiple country qualifiers present");
  }
  if (num_lat_lon > 1) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MultipleSourceQualifiers, "Multiple lat_lon qualifiers present");
  }
  if (num_fwd_primer_seq > 1) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MultipleSourceQualifiers, "Multiple fwd_primer_seq qualifiers present");
  }
  if (num_rev_primer_seq > 1) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MultipleSourceQualifiers, "Multiple rev_primer_seq qualifiers present");
  }
  if (num_fwd_primer_name > 1) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MultipleSourceQualifiers, "Multiple fwd_primer_name qualifiers present");
  }
  if (num_rev_primer_name > 1) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MultipleSourceQualifiers, "Multiple rev_primer_name qualifiers present");
  }

  if (countryname != NULL && lat_lon != NULL) {
    csp = GetLatLonCountryData ();
    if (csp != NULL) {
      NewerValidateCountryLatLon (vsp, gcp, countryname, lat_lon);
    }
  }

  if (has_pcr_name) {
    if ((! has_fwd_pcr_seq) || (! has_rev_pcr_seq)) {
      /*
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadPCRPrimerSequence, "PCR primer has name but not both sequences");
      */
    }
  } else if (has_fwd_pcr_seq || has_rev_pcr_seq) {
    if (! (has_fwd_pcr_seq && has_rev_pcr_seq)) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadPCRPrimerSequence, "PCR primer does not have both sequences");
    }
  }

  pset = ParsePCRSet (biop);
  if (pset != NULL) {
    pset = ValNodeSort (pset, SortVnpByPCRSetSeq);
    primer_len_before = ValNodeLen (pset);
    pset = UniqueVnpByPCRSetSeq (pset);
    primer_len_after = ValNodeLen (pset);
    FreePCRSet (pset);
    if (primer_len_before != primer_len_after) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_DuplicatePCRPrimerSequence,
                "PCR primer sequence has duplicates");
    }
  }

  for (prp = biop->pcr_primers; prp != NULL; prp = prp->next) {

    for (ppp = prp->forward; ppp != NULL; ppp = ppp->next) {
      if (StringDoesHaveText (ppp->seq) && (! PrimerSeqIsValid (vsp, ppp->seq, &badch))) {
        if (badch < ' ' || badch > '~') {
          badch = '?';
        }
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadPCRPrimerSequence,
                  "PCR forward primer sequence format is incorrect, first bad character is '%c'", (char) badch);
      }
      if (StringLen (ppp->name) > 10 && PrimerSeqIsValid (vsp, ppp->name, &badch)) {
        if (badch < ' ' || badch > '~') {
          badch = '?';
        }
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadPCRPrimerName, "PCR forward primer name appears to be a sequence");
      }
    }

    for (ppp = prp->reverse; ppp != NULL; ppp = ppp->next) {
      if (StringDoesHaveText (ppp->seq) && (! PrimerSeqIsValid (vsp, ppp->seq, &badch))) {
        if (badch < ' ' || badch > '~') {
          badch = '?';
        }
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadPCRPrimerSequence,
                  "PCR reverse primer sequence format is incorrect, first bad character is '%c'", (char) badch);
      }
      if (StringLen (ppp->name) > 10 && PrimerSeqIsValid (vsp, ppp->name, &badch)) {
        if (badch < ' ' || badch > '~') {
          badch = '?';
        }
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadPCRPrimerName, "PCR reverse primer name appears to be a sequence");
      }
    }
  }

  if (germline && rearranged) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Germline and rearranged should not both be present");
  }
  if (is_transgenic && is_env_sample) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Transgenic and environmental sample should not both be present");
  }
  if (is_metagenomic && (! is_env_sample)) {
    ValidErr (vsp, SEV_REJECT, ERR_SEQ_DESCR_BioSourceInconsistency, "Metagenomic should also have environmental sample annotated");
  }
  if (is_sex && is_mating_type) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Sex and mating type should not both be present");
  }

  if (biop->org != NULL
      && biop->org->orgname != NULL
      && StringISearch (biop->org->orgname->lineage, "metagenomes") != NULL
      && !is_metagenomic) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "If metagenomes appears in lineage, BioSource should have metagenomic qualifier");
  }
  if (chromcount > 1) {
    if (chromconf) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MultipleChromosomes, "Multiple conflicting chromosome qualifiers");
    } else {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MultipleChromosomes, "Multiple identical chromosome qualifiers");
    }
  }
  orp = biop->org;
  if (orp != NULL) {
    /*
    if (StringICmp (orp->taxname, "Human immunodeficiency virus") == 0 ||
        StringICmp (orp->taxname, "Human immunodeficiency virus 1") == 0 ||
        StringICmp (orp->taxname, "Human immunodeficiency virus 2") == 0) {
      ValidateLocationForHIV (vsp, biop);
    } else */
    if (StringICmp (orp->taxname, "uncultured bacterium") == 0) {
      bsp = NULL;
      if (sfp != NULL) {
        bsp = BioseqFindFromSeqLoc (sfp->location);
      } else if (sdp != NULL && sdp->extended != 0) {
        ovp = (ObjValNodePtr) sdp;
        if (ovp->idx.parenttype == OBJ_BIOSEQ) {
          bsp = (BioseqPtr) ovp->idx.parentptr;
        } else if (ovp->idx.parenttype == OBJ_BIOSEQSET) {
          bssp = (BioseqSetPtr) ovp->idx.parentptr;
          if (bssp != NULL) {
            sep = bssp->seqentry;
            if (sep != NULL) {
              sep = FindNthBioseq (sep, 1);
              if (sep != NULL && IS_Bioseq (sep)) {
                bsp = (BioseqPtr) sep->data.ptrvalue;
              }
            }
          }
        }
      }
      if (bsp != NULL && bsp->length >= 10000) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Uncultured bacterium sequence length is suspiciously high");
      }
    }
    if (StringNICmp (orp->taxname, "uncultured ", 11) == 0) {
      if (! is_env_sample) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Uncultured should also have /environmental_sample");
      }
    }
    str = orp->taxname;
    if (StringDoesHaveText (str)) {
      if (UnbalancedParentheses (str)) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_UnbalancedParentheses,
                    "Unbalanced parentheses in taxname '%s'", str);
      }
      if (StringHasSgml (vsp, str)) {
        ValidErr (vsp, SEV_WARNING, ERR_GENERIC_SgmlPresentInText, "taxname %s has SGML", str);
      }
    }
  }

  for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
    if (ssp->subtype == SUBSRC_germline ||
              ssp->subtype == SUBSRC_rearranged ||
              ssp->subtype == SUBSRC_transgenic ||
              ssp->subtype == SUBSRC_environmental_sample ||
              ssp->subtype == SUBSRC_metagenomic) continue;
    str = ssp->name;
    if (StringHasNoText (str)) continue;
    if (UnbalancedParentheses (str)) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_UnbalancedParentheses,
                  "Unbalanced parentheses in subsource '%s'", str);
    }
    if (StringHasSgml (vsp, str)) {
      ValidErr (vsp, SEV_WARNING, ERR_GENERIC_SgmlPresentInText, "subsource %s has SGML", str);
    }
  }

  if (orp == NULL || (StringHasNoText (orp->taxname) && StringHasNoText (orp->common))) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_NoOrgFound, "No organism name has been applied to this Bioseq.  Other qualifiers may exist.");
  }
  if (orp == NULL) {
    if (is_env_sample && (! is_iso_source) && (! is_specific_host)) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Environmental sample should also have isolation source or specific host annotated");
    }
    return;
  }
  onp = orp->orgname;
  if (onp == NULL || StringHasNoText (onp->lineage)) {
    if (! vsp->seqSubmitParent && vsp->indexerVersion) { /* suppress when validator run from tbl2asn or when not indexer version */
      sev = SEV_ERROR;
      if (vsp->is_refseq_in_sep) {
        for (db = orp->db; db != NULL; db = db->next) {
          dbt = (DbtagPtr) db->data.ptrvalue;
          if (dbt != NULL) {
            if (StringICmp (dbt->db, "taxon") == 0) {
              sev = SEV_REJECT;
            }
          }
        }
      }
      if (vsp->is_embl_ddbj_in_sep) {
        sev = SEV_WARNING;
      }
      ValidErr (vsp, sev, ERR_SEQ_DESCR_MissingLineage, "No lineage for this BioSource.");
    }
  } else {
    if (biop->genome == GENOME_kinetoplast) {
      if (StringStr (onp->lineage, "Kinetoplastida") == 0) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadOrganelle, "Only Kinetoplastida have kinetoplasts");
      }
    } else if (biop->genome == GENOME_nucleomorph) {
      if (StringStr (onp->lineage, "Chlorarachniophyceae") == 0 && StringStr (onp->lineage, "Cryptophyta") == 0) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadOrganelle, "Only Chlorarachniophyceae and Cryptophyta have nucleomorphs");
      }
    } else if (biop->genome == GENOME_macronuclear) {
      if (StringStr (onp->lineage, "Ciliophora") == 0) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadOrganelle, "Only Ciliophora have macronuclear locations");
      }
    }

    /* warn if bacteria has organelle location */
    if (StringCmp (onp->div, "BCT") == 0 || StringCmp (onp->div, "VRL") == 0) {
      if (biop->genome == GENOME_unknown
          || biop->genome == GENOME_genomic
          || biop->genome == GENOME_plasmid
          || biop->genome == GENOME_chromosome
          || (biop->genome == GENOME_proviral && StringCmp (onp->div, "VRL") == 0)) {
        /* it's ok */
      } else {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Bacterial or viral source should not have organelle location");
      }
    }

    if (StringCmp (onp->div, "ENV") == 0 && (! is_env_sample)) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_BioSourceInconsistency, "BioSource with ENV division is missing environmental sample subsource");
    }
  }
  for (db = orp->db; db != NULL; db = db->next) {
    dbt = (DbtagPtr) db->data.ptrvalue;
    if (dbt != NULL) {
      if (last_db != NULL) {
        if (StringICmp (dbt->db, last_db) == 0) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceDbTagConflict, "BioSource uses db %s multiple times", last_db);
        }
      }
      last_db = dbt->db;
    }
  }
  if (onp != NULL) {
    omp = onp->mod;
    varietyOK = FALSE;
    while (omp != NULL) {
      if (omp->subtype == 0 || omp->subtype == 1) {
        ValidErr (vsp, SEV_REJECT, ERR_SEQ_DESCR_BadOrgMod, "Unknown orgmod subtype %d", (int) (omp->subtype));
      } else if (omp->subtype == ORGMOD_strain) {
        if (has_strain) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadOrgMod, "Multiple strain qualifiers on the same BioSource");
        }
        has_strain = TRUE;
      } else if (omp->subtype == ORGMOD_variety) {
        if ((StringHasNoText (onp->div) || StringICmp (onp->div, "PLN") != 0) &&
            StringStr (onp->lineage, "Cyanobacteria") == 0 &&
            StringStr (onp->lineage, "Myxogastria") == 0 &&
            StringStr (onp->lineage, "Oomycetes") == 0) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadOrgMod, "Orgmod variety should only be in plants, fungi, or cyanobacteria");
        }
        varietyOK = ValidateOrgModInTaxName (vsp, omp, orp->taxname, varietyOK);
      } else if (omp->subtype == ORGMOD_nat_host) {
        is_specific_host = TRUE;
        if (StringICmp (omp->subname, orp->taxname) == 0) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadOrgMod, "Specific host is identical to taxname");
        }
      } else if (omp->subtype == ORGMOD_other) {
        ValidateSourceQualTags (vsp, gcp, biop, omp->subname);
      } else if (omp->subtype == ORGMOD_biovar ||
                 omp->subtype == ORGMOD_forma ||
                 omp->subtype == ORGMOD_forma_specialis ||
                 omp->subtype == ORGMOD_sub_species ||
                 omp->subtype == ORGMOD_pathovar) {
        ValidateOrgModInTaxName (vsp, omp, orp->taxname, varietyOK);
      } else if (omp->subtype == ORGMOD_specimen_voucher) {
        num_specimen_voucher++;
        ValidateOrgModVoucher (vsp, omp);
      } else if (omp->subtype == ORGMOD_culture_collection) {
        num_culture_collection++;
        ValidateOrgModVoucher (vsp, omp);
      } else if (omp->subtype == ORGMOD_bio_material) {
        num_bio_material++;
        ValidateOrgModVoucher (vsp, omp);
      } else if (omp->subtype == ORGMOD_metagenome_source) {
        has_metagenome_source = TRUE;
      } else if (omp->subtype == ORGMOD_common) {
        if (StringICmp (omp->subname, orp->common) == 0) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadOrgMod, "OrgMod common is identical to Org-ref common");
        }
      } else if (omp->subtype == ORGMOD_synonym) {
        synonym = omp->subname;
      } else if (omp->subtype == ORGMOD_gb_synonym) {
        gb_synonym = omp->subname;
      }

      if (isViral && IsUnexpectedViralOrgModQualifier(omp->subtype)) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Virus has unexpected %s qualifier", GetOrgModQualName (omp->subtype));
      }

      omp = omp->next;
    }

    for (omp = onp->mod; omp != NULL; omp = omp->next) {
      if (omp->subtype != ORGMOD_specimen_voucher &&
          omp->subtype != ORGMOD_culture_collection &&
          omp->subtype != ORGMOD_bio_material) continue;
      nxtomp = omp->next;
      if (nxtomp == NULL) continue;
      inst1 = NULL;
      inst2 = NULL;
      id1 = NULL;
      id2 = NULL;
      coll1 = NULL;
      coll2 = NULL;
      StringNCpy_0 (buf1, omp->subname, sizeof (buf1));
      StringNCpy_0 (buf2, nxtomp->subname, sizeof (buf2));
      if (StringChr (buf1, ':') == NULL || StringChr (buf2, ':') == NULL) continue;
      if (! ParseStructuredVoucher (buf1, &inst1, &id1)) continue;
      if (! ParseStructuredVoucher (buf2, &inst2, &id2)) continue;
      if (inst1 == NULL || inst2 == NULL) continue;
      if (StringNICmp (inst1, "personal", 8) == 0) continue;
      if (StringNICmp (inst2, "personal", 8) == 0) continue;
      coll1 = StringChr (inst1, ':');
      if (coll1 != NULL) {
        *coll1 = '\0';
        coll1++;
      }
      coll2 = StringChr (inst2, ':');
      if (coll2 != NULL) {
        *coll2 = '\0';
        coll2++;
      }
      if (StringICmp (inst1, inst2) != 0) continue;
      if (omp->subtype != nxtomp->subtype) continue;
      if (StringCmp (coll1, "DNA") == 0 || StringCmp (coll2, "DNA") == 0) continue;
      if (coll1 != NULL && coll2 != NULL && StringICmp (coll1, coll2) == 0) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MultipleSourceVouchers, "Multiple vouchers with same institution:collection");
      } else {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MultipleSourceVouchers, "Multiple vouchers with same institution");
      }
    }
  }

  if (onp != NULL) {
    for (omp = onp->mod; omp != NULL; omp = omp->next) {
      str = omp->subname;
      if (StringHasNoText (str)) continue;
      if (UnbalancedParentheses (str)) {
        if (omp->subtype == ORGMOD_old_name) continue;
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_UnbalancedParentheses,
                  "Unbalanced parentheses in orgmod '%s'", str);
      }
      if (StringHasSgml (vsp, str)) {
        ValidErr (vsp, SEV_WARNING, ERR_GENERIC_SgmlPresentInText, "orgmod %s has SGML", str);
      }
    }
  }

  /*
  if (num_bio_material > 1) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MultipleSourceQualifiers, "Multiple bio_material qualifiers present");
  }
  if (num_culture_collection > 1) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MultipleSourceQualifiers, "Multiple culture_collection qualifiers present");
  }
  if (num_specimen_voucher > 1) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MultipleSourceQualifiers, "Multiple specimen_voucher qualifiers present");
  }
  */
  if (is_env_sample && has_strain) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_BioSourceInconsistency, "Strain should not be present in an environmental sample");
  }
  if (is_env_sample && (! is_iso_source) && (! is_specific_host)) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Environmental sample should also have isolation source or specific host annotated");
  }
  if (has_metagenome_source && (! is_metagenomic)) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_BioSourceInconsistency, "Metagenome source should also have metagenomic qualifier");
  }
  if (StringDoesHaveText (synonym) && StringDoesHaveText (gb_synonym)) {
    if (StringICmp (synonym, gb_synonym) == 0) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "OrgMod synonym is identical to OrgMod gb_synonym");
    }
  }

  for (db = orp->db; db != NULL; db = db->next) {
    id = -1;
    dbt = (DbtagPtr) db->data.ptrvalue;
    if (dbt != NULL && dbt->db != NULL) {

      if (DbxrefIsValid (dbt->db, &is_rf, &is_sc, &is_bc, &good)) {
        if (is_bc) {
          if (StringHasNoText (good)) {
            good = "?";
          }
          if (is_sc) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_IllegalDbXref,
                      "Illegal db_xref type %s, legal capitalization is %s", dbt->db, good);
          } else {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_IllegalDbXref,
                      "Illegal db_xref type %s, legal capitalization is %s, but should not be used on an OrgRef",
                      dbt->db, good);
          }
        } else if (is_rf) {
          if (vsp->is_refseq_in_sep || vsp->is_gps_in_sep) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_IllegalDbXref,
                      "RefSeq-specific db_xref type %s should not used on an OrgRef", dbt->db);
          } else {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_IllegalDbXref,
                      "RefSeq-specific db_xref type %s should not used on a non-RefSeq OrgRef", dbt->db);
          }
        } else if (is_sc) {
        } else {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_IllegalDbXref,
                    "db_xref type %s should not used on an OrgRef", dbt->db);
        }
      } else {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_IllegalDbXref, "Illegal db_xref type %s", dbt->db);
      }

      if (StringDoesHaveText (dbt->db)) {
        if (StringHasSgml (vsp, dbt->db)) {
          ValidErr (vsp, SEV_WARNING, ERR_GENERIC_SgmlPresentInText, "dbxref database %s has SGML", dbt->db);
        }
      }

      oip = dbt->tag;
      if (oip != NULL && StringDoesHaveText (oip->str)) {
        if (StringHasSgml (vsp, oip->str)) {
          ValidErr (vsp, SEV_WARNING, ERR_GENERIC_SgmlPresentInText, "dbxref value %s has SGML", oip->str);
        }
      }

      /*
      dbxerr = NULL;
      dbvalid = IsDbxrefValid (dbt->db, NULL, orp, FALSE, &dbxerr);
      if (dbxerr != NULL) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_IllegalDbXref, dbxerr);
        dbxerr = MemFree (dbxerr);
      }
      */
    }
  }

  if (GetAppProperty ("InternalNcbiSequin") == NULL) return;

  for (db = orp->db; db != NULL; db = db->next) {
    dbt = (DbtagPtr) db->data.ptrvalue;
    if (dbt != NULL) {
      if (StringICmp (dbt->db, "taxon") == 0)
        return;
    }
  }
  if (! vsp->seqSubmitParent) { /* suppress when validator run from tbl2asn */
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_NoTaxonID, "BioSource is missing taxon ID");
  }
}

static Boolean IsXr (ValNodePtr sdp)

{
  BioseqPtr      bsp;
  ObjValNodePtr  ovp;
  SeqIdPtr       sip;
  TextSeqIdPtr   tsip;

  if (sdp->extended == 0) return FALSE;
  ovp = (ObjValNodePtr) sdp;
  if (ovp->idx.parenttype != OBJ_BIOSEQ) return FALSE;
  bsp = (BioseqPtr) ovp->idx.parentptr;
  if (bsp == NULL) return FALSE;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice != SEQID_OTHER) continue;
    tsip = (TextSeqIdPtr) sip->data.ptrvalue;
    if (tsip == NULL) continue;
    if (StringNICmp (tsip->accession, "XR_", 3) == 0) return TRUE;
  }
  return FALSE;
}

static Boolean IsSynthetic (BioseqPtr bsp)

{
  BioSourcePtr       biop;
  SeqMgrDescContext  dcontext;
  OrgNamePtr         onp;
  OrgRefPtr          orp;
  SeqDescrPtr        sdp;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp == NULL) return FALSE;
  biop = (BioSourcePtr) sdp->data.ptrvalue;
  if (biop == NULL) return FALSE;
  if (biop->origin == 5) return TRUE;
  orp = biop->org;
  if (orp == NULL) return FALSE;
  onp = orp->orgname;
  if (onp == NULL) return FALSE;
  if (StringICmp (onp->div, "SYN") == 0) return TRUE;
  return FALSE;
}

static Boolean IsMicroRNA (BioseqPtr bsp)

{
  SeqMgrFeatContext  fcontext;
  RnaRefPtr          rrp;
  SeqFeatPtr         sfp;
  CharPtr            str;

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_otherRNA, &fcontext);
  while (sfp != NULL) {
    if (sfp->data.choice == SEQFEAT_RNA) {
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp != NULL && rrp->ext.choice == 1) {
        str = (CharPtr) rrp->ext.value.ptrvalue;
        if (StringStr (str, "microRNA") != NULL) return TRUE;
      }
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, FEATDEF_otherRNA, &fcontext);
  }
  return FALSE;
}

static Boolean IsOtherDNA (BioseqPtr bsp)

{
  SeqMgrDescContext  dcontext;
  MolInfoPtr         mip;
  SeqDescrPtr        sdp;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
  if (sdp == NULL) return FALSE;
  mip = (MolInfoPtr) sdp->data.ptrvalue;
  if (mip == NULL) return FALSE;
  if (mip->biomol == 255) return TRUE;
  return FALSE;
}

static Boolean StringHasPMID (CharPtr str)

{
  Char     ch;
  Int2     numdigits = 0;
  CharPtr  ptr;

  if (StringHasNoText (str)) return FALSE;

  ptr = StringStr (str, "(PMID ");
  if (ptr == NULL) return FALSE;

  ptr += 6;
  ch = *ptr;
  while (ch != '\0') {
    if (ch == ')') {
      if (numdigits > 0) return TRUE;
      return FALSE;
    } else if (IS_DIGIT (ch)) {
      numdigits++;
    }
    ptr++;
    ch = *ptr;
  }

  return FALSE;
}


static Boolean HasStructuredCommentPrefix (UserObjectPtr uop)
{
  UserFieldPtr ufp;

  if (uop == NULL) {
    return FALSE;
  }
  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    if (ufp->label != NULL && StringCmp (ufp->label->str, "StructuredCommentPrefix") == 0) {
      return TRUE;
    }
  }
  return FALSE;
}


static Boolean ValidateSeqDescrCommon (ValNodePtr sdp, BioseqValidStrPtr bvsp, ValidStructPtr vsp, Uint4 descitemid)
{
  ValNodePtr      vnp, vnp2;
  OrgRefPtr       this_org = NULL, that_org = NULL;
  int             tmpval;
  Char            buf1[20], buf2[20], ch;
  EMBLBlockPtr    ebp;
  GBBlockPtr      gbp;
  ValNodePtr      keywords = NULL;
  PubdescPtr      pdp;
  MolInfoPtr      mip;
  ObjectIdPtr     oip;
  Uint2           olditemtype = 0;
  Uint4           olditemid = 0;
  BioSourcePtr    biop;
  GatherContextPtr gcp = NULL;
  CharPtr         str, ptr;
  SeqFeatPtr      sfp;
  Boolean         tpa_exp;
  Boolean         tpa_inf;
  UserObjectPtr   uop;
  BioseqPtr       bsp;
  DatePtr         dp;
  size_t          len;
  SeqMgrFeatContext  fcontext;
  Int2            baddate;
  static char    *badmod = "Inconsistent GIBB-mod [%d] and [%d]";

  vsp->sfp = NULL;
  vnp = sdp;
  vsp->descr = vnp;

  if (descitemid > 0) {
    gcp = vsp->gcp;
    if (gcp != NULL) {
      olditemid = gcp->itemID;
      olditemtype = gcp->thistype;
      gcp->itemID = descitemid;
      gcp->thistype = OBJ_SEQDESC;
    }
  }

  switch (vnp->choice) {
  case Seq_descr_mol_type:
    tmpval = (int) (vnp->data.intvalue);
    switch (tmpval) {
    case 8:                    /* peptide */
      if (!bvsp->is_aa)
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_InvalidForType, "Nucleic acid with GIBB-mol = peptide");
      break;
    case 0:                    /* unknown */
    case 255:                  /* other */
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_InvalidForType, "GIBB-mol unknown or other used");
      break;
    default:                   /* the rest are nucleic acid */
      if (bvsp->is_aa) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_InvalidForType, "GIBB-mol [%d] used on protein", tmpval);
      } else {
        if (bvsp->last_na_mol) {
          if (bvsp->last_na_mol != (int) vnp->data.intvalue) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Inconsistent GIBB-mol [%d] and [%d]", bvsp->last_na_mol, tmpval);
          }
        } else
          bvsp->last_na_mol = tmpval;
      }
      break;
    }
    break;
  case Seq_descr_modif:
    for (vnp2 = (ValNodePtr) (vnp->data.ptrvalue); vnp2 != NULL; vnp2 = vnp2->next) {
      tmpval = (int) (vnp2->data.intvalue);
      switch (tmpval) {
      case 0:                  /* dna */
      case 1:                  /* rna */
        if (bvsp->is_aa) {      /* only temporarily on 0 */
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_InvalidForType, "Nucleic acid GIBB-mod [%d] on protein", tmpval);
        } else if (bvsp->last_na_mod) {
          if (tmpval != bvsp->last_na_mod) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, badmod, bvsp->last_na_mod, tmpval);
          }
        } else
          bvsp->last_na_mod = tmpval;
        break;
      case 4:                  /* mitochondria */
      case 5:                  /* chloroplast */
      case 6:                  /* kinetoplast */
      case 7:                  /* cyanelle */
      case 18:                 /* macronuclear */
        if (bvsp->last_organelle) {
          if (tmpval != bvsp->last_organelle) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, badmod, bvsp->last_organelle, tmpval);
          }
        } else
          bvsp->last_organelle = tmpval;
        break;
      case 10:                 /* partial */
      case 11:                 /* complete */
        if (bvsp->last_partialness) {
          if (tmpval != bvsp->last_partialness) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, badmod, bvsp->last_partialness, tmpval);
          }
        } else
          bvsp->last_partialness = tmpval;
        if ((bvsp->last_left_right) && (tmpval == 11)) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, badmod, bvsp->last_left_right, tmpval);
        }
        break;
      case 16:                 /* no left */
      case 17:                 /* no right */
        if (bvsp->last_partialness == 11) {     /* complete */
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, badmod, bvsp->last_partialness, tmpval);
        }
        bvsp->last_left_right = tmpval;
        break;
      case 255:                /* other */
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Unknown, "GIBB-mod = other used");
        break;
      default:
        break;

      }
    }
    break;
  case Seq_descr_method:
    if (!bvsp->is_aa) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_InvalidForType, "Nucleic acid with protein sequence method");
    }
    break;
  /*
  case Seq_descr_comment:
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_MissingText, "Comment descriptor needs text");
    }
    if (SerialNumberInString (str)) {
      ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_SerialInComment,
                "Comment may refer to reference by serial number - attach reference specific comments to the reference REMARK instead.");
    }
    if (StringLooksLikeFakeStructuredComment (str)) {
      ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_FakeStructuredComment,
                "Comment may be formatted to look like a structured comment.");
    }
    for (vnp2 = vnp->next; vnp2 != NULL; vnp2 = vnp2->next) {
      if (vnp2->choice == Seq_descr_comment) {
        ptr = (CharPtr) vnp2->data.ptrvalue;
        if (StringDoesHaveText (ptr) && StringICmp (str, ptr) == 0) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MultipleComments, "Undesired multiple comment descriptors, identical text");
        }
      }
    }
    break;
  */
  case Seq_descr_genbank:
    if (bvsp->last_gb != NULL)
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Multiple GenBank blocks");
    else
      bvsp->last_gb = vnp;
    if (vnp != NULL) {
      gbp = (GBBlockPtr) vnp->data.ptrvalue;
      if (gbp != NULL) {
        keywords = gbp->keywords;
      }
    }
    break;
  case Seq_descr_embl:
    if (bvsp->last_embl != NULL)
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Multiple EMBL blocks");
    else
      bvsp->last_embl = vnp;
    if (vnp != NULL) {
      ebp = (EMBLBlockPtr) vnp->data.ptrvalue;
      if (ebp != NULL) {
        keywords = ebp->keywords;
      }
    }
    break;
  case Seq_descr_pir:
    if (bvsp->last_pir != NULL)
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Multiple PIR blocks");
    else
      bvsp->last_pir = vnp;
    break;
  case Seq_descr_sp:
    if (bvsp->last_sp != NULL)
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Multiple SWISS-PROT blocks");
    else
      bvsp->last_sp = vnp;
    break;
  case Seq_descr_pdb:
    if (bvsp->last_pdb != NULL)
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Multiple PDB blocks");
    else
      bvsp->last_pdb = vnp;
    break;
  case Seq_descr_prf:
    if (bvsp->last_prf != NULL)
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Multiple PRF blocks");
    else
      bvsp->last_prf = vnp;
    break;
  case Seq_descr_create_date:
    dp = (DatePtr) vnp->data.ptrvalue;
    if (DateIsBad (dp, TRUE, &baddate)) {
      PrintBadDateError (vsp, baddate, SEV_ERROR, ERR_GENERIC_BadDate, "Create date has error");
    }
    if (bvsp->last_create != NULL) {
      tmpval = (int) DateMatch ((DatePtr) vnp->data.ptrvalue, (DatePtr) (bvsp->last_create->data.ptrvalue), FALSE);
      if (tmpval) {
        DatePrint ((DatePtr) (vnp->data.ptrvalue), buf1);
        DatePrint ((DatePtr) (bvsp->last_create->data.ptrvalue), buf2);
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_Inconsistent, "Inconsistent create_dates [%s] and [%s]", buf1, buf2);
      }
    } else
      bvsp->last_create = vnp;
    if (bvsp->last_update != NULL) {
      tmpval = (int) DateMatch ((DatePtr) vnp->data.ptrvalue, (DatePtr) (bvsp->last_update->data.ptrvalue), FALSE);
      if (tmpval == 1) {
        DatePrint ((DatePtr) (vnp->data.ptrvalue), buf1);
        DatePrint ((DatePtr) (bvsp->last_update->data.ptrvalue), buf2);
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_Inconsistent, "Inconsistent create_date [%s] and update_date [%s]", buf1, buf2);
      }
    }
    break;
  case Seq_descr_update_date:
    dp = (DatePtr) vnp->data.ptrvalue;
    if (DateIsBad (dp, TRUE, &baddate)) {
      PrintBadDateError (vsp, baddate, SEV_ERROR, ERR_GENERIC_BadDate, "Update date has error");
    }
    if (bvsp->last_create != NULL) {
      tmpval = (int) DateMatch ((DatePtr) bvsp->last_create->data.ptrvalue, (DatePtr) (vnp->data.ptrvalue), FALSE);
      if (tmpval == 1) {
        DatePrint ((DatePtr) (bvsp->last_create->data.ptrvalue), buf1);
        DatePrint ((DatePtr) (vnp->data.ptrvalue), buf2);
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_Inconsistent, "Inconsistent create_date [%s] and update_date [%s]", buf1, buf2);
      }
    }
    if (bvsp->last_update == NULL)
      bvsp->last_update = vnp;
    break;
  case Seq_descr_source:
    biop = (BioSourcePtr) vnp->data.ptrvalue;
    bsp = bvsp->bsp;
    if (biop != NULL && biop->is_focus && bsp != NULL) {
      if (ISA_aa (bsp->mol) || bsp->repr == Seq_repr_seg || SeqMgrGetParentOfPart (bsp, NULL) != NULL) {
        /* skip proteins, segmented bioseqs, or segmented parts */
      } else {
        sfp = SeqMgrGetNextFeature (bvsp->bsp, NULL, SEQFEAT_BIOSRC, 0, &fcontext);
        if (sfp == NULL) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_UnnecessaryBioSourceFocus, "BioSource descriptor has focus, but no BioSource feature");
        }
      }
    }
    if (biop != NULL && biop->origin == 5) {
      bsp = bvsp->bsp;
      if (! IsOtherDNA (bsp) && !ISA_aa (bsp->mol)) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_InvalidForType, "Molinfo-biomol other should be used if Biosource-location is synthetic");
      }
    }
    /* ValidateBioSource (vsp, gcp, biop, NULL, vnp); */
    if (biop != NULL) {
      this_org = biop->org;
    }
    /* fall into Seq_descr_org */
  case Seq_descr_org:
    if (this_org == NULL)
      this_org = (OrgRefPtr) (vnp->data.ptrvalue);
    if (bvsp->last_org != NULL) {
      if ((this_org->taxname != NULL) && (bvsp->last_org->taxname != NULL)) {
        if (StringCmp (this_org->taxname, bvsp->last_org->taxname)) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Inconsistent taxnames [%s] and [%s]", this_org->taxname, bvsp->last_org->taxname);
        }
      }
    } else
      bvsp->last_org = this_org;

    for (vnp2 = vnp->next; vnp2 != NULL; vnp2 = vnp2->next) {
      if (vnp2->choice == Seq_descr_source || vnp2->choice == Seq_descr_org) {
        that_org = NULL;
        if (vnp2->choice == Seq_descr_source) {
          that_org = ((BioSourcePtr) (vnp2->data.ptrvalue))->org;
        }
        if (that_org == NULL) {
          that_org = (OrgRefPtr) (vnp2->data.ptrvalue);
        }
        if (that_org != NULL) {
          if ((this_org->taxname != NULL) && (that_org->taxname != NULL) && StringCmp (this_org->taxname, that_org->taxname) == 0) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_MultipleBioSources, "Undesired multiple source descriptors");
          }
        }
      }
    }
    break;
  case Seq_descr_title:
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_MissingText, "Title descriptor needs text");
    }
    for (vnp2 = vnp->next; vnp2 != NULL; vnp2 = vnp2->next) {
      if (vnp2->choice == Seq_descr_title) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_MultipleTitles, "Undesired multiple title descriptors");
        break;
      }
    }
    len = StringLen (str);
    if (len > 4) {
      ch = str [len - 1];
      while (ch == ' ' && len > 4) {
        len--;
        ch = str [len - 1];
      }
      if (ch == '.' && len > 4) {
        len--;
        ch = str [len - 1];
      }
      if (ch == '.' || ch == ',' || ch == ';' || ch == ':') {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadPunctuation, "Title descriptor ends in bad punctuation");
      }
    }
    if (StringHasPMID (str)) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_TitleHasPMID, "Title descriptor has internal PMID");
    }
    break;
  case Seq_descr_name:
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_MissingText, "Name descriptor needs text");
    }
    for (vnp2 = vnp->next; vnp2 != NULL; vnp2 = vnp2->next) {
      if (vnp2->choice == Seq_descr_name) {
        ptr = (CharPtr) vnp2->data.ptrvalue;
        if (StringDoesHaveText (ptr) && StringICmp (str, ptr) == 0) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MultipleNames, "Undesired multiple name descriptors, identical text");
        } else {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MultipleNames, "Undesired multiple name descriptors, different text");
        }
      }
    }
    break;
  case Seq_descr_region:
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_MissingText, "Region descriptor needs text");
    }
    break;
  case Seq_descr_user:
    uop = (UserObjectPtr) vnp->data.ptrvalue;
    if (uop != NULL) {
      oip = uop->type;
      if (oip != NULL) {
        if (StringCmp (oip->str, "StructuredComment") == 0) {
          if (uop->data == NULL) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_UserObjectProblem, "Structured Comment user object descriptor is empty");
          }
          if (!HasStructuredCommentPrefix (uop)) {
            ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_StructuredCommentPrefixOrSuffixMissing, "Structured Comment lacks prefix");
          }
        }
      }
    }
    break;
  case Seq_descr_pub:
    bvsp->got_a_pub = TRUE;
    pdp = (PubdescPtr) vnp->data.ptrvalue;
    /*
       ValidatePubdesc (vsp, pdp);
     */
    break;
  case Seq_descr_molinfo:
    mip = (MolInfoPtr) vnp->data.ptrvalue;
    if (mip != NULL) {
      switch (mip->biomol) {
      case MOLECULE_TYPE_PEPTIDE:      /* peptide */
        if (!bvsp->is_aa) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_InvalidForType, "Nucleic acid with Molinfo-biomol = peptide");
        }
        break;
      case MOLECULE_TYPE_OTHER_GENETIC_MATERIAL:
        if (! bvsp->is_artificial) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_InvalidForType, "Molinfo-biomol = other genetic");
        }
        break;
      case 0:                  /* unknown */
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_InvalidForType, "Molinfo-biomol unknown used");
        break;
      case 255:                /* other */
        if (! IsXr (vnp)) {
          bsp = bvsp->bsp;
          if (! IsSynthetic (bsp)) {
            if (! IsMicroRNA (bsp)) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_InvalidForType, "Molinfo-biomol other used");
            }
          }
        }
        break;
      default:                 /* the rest are nucleic acid */
        if (bvsp->is_aa) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_InvalidForType, "Molinfo-biomol [%d] used on protein", (int) mip->biomol);
        } else {
          if (bvsp->last_biomol) {
            if (bvsp->last_biomol != (int) mip->biomol) {
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Inconsistent Molinfo-biomol [%d] and [%d]", bvsp->last_biomol, (int) mip->biomol);
            }
          } else {
            bvsp->last_biomol = (int) mip->biomol;
          }
        }
        break;
      }

      if (bvsp->is_syn_constr) {
        if (mip->biomol != MOLECULE_TYPE_OTHER_GENETIC_MATERIAL && !bvsp->is_aa) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_InvalidForType, "synthetic construct should have other-genetic");
        }
        if (! bvsp->is_artificial) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_InvalidForType, "synthetic construct should have artificial origin");
        }
      } else if (bvsp->is_artificial) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_InvalidForType, "artificial origin should have other-genetic and synthetic construct");
      }
      if (bvsp->is_artificial) {
        if (mip->biomol != MOLECULE_TYPE_OTHER_GENETIC_MATERIAL && mip->biomol != MOLECULE_TYPE_PEPTIDE) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_InvalidForType, "artificial origin should have other-genetic");
        }
      }
      if (!bvsp->is_aa) {
        switch (mip->tech) {
        case MI_TECH_concept_trans:
        case MI_TECH_seq_pept:
        case MI_TECH_both:
        case MI_TECH_seq_pept_overlap:
        case MI_TECH_seq_pept_homol:
        case MI_TECH_concept_trans_a:
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_InvalidForType, "Nucleic acid with protein sequence method");
          break;
        default:
          break;
        }
      } else {
        switch (mip->tech) {
        case MI_TECH_est:
        case MI_TECH_sts:
        case MI_TECH_genemap:
        case MI_TECH_physmap:
        case MI_TECH_htgs_1:
        case MI_TECH_htgs_2:
        case MI_TECH_htgs_3:
        case MI_TECH_fli_cdna:
        case MI_TECH_htgs_0:
        case MI_TECH_htc:
        case MI_TECH_wgs:
        case MI_TECH_barcode:
        case MI_TECH_composite_wgs_htgs:
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_InvalidForType, "Protein with nucleic acid sequence method");
          break;
        default:
          break;
        }
      }
      if (bvsp->last_tech) {
        if (bvsp->last_tech != (int) mip->tech) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Inconsistent Molinfo-tech [%d] and [%d]", bvsp->last_tech, (int) mip->tech);
        }
      } else {
        bvsp->last_tech = (int) mip->tech;
      }
      if (bvsp->last_completeness) {
        if (bvsp->last_completeness != (int) mip->completeness) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Inconsistent Molinfo-completeness [%d] and [%d]",
                    bvsp->last_completeness, (int) mip->completeness);
        }
      } else {
        bvsp->last_completeness = (int) mip->completeness;
      }
    }
    break;
  default:
    break;
  }

  if (keywords != NULL) {
    tpa_exp = FALSE;
    tpa_inf = FALSE;
    for (vnp = keywords; vnp != NULL; vnp = vnp->next) {
      if (StringICmp ((CharPtr) vnp->data.ptrvalue, "TPA:experimental") == 0) {
        tpa_exp = TRUE;
      } else if (StringICmp ((CharPtr) vnp->data.ptrvalue, "TPA:inferential") == 0) {
        tpa_inf = TRUE;
      }
    }
    if (tpa_exp && tpa_inf) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "TPA:experimental and TPA:inferential should not both be in the same set of keywords");
    }
  }

  if (gcp != NULL) {
    gcp->itemID = olditemid;
    gcp->thistype = olditemtype;
  }

  return TRUE;
}

static Boolean LIBCALLBACK ValidateSeqDescrIndexed (ValNodePtr sdp, SeqMgrDescContextPtr context)
{
  ValidStructPtr  vsp;
  BioseqValidStrPtr bvsp;

  bvsp = (BioseqValidStrPtr) context->userdata;
  vsp = bvsp->vsp;

  return ValidateSeqDescrCommon (sdp, bvsp, vsp, context->itemID);
}

static void ValidateSeqDescrContext (GatherContextPtr gcp)
{
  ValidStructPtr  vsp;
  BioseqValidStrPtr bvsp;
  ValNodePtr      sdp;

  bvsp = (BioseqValidStrPtr) (gcp->userdata);
  vsp = bvsp->vsp;
  sdp = (ValNodePtr) (gcp->thisitem);

  ValidateSeqDescrCommon (sdp, bvsp, vsp, 0);
}

/*****************************************************************************
*
*   ValidateBioseqContextGather(gcp)
*      Gather callback for validating context on a Bioseq
*
*****************************************************************************/
static Boolean DifferentDbxrefs (ValNodePtr dbxref1, ValNodePtr dbxref2)
{
  DbtagPtr        dbt1, dbt2;
  ObjectIdPtr     oip1, oip2;

  if (dbxref1 == NULL || dbxref2 == NULL)
    return FALSE;
  dbt1 = (DbtagPtr) dbxref1->data.ptrvalue;
  dbt2 = (DbtagPtr) dbxref2->data.ptrvalue;
  if (dbt1 == NULL || dbt2 == NULL)
    return FALSE;
  if (StringICmp (dbt1->db, dbt2->db) != 0)
    return TRUE;
  oip1 = dbt1->tag;
  oip2 = dbt2->tag;
  if (oip1 == NULL || oip2 == NULL)
    return FALSE;
  if (oip1->str == NULL && oip2->str == NULL) {
    if (oip1->id != oip2->id)
      return TRUE;
  } else {
    if (StringICmp (oip1->str, oip2->str) != 0)
      return TRUE;
  }
  return FALSE;
}

static Boolean FlybaseDbxrefs (ValNodePtr vnp)

{
  DbtagPtr  dbt;

  while (vnp != NULL) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt != NULL) {
      if (StringCmp (dbt->db, "FLYBASE") == 0 || StringCmp (dbt->db, "FlyBase") == 0) {
        return TRUE;
      }
    }
    vnp = vnp->next;
  }
  return FALSE;
}

static Boolean GPSorNTorNCorNGorNW (SeqEntryPtr sep, SeqLocPtr location)
{
  BioseqPtr       bsp;
  BioseqSetPtr    bssp;
  SeqIdPtr        sip;
  TextSeqIdPtr    tsip;

  if (sep != NULL && IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && bssp->_class == BioseqseqSet_class_gen_prod_set) {
      return TRUE;
    }
  }
  bsp = BioseqFindFromSeqLoc (location);
  if (bsp != NULL) {
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_OTHER) {
        tsip = (TextSeqIdPtr) sip->data.ptrvalue;
        if (tsip != NULL && tsip->accession != NULL) {
          if (StringNICmp (tsip->accession, "NT_", 3) == 0) {
            return TRUE;
          } else if (StringNICmp (tsip->accession, "NC_", 3) == 0) {
            return TRUE;
          } else if (StringNICmp (tsip->accession, "NG_", 3) == 0) {
            return TRUE;
          } else if (StringNICmp (tsip->accession, "NW_", 3) == 0) {
            return TRUE;
          }
        }
      }
    }
  }
  return FALSE;
}

static Boolean IsGenBankAccn (SeqEntryPtr sep, SeqLocPtr location)
{
  BioseqPtr  bsp;
  SeqIdPtr   sip;

  bsp = BioseqFindFromSeqLoc (location);
  if (bsp != NULL) {
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_GENBANK) return TRUE;
    }
  }
  return FALSE;
}

static Boolean IsEMBLAccn (SeqEntryPtr sep, SeqLocPtr location)
{
  BioseqPtr  bsp;
  SeqIdPtr   sip;

  bsp = BioseqFindFromSeqLoc (location);
  if (bsp != NULL) {
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_EMBL) return TRUE;
    }
  }
  return FALSE;
}

static Boolean IsGeneralAccn (SeqEntryPtr sep, SeqLocPtr location)
{
  BioseqPtr  bsp;
  DbtagPtr   dbt;
  SeqIdPtr   sip;

  bsp = BioseqFindFromSeqLoc (location);
  if (bsp != NULL) {
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      if (sip->choice != SEQID_GENERAL) continue;
      dbt = (DbtagPtr) sip->data.ptrvalue;
      if (dbt == NULL) continue;
      if (IsSkippableDbtag(dbt)) continue;
      return TRUE;
    }
  }
  return FALSE;
}

static Boolean NGorNT (SeqEntryPtr sep, SeqLocPtr location, BoolPtr is_nc)
{
  BioseqPtr       bsp;
  SeqIdPtr        sip;
  TextSeqIdPtr    tsip;

  if (is_nc != NULL) {
    *is_nc = FALSE;
  }
  bsp = BioseqFindFromSeqLoc (location);
  if (bsp != NULL) {
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_OTHER) {
        tsip = (TextSeqIdPtr) sip->data.ptrvalue;
        if (tsip != NULL && tsip->accession != NULL) {
          if (StringNICmp (tsip->accession, "NT_", 3) == 0) {
            return TRUE;
          } else if (StringNICmp (tsip->accession, "NG_", 3) == 0) {
            return TRUE;
          } else if (StringNICmp (tsip->accession, "NW_", 3) == 0) {
            return TRUE;
          } else if (StringNICmp (tsip->accession, "NC_", 3) == 0 && is_nc != NULL) {
            *is_nc = TRUE;
          }
        }
      }
    }
  }
  return FALSE;
}

static Boolean GPSorRefSeq (SeqEntryPtr sep, SeqLocPtr location)
{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  SeqIdPtr      sip;

  if (sep != NULL && IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && bssp->_class == BioseqseqSet_class_gen_prod_set) {
      return TRUE;
    }
  }
  bsp = BioseqFindFromSeqLoc (location);
  if (bsp != NULL) {
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_OTHER) {
        return TRUE;
      }
    }
  }
  return FALSE;
}

static Boolean IsNCorNT (SeqEntryPtr sep, SeqLocPtr location)
{
  BioseqPtr       bsp;
  SeqIdPtr        sip;
  TextSeqIdPtr    tsip;

  bsp = BioseqFindFromSeqLoc (location);
  if (bsp != NULL) {
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_OTHER) {
        tsip = (TextSeqIdPtr) sip->data.ptrvalue;
        if (tsip != NULL && tsip->accession != NULL) {
          if (StringNICmp (tsip->accession, "NC_", 3) == 0) {
            return TRUE;
          } else if (StringNICmp (tsip->accession, "NT_", 3) == 0) {
            return TRUE;
          }
        }
      }
    }
  }
  return FALSE;
}

static Boolean IsNCorNTorNW (SeqEntryPtr sep, SeqLocPtr location)
{
  BioseqPtr       bsp;
  SeqIdPtr        sip;
  TextSeqIdPtr    tsip;

  bsp = BioseqFindFromSeqLoc (location);
  if (bsp != NULL) {
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_OTHER) {
        tsip = (TextSeqIdPtr) sip->data.ptrvalue;
        if (tsip != NULL && tsip->accession != NULL) {
          if (StringNICmp (tsip->accession, "NC_", 3) == 0) {
            return TRUE;
          } else if (StringNICmp (tsip->accession, "NT_", 3) == 0) {
            return TRUE;
          } else if (StringNICmp (tsip->accession, "NW_", 3) == 0) {
            return TRUE;
          }
        }
      }
    }
  }
  return FALSE;
}

static Boolean NotPeptideException (SeqFeatPtr sfp, SeqFeatPtr last)
{
  if (sfp != NULL && sfp->excpt) {
    if (StringISearch (sfp->except_text, "alternative processing") != NULL)
      return FALSE;
  }
  if (last != NULL && last->excpt) {
    if (StringISearch (last->except_text, "alternative processing") != NULL)
      return FALSE;
  }
  return TRUE;
}

static Boolean DescsSame (AnnotDescrPtr adp1, AnnotDescrPtr adp2)

{
  if (adp1 == NULL || adp2 == NULL) return TRUE;
  if (adp1->choice != adp2->choice) return FALSE;
  if (adp1->choice == Annot_descr_name || adp1->choice == Annot_descr_title) {
    if (StringICmp ((CharPtr) adp1->data.ptrvalue, (CharPtr) adp2->data.ptrvalue) == 0) return TRUE;
  }
  return FALSE;
}

typedef struct gmcdata {
  SeqFeatPtr  gene;
  SeqFeatPtr  feat;
} GmcData, PNTR GmcDataPtr;

static int LIBCALLBACK SortGmcByGenePtr (
  VoidPtr vp1,
  VoidPtr vp2
)

{
  GmcDataPtr gdp1, gdp2;

  if (vp1 == NULL || vp2 == NULL) return 0;
  gdp1 = (GmcDataPtr) vp1;
  gdp2 = (GmcDataPtr) vp2;
  if (gdp1 == NULL || gdp2 == NULL) return 0;

  if (gdp1->gene > gdp2->gene) return -1;
  if (gdp1->gene < gdp2->gene) return 1;

  if (gdp1->feat > gdp2->feat) return -1;
  if (gdp1->feat < gdp2->feat) return 1;

  return 0;
}

static void ValidateLocusTagGeneral (ValidStructPtr vsp, BioseqPtr bsp)

{
  DbtagPtr           dbt;
  SeqMgrFeatContext  fcontext;
  GatherContextPtr   gcp;
  SeqFeatPtr         gene;
  GeneRefPtr         grp;
  ObjectIdPtr        oip;
  Uint2              olditemtype = 0;
  Uint4              olditemid = 0;
  BioseqPtr          prod;
  CharPtr            ptr;
  SeqFeatPtr         sfp;
  SeqIdPtr           sip;
  Char               tmp [64];

  if (vsp == NULL || bsp == NULL) return;
  if (! ISA_na (bsp->mol)) return;

  gcp = vsp->gcp;
  if (gcp != NULL) {
    olditemid = gcp->itemID;
    olditemtype = gcp->thistype;
  }

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  while (sfp != NULL) {
    if (sfp->idx.subtype == FEATDEF_CDS || sfp->idx.subtype == FEATDEF_mRNA) {
      grp = SeqMgrGetGeneXref (sfp);
      if (! SeqMgrGeneIsSuppressed (grp)) {
        if (grp == NULL) {
          gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
          if (gene != NULL) {
            grp = (GeneRefPtr) sfp->data.value.ptrvalue;
          }
        }
        if (grp != NULL && StringDoesHaveText (grp->locus_tag)) {
          prod = BioseqFindFromSeqLoc (sfp->product);
          if (prod != NULL) {
            for (sip = prod->id; sip != NULL; sip = sip->next) {
              if (sip->choice != SEQID_GENERAL) continue;
              dbt = (DbtagPtr) sip->data.ptrvalue;
              if (dbt == NULL) continue;
              if (IsSkippableDbtag(dbt)) continue;
              oip = dbt->tag;
              if (oip == NULL) continue;
              if (StringHasNoText (oip->str)) continue;
              StringNCpy_0 (tmp, oip->str, sizeof (tmp));
              ptr = StringChr (tmp, '-');
              if (ptr != NULL) {
                *ptr = '\0';
              }
              if (StringICmp (grp->locus_tag, tmp) != 0) {
                if (gcp != NULL) {
                  gcp->itemID = sfp->idx.itemID;
                  gcp->thistype = OBJ_SEQFEAT;
                }
                vsp->descr = NULL;
                vsp->sfp = sfp;
                ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_LocusTagProductMismatch, "Gene locus_tag does not match general ID of product");
              }
            }
          }
        }
      }
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
  }

  if (gcp != NULL) {
    gcp->itemID = olditemid;
    gcp->thistype = olditemtype;
  }
}

static Boolean ReplaceQualsDiffer (GBQualPtr sfpqual, GBQualPtr lastqual)

{
  if (sfpqual == NULL || lastqual == NULL) return FALSE;

  while (sfpqual != NULL && StringICmp (sfpqual->qual, "replace") != 0) {
    sfpqual = sfpqual->next;
  }
  while (lastqual != NULL && StringICmp (lastqual->qual, "replace") != 0) {
    lastqual = lastqual->next;
  }
  if (sfpqual == NULL || lastqual == NULL) return FALSE;

  if (StringICmp (sfpqual->val, lastqual->val) != 0) return TRUE;

  return FALSE;
}

static Boolean GBQualsDiffer (GBQualPtr sfpqual, GBQualPtr lastqual)

{
  if (sfpqual == NULL || lastqual == NULL) return FALSE;

  /* depends upon sorted order of gbquals imposed by BasicSeqEntryCleanup */

  while (sfpqual != NULL && lastqual != NULL) {
    if (StringICmp (sfpqual->qual, lastqual->qual) != 0) return TRUE;
    if (StringICmp (sfpqual->val, lastqual->val) != 0) return TRUE;
    sfpqual = sfpqual->next;
    lastqual = lastqual->next;
  }

  if (sfpqual != NULL || lastqual != NULL) return TRUE;

  return FALSE;
}

static CharPtr MakePubLabelString (PubdescPtr pdp)

{
  Char        buf [521];
  CitGenPtr   cgp;
  ValNodePtr  vnp;

  if (pdp == NULL) return NULL;

  vnp = pdp->pub;

  /* skip over just serial number */

  if (vnp != NULL && vnp->choice == PUB_Gen && vnp->next != NULL) {
    cgp = (CitGenPtr) vnp->data.ptrvalue;
    if (cgp != NULL) {
      if (StringNICmp ("BackBone id_pub", cgp->cit, 15) != 0) {
        if (cgp->cit == NULL && cgp->journal == NULL && cgp->date == NULL && cgp->serial_number) {
          vnp = vnp->next;
        }
      }
    }
  }

  if (PubLabelUnique (vnp, buf, sizeof (buf) - 1, OM_LABEL_CONTENT, TRUE) > 0) {
    return StringSaveNoNull (buf);
  }

  return NULL;
}

static CharPtr ValGetAuthorsPlusConsortium (
  AuthListPtr alp
)

{
  CharPtr  consortium;
  CharPtr  str;
  CharPtr  tmp;

  consortium = NULL;
  str = GetAuthorsString (GENBANK_FMT, alp, &consortium, NULL, NULL);
  if (str == NULL) return consortium;
  if (consortium == NULL) return str;
  tmp = MemNew (StringLen (str) + StringLen (consortium) + 5);
  if (tmp == NULL) return NULL;
  StringCpy (tmp, str);
  StringCat (tmp, "; ");
  StringCat (tmp, consortium);
  MemFree (str);
  MemFree (consortium);
  return tmp;
}

static Boolean IsIdenticalPublication (PubdescPtr pdp1, PubdescPtr pdp2)

{
  AuthListPtr  alp1, alp2;
  Boolean      rsult = TRUE;
  CharPtr      str1, str2;

  if (pdp1 == NULL || pdp2 == NULL) return FALSE;

  str1 = MakePubLabelString (pdp1);
  str2 = MakePubLabelString (pdp2);
  if (StringDoesHaveText (str1) && StringDoesHaveText (str2)) {
    if (StringICmp (str1, str2) != 0) {
       rsult = FALSE;
    }
  }
  MemFree (str1);
  MemFree (str2);
  if (! rsult) return rsult;

  alp1 = GetAuthListPtr (pdp1, NULL);
  alp2 = GetAuthListPtr (pdp2, NULL);
  if (alp1 != NULL && alp2 != NULL) {
    str1 = ValGetAuthorsPlusConsortium (alp1);
    str2 = ValGetAuthorsPlusConsortium (alp2);
    if (StringDoesHaveText (str1) && StringDoesHaveText (str2)) {
      if (StringICmp (str1, str2) != 0) {
         rsult = FALSE;
      }
    }
    MemFree (str1);
    MemFree (str2);
  }

  return rsult;
}

static Boolean IsIdenticalBioSource (BioSourcePtr biop1, BioSourcePtr biop2)

{
  DbtagPtr      dbt1, dbt2;
  ObjectIdPtr   oip1, oip2;
  OrgModPtr     omp1, omp2;
  OrgNamePtr    onp1, onp2;
  OrgRefPtr     orp1, orp2;
  SubSourcePtr  ssp1, ssp2;
  ValNodePtr    vnp1, vnp2;

  if (biop1 == NULL || biop2 == NULL) return FALSE;

  if (biop1->is_focus != biop2->is_focus) return FALSE;

  orp1 = biop1->org;
  orp2 = biop2->org;
  if (orp1 == NULL || orp2 == NULL) return FALSE;
  if (StringICmp (orp1->taxname, orp2->taxname) != 0) return FALSE;

  onp1 = orp1->orgname;
  onp2 = orp2->orgname;
  if (onp1 == NULL || onp2 == NULL) return FALSE;

  omp1 = onp1->mod;
  omp2 = onp2->mod;
  while (omp1 != NULL && omp2 != NULL) {
    if (omp1->subtype != omp2->subtype) return FALSE;
    if (StringICmp (omp1->subname, omp2->subname) != 0) return FALSE;
    omp1 = omp1->next;
    omp2 = omp2->next;
  }
  if (omp1 != NULL || omp2 != NULL) return FALSE;

  ssp1 = biop1->subtype;
  ssp2 = biop2->subtype;
  while (ssp1 != NULL && ssp2 != NULL) {
    if (ssp1->subtype != ssp2->subtype) return FALSE;
    if (StringICmp(ssp1->name, ssp2->name) != 0) return FALSE;
    ssp1 = ssp1->next;
    ssp2 = ssp2->next;
  }
  if (ssp1 != NULL || ssp2 != NULL) return FALSE;

  vnp1 = orp1->db;
  vnp2 = orp2->db;
  while (vnp1 != NULL && vnp2 != NULL) {
    dbt1 = (DbtagPtr) vnp1->data.ptrvalue;
    dbt2 = (DbtagPtr) vnp2->data.ptrvalue;

    if ((dbt1 != NULL) && (dbt2 != NULL)) {
      if (StringCmp (dbt1->db, dbt2->db) != 0) return FALSE;

      oip1 = dbt1->tag;
      oip2 = dbt2->tag;
      if ((oip1 != NULL) && (oip2 != NULL)) {
        if (oip1->str != NULL) {
          if (StringICmp(oip1->str, oip2->str) != 0) return FALSE;
        } else  {
          if (oip1->id != oip2->id) return FALSE;
        }
      }
      else if (oip1 != NULL)
        return FALSE;
      else if (oip2 != NULL)
        return FALSE;
    }
    else if (dbt1 != NULL)
      return FALSE;
    else if (dbt2 != NULL)
      return FALSE;

    vnp1 = vnp1->next;
    vnp2 = vnp2->next;
  }
  if (vnp1 != NULL || vnp2 != NULL) return FALSE;

  return TRUE;
}

typedef struct lpdata {
  Int2        count;
  SeqFeatPtr  cds;
  SeqFeatPtr  mrna;
  Char        firstid [64];
  Boolean     products_unique;
  Boolean     featid_matched;
} LpData, PNTR LpDataPtr;

static Boolean IdXrefsAreReciprocal (
  SeqFeatPtr cds,
  SeqFeatPtr mrna
)

{
  SeqFeatXrefPtr  xref;
  Boolean         match1 = FALSE, match2 = FALSE;
  SeqFeatPtr      matchsfp;

  if (cds == NULL || mrna == NULL) return FALSE;
  if (cds->id.choice != 3 || mrna->id.choice != 3) return FALSE;

  for (xref = cds->xref; xref != NULL; xref = xref->next) {
    if (xref->id.choice != 0) {
      matchsfp = SeqMgrGetFeatureByFeatID (cds->idx.entityID, NULL, NULL, xref, NULL);
      if (matchsfp == mrna) {
        match1 = TRUE;
      }
    }
  }

  for (xref = mrna->xref; xref != NULL; xref = xref->next) {
    if (xref->id.choice != 0) {
      matchsfp = SeqMgrGetFeatureByFeatID (mrna->idx.entityID, NULL, NULL, xref, NULL);
      if (matchsfp == cds) {
        match2 = TRUE;
      }
    }
  }

  if (match1 && match2) return TRUE;
  return FALSE;
}

static Int2 IdXrefsNotReciprocal (
  SeqFeatPtr cds,
  SeqFeatPtr mrna
)

{
  Int4            giu = 0, gip = 0;
  SeqFeatPtr      matchsfp;
  ObjectIdPtr     oip;
  SeqIdPtr        sip;
  CharPtr         tmp;
  UserFieldPtr    ufp;
  UserObjectPtr   uop;
  SeqFeatXrefPtr  xref;

  if (cds == NULL || mrna == NULL) return 0;
  if (cds->id.choice != 3 || mrna->id.choice != 3) return 0;

  for (xref = cds->xref; xref != NULL; xref = xref->next) {
    if (xref->id.choice != 0) {
      matchsfp = SeqMgrGetFeatureByFeatID (cds->idx.entityID, NULL, NULL, xref, NULL);
      if (matchsfp != NULL && matchsfp->idx.subtype == FEATDEF_mRNA && matchsfp != mrna) {
        return 1;
      }
    }
  }

  for (xref = mrna->xref; xref != NULL; xref = xref->next) {
    if (xref->id.choice != 0) {
      matchsfp = SeqMgrGetFeatureByFeatID (mrna->idx.entityID, NULL, NULL, xref, NULL);
      if (matchsfp != NULL && matchsfp->idx.subtype == FEATDEF_CDS && matchsfp != cds) {
        return 1;
      }
    }
  }

  if (cds->product == NULL) return 0;
  if (mrna->ext == NULL) return 0;
  uop = FindUopByTag (mrna->ext, "MrnaProteinLink");
  if (uop == NULL) return 0;
  sip = SeqLocId (cds->product);
  if (sip == NULL) return 0;
  if (sip->choice == SEQID_GI) {
    gip = (Int4) sip->data.intvalue;
  } else {
    gip = GetGIForSeqId (sip);
  }
  if (gip == 0) return 0;
  ufp = uop->data;
  if (ufp == NULL || ufp->choice != 1) return 0;
  oip = ufp->label;
  if (oip == NULL || StringICmp (oip->str, "protein seqID") != 0) return 0;
  tmp = (CharPtr) ufp->data.ptrvalue;
  if (StringHasNoText (tmp)) return 0;
  sip = MakeSeqID (tmp);
  if (sip == NULL) return 0;
  if (sip->choice == SEQID_GI) {
    giu = (Int4) sip->data.intvalue;
  } else {
    giu = GetGIForSeqId (sip);
  }
  SeqIdFree (sip);
  if (giu == 0) return 0;
  if (gip != giu) return 2;

  return 0;
}

static Boolean LIBCALLBACK FindSingleMrnaProc (
  SeqFeatPtr sfp,
  SeqMgrFeatContextPtr context
)

{
  Char        buf [64];
  SeqFeatPtr  cds;
  LpDataPtr   ldp;
  SeqIdPtr    sip;
  VvmDataPtr  vdp;

  ldp = (LpDataPtr) context->userdata;
  if (ldp == NULL) return TRUE;
  cds = ldp->cds;
  if (cds == NULL) return TRUE;

  if (sfp->product) {
    if (StringHasNoText (ldp->firstid)) {
      sip = SeqLocId (sfp->product);
      SeqIdWrite (sip, ldp->firstid, PRINTID_FASTA_LONG, sizeof (ldp->firstid) - 1);
    } else {
      sip = SeqLocId (sfp->product);
      SeqIdWrite (sip, buf, PRINTID_FASTA_LONG, sizeof (buf) - 1);
      if (StringCmp (ldp->firstid, buf) == 0) {
        ldp->products_unique = FALSE;
      }
    }
  }

  vdp = (VvmDataPtr) sfp->idx.scratch;
  if (vdp != NULL && vdp->accounted_for) return TRUE;

  (ldp->count)++;
  ldp->mrna = sfp;

  if (IdXrefsAreReciprocal (cds, sfp)) {
    ldp->featid_matched = TRUE;
  }

  return TRUE;
}


static Boolean MarkMrnasFromCDSXrefs (SeqFeatPtr cds, LpDataPtr lp)
{
  SeqFeatXrefPtr     xref;
  Boolean            has_xref = FALSE;
  BioseqExtraPtr     bspextra;
  SMFidItemPtr PNTR  array;
  Char               buf [32];
  CharPtr            featid = NULL;
  SMFeatItemPtr      feat;
  SMFidItemPtr       item;
  Int4               L;
  Int4               mid;
  Int4               num;
  ObjectIdPtr        oip;
  ObjMgrDataPtr      omdp;
  Int4               R;
  SeqFeatPtr         sfp;
  SeqMgrFeatContext  context;

  if (cds == NULL) {
    return FALSE;
  }

  omdp = ObjMgrGetData (cds->idx.entityID);
  if (omdp == NULL) {
    return FALSE;
  }

  bspextra = (BioseqExtraPtr) omdp->extradata;
  if (bspextra == NULL) return FALSE;
  array = bspextra->featsByFeatID;
  num = bspextra->numfids;
  if (array == NULL || num < 1) return FALSE;

  context.userdata = lp;

  for (xref = cds->xref; xref != NULL; xref = xref->next) {
    if (xref != NULL && xref->id.choice == 3) {
      featid = NULL;
      oip = (ObjectIdPtr) xref->id.value.ptrvalue;
      if (oip != NULL) {
        if (StringDoesHaveText (oip->str)) {
          featid = oip->str;
        } else {
          sprintf (buf, "%ld", (long) oip->id);
          featid = buf;
        }
      }
      if (StringHasNoText (featid)) continue;

      L = 0;
      R = num - 1;
      while (L < R) {
        mid = (L + R) / 2;
        item = array [mid];
        if (item != NULL && StringICmp (item->fid, featid) < 0) {
          L = mid + 1;
        } else {
          R = mid;
        }
      }

      while (R > 0 && StringICmp (array[R - 1]->fid, featid) == 0) {
        R--;
      }

      while (R < num && StringICmp (array[R]->fid, featid) == 0) {
        item = array [R];
        feat = item->feat;
        if (feat != NULL
            && !feat->ignore
            && feat->sfp != NULL
            && feat->sfp->idx.subtype == FEATDEF_mRNA
            && IdXrefsAreReciprocal(cds, feat->sfp)) {
          has_xref = TRUE;
          sfp = feat->sfp;
          context.entityID = sfp->idx.entityID;
          context.itemID = feat->itemID;
          context.sfp = sfp;
          context.sap = feat->sap;
          context.bsp = feat->bsp;
          context.label = feat->label;
          context.left = feat->left;
          context.right = feat->right;
          context.dnaStop = feat->dnaStop;
          context.partialL = feat->partialL;
          context.partialR = feat->partialR;
          context.farloc = feat->farloc;
          context.bad_order = feat->bad_order;
          context.mixed_strand = feat->mixed_strand;
          context.strand = feat->strand;
          context.seqfeattype = sfp->data.choice;;
          context.featdeftype = feat->subtype;
          context.numivals = feat->numivals;
          context.ivals = feat->ivals;
          context.omdp = (Pointer) omdp;
          context.index = R + 1;
          FindSingleMrnaProc (sfp, &context);
        }
        R++;
      }
    }
  }
  return has_xref;
}


/*
static Boolean LIBCALLBACK DummyCM121Proc (
  SeqFeatPtr sfp,
  SeqMgrFeatContextPtr context
)

{
  return TRUE;
}
*/

static void ValidateCDSmRNAmatch (
  ValidStructPtr vsp,
  BioseqPtr bsp,
  Int2 numgene,
  Int2 numcds,
  Int2 nummrna,
  Boolean suppress_duplicate_messages
)

{
  BioSourcePtr       biop;
  ValNodePtr         cdshead = NULL;
  ValNodePtr         cdstail = NULL;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext, rcontext;
  GatherContextPtr   gcp;
  GmcDataPtr         gdp, head;
  SeqFeatPtr         gene;
  Boolean            goOn, pseudo;
  GeneRefPtr         grp;
  Int2               i, j, k, numfeats, tmpnumcds, tmpnummrna, count;
  Boolean            is_genbank = FALSE;
  LpData             ld;
  Int2               num_no_mrna = 0;
  Int4               num_repeat_regions;
  Uint2              olditemtype = 0;
  Uint4              olditemid = 0;
  OrgNamePtr         onp;
  OrgRefPtr          orp;
  Int2               recip;
  VoidPtr            repeat_region_array;
  SeqFeatPtr         rpt_region;
  SeqDescrPtr        sdp;
  ErrSev             sev = /* SEV_INFO */ SEV_WARNING;
  SeqFeatPtr         sfp;
  SeqIdPtr           sip;
  VvmDataPtr         vdp;
  ValNodePtr         vnp;

  if (vsp == NULL || bsp == NULL) return;
  if (! ISA_na (bsp->mol)) return;

  gcp = vsp->gcp;
  if (gcp != NULL) {
    olditemid = gcp->itemID;
    olditemtype = gcp->thistype;
  }

  /*
  if (GetAppProperty ("ValidateCDSmRNAoneToOne") != NULL) {
    cdsMrnaOneToOne = TRUE;
  }
  */

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_OTHER) {
      sev = SEV_WARNING;
    } else if (sip->choice == SEQID_GENBANK) {
      is_genbank = TRUE;
    }
  }

  if (is_genbank) {
    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
    if (sdp != NULL) {
      biop = (BioSourcePtr) sdp->data.ptrvalue;
      if (biop != NULL) {
        orp = biop->org;
        if (orp != NULL) {
          onp = orp->orgname;
          if (onp != NULL) {
            if (StringDoesHaveText (onp->div) &&
                StringCmp (onp->div, "BCT") != 0 &&
                StringCmp (onp->div, "VRL") != 0) {
              is_genbank = FALSE;
            }
          }
        }
      }
    }
  }

  repeat_region_array = SeqMgrBuildFeatureIndex (bsp, &num_repeat_regions, 0, FEATDEF_repeat_region);

  if (numgene > 0 && numcds > 0 && nummrna > 0) {
    numfeats = numcds + nummrna;
    head = (GmcDataPtr) MemNew (sizeof (GmcData) * (size_t) (numfeats + 1));
    if (head != NULL) {
      gdp = head;
      sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
      while (sfp != NULL) {
        if (sfp->idx.subtype == FEATDEF_CDS || sfp->idx.subtype == FEATDEF_mRNA) {
          gdp->feat = sfp;
          grp = SeqMgrGetGeneXref (sfp);
          if (grp == NULL) {
            gdp->gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
          } else if (! SeqMgrGeneIsSuppressed (grp)) {
            if (StringDoesHaveText (grp->locus_tag)) {
              gdp->gene = SeqMgrGetGeneByLocusTag (bsp, grp->locus_tag, NULL);
            } else if (StringDoesHaveText (grp->locus)) {
              gdp->gene = SeqMgrGetFeatureByLabel (bsp, grp->locus, SEQFEAT_GENE, 0, NULL);
            }
          }
          gdp++;
        }
        sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
      }
      HeapSort (head, (size_t) numfeats, sizeof (GmcData), SortGmcByGenePtr);
      for (i = 0; i < numfeats; i += j) {
        gene = head [i].gene;
        for (j = 1; i + j < numfeats && gene == head [i + j].gene; j++) continue;
        if (j > 1 && gene != NULL) {
          /* is alt splicing */
          tmpnumcds = 0;
          tmpnummrna = 0;
          for (k = 0; k < j; k++) {
            sfp = head [i + k].feat;
            if (sfp == NULL) continue;
            if (sfp->idx.subtype == FEATDEF_CDS) {
              tmpnumcds++;
            }
            if (sfp->idx.subtype == FEATDEF_mRNA) {
              tmpnummrna++;
            }
          }
          if (tmpnumcds > 0 && tmpnummrna > 1 && tmpnumcds != tmpnummrna && (! is_genbank)) {

            if (gcp != NULL) {
              gcp->itemID = gene->idx.itemID;
              gcp->thistype = OBJ_SEQFEAT;
            }
            vsp->descr = NULL;
            vsp->sfp = gene;
            ValidErr (vsp, sev, ERR_SEQ_FEAT_CDSmRNAmismatch, "mRNA count (%d) does not match CDS (%d) count for gene",
                      (int) tmpnummrna, (int) tmpnumcds);
          }
        }
      }
    }
    MemFree (head);
  }

  /* loop through CDS features, finding single unused mRNA partner */

  goOn = TRUE;
  while (goOn) {
    goOn = FALSE;
    sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
    while (sfp != NULL) {
      vdp = (VvmDataPtr) sfp->idx.scratch;
      if (vdp != NULL && (! vdp->accounted_for)) {
        vdp->num_mrnas = 0;
        ld.count = 0;
        ld.cds = sfp;
        ld.mrna = NULL;
        ld.firstid [0] = '\0';
        ld.products_unique = TRUE;
        ld.featid_matched = FALSE;

        if (sfp->excpt &&
          (StringISearch (sfp->except_text, "ribosomal slippage") != NULL ||
            StringISearch (sfp->except_text, "trans-splicing") != NULL)) {
          count = SeqMgrGetAllOverlappingFeatures (sfp->location, FEATDEF_mRNA, NULL, 0,
                                                  LOCATION_SUBSET, (Pointer) &ld, FindSingleMrnaProc);
        } else {
          count = SeqMgrGetAllOverlappingFeatures (sfp->location, FEATDEF_mRNA, NULL, 0,
                                                  CHECK_INTERVALS, (Pointer) &ld, FindSingleMrnaProc);
        }

        if (!ld.featid_matched) {
          MarkMrnasFromCDSXrefs (sfp, &ld);
        }

        if (ld.count == 1 && ld.mrna != NULL) {
          vdp->accounted_for = TRUE;
          vdp->num_mrnas = ld.count;
          vdp->featid_matched = ld.featid_matched;
          vdp = (VvmDataPtr) ld.mrna->idx.scratch;
          if (vdp != NULL) {
            vdp->accounted_for = TRUE;
            goOn = TRUE;
            recip = IdXrefsNotReciprocal (sfp, ld.mrna);
            if (recip == 1) {
              if (gcp != NULL) {
                gcp->itemID = sfp->idx.itemID;
                gcp->thistype = OBJ_SEQFEAT;
              }
              vsp->descr = NULL;
              vsp->sfp = sfp;
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_SeqFeatXrefNotReciprocal, "CDS/mRNA unambiguous pair have erroneous cross-references");
            } else if (recip == 2) {
              if (gcp != NULL) {
                gcp->itemID = ld.mrna->idx.itemID;
                gcp->thistype = OBJ_SEQFEAT;
              }
              vsp->descr = NULL;
              vsp->sfp = ld.mrna;
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_SeqFeatXrefProblem, "MrnaProteinLink inconsistent with feature ID cross-references");
            }
          }
          if (SeqLocAinB (sfp->location, ld.mrna->location) < 0) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_CDSmRNAXrefLocationProblem, "CDS not contained within cross-referenced mRNA");
          }
        } else {
          vdp->num_mrnas = ld.count;
          vdp->products_unique = ld.products_unique;
          vdp->featid_matched = ld.featid_matched;
        }
      }
      sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &fcontext);
    }
  }

  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
  while (sfp != NULL && (! is_genbank)) {
    vdp = (VvmDataPtr) sfp->idx.scratch;
    if (vdp != NULL) {
      count = vdp->num_mrnas;
      /*
      count = SeqMgrGetAllOverlappingFeatures (sfp->location, FEATDEF_mRNA, NULL, 0,
                                                 CHECK_INTERVALS, NULL, DummyCM121Proc);
      */
      if (count > 1) {
        if (gcp != NULL) {
          gcp->itemID = sfp->idx.itemID;
          gcp->thistype = OBJ_SEQFEAT;
        }
        vsp->descr = NULL;
        vsp->sfp = sfp;
        if (vdp->featid_matched) {
          /* presence of reciprocal link suppresses warnings */
        } else if (vdp->products_unique) {
          /*
          if (! suppress_duplicate_messages) {
            ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_CDSwithMultipleMRNAs,
                      "CDS overlapped by %d mRNAs, but product locations are unique", (int) count);
          }
          */
          ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_CDSwithMultipleMRNAs,
                    "CDS overlapped by %d mRNAs, but product locations are unique", (int) count);
        } else {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_CDSwithMultipleMRNAs, "CDS overlapped by %d mRNAs", (int) count);
        }
      } else if (count == 0 && numgene > 0 && numcds > 0 && nummrna > 0) {
        pseudo = sfp->pseudo;
        if (! pseudo) {
          grp = SeqMgrGetGeneXref (sfp);
          if (grp != NULL) {
            pseudo = grp->pseudo;
          } else {
            gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
            if (gene != NULL) {
              pseudo = gene->pseudo;
              if (! pseudo) {
                grp = (GeneRefPtr) gene->data.value.ptrvalue;
                if (grp != NULL) {
                  pseudo = grp->pseudo;
                }
              }
            }
          }
        }
        if (! pseudo) {
          rpt_region = SeqMgrGetOverlappingFeature (sfp->location, 0, repeat_region_array, num_repeat_regions,
                                                    NULL, CONTAINED_WITHIN, &rcontext);
          if (rpt_region == NULL) {
            if (StringStr (sfp->except_text, "rearrangement required for product") == NULL) {
              /*
              if (gcp != NULL) {
                gcp->itemID = sfp->idx.itemID;
                gcp->thistype = OBJ_SEQFEAT;
              }
              vsp->descr = NULL;
              vsp->sfp = sfp;
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_CDSwithNoMRNAOverlap, "CDS overlapped by 0 mRNAs");
              */
              vnp = ValNodeAddPointer (&cdstail, 0, (Pointer) sfp);
              if (cdshead == NULL) {
                cdshead = vnp;
              }
              cdstail = vnp;
              num_no_mrna++;
            }
          }
        }
      }
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &fcontext);
  }

  MemFree (repeat_region_array);

  if (num_no_mrna > 0) {
    if (num_no_mrna >= 10) {
      if (gcp != NULL) {
        gcp->itemID = olditemid;
        gcp->thistype = olditemtype;
      }
      vsp->descr = NULL;
      vsp->sfp = NULL;
      vsp->bsp = bsp;
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_CDSwithNoMRNAOverlap,
                "%d out of %d CDSs overlapped by 0 mRNAs", (int) num_no_mrna, (int) numcds);
    } else {
      for (vnp = cdshead; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp == NULL) continue;
        if (gcp != NULL) {
          gcp->itemID = sfp->idx.itemID;
          gcp->thistype = OBJ_SEQFEAT;
        }
        vsp->descr = NULL;
        vsp->sfp = sfp;
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_CDSwithNoMRNAOverlap, "CDS overlapped by 0 mRNAs");
      }
    }
  }

  ValNodeFree (cdshead);

  if (gcp != NULL) {
    gcp->itemID = olditemid;
    gcp->thistype = olditemtype;
  }
}

static Boolean HaveUniqueFeatIDXrefs (SeqFeatXrefPtr xref1, SeqFeatXrefPtr xref2)

{
  ObjectIdPtr  oip1 = NULL, oip2 = NULL;

  while (xref1 != NULL) {
    if (xref1->id.choice == 3) {
      oip1 = (ObjectIdPtr) xref1->id.value.ptrvalue;
    }
    xref1 = xref1->next;
  }

  while (xref2 != NULL) {
    if (xref2->id.choice == 3) {
      oip2 = (ObjectIdPtr) xref2->id.value.ptrvalue;
    }
    xref2 = xref2->next;
  }

  if (oip1 == NULL || oip2 == NULL) return FALSE;
  if (oip1->str == NULL && oip2->str == NULL) {
    if (oip1->id != oip2->id && oip1->id > 0 && oip2->id > 0) return TRUE;
  }

  return FALSE;
}

#define LEFT_RIBOSOMAL_SUBUNIT  1
#define INTERNAL_SPACER_1        2
#define MIDDLE_RIBOSOMAL_SUBUNIT 3
#define INTERNAL_SPACER_2        4
#define RIGHT_RIBOSOMAL_SUBUNIT  5
#define INTERNAL_SPACER_X        6

static Int2 WhichRNA (SeqFeatPtr sfp)

{
  RnaRefPtr  rrp;
  CharPtr    str;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA) return 0;
  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (rrp == NULL) return 0;
  str = GetRNARefProductString (rrp, NULL);
  if (StringHasNoText (str)) return 0;
  if (rrp->type == RNA_TYPE_rRNA) {
    if (StringNICmp (str, "small ", 6) == 0) return LEFT_RIBOSOMAL_SUBUNIT;
    if (StringNICmp (str, "18S ", 4) == 0) return LEFT_RIBOSOMAL_SUBUNIT;
    if (StringNICmp (str, "16S ", 4) == 0) return LEFT_RIBOSOMAL_SUBUNIT;
    if (StringNICmp (str, "5.8S ", 5) == 0) return MIDDLE_RIBOSOMAL_SUBUNIT;
    if (StringNICmp (str, "large ", 6) == 0) return RIGHT_RIBOSOMAL_SUBUNIT;
    if (StringNICmp (str, "26S ", 4) == 0) return RIGHT_RIBOSOMAL_SUBUNIT;
    if (StringNICmp (str, "28S ", 4) == 0) return RIGHT_RIBOSOMAL_SUBUNIT;
    if (StringNICmp (str, "23S ", 4) == 0) return RIGHT_RIBOSOMAL_SUBUNIT;
    /* variant spellings */
    if (StringNICmp (str, "18 ", 3) == 0) return LEFT_RIBOSOMAL_SUBUNIT;
    if (StringNICmp (str, "16 ", 3) == 0) return LEFT_RIBOSOMAL_SUBUNIT;
    if (StringNICmp (str, "5.8 ", 4) == 0) return MIDDLE_RIBOSOMAL_SUBUNIT;
    if (StringNICmp (str, "26 ", 3) == 0) return RIGHT_RIBOSOMAL_SUBUNIT;
    if (StringNICmp (str, "28 ", 3) == 0) return RIGHT_RIBOSOMAL_SUBUNIT;
    if (StringNICmp (str, "23 ", 3) == 0) return RIGHT_RIBOSOMAL_SUBUNIT;
  }
  if (rrp->type == RNA_TYPE_misc_RNA) {
    if (StringICmp (str, "internal transcribed spacer 1") == 0) return INTERNAL_SPACER_1;
    if (StringICmp (str, "internal transcribed spacer 2") == 0) return INTERNAL_SPACER_2;
    /* variant spellings */
    if (StringICmp (str, "internal transcribed spacer1") == 0) return INTERNAL_SPACER_1;
    if (StringICmp (str, "internal transcribed spacer2") == 0) return INTERNAL_SPACER_2;
    if (StringICmp (str, "internal transcribed spacer") == 0) return INTERNAL_SPACER_X;
    if (StringICmp (str, "ITS") == 0) return INTERNAL_SPACER_X;
    if (StringICmp (str, "16S-23S ribosomal RNA intergenic spacer") == 0) return INTERNAL_SPACER_X;
    if (StringICmp (str, "16S-23S intergenic spacer") == 0) return INTERNAL_SPACER_X;
    if (StringICmp (str, "intergenic spacer") == 0) return INTERNAL_SPACER_X;
  }
  return 0;
}

static Boolean CDSsLinkedToDifferentMRNAs (SeqFeatPtr sfp, SeqFeatPtr last)

{
  SeqFeatPtr      mrna1 = NULL, mrna2 = NULL;
  SeqFeatXrefPtr  xref;

  if (sfp == NULL || last == NULL) return FALSE;
  if (sfp->idx.subtype != FEATDEF_CDS || last->idx.subtype != FEATDEF_CDS) return FALSE;

  for (xref = sfp->xref; xref != NULL && mrna1 == NULL; xref = xref->next) {
    if (xref->id.choice != 0) {
      mrna1 = SeqMgrGetFeatureByFeatID (sfp->idx.entityID, NULL, NULL, xref, NULL);
      if (mrna1 != NULL && mrna1->idx.subtype != FEATDEF_mRNA) {
        mrna1 = NULL;
      }
    }
  }

  for (xref = last->xref; xref != NULL && mrna2 == NULL; xref = xref->next) {
    if (xref->id.choice != 0) {
      mrna2 = SeqMgrGetFeatureByFeatID (last->idx.entityID, NULL, NULL, xref, NULL);
      if (mrna2 != NULL && mrna2->idx.subtype != FEATDEF_mRNA) {
        mrna2 = NULL;
      }
    }
  }

  if (mrna1 != NULL && mrna2 != NULL && mrna1 != mrna2) return TRUE;

  return FALSE;
}

static Boolean MRNAsLinkedToDifferentCDSs (SeqFeatPtr sfp, SeqFeatPtr last)

{
  SeqFeatPtr      cds1 = NULL, cds2 = NULL;
  CdRegionPtr     crp1, crp2;
  SeqFeatXrefPtr  xref;

  if (sfp == NULL || last == NULL) return FALSE;
  if (sfp->idx.subtype != FEATDEF_mRNA || last->idx.subtype != FEATDEF_mRNA) return FALSE;

  for (xref = sfp->xref; xref != NULL && cds1 == NULL; xref = xref->next) {
    if (xref->id.choice != 0) {
      cds1 = SeqMgrGetFeatureByFeatID (sfp->idx.entityID, NULL, NULL, xref, NULL);
      if (cds1 != NULL && cds1->idx.subtype != FEATDEF_CDS) {
        cds1 = NULL;
      }
    }
  }

  for (xref = last->xref; xref != NULL && cds2 == NULL; xref = xref->next) {
    if (xref->id.choice != 0) {
      cds2 = SeqMgrGetFeatureByFeatID (last->idx.entityID, NULL, NULL, xref, NULL);
      if (cds2 != NULL && cds2->idx.subtype != FEATDEF_CDS) {
        cds2 = NULL;
      }
    }
  }

  if (cds1 == NULL || cds2 == NULL || cds1 == cds2) return FALSE;

  crp1 = (CdRegionPtr) cds1->data.value.ptrvalue;
  crp2 = (CdRegionPtr) cds2->data.value.ptrvalue;
  if (crp1 == NULL || crp2 == NULL) return FALSE;

  if (SeqLocCompare (cds1->location, cds2->location) != SLC_A_EQ_B) return TRUE;

  if (crp1->frame < 2 && crp2->frame < 2) return FALSE;
  if (crp1->frame != crp2->frame) return TRUE;

  return FALSE;
}

static Boolean BaseRangeIsVirtual (BioseqPtr bsp, Int4 left, Int4 right)

{
  Uint1        res;
  StreamCache  sc;

  if (! StreamCacheSetup (bsp, NULL, EXPAND_GAPS_TO_DASHES, &sc)) return FALSE;

  StreamCacheSetPosition (&sc, left - 1);
  res = StreamCacheGetResidue (&sc);
  if (res == '-') return FALSE;
  res = StreamCacheGetResidue (&sc);
  if (res != '-') return FALSE;

  StreamCacheSetPosition (&sc, right - 1);
  res = StreamCacheGetResidue (&sc);
  if (res != '-') return FALSE;
  res = StreamCacheGetResidue (&sc);
  if (res == '-') return FALSE;

  return TRUE;
}

static Boolean LIBCALLBACK GetFeatsInGaps (
  SeqFeatPtr sfp,
  SeqMgrFeatContextPtr fcontext
)

{
  BioseqPtr         bsp;
  Int4              dashes;
  Int2              first = 0;
  GatherContextPtr  gcp;
  Int2              last = 0;
  Int4              len;
  SeqLocPtr         loc;
  Int2              localfirst;
  Int2              locallast;
  ErrSev            logsev;
  ErrSev            msgsev;
  Boolean           needToStream = TRUE;
  Int4              Ns;
  Uint2             olditemtype = 0;
  Uint4             olditemid = 0;
  Int4              plusses;
  Int2              prefix = 0;
  Int2              suffix = 0;
  Int4              realBases;
  Int2              res;
  StreamCache       sc;
  SeqIntPtr         sintp;
  SeqLocPtr         slp;
  Boolean           startsOrEndsInGap = FALSE;
  ValidStructPtr    vsp;

  if (sfp == NULL || fcontext == NULL) return FALSE;
  vsp = (ValidStructPtr) fcontext->userdata;
  if (vsp == NULL) return FALSE;
  gcp = vsp->gcp;
  if (gcp == NULL) return FALSE;

  if (sfp->idx.subtype == FEATDEF_gap) return TRUE;
  loc = sfp->location;
  if (loc == NULL) return TRUE;

  olditemid = gcp->itemID;
  olditemtype = gcp->thistype;

  gcp->itemID = fcontext->itemID;
  gcp->thistype = OBJ_SEQFEAT;
  vsp->sfp = sfp;


  dashes = 0;
  plusses = 0;
  Ns = 0;
  realBases = 0;

  msgsev = ErrSetMessageLevel (SEV_MAX);
  logsev = ErrSetLogLevel (SEV_MAX);

  /* special check for single interval misc_features that may exactly cover a gap */
  if (loc->choice == SEQLOC_INT && sfp->idx.subtype == FEATDEF_misc_feature) {
    sintp = (SeqIntPtr) loc->data.ptrvalue;
    if (sintp != NULL) {
      bsp = BioseqFind (sintp->id);
      if (bsp != NULL && sintp->from > 0 && sintp->to < bsp->length - 1) {
        len = SeqLocLen (loc);
        if (StreamCacheSetup (bsp, NULL, EXPAND_GAPS_TO_DASHES | KNOWN_GAP_AS_PLUS, &sc)) {
          StreamCacheSetPosition (&sc, sintp->from - 1);
          prefix = StreamCacheGetResidue (&sc);
          while ((res = StreamCacheGetResidue (&sc)) != '\0' && len > 0) {
            if (IS_LOWER (res)) {
              res = TO_UPPER (res);
            }
            if (first == 0) {
              first = res;
            }
            last = res;
            if (res == '-') {
              dashes++;
            } else if (res == '+') {
              plusses++;
            } else if (res == 'N') {
              Ns++;
            } else if (res != 0) {
              realBases++;
            }
            len--;
          }
          suffix = StreamCacheGetResidue (&sc);
          needToStream = FALSE;
        }
      }
    }
  }

  /*
  if (needToStream && StreamCacheSetup (NULL, loc, EXPAND_GAPS_TO_DASHES | KNOWN_GAP_AS_PLUS, &sc)) {
    while ((res = StreamCacheGetResidue (&sc)) != '\0') {
      if (IS_LOWER (res)) {
        res = TO_UPPER (res);
      }
      if (first == 0) {
        first = res;
      }
      last = res;
      if (res == '-') {
        dashes++;
      } else if (res == '+') {
        plusses++;
      } else if (res == 'N') {
        Ns++;
      } else if (res != 0) {
        realBases++;
      }
    }
  }
  */

  if (needToStream) {
    for (slp = SeqLocFindNext (loc, NULL); slp != NULL; slp = SeqLocFindNext (loc, slp)) {
      if (StreamCacheSetup (NULL, slp, EXPAND_GAPS_TO_DASHES | KNOWN_GAP_AS_PLUS, &sc)) {
        localfirst = 0;
        locallast = 0;
        while ((res = StreamCacheGetResidue (&sc)) != '\0') {
          if (IS_LOWER (res)) {
            res = TO_UPPER (res);
          }
          if (first == 0) {
            first = res;
          }
          if (localfirst == 0) {
            localfirst = res;
          }
          last = res;
          locallast = res;
          if (res == '-') {
            dashes++;
          } else if (res == '+') {
            plusses++;
          } else if (res == 'N') {
            Ns++;
          } else if (res != 0) {
            realBases++;
          }
        }
        if (localfirst == '-' || localfirst == '+' || locallast == '-' || locallast == '+') {
          startsOrEndsInGap = TRUE;
        }
      }
    }
  }

  ErrSetLogLevel (logsev);
  ErrSetMessageLevel (msgsev);

  if (dashes == 0 && plusses == 0 && Ns == 0) {
    /* ignore features that do not cover any gap characters */
  } else if (first == '-' || first == '+' || last == '-' || last == '+') {
    if (realBases > 0) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FeatureBeginsOrEndsInGap, "Feature begins or ends in gap");
    } else if (IS_ALPHA (prefix) && IS_ALPHA (suffix)) {
      /* ignore (misc_) features that exactly cover the gap */
    } else {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FeatureInsideGap, "Feature inside sequence gap");
    }
  } else if (realBases == 0 && dashes == 0 && plusses == 0 && Ns >= 50) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FeatureInsideGap, "Feature inside gap of Ns");
  } else if ((sfp->data.choice == SEQFEAT_CDREGION || sfp->data.choice == SEQFEAT_RNA) && dashes > 0) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FeatureCrossesGap, "Feature crosses gap of unknown length");
  } else if (startsOrEndsInGap) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_IntervalBeginsOrEndsInGap, "Internal interval begins or ends in gap");
  }

  gcp->itemID = olditemid;
  gcp->thistype = olditemtype;
  vsp->sfp = NULL;

  return TRUE;
}

static void CheckBioseqForFeatsInGap (
  BioseqPtr bsp,
  ValidStructPtr vsp
)

{
  SeqMgrFeatContext  fcontext;
  SeqFeatPtr         sfp;
  SeqIdPtr           sip;

  if (bsp == NULL || ISA_aa (bsp->mol)) return;
  sip = SeqIdFindBest (bsp->id, 0);
  if (sip == NULL) return;

  for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext)) {
    fcontext.userdata = (Pointer) vsp;
    GetFeatsInGaps (sfp, &fcontext);
  }
}

static Boolean ReportGeneCollision (GeneRefPtr grp, GeneRefPtr lastgrp)

{
  if (grp == NULL || lastgrp == NULL) return TRUE;

  if (StringDoesHaveText (grp->locus) && StringDoesHaveText (lastgrp->locus)) {
    if (StringICmp (grp->locus, lastgrp->locus) == 0) return TRUE;
  }

  if (StringDoesHaveText (grp->locus_tag) && StringDoesHaveText (lastgrp->locus_tag)) {
    if (StringICmp (grp->locus_tag, lastgrp->locus_tag) == 0) return TRUE;
  }

  if (StringDoesHaveText (grp->desc) && StringDoesHaveText (lastgrp->desc)) {
    if (StringICmp (grp->desc, lastgrp->desc) == 0) return FALSE;
  }

  return TRUE;
}

static Boolean FeatureSequencesIdentical (SeqFeatPtr sfp, SeqFeatPtr lastsfp)

{
  Boolean  rsult = FALSE;
  CharPtr  tmp1, tmp2;

  if (sfp == NULL || lastsfp == NULL) return rsult;

  tmp1 = GetSequenceByFeature (sfp);
  tmp2 = GetSequenceByFeature (lastsfp);

  if (tmp1 != NULL && tmp2 != NULL) {
    if (StringCmp (tmp1, tmp2) == 0) {
      rsult = TRUE;
    }
  }

  MemFree (tmp1);
  MemFree (tmp2);

  return rsult;
}

static Boolean GeneXrefsDifferent (SeqFeatPtr sfp, SeqFeatPtr lastsfp)

{
  SeqFeatPtr  gene, lastgene;

  if (sfp == NULL || lastsfp == NULL) return FALSE;

  gene = GetGeneForFeature (sfp);
  lastgene = GetGeneForFeature (lastsfp);
  if (gene == NULL || lastgene == NULL) return FALSE;

  if (gene != lastgene) return TRUE;

  return FALSE;
}


static void CheckForGenesInconsistent (BioseqPtr bsp, ValidStructPtr vsp)
{
  SeqMgrFeatContext  fcontext;
  SeqFeatPtr mrna, genomic_gene = NULL, mrna_gene = NULL;
  GeneRefPtr genomicgrp, grp, found_grp;
  BioseqPtr  genomic_bsp;
  Boolean    is_error = FALSE;

  mrna = SeqMgrGetRNAgivenProduct (bsp, &fcontext);
  if (mrna != NULL) {
    genomicgrp = SeqMgrGetGeneXref (mrna);
    if (genomicgrp == NULL) {
      genomic_gene = SeqMgrGetOverlappingGene (mrna->location, NULL);
      if (genomic_gene != NULL) {
        genomicgrp = (GeneRefPtr) genomic_gene->data.value.ptrvalue;
      }
    }
    if (genomicgrp != NULL
        && (mrna_gene = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_GENE, 0, &fcontext)) != NULL
        && (grp = (GeneRefPtr) mrna_gene->data.value.ptrvalue) != NULL) {
      if (StringCmp (grp->locus, genomicgrp->locus) != 0 ||
          StringCmp (grp->allele, genomicgrp->allele) != 0 ||
          StringCmp (grp->desc, genomicgrp->desc) != 0 ||
          StringCmp (grp->locus_tag, genomicgrp->locus_tag) != 0) {
        is_error = TRUE;
        if (genomic_gene == NULL
            && ((StringHasNoText (genomicgrp->desc) && StringDoesHaveText(grp->desc))
                || (StringHasNoText (genomicgrp->allele) && StringDoesHaveText(grp->allele)))) {
          genomic_bsp = BioseqFindFromSeqLoc (mrna->location);
          if (StringDoesHaveText (genomicgrp->locus_tag)) {
            genomic_gene = SeqMgrGetGeneByLocusTag (genomic_bsp, genomicgrp->locus_tag, &fcontext);
            if (genomic_gene != NULL && (found_grp = (GeneRefPtr)genomic_gene->data.value.ptrvalue)
              && StringCmp (found_grp->locus, genomicgrp->locus) != 0) {
              genomic_gene = NULL;
            }
          } else if (StringDoesHaveText (genomicgrp->locus)) {
            genomic_gene = SeqMgrGetFeatureByLabel (genomic_bsp, genomicgrp->locus, SEQFEAT_GENE, 0, &fcontext);
            if (genomic_gene != NULL
                && (found_grp = (GeneRefPtr)genomic_gene->data.value.ptrvalue) != NULL
                && StringDoesHaveText (found_grp->locus_tag)) {
              genomic_gene = NULL;
            }
          }
          if (genomic_gene != NULL && (genomicgrp = (GeneRefPtr) genomic_gene->data.value.ptrvalue) != NULL
              && StringCmp (grp->locus, genomicgrp->locus) == 0
              && StringCmp (grp->allele, genomicgrp->allele) == 0
              && StringCmp (grp->desc, genomicgrp->desc) == 0
              && StringCmp (grp->locus_tag, genomicgrp->locus_tag) == 0) {
            is_error = FALSE;
          }
        }
        if (is_error) {
          vsp->descr = NULL;
          vsp->sfp = mrna_gene;
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_GenesInconsistent, "Gene on mRNA bioseq does not match gene on genomic bioseq");
        }
      }
    }
  }
}

static void CheckForNonViralComplete (BioseqPtr bsp, ValidStructPtr vsp, GatherContextPtr gcp)

{
  BioSourcePtr       biop = NULL;
  SeqMgrDescContext  dcontext;
  MolInfoPtr         mip = NULL;
  Uint2              olditemtype = 0;
  Uint4              olditemid = 0;
  OrgNamePtr         onp;
  OrgRefPtr          orp;
  ObjValNodePtr      ovp;
  SeqDescrPtr        sdp;
  CharPtr            title = NULL;
  SubSourcePtr       ssp;

  if (bsp == NULL || vsp == NULL) return;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
  if (sdp != NULL) {
    mip = (MolInfoPtr) sdp->data.ptrvalue;
  }
  if (mip == NULL) return;
  if (mip->biomol != MOLECULE_TYPE_GENOMIC || mip->completeness != 1) return;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_title, &dcontext);
  if (sdp != NULL) {
    title = (CharPtr) sdp->data.ptrvalue;
  }
  if (title == NULL) return;
  if (StringStr (title, "complete genome") == NULL) return;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp != NULL) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
  }
  if (biop == NULL) return;
  if (biop->genome != GENOME_genomic && biop->genome != GENOME_unknown) return;

  orp = biop->org;
  if (orp == NULL) return;
  onp = orp->orgname;
  if (onp == NULL) return;
  if (StringNICmp (onp->lineage, "Viruses; ", 9) == 0) return;
  if (StringNICmp (onp->lineage, "Viroids; ", 9) == 0) return;
  if (StringICmp (onp->lineage, "Viruses") == 0 && StringICmp (onp->div, "PHG") == 0) return;

  for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
    if (ssp->subtype == SUBSRC_endogenous_virus_name) {
      return;
    }
  }

  if (gcp != NULL) {
    olditemid = gcp->itemID;
    olditemtype = gcp->thistype;
  }

  if (sdp != NULL && sdp->extended != 0) {
    ovp = (ObjValNodePtr) sdp;
    if (ovp != NULL && gcp != NULL) {
      gcp->itemID = ovp->idx.itemID;
      gcp->thistype = ovp->idx.itemtype;
    }
  }

  ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceNeedsChromosome, "Non-viral complete genome not labeled as chromosome");

  if (gcp != NULL) {
    gcp->itemID = olditemid;
    gcp->thistype = olditemtype;
  }
}


static void LookForViralMolInfoInLineage (BioseqPtr bsp, ValidStructPtr vsp, GatherContextPtr gcp)
{
  SeqDescPtr        sdp;
  SeqMgrDescContext context;
  BioSourcePtr      biop;
  CharPtr           lineage;
  ObjValNodePtr     ovp;
  Uint2             olditemtype = 0;
  Uint4             olditemid = 0;

  if (bsp == NULL || ISA_aa (bsp->mol) || vsp == NULL) return;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
  if (sdp == NULL || (biop = (BioSourcePtr) sdp->data.ptrvalue) == NULL
      || biop->org == NULL || biop->org->orgname == NULL
      || (lineage = biop->org->orgname->lineage) == NULL
      || StringNICmp (biop->org->orgname->lineage, "Viruses; ", 9) != 0) {
    return;
  }

  if (gcp != NULL) {
    olditemid = gcp->itemID;
    olditemtype = gcp->thistype;
  }

  if (sdp != NULL && sdp->extended != 0) {
    ovp = (ObjValNodePtr) sdp;
    if (ovp != NULL && gcp != NULL) {
      gcp->itemID = ovp->idx.itemID;
      gcp->thistype = ovp->idx.itemtype;
    }
  }

  vsp->bsp = bsp;
  vsp->descr = sdp;

  if ((StringSearch (lineage, " ssRNA viruses; ") != NULL
       || StringSearch (lineage, " ssRNA negative-strand viruses; ") != NULL
       || StringSearch (lineage, " ssRNA positive-strand viruses, no DNA stage; ") != NULL
       || StringSearch (lineage, " unassigned ssRNA viruses; ") != NULL)
      && (bsp->mol != Seq_mol_rna)) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MolInfoConflictsWithBioSource, "Taxonomy indicates single-stranded RNA, sequence does not agree.");
  }

  if (StringSearch (lineage, " dsRNA viruses; ") != NULL
      && (bsp->mol != Seq_mol_rna)) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MolInfoConflictsWithBioSource, "Taxonomy indicates double-stranded RNA, sequence does not agree.");
  }

  if (StringSearch (lineage, " ssDNA viruses; ") != NULL
      && (bsp->mol != Seq_mol_dna)) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MolInfoConflictsWithBioSource, "Taxonomy indicates single-stranded DNA, sequence does not agree.");
  }

  if (StringSearch (lineage, " dsDNA viruses; ") != NULL
      && (bsp->mol != Seq_mol_dna)) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MolInfoConflictsWithBioSource, "Taxonomy indicates double-stranded DNA, sequence does not agree.");
  }

  if (gcp != NULL) {
    gcp->itemID = olditemid;
    gcp->thistype = olditemtype;
  }
}


static Boolean ValidateBioseqContextIndexed (BioseqPtr bsp, BioseqValidStrPtr bvsp)

{
  ValidStructPtr     vsp;
  ObjMgrDataPtr      omdp;
  SeqSubmitPtr       ssp;
  SubmitBlockPtr     sbp;
  GatherContextPtr   gcp;
  SeqFeatPtr         sfp, lastsfp;
  SeqMgrFeatContext  fcontext;
  Uint2              featdeftype = 0;
  Boolean            firstCDS;
  GeneRefPtr         grp, lastgrp;
  SeqFeatPtr         last = NULL;
  Boolean            leave;
  CharPtr            label = NULL;
  CharPtr            comment = NULL;
  Int4               left = 0;
  Boolean            partialL = FALSE;
  Boolean            partialR = FALSE;
  Int4               right = 0;
  Uint1              strand = 0;
  Int2               numivals = 0;
  Int4Ptr            ivals = NULL;
  Boolean            ivalssame;
  SeqAnnotPtr        sap = NULL;
  Uint2              olditemtype = 0;
  Uint4              olditemid = 0;
  CharPtr            lastLabel;
  CharPtr            message;
  Int2               i;
  Boolean            isCuratedFlybase = FALSE;
  Boolean            isDrosophila = FALSE;
  Boolean            isGenBankAccn = FALSE;
  Boolean            isGeneralAccn = FALSE;
  Boolean            isGPSorNTorNCorNGorNW = FALSE;
  Boolean            isViral = FALSE;
  Int2               j;
  CdRegionPtr        crp;
  Uint1              frame = 0;
  Boolean            samelabel;
  int                severity;
  int                overlapPepSev;
  BioSourcePtr       biop = NULL, lastbiop;
  OrgRefPtr          orp = NULL;
  OrgNamePtr         onp = NULL;
  Int4               fiveUTRright;
  Int4               cdsRight;
  Int4               threeUTRright;
  Int4               cdscount, genecount, utr5count, utr3count;
  SeqFeatPtr         gene, cdsgene, utr5gene, utr3gene;
  PubdescPtr         pdp = NULL, lastpdp;
  SeqDescrPtr        sdp;
  SeqMgrDescContext  dcontext;
  Boolean            showBadFullSource;
  Int2               numBadFullSource;
  SubSourcePtr       sbsp;
  Int2               numgene, numcds, nummrna, numcdsproducts, nummrnaproducts,
                     numcdspseudo, nummrnapseudo, numrearrangedcds, lastrnatype,
                     thisrnatype;
  Boolean            cds_products_unique = TRUE, mrna_products_unique = TRUE,
                     suppress_duplicate_messages = FALSE, pseudo;
  SeqIdPtr           sip;
  Char               buf [96];
  SeqFeatXrefPtr     xref = NULL;
  CharPtr            except_text = NULL;
  ValNodePtr         vnp, cds_prod_head = NULL, mrna_prod_head = NULL,
                     lastcdsprod = NULL, lastmrnaprod = NULL;
  StreamCache        sc;
  Int2               res;
  Int4               dashes;
  Int4               Ns;
  Int4               realBases;
  Int4               estimated_length;
  Int4               loclen;
  GBQualPtr          gbq;
  long int           val;
  SeqLocPtr          slp;
  MolInfoPtr         mip = NULL;
  SeqFeatPtr         cds;
  BioseqPtr          nbsp;
  Boolean            last_reported;
  Boolean            found_overlapping_peptide;

  gcp = bvsp->gcp;
  vsp = bvsp->vsp;
  vsp->descr = NULL;
  vsp->sfp = NULL;
  vsp->gcp = gcp;               /* needed for ValidErr */

  vsp->rrna_array = SeqMgrBuildFeatureIndex (bsp, &(vsp->numrrna), 0, FEATDEF_rRNA);
  vsp->trna_array = SeqMgrBuildFeatureIndex (bsp, &(vsp->numtrna), 0, FEATDEF_tRNA);

  SeqMgrExploreFeatures (bsp, (Pointer) bvsp, ValidateSeqFeatIndexed, NULL, NULL, NULL);

  vsp->rrna_array = MemFree (vsp->rrna_array);
  vsp->trna_array = MemFree (vsp->trna_array);

  overlapPepSev = SEV_WARNING;
  if (GetAppProperty ("SpliceValidateAsError") != NULL) {
    overlapPepSev = SEV_ERROR;
  }

  if (gcp != NULL) {
    olditemid = gcp->itemID;
    olditemtype = gcp->thistype;
  }

  numgene = 0;
  numcds = 0;
  nummrna = 0;
  numcdsproducts = 0;
  nummrnaproducts = 0;
  numcdspseudo = 0;
  nummrnapseudo = 0;
  numrearrangedcds = 0;

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  while (sfp != NULL) {
    switch (sfp->idx.subtype) {
      case FEATDEF_GENE :
        numgene++;
        break;
      case FEATDEF_CDS :
        numcds++;
        if (StringStr (sfp->except_text, "rearrangement required for product") != NULL) {
          numrearrangedcds++;
        }
        if (sfp->product != NULL) {
          numcdsproducts++;
          sip = SeqLocId (sfp->product);
          SeqIdWrite (sip, buf, PRINTID_FASTA_LONG, sizeof (buf) - 1);
          if (StringDoesHaveText (buf)) {
            vnp = ValNodeCopyStr (&lastcdsprod, 0, buf);
            if (cds_prod_head == NULL) {
              cds_prod_head = vnp;
            }
            lastcdsprod = vnp;
          }
        } else {
          pseudo = sfp->pseudo;
          if (! pseudo) {
            grp = SeqMgrGetGeneXref (sfp);
            if (grp != NULL) {
              pseudo = grp->pseudo;
            } else {
              gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
              if (gene != NULL) {
                pseudo = gene->pseudo;
                if (! pseudo) {
                  grp = (GeneRefPtr) gene->data.value.ptrvalue;
                  if (grp != NULL) {
                    pseudo = grp->pseudo;
                  }
                }
              }
            }
          }
          if (pseudo) {
            numcdspseudo++;
          }
        }
        break;
      case FEATDEF_mRNA :
        nummrna++;
        if (sfp->product != NULL) {
          nummrnaproducts++;
          sip = SeqLocId (sfp->product);
          SeqIdWrite (sip, buf, PRINTID_FASTA_LONG, sizeof (buf) - 1);
          if (StringDoesHaveText (buf)) {
            vnp = ValNodeCopyStr (&lastmrnaprod, 0, buf);
            if (mrna_prod_head == NULL) {
              mrna_prod_head = vnp;
            }
            lastmrnaprod = vnp;
          }
        } else {
          pseudo = sfp->pseudo;
          if (! pseudo) {
            grp = SeqMgrGetGeneXref (sfp);
            if (grp != NULL) {
              pseudo = grp->pseudo;
            } else {
              gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
              if (gene != NULL) {
                pseudo = gene->pseudo;
                if (! pseudo) {
                  grp = (GeneRefPtr) gene->data.value.ptrvalue;
                  if (grp != NULL) {
                    pseudo = grp->pseudo;
                  }
                }
              }
            }
          }
          if (pseudo) {
            nummrnapseudo++;
          }
        }
        break;
      default :
        break;
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
  }

  vsp->descr = NULL;
  vsp->sfp = NULL;
  vsp->bssp = NULL;
  vsp->bsp = bsp;

  if (cds_prod_head != NULL) {
    cds_prod_head = ValNodeSort (cds_prod_head, SortVnpByString);
    cds_products_unique = StringListIsUnique (cds_prod_head);
  }
  if (mrna_prod_head != NULL) {
    mrna_prod_head = ValNodeSort (mrna_prod_head, SortVnpByString);
    mrna_products_unique = StringListIsUnique (mrna_prod_head);;
  }

  if (numcds > 0 && nummrna > 1) {
    if (numcdsproducts + numcdspseudo == numcds &&
        (nummrnaproducts + nummrnapseudo == nummrna || nummrnaproducts == 0) &&
        cds_products_unique && mrna_products_unique) {
      suppress_duplicate_messages = TRUE;
    }
    if (numcdsproducts > 0 && numcdsproducts + numcdspseudo != numcds && numcdsproducts + numcdspseudo + numrearrangedcds != numcds) {
      if (gcp != NULL) {
        gcp->itemID = olditemid;
        gcp->thistype = olditemtype;
      }
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FeatureProductInconsistency, "%d CDS features have %d product references",
                (int) numcds, (int) numcdsproducts);
    }
    if (numcdsproducts > 0 && (! cds_products_unique)) {
      if (gcp != NULL) {
        gcp->itemID = olditemid;
        gcp->thistype = olditemtype;
      }
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FeatureProductInconsistency, "CDS products are not unique");
    }
    if (nummrnaproducts > 0 && nummrnaproducts + nummrnapseudo != nummrna) {
      if (gcp != NULL) {
        gcp->itemID = olditemid;
        gcp->thistype = olditemtype;
      }
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FeatureProductInconsistency, "%d mRNA features have %d product references",
                (int) nummrna, (int) nummrnaproducts);
    }
    if (nummrnaproducts > 0 && (! mrna_products_unique)) {
      if (gcp != NULL) {
        gcp->itemID = olditemid;
        gcp->thistype = olditemtype;
      }
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FeatureProductInconsistency, "mRNA products are not unique");
    }
    /*
    if (numcds > nummrna) {
      if (gcp != NULL) {
        gcp->itemID = olditemid;
        gcp->thistype = olditemtype;
      }
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FeatureProductInconsistency, "%d CDS features and only %d mRNA features",
                (int) numcds, (int) nummrna);
    } else if (numcds < nummrna) {
      if (gcp != NULL) {
        gcp->itemID = olditemid;
        gcp->thistype = olditemtype;
      }
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FeatureProductInconsistency, "%d mRNA features and only %d CDS features",
                (int) nummrna, (int) numcds);
    }
    */
  }

  ValNodeFreeData (cds_prod_head);
  ValNodeFreeData (mrna_prod_head);

  /*
  SeqEntryToBioSource (vsp->sep, NULL, NULL, 0, &biop);
  */
  BioseqToGeneticCode (bsp, NULL, NULL, NULL, NULL, 0, &biop);
  if (biop != NULL) {
    orp = biop->org;
    if (orp != NULL) {
      /* curated fly source still has duplicate features */
      if (StringICmp (orp->taxname, "Drosophila melanogaster") == 0) {
        isDrosophila = TRUE;
      }
      onp = orp->orgname;
      if (onp != NULL) {
        if (StringNICmp (onp->lineage, "Viruses; ", 9) == 0) {
          isViral = TRUE;
        }
      }
    }
  }
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
  if (sdp != NULL && sdp->choice == Seq_descr_molinfo) {
    mip = (MolInfoPtr) sdp->data.ptrvalue;
  }

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  last_reported = FALSE;
  while (sfp != NULL) {
    HasFeatId(sfp, 932);
    leave = TRUE;
    if (last != NULL) {
      ivalssame = FALSE;
      if (fcontext.left == left && fcontext.right == right && fcontext.featdeftype == featdeftype) {
        if ((fcontext.strand == Seq_strand_minus && strand == Seq_strand_minus)
            || (fcontext.strand != Seq_strand_minus && strand != Seq_strand_minus)) {
          ivalssame = TRUE;
          if (fcontext.numivals != numivals || fcontext.ivals == NULL || ivals == NULL) {
            ivalssame = FALSE;
          } else {
            for (i = 0, j = 0; i < numivals; i++, j += 2) {
              if (fcontext.ivals[j] != ivals[j]) {
                ivalssame = FALSE;
              }
              if (fcontext.ivals[j + 1] != ivals[j + 1]) {
                ivalssame = FALSE;
              }
            }
          }
          if (ivalssame &&      /* StringICmp (fcontext.label, label) == 0 && */
              (fcontext.sap == sap || (fcontext.sap->desc == NULL && sap->desc == NULL) || DescsSame (fcontext.sap->desc, sap->desc))) {
            if (gcp != NULL) {
              gcp->itemID = fcontext.itemID;
              gcp->thistype = OBJ_SEQFEAT;
            }
            vsp->descr = NULL;
            vsp->sfp = sfp;
            severity = SEV_ERROR;
            samelabel = TRUE;
            if (StringICmp (fcontext.label, label) != 0 || StringICmp (sfp->comment, comment) != 0) {
              samelabel = FALSE;
            }
            if (GBQualsDiffer (sfp->qual, last->qual)) {
              samelabel = FALSE;
            }
            if (featdeftype == FEATDEF_PUB ||
                featdeftype == FEATDEF_REGION || featdeftype == FEATDEF_misc_feature || featdeftype == FEATDEF_STS || featdeftype == FEATDEF_variation) {
              severity = SEV_WARNING;
            } else {
              if (isGPSorNTorNCorNGorNW || GPSorNTorNCorNGorNW (vsp->sep, sfp->location)) {
                isGPSorNTorNCorNGorNW = TRUE;
                if (! isCuratedFlybase) {
                  if (isDrosophila) {
                    isCuratedFlybase = TRUE;
                  }
                }
                if (isCuratedFlybase) {
                  /* curated fly source still has duplicate features */
                  severity = SEV_WARNING;
                }
              } else if (isGenBankAccn || IsGenBankAccn (vsp->sep, sfp->location)) {
                isGenBankAccn = TRUE;
                if (! isCuratedFlybase) {
                  if (isDrosophila) {
                    isCuratedFlybase = TRUE;
                  }
                }
                if (isCuratedFlybase) {
                  /* curated fly source still has duplicate features */
                  severity = SEV_WARNING;
                }
              } else if (isGeneralAccn || IsGeneralAccn (vsp->sep, sfp->location)) {
                isGeneralAccn = TRUE;
                if (! isCuratedFlybase) {
                  if (isDrosophila) {
                    isCuratedFlybase = TRUE;
                  }
                }
                if (isCuratedFlybase) {
                  /* curated fly source still has duplicate features */
                  severity = SEV_WARNING;
                }
              } else {
                severity = SEV_WARNING;
              }
            }
            /* if different CDS frames, lower to warning */
            if (sfp->data.choice == SEQFEAT_CDREGION) {
              crp = (CdRegionPtr) sfp->data.value.ptrvalue;
              if (crp != NULL) {
                if (frame > 1 || crp->frame > 1) {
                  if (frame != crp->frame) {
                    severity = SEV_WARNING;
                    if (! samelabel) {
                      if (numivals == 1 && left == 0 && right == bsp->length - 1 && partialL && partialR) {
                        /* skip full length partial CDS features in different frames with different products */
                        severity = SEV_NONE;
                      }
                    }
                  }
                }
              }
            }
            if (isGPSorNTorNCorNGorNW || GPSorNTorNCorNGorNW (vsp->sep, sfp->location)) {
              isGPSorNTorNCorNGorNW = TRUE;
              severity = SEV_WARNING;
            }
            if (FlybaseDbxrefs (last->dbxref) || FlybaseDbxrefs (sfp->dbxref)) {
              severity = SEV_ERROR;
            }
            if (featdeftype == FEATDEF_repeat_region) {
              severity = SEV_WARNING;
            }
            if (featdeftype == FEATDEF_SITE || featdeftype == FEATDEF_BOND) {
              severity = SEV_WARNING;
            }
            if (severity == SEV_NONE) {
              /* skip full length partial CDS features in different frames with different products */
            } else if (featdeftype == FEATDEF_REGION && DifferentDbxrefs (last->dbxref, sfp->dbxref)) {
              /* do not report if both have dbxrefs and they are different */
            } else if (featdeftype == FEATDEF_variation && ReplaceQualsDiffer (sfp->qual, last->qual)) {
              /* do not report if both have replace quals and they are different */
            } else if (CDSsLinkedToDifferentMRNAs (sfp, last)) {
              /* do not report if CDSs are linked to two different mRNAs */
            } else if (MRNAsLinkedToDifferentCDSs (sfp, last)) {
              /* do not report if mRNAs are linked to two different CDSs */
            } else if (fcontext.sap == sap) {
              if (samelabel) {
                if (GeneXrefsDifferent (sfp, last)) {
                  severity = SEV_WARNING;
                }
                ValidErr (vsp, severity, ERR_SEQ_FEAT_FeatContentDup, "Duplicate feature");
              } else if (featdeftype != FEATDEF_PUB) {
                if (fcontext.partialL != partialL || fcontext.partialR != partialR) {
                  /* do not report if partial flags are different */
                } else {
                  if (suppress_duplicate_messages && (featdeftype == FEATDEF_CDS || featdeftype == FEATDEF_mRNA) && HaveUniqueFeatIDXrefs (xref, sfp->xref)) {
                    /* do not report CDS or mRNA if every one has a unique product and unique featID xrefs */
                  } else if (featdeftype == FEATDEF_GENE &&
                             StringStr (sfp->except_text, "dicistronic gene") != NULL &&
                             StringStr (except_text, "dicistronic gene") != NULL &&
                             isCuratedFlybase) {
                    /* do not report genes marked dicistronic */
                  } else if (vsp->is_small_genome_set && SeqLocCompare (sfp->location, last->location) != SLC_A_EQ_B) {
                    /* do not report trans-spliced features that really are different on far components */
                  } else {
                    if (featdeftype == FEATDEF_GENE && isViral && (sfp->partial || last->partial)) {
                      severity = SEV_WARNING;
                    }
                    if (featdeftype == FEATDEF_CDS && isViral && (sfp->partial || last->partial)) {
                      severity = SEV_WARNING;
                    }
                    if (featdeftype == FEATDEF_mRNA && (sfp->partial || last->partial)) {
                      severity = SEV_WARNING;
                    }
                    if (featdeftype == FEATDEF_GENE && (sfp->partial || last->partial)) {
                      severity = SEV_WARNING;
                    }
                    if (featdeftype == FEATDEF_GENE && sfp->pseudo && last->pseudo) {
                      severity = SEV_WARNING;
                    }
                    if (featdeftype == FEATDEF_GENE && isViral) {
                      severity = SEV_WARNING;
                    }
                    if (severity == SEV_ERROR && featdeftype == FEATDEF_CDS && mip != NULL) {
                      if (mip->tech >= MI_TECH_htgs_1 && mip->tech <= MI_TECH_htgs_3) {
                        severity = SEV_WARNING;
                      }
                    }
                    if (fcontext.seqfeattype == SEQFEAT_IMP) {
                      severity = SEV_WARNING;
                    }
                    ValidErr (vsp, severity, ERR_SEQ_FEAT_DuplicateFeat, "Features have identical intervals, but labels differ");
                  }
                }
              }
            } else {
              if (samelabel) {
                if (GeneXrefsDifferent (sfp, last)) {
                  severity = SEV_WARNING;
                }
                ValidErr (vsp, severity, ERR_SEQ_FEAT_FeatContentDup, "Duplicate feature (packaged in different feature table)");
              } else if (featdeftype != FEATDEF_PUB) {
                if (suppress_duplicate_messages && (featdeftype == FEATDEF_CDS || featdeftype == FEATDEF_mRNA) && HaveUniqueFeatIDXrefs (xref, sfp->xref)) {
                  /* do not report CDS or mRNA if every one has a unique product and unique featID xrefs */
                } else {
                  ValidErr (vsp, /* severity */ SEV_WARNING, ERR_SEQ_FEAT_DuplicateFeat, "Features have identical intervals, but labels differ (packaged in different feature table)");
                }
              }
            }
            vsp->sfp = NULL;
            if (gcp != NULL) {
              gcp->itemID = olditemid;
              gcp->thistype = olditemtype;
            }
          }
        }
      }
      found_overlapping_peptide = FALSE;
      if (fcontext.featdeftype == FEATDEF_mat_peptide_aa ||
          fcontext.featdeftype == FEATDEF_sig_peptide_aa || fcontext.featdeftype == FEATDEF_transit_peptide_aa) {
        if (featdeftype == FEATDEF_mat_peptide_aa || featdeftype == FEATDEF_sig_peptide_aa || featdeftype == FEATDEF_transit_peptide_aa) {
          if (fcontext.left <= right && NotPeptideException (sfp, last)) {
            if (gcp != NULL) {
              gcp->itemID = fcontext.itemID;
              gcp->thistype = OBJ_SEQFEAT;
            }
            buf [0] = '\0';
            cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
            if (cds != NULL) {
              nbsp = BioseqFindFromSeqLoc (cds->location);
              if (nbsp != NULL) {
                SeqIdWrite (nbsp->id, buf, PRINTID_FASTA_LONG, sizeof (buf) - 1);
              }
            }
            vsp->descr = NULL;
            if (StringDoesHaveText (buf)) {
              if (!last_reported) {
                vsp->sfp = last;
                ValidErr (vsp, overlapPepSev, ERR_SEQ_FEAT_OverlappingPeptideFeat,
                          "Signal, Transit, or Mature peptide features overlap (parent CDS is on %s)", buf);
              }
              vsp->sfp = sfp;
              ValidErr (vsp, overlapPepSev, ERR_SEQ_FEAT_OverlappingPeptideFeat,
                        "Signal, Transit, or Mature peptide features overlap (parent CDS is on %s)", buf);
            } else {
              if (!last_reported) {
                vsp->sfp = last;
                ValidErr (vsp, overlapPepSev, ERR_SEQ_FEAT_OverlappingPeptideFeat, "Signal, Transit, or Mature peptide features overlap");
              }
              vsp->sfp = sfp;
              ValidErr (vsp, overlapPepSev, ERR_SEQ_FEAT_OverlappingPeptideFeat, "Signal, Transit, or Mature peptide features overlap");
            }
            found_overlapping_peptide = TRUE;
            vsp->sfp = NULL;
            if (gcp != NULL) {
              gcp->itemID = olditemid;
              gcp->thistype = olditemtype;
            }
          }
        }
      }
      last_reported = found_overlapping_peptide;
    }
    if (leave) {
      last = sfp;
      left = fcontext.left;
      right = fcontext.right;
      label = fcontext.label;
      comment = sfp->comment;
      strand = fcontext.strand;
      partialL = fcontext.partialL;
      partialR = fcontext.partialR;
      featdeftype = fcontext.featdeftype;
      numivals = fcontext.numivals;
      ivals = fcontext.ivals;
      sap = fcontext.sap;
      xref = sfp->xref;
      except_text = sfp->except_text;
      frame = 0;
      if (sfp->data.choice == SEQFEAT_CDREGION) {
        crp = (CdRegionPtr) sfp->data.value.ptrvalue;
        if (crp != NULL) {
          frame = crp->frame;
        }
      }
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
  }

  lastLabel = NULL;
  lastsfp = NULL;
  lastgrp = NULL;
  grp = NULL;
  sfp = SeqMgrGetNextFeatureByLabel (bsp, NULL, SEQFEAT_GENE, 0, &fcontext);
  while (sfp != NULL) {
    grp = (GeneRefPtr) sfp->data.value.ptrvalue;
    label = fcontext.label;
    if (lastLabel != NULL) {
      message = NULL;
      if (StringCmp (lastLabel, label) == 0) {
        message = "Colliding names in gene features";
      } else if (StringICmp (lastLabel, label) == 0) {
        message = "Colliding names (with different capitalization) in gene features";
      }
      if (message != NULL && (ReportGeneCollision (grp, lastgrp))) {
        if (gcp != NULL) {
          gcp->itemID = fcontext.itemID;
          gcp->thistype = OBJ_SEQFEAT;
        }
        vsp->descr = NULL;
        vsp->sfp = sfp;

        if (lastsfp != NULL && SeqLocCompare (sfp->location, lastsfp->location) == SLC_A_EQ_B) {
          ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_MultiplyAnnotatedGenes, "%s, but feature locations are identical", message);
        } else if (vsp->is_small_genome_set && StringISearch (lastsfp->except_text, "trans-splicing") != NULL && StringISearch (sfp->except_text, "trans-splicing") != NULL) {
          /* suppress for trans-spliced genes on small genome set */
        } else if (FeatureSequencesIdentical (sfp, lastsfp)) {
          ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_ReplicatedGeneSequence, "%s, but underlying sequences are identical", message);
        } else {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_CollidingGeneNames, "%s", message);
        }
        vsp->sfp = NULL;
        if (gcp != NULL) {
          gcp->itemID = olditemid;
          gcp->thistype = olditemtype;
        }
      }
    }
    lastLabel = label;
    lastsfp = sfp;
    lastgrp = grp;
    sfp = SeqMgrGetNextFeatureByLabel (bsp, sfp, SEQFEAT_GENE, 0, &fcontext);
  }

  lastLabel = NULL;
  sfp = SeqMgrGetNextGeneByLocusTag (bsp, NULL, &fcontext);
  while (sfp != NULL) {
    label = NULL;
    if (sfp->data.choice == SEQFEAT_GENE) {
      grp = (GeneRefPtr) sfp->data.value.ptrvalue;
      if (grp != NULL) {
        label = grp->locus_tag;
      }
    }
    if (lastLabel != NULL) {
      message = NULL;
      if (StringCmp (lastLabel, label) == 0) {
        message = "Colliding locus_tags in gene features";
      } else if (StringICmp (lastLabel, label) == 0) {
        message = "Colliding locus_tags (with different capitalization) in gene features";
      }
      if (message != NULL) {
        if (gcp != NULL) {
          gcp->itemID = fcontext.itemID;
          gcp->thistype = OBJ_SEQFEAT;
        }
        vsp->descr = NULL;
        vsp->sfp = sfp;
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_CollidingLocusTags, "%s", message);
        vsp->sfp = NULL;
        if (gcp != NULL) {
          gcp->itemID = olditemid;
          gcp->thistype = olditemtype;
        }
      }
    }
    lastLabel = label;
    sfp = SeqMgrGetNextGeneByLocusTag (bsp, sfp, &fcontext);
  }

  /* do UTR vs. CDS check on genomic if only one CDS, still need separate minus strand logic */
  cdscount = 0;
  genecount = 0;
  utr5count = 0;
  utr3count = 0;
  cdsgene = NULL;
  utr5gene = NULL;
  utr3gene = NULL;
  strand = Seq_strand_plus;
  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  while (sfp != NULL /* && cdscount < 2 && genecount < 2 */) {
    if (sfp->idx.subtype == FEATDEF_CDS) {
      strand = fcontext.strand;
      cdscount++;
      cdsgene = SeqMgrGetOverlappingGene (sfp->location, NULL);
    } else if (sfp->idx.subtype == FEATDEF_GENE) {
      genecount++;
    } else if (sfp->idx.subtype == FEATDEF_5UTR) {
      utr5count++;
      utr5gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
    } else if (sfp->idx.subtype == FEATDEF_3UTR) {
      utr3count++;
      utr3gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
  }
  if (bvsp->is_mrna || cdscount == 1 && genecount < 2) {
    if (bvsp->is_mrna) {
      strand = Seq_strand_plus;
    }
    fiveUTRright = 0;
    cdsRight = 0;
    threeUTRright = 0;
    firstCDS = TRUE;

    if (strand == Seq_strand_minus) {

      sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
      while (sfp != NULL) {
        if (gcp != NULL) {
          gcp->itemID = fcontext.itemID;
          gcp->thistype = OBJ_SEQFEAT;
        }
        vsp->descr = NULL;
        vsp->sfp = sfp;
        if (sfp->idx.subtype == FEATDEF_3UTR && utr3count < 2) {
          if (fcontext.strand != Seq_strand_minus) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UTRdoesNotAbutCDS, "3'UTR is not on minus strand");
          } else if (threeUTRright > 0) {
            if (threeUTRright + 1 != fcontext.left) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UTRdoesNotAbutCDS, "Previous 3'UTR does not abut next 3'UTR");
            }
          }
          threeUTRright = fcontext.right;
        } else if (sfp->idx.subtype == FEATDEF_CDS) {
          cdsRight = fcontext.right;
          if (threeUTRright > 0 && firstCDS) {
            if (threeUTRright + 1 != fcontext.left) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UTRdoesNotAbutCDS, "CDS does not abut 3'UTR");
            }
          }
          firstCDS = FALSE;
        } else if (sfp->idx.subtype == FEATDEF_5UTR && utr5count < 2) {
          if (fcontext.strand != Seq_strand_minus) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UTRdoesNotAbutCDS, "5'UTR is not on minus strand");
          } else if (cdsRight > 0) {
            if (cdsRight + 1 != fcontext.left) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UTRdoesNotAbutCDS, "5'UTR does not abut CDS");
            }
          }
          threeUTRright = fcontext.right;
        }
        sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
      }

    } else {

      sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
      while (sfp != NULL) {
        if (gcp != NULL) {
          gcp->itemID = fcontext.itemID;
          gcp->thistype = OBJ_SEQFEAT;
        }
        vsp->descr = NULL;
        vsp->sfp = sfp;
        if (sfp->idx.subtype == FEATDEF_5UTR && utr5count < 2) {
          if (fcontext.strand == Seq_strand_minus) {
            if (genecount > 1 && cdsgene != NULL && utr5gene != NULL && cdsgene != utr5gene) {
              /* ignore */
            } else {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UTRdoesNotAbutCDS, "5'UTR is not on plus strand");
            }
          }
          fiveUTRright = fcontext.right;
        } else if (sfp->idx.subtype == FEATDEF_CDS) {
          cdsRight = fcontext.right;
          if (fiveUTRright > 0 && firstCDS) {
            if (fiveUTRright + 1 != fcontext.left) {
              if (genecount > 1 && cdsgene != NULL && utr5gene != NULL && cdsgene != utr5gene) {
                /* ignore */
              } else {
                ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UTRdoesNotAbutCDS, "5'UTR does not abut CDS");
              }
            }
          }
          firstCDS = FALSE;
        } else if (sfp->idx.subtype == FEATDEF_3UTR && utr3count < 2) {
          if (fcontext.strand == Seq_strand_minus) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UTRdoesNotAbutCDS, "3'UTR is not on plus strand");
          } else if (threeUTRright > 0) {
            if (threeUTRright + 1 != fcontext.left) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UTRdoesNotAbutCDS, "Previous 3'UTR does not abut next 3'UTR");
            }
          } else if (cdsRight > 0) {
            if (cdsRight + 1 != fcontext.left) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UTRdoesNotAbutCDS, "CDS does not abut 3'UTR");
            }
          }
          if (bvsp->is_mrna && cdscount == 1 && utr3count == 1 && fcontext.right != bsp->length - 1) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UTRdoesNotExtendToEnd, "3'UTR does not extend to end of mRNA");
          }
          threeUTRright = fcontext.right;
        }
        sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
      }
    }
  }

  if (! bvsp->is_mrna) {
    last = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_RNA, 0, &fcontext);
    if (last != NULL) {
      lastrnatype = WhichRNA (last);
      left = fcontext.left;
      right = fcontext.right;
      strand = fcontext.strand;
      partialL = fcontext.partialL;
      partialR = fcontext.partialR;
      sfp = SeqMgrGetNextFeature (bsp, last, SEQFEAT_RNA, 0, &fcontext);
      while (sfp != NULL) {
        thisrnatype = WhichRNA (sfp);
        if (fcontext.strand == strand || (strand != Seq_strand_minus && fcontext.strand != Seq_strand_minus)) {
          if (lastrnatype != 0 && thisrnatype != 0) {
            if (right + 1 < fcontext.left) {
              /* gap */
              if (BaseRangeIsVirtual (bsp, right + 1, fcontext.left)) {
                /* ignore if abuts gap */
              } else if (strand == Seq_strand_minus) {
                if ((lastrnatype == RIGHT_RIBOSOMAL_SUBUNIT && (thisrnatype == INTERNAL_SPACER_2 || thisrnatype == INTERNAL_SPACER_X)) ||
                    (lastrnatype == INTERNAL_SPACER_2 && thisrnatype == MIDDLE_RIBOSOMAL_SUBUNIT) ||
                    (lastrnatype == MIDDLE_RIBOSOMAL_SUBUNIT && thisrnatype == INTERNAL_SPACER_1) ||
                    (lastrnatype == INTERNAL_SPACER_1 && thisrnatype == LEFT_RIBOSOMAL_SUBUNIT) ||
                    (lastrnatype == INTERNAL_SPACER_X && thisrnatype == LEFT_RIBOSOMAL_SUBUNIT)) {
                  if (gcp != NULL) {
                    gcp->itemID = fcontext.itemID;
                    gcp->thistype = OBJ_SEQFEAT;
                  }
                  vsp->descr = NULL;
                  vsp->sfp = sfp;
                  ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_ITSdoesNotAbutRRNA, "ITS does not abut adjacent rRNA component");
                }
              } else {
                if ((lastrnatype == LEFT_RIBOSOMAL_SUBUNIT && (thisrnatype == INTERNAL_SPACER_1 || thisrnatype == INTERNAL_SPACER_X)) ||
                    (lastrnatype == MIDDLE_RIBOSOMAL_SUBUNIT && thisrnatype == INTERNAL_SPACER_2) ||
                    (lastrnatype == INTERNAL_SPACER_1 && thisrnatype == MIDDLE_RIBOSOMAL_SUBUNIT) ||
                    (lastrnatype == INTERNAL_SPACER_2 && thisrnatype == RIGHT_RIBOSOMAL_SUBUNIT) ||
                    (lastrnatype == INTERNAL_SPACER_X && thisrnatype == RIGHT_RIBOSOMAL_SUBUNIT)) {
                  if (gcp != NULL) {
                    gcp->itemID = fcontext.itemID;
                    gcp->thistype = OBJ_SEQFEAT;
                  }
                  vsp->descr = NULL;
                  vsp->sfp = sfp;
                  ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_ITSdoesNotAbutRRNA, "ITS does not abut adjacent rRNA component");
                }
              }
            } else if (right + 1 > fcontext.left) {
              /* overlaps */
              if (strand == Seq_strand_minus) {
                if ((lastrnatype == RIGHT_RIBOSOMAL_SUBUNIT && (thisrnatype == INTERNAL_SPACER_2 || thisrnatype == INTERNAL_SPACER_X)) ||
                    (lastrnatype == INTERNAL_SPACER_2 && thisrnatype == MIDDLE_RIBOSOMAL_SUBUNIT) ||
                    (lastrnatype == MIDDLE_RIBOSOMAL_SUBUNIT && thisrnatype == INTERNAL_SPACER_1) ||
                    (lastrnatype == INTERNAL_SPACER_1 && thisrnatype == LEFT_RIBOSOMAL_SUBUNIT) ||
                    (lastrnatype == INTERNAL_SPACER_X && thisrnatype == LEFT_RIBOSOMAL_SUBUNIT)) {
                  if (gcp != NULL) {
                    gcp->itemID = fcontext.itemID;
                    gcp->thistype = OBJ_SEQFEAT;
                  }
                  vsp->descr = NULL;
                  vsp->sfp = sfp;
                  ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadRRNAcomponentOverlap, "ITS overlaps adjacent rRNA component");
                }
              } else {
                if ((lastrnatype == LEFT_RIBOSOMAL_SUBUNIT && (thisrnatype == INTERNAL_SPACER_1 || thisrnatype == INTERNAL_SPACER_X)) ||
                    (lastrnatype == MIDDLE_RIBOSOMAL_SUBUNIT && thisrnatype == INTERNAL_SPACER_2) ||
                    (lastrnatype == INTERNAL_SPACER_1 && thisrnatype == MIDDLE_RIBOSOMAL_SUBUNIT) ||
                    (lastrnatype == INTERNAL_SPACER_2 && thisrnatype == RIGHT_RIBOSOMAL_SUBUNIT) ||
                    (lastrnatype == INTERNAL_SPACER_X && thisrnatype == RIGHT_RIBOSOMAL_SUBUNIT)) {
                  if (gcp != NULL) {
                    gcp->itemID = fcontext.itemID;
                    gcp->thistype = OBJ_SEQFEAT;
                  }
                  vsp->descr = NULL;
                  vsp->sfp = sfp;
                  ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadRRNAcomponentOverlap, "ITS overlaps adjacent rRNA component");
                }
              }
            } else {
              /* abuts */
              if (strand == Seq_strand_minus) {
                if (lastrnatype == thisrnatype && partialL && fcontext.partialR && bsp->repr == Seq_repr_seg) {
                  /* okay in segmented set */
                } else if ((lastrnatype == RIGHT_RIBOSOMAL_SUBUNIT && (thisrnatype != INTERNAL_SPACER_2 && thisrnatype != INTERNAL_SPACER_X)) ||
                    (lastrnatype == INTERNAL_SPACER_2 && thisrnatype != MIDDLE_RIBOSOMAL_SUBUNIT) ||
                    (lastrnatype == MIDDLE_RIBOSOMAL_SUBUNIT && thisrnatype != INTERNAL_SPACER_1) ||
                    (lastrnatype == INTERNAL_SPACER_1 && thisrnatype != LEFT_RIBOSOMAL_SUBUNIT) ||
                    (lastrnatype == INTERNAL_SPACER_X && thisrnatype != LEFT_RIBOSOMAL_SUBUNIT)) {
                  if (gcp != NULL) {
                    gcp->itemID = fcontext.itemID;
                    gcp->thistype = OBJ_SEQFEAT;
                  }
                  vsp->descr = NULL;
                  vsp->sfp = sfp;
                  ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadRRNAcomponentOrder, "Problem with order of abutting rRNA components");
                }
              } else {
                if (lastrnatype == thisrnatype && partialR && fcontext.partialL && bsp->repr == Seq_repr_seg) {
                  /* okay in segmented set */
                } else if ((lastrnatype == LEFT_RIBOSOMAL_SUBUNIT && (thisrnatype != INTERNAL_SPACER_1 && thisrnatype != INTERNAL_SPACER_X)) ||
                    (lastrnatype == MIDDLE_RIBOSOMAL_SUBUNIT && thisrnatype != INTERNAL_SPACER_2) ||
                    (lastrnatype == INTERNAL_SPACER_1 && thisrnatype != MIDDLE_RIBOSOMAL_SUBUNIT) ||
                    (lastrnatype == INTERNAL_SPACER_2 && thisrnatype != RIGHT_RIBOSOMAL_SUBUNIT) ||
                    (lastrnatype == INTERNAL_SPACER_X && thisrnatype != RIGHT_RIBOSOMAL_SUBUNIT)) {
                  if (gcp != NULL) {
                    gcp->itemID = fcontext.itemID;
                    gcp->thistype = OBJ_SEQFEAT;
                  }
                  vsp->descr = NULL;
                  vsp->sfp = sfp;
                  ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadRRNAcomponentOrder, "Problem with order of abutting rRNA components");
                }
              }
            }
          }
        }
        last = sfp;
        left = fcontext.left;
        right = fcontext.right;
        strand = fcontext.strand;
        partialL = fcontext.partialL;
        partialR = fcontext.partialR;
        lastrnatype = thisrnatype;
        sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_RNA, 0, &fcontext);
      }
    }
  }

  vsp->sfp = NULL;
  if (gcp != NULL) {
    gcp->itemID = olditemid;
    gcp->thistype = olditemtype;
  }

  CheckForGenesInconsistent(bsp, vsp);

  if (ISA_na (bsp->mol)) {
    sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_gap, &fcontext);
    while (sfp != NULL) {
      estimated_length = 0;
      for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
        if (StringICmp (gbq->qual, "estimated_length") != 0) continue;
        if (StringHasNoText (gbq->val)) continue;
        if (StringICmp (gbq->val, "unknown") == 0) continue;
        if (sscanf (gbq->val, "%ld", &val) == 1) {
          estimated_length = val;
        }
      }
      if (StreamCacheSetup (NULL, sfp->location, EXPAND_GAPS_TO_DASHES, &sc)) {
        dashes = 0;
        Ns = 0;
        realBases = 0;
        while ((res = StreamCacheGetResidue (&sc)) != '\0') {
          if (IS_LOWER (res)) {
            res = TO_UPPER (res);
          }
          if (res == '-') {
            dashes++;
          } else if (res == 'N') {
            Ns++;
          } else {
            realBases++;
          }
        }
        if (gcp != NULL) {
          gcp->itemID = fcontext.itemID;
          gcp->thistype = OBJ_SEQFEAT;
        }
        vsp->descr = NULL;
        vsp->sfp = sfp;
        loclen = SeqLocLen (sfp->location);
        if (estimated_length > 0 && estimated_length != loclen) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_GapFeatureProblem, "Gap feature estimated_length %ld does not match %ld feature length",
                    (long) estimated_length, (long) loclen);
        } else if (realBases > 0 && Ns > 0) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_GapFeatureProblem, "Gap feature over %ld real bases and %ld Ns", (long) realBases, (long) Ns);
        } else if (realBases > 0) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_GapFeatureProblem, "Gap feature over %ld real bases", (long) realBases);
        } else if (Ns > 0) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_GapFeatureProblem, "Gap feature over %ld Ns", (long) Ns);
        } else if (estimated_length > 0 && dashes != estimated_length) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_GapFeatureProblem, "Gap feature estimated_length %ld does not match %ld gap characters",
                    (long) estimated_length, (long) dashes);
        }
      }
      sfp = SeqMgrGetNextFeature (bsp, sfp, 0, FEATDEF_gap, &fcontext);
    }
  }
  if (gcp != NULL) {
    gcp->itemID = olditemid;
    gcp->thistype = olditemtype;
  }
  vsp->descr = NULL;
  vsp->sfp = NULL;

  if (ISA_aa (bsp->mol)) {
    sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
    while (sfp != NULL) {
      slp = SeqLocFindNext (sfp->location, NULL);
      while (slp != NULL) {
        if (SeqLocStrand (slp) == Seq_strand_minus) {
          if (gcp != NULL) {
            gcp->itemID = fcontext.itemID;
            gcp->thistype = OBJ_SEQFEAT;
          }
          vsp->descr = NULL;
          vsp->sfp = sfp;
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_MinusStrandProtein, "Feature on protein indicates negative strand");
        }
        slp = SeqLocFindNext (sfp->location, slp);
      }
      sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
    }
  }
  if (gcp != NULL) {
    gcp->itemID = olditemid;
    gcp->thistype = olditemtype;
  }
  vsp->descr = NULL;
  vsp->sfp = NULL;

  lastbiop = NULL;
  lastsfp = NULL;
  numBadFullSource = 0;
  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_BIOSRC, 0, &fcontext);
  if (sfp != NULL) {
    if (fcontext.left == 0 && fcontext.right == bsp->length - 1 && fcontext.numivals == 1) {
      showBadFullSource = TRUE;
      sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
      if (sdp != NULL) {
        biop = (BioSourcePtr) sdp->data.ptrvalue;
        if (biop != NULL) {
          if (biop->is_focus) {
            showBadFullSource = FALSE;
          }
          for (sbsp = biop->subtype; sbsp != NULL; sbsp = sbsp->next) {
            if (sbsp->subtype == SUBSRC_transgenic) {
              showBadFullSource = FALSE;
            }
          }
        }
      }
      if (showBadFullSource) {
        if (gcp != NULL) {
          gcp->itemID = fcontext.itemID;
          gcp->thistype = OBJ_SEQFEAT;
        }
        vsp->descr = NULL;
        vsp->sfp = sfp;
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadFullLengthFeature, "Source feature is full length, should be descriptor");
        vsp->sfp = NULL;
        if (gcp != NULL) {
          gcp->itemID = olditemid;
          gcp->thistype = olditemtype;
        }
      }
    }
  }
  /* and fall through to continue testing first and remaining source features */
  while (sfp != NULL) {
    if (fcontext.left == 0 && fcontext.right == bsp->length - 1 && fcontext.numivals == 1) {
      numBadFullSource++;
      if (numBadFullSource > 1) {
        if (gcp != NULL) {
          gcp->itemID = fcontext.itemID;
          gcp->thistype = OBJ_SEQFEAT;
        }
        vsp->descr = NULL;
        vsp->sfp = sfp;
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadFullLengthFeature, "Multiple full-length source features, should only be one if descriptor is transgenic");
        vsp->sfp = NULL;
        if (gcp != NULL) {
          gcp->itemID = olditemid;
          gcp->thistype = olditemtype;
        }
      }
    }
    biop = (BioSourcePtr) sfp->data.value.ptrvalue;
    if (biop != NULL && lastbiop != NULL) {
      if (lastsfp != NULL) {
        if (StringDoesHaveText (lastsfp->comment) && StringDoesHaveText (sfp->comment) && StringICmp (lastsfp->comment, sfp->comment) != 0) {
          /* different comments, so ignore */
        } else if (IsIdenticalBioSource (biop, lastbiop) && (! bvsp->is_synthetic) && (!bvsp->is_artificial)) {
          if (gcp != NULL) {
            gcp->itemID = fcontext.itemID;
            gcp->thistype = OBJ_SEQFEAT;
          }
          vsp->descr = NULL;
          vsp->sfp = sfp;
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_MultipleEquivBioSources, "Multiple equivalent source features should be combined into one multi-interval feature");
          vsp->sfp = NULL;
          if (gcp != NULL) {
            gcp->itemID = olditemid;
            gcp->thistype = olditemtype;
          }
        }
      }
    }
    lastbiop = biop;
    lastsfp = sfp;
    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_BIOSRC, 0, &fcontext);
  }

  lastpdp = NULL;
  lastsfp = NULL;
  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_PUB, 0, &fcontext);
  if (sfp != NULL) {
    if (fcontext.left == 0 && fcontext.right == bsp->length - 1 && fcontext.numivals == 1) {
      if (gcp != NULL) {
        gcp->itemID = fcontext.itemID;
        gcp->thistype = OBJ_SEQFEAT;
      }
      vsp->descr = NULL;
      vsp->sfp = sfp;
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadFullLengthFeature, "Publication feature is full length, should be descriptor");
      vsp->sfp = NULL;
      if (gcp != NULL) {
        gcp->itemID = olditemid;
        gcp->thistype = olditemtype;
      }
    }
  }
  while (sfp != NULL) {
    pdp = (PubdescPtr) sfp->data.value.ptrvalue;
    if (pdp != NULL && lastpdp != NULL) {
      if (lastsfp != NULL) {
        if (StringDoesHaveText (lastsfp->comment) && StringDoesHaveText (sfp->comment) && StringICmp (lastsfp->comment, sfp->comment) != 0) {
          /* different comments, so ignore */
        } else if (IsIdenticalPublication (pdp, lastpdp)) {
          if (gcp != NULL) {
            gcp->itemID = fcontext.itemID;
            gcp->thistype = OBJ_SEQFEAT;
          }
          vsp->descr = NULL;
          vsp->sfp = sfp;
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_MultipleEquivPublications, "Multiple equivalent publication features should be combined into one multi-interval feature");
          vsp->sfp = NULL;
          if (gcp != NULL) {
            gcp->itemID = olditemid;
            gcp->thistype = olditemtype;
          }
        }
      }
    }
    lastpdp = pdp;
    lastsfp = sfp;
    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_PUB, 0, &fcontext);
  }

  SeqMgrExploreDescriptors (bsp, (Pointer) bvsp, ValidateSeqDescrIndexed, NULL);

  omdp = ObjMgrGetData (gcp->entityID);
  if (omdp != NULL && omdp->datatype == OBJ_SEQSUB) {
    ssp = (SeqSubmitPtr) omdp->dataptr;
    if (ssp != NULL) {
      sbp = ssp->sub;
      if (sbp != NULL) {
        bvsp->got_a_pub = TRUE;
      }
    }
  }

  ValidateCDSmRNAmatch (vsp, bsp, numgene, numcds, nummrna, suppress_duplicate_messages);

  if (vsp->locusTagGeneralMatch) {
    ValidateLocusTagGeneral (vsp, bsp);
  }

  if (ISA_na (bsp->mol) && SeqMgrGetParentOfPart (bsp, NULL) == NULL) {
    LookForMultipleUnpubPubs (vsp, gcp, bsp);
  }

  if (bsp->repr == Seq_repr_delta && ISA_na (bsp->mol)) {
    CheckBioseqForFeatsInGap (bsp, vsp);
  }

  CheckForNonViralComplete (bsp, vsp, gcp);

  LookForViralMolInfoInLineage (bsp, vsp, gcp);

  return TRUE;
}

static Boolean ValidateBioseqContextGather (GatherContextPtr gcp)
{
  ValidStructPtr  vsp;
  BioseqValidStrPtr bvsp;
  CitSubPtr       csp;

  bvsp = (BioseqValidStrPtr) (gcp->userdata);
  vsp = bvsp->vsp;
  vsp->descr = NULL;
  vsp->sfp = NULL;
  vsp->gcp = gcp;               /* needed for ValidErr */

  switch (gcp->thistype) {
  case OBJ_SEQFEAT:
    ValidateSeqFeatContext (gcp);
    break;
  case OBJ_SEQDESC:
    ValidateSeqDescrContext (gcp);
    break;
  case OBJ_SEQSUB_CIT:
    bvsp->got_a_pub = TRUE;
    csp = (CitSubPtr) gcp->thisitem;
    ValidateCitSub (vsp, csp);
    break;
  default:
    break;
  }
  return TRUE;
}


static ValNodePtr ListFeaturesContainedInLocation (BioseqPtr bsp, SeqLocPtr slp, Uint1 seqfeatChoice, Uint1 featdefChoice)
{
  ValNodePtr        feat_list = NULL;
  SeqMgrFeatContext fcontext;
  SeqFeatPtr        sfp;
  Int4              loc_left, loc_right, tmp;
  Int4              cmp;

  if (bsp == NULL || slp == NULL) return NULL;

  loc_left = SeqLocStart (slp);
  loc_right = SeqLocStop (slp);
  if (loc_left > loc_right) {
    tmp = loc_left;
    loc_left = loc_right;
    loc_right = tmp;
  }
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, seqfeatChoice, featdefChoice, &fcontext);
       sfp != NULL && fcontext.left <= loc_right;
       sfp = SeqMgrGetNextFeature (bsp, sfp, seqfeatChoice, featdefChoice, &fcontext))
  {
    cmp = SeqLocCompare (sfp->location, slp);
    if (cmp == SLC_A_EQ_B || cmp == SLC_A_IN_B)
    {
      ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, sfp);
    }
  }
  return feat_list;
}


typedef struct multigeneoverlap {
  SeqFeatPtr gene;
  Int4       left;
  Int4       right;
  Boolean    reported;
} MultiGeneOverlapData, PNTR MultiGeneOverlapPtr;

static MultiGeneOverlapPtr MultiGeneOverlapNew (SeqFeatPtr gene, Int4 left, Int4 right)
{
  MultiGeneOverlapPtr m;

  m = (MultiGeneOverlapPtr) MemNew (sizeof (MultiGeneOverlapData));
  m->gene = gene;
  m->left = left;
  m->right = right;
  m->reported = FALSE;
  return m;
}

static MultiGeneOverlapPtr MultiGeneOverlapFree (MultiGeneOverlapPtr m)
{
  m = MemFree (m);
  return m;
}


static void ReportContainedGenes (SeqFeatPtr gene, ValNodePtr contained_list, ValidStructPtr vsp)
{
  ValNodePtr vnp;
  Int4       cmp, num_overlap = 0;
  MultiGeneOverlapPtr g;

  if (gene == NULL || contained_list == NULL || contained_list->next == NULL) {
    return;
  }
  for (vnp = contained_list; vnp != NULL; vnp = vnp->next) {
    g = (MultiGeneOverlapPtr) vnp->data.ptrvalue;
    cmp = SeqLocCompare (gene->location, g->gene->location);
    if (cmp == SLC_A_EQ_B || cmp == SLC_B_IN_A) {
      num_overlap++;
    }
  }
  if (num_overlap > 1) {
    vsp->descr = NULL;
    vsp->sfp = gene;

    vsp->gcp->entityID = gene->idx.entityID;
    vsp->gcp->itemID = gene->idx.itemID;
    vsp->gcp->thistype = OBJ_SEQFEAT;
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_MultipleGeneOverlap, "Gene contains %d other genes", num_overlap);
  }
}


static void AddMultiGeneOverlapToList (ValNodePtr PNTR list, MultiGeneOverlapPtr m, ValidStructPtr vsp)
{
  ValNodePtr vnp;
  MultiGeneOverlapPtr g;

  if (list == NULL || m == NULL || vsp == NULL) {
    return;
  }

  if (*list == NULL) {
    ValNodeAddPointer (list, 0, m);
    return;
  }

  /* we're examining the interval for the first gene on the list */
  g = (MultiGeneOverlapPtr) (*list)->data.ptrvalue;

  /* if the new gene is outside the interval for the first gene, it's time to remove the first gene
   */
  while (g != NULL && g->right < m->left) {
    /* first, figure out if this gene contains other genes */
    ReportContainedGenes (g->gene, (*list)->next, vsp);
    /* now remove it from the list */
    vnp = *list;
    *list = (*list)->next;
    vnp->next = NULL;
    vnp = ValNodeFreeData (vnp);
    if ((*list) == NULL) {
      g = NULL;
    } else {
      g = (*list)->data.ptrvalue;
    }
  }

  /* now add this gene to list */
  ValNodeAddPointer (list, 0, m);
}


static void NewFindMultiGeneOverlaps (BioseqPtr bsp, ValidStructPtr vsp)
{
  SeqMgrFeatContext context;
  SeqFeatPtr        gene;
  ValNodePtr        gene_list = NULL, vnp;
  MultiGeneOverlapPtr m;

  vsp->bsp = bsp;

  for (gene = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_GENE, 0, &context);
       gene != NULL;
       gene = SeqMgrGetNextFeature (bsp, gene, SEQFEAT_GENE, 0, &context)) {
    m = MultiGeneOverlapNew (gene, context.left, context.right);
    AddMultiGeneOverlapToList (&gene_list, m, vsp);
  }
  for (vnp = gene_list; vnp != NULL && vnp->next != NULL; vnp = vnp->next) {
    m = (MultiGeneOverlapPtr) vnp->data.ptrvalue;
    ReportContainedGenes (m->gene, vnp->next, vsp);
  }
  gene_list = ValNodeFreeData (gene_list);
}


/*****************************************************************************
*   FindMultiGeneOverlaps (BioseqPtr bsp, ValidStructPtr vsp)
*
*      This function reports genes that overlap two or more other genes.
*****************************************************************************/
static void FindMultiGeneOverlaps (BioseqPtr bsp, ValidStructPtr vsp)
{
  GatherContextPtr  gcp;
  Uint2           oldEntityID, oldItemtype;
  Uint4           oldItemID;

  if (bsp == NULL || vsp == NULL || vsp->gcp == NULL) {
    return;
  }

  gcp = vsp->gcp;

  oldEntityID = gcp->entityID;
  oldItemID = gcp->itemID;
  oldItemtype = gcp->thistype;

  NewFindMultiGeneOverlaps (bsp, vsp);

  gcp->entityID = oldEntityID;
  gcp->itemID = oldItemID;
  gcp->thistype = oldItemtype;

}

static void ValidateTSASequenceForNs (BioseqPtr bsp, ValidStructPtr vsp)
{
  Int4              total = 0, max_stretch = 0;
  GatherContextPtr  gcp;
  ErrSev            logsev;
  ErrSev            msgsev;
  Uint2             oldEntityID, oldItemtype;
  Uint4             oldItemID;
  Int4              percent_N, allowed_percentN = 10;
  SeqFeat           sf;
  SeqInt            si;
  CharPtr           str;
  ValNode           vn;

  if (ISA_aa (bsp->mol)) {
    return;
  }
  gcp = vsp->gcp;

  oldEntityID = gcp->entityID;
  oldItemID = gcp->itemID;
  oldItemtype = gcp->thistype;

  if (IsTSA (bsp)) {
    msgsev = ErrSetMessageLevel (SEV_MAX);
    logsev = ErrSetLogLevel (SEV_MAX);

    CountNsInSequence (bsp, &total, &max_stretch, FALSE);

    ErrSetLogLevel (logsev);
    ErrSetMessageLevel (msgsev);

    percent_N = (total * 100) / bsp->length;
    if (percent_N > allowed_percentN) {
      vsp->bsp = bsp;
      vsp->descr = NULL;
      vsp->sfp = NULL;
      gcp->entityID = bsp->idx.entityID;
      gcp->itemID = bsp->idx.itemID;
      gcp->thistype = OBJ_BIOSEQ;

      ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_HighNContentPercent, "Sequence contains %d percent Ns", percent_N);
    }
    if (max_stretch > 15) {
      vsp->bsp = bsp;
      vsp->descr = NULL;
      vsp->sfp = NULL;
      gcp->entityID = bsp->idx.entityID;
      gcp->itemID = bsp->idx.itemID;
      gcp->thistype = OBJ_BIOSEQ;

      ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_HighNContentStretch, "Sequence has a stretch of %d Ns", max_stretch);
    } else if (bsp->length > 10) {
      vsp->bsp = bsp;
      vsp->descr = NULL;
      vsp->sfp = NULL;
      gcp->entityID = bsp->idx.entityID;
      gcp->itemID = bsp->idx.itemID;
      gcp->thistype = OBJ_BIOSEQ;

      MemSet ((Pointer) &sf, 0, sizeof (SeqFeat));
      MemSet ((Pointer) &si, 0, sizeof (SeqInt));
      MemSet ((Pointer) &vn, 0, sizeof (ValNode));
      sf.location = &vn;
      vn.choice = SEQLOC_INT;
      vn.data.ptrvalue = (Pointer) &si;
      si.id = bsp->id;
      si.from = 0;
      si.to = 9;
      str = GetSequenceByFeature (&sf);
      if (StringStr (str, "NNNNN") != NULL) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_HighNContentStretch, "Sequence has a stretch of at least 5 Ns within the first 10 bases");
      }
      MemFree (str);
      si.from = bsp->length - 10;
      si.to = bsp->length - 1;
      str = GetSequenceByFeature (&sf);
      if (StringStr (str, "NNNNN") != NULL) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_HighNContentStretch, "Sequence has a stretch of at least 5 Ns within the last 10 bases");
      }
      MemFree (str);
    }
  } else {
    msgsev = ErrSetMessageLevel (SEV_MAX);
    logsev = ErrSetLogLevel (SEV_MAX);

    CountNsInSequence (bsp, &total, &max_stretch, FALSE);

    ErrSetLogLevel (logsev);
    ErrSetMessageLevel (msgsev);

    percent_N = (total * 100) / bsp->length;
    if (percent_N > 50) {
      vsp->bsp = bsp;
      vsp->descr = NULL;
      vsp->sfp = NULL;
      gcp->entityID = bsp->idx.entityID;
      gcp->itemID = bsp->idx.itemID;
      gcp->thistype = OBJ_BIOSEQ;

      ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_HighNContentPercent, "Sequence contains %d percent Ns", percent_N);
    }
  }
  gcp->entityID = oldEntityID;
  gcp->itemID = oldItemID;
  gcp->thistype = oldItemtype;
}


static void ValidateRefSeqTitle (BioseqPtr bsp, ValidStructPtr vsp, Boolean is_virus)
{
  SeqDescrPtr       sdp;
  SeqMgrDescContext dcontext;
  SeqMgrFeatContext fcontext;
  CharPtr           taxname = NULL, title;
  BioSourcePtr      biop;
  SeqFeatPtr        cds, src_feat;
  size_t            len, tlen;

  if (bsp == NULL || vsp == NULL) {
    return;
  }

  if (is_virus) return;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp != NULL
      && (biop = (BioSourcePtr) sdp->data.ptrvalue) != NULL
      && biop->org != NULL) {
    taxname = biop->org->taxname;
  }
  if (ISA_aa(bsp->mol)) {
    cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
    if (cds != NULL) {
      src_feat = SeqMgrGetOverlappingSource (cds->location, &fcontext);
      if (src_feat != NULL
          && (biop = (BioSourcePtr) src_feat->data.value.ptrvalue) != NULL
          && biop->org != NULL) {
        taxname = biop->org->taxname;
      }
    }
  }

  if (StringDoesHaveText (taxname)) {
    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_title, &dcontext);
    if (sdp != NULL) {
      title = (CharPtr) sdp->data.ptrvalue;
      if (StringDoesHaveText (title)) {
        if (StringNCmp (title, "PREDICTED: ", 11) == 0) {
          title += 11;
        }
        len = StringLen (taxname);
        tlen = StringLen (title);
        if (ISA_na (bsp->mol)) {
          if (tlen < len || StringNICmp (title, taxname, len) != 0) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_NoOrganismInTitle, "RefSeq nucleotide title does not start with organism name");
          }
        } else if (ISA_aa (bsp->mol)) {
          if (tlen < len + 3 ||
              StringNICmp (title + tlen - len - 1, taxname, len) != 0 ||
              title [tlen - len - 2] != '[' ||
              title [tlen - 1] != ']') {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_NoOrganismInTitle, "RefSeq protein title does not end with organism name");
          }
        }
      }
    }
  }
}


static Boolean EndsWithSuffixPlusFieldValue (CharPtr str, CharPtr suffix, CharPtr val)
{
  CharPtr cp, last_word;

  cp = StringSearch (str, suffix);
  if (cp == NULL) {
    return FALSE;
  }
  last_word = StringRChr (str, ' ');
  if (last_word == NULL || last_word < cp) {
    return FALSE;
  }
  if (StringCmp (last_word + 1, val) == 0) {
    return TRUE;
  } else {
    return FALSE;
  }

}


static void ValidateBarcodeIndexNumber (CharPtr bin, BioseqPtr bsp, ValidStructPtr vsp)
{
  SeqDescPtr        sdp;
  SeqMgrDescContext context;
  BioSourcePtr      biop;
  Int4              bin_len;

  if (StringHasNoText (bin) || bsp == NULL || vsp == NULL) {
    return;
  }

  bin_len = StringLen (bin);

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
  if (sdp == NULL || (biop = (BioSourcePtr) sdp->data.ptrvalue) == NULL || biop->org == NULL) {
    return;
  }
  /* only check if name contains "sp." or "bacterium" */
  if (StringISearch (biop->org->taxname, "sp.") == NULL && StringISearch (biop->org->taxname, "bacterium") == NULL) {
    return;
  }
  /* only check if name contains BOLD */
  if (StringSearch (biop->org->taxname, "BOLD") == NULL) {
    return;
  }
  if (!EndsWithSuffixPlusFieldValue(biop->org->taxname, "sp. ", bin)
      && !EndsWithSuffixPlusFieldValue(biop->org->taxname, "bacterium ", bin)) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_BadStrucCommInvalidFieldValue, "Organism name should end with sp. plus Barcode Index Number (%s)", bin);
  }
}


static void ValidateStructuredCommentsInContext (BioseqPtr bsp, ValidStructPtr vsp)
{
  SeqDescPtr    sdp;
  SeqMgrDescContext dcontext;
  UserObjectPtr uop;
  ObjectIdPtr   oip;
  UserFieldPtr  curr;

  /* validate structured comments in context */
  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &dcontext);
       sdp != NULL;
       sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &dcontext))
  {
    uop = sdp->data.ptrvalue;
    if (uop != NULL && uop->type != NULL && StringICmp (uop->type->str, "StructuredComment") == 0) 
    {
      for (curr = uop->data; curr != NULL; curr = curr->next) 
      {
        if (curr->choice != 1) continue;
        oip = curr->label;
        if (oip == NULL || StringCmp (oip->str, "Barcode Index Number") != 0) continue;
        ValidateBarcodeIndexNumber ((CharPtr) curr->data.ptrvalue, bsp, vsp);
      }
    }
  }
}


/*****************************************************************************
*
*   ValidateBioseqContext(gcp)
*      Validate one Bioseq for descriptors, features, and context
*      This is done as a second Gather, focussed on the Bioseq in
*        question.
*
*****************************************************************************/
static void ValidateBioseqContext (GatherContextPtr gcp)
{
  size_t          acclen;
  ValidStructPtr  vsp;
  BioseqPtr       bsp;
  GatherScope     gs;
  BioseqValidStr  bvs;
  SeqFeatPtr      sfp;
  ValNode         fake_whole;
  SeqIdPtr        sip;
  ValNodePtr      vnp = NULL;
  MolInfoPtr      mip = NULL;
  SeqMgrDescContext dcontext;
  SeqMgrFeatContext fcontext;
  BioseqContextPtr bcp;
  Uint2           oldEntityID, oldItemtype;
  Uint4           oldItemID;
  Uint2           mipEntityID = 0, mipItemtype = 0;
  Uint4           mipItemID = 0;
  ObjMgrDataPtr   omdp;
  BioseqPtr       parent;
  PatentSeqIdPtr  psip;
  IdPatPtr        ipp;
  Boolean         isPDB = FALSE;
  Boolean         is_wgs = FALSE;
  Boolean         is_gb = FALSE;
  Boolean         is_ac = FALSE;
  Boolean         is_ch_or_cm = FALSE;
  Boolean         is_nc = FALSE;
  Boolean         is_local = FALSE;
  Boolean         is_local_only = TRUE;
  Boolean         is_organelle = FALSE;
  Boolean         is_plasmid = FALSE;
  Boolean         is_prokaryote = FALSE;
  Boolean         is_refseq = FALSE;
  Boolean         is_neg_strand_virus = FALSE;
  Boolean         is_ambisense_virus = FALSE;
  Boolean         is_synthetic = FALSE;
  Boolean         is_transgenic = FALSE;
  Boolean         is_virus = FALSE;
  Boolean         has_cds = FALSE;
  Boolean         has_chromosome = FALSE;
  ErrSev          sev;
  SubSourcePtr    ssp;
  CharPtr         str;
  CharPtr         taxname = NULL;
  TextSeqIdPtr    tsip;
  BioSourcePtr    biop = NULL;
  OrgRefPtr       orp;
  OrgNamePtr      onp;
  OrgModPtr       omp;
  /*
  Char            buf1[255];
  */

  vsp = (ValidStructPtr) (gcp->userdata);
  bsp = (BioseqPtr) (gcp->thisitem);
  vsp->bsp = bsp;
  vsp->descr = NULL;
  vsp->sfp = NULL;
  vsp->bssp = (BioseqSetPtr) (gcp->parentitem);

  MemSet (&gs, 0, sizeof (GatherScope));
  fake_whole.choice = SEQLOC_WHOLE;
  sip = SeqIdFindBest (bsp->id, 0);

  fake_whole.data.ptrvalue = sip;

  fake_whole.next = NULL;
  gs.target = &fake_whole;
  gs.get_feats_location = TRUE;
  gs.nointervals = TRUE;
  MemSet ((Pointer) (gs.ignore), (int) TRUE, (size_t) (sizeof (Boolean) * OBJ_MAX));
  gs.ignore[OBJ_SEQDESC] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  gs.ignore[OBJ_SUBMIT_BLOCK] = FALSE;
  gs.ignore[OBJ_SEQSUB_CIT] = FALSE;

  gs.scope = vsp->sep;

  MemSet (&bvs, 0, sizeof (BioseqValidStr));
  bvs.vsp = vsp;

  /* now looking for molinfo on every bioseq (okay on segset) */
  if (bsp != NULL) {
    vnp = NULL;
    if (vsp->useSeqMgrIndexes) {
      vnp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
      if (vnp != NULL) {
        mipEntityID = dcontext.entityID;
        mipItemID = dcontext.itemID;
        mipItemtype = OBJ_SEQDESC;
      }
    } else {
      bcp = BioseqContextNew (bsp);
      vnp = BioseqContextGetSeqDescr (bcp, Seq_descr_molinfo, NULL, NULL);
      BioseqContextFree (bcp);
    }
    if (vnp != NULL) {
      mip = (MolInfoPtr) vnp->data.ptrvalue;
    }

    vnp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
    if (vnp != NULL) {
      biop = (BioSourcePtr) vnp->data.ptrvalue;
      if (biop != NULL) {
        orp = biop->org;
        if (orp != NULL) {
          taxname = orp->taxname;
          if (StringICmp (orp->taxname, "Human immunodeficiency virus") == 0 ||
              StringICmp (orp->taxname, "Human immunodeficiency virus 1") == 0 ||
              StringICmp (orp->taxname, "Human immunodeficiency virus 2") == 0) {
            ValidateLocationForHIV (vsp, biop, bsp);
          } else if (StringICmp (orp->taxname, "synthetic construct") == 0) {
            bvs.is_syn_constr = TRUE;
          } else if (StringISearch (orp->taxname, "vector") != NULL) {
            bvs.is_syn_constr = TRUE;
          }
          onp = orp->orgname;
          if (onp != NULL) {
            if (StringICmp (onp->div, "SYN") == 0) {
              bvs.is_syn_constr = TRUE;
              is_synthetic = TRUE;
            }
            if (StringISearch (onp->lineage, "artificial sequences") != NULL) {
              bvs.is_syn_constr = TRUE;
            }
            if (StringISearch (onp->lineage, "negative-strand viruses") != NULL) {
              is_neg_strand_virus = TRUE;
            }
            if (StringISearch (onp->lineage, "Arenavirus") != NULL ||
                StringISearch (onp->lineage, "Arenaviridae") != NULL ||
                StringISearch (onp->lineage, "Phlebovirus") != NULL ||
                StringISearch (onp->lineage, "Tospovirus") != NULL ||
                StringISearch (onp->lineage, "Tenuivirus") != NULL) {
              is_ambisense_virus = TRUE;
            }
            if (StringNICmp (onp->lineage, "Viruses; ", 9) == 0 ||
                StringNICmp (onp->lineage, "Bacteria; ", 10) == 0 ||
                StringNICmp (onp->lineage, "Archaea; ", 9) == 0 ||
                StringCmp (onp->div, "BCT") == 0 ||
                StringCmp (onp->div, "VRL") == 0) {
              is_prokaryote = TRUE;
            }
            if (StringNICmp (onp->lineage, "Viruses; ", 9) == 0) {
              is_virus = TRUE;
            }
            is_organelle = IsLocationOrganelle (biop->genome);
            is_plasmid = (Boolean) (biop->genome == GENOME_plasmid);
            for (omp = onp->mod; omp != NULL; omp = omp->next) {
              if (omp->subtype == ORGMOD_other) {
                if (mip != NULL && (StringICmp (omp->subname, "cRNA") == 0)) {
                  oldEntityID = gcp->entityID;
                  oldItemID = gcp->itemID;
                  oldItemtype = gcp->thistype;

                  gcp->entityID = dcontext.entityID;
                  gcp->itemID = dcontext.itemID;
                  gcp->thistype = OBJ_SEQDESC;

                  if (mip->biomol == MOLECULE_TYPE_CRNA) {
                    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "cRNA note redundant with molecule type");
                  } else {
                    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "cRNA note conflicts with molecule type");
                  }

                  gcp->entityID = oldEntityID;
                  gcp->itemID = oldItemID;
                  gcp->thistype = oldItemtype;
                }
              }
            }
            if (mip != NULL && mip->biomol == MOLECULE_TYPE_GENOMIC && bsp->mol == MOLECULE_CLASS_DNA) {
              if (StringNICmp (onp->lineage, "Viruses; ", 9) == 0 && StringISearch (onp->lineage, "no DNA stage") != NULL) {
                ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Genomic DNA viral lineage indicates no DNA stage");
              }
            }
          }
        }
        if (biop->origin == ORG_ARTIFICIAL) {
          bvs.is_artificial = TRUE;
        } else if (biop->origin == ORG_SYNTHETIC) {
          bvs.is_synthetic = TRUE;
        }
        if (biop->origin == ORG_MUT || biop->origin == ORG_ARTIFICIAL || biop->origin == ORG_SYNTHETIC) {
          is_synthetic = TRUE;
        }
        for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
          if (ssp->subtype == SUBSRC_transgenic) {
            is_transgenic = TRUE;
          } else if (ssp->subtype == SUBSRC_chromosome) {
            if (StringDoesHaveText (ssp->name)) {
              has_chromosome = TRUE;
            }
          } else if (ssp->subtype == SUBSRC_other) {
            if (mip != NULL && (StringICmp (ssp->name, "cRNA") == 0)) {
              oldEntityID = gcp->entityID;
              oldItemID = gcp->itemID;
              oldItemtype = gcp->thistype;

              gcp->entityID = dcontext.entityID;
              gcp->itemID = dcontext.itemID;
              gcp->thistype = OBJ_SEQDESC;

              if (mip->biomol == MOLECULE_TYPE_CRNA) {
                ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "cRNA note redundant with molecule type");
              } else {
                ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "cRNA note conflicts with molecule type");
              }

              gcp->entityID = oldEntityID;
              gcp->itemID = oldItemID;
              gcp->thistype = oldItemtype;
            }
          }
        }
        if (is_transgenic && ISA_na (bsp->mol)) {
          if (SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_BIOSRC, 0, &fcontext) == NULL) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_TransgenicProblem, "Transgenic source descriptor requires presence of source feature");
          }
        }
      }
    }

    if (mip != NULL && mip->tech == MI_TECH_tsa && bsp->mol == MOLECULE_CLASS_DNA) {
      oldEntityID = gcp->entityID;
      oldItemID = gcp->itemID;
      oldItemtype = gcp->thistype;

      gcp->entityID = mipEntityID;
      gcp->itemID = mipItemID;
      gcp->thistype = mipItemtype;

      ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_ConflictingBiomolTech, "TSA sequence should not be DNA");

      gcp->entityID = oldEntityID;
      gcp->itemID = oldItemID;
      gcp->thistype = oldItemtype;
    }
    if (BioseqHasKeyword(bsp, "BARCODE") && BioseqHasKeyword(bsp, "UNVERIFIED")) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadKeyword, "Sequence has both BARCODE and UNVERIFIED keywords");
    }
  }

  if (is_neg_strand_virus && mip != NULL) {
    oldEntityID = gcp->entityID;
    oldItemID = gcp->itemID;
    oldItemtype = gcp->thistype;

    gcp->entityID = mipEntityID;
    gcp->itemID = mipItemID;
    gcp->thistype = mipItemtype;

    sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
    while (sfp != NULL) {
      has_cds = TRUE;
      if (SeqLocStrand (sfp->location) == Seq_strand_minus) {
        if (mip->biomol != MOLECULE_TYPE_GENOMIC) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Negative-strand virus with minus strand CDS should be genomic");
        }
      } else {
        if (mip->biomol != MOLECULE_TYPE_MRNA && mip->biomol != MOLECULE_TYPE_CRNA && (! is_ambisense_virus) && (! is_synthetic)) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Negative-strand virus with plus strand CDS should be mRNA or cRNA");
        }
      }
      sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &fcontext);
    }
    if (! has_cds) {
      sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_misc_feature, &fcontext);
      while (sfp != NULL) {
        if (StringISearch (sfp->comment, "nonfunctional") != NULL) {
          if (SeqLocStrand (sfp->location) == Seq_strand_minus) {
            if (mip->biomol != MOLECULE_TYPE_GENOMIC) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Negative-strand virus with nonfunctional minus strand misc_feature should be genomic");
            }
          } else {
            if (mip->biomol != MOLECULE_TYPE_MRNA && mip->biomol != MOLECULE_TYPE_CRNA && (! is_ambisense_virus)) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Negative-strand virus with nonfunctional plus strand misc_feature should be mRNA or cRNA");
            }
          }
        }
        sfp = SeqMgrGetNextFeature (bsp, sfp, 0, FEATDEF_misc_feature, &fcontext);
      }
    }

    gcp->entityID = oldEntityID;
    gcp->itemID = oldItemID;
    gcp->thistype = oldItemtype;
  }

  bvs.is_mrna = FALSE;
  bvs.is_prerna = FALSE;
  if (bsp != NULL && ISA_na (bsp->mol)) {
    if (mip != NULL) {
      if (mip->biomol == MOLECULE_TYPE_GENOMIC && mip->completeness == 1) {
        sev = SEV_ERROR;
        if (mip->tech == MI_TECH_htgs_3) {
          sev = SEV_WARNING;
        }
        for (sip = bsp->id; sip != NULL; sip = sip->next) {
          if (sip->choice == SEQID_GENBANK) {
            is_gb = TRUE;
          }
        }
        if (is_gb) {
          vnp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_title, &dcontext);
          if (vnp != NULL) {
            str = (CharPtr) vnp->data.ptrvalue;
            if (! StringHasNoText (str)) {
              if (StringISearch (str, "complete sequence") == NULL &&
                  StringISearch (str, "complete genome") == NULL) {

                oldEntityID = gcp->entityID;
                oldItemID = gcp->itemID;
                oldItemtype = gcp->thistype;

                gcp->entityID = mipEntityID;
                gcp->itemID = mipItemID;
                gcp->thistype = mipItemtype;

                if (bsp->topology == 2) {
                  ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_CompleteCircleProblem, "Circular topology has complete flag set, but title should say complete sequence or complete genome");
                } else {
                  ValidErr (vsp, sev, ERR_SEQ_DESCR_UnwantedCompleteFlag, "Suspicious use of complete");
                }

                gcp->entityID = oldEntityID;
                gcp->itemID = oldItemID;
                gcp->thistype = oldItemtype;
              }
            }
          }
        }
      } else if (mip->biomol == MOLECULE_TYPE_MRNA) {
        bvs.is_mrna = TRUE;
      } else if (mip->biomol == MOLECULE_TYPE_PRE_MRNA) {
        bvs.is_prerna = TRUE;
      }
      if (mip->biomol >= MOLECULE_TYPE_PRE_MRNA && mip->biomol <= MOLECULE_TYPE_SCRNA && bsp->mol == Seq_mol_dna) {
        /* - this is how we indicate an mRNA sequenced from a cDNA, so no error
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_RnaDnaConflict, "MolInfo says RNA, Bioseq says DNA");
        */
      }
    } else if (bsp->mol == Seq_mol_rna) {
      bvs.is_mrna = TRUE;       /* if no molinfo, assume rna is mrna */
    }
  }

  if (mip != NULL) {

    oldEntityID = gcp->entityID;
    oldItemID = gcp->itemID;
    oldItemtype = gcp->thistype;

    gcp->entityID = mipEntityID;
    gcp->itemID = mipItemID;
    gcp->thistype = mipItemtype;

    if (mip->tech == MI_TECH_sts ||
        mip->tech == MI_TECH_survey ||
        mip->tech == MI_TECH_wgs ||
        mip->tech == MI_TECH_htgs_0 || mip->tech == MI_TECH_htgs_1 ||
        mip->tech == MI_TECH_htgs_2 || mip->tech == MI_TECH_htgs_3 ||
        mip->tech == MI_TECH_composite_wgs_htgs) {
      if (mip->tech == MI_TECH_sts && bsp->mol == Seq_mol_rna && mip->biomol == MOLECULE_TYPE_MRNA) {
        /* there are some STS sequences derived from cDNAs, so do not report these */
      } else if (mip->biomol != MOLECULE_TYPE_GENOMIC) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_ConflictingBiomolTech, "HTGS/STS/GSS/WGS sequence should be genomic");
      } else if (bsp == NULL || (bsp->mol != Seq_mol_dna && bsp->mol != Seq_mol_na)) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_ConflictingBiomolTech, "HTGS/STS/GSS/WGS sequence should not be RNA");
      }
    } else if (mip->tech == MI_TECH_est && mip->biomol != MOLECULE_TYPE_MRNA) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_ConflictingBiomolTech, "EST sequence should be mRNA");
    }

    gcp->entityID = oldEntityID;
    gcp->itemID = oldItemID;
    gcp->thistype = oldItemtype;
  }

  if (ISA_aa (bsp->mol)) {
    bvs.is_aa = TRUE;
    /* check proteins in nuc-prot set have a CdRegion */
    if (vsp->bssp != NULL) {
      if (vsp->bssp->_class == 1) {     /* in a nuc-prot set */
        if (vsp->useSeqMgrIndexes) {
          sfp = SeqMgrGetCDSgivenProduct (bsp, NULL);
          if (sfp == NULL) {
            sfp = SeqMgrGetPROTgivenProduct (bsp, NULL); /* now instantiating and indexing products of protein processing */
          }
        } else {
          sfp = SeqEntryGetSeqFeat (vsp->sep, 3, NULL, NULL, 1, bsp);
        }
        if (sfp == NULL)        /* no CdRegion points to this bsp */
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_PKG_NoCdRegionPtr, "No CdRegion in nuc-prot set points to this protein");
      }
    }
  }

  if (vsp->useSeqMgrIndexes) {
    bvs.gcp = gcp;
    bvs.bsp = bsp;
    ValidateBioseqContextIndexed (bsp, &bvs);
  } else {
    GatherSeqEntry (vsp->sep, &bvs, ValidateBioseqContextGather, &gs);
  }

  vsp->gcp = gcp;               /* reset the gcp pointer changed in previous gather */
  vsp->descr = NULL;
  vsp->sfp = NULL;

  if ((!bvs.got_a_pub) && (!vsp->suppress_no_pubs) && (! vsp->seqSubmitParent)) {
    omdp = NULL;
    if (gcp != NULL) {
      omdp = ObjMgrGetData (gcp->entityID);
    }
    if (omdp == NULL || omdp->datatype != OBJ_SEQSUB) {
      sev = SEV_ERROR;
      if (!IsNoncuratedRefSeq (bsp, &sev)) {
        if (! IsWgsIntermediate (vsp->sep)) {
          ValidErr (vsp, sev, ERR_SEQ_DESCR_NoPubFound, "No publications refer to this Bioseq.");
        }
      }
    }
  }

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_LOCAL) {
      is_local = TRUE;
    } else {
      is_local_only = FALSE;
    }
    if (sip->choice == SEQID_PATENT) {
      psip = (PatentSeqIdPtr) sip->data.ptrvalue;
      if (psip != NULL) {
        ipp = psip->cit;
        if (ipp != NULL && StringICmp (ipp->country, "US") == 0)
          return;
      }
      return;
    } else if (sip->choice == SEQID_PDB) {
      isPDB = TRUE;
    } else if (sip->choice == SEQID_GENBANK ||
               sip->choice == SEQID_EMBL ||
               sip->choice == SEQID_DDBJ) {
      is_gb = TRUE;
      tsip = (TextSeqIdPtr) sip->data.ptrvalue;
      if (tsip != NULL && tsip->accession != NULL) {
        acclen = StringLen (tsip->accession);
        if (acclen == 12) {
          is_wgs = TRUE;
        } else if (acclen == 13) {
          is_wgs = TRUE;
        /*
        } else if (StringNCmp (tsip->accession, "CH", 2) == 0 ||
                   StringNCmp (tsip->accession, "CM", 2) == 0) {
          is_ch_or_cm = TRUE;
        */
        } else if (WHICH_db_accession (tsip->accession) == ACCN_NCBI_SEGSET) {
          /* NOTE '==' is appropriate here, rather than '|', because we
           * really do only want to suppress if the type is exactly ACCN_NCBI_SEGSET
           * and NOT if the type is ACCN_NCBI_SEGSET | ACCN_AMBIGOUS_MOL (prefix is AH)
           */
          is_ch_or_cm = TRUE;
        }
      }
    } else if (sip->choice == SEQID_OTHER) {
      is_refseq = TRUE;
      tsip = (TextSeqIdPtr) sip->data.ptrvalue;
      if (tsip != NULL && tsip->accession != NULL) {
        if (StringNCmp (tsip->accession, "NM_", 3) == 0 ||
            StringNCmp (tsip->accession, "NP_", 3) == 0 ||
            StringNCmp (tsip->accession, "NG_", 3) == 0 ||
            StringNCmp (tsip->accession, "NR_", 3) == 0) {
          is_gb = TRUE;
        } else if (StringNCmp (tsip->accession, "NC_", 3) == 0) {
          is_nc = TRUE;
        } else if (StringNCmp (tsip->accession, "AC_", 3) == 0) {
          is_ac = TRUE;
        }
      }
    }
  }
  if (is_nc || is_ac) {
    if (! is_prokaryote && ! is_organelle && ! has_chromosome && ! is_plasmid) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_MissingChromosome, "Missing chromosome qualifier on NC or AC RefSeq record");
    }
  }
  if (! is_local) {
    is_local_only = FALSE;
  }
  if (is_wgs) {
    if (mip == NULL || mip->tech != MI_TECH_wgs) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "WGS accession should have Mol-info.tech of wgs");
    }
  } else if (mip != NULL && mip->tech == MI_TECH_wgs && is_gb) {
    if (is_ch_or_cm || is_local_only) {
      /* skip warning if CH or CM or SEQID_LOCAL only */
    } else {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Mol-info.tech of wgs should have WGS accession");
    }
  }
  if (is_nc) {
    if (mip != NULL && mip->biomol != MOLECULE_TYPE_GENOMIC && mip->biomol != MOLECULE_TYPE_CRNA) {
      if (ISA_na (bsp->mol)) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "NC nucleotide should be genomic or cRNA");
      }
    }
  }

  if (bsp != NULL && is_refseq) {
    ValidateRefSeqTitle (bsp, vsp, is_virus);
  }

  if ((!bvs.last_org) && (!vsp->suppress_no_biosrc))
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_NoOrgFound, "No organism name has been applied to this Bioseq.  Other qualifiers may exist.");


  if ((bvs.is_aa) && (bvs.num_full_length_prot_ref == 0) && (!isPDB) && (bsp->repr != Seq_repr_virtual)) {
    parent = SeqMgrGetParentOfPart (bsp, NULL);
    if (parent == NULL || SeqMgrGetBestProteinFeature (bsp, NULL) == NULL) {

      oldEntityID = gcp->entityID;
      oldItemID = gcp->itemID;
      oldItemtype = gcp->thistype;

      if (SeqMgrGetCDSgivenProduct (bsp, &fcontext) != NULL) {
        gcp->entityID = fcontext.entityID;
        gcp->itemID = fcontext.itemID;
        gcp->thistype = OBJ_SEQFEAT;
      }

      ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_NoProtRefFound, "No full length Prot-ref feature applied to this Bioseq");

      gcp->entityID = oldEntityID;
      gcp->itemID = oldItemID;
      gcp->thistype = oldItemtype;
    }
  } else if (bvs.is_aa && bvs.num_full_length_prot_ref > 1) {
    if (bvs.num_justprot > 1 ||
        bvs.num_preprot > 1 ||
        bvs.num_matpep > 1 ||
        bvs.num_sigpep > 1 ||
        bvs.num_transpep > 1) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_MultipleProtRefs, "%d full-length protein features present on protein",
                (int) bvs.num_full_length_prot_ref);
    }
  }

  /* now flag missing molinfo even if not in Sequin */
  if (mip == NULL && (!isPDB)) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_NoMolInfoFound, "No Mol-info applies to this Bioseq");
  }

#if 0 /* temporarily suppress */
  /* if tech is TSA, must have assembly */
  if (mip != NULL && mip->tech == MI_TECH_tsa
      && (bsp->hist == NULL || bsp->hist->assembly == NULL)) {
    SeqIdWrite (bsp->id, buf1, PRINTID_FASTA_SHORT, sizeof (buf1) - 1);
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_TSAHistAssemblyMissing, "TSA record %s should have Seq-hist.assembly", buf1);
  }
#endif

  /* look for genes that overlap two other genes */
  FindMultiGeneOverlaps (bsp, vsp);

  /* TSA checks */
  ValidateTSASequenceForNs (bsp, vsp);

  /* validate structured comments in context */
  ValidateStructuredCommentsInContext (bsp, vsp);
}

/*****************************************************************************
*
*   ValidateSeqFeat(gcp)
*
*****************************************************************************/
static Boolean EmptyOrNullString (CharPtr str)
{
  Char            ch;

  if (str == NULL)
    return TRUE;
  ch = *str;
  while (ch != '\0') {
    if (ch > ' ' && ch <= '~')
      return FALSE;
    str++;
    ch = *str;
  }
  return TRUE;
}

static void CheckPeptideOnCodonBoundary (ValidStructPtr vsp, GatherContextPtr gcp, SeqFeatPtr sfp, CharPtr key)
{
  SeqFeatPtr      cds;
  CdRegionPtr     crp;
  SeqLocPtr       first = NULL, last = NULL, slp = NULL;
  Boolean         partial5, partial3;
  Int4            pos1, pos2, adjust = 0, mod1, mod2;

  if (SeqLocStop (sfp->location) == 2150166) {
    mod1 = 0;
  }
  cds = SeqMgrGetOverlappingCDS (sfp->location, NULL);
  if (cds == NULL)
    return;
  crp = (CdRegionPtr) cds->data.value.ptrvalue;
  if (crp == NULL)
    return;
  if (crp->frame == 2) {
    adjust = 1;
  } else if (crp->frame == 3) {
    adjust = 2;
  }

  while ((slp = SeqLocFindNext (sfp->location, slp)) != NULL) {
    last = slp;
    if (first == NULL) {
      first = slp;
    }
  }
  if (first == NULL || last == NULL)
    return;

  pos1 = GetOffsetInLoc (first, cds->location, SEQLOC_START) - adjust;
  pos2 = GetOffsetInLoc (last, cds->location, SEQLOC_STOP) - adjust;
  mod1 = pos1 % 3;
  mod2 = pos2 % 3;

  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  if (partial5) {
    mod1 = 0;
  }
  if (partial3) {
    mod2 = 2;
  }

  if (pos1 < 0 && pos2 < 0 && StringICmp (key, "sig_peptide") == 0) {
    /* ignore special case of sig_peptide completely before codon_start of CDS */
  } else if (mod1 != 0 && mod2 != 2) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PeptideFeatOutOfFrame, "Start and stop of %s are out of frame with CDS codons", key);
  } else if (mod1 != 0) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PeptideFeatOutOfFrame, "Start of %s is out of frame with CDS codons", key);
  } else if (mod2 != 2) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PeptideFeatOutOfFrame, "Stop of %s is out of frame with CDS codons", key);
  }
}

static CharPtr  legal_repeat_types[] = {
  "tandem", "inverted", "flanking", "terminal",
  "direct", "dispersed", "other", NULL
};

static CharPtr legal_cons_splice_strings [] = {
  "(5'site:YES, 3'site:YES)",
  "(5'site:YES, 3'site:NO)",
  "(5'site:YES, 3'site:ABSENT)",
  "(5'site:NO, 3'site:YES)",
  "(5'site:NO, 3'site:NO)",
  "(5'site:NO, 3'site:ABSENT)",
  "(5'site:ABSENT, 3'site:YES)",
  "(5'site:ABSENT, 3'site:NO)",
  "(5'site:ABSENT, 3'site:ABSENT)",
  NULL
};

static CharPtr legal_mobile_element_strings [] = {
  "transposon",
  "retrotransposon",
  "integron",
  "insertion sequence",
  "non-LTR retrotransposon",
  "SINE",
  "MITE",
  "LINE",
  "other",
  NULL
};

NLM_EXTERN Boolean LookForECnumberPattern (CharPtr str)

{
  Char     ch;
  Boolean  is_ambig;
  Int2     numdashes;
  Int2     numdigits;
  Int2     numperiods;
  CharPtr  ptr;

  if (StringHasNoText (str)) return FALSE;

  is_ambig = FALSE;
  numperiods = 0;
  numdigits = 0;
  numdashes = 0;

  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    if (IS_DIGIT (ch)) {
      numdigits++;
      if (is_ambig) {
        is_ambig = FALSE;
        numperiods = 0;
        numdashes = 0;
      }
      ptr++;
      ch = *ptr;
    } else if (ch == '-') {
      numdashes++;
      is_ambig = TRUE;
      ptr++;
      ch = *ptr;
    } else if (ch == 'n') {
      numdashes++;
      is_ambig = TRUE;
      ptr++;
      ch = *ptr;
    } else if (ch == '.') {
      numperiods++;
      if (numdigits > 0 && numdashes > 0) {
        is_ambig = FALSE;
        numperiods = 0;
        numdigits = 0;
        numdashes = 0;
      } else if (numdigits == 0 && numdashes == 0) {
        is_ambig = FALSE;
        numperiods = 0;
        numdigits = 0;
        numdashes = 0;
      } else if (numdashes > 1) {
        is_ambig = FALSE;
        numperiods = 0;
        numdigits = 0;
        numdashes = 0;
      }
      numdigits = 0;
      numdashes = 0;
      ptr++;
      ch = *ptr;
    } else {
      if (numperiods == 3) {
        if (numdigits > 0 && numdashes > 0) {
        is_ambig = FALSE;
        numperiods = 0;
        numdigits = 0;
        numdashes = 0;
        } else if (numdigits > 0 || numdashes == 1) return TRUE;
      }
      ptr++;
      ch = *ptr;
      is_ambig = FALSE;
      numperiods = 0;
      numdigits = 0;
      numdashes = 0;
    }
  }

  if (numperiods == 3) {
    if (numdigits > 0 && numdashes > 0) return FALSE;
    if (numdigits > 0 || numdashes == 1) return TRUE;
  }

  return FALSE;
}

NLM_EXTERN Boolean ValidateECnumber (CharPtr str)

{
  Char     ch;
  Boolean  is_ambig;
  Int2     numdashes;
  Int2     numdigits;
  Int2     numperiods;
  CharPtr  ptr;

  if (StringHasNoText (str)) return FALSE;

  is_ambig = FALSE;
  numperiods = 0;
  numdigits = 0;
  numdashes = 0;

  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    if (IS_DIGIT (ch)) {
      numdigits++;
      if (is_ambig) return FALSE;
      ptr++;
      ch = *ptr;
    } else if (ch == '-') {
      numdashes++;
      is_ambig = TRUE;
      ptr++;
      ch = *ptr;
    } else if (ch == 'n') {
      if (numperiods == 3 && numdigits == 0 && IS_DIGIT(*(ptr + 1))) {
        /* allow/ignore n in first position of fourth number to not mean ambiguous, if followed by digit */
      } else {
        numdashes++;
        is_ambig = TRUE;
      }
      ptr++;
      ch = *ptr;
    } else if (ch == '.') {
      numperiods++;
      if (numdigits > 0 && numdashes > 0) return FALSE;
      if (numdigits == 0 && numdashes == 0) return FALSE;
      if (numdashes > 1) return FALSE;
      numdigits = 0;
      numdashes = 0;
      ptr++;
      ch = *ptr;
    } else {
      ptr++;
      ch = *ptr;
    }
  }

  if (numperiods == 3) {
    if (numdigits > 0 && numdashes > 0) return FALSE;
    if (numdigits > 0 || numdashes == 1) return TRUE;
  }

  return FALSE;
}

NLM_EXTERN void ECNumberFSAFreeAll (void)

{
  CtrySetPtr  ctsp;
  TextFsaPtr  fsa;

  fsa = (TextFsaPtr) GetAppProperty ("SpecificECNumberFSA");
  if (fsa != NULL) {
    SetAppProperty ("SpecificECNumberFSA", NULL);
    TextFsaFree (fsa);
  }

  fsa = (TextFsaPtr) GetAppProperty ("AmbiguousECNumberFSA");
  if (fsa != NULL) {
    SetAppProperty ("AmbiguousECNumberFSA", NULL);
    TextFsaFree (fsa);
  }

  fsa = (TextFsaPtr) GetAppProperty ("DeletedECNumberFSA");
  if (fsa != NULL) {
    SetAppProperty ("DeletedECNumberFSA", NULL);
    TextFsaFree (fsa);
  }

  fsa = (TextFsaPtr) GetAppProperty ("ReplacedEECNumberFSA");
  if (fsa != NULL) {
    SetAppProperty ("ReplacedEECNumberFSA", NULL);
    TextFsaFree (fsa);
  }

  fsa = (TextFsaPtr) GetAppProperty ("BodiesOfWaterFSA");
  if (fsa != NULL) {
    SetAppProperty ("BodiesOfWaterFSA", NULL);
    TextFsaFree (fsa);
  }

  ctsp = (CtrySetPtr) GetAppProperty ("CountryLatLonData");
  if (ctsp != NULL) {
    SetAppProperty ("CountryLatLonData", NULL);
    FreeLatLonCountryData (ctsp);
  }

  ctsp = (CtrySetPtr) GetAppProperty ("WaterLatLonData");
  if (ctsp != NULL) {
    SetAppProperty ("WaterLatLonData", NULL);
    FreeLatLonCountryData (ctsp);
  }

  ic_code_data = MemFree (ic_code_data);
  ic_code_list = ValNodeFreeData (ic_code_list);
}

static TextFsaPtr GetECNumberFSA (CharPtr prop, CharPtr file, CharPtr PNTR local, size_t numitems, Boolean trimAtTab)

{
  FileCache   fc;
  FILE        *fp = NULL;
  TextFsaPtr  fsa;
  Int2        i;
  Char        line [512];
  Char        path [PATH_MAX];
  CharPtr     ptr;
  ErrSev      sev;
  CharPtr     str;
  Char        tmp [128];

  fsa = (TextFsaPtr) GetAppProperty (prop);
  if (fsa != NULL) return fsa;

  if (FindPath ("ncbi", "ncbi", "data", path, sizeof (path))) {
    FileBuildPath (path, NULL, file);
    sev = ErrSetMessageLevel (SEV_ERROR);
    fp = FileOpen (path, "r");
    ErrSetMessageLevel (sev);
  }

  fsa = TextFsaNew ();
  if (fsa != NULL) {
    if (fp != NULL) {
      FileCacheSetup (&fc, fp);

      str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
      while (str != NULL) {
        if (StringDoesHaveText (str)) {
          if (trimAtTab) {
            ptr = StringChr (str, '\t');
            if (ptr != NULL) {
              *ptr = '\0';
            }
          }
          if (StringLen (str) + 3 < sizeof (tmp)) {
            StringCpy (tmp, " ");
            StringCat (tmp, str);
            StringCat (tmp, " ");
            TextFsaAdd (fsa, tmp);
          }
        }
        str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
      }

    } else if (local != NULL) {
      for (i = 0; /* local [i] != NULL */ i < numitems; i++) {
        str = local [i];
        if (StringDoesHaveText (str)) {
          if (StringLen (str) + 3 < sizeof (tmp)) {
            StringCpy (tmp, " ");
            StringCat (tmp, str);
            StringCat (tmp, " ");
            TextFsaAdd (fsa, tmp);
          }
        }
      }
    }
  }

  if (fp != NULL) {
    FileClose (fp);
  }

  SetAppProperty (prop, (Pointer) fsa);

  return fsa;
}

static TextFsaPtr GetSpecificECNumberFSA (void)

{
  return (GetECNumberFSA ("SpecificECNumberFSA", "ecnum_specific.txt", (CharPtr PNTR) kECNum_specific, sizeof (kECNum_specific) / sizeof (char*), TRUE));
}

static TextFsaPtr GetAmbiguousECNumberFSA (void)

{
  return (GetECNumberFSA ("AmbiguousECNumberFSA", "ecnum_ambiguous.txt", (CharPtr PNTR) kECNum_ambiguous, sizeof (kECNum_ambiguous) / sizeof (char*), TRUE));
}

static TextFsaPtr GetDeletedECNumberFSA (void)

{
  return (GetECNumberFSA ("DeletedECNumberFSA", "ecnum_deleted.txt", (CharPtr PNTR) kECNum_deleted, sizeof (kECNum_deleted) / sizeof (char*), TRUE));
}

static TextFsaPtr GetReplacedECNumberFSA (void)

{
  return (GetECNumberFSA ("ReplacedEECNumberFSA", "ecnum_replaced.txt", (CharPtr PNTR) kECNum_replaced, sizeof (kECNum_replaced) / sizeof (char*), TRUE));
}

NLM_EXTERN Boolean ECnumberNotInList (CharPtr str)

{
  Char        ch;
  TextFsaPtr  fsa;
  ValNodePtr  matches;
  CharPtr     ptr;
  Int4        state;

  fsa = GetSpecificECNumberFSA ();
  if (fsa == NULL) return FALSE;

  state = 0;
  matches = NULL;
  state = TextFsaNext (fsa, state, ' ', &matches);
  for (ptr = str, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
    state = TextFsaNext (fsa, state, ch, &matches);
  }
  state = TextFsaNext (fsa, state, ' ', &matches);
  if (matches != NULL) return FALSE;

  fsa = GetAmbiguousECNumberFSA ();
  if (fsa == NULL) return FALSE;

  state = 0;
  matches = NULL;
  state = TextFsaNext (fsa, state, ' ', &matches);
  for (ptr = str, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
    state = TextFsaNext (fsa, state, ch, &matches);
  }
  state = TextFsaNext (fsa, state, ' ', &matches);
  if (matches != NULL) return FALSE;

  return TRUE;
}

NLM_EXTERN Boolean ECnumberWasDeleted (CharPtr str)

{
  Char        ch;
  TextFsaPtr  fsa;
  ValNodePtr  matches;
  CharPtr     ptr;
  Int4        state;

  fsa = GetDeletedECNumberFSA ();
  if (fsa == NULL) return FALSE;

  state = 0;
  matches = NULL;
  state = TextFsaNext (fsa, state, ' ', &matches);
  for (ptr = str, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
    state = TextFsaNext (fsa, state, ch, &matches);
  }
  state = TextFsaNext (fsa, state, ' ', &matches);
  if (matches != NULL) return TRUE;

  return FALSE;
}

NLM_EXTERN Boolean ECnumberWasReplaced (CharPtr str)

{
  Char        ch;
  TextFsaPtr  fsa;
  ValNodePtr  matches;
  CharPtr     ptr;
  Int4        state;

  fsa = GetReplacedECNumberFSA ();
  if (fsa == NULL) return FALSE;

  state = 0;
  matches = NULL;
  state = TextFsaNext (fsa, state, ' ', &matches);
  for (ptr = str, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
    state = TextFsaNext (fsa, state, ch, &matches);
  }
  state = TextFsaNext (fsa, state, ' ', &matches);
  if (matches != NULL) return TRUE;

  return FALSE;
}

static Boolean RptUnitIsBaseRange (CharPtr str, Int4Ptr fromP, Int4Ptr toP)

{
  CharPtr   ptr;
  Char      tmp [32];
  long int  val;

  if (StringLen (str) > 25) return FALSE;
  StringNCpy_0 (tmp, str, sizeof (tmp));
  ptr = StringStr (tmp, "..");
  if (ptr == NULL) return FALSE;
  *ptr = '\0';
  if (StringHasNoText (tmp)) return FALSE;
  if (sscanf (tmp, "%ld", &val) != 1 || val < 1) return FALSE;
  if (fromP != NULL) {
    *fromP = val - 1;
  }
  ptr += 2;
  if (StringHasNoText (ptr)) return FALSE;
  if (sscanf (ptr, "%ld", &val) != 1 || val < 1) return FALSE;
  if (toP != NULL) {
    *toP = val - 1;
  }
  return TRUE;
}


static void ValidateRptUnit (ValidStructPtr vsp, GatherContextPtr gcp, SeqFeatPtr sfp, GBQualPtr gbqual, Int2 qual, CharPtr key)

{
  Boolean            badchars, found, just_nuc_letters, multi_rpt_unit;
  Char               ch;
  SeqMgrFeatContext  context;
  Int4               from = -1, to = -1, ffrom, fto, ftmp;
  CharPtr            ptr, tmp;

  if (vsp == NULL || gcp == NULL || sfp == NULL || gbqual == NULL || gbqual->val == NULL || key == NULL) return;

  found = FALSE;
  multi_rpt_unit = TRUE;
  for (ptr = gbqual->val, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
    if (ch <= ' ') {
      found = TRUE;
    } else if (ch == '(' || ch == ')' || ch == ',' || ch == '.' || IS_DIGIT (ch)) {
    } else {
      multi_rpt_unit = FALSE;
    }
  }
  /*
  if (found) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "Illegal value for qualifier %s", gbqual->qual);
  } else if ((!multi_rpt_unit) && StringLen (gbqual->val) > 48) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "Illegal value for qualifier %s", gbqual->qual);
  }
  */
  if (StringICmp (key,"repeat_region") == 0 && qual == GBQUAL_rpt_unit_seq && ! multi_rpt_unit) {
    if (StringLen (gbqual->val) <= SeqLocLen (sfp->location)) {
      just_nuc_letters = TRUE;
      for (ptr = gbqual->val, ch = *ptr; ch != '\0' && just_nuc_letters; ptr++, ch = *ptr) {
        if (StringChr ("ACGTNacgtn", ch) == NULL) {
          just_nuc_letters = FALSE;
        }
      }
      if (just_nuc_letters) {
        tmp = GetSequenceByFeature (sfp);
        if (tmp != NULL) {
          if (StringISearch (tmp, gbqual->val) == NULL) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_InvalidQualifierValue, "repeat_region /rpt_unit and underlying sequence do not match");
          }
          MemFree (tmp);
        }
      } else {
        /* illegal character test is now handled better at the end of this function */
        /*
        ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_InvalidQualifierValue, "rpt_unit_seq qualifier contains invalid characters");
        */
      }
    } else {
      ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_InvalidQualifierValue, "Length of rpt_unit_seq is greater than feature length");
    }
  }

  if (qual == GBQUAL_rpt_unit_range) {
    if (RptUnitIsBaseRange (gbqual->val, &from, &to)) {
      if (SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL, 0, 0, sfp, &context) == sfp) {
        if (from < context.left || from > context.right || to < context.left || to > context.right) {
          /* could be segmented sequence */
          ffrom = SeqLocStart (sfp->location);
          fto = SeqLocStop (sfp->location);
          if (ffrom > fto) {
            ftmp = ffrom;
            ffrom = fto;
            fto = ftmp;
          }
          if (from < ffrom || from > fto || to < ffrom || to > fto) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_RptUnitRangeProblem, "/rpt_unit_range is not within sequence length");
          }
        }
      }
    } else {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_InvalidQualifierValue, "/rpt_unit_range is not a base range");
    }
  }

  if (qual == GBQUAL_rpt_unit_seq) {
    badchars = FALSE;
    for (ptr = gbqual->val, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
      if (ch <= ' ') {
        badchars = TRUE;
      } else if (ch == '(' || ch == ')' || IS_DIGIT (ch) || IS_ALPHA (ch)) {
      } else if (ch == ',' || ch == ';') {
      } else {
        badchars = TRUE;
      }
    }
    if (badchars) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_InvalidQualifierValue, "/rpt_unit_seq has illegal characters");
    }
  }
}


static Boolean IsGbIndexQualPairValid (Int2 index, Int2 val)
{
  Int2    i;
  Boolean found = FALSE;

  for (i = 0; i < ParFlat_GBFeat[index].opt_num && !found; i++) {
    if (ParFlat_GBFeat[index].opt_qual[i] == val) {
      found = TRUE;
    }
  }
  for (i = 0; i < ParFlat_GBFeat[index].mand_num && !found; i++) {
    if (ParFlat_GBFeat[index].mand_qual[i] == val) {
      found = TRUE;
    }
  }
  return found;
}


NLM_EXTERN CharPtr GetGBFeatKeyForFeature (SeqFeatPtr sfp)
{
  CharPtr key = NULL;
  ImpFeatPtr ifp;

  if (sfp == NULL) {
    return NULL;
  }

  if (sfp->data.choice == SEQFEAT_IMP) {
    ifp = (ImpFeatPtr) sfp->data.value.ptrvalue;
    if (StringCmp (ifp->key, "-") == 0) {
      key = StringSave ("misc_feature");
    } else {
      key = StringSaveNoNull (ifp->key);
    }
  } else {
    key = StringSaveNoNull (FeatDefTypeLabel (sfp));
    if (StringCmp (key, "Gene") == 0) {
      *key = 'g';
    } else if (StringCmp (key, "preRNA") == 0) {
      key = MemFree (key);
      key = StringSave ("precursor_RNA");
    }
  }
  return key;
}


NLM_EXTERN Boolean ShouldSuppressGBQual(Uint1 subtype, CharPtr qual_name)
{
  if (StringHasNoText (qual_name)) {
    return FALSE;
  }

  /* always suppress experiment and inference quals */
  if (StringCmp (qual_name, "experiment") == 0 || StringCmp (qual_name, "inference") == 0) {
    return TRUE;
  }
  
  if (subtype == FEATDEF_ncRNA) {
    if (StringCmp (qual_name, "product") == 0
        || StringCmp (qual_name, "ncRNA_class") == 0) {
      return TRUE;
    }
  } else if (subtype == FEATDEF_tmRNA) {
    if (StringCmp (qual_name, "product") == 0
        || StringCmp (qual_name, "tag_peptide") == 0) {
      return TRUE;
    }
  } else if (subtype == FEATDEF_otherRNA) {
    if (StringCmp (qual_name, "product") == 0) {
      return TRUE;
    }
  }

  return FALSE;
}
    

NLM_EXTERN Boolean ShouldBeAGBQual (Uint1 subtype, Int2 qual, Boolean allowProductGBQual)

{
  if (qual < 0) return FALSE;
  if (allowProductGBQual && qual == GBQUAL_product) return TRUE;
  if (qual == GBQUAL_citation ||
      qual == GBQUAL_db_xref ||
      qual == GBQUAL_evidence ||
      qual == GBQUAL_exception ||
      qual == GBQUAL_gene ||
      qual == GBQUAL_gene_synonym ||
      qual == GBQUAL_insertion_seq ||
      qual == GBQUAL_label ||
      qual == GBQUAL_locus_tag ||
      qual == GBQUAL_non_functional ||
      qual == GBQUAL_note ||
      qual == GBQUAL_partial ||
      qual == GBQUAL_product ||
      qual == GBQUAL_pseudo ||
      qual == GBQUAL_pseudogene ||
      qual == GBQUAL_rpt_unit ||
      qual == GBQUAL_transposon ||
      qual == GBQUAL_experiment ||
      qual == GBQUAL_trans_splicing ||
      qual == GBQUAL_ribosomal_slippage ||
      qual == GBQUAL_standard_name ||
      qual == GBQUAL_inference) 
  {
    return FALSE;
  }
  if (subtype == FEATDEF_CDS) 
  {
    if (qual == GBQUAL_codon_start 
        || qual == GBQUAL_codon 
        || qual == GBQUAL_EC_number
        || qual == GBQUAL_gdb_xref
        || qual == GBQUAL_number
        || qual == GBQUAL_protein_id
        || qual == GBQUAL_transl_except
        || qual == GBQUAL_transl_table
        || qual == GBQUAL_translation
        || qual == GBQUAL_allele
        || qual == GBQUAL_function
        || qual == GBQUAL_old_locus_tag)
    {
      return FALSE;
    }
  }
  if (qual == GBQUAL_map && subtype != FEATDEF_ANY && subtype != FEATDEF_repeat_region && subtype != FEATDEF_gap) return FALSE;
  if (qual == GBQUAL_operon && subtype != FEATDEF_ANY && subtype != FEATDEF_operon) return FALSE;
  if (Nlm_GetAppProperty ("SequinUseEMBLFeatures") == NULL) 
  {
    if (qual == GBQUAL_usedin) 
    {
      return FALSE;
    }
  }

  if (qual > -1 && ShouldSuppressGBQual (subtype, ParFlat_GBQual_names [qual].name)) {
    return FALSE;
  }

  return TRUE;
}


static CharPtr sWrongQualReasons[] = {
  "conflicting codon_start values",
  "codon_start value should be 1, 2, or 3" 
};

typedef enum {
  eWrongQualReason_conflicting_codon_start = 0,
  eWrongQualReason_bad_codon_start_value
} EWrongQualReason;

/* 
 * Return values:
 * 1: yes
 * 0: no
 * -1: don't know
 * 2: no for special reasons
 */
NLM_EXTERN Int4 IsQualValidForFeature (GBQualPtr gbqual, SeqFeatPtr sfp)
{
  CharPtr     key = NULL;
  Int2        val;
  Int4        rval = -1;
  Int2        index;
  CdRegionPtr crp;

  if (sfp == NULL || gbqual == NULL) {
    return -1;
  }

  key = GetGBFeatKeyForFeature (sfp);
  index = GBFeatKeyNameValid (&key, FALSE);
  key = MemFree (key);

  if (index == -1) {
    /* unknown */
    rval = -1;
  } else if (StringCmp (gbqual->qual, "gsdb_id") == 0) {
    /* force good */
    rval = 1;
  } else if (sfp->data.choice == SEQFEAT_GENE &&
             (StringCmp (gbqual->qual, "gen_map") == 0 ||
              StringCmp (gbqual->qual, "cyt_map") == 0 ||
              StringCmp (gbqual->qual, "rad_map") == 0)) {
    rval = 1;
  } else if (sfp->data.choice == SEQFEAT_CDREGION
             && StringCmp (gbqual->qual, "orig_transcript_id") == 0) {
    rval = 1;
  } else if (sfp->data.choice == SEQFEAT_RNA && 
             (StringCmp (gbqual->qual, "orig_protein_id") == 0 ||
              StringCmp (gbqual->qual, "orig_transcript_id") == 0)) {
    rval = 1;
  } else if ((val = GBQualNameValid (gbqual->qual)) == -1) {
    rval = -1;
  } else if (sfp->data.choice == SEQFEAT_CDREGION
             && val == GBQUAL_codon_start) {
    crp = (CdRegionPtr) sfp->data.value.ptrvalue;
    if (crp != NULL) {
      if (crp->frame > 0) {
        rval = eWrongQualReason_conflicting_codon_start + 2;
      } else {
        rval = eWrongQualReason_bad_codon_start_value + 2;
      }
    }
  } else if (IsGbIndexQualPairValid (index, val)) {
    rval = 1;
  } else {
    rval = 0;
  }
  return rval;
}


static void ValidateImpFeat (ValidStructPtr vsp, GatherContextPtr gcp, SeqFeatPtr sfp, ImpFeatPtr ifp)

{
  Int2            adv;
  BioseqPtr       bsp;
  Char            ch;
  Boolean         failed;
  Boolean         found;
  IntFuzzPtr      fuzz;
  GBQualPtr       gbqual;
  SeqMgrFeatContext gcontext;
  SeqFeatPtr      gene;
  GeneRefPtr      grp;
  Int2            i;
  Int2            index;
  Boolean         just_nuc_letters;
  Boolean         just_prt_letters;
  CharPtr         key;
  size_t          len;
  Boolean         multi_compare;
  Boolean         no_white_space;
  Boolean         only_digits;
  CharPtr         ptr;
  Int2            qual;
  Char            range[32];
  ErrSev          sev;
  SeqIntPtr       sint;
  SeqIdPtr        sip;
  SeqLocPtr       slp;
  SeqPntPtr       spp;
  CharPtr         str;
  CharPtr         tmp;
  Int2            val;
  Int4            qvalid;

  if (vsp == NULL || gcp == NULL || sfp == NULL || ifp == NULL)
    return;
  if (StringCmp (ifp->key, "-") == 0) {
    key = StringSave ("misc_feature");
  } else {
    key = StringSaveNoNull (ifp->key);
  }
  index = GBFeatKeyNameValid (&key, FALSE);
  if (index == -1) {
    if (key != NULL) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UnknownImpFeatKey, "Unknown feature key %s", key);
    } else {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UnknownImpFeatKey, "NULL feature key");
    }
  } else if (StringICmp (key, "virion") == 0 ||
             StringICmp (key, "mutation") == 0 ||
             StringICmp (key, "allele") == 0 ||
             StringICmp (key, "Import") == 0) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_UnknownImpFeatKey, "Feature key %s is no longer legal", key);
  } else if (StringICmp (key, "polyA_site") == 0) {
    sev = SEV_WARNING;
    if (vsp->is_refseq_in_sep) {
      sev = SEV_ERROR;
    }
    if (SeqLocStart (sfp->location) != SeqLocStop (sfp->location)) {
      ValidErr (vsp, sev, ERR_SEQ_FEAT_PolyAsiteNotPoint, "PolyA_site should be a single point");
    }
  } else if (StringICmp (key, "polyA_signal") == 0) {
    sev = SEV_WARNING;
    if (vsp->is_refseq_in_sep) {
      sev = SEV_ERROR;
    }
    if (SeqLocStart (sfp->location) == SeqLocStop (sfp->location)) {
      ValidErr (vsp, sev, ERR_SEQ_FEAT_PolyAsignalNotRange, "PolyA_signal should be a range");
    }
  } else if (StringICmp (key, "mat_peptide") == 0 ||
             StringICmp (key, "sig_peptide") == 0 ||
             StringICmp (key, "transit_peptide") == 0) {
    sev = SEV_WARNING;
    if (vsp->is_refseq_in_sep) {
      sev = SEV_ERROR;
    }
    ValidErr (vsp, sev, ERR_SEQ_FEAT_InvalidForType, "Peptide processing feature should be converted to the appropriate protein feature subtype");
    CheckPeptideOnCodonBoundary (vsp, gcp, sfp, key);
  } else if (StringICmp (key, "preprotein") == 0 ||
             StringICmp (key, "proprotein") == 0) {
    sev = SEV_WARNING;
    if (vsp->is_refseq_in_sep) {
      sev = SEV_ERROR;
    }
    ValidErr (vsp, sev, ERR_SEQ_FEAT_InvalidForType, "Peptide processing feature should be converted to the appropriate protein feature subtype");
  } else if (StringICmp (key, "mRNA") == 0 ||
             StringICmp (key, "tRNA") == 0 ||
             StringICmp (key, "rRNA") == 0 ||
             StringICmp (key, "snRNA") == 0 ||
             StringICmp (key, "scRNA") == 0 ||
             StringICmp (key, "snoRNA") == 0 ||
             StringICmp (key, "misc_RNA") == 0 ||
             StringICmp (key, "precursor_RNA") == 0) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidForType,
              "RNA feature should be converted to the appropriate RNA feature subtype, location should be converted manually");
  } else if (StringICmp (key, "CDS") == 0) {
    failed = TRUE;              /* impfeat CDS must be pseudo; fail if not */
    if (sfp->pseudo) {
      failed = FALSE;
    } else {
      grp = SeqMgrGetGeneXref (sfp);
      if (grp != NULL && grp->pseudo) {
        failed = FALSE;
      } else {
        gene = SeqMgrGetOverlappingGene (sfp->location, &gcontext);
        if (gene != NULL) {
          if (gene->pseudo) {
            failed = FALSE;
          } else {
            grp = (GeneRefPtr) gene->data.value.ptrvalue;
            if (grp != NULL && grp->pseudo) {
              failed = FALSE;
            }
          }
        }
      }
    }
    for (gbqual = sfp->qual; gbqual != NULL; gbqual = gbqual->next) {
      if (StringCmp (gbqual->qual, "translation") == 0) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_ImpCDShasTranslation, "ImpFeat CDS with /translation found");
      }
    }
    if (failed) {
      ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_ImpCDSnotPseudo, "ImpFeat CDS should be pseudo");
    }
  } else if (StringICmp (key, "misc_feature") == 0) {
    for (gbqual = sfp->qual; gbqual != NULL; gbqual = gbqual->next) {
      if (StringCmp (gbqual->qual, "standard_name") == 0) {
        if (StringCmp (gbqual->val, "Vector Contamination") == 0) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_VectorContamination, "Vector Contamination region should be trimmed from sequence");
        }
      }
    }
    if (StringHasNoText(sfp->comment) && sfp->qual == NULL && sfp->dbxref == NULL) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_NeedsNote, "A note or other qualifier is required for a misc_feature");
    }
  }
  for (gbqual = sfp->qual; gbqual != NULL; gbqual = gbqual->next) {
    qvalid = IsQualValidForFeature (gbqual, sfp);
    if (qvalid == 0) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_WrongQualOnImpFeat, "Wrong qualifier %s for feature %s", gbqual->qual, key);
    }

    if (StringCmp (gbqual->qual, "gsdb_id") == 0) {
      continue;
    }
    val = GBQualNameValid (gbqual->qual);
    if (val == -1) {
      if (gbqual->qual != NULL) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UnknownImpFeatQual, "Unknown qualifier %s", gbqual->qual);
      } else {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UnknownImpFeatQual, "NULL qualifier");
      }
    } else if (index != -1) {
      if (gbqual->val != NULL) {
        if (val == GBQUAL_rpt_type) {
          failed = FALSE;
          tmp = StringSave (gbqual->val);
          str = tmp;
          if (*str == '(') {
            str++;
          }
          while (!StringHasNoText (str)) {
            ptr = StringChr (str, ',');
            if (ptr == NULL) {
              ptr = StringChr (str, ')');
            }
            if (ptr != NULL) {
              *ptr = '\0';
              ptr++;
            }
            found = FALSE;
            for (i = 0; legal_repeat_types[i] != NULL; i++) {
              if (StringICmp (str, legal_repeat_types[i]) == 0) {
                found = TRUE;
                break;
              }
            }
            if (!found) {
              failed = TRUE;
            }
            str = ptr;
          }
          MemFree (tmp);
          if (failed) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "%s is not a legal value for qualifier %s", gbqual->val, gbqual->qual);
          }
        } else if (val == GBQUAL_rpt_unit || val == GBQUAL_rpt_unit_range || val == GBQUAL_rpt_unit_seq) {
          ValidateRptUnit (vsp, gcp, sfp, gbqual, val, key);
        } else if (val == GBQUAL_label) {
          no_white_space = TRUE;
          only_digits = TRUE;
          for (ptr = gbqual->val, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
            if (IS_WHITESP (ch)) {
              no_white_space = FALSE;
            }
            if (! IS_DIGIT (ch)) {
              only_digits = FALSE;
            }
          }
          if (only_digits || (! no_white_space)) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "Illegal value for qualifier %s", gbqual->qual);
          }
        } else if (val == GBQUAL_replace) {
          bsp = BioseqFindFromSeqLoc (sfp->location);
          if (bsp != NULL) {
            if (ISA_na (bsp->mol)) {
              if (StringICmp (key, "variation") == 0) {
                just_nuc_letters = TRUE;
                for (ptr = gbqual->val, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
                  if (StringChr ("acgt", ch) == NULL) {
                    just_nuc_letters = FALSE;
                  }
                }
                if (!just_nuc_letters) {
                  ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue,
                            "%s is not a legal value for qualifier %s - should only be composed of acgt unambiguous nucleotide bases",
                            gbqual->val, gbqual->qual);
                }
              } else {
                just_nuc_letters = TRUE;
                for (ptr = gbqual->val, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
                  if (StringChr ("acgtmrwsykvhdbn", ch) == NULL) {
                    just_nuc_letters = FALSE;
                  }
                }
                if (!just_nuc_letters) {
                  ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue,
                            "%s is not a legal value for qualifier %s - should only be composed of acgtmrwsykvhdbn nucleotide bases",
                            gbqual->val, gbqual->qual);
                }
              }
            } else if (ISA_aa (bsp->mol)) {
              just_prt_letters = TRUE;
              for (ptr = gbqual->val, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
                if (StringChr ("acdefghiklmnpqrstuvwy*", ch) == NULL) {
                  just_prt_letters = FALSE;
                }
              }
              if (!just_prt_letters) {
                ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue,
                          "%s is not a legal value for qualifier %s - should only be composed of acdefghiklmnpqrstuvwy* amino acids",
                          gbqual->val, gbqual->qual);
              }
            }
            slp = sfp->location;
            fuzz = NULL;
            if (slp != NULL && slp->choice == SEQLOC_PNT) {
              spp = (SeqPntPtr) slp->data.ptrvalue;
              if (spp != NULL) {
                fuzz = spp->fuzz;
              }
            }
            if (slp != NULL && StringLen (gbqual->val) == SeqLocLen (slp) && fuzz == NULL) {
              tmp = GetSequenceByFeature (sfp);
              if (tmp != NULL) {
                if (StringICmp (tmp, gbqual->val) == 0) {
                  ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_SuspiciousQualifierValue, "/replace already matches underlying sequence (%s)", gbqual->val);
                }
                MemFree (tmp);
              }
            }
          }
        } else if (val == GBQUAL_cons_splice) {
          found = FALSE;
          for (i = 0; legal_cons_splice_strings[i] != NULL; i++) {
            if (StringICmp (gbqual->val, legal_cons_splice_strings[i]) == 0) {
              found = TRUE;
              break;
            }
          }
          if (!found) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "%s is not a legal value for qualifier %s", gbqual->val, gbqual->qual);
          }
        } else if (val == GBQUAL_mobile_element_type) {
          found = FALSE;
          str = NULL;
          for (i = 0; legal_mobile_element_strings[i] != NULL; i++) {
            ptr = legal_mobile_element_strings[i];
            len = StringLen (ptr);
            if (StringNICmp (gbqual->val, ptr, len) == 0) {
              found = TRUE;
              str = gbqual->val + len;
              break;
            }
          }
          if (found) {
            if (StringDoesHaveText (str) && (str [0] != ':' || str [1] == '\0')) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_InvalidQualifierValue, "%s is not a legal value for qualifier %s", gbqual->val, gbqual->qual);
            } else if (StringNICmp (gbqual->val, "other", 5) == 0) {
              if (str [0] != ':' || str [1] == '\0') {
                ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_InvalidQualifierValue, "%s is not a legal value for qualifier %s", gbqual->val, gbqual->qual);
              }
            }
          }
          if (!found) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_InvalidQualifierValue, "%s is not a legal value for qualifier %s", gbqual->val, gbqual->qual);
          }
        } else if (val == GBQUAL_frequency) {
          if (StringCmp (gbqual->val, "1") == 0 || StringCmp (gbqual->val, "1.0") == 0 || StringCmp (gbqual->val, "1.00") == 0) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_InvalidQualifierValue, "%s is a suspicious value for qualifier %s", gbqual->val, gbqual->qual);
          }
        } else if (val == GBQUAL_compare) {
          multi_compare = FALSE;
          ptr = gbqual->val;
          ch = *ptr;
          if (ch == '(') {
            multi_compare = TRUE;
          }
          if (! multi_compare) {
            adv = ValidateAccnDotVer (gbqual->val);
            if (adv == -5) {
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "%s accession missing version for qualifier %s", gbqual->val, gbqual->qual);
            } else if (adv == -6) {
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "%s accession has bad version for qualifier %s", gbqual->val, gbqual->qual);
            } else if (adv != 0) {
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "%s is not a legal accession for qualifier %s", gbqual->val, gbqual->qual);
            } else if (StringChr (gbqual->val, '_') != NULL) {
              if (vsp->is_insd_in_sep) {
                ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "RefSeq accession %s cannot be used for qualifier %s", gbqual->val, gbqual->qual);
              }
            }
          }
        }
      }
    }
  }
  if (index != -1 && ParFlat_GBFeat[index].mand_num > 0) {
    for (i = 0; i < ParFlat_GBFeat[index].mand_num; i++) {
      found = FALSE;
      qual = ParFlat_GBFeat[index].mand_qual[i];
      for (gbqual = sfp->qual; gbqual != NULL; gbqual = gbqual->next) {
        val = GBQualNameValid (gbqual->qual);
        if (qual == val) {
          found = TRUE;
          break;
        }
      }
      if (!found) {
        if (qual == GBQUAL_citation) {
          if (sfp->cit != NULL) {
            found = TRUE;
          } else if (! StringHasNoText (sfp->comment)) {
            /* RefSeq allows conflict with accession in comment instead of sfp->cit */
            if (StringICmp (key, "conflict") == 0) {
              bsp = BioseqFindFromSeqLoc (sfp->location);
              if (bsp != NULL) {
                for (sip = bsp->id; sip != NULL; sip = sip->next) {
                  if (sip->choice == SEQID_OTHER) {
                    found = TRUE;
                  }
                }
              }
            }
          }
        }
      }
      if (!found) {
        if (StringICmp (key, "conflict") == 0 || StringICmp (key, "old_sequence") == 0) {
          /* compare qualifier can now substitute for citation qualifier for conflict and old_sequence */
          for (gbqual = sfp->qual; gbqual != NULL; gbqual = gbqual->next) {
            if (StringICmp (gbqual->qual, "compare") == 0 && StringDoesHaveText (gbqual->val)) {
              found = TRUE;
            }
          }
        }
      }
      if (!found) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_MissingQualOnImpFeat, "Missing qualifier %s for feature %s", ParFlat_GBQual_names[qual].name, key);
      }
    }
  }
  if (!StringHasNoText (ifp->loc)) {
    slp = sfp->location;
    if (StringStr (ifp->loc, "one-of") != NULL) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_ImpFeatBadLoc, "ImpFeat loc %s has obsolete 'one-of' text for feature %s", ifp->loc, key);
    } else if (slp != NULL && slp->choice == SEQLOC_INT) {
      sint = (SeqIntPtr) slp->data.ptrvalue;
      if (sint != NULL && sint->strand != Seq_strand_minus) {
        sprintf (range, "%ld..%ld", (long) (sint->from + 1), (long) (sint->to + 1));
        if (StringCmp (ifp->loc, range) != 0) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_ImpFeatBadLoc, "ImpFeat loc %s does not equal feature location %s for feature %s", ifp->loc, range, key);
        }
      }
    }
  }
  MemFree (key);
}

static void ValidateNonImpFeat (ValidStructPtr vsp, GatherContextPtr gcp, SeqFeatPtr sfp)
{
  Int2       adv;
  BioseqPtr  bsp;
  Char       ch;
  Boolean    failed;
  Boolean    found;
  GBQualPtr  gbqual;
  Int2       i;
  Int2       index;
  CharPtr    key;
  Boolean    multi_compare;
  Boolean    no_white_space;
  Boolean    only_digits;
  CharPtr    ptr;
  Int2       qual;
  RNAGenPtr  rgp;
  RnaRefPtr  rrp;
  ErrSev     sev;
  SeqIdPtr   sip;
  CharPtr    str;
  CharPtr    tmp;
  Int2       val;
  Int4       qvalid;

  if (vsp == NULL || gcp == NULL || sfp == NULL)
    return;
  key = StringSaveNoNull (FeatDefTypeLabel (sfp));
  if (StringCmp (key, "Gene") == 0) {
    *key = 'g';
  } else if (StringCmp (key, "preRNA") == 0) {
    key = MemFree (key);
    key = StringSave ("precursor_RNA");
  }
  index = GBFeatKeyNameValid (&key, FALSE);
  for (gbqual = sfp->qual; gbqual != NULL; gbqual = gbqual->next) {
    qvalid = IsQualValidForFeature (gbqual, sfp);
    if (qvalid == 0) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_WrongQualOnFeature, "Wrong qualifier %s for feature %s", gbqual->qual, key);
    } else if (qvalid > 1) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_WrongQualOnFeature, sWrongQualReasons[qvalid - 2]);
    }
    
    if (StringCmp (gbqual->qual, "gsdb_id") == 0) {
      continue;
    }
    val = GBQualNameValid (gbqual->qual);
    if (val == -1) {
      if (gbqual->qual != NULL) {
        if (sfp->data.choice == SEQFEAT_GENE) {
          if (StringCmp (gbqual->qual, "gen_map") == 0) continue;
          if (StringCmp (gbqual->qual, "cyt_map") == 0) continue;
          if (StringCmp (gbqual->qual, "rad_map") == 0) continue;
        }
        if (sfp->data.choice == SEQFEAT_CDREGION) {
          if (StringCmp (gbqual->qual, "orig_transcript_id") == 0) continue;
        }
        if (sfp->data.choice == SEQFEAT_RNA) {
          if (StringCmp (gbqual->qual, "orig_protein_id") == 0) continue;
          if (StringCmp (gbqual->qual, "orig_transcript_id") == 0) continue;
        }
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UnknownFeatureQual, "Unknown qualifier %s", gbqual->qual);
      } else {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UnknownFeatureQual, "NULL qualifier");
      }
    } else if (index != -1) {
      if (gbqual->val != NULL) {
        if (val == GBQUAL_rpt_type) {
          failed = FALSE;
          tmp = StringSave (gbqual->val);
          str = tmp;
          if (*str == '(') {
            str++;
          }
          while (!StringHasNoText (str)) {
            ptr = StringChr (str, ',');
            if (ptr == NULL) {
              ptr = StringChr (str, ')');
            }
            if (ptr != NULL) {
              *ptr = '\0';
              ptr++;
            }
            found = FALSE;
            for (i = 0; legal_repeat_types[i] != NULL; i++) {
              if (StringICmp (str, legal_repeat_types[i]) == 0) {
                found = TRUE;
                break;
              }
            }
            if (!found) {
              failed = TRUE;
            }
            str = ptr;
          }
          MemFree (tmp);
          if (failed) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "%s is not a legal value for qualifier %s", gbqual->val, gbqual->qual);
          }
        } else if (val == GBQUAL_rpt_unit || val == GBQUAL_rpt_unit_range || val == GBQUAL_rpt_unit_seq) {
          ValidateRptUnit (vsp, gcp, sfp, gbqual, val, key);
        } else if (val == GBQUAL_label) {
          no_white_space = TRUE;
          only_digits = TRUE;
          for (ptr = gbqual->val, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
            if (IS_WHITESP (ch)) {
              no_white_space = FALSE;
            }
            if (! IS_DIGIT (ch)) {
              only_digits = FALSE;
            }
          }
          if (only_digits || (! no_white_space)) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "Illegal value for qualifier %s", gbqual->qual);
          }
        } else if (val == GBQUAL_cons_splice) {
          found = FALSE;
          for (i = 0; legal_cons_splice_strings[i] != NULL; i++) {
            if (StringICmp (gbqual->val, legal_cons_splice_strings[i]) == 0) {
              found = TRUE;
              break;
            }
          }
          if (!found) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "%s is not a legal value for qualifier %s", gbqual->val, gbqual->qual);
          }
        } else if (val == GBQUAL_compare) {
          multi_compare = FALSE;
          ptr = gbqual->val;
          ch = *ptr;
          if (ch == '(') {
            multi_compare = TRUE;
          }
          if (! multi_compare) {
            adv = ValidateAccnDotVer (gbqual->val);
            if (adv == -5) {
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "%s accession missing version for qualifier %s", gbqual->val, gbqual->qual);
            } else if (adv == -6) {
             ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "%s accession has bad version for qualifier %s", gbqual->val, gbqual->qual);
            } else if (adv != 0) {
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "%s is not a legal accession for qualifier %s", gbqual->val, gbqual->qual);
            } else if (StringChr (gbqual->val, '_') != NULL) {
              if (vsp->is_insd_in_sep) {
                ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "RefSeq accession %s cannot be used for qualifier %s", gbqual->val, gbqual->qual);
              }
            }
          }
        }
      }
    }
  }
  if (index != -1 && ParFlat_GBFeat[index].mand_num > 0) {
    for (i = 0; i < ParFlat_GBFeat[index].mand_num; i++) {
      sev = SEV_WARNING;
      found = FALSE;
      qual = ParFlat_GBFeat[index].mand_qual[i];
      for (gbqual = sfp->qual; gbqual != NULL; gbqual = gbqual->next) {
        val = GBQualNameValid (gbqual->qual);
        if (qual == val) {
          found = TRUE;
          break;
        }
      }
      if (!found) {
        if (qual == GBQUAL_citation) {
          if (sfp->cit != NULL) {
            found = TRUE;
          } else if (! StringHasNoText (sfp->comment)) {
            /* RefSeq allows conflict with accession in comment instead of sfp->cit */
            if (StringICmp (key, "conflict") == 0) {
              bsp = BioseqFindFromSeqLoc (sfp->location);
              if (bsp != NULL) {
                for (sip = bsp->id; sip != NULL; sip = sip->next) {
                  if (sip->choice == SEQID_OTHER) {
                    found = TRUE;
                  }
                }
              }
            }
          }
        }
      }
      if (!found) {
        if (StringICmp (key, "conflict") == 0 || StringICmp (key, "old_sequence") == 0) {
          /* compare qualifier can now substitute for citation qualifier for conflict and old_sequence */
          for (gbqual = sfp->qual; gbqual != NULL; gbqual = gbqual->next) {
            if (StringICmp (gbqual->qual, "compare") == 0 && StringDoesHaveText (gbqual->val)) {
              found = TRUE;
            }
          }
        }
      }
      if (!found) {
        if (qual == GBQUAL_ncRNA_class) {
          sev = SEV_ERROR;
          if (sfp->data.choice == SEQFEAT_RNA) {
            rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
            if (rrp != NULL) {
              if (rrp->ext.choice == 3) {
                rgp = (RNAGenPtr) rrp->ext.value.ptrvalue;
                if (rgp != NULL) {
                  if (StringDoesHaveText (rgp->_class)) {
                    found = TRUE;
                  }
                }
              }
            }
          }
        }
      }
      if (!found) {
        ValidErr (vsp, sev, ERR_SEQ_FEAT_MissingQualOnFeature,
                  "Missing qualifier %s for feature %s", ParFlat_GBQual_names[qual].name, key);
      }
    }
  }
  if (StringICmp (key, "mat_peptide") == 0 ||
      StringICmp (key, "sig_peptide") == 0 ||
      StringICmp (key, "transit_peptide") == 0) {
    bsp = BioseqFindFromSeqLoc (sfp->location);
    if (bsp != NULL) {
      if (ISA_na (bsp->mol)) {
        sev = SEV_WARNING;
        if (vsp->is_refseq_in_sep) {
          sev = SEV_ERROR;
        }
        ValidErr (vsp, sev, ERR_SEQ_FEAT_InvalidForType, "Peptide processing feature should be remapped to the appropriate protein bioseq");
        CheckPeptideOnCodonBoundary (vsp, gcp, sfp, key);
      }
    }
  } else if (StringICmp (key, "preprotein") == 0 ||
      StringICmp (key, "proprotein") == 0) {
    bsp = BioseqFindFromSeqLoc (sfp->location);
    if (bsp != NULL) {
      if (ISA_na (bsp->mol)) {
        sev = SEV_WARNING;
        if (vsp->is_refseq_in_sep) {
          sev = SEV_ERROR;
        }
        ValidErr (vsp, sev, ERR_SEQ_FEAT_InvalidForType, "Peptide processing feature should be remapped to the appropriate protein bioseq");
      }
    }
  }
  MemFree (key);
}

/* PartialAtSpliceSiteOrGap uses code taken from SpliceCheckEx */
static Boolean PartialAtSpliceSiteOrGap (ValidStructPtr vsp, SeqLocPtr head, Uint2 slpTag, BoolPtr isgapP, BoolPtr badseqP)
{
  BioseqPtr       bsp;
  Int2            residue1, residue2;
  Boolean         rsult = FALSE;
  SeqIdPtr        sip;
  SeqLocPtr       slp = NULL, first = NULL, last = NULL;
  /*
  SeqPortPtr      spp = NULL;
  */
  Uint1           strand;
  Int4            strt, stp, donor, acceptor, len;
  StreamCache     sc;
  SeqInt          sint;
  ValNode         vn;

  if (isgapP != NULL) {
    *isgapP = FALSE;
  }
  if (badseqP != NULL) {
    *badseqP = FALSE;
  }
  if (slpTag != SLP_NOSTART && slpTag != SLP_NOSTOP)
    return FALSE;
  while ((slp = SeqLocFindPart (head, slp, EQUIV_IS_ONE)) != NULL) {
    if (first == NULL) {
      first = slp;
    }
    last = slp;
  }
  if (first == NULL)
    return FALSE;

  strand = SeqLocStrand (first);
  if (SeqLocStrand (last) != strand)
    return FALSE;

  if (slpTag == SLP_NOSTART) {
    slp = first;
  } else {
    slp = last;
  }
  sip = SeqLocId (slp);
  if (sip == NULL)
    return FALSE;

  bsp = NULL;
  if (sip != NULL && (sip->choice != SEQID_GI || sip->data.intvalue > 0)) {
    bsp = BioseqLockById (sip);
  }
  if (bsp == NULL)
    return FALSE;
  len = bsp->length;

  acceptor = SeqLocStart (slp);
  donor = SeqLocStop (slp);

  if (acceptor < 0 || acceptor >= len || donor < 0 || donor >= len) {
    /*
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_Range,
              "Unable to check splice consensus because feature outside range of sequence");
    */
    BioseqUnlock (bsp);
    return FALSE;
  }

  if (strand != Seq_strand_minus) {
    if (! StreamCacheSetup (bsp, NULL, EXPAND_GAPS_TO_DASHES, &sc)) {
      BioseqUnlock (bsp);
      return FALSE;
    }
  } else {
    sint.from = 0;
    sint.to = len - 1;
    sint.strand = strand;
    sint.id = sip;
    vn.choice = SEQLOC_INT;
    vn.data.ptrvalue = (Pointer) &sint;
    vn.next = NULL;
    if (! StreamCacheSetup (NULL, &vn, EXPAND_GAPS_TO_DASHES, &sc)) {
      BioseqUnlock (bsp);
      return FALSE;
    }
  }
  /* spp = SeqPortNew (bsp, 0, -1, strand, Seq_code_ncbi4na); */
  BioseqUnlock (bsp);
  /*
  if (spp == NULL)
    return FALSE;
  */

  if (strand != Seq_strand_minus) {
    strt = acceptor;
    stp = donor;
  } else {
    strt = donor;
    donor = acceptor;
    acceptor = strt;
    stp = len - donor - 1;
    strt = len - acceptor - 1;
  }

  if (slpTag == SLP_NOSTOP && stp < len - 2) {
    StreamCacheSetPosition (&sc, stp + 1);
    residue1 = StreamCacheGetResidue (&sc);
    residue2 = StreamCacheGetResidue (&sc);
    /*
    SeqPortSeek (spp, (stp + 1), SEEK_SET);
    residue1 = SeqPortGetResidue (spp);
    residue2 = SeqPortGetResidue (spp);
    */
    if (residue1 == '-' && residue2 == '-') {
      if (isgapP != NULL) {
        *isgapP = TRUE;
      }
      rsult = TRUE;
    } else if (IS_residue (residue1) && IS_residue (residue2)) {
      if ((residue1 == 'G') && (residue2  == 'T')) {
        rsult = TRUE;
      } else if ((residue1 == 'G') && (residue2 == 'C')) {
        rsult = TRUE;
      }
    } else if (badseqP != NULL) {
      *badseqP = TRUE;
    }
  } else if (slpTag == SLP_NOSTART && strt > 1) {
    StreamCacheSetPosition (&sc, strt - 2);
    residue1 = StreamCacheGetResidue (&sc);
    residue2 = StreamCacheGetResidue (&sc);
    /*
    SeqPortSeek (spp, (strt - 2), SEEK_SET);
    residue1 = SeqPortGetResidue (spp);
    residue2 = SeqPortGetResidue (spp);
    */
    if (residue1 == '-' && residue2 == '-') {
      if (isgapP != NULL) {
        *isgapP = TRUE;
      }
      rsult = TRUE;
    } else if (IS_residue (residue1) && IS_residue (residue2) && IS_ALPHA ((Char) residue1) && IS_ALPHA ((Char) residue2)) {
      if ((residue1 == 'A') && (residue2 == 'G')) {
        rsult = TRUE;
      }
    } else if (badseqP != NULL) {
      *badseqP = TRUE;
    }
  }

  /* spp = SeqPortFree (spp); */
  return rsult;
}

static Boolean PartialAtGapOrNs (ValidStructPtr vsp, SeqLocPtr head, Uint2 slpTag)

{
  BioseqPtr       bsp;
  Int2            residue;
  Boolean         rsult = FALSE;
  SeqIdPtr        sip;
  SeqLocPtr       slp = NULL, first = NULL, last = NULL;
  Uint1           strand;
  Int4            strt, stp, donor, acceptor, len;
  StreamCache     sc;
  SeqInt          sint;
  ValNode         vn;

  if (slpTag != SLP_NOSTART && slpTag != SLP_NOSTOP) return FALSE;

  while ((slp = SeqLocFindPart (head, slp, EQUIV_IS_ONE)) != NULL) {
    if (first == NULL) {
      first = slp;
    }
    last = slp;
  }
  if (first == NULL) return FALSE;

  strand = SeqLocStrand (first);
  if (SeqLocStrand (last) != strand) return FALSE;

  if (slpTag == SLP_NOSTART) {
    slp = first;
  } else {
    slp = last;
  }
  sip = SeqLocId (slp);
  if (sip == NULL) return FALSE;

  bsp = NULL;
  if (sip != NULL && (sip->choice != SEQID_GI || sip->data.intvalue > 0)) {
    bsp = BioseqLockById (sip);
  }
  if (bsp == NULL) return FALSE;
  len = bsp->length;

  acceptor = SeqLocStart (slp);
  donor = SeqLocStop (slp);

  if (acceptor < 0 || acceptor >= len || donor < 0 || donor >= len) {
    BioseqUnlock (bsp);
    return FALSE;
  }

  if (strand != Seq_strand_minus) {
    if (! StreamCacheSetup (bsp, NULL, EXPAND_GAPS_TO_DASHES, &sc)) {
      BioseqUnlock (bsp);
      return FALSE;
    }
  } else {
    sint.from = 0;
    sint.to = len - 1;
    sint.strand = strand;
    sint.id = sip;
    vn.choice = SEQLOC_INT;
    vn.data.ptrvalue = (Pointer) &sint;
    vn.next = NULL;
    if (! StreamCacheSetup (NULL, &vn, EXPAND_GAPS_TO_DASHES, &sc)) {
      BioseqUnlock (bsp);
      return FALSE;
    }
  }
  BioseqUnlock (bsp);

  if (strand != Seq_strand_minus) {
    strt = acceptor;
    stp = donor;
  } else {
    strt = donor;
    donor = acceptor;
    acceptor = strt;
    stp = len - donor - 1;
    strt = len - acceptor - 1;
  }

  if (slpTag == SLP_NOSTOP && stp < len - 2) {
    StreamCacheSetPosition (&sc, stp + 1);
    residue = StreamCacheGetResidue (&sc);
    if (residue == '-' || residue == 'N') {
      rsult = TRUE;
    }
  } else if (slpTag == SLP_NOSTART && strt > 1) {
    StreamCacheSetPosition (&sc, strt - 1);
    residue = StreamCacheGetResidue (&sc);
    if (residue == '-' || residue == 'N') {
      rsult = TRUE;
    }
  }

  return rsult;
}


#if 0
static void CheckTrnaCodons (ValidStructPtr vsp, GatherContextPtr gcp, SeqFeatPtr sfp, tRNAPtr trp)
{
  Uint1           aa = 0;
  BioseqPtr       bsp;
  Int2            code = 0;
  CharPtr         codes = NULL;
  Uint1           codon [4];
  Uint1           from;
  CharPtr         gen_code_name = NULL;
  GeneticCodePtr  gncp;
  Uint2           idx;
  Int2            j;
  Int2            k;
  ErrSev          sev = SEV_ERROR;
  SeqMapTablePtr  smtp;
  Uint1           taa;
  CharPtr         three_letter_aa = NULL;
  ValNodePtr      vnp;

  if (vsp == NULL || gcp == NULL || sfp == NULL || trp == NULL)
    return;

  aa = 0;
  if (trp->aatype == 2) {
    aa = trp->aa;
  } else {
    from = 0;
    switch (trp->aatype) {
    case 0:
      from = 0;
      break;
    case 1:
      from = Seq_code_iupacaa;
      break;
    case 2:
      from = Seq_code_ncbieaa;
      break;
    case 3:
      from = Seq_code_ncbi8aa;
      break;
    case 4:
      from = Seq_code_ncbistdaa;
      break;
    default:
      break;
    }
    smtp = SeqMapTableFind (Seq_code_ncbieaa, from);
    if (smtp != NULL) {
      aa = SeqMapTableConvert (smtp, trp->aa);
    }
  }

  for (j = 0; j < 6; j++) {
    if (trp->codon[j] < 64) {
      if (codes == NULL) {
        bsp = GetBioseqGivenSeqLoc (sfp->location, gcp->entityID);
        /*
        sep = GetBestTopParentForData (gcp->entityID, bsp);
        code = SeqEntryToGeneticCode (sep, NULL, NULL, 0);
        */
        BioseqToGeneticCode (bsp, &code, NULL, NULL, NULL, 0, NULL);
        gncp = GeneticCodeFind (code, NULL);
        if (gncp == NULL) {
          gncp = GeneticCodeFind (1, NULL);
          code = 1;
        }
        if (gncp == NULL)
          return;
        for (vnp = (ValNodePtr) gncp->data.ptrvalue; vnp != NULL; vnp = vnp->next) {
          if (vnp->choice == 3) {
            codes = (CharPtr) vnp->data.ptrvalue;
          }
        }
      }
      if (codes == NULL)
        return;
      taa = codes[trp->codon[j]];
      if (aa > 0 && aa != 255) {
        if (taa != aa) {
          if (aa == 'U' || aa == 'O') {
            sev = SEV_WARNING;
          }
          if (aa == 'U' && taa == '*' && trp->codon [j] == 14) {
            /* selenocysteine normally uses TGA (14), so ignore without requiring exception in record */
          } else if (aa == 'O' && taa == '*' && trp->codon [j] == 11) {
            /* pyrrolysine normally uses TAG (11) in archaebacteria, so ignore without requiring exception in record */

            /* TAA (10) is not yet known to be used for an exceptional amino acid */
          } else if (StringISearch (sfp->except_text, "modified codon recognition") == NULL) {
            codon [0] = '\0';
            if (CodonForIndex (trp->codon [j], Seq_code_iupacna, codon)) {
              for (k = 0; k < 3; k++) {
                if (codon [k] == 'T') {
                  codon [k] = 'U';
                }
              }
              codon [3] = '\0';
            } else {
              StringCpy ((CharPtr) codon, "?");
            }
            three_letter_aa = Get3LetterSymbol (NULL, Seq_code_ncbieaa, NULL, aa);
            if (StringHasNoText (three_letter_aa)) {
              three_letter_aa = "?";
            }
            for (vnp = genetic_code_name_list; vnp != NULL; vnp = vnp->next) {
              if (vnp->choice != (Uint1) code) continue;
              gen_code_name = (CharPtr) vnp->data.ptrvalue;
              break;
            }
            if (StringHasNoText (gen_code_name)) {
              gen_code_name = "?";
            }
            ValidErr (vsp, sev, ERR_SEQ_FEAT_TrnaCodonWrong,
                      "Codon recognized by tRNA (%s) does not match amino acid (%c/%s) specified by genetic code (%d/%s)",
                      (char *) codon, (char) aa, (char *) three_letter_aa, (int) code, (char *) gen_code_name);
          }
        }
      }
    } else if (trp->codon [j] < 255) {
      ValidErr (vsp, sev, ERR_SEQ_FEAT_BadTrnaCodon, "tRNA codon value %d is greater than maximum 63", (int) trp->codon [j]);
    }
  }

  if (sfp->pseudo) return;

  if (aa > 0 && aa != 255) {
    /* - no gaps now that O and J are added
    if (aa <= 74) {
      shift = 0;
    } else if (aa > 79) {
      shift = 2;
    } else {
      shift = 1;
    }
    */
    if (aa != '*') {
      idx = aa - (64 /* + shift */);
    } else {
      idx = 25; /* termination */
    }
    if (idx > 0 && idx < 28) {
      /* valid trna amino acid */
    } else {
      ValidErr (vsp, sev, ERR_SEQ_FEAT_BadTrnaAA, "Invalid tRNA amino acid");
    }
  } else {
    ValidErr (vsp, sev, ERR_SEQ_FEAT_BadTrnaAA, "Missing tRNA amino acid");
  }
}
#endif

static Boolean TwoListsHaveCommonItem (
  ValNodePtr list1,
  ValNodePtr list2
)

{
  CharPtr     str1;
  CharPtr     str2;
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  for (vnp1 = list1; vnp1 != NULL; vnp1 = vnp1->next) {
    str1 = (CharPtr) vnp1->data.ptrvalue;
    if (StringHasNoText (str1)) continue;
    for (vnp2 = list2; vnp2 != NULL; vnp2 = vnp2->next) {
      str2 = (CharPtr) vnp2->data.ptrvalue;
      if (StringHasNoText (str2)) continue;
      if (StringICmp (str1, str2) == 0) return TRUE;
    }
  }

  return FALSE;
}


static void CheckTrnaCodons (
  ValidStructPtr vsp,
  GatherContextPtr gcp,
  SeqFeatPtr sfp,
  tRNAPtr trp
)

{
  Uint1           aa = 0;
  Uint1           anticodon [4];
  Char            ch;
  Int2            code = 0;
  CharPtr         codes = NULL;
  Uint1           codon [4];
  CharPtr         complementBase = " TVGH  CD  M KN   YSAABW R ";
  CharPtr         gen_code_name = NULL;
  Int2            i;
  Uint2           idx;
  Uint1           index;
  Int2            j;
  Int2            k;
  Uint1           letterToComp [256];
  Char            lttr;
  Boolean         okay;
  ValNodePtr      possibles = NULL;
  ValNodePtr      recognizes = NULL;
  StreamCache     sc;
  ErrSev          sev = SEV_ERROR;
  SeqLocPtr       slp;
  CharPtr         str;
  Uint1           taa;
  CharPtr         three_letter_aa = NULL;
  ValNodePtr      vnp;
  CharPtr         wobble = NULL;
  Boolean         rna_editing = FALSE;

  if (vsp == NULL || gcp == NULL || sfp == NULL || trp == NULL) return;

  anticodon [0] = '\0';

  /* extract indicated amino acid */

  aa = GetAaFromtRNA (trp);

  three_letter_aa = Get3LetterSymbol (NULL, Seq_code_ncbieaa, NULL, aa);
  if (StringHasNoText (three_letter_aa)) {
    three_letter_aa = "?";
  }

  /* find genetic code table */
  codes = GetCodesFortRNA(sfp, &code);

  if (codes == NULL) return;

  for (vnp = genetic_code_name_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice != (Uint1) code) continue;
    gen_code_name = (CharPtr) vnp->data.ptrvalue;
    break;
  }
  if (StringHasNoText (gen_code_name)) {
    gen_code_name = "?";
  }

  /* set up nucleotide complementation lookup table */

  for (i = 0; i < 256; i++) {
    letterToComp [i] = '\0';
  }
  for (ch = 'A', i = 1; ch <= 'Z'; ch++, i++) {
    lttr = complementBase [i];
    if (lttr != ' ') {
      letterToComp [(int) (Uint1) ch] = lttr;
    }
  }
  for (ch = 'a', i = 1; ch <= 'z'; ch++, i++) {
    lttr = complementBase [i];
    if (lttr != ' ') {
      letterToComp [(int) (Uint1) ch] = lttr;
    }
  }

  if (StringCmp (sfp->except_text, "RNA editing") == 0) {
    rna_editing = TRUE;
  }

  /* loop through codon_recognized array */

  for (j = 0; j < 6; j++) {
    index = (Uint1) trp->codon [j];

    if (index == 255) continue;

    if (index >= 64) {
      ValidErr (vsp, sev, ERR_SEQ_FEAT_BadTrnaCodon,
                "tRNA codon value %d is greater than maximum 63",
                (int) (index));
      continue;
    }
    if (rna_editing) {
      continue;
    }

    taa = codes [index];

    codon [0] = '\0';
    if (CodonForIndex (index, Seq_code_iupacna, codon)) {
      /*
      for (k = 0; k < 3; k++) {
        if (codon [k] == 'T') {
          codon [k] = 'U';
        }
      }
      */
      codon [3] = '\0';
    } else {
      StringCpy ((CharPtr) codon, "?");
    }

    /* save codon recognized and translated amino acid for anticodon reality check */

    ValNodeCopyStr (&recognizes, taa, (CharPtr) codon);

    if (aa == 0 || aa == 255) continue;

    /* only report if encoded amino acid does not match indicated amino acid */

    if (taa == aa) continue;

    if (aa == 'U' || aa == 'O') {
      sev = SEV_WARNING;
    }

    /* selenocysteine normally uses TGA (14), so ignore without requiring exception in record */
    if (aa == 'U' && taa == '*' && index == 14) continue;

    /* pyrrolysine normally uses TAG (11) in archaebacteria, ignore without requiring exception */
    if (aa == 'O' && taa == '*' && index == 11) continue;

    /* TAA (10) is not yet known to be used for an exceptional amino acid, but the night is young */

    /* ignore if modified codon recognition exception is present */
    if (StringISearch (sfp->except_text, "modified codon recognition") != NULL) continue;

    for (k = 0; k < 3; k++) {
      if (codon [k] == 'T') {
        codon [k] = 'U';
      }
    }
    codon [3] = '\0';

    ValidErr (vsp, sev, ERR_SEQ_FEAT_TrnaCodonWrong,
              "Codon recognized by tRNA (%s) does not match amino acid (%c/%s) specified by genetic code (%d/%s)",
              (char *) codon, (char) aa, (char *) three_letter_aa, (int) code, (char *) gen_code_name);
  }

  /* see if anticodon is compatible with codons recognized and amino acid */

  slp = trp->anticodon;
  if (slp != NULL && SeqLocLen (slp) == 3) {

    /* read sequence under anticodon */

    if (StreamCacheSetup (NULL, slp, 0, &sc)) {
      for (i = 0; i < 3; i++) {
        ch = (Char) StreamCacheGetResidue (&sc);
        anticodon [i] = ch;
      }
      anticodon [3] = '\0';

      /* reverse complement non-wobble bases */

      codon [0] = letterToComp [(int) (Uint1) anticodon [2]];
      codon [1] = letterToComp [(int) (Uint1) anticodon [1]];
      codon [3] = '\0';

      /* expand wobble base to known binding partners */

      ch = anticodon [0];
      switch (ch) {
        case 'A' :
          wobble = "ACT";
          break;
        case 'C' :
          wobble = "G";
          break;
        case 'G' :
          wobble = "CT";
          break;
        case 'T' :
          wobble = "AG";
          break;
        default :
          break;
      }

      if (wobble != NULL) {
        for (i = 0; wobble [i] != '\0'; i++) {
          codon [2] = wobble [i];
          index = IndexForCodon (codon, Seq_code_iupacna);
          if (index < 64) {
            taa = codes [index];

            /* save possible codon recognized and translated amino acid */

            ValNodeCopyStr (&possibles, taa, (CharPtr) codon);
          }
        }
      }
    }
  }

  for (k = 0; k < 3; k++) {
    if (anticodon [k] == 'T') {
      anticodon [k] = 'U';
    }
  }
  anticodon [3] = '\0';

  if (StringHasNoText ((CharPtr) anticodon)) {
    StringCpy ((CharPtr) anticodon, "?");
  }

  /* check that codons predicted from anticodon can transfer indicated amino acid */

  if (possibles != NULL) {
    okay = FALSE;
    for (vnp = possibles; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      taa = vnp->choice;
      if (taa == aa) {
        okay = TRUE;
      }
    }
    if (! okay) {
      if (aa == 'U' && StringCmp ((CharPtr) anticodon, "UCA") == 0) {
        /* ignore TGA codon for selenocysteine */
      } else if (aa == 'O' && StringCmp ((CharPtr) anticodon, "CUA") == 0) {
        /* ignore TAG codon for pyrrolysine */
      } else if (StringISearch (sfp->except_text, "modified codon recognition") == NULL &&
                 StringISearch (sfp->except_text, "RNA editing") == NULL) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadAnticodonAA,
                  "Codons predicted from anticodon (%s) cannot produce amino acid (%c/%s)",
                  (char *) anticodon, (char) aa, (char *) three_letter_aa);
      }
    }
  }

  /* check that codons recognized match codons predicted from anticodon */

  if (recognizes != NULL && possibles != NULL) {
    okay = FALSE;
    if (TwoListsHaveCommonItem (recognizes, possibles)) {
      okay = TRUE;
    }
    if (! okay) {
      if (StringISearch (sfp->except_text, "RNA editing") == NULL) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadAnticodonCodon,
                  "Codon recognized cannot be produced from anticodon (%s)",
                  (char*) anticodon);
      }
    }
  }

  ValNodeFreeData (recognizes);
  ValNodeFreeData (possibles);

  if (sfp->pseudo) return;

  if (aa == 0 || aa == 255) {
    ValidErr (vsp, sev, ERR_SEQ_FEAT_BadTrnaAA, "Missing tRNA amino acid");
    return;
  }

  /* verify that legal amino acid is indicated */

  /* - no gaps now that O and J are added
  if (aa <= 74) {
    shift = 0;
  } else if (aa > 79) {
    shift = 2;
  } else {
    shift = 1;
  }
  */
  if (aa != '*') {
    idx = aa - (64 /* + shift */);
  } else {
    idx = 25; /* termination */
  }
  if (idx == 0 || idx >= 28) {
    ValidErr (vsp, sev, ERR_SEQ_FEAT_BadTrnaAA, "Invalid tRNA amino acid");
  }
}

static void CheckRnaProductType (ValidStructPtr vsp, GatherContextPtr gcp, SeqFeatPtr sfp, RnaRefPtr rrp)

{
  BioseqPtr          bsp;
  SeqMgrDescContext  context;
  MolInfoPtr         mip;
  SeqDescrPtr        sdp;

  if (vsp == NULL || gcp == NULL || sfp == NULL || rrp == NULL) return;
  bsp = BioseqFindFromSeqLoc (sfp->product);
  if (bsp == NULL) return;
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &context);
  if (sdp == NULL) return;
  mip = (MolInfoPtr) sdp->data.ptrvalue;
  if (mip == NULL) return;
  switch (rrp->type) {
    case 2 : /* mRNA */
      if (mip->biomol == MOLECULE_TYPE_MRNA) return;
      break;
    case 3 : /* tRNA */
      if (mip->biomol == MOLECULE_TYPE_TRNA) return;
      break;
    case 4 : /* rRNA */
      if (mip->biomol == MOLECULE_TYPE_RRNA) return;
      break;
    default :
      return;
  }
  ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_RnaProductMismatch, "Type of RNA does not match MolInfo of product Bioseq");
}

static BioseqSetPtr GetParentNPS (BioseqPtr bsp)
{
  BioseqSetPtr    bssp;

  if (bsp == NULL)
    return NULL;
  if (bsp->idx.parenttype != OBJ_BIOSEQSET)
    return NULL;
  bssp = (BioseqSetPtr) bsp->idx.parentptr;
  while (bssp != NULL && bssp->_class != BioseqseqSet_class_nuc_prot && bssp->idx.parenttype == OBJ_BIOSEQSET) {
    bssp = (BioseqSetPtr) bssp->idx.parentptr;
  }
  if (bssp != NULL && bssp->_class == BioseqseqSet_class_nuc_prot)
    return bssp;
  return NULL;
}

static Boolean NucAndProtNotInNPS (BioseqPtr nuc, BioseqPtr prot)
{
  BioseqSetPtr    bssp;

  if (nuc == NULL || prot == NULL)
    return FALSE;
  bssp = GetParentNPS (nuc);
  if (bssp == NULL)
    return TRUE;
  if (GetParentNPS (prot) != bssp)
    return TRUE;
  return FALSE;
}

static Boolean CDS5primePartialTest (
  SeqFeatPtr sfp
)

{
  BioseqPtr  nbsp;
  SeqLocPtr  slp = NULL;

  if (sfp == NULL) return FALSE;
  nbsp = BioseqFindFromSeqLoc (sfp->location);
  if (nbsp != NULL) {
    slp = SeqLocFindNext (sfp->location, NULL);
    if (slp != NULL) {
      if (SeqLocStrand (slp) == Seq_strand_minus) {
        if (SeqLocStop (slp) == nbsp->length - 1) {
          return TRUE;
        }
      } else {
        if (SeqLocStart (slp) == 0) {
          return TRUE;
        }
      }
    }
  }
  return FALSE;
}

static Boolean CDS3primePartialTest (
  SeqFeatPtr sfp
)

{
  BioseqPtr  nbsp;
  SeqLocPtr  last = NULL;
  SeqLocPtr  slp = NULL;

  if (sfp == NULL) return FALSE;
  nbsp = BioseqFindFromSeqLoc (sfp->location);
  if (nbsp != NULL) {
    last = NULL;
    slp = SeqLocFindNext (sfp->location, NULL);
    while (slp != NULL) {
      last = slp;
      slp = SeqLocFindNext (sfp->location, last);
    }
    if (last != NULL) {
      if (SeqLocStrand (last) == Seq_strand_minus) {
        if (SeqLocStart (last) == 0) {
          return TRUE;
        }
      } else {
        if (SeqLocStop (last) == nbsp->length - 1) {
          return TRUE;
        }
      }
    }
  }
  return FALSE;
}

static void CheckCDSPartial (ValidStructPtr vsp, SeqFeatPtr sfp)

{
  BioseqPtr          bsp;
  SeqMgrDescContext  context;
  MolInfoPtr         mip;
  Boolean            partial5;
  Boolean            partial3;
  SeqDescrPtr        sdp;
  ErrSev             sev;
  Boolean            need_to_unlock = FALSE;

  if (vsp == NULL || sfp == NULL) return;
  if (sfp->product == NULL) return;
  if (!vsp->useSeqMgrIndexes) return;
  bsp = BioseqFindFromSeqLoc (sfp->product);
  if (bsp == NULL && vsp->farFetchCDSproducts) {
    bsp = BioseqLockById (SeqLocId(sfp->product));
    if (bsp != NULL) {
      need_to_unlock = TRUE;
    }
  }
  if (bsp == NULL) return;
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &context);
  if (sdp == NULL) {
    if (need_to_unlock) {
      BioseqUnlock(bsp);
    }
    return;
  }
  mip = (MolInfoPtr) sdp->data.ptrvalue;
  if (mip == NULL) {
    if (need_to_unlock) {
      BioseqUnlock (bsp);
    }
    return;
  }
  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);

  switch (mip->completeness) {
    case 0 : /* unknown */
      break;
    case 1 : /* complete */
      if (partial5 || partial3) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_PartialProblem, "CDS is partial but protein is complete");
      }
      break;
    case 2 : /* partial */
      break;
    case 3 : /* no-left */
      if (! partial5) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_PartialProblem, "CDS is 5' complete but protein is NH2 partial");
      }
      if (partial3) {
        sev = SEV_ERROR;
        if (CDS3primePartialTest (sfp)) {
          sev = SEV_WARNING;
        }
        ValidErr (vsp, sev, ERR_SEQ_FEAT_PartialProblem, "CDS is 3' partial but protein is NH2 partial");
      }
      break;
    case 4 : /* no-right */
      if (! partial3) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_PartialProblem, "CDS is 3' complete but protein is CO2 partial");
      }
      if (partial5) {
        sev = SEV_ERROR;
        if (CDS5primePartialTest (sfp)) {
          sev = SEV_WARNING;
        }
        ValidErr (vsp, sev, ERR_SEQ_FEAT_PartialProblem, "CDS is 5' partial but protein is CO2 partial");
      }
      break;
    case 5 : /* no-ends */
      if (partial5 && partial3) {
      } else if (partial5) {
        sev = SEV_ERROR;
        if (CDS5primePartialTest (sfp)) {
          sev = SEV_WARNING;
        }
        ValidErr (vsp, sev, ERR_SEQ_FEAT_PartialProblem, "CDS is 5' partial but protein has neither end");
      } else if (partial3) {
        sev = SEV_ERROR;
        if (CDS3primePartialTest (sfp)) {
          sev = SEV_WARNING;
        }
        ValidErr (vsp, sev, ERR_SEQ_FEAT_PartialProblem, "CDS is 3' partial but protein has neither end");
      } else {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_PartialProblem, "CDS is complete but protein has neither end");
      }
      break;
    case 6 : /* has-left */
      break;
    case 7 : /* has-right */
      break;
    default :
      break;
  }
  if (need_to_unlock) {
    BioseqUnlock (bsp);
  }
}

static void CheckForCommonCDSProduct (ValidStructPtr vsp, SeqFeatPtr sfp)
{
  BioseqPtr       bsp;
  BioseqSetPtr    bssp;
  SeqFeatPtr      cds;
  CdRegionPtr     crp;
  SeqFeatPtr      gene;
  GeneRefPtr      grp;
  Boolean         is_nc = FALSE;
  Boolean         is_nc_gps = FALSE;
  Boolean         is_nt = FALSE;
  Boolean         is_nw = FALSE;
  BioseqPtr       nuc;
  SeqEntryPtr     sep;
  SeqIdPtr        sip;
  TextSeqIdPtr    tsip;

  if (sfp == NULL || sfp->pseudo)
    return;
  if (!vsp->useSeqMgrIndexes)
    return;
  crp = (CdRegionPtr) sfp->data.value.ptrvalue;
  if (crp != NULL && crp->orf)
    return;
  
  grp = SeqMgrGetGeneXref (sfp);
  if (grp == NULL || (!SeqMgrGeneIsSuppressed (grp))) {
    gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
    if (gene != NULL) {
      if (gene->pseudo) return;
      grp = (GeneRefPtr) gene->data.value.ptrvalue;
      if (grp != NULL && grp->pseudo) return;
    }
  }
  if (sfp->product == NULL) return;
  bsp = BioseqFindFromSeqLoc (sfp->product);
  if (bsp == NULL) {
    sip = SeqLocId (sfp->product);
    /* okay to have far RefSeq product... */
    if (sip == NULL || sip->choice != SEQID_OTHER) {
      sep = vsp->sep;
      if (sep != NULL && IS_Bioseq_set (sep)) {
        bssp = (BioseqSetPtr) sep->data.ptrvalue;
        /* but only if genomic product set */
        if (bssp != NULL && bssp->_class == BioseqseqSet_class_gen_prod_set)
          return;
        if (bssp != NULL && bssp->_class == BioseqseqSet_class_genbank) {
          sep = bssp->seq_set;
          if (sep != NULL && IS_Bioseq_set (sep)) {
            bssp = (BioseqSetPtr) sep->data.ptrvalue;
            if (bssp != NULL && bssp->_class == BioseqseqSet_class_gen_prod_set)
              return;
          }
        }
      }
      /* or just a bioseq */
      if (sep != NULL && IS_Bioseq (sep))
        return;
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_MissingCDSproduct, "Unable to find product Bioseq from CDS feature");
    }
    return;
  }
  nuc = BioseqFindFromSeqLoc (sfp->location);
  if (nuc != NULL) {
    for (sip = nuc->id; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_OTHER) {
        tsip = (TextSeqIdPtr) sip->data.ptrvalue;
        if (tsip != NULL && tsip->accession != NULL) {
          if (StringNICmp (tsip->accession, "NT_", 3) == 0) {
            is_nt = TRUE;
          } else if (StringNICmp (tsip->accession, "NC_", 3) == 0) {
            is_nc = TRUE;
          } else if (StringNICmp (tsip->accession, "NW_", 3) == 0) {
            is_nw = TRUE;
          }
        }
      }
    }
    if (/* (is_nc || is_nw) && */ nuc->idx.parenttype == OBJ_BIOSEQSET) {
      bssp = (BioseqSetPtr) nuc->idx.parentptr;
      if (bssp != NULL) {
        if (bssp->_class == BioseqseqSet_class_gen_prod_set) {
          is_nc_gps = TRUE;
        }
      }
    }
    if (NucAndProtNotInNPS (nuc, bsp) && (! is_nt) && (! is_nc_gps)) {
      if (vsp->is_small_genome_set) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_CDSproductPackagingProblem, "Protein product not packaged in nuc-prot set with nucleotide in small genome set");
      } else {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_CDSproductPackagingProblem, "Protein product not packaged in nuc-prot set with nucleotide");
      }
    }
  }
  cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
  if (cds == NULL) return;
  if (cds != sfp) {
    /* if genomic product set, with one cds on contig and one on cdna, do not report */
    sep = vsp->sep;
    if (sep != NULL && IS_Bioseq_set (sep)) {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      if (bssp != NULL && bssp->_class == BioseqseqSet_class_gen_prod_set) {
        /* feature packaging test will do final contig vs. cdna check */
        if (BioseqFindFromSeqLoc (cds->location) != BioseqFindFromSeqLoc (sfp->location))
          return;
      }
      if (bssp != NULL && bssp->_class == BioseqseqSet_class_genbank) {
        sep = bssp->seq_set;
        if (sep != NULL && IS_Bioseq_set (sep)) {
          bssp = (BioseqSetPtr) sep->data.ptrvalue;
          if (bssp != NULL && bssp->_class == BioseqseqSet_class_gen_prod_set)
            if (BioseqFindFromSeqLoc (cds->location) != BioseqFindFromSeqLoc (sfp->location))
              return;
        }
      }
    }

    ValidErr (vsp, SEV_REJECT, ERR_SEQ_FEAT_MultipleCDSproducts, "Same product Bioseq from multiple CDS features");
  }
}

static void CheckForCommonMRNAProduct (ValidStructPtr vsp, SeqFeatPtr sfp)
{
  BioseqPtr       bsp;
  BioseqSetPtr    bssp;
  SeqFeatPtr      gene;
  GeneRefPtr      grp;
  SeqFeatPtr      mrna;
  SeqEntryPtr     oldscope;
  SeqEntryPtr     sep;
  SeqIdPtr        sip;

  if (sfp == NULL || sfp->pseudo)
    return;
  if (!vsp->useSeqMgrIndexes)
    return;
  grp = SeqMgrGetGeneXref (sfp);
  if (grp == NULL || (!SeqMgrGeneIsSuppressed (grp))) {
    gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
    if (gene == NULL || gene->pseudo)
      return;
    grp = (GeneRefPtr) gene->data.value.ptrvalue;
    if (grp != NULL && grp->pseudo)
      return;
  }
  if (sfp->product == NULL)
    return;

  oldscope = SeqEntrySetScope (vsp->sep);
  bsp = BioseqFindFromSeqLoc (sfp->product);
  SeqEntrySetScope (oldscope);
  if (bsp == NULL) {
    sip = SeqLocId (sfp->product);
    if (sip != NULL && sip->choice == SEQID_LOCAL) {
      sep = vsp->sep;
      if (sep != NULL && IS_Bioseq_set (sep)) {
        bssp = (BioseqSetPtr) sep->data.ptrvalue;
        if (bssp != NULL) {
          if (bssp->_class == BioseqseqSet_class_gen_prod_set ||
              bssp->_class == BioseqseqSet_class_other) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_MissingMRNAproduct,
            "Product Bioseq of mRNA feature is not packaged in the record");
          }
        }
      }
    }
    return;
  }

  mrna = SeqMgrGetRNAgivenProduct (bsp, NULL);
  if (mrna == NULL)
    return;
  if (mrna != sfp) {
    ValidErr (vsp, SEV_REJECT, ERR_SEQ_FEAT_MultipleMRNAproducts, "Same product Bioseq from multiple mRNA features");
  }
}

static void CheckForBadGeneOverlap (ValidStructPtr vsp, SeqFeatPtr sfp)
{
  SeqMgrFeatContext fcontext;
  SeqFeatPtr      gene, operon;
  GeneRefPtr      grp;
  ErrSev          sev = /* SEV_ERROR */ SEV_WARNING;

  if (sfp == NULL)
    return;
  grp = SeqMgrGetGeneXref (sfp);
  if (grp != NULL)
    return;
  gene = SeqMgrGetOverlappingGene (sfp->location, &fcontext);
  if (gene != NULL)
    return;
  gene = SeqMgrGetOverlappingFeature (sfp->location, FEATDEF_GENE, NULL, 0, NULL, SIMPLE_OVERLAP, &fcontext);
  if (gene == NULL)
    return;
  if (IsNCorNT (vsp->sep, sfp->location)) {
    sev = SEV_WARNING;
  }
  if (sfp->data.choice == SEQFEAT_CDREGION) {
    ValidErr (vsp, sev, ERR_SEQ_FEAT_CDSgeneRange, "gene overlaps CDS but does not completely contain it");
  } else if (sfp->data.choice == SEQFEAT_RNA) {
    operon = SeqMgrGetOverlappingOperon (sfp->location, &fcontext);
    if (operon != NULL)
      return;
    ValidErr (vsp, sev, ERR_SEQ_FEAT_mRNAgeneRange, "gene overlaps mRNA but does not completely contain it");
  }
}

static void CheckForBadMRNAOverlap (ValidStructPtr vsp, SeqFeatPtr sfp)

{
  BioseqPtr          bsp;
  SeqMgrFeatContext  fcontext;
  SeqFeatPtr         gene = NULL;
  GeneRefPtr         grp;
  SeqFeatPtr         mrna;
  Boolean            pseudo = FALSE;
  ErrSev             sev = /* SEV_ERROR */ SEV_WARNING;

  if (sfp == NULL)
    return;

  if (sfp->pseudo) {
    pseudo = TRUE;
  }
  grp = SeqMgrGetGeneXref (sfp);
  if (grp != NULL) {
    if (SeqMgrGeneIsSuppressed (grp)) {
    } else {
      if (grp->pseudo) return;
      bsp = BioseqFindFromSeqLoc (sfp->location);
      if (bsp != NULL) {
        if (StringDoesHaveText (grp->locus_tag)) {
          gene = SeqMgrGetGeneByLocusTag (bsp, grp->locus_tag, &fcontext);
        } else if (StringDoesHaveText (grp->locus)) {
          gene = SeqMgrGetFeatureByLabel (bsp, grp->locus_tag, SEQFEAT_GENE, 0, &fcontext);
        }
        if (gene != NULL) {
          grp = (GeneRefPtr) gene->data.value.ptrvalue;
          if (grp != NULL && grp->pseudo) {
            pseudo = TRUE;
          }
        }
      }
    }
  }

  mrna = SeqMgrGetOverlappingFeature (sfp->location, FEATDEF_mRNA, NULL, 0, NULL, SIMPLE_OVERLAP, &fcontext);
  if (mrna == NULL)
    return;
  mrna = SeqMgrGetOverlappingFeature (sfp->location, FEATDEF_mRNA, NULL, 0, NULL, CHECK_INTERVALS, &fcontext);
  if (mrna != NULL)
    return;
  mrna = SeqMgrGetOverlappingFeature (sfp->location, FEATDEF_mRNA, NULL, 0, NULL, INTERVAL_OVERLAP, &fcontext);
  if (mrna == NULL)
    return;
  if (IsNCorNTorNW (vsp->sep, sfp->location)) {
    sev = SEV_WARNING;
  }
  if (sfp->excpt) {
    sev = SEV_WARNING;
  }
  mrna = SeqMgrGetOverlappingFeature (sfp->location, FEATDEF_mRNA, NULL, 0, NULL, LOCATION_SUBSET, &fcontext);
  if (mrna != NULL) {
    if (StringISearch (sfp->except_text, "ribosomal slippage") == NULL && StringISearch (sfp->except_text, "trans-splicing") == NULL) {
      if (pseudo) {
        ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_PseudoCDSmRNArange, "mRNA contains CDS but internal intron-exon boundaries do not match");
      } else {
        ValidErr (vsp, sev, ERR_SEQ_FEAT_CDSmRNArange, "mRNA contains CDS but internal intron-exon boundaries do not match");
      }
    }
  } else {
    if (pseudo) {
      ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_PseudoCDSmRNArange, "mRNA overlaps or contains CDS but does not completely contain intervals");
    } else {
      ValidErr (vsp, sev, ERR_SEQ_FEAT_CDSmRNArange, "mRNA overlaps or contains CDS but does not completely contain intervals");
    }
  }
}

/*
static void CheckForBothStrands (ValidStructPtr vsp, SeqFeatPtr sfp)
{
  Boolean    bothstrands = FALSE, bothreverse = FALSE;
  SeqLocPtr  location, slp = NULL;
  Uint1      strand;

  if (sfp == NULL)
    return;
  location = sfp->location;
  if (location == NULL)
    return;
  while ((slp = SeqLocFindNext (location, slp)) != NULL) {
    strand = SeqLocStrand (slp);
    if (strand == Seq_strand_both) {
      bothstrands = TRUE;
    } else if (strand == Seq_strand_both_rev) {
      bothreverse = TRUE;
    }
  }
  if (bothstrands && bothreverse) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_BothStrands, "mRNA or CDS may not be on both (forward and reverse) strands");
  } else if (bothstrands) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_BothStrands, "mRNA or CDS may not be on both (forward) strands");
  } else if (bothreverse) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_BothStrands, "mRNA or CDS may not be on both (reverse) strands");
  }
}
*/

static void CheckForBothOrBothRev (ValidStructPtr vsp, SeqFeatPtr sfp)

{
  Boolean    bothstrands = FALSE, bothreverse = FALSE, iswhole = FALSE;
  SeqLocPtr  location, slp = NULL;
  CharPtr    prefix = "Feature";
  ErrSev     sev = SEV_WARNING;
  CharPtr    suffix = "";
  Uint1      strand;

  if (sfp == NULL) return;
  location = sfp->location;
  if (location == NULL) return;

  if (sfp->idx.subtype == FEATDEF_CDS) {
    sev = SEV_ERROR;
    prefix = "CDS";
  } else if (sfp->idx.subtype == FEATDEF_mRNA) {
    sev = SEV_ERROR;
    prefix = "mRNA";
  }

  while ((slp = SeqLocFindNext (location, slp)) != NULL) {
    if (slp->choice == SEQLOC_WHOLE) {
      iswhole = TRUE;
    } else {
      strand = SeqLocStrand (slp);
      if (strand == Seq_strand_both) {
        bothstrands = TRUE;
      } else if (strand == Seq_strand_both_rev) {
        bothreverse = TRUE;
      }
    }
  }
  if (bothstrands && bothreverse) {
    suffix = "(forward and reverse)";
  } else if (bothstrands) {
    suffix = "(forward)";
  } else if (bothreverse) {
    suffix = "(reverse)";
  }
  if (bothstrands || bothreverse) {
    ValidErr (vsp, sev, ERR_SEQ_FEAT_BothStrands, "%s may not be on both %s strands", prefix, suffix);
  }
  if (iswhole) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_WholeLocation, "%s may not have whole location", prefix);
  }
}

static Boolean OverlappingGeneIsPseudo (SeqFeatPtr sfp)
{
  SeqFeatPtr      gene;
  GeneRefPtr      grp;

  if (sfp == NULL)
    return FALSE;
  grp = SeqMgrGetGeneXref (sfp);
  if (grp != NULL) {
    if (grp->pseudo)
      return TRUE;
    return FALSE;
  }
  gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
  if (gene != NULL) {
    if (gene->pseudo)
      return TRUE;
    grp = (GeneRefPtr) gene->data.value.ptrvalue;
    if (grp != NULL) {
      if (grp->pseudo)
        return TRUE;
    }
  }
  return FALSE;
}

static void CheckForIllegalDbxref (ValidStructPtr vsp, GatherContextPtr gcp, ValNodePtr dbxref)

{
  DbtagPtr     db;
  CharPtr      good;
  Int4         id;
  Boolean      is_bc;
  Boolean      is_rf;
  Boolean      is_sc;
  ObjectIdPtr  oip;
  ValNodePtr   vnp;

  for (vnp = dbxref; vnp != NULL; vnp = vnp->next) {
    id = -1;
    db = (DbtagPtr) vnp->data.ptrvalue;
    if (db != NULL && db->db != NULL) {

      if (DbxrefIsValid (db->db, &is_rf, &is_sc, &is_bc, &good)) {
        if (is_bc) {
          if (StringHasNoText (good)) {
            good = "?";
          }
          if (is_sc && StringICmp (db->db, "taxon") == 0) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_IllegalDbXref,
                      "Illegal db_xref type %s, legal capitalization is %s, but should only be used on an OrgRef",
                      db->db, good);
          } else {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_IllegalDbXref,
                      "Illegal db_xref type %s, legal capitalization is %s",
                      db->db, good);
          }
        } else if (is_rf) {
          if (vsp->is_refseq_in_sep || vsp->is_gps_in_sep) {
          } else {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_IllegalDbXref,
                      "db_xref type %s is only legal for RefSeq", db->db);
          }
        } else if (is_sc && StringICmp (db->db, "taxon") == 0) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_IllegalDbXref,
                    "db_xref type %s should only be used on an OrgRef", db->db);
        } else {
        }
      } else {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_IllegalDbXref, "Illegal db_xref type %s", db->db);
      }

      if (StringDoesHaveText (db->db)) {
        if (StringHasSgml (vsp, db->db)) {
          ValidErr (vsp, SEV_WARNING, ERR_GENERIC_SgmlPresentInText, "dbxref database %s has SGML", db->db);
        }
      }

      oip = db->tag;
      if (oip != NULL && StringDoesHaveText (oip->str)) {
        if (StringHasSgml (vsp, oip->str)) {
          ValidErr (vsp, SEV_WARNING, ERR_GENERIC_SgmlPresentInText, "dbxref value %s has SGML", oip->str);
        }
      }

      /*
      dbxerr = NULL;
      valid = IsDbxrefValid (db->db, sfp, NULL,
                             GPSorRefSeq (vsp->sep, sfp->location),
                             &dbxerr);
      if (dbxerr != NULL) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_IllegalDbXref, dbxerr);
        dbxerr = MemFree (dbxerr);
      }
      */
    }
  }
}

static CharPtr plastidtxt [] = {
  "",
  "",
  "chloroplast",
  "chromoplast",
  "",
  "",
  "plastid",
  "",
  "",
  "",
  "",
  "",
  "cyanelle",
  "",
  "",
  "",
  "apicoplast",
  "leucoplast",
  "proplastid",
  "",
  ""
};

static CharPtr legal_exception_strings [] = {
  "RNA editing",
  "reasons given in citation",
  "rearrangement required for product",
  "ribosomal slippage",
  "trans-splicing",
  "alternative processing",
  "artificial frameshift",
  "nonconsensus splice site",
  "modified codon recognition",
  "alternative start codon",
  "dicistronic gene",
  "transcribed product replaced",
  "translated product replaced",
  "transcribed pseudogene",
  "annotated by transcript or proteomic data",
  "heterogeneous population sequenced",
  "low-quality sequence region",
  "unextendable partial coding region",
  "artificial location",
  NULL
};

static CharPtr refseq_exception_strings [] = {
  "unclassified transcription discrepancy",
  "unclassified translation discrepancy",
  "mismatches in transcription",
  "mismatches in translation",
  "adjusted for low-quality genome",
  NULL
};

static void ValidateExceptText (ValidStructPtr vsp, GatherContextPtr gcp, SeqFeatPtr sfp)

{
  Boolean    art_loc_except = FALSE;
  Boolean    found;
  GBQualPtr  gbq;
  Int2       i;
  Boolean    has_inference = FALSE;
  CharPtr    ptr;
  Boolean    reasons_given_except = FALSE;
  Boolean    redundant_with_comment = FALSE;
  Boolean    refseq_except = FALSE;
  ErrSev     sev = SEV_ERROR;
  Boolean    trans_prot_except = FALSE;
  CharPtr    str;
  CharPtr    tmp;

  str = StringSave (sfp->except_text);
  if (str == NULL) return;
  tmp = str;
  while (! StringHasNoText (tmp)) {
    ptr = StringChr (tmp, ',');
    if (ptr != NULL) {
      *ptr = '\0';
      ptr++;
    }
    TrimSpacesAroundString (tmp);
    found = FALSE;
    for (i = 0; legal_exception_strings[i] != NULL; i++) {
      if (StringICmp (tmp, legal_exception_strings[i]) == 0) {
        found = TRUE;
        if (StringICmp (tmp, "reasons given in citation") == 0) {
          reasons_given_except = TRUE;
        } else if (StringICmp (tmp, "annotated by transcript or proteomic data") == 0) {
          trans_prot_except = TRUE;
        } else if (StringICmp (tmp, "artificial location") == 0) {
          art_loc_except = TRUE;
        }
        break;
      }
    }
    if (!found) {
      if (GPSorRefSeq (vsp->sep, sfp->location)) {
        for (i = 0; refseq_exception_strings[i] != NULL; i++) {
          if (StringICmp (tmp, refseq_exception_strings[i]) == 0) {
            found = TRUE;
            refseq_except = TRUE;
            break;
          }
        }
      }
      if (! found) {
        if (IsNCorNT (vsp->sep, sfp->location)) {
          sev = SEV_WARNING;
        }
        ValidErr (vsp, sev, ERR_SEQ_FEAT_ExceptionProblem, "%s is not a legal exception explanation", tmp);
      }
    }
    if (sfp->comment != NULL && StringISearch (sfp->comment, tmp) != NULL) {
      if (StringICmp (tmp, "ribosomal slippage") != 0 &&
          StringICmp (tmp, "trans-splicing") != 0 &&
          StringICmp (tmp, "artificial location") != 0) {
        redundant_with_comment = TRUE;
      } else if (StringICmp (sfp->comment, tmp) == 0) {
        redundant_with_comment = TRUE;
      }
    }
    tmp = ptr;
  }
  MemFree (str);
  if (redundant_with_comment) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_ExceptionProblem, "Exception explanation text is also found in feature comment");
  }
  if (refseq_except) {
    found = FALSE;
    for (i = 0; refseq_exception_strings[i] != NULL; i++) {
      if (StringICmp (sfp->except_text, refseq_exception_strings[i]) == 0) {
        found = TRUE;
        refseq_except = TRUE;
        break;
      }
    }
    if (! found) {
      if (! vsp->is_gpipe_in_sep) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_ExceptionProblem, "Genome processing exception should not be combined with other explanations");
      }
    }
  }
  if (reasons_given_except && sfp->cit == NULL) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_ExceptionProblem, "Reasons given in citation exception does not have the required citation");
  }
  if (trans_prot_except) {
    for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
      if (StringICmp (gbq->qual, "inference") == 0) {
        has_inference = TRUE;
      }
    }
    if (! has_inference) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_ExceptionProblem,
                "Annotated by transcript or proteomic data exception does not have the required inference qualifier");
    }
  }
  if (art_loc_except) {
    if (! vsp->is_embl_ddbj_in_sep) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_ExceptionProblem,
                "Artificial location should only be used directly by EMBL or DDBJ records");
    }
  }
}

typedef struct samecds {
  Boolean               found;
  SeqMgrFeatContextPtr  gcontext;
  Uint2                 slpTag;
  Uint1                 subtype;
  Boolean               bypassGeneTest;
} SameCds, PNTR SameCdsPtr;

static Boolean LIBCALLBACK FindSameCDS (SeqFeatPtr sfp, SeqMgrFeatContextPtr ccontext)

{
  SeqMgrFeatContextPtr  gcontext;
  Int2                  i;
  SameCdsPtr            same;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return TRUE;
  same = (SameCdsPtr) ccontext->userdata;
  gcontext = same->gcontext;
  if (gcontext == NULL || gcontext->sfp == NULL ||
      gcontext->ivals == NULL || ccontext->ivals == NULL) return TRUE;
  if (gcontext->strand == ccontext->strand ||
      (ccontext->strand == Seq_strand_unknown && gcontext->strand != Seq_strand_minus) ||
      (gcontext->strand == Seq_strand_unknown && ccontext->strand != Seq_strand_minus) ||
      ccontext->strand == Seq_strand_both) {
    /* test for strands from SeqMgrGetBestOverlappingFeat, keep going if okay */
  } else {
    return TRUE;
  }
  if (same->subtype == FEATDEF_GENE) {
    if (gcontext->left == ccontext->left &&
        gcontext->right == ccontext->right) {
      same->found = TRUE;
      return FALSE;
    }
  } else if (same->subtype == FEATDEF_mRNA) {
    if (gcontext->left == ccontext->left &&
        gcontext->right == ccontext->right &&
        gcontext->numivals == ccontext->numivals) {
      if (SeqLocAinB (sfp->location, gcontext->sfp->location) >= 0) {
        if (gcontext->numivals == 1) {
          same->found = TRUE;
          return FALSE;
        } else {
          for (i = 0; i < gcontext->numivals; i++) {
            if (gcontext->ivals [2 * i] != ccontext->ivals [2 * i]) return TRUE;
            if (gcontext->ivals [2 * i + 1] != ccontext->ivals [2 * i + 1]) return TRUE;
          }
          same->found = TRUE;
          return FALSE;
        }
      }
    } else if (SeqLocAinB (sfp->location, gcontext->sfp->location) > 0) {

      if (ccontext->strand == Seq_strand_minus || gcontext->strand == Seq_strand_minus) {
        if (same->slpTag == SLP_NOSTART && gcontext->partialL) {
          if (gcontext->right == ccontext->right) {
            same->found = TRUE;
            return FALSE;
          }
          if (gcontext->right > ccontext->right) {
            same->bypassGeneTest = TRUE;
            return FALSE;
          }
        } else if (same->slpTag == SLP_NOSTOP && gcontext->partialR) {
          if (gcontext->left == ccontext->left) {
            same->found = TRUE;
            return FALSE;
          }
          if (gcontext->left < ccontext->left) {
            same->bypassGeneTest = TRUE;
            return FALSE;
          }
        }

      } else {

        if (same->slpTag == SLP_NOSTART && gcontext->partialL) {
          if (gcontext->left == ccontext->left) {
            same->found = TRUE;
            return FALSE;
          }
          if (gcontext->left < ccontext->left) {
            same->bypassGeneTest = TRUE;
            return FALSE;
          }
        } else if (same->slpTag == SLP_NOSTOP && gcontext->partialR) {
          if (gcontext->right == ccontext->right) {
            same->found = TRUE;
            return FALSE;
          }
          if (gcontext->right > ccontext->right) {
            same->bypassGeneTest = TRUE;
            return FALSE;
          }
        }
      }
    }
  }
  return TRUE;
}

static Boolean SameAsCDS (SeqFeatPtr sfp, Uint2 slpTag, BoolPtr bypassGeneTestP)

{
  BioseqPtr          bsp;
  SeqMgrFeatContext  ccontext;
  SeqFeatPtr         cds;
  Boolean            cdsFilt [SEQFEAT_MAX];
  SeqMgrFeatContext  gcontext;
  SameCds            same;
  VvmDataPtr         vdp;
  SeqFeatXrefPtr     xref;

  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) return FALSE;
  if (SeqMgrGetDesiredFeature (0, bsp, 0, 0, sfp, &gcontext) != sfp) return FALSE;
  same.found = FALSE;
  same.gcontext = &gcontext;
  same.slpTag = slpTag;
  same.subtype = sfp->idx.subtype;
  same.bypassGeneTest = FALSE;

  vdp = (VvmDataPtr) sfp->idx.scratch;
  if (vdp != NULL && vdp->nearbycds != NULL) {
    cds = SeqMgrGetDesiredFeature (0, bsp, 0, 0, vdp->nearbycds, &ccontext);
    if (cds != NULL && cds->idx.subtype == FEATDEF_CDS && cds == vdp->nearbycds) {
      ccontext.userdata = (Pointer) &same;
      FindSameCDS (cds, &ccontext);
      if (same.found) {
        if (bypassGeneTestP != NULL) {
          *bypassGeneTestP = same.bypassGeneTest;
        }
        return same.found;
      }
      same.bypassGeneTest = FALSE;
    }
  }

  for (xref = sfp->xref; xref != NULL; xref = xref->next) {
    if (xref->id.choice == 0) continue;
    cds = SeqMgrGetFeatureByFeatID (sfp->idx.entityID, NULL, NULL, xref, &ccontext);
    if (cds == NULL || cds->idx.subtype != FEATDEF_CDS) continue;
    ccontext.userdata = (Pointer) &same;
    FindSameCDS (cds, &ccontext);
    if (same.found) {
      if (bypassGeneTestP != NULL) {
        *bypassGeneTestP = same.bypassGeneTest;
      }
      return same.found;
    }
    same.bypassGeneTest = FALSE;
  }

  MemSet ((Pointer) &cdsFilt, 0, sizeof (cdsFilt));
  cdsFilt [SEQFEAT_CDREGION] = TRUE;
  SeqMgrExploreFeatures (bsp, (Pointer) &same, FindSameCDS, sfp->location, cdsFilt, NULL);
  if (bypassGeneTestP != NULL) {
    *bypassGeneTestP = same.bypassGeneTest;
  }
  return same.found;
}

static Boolean LIBCALLBACK FindSameMRNA (SeqFeatPtr sfp, SeqMgrFeatContextPtr ccontext)

{
  SeqMgrFeatContextPtr  gcontext;
  Int2                  i;
  SameCdsPtr            same;

  if (sfp == NULL || sfp->idx.subtype != FEATDEF_mRNA) return TRUE;
  same = (SameCdsPtr) ccontext->userdata;
  gcontext = same->gcontext;
  if (gcontext == NULL || gcontext->sfp == NULL ||
      gcontext->ivals == NULL || ccontext->ivals == NULL) return TRUE;
  if (gcontext->strand == ccontext->strand ||
      (ccontext->strand == Seq_strand_unknown && gcontext->strand != Seq_strand_minus) ||
      (gcontext->strand == Seq_strand_unknown && ccontext->strand != Seq_strand_minus) ||
      ccontext->strand == Seq_strand_both) {
    /* test for strands from SeqMgrGetBestOverlappingFeat, keep going if okay */
  } else {
    return TRUE;
  }
  if (same->subtype == FEATDEF_GENE) {
    if (gcontext->left == ccontext->left &&
        gcontext->right == ccontext->right) {
      same->found = TRUE;
      return FALSE;
    }
  } else if (same->subtype == FEATDEF_CDS) {
    if (gcontext->left == ccontext->left &&
        gcontext->right == ccontext->right &&
        gcontext->numivals == ccontext->numivals) {
      if (SeqLocAinB (gcontext->sfp->location, sfp->location) >= 0) {
        if (gcontext->numivals == 1) {
          same->found = TRUE;
          return FALSE;
        } else {
          for (i = 0; i < gcontext->numivals; i++) {
            if (gcontext->ivals [2 * i] != ccontext->ivals [2 * i]) return TRUE;
            if (gcontext->ivals [2 * i + 1] != ccontext->ivals [2 * i + 1]) return TRUE;
          }
          same->found = TRUE;
          return FALSE;
        }
      }
    }
  } else if (same->subtype == FEATDEF_exon) {
    if (ccontext->strand == Seq_strand_minus || gcontext->strand == Seq_strand_minus) {
      if (same->slpTag == SLP_NOSTART && ccontext->partialL) {
        if (gcontext->right == ccontext->right) {
          same->found = TRUE;
          return FALSE;
        }
      } else if (same->slpTag == SLP_NOSTOP && ccontext->partialR) {
        if (gcontext->left == ccontext->left) {
          same->found = TRUE;
          return FALSE;
        }
      }

    } else {

      if (same->slpTag == SLP_NOSTART && ccontext->partialL) {
        if (gcontext->left == ccontext->left) {
          same->found = TRUE;
          return FALSE;
        }
      } else if (same->slpTag == SLP_NOSTOP && ccontext->partialR) {
        if (gcontext->right == ccontext->right) {
          same->found = TRUE;
          return FALSE;
        }
      }
    }
  }
  return TRUE;
}

static Boolean SameAsMRNA (SeqFeatPtr sfp, Uint2 slpTag)

{
  BioseqPtr          bsp;
  SeqMgrFeatContext  mcontext;
  SeqFeatPtr         mrna;
  Boolean            mrnaFilt [FEATDEF_MAX];
  SeqMgrFeatContext  gcontext;
  SameCds            same;
  VvmDataPtr         vdp;
  SeqFeatXrefPtr     xref;

  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) return FALSE;
  if (SeqMgrGetDesiredFeature (0, bsp, 0, 0, sfp, &gcontext) != sfp) return FALSE;
  same.found = FALSE;
  same.gcontext = &gcontext;
  same.slpTag = slpTag;
  same.subtype = sfp->idx.subtype;

  vdp = (VvmDataPtr) sfp->idx.scratch;
  if (vdp != NULL && vdp->nearbymrna != NULL) {
    mrna = SeqMgrGetDesiredFeature (0, bsp, 0, 0, vdp->nearbymrna, &mcontext);
    if (mrna != NULL && mrna->idx.subtype == FEATDEF_mRNA && mrna == vdp->nearbymrna) {
      mcontext.userdata = (Pointer) &same;
      FindSameMRNA (mrna, &mcontext);
      if (same.found) {
        return same.found;
      }
    }
  }

  for (xref = sfp->xref; xref != NULL; xref = xref->next) {
    if (xref->id.choice == 0) continue;
    mrna = SeqMgrGetFeatureByFeatID (sfp->idx.entityID, NULL, NULL, xref, &mcontext);
    if (mrna == NULL || mrna->idx.subtype != FEATDEF_mRNA) continue;
    mcontext.userdata = (Pointer) &same;
    FindSameMRNA (mrna, &mcontext);
    if (same.found) {
      return same.found;
    }
  }

  MemSet ((Pointer) &mrnaFilt, 0, sizeof (mrnaFilt));
  mrnaFilt [FEATDEF_mRNA] = TRUE;
  SeqMgrExploreFeatures (bsp, (Pointer) &same, FindSameMRNA, sfp->location, NULL, mrnaFilt);
  return same.found;
}

static Boolean TestSameGene (SeqMgrFeatContextPtr ccontext, SeqMgrFeatContextPtr gcontext)

{
  if (ccontext == NULL || ccontext->sfp == NULL ||
      gcontext == NULL || gcontext->sfp == NULL ||
      gcontext->ivals == NULL || ccontext->ivals == NULL) return FALSE;
  if (gcontext->strand == ccontext->strand ||
      (ccontext->strand == Seq_strand_unknown && gcontext->strand != Seq_strand_minus) ||
      (gcontext->strand == Seq_strand_unknown && ccontext->strand != Seq_strand_minus) ||
      ccontext->strand == Seq_strand_both) {
    /* test for strands from SeqMgrGetBestOverlappingFeat, keep going if okay */
  } else {
    return FALSE;
  }
  if (gcontext->left == ccontext->left &&
      gcontext->right == ccontext->right) {
    return TRUE;
  }
  return FALSE;
}

static Boolean SameAsGene (SeqFeatPtr sfp)

{
  BioseqPtr          bsp;
  SeqMgrFeatContext  ccontext;
  SeqFeatPtr         gene;
  SeqMgrFeatContext  gcontext;
  GeneRefPtr         grp;

  grp = SeqMgrGetGeneXref (sfp);
  if (grp != NULL) return FALSE;
  gene = SeqMgrGetOverlappingGene (sfp->location, &gcontext);
  if (gene == NULL) return FALSE;
  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) return FALSE;
  if (SeqMgrGetDesiredFeature (0, bsp, 0, 0, sfp, &ccontext) != sfp) return FALSE;
  return TestSameGene (&gcontext, &ccontext);
}

static Boolean SplicingNotExpected (SeqFeatPtr sfp)

{
  BioSourcePtr       biop;
  BioseqPtr          bsp;
  SeqMgrDescContext  dcontext;
  OrgNamePtr         onp;
  OrgRefPtr          orp;
  SeqDescrPtr        sdp;

  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) return FALSE;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp == NULL) return FALSE;
  biop = (BioSourcePtr) sdp->data.ptrvalue;
  if (biop == NULL) return FALSE;
  orp = biop->org;
  if (orp == NULL) return FALSE;
  onp = orp->orgname;
  if (onp == NULL) return FALSE;

  if (StringCmp (onp->div, "BCT") == 0) return TRUE;
  if (StringCmp (onp->div, "VRL") == 0) return TRUE;
  if (StringNICmp (onp->lineage, "Bacteria; ", 10) == 0) return TRUE;
  if (StringNICmp (onp->lineage, "Archaea; ", 9) == 0) return TRUE;

  return FALSE;
}

static Boolean RareConsensusNotExpected (SeqFeatPtr sfp)

{
  BioSourcePtr       biop;
  BioseqPtr          bsp;
  SeqMgrDescContext  dcontext;
  OrgNamePtr         onp;
  OrgRefPtr          orp;
  SeqDescrPtr        sdp;

  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) return FALSE;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp == NULL) return FALSE;
  biop = (BioSourcePtr) sdp->data.ptrvalue;
  if (biop == NULL) return FALSE;
  orp = biop->org;
  if (orp == NULL) return FALSE;
  onp = orp->orgname;
  if (onp == NULL) return FALSE;

  if (StringCmp (onp->div, "PLN") != 0) return TRUE;
  if (StringNICmp (onp->lineage, "Eukaryota; Viridiplantae; ", 26) != 0) return TRUE;

  return FALSE;
}

static Boolean HasUnderscore (CharPtr str)
{
  if (StringChr(str, '_') != NULL)
    return TRUE;
  else
    return FALSE;
}

static Boolean IsUpperCaseChar (Char ch)
{
  if (StringChr("ABCDEFGHIJKLMNOPQRSTUVWXYZ",ch) != NULL)
    return TRUE;
  else
    return FALSE;
}

/*
static Boolean IsNumericChar (Char ch)
{
  if (StringChr("0123456789",ch) != NULL)
    return TRUE;
  else
    return FALSE;
}
*/

NLM_EXTERN Boolean IsNuclAcc (CharPtr name)

{
  if (!IsUpperCaseChar (name[0]))
    return FALSE;

  if (!HasUnderscore (name))
    return FALSE;

  return TRUE;
}

static Boolean IsCddFeat (
  SeqFeatPtr sfp
)

{
  DbtagPtr    dbt;
  ValNodePtr  vnp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_REGION) return FALSE;

  for (vnp = sfp->dbxref; vnp != NULL; vnp = vnp->next) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt == NULL) continue;
    if (StringCmp (dbt->db, "CDD") == 0 || StringCmp (dbt->db, "cdd") == 0) return TRUE;
  }

  return FALSE;
}

/* derived from ValidateSeqLoc */
static void ValidateAnticodon (ValidStructPtr vsp, SeqLocPtr slp)
{
  SeqLocPtr       tmp;
  Boolean         retval = TRUE, tmpval, mixed_strand = FALSE, unmarked_strand = FALSE,
                  ordered = TRUE, adjacent = FALSE, circular = FALSE, exception = FALSE;
  CharPtr         ctmp;
  Uint1           strand1 = Seq_strand_other, strand2 = Seq_strand_other;
  Int4            from1 = -1, from2 = -1, to1 = -1, to2 = -1;
  SeqIntPtr       sip;
  SeqPntPtr       spp;
  SeqIdPtr        id1 = NULL, id2 = NULL;
  BioseqPtr       bsp;
  SeqFeatPtr      sfp = NULL;

  if (slp == NULL)
    return;

  sfp = vsp->sfp;

  bsp = BioseqFindFromSeqLoc (slp);
  if (bsp != NULL && bsp->topology == 2) {
    circular = TRUE;
  }

  tmp = NULL;

  for (tmp = SeqLocFindNext (slp, NULL); tmp != NULL; tmp = SeqLocFindNext (slp, tmp)) {
    tmpval = TRUE;
    switch (tmp->choice) {
    case SEQLOC_INT:
      sip = (SeqIntPtr) (tmp->data.ptrvalue);
      if (sip == NULL) continue;
      strand2 = sip->strand;
      id2 = sip->id;
      from2 = sip->from;
      to2 = sip->to;
      tmpval = SeqIntCheck (sip);
      break;
    case SEQLOC_PNT:
      spp = (SeqPntPtr) (tmp->data.ptrvalue);
      if (spp == NULL) continue;
      strand2 = spp->strand;
      id2 = spp->id;
      from2 = spp->point;
      to2 = spp->point;
      tmpval = SeqPntCheck (spp);
      break;
    case SEQLOC_NULL:
      continue;
    default:
      break;
    }

    if (id1 != NULL && id2 != NULL) {
      if (SeqIdForSameBioseq (id1, id2)) {
        if ((ordered) /* && (! circular) */) {
          if (strand2 == Seq_strand_minus) {
            if (to1 < to2)
              ordered = FALSE;
            if (to2 + 1 == from1)
              adjacent = TRUE;
          } else {
            if (to1 > to2)
              ordered = FALSE;
            if (to1 + 1 == from2)
              adjacent = TRUE;
          }
        }
        if (strand1 == strand2 && from1 == from2 && to1 == to2) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_DuplicateInterval, "Duplicate anticodon exons in location");
        }
      }
    }

    if (!tmpval) {
      retval = FALSE;
      ctmp = SeqLocPrint (tmp);
      if (ctmp != NULL && StringLen (ctmp) > 800) {
        StringCpy (ctmp + 797, "...");
      }
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_Range, "Anticodon location [%s] out of range", ctmp);
      MemFree (ctmp);
    }

    if ((strand1 != Seq_strand_other) && (strand2 != Seq_strand_other)) {
      if (SeqIdForSameBioseq (id1, id2)) {
        if (strand1 != strand2) {
          if (strand1 == Seq_strand_plus && strand2 == Seq_strand_unknown) {
            unmarked_strand = TRUE;
          } else if (strand1 == Seq_strand_unknown && strand2 == Seq_strand_plus) {
            unmarked_strand = TRUE;
          } else {
            mixed_strand = TRUE;
          }
        }
      }
    }

    from1 = from2;
    to1 = to2;
    id1 = id2;
    strand1 = strand2;
  }

  if (sfp != NULL && sfp->excpt) {
    exception = TRUE;
  }

  if (adjacent) {
    ctmp = SeqLocPrint (slp);
    /*
    if (exception) {
      sev = SEV_WARNING;
    } else {
      sev = SEV_ERROR;
    }
    */
    if (ctmp != NULL && StringLen (ctmp) > 800) {
      StringCpy (ctmp + 797, "...");
    }
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_AbuttingIntervals, "Adjacent intervals in Anticodon [%s]", ctmp);
    MemFree (ctmp);
  }

  if (sfp != NULL) {
    strand1 = SeqLocStrand (sfp->location);
    strand2 = SeqLocStrand (slp);
    if (strand1 == Seq_strand_minus && strand2 != Seq_strand_minus) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_BadAnticodonStrand, "Anticodon should be on minus strand");
    } else if (strand1 != Seq_strand_minus && strand2 == Seq_strand_minus) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_BadAnticodonStrand, "Anticodon should be on plus strand");
    }
  }

  if (exception) {
    /* trans splicing exception turns off both mixed_strand and out_of_order messages */
    if (StringISearch (sfp->except_text, "trans-splicing") != NULL) {
      return;
    }
  }

  if (mixed_strand || unmarked_strand || (!ordered)) {
    ctmp = SeqLocPrint (slp);
    if (ctmp != NULL && StringLen (ctmp) > 800) {
      StringCpy (ctmp + 797, "...");
    }
    if (mixed_strand) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_MixedStrand, "Mixed strands in Anticodon [%s]", ctmp);
    } else if (unmarked_strand) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_MixedStrand, "Mixed plus and unknown strands in Anticodon [%s]", ctmp);
    }
    if (!ordered)
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_SeqLocOrder, "Intervals out of order in Anticodon [%s]", ctmp);
    MemFree (ctmp);
    return;
  }

  /* newer check for intervals out of order on segmented bioseq */

  if (bsp == NULL || bsp->repr != Seq_repr_seg) return;

  if (SeqLocBadSortOrder (bsp, slp)) {
    ctmp = SeqLocPrint (slp);
    if (ctmp != NULL && StringLen (ctmp) > 800) {
      StringCpy (ctmp + 797, "...");
    }
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_SeqLocOrder, "Intervals out of order in Anticodon [%s]", ctmp);
    MemFree (ctmp);
  }

  /* newer check for mixed strand on segmented bioseq */

  if (SeqLocMixedStrands (bsp, slp)) {
    ctmp = SeqLocPrint (slp);
    if (ctmp != NULL && StringLen (ctmp) > 800) {
      StringCpy (ctmp + 797, "...");
    }
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_MixedStrand, "Mixed strands in Anticodon [%s]", ctmp);
    MemFree (ctmp);
  }
}

static Boolean JustQuotes (CharPtr str)

{
  Char  ch;

  if (str == NULL) return FALSE;

  ch = *str;
  while (ch != '\0') {
    if (ch != '"' && ch != ' ') return FALSE;
    str++;
    ch = *str;
  }

  return TRUE;
}

typedef struct dummysmfedata {
  Int4        max;
  Int4        num_at_max;
  Int4        num_trans_spliced;
  Boolean     equivalent_genes;
  GeneRefPtr  grp_at_max;
} DummySmfeData, PNTR DummySmfePtr;

static Boolean LIBCALLBACK DummySMFEProc (
  SeqFeatPtr sfp,
  SeqMgrFeatContextPtr context
)


{
  DummySmfePtr  dsp;
  GeneRefPtr    grp, grpx;
  Int4          len;
  Boolean       redundantgenexref = FALSE;
  CharPtr       syn1, syn2;

  if (sfp == NULL || context == NULL) return TRUE;
  dsp = context->userdata;
  if (dsp == NULL) return TRUE;
  if (sfp->data.choice != SEQFEAT_GENE) return TRUE;
  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (grp == NULL) return TRUE;

  len = SeqLocLen (sfp->location);
  if (len < dsp->max) {
    dsp->max = len;
    dsp->num_at_max = 1;
    dsp->num_trans_spliced = 0;
    if (StringISearch (sfp->except_text, "trans-splicing") != NULL) {
      (dsp->num_trans_spliced)++;
    }
    dsp->equivalent_genes = FALSE;
    dsp->grp_at_max = grp;
  } else if (len == dsp->max) {
    (dsp->num_at_max)++;
    if (StringISearch (sfp->except_text, "trans-splicing") != NULL) {
      (dsp->num_trans_spliced)++;
    }
    grpx = dsp->grp_at_max;
    if (grpx != NULL) {
      redundantgenexref = FALSE;
      if (StringDoesHaveText (grp->locus_tag) && StringDoesHaveText (grpx->locus_tag)) {
        if (StringICmp (grp->locus_tag, grpx->locus_tag) == 0) {
          redundantgenexref = TRUE;
        }
      } else if (StringDoesHaveText (grp->locus) && StringDoesHaveText (grpx->locus)) {
        if (StringICmp (grp->locus, grpx->locus) == 0) {
          redundantgenexref = TRUE;
        }
      } else if (grp->syn != NULL && grpx->syn != NULL) {
        syn1 = (CharPtr) grp->syn->data.ptrvalue;
        syn2 = (CharPtr) grpx->syn->data.ptrvalue;
        if (StringDoesHaveText (syn1) && StringDoesHaveText (syn2)) {
          if (StringICmp (syn1, syn2) == 0) {
            redundantgenexref = TRUE;
          }
        }
      }
    }
    if (redundantgenexref) {
      dsp->equivalent_genes = TRUE;
    }
  }

  return TRUE;
}

typedef struct govstruc {
  CharPtr  term;
  CharPtr  goid;
  CharPtr  evidence;
  Int4     pmid;
} GovStruc, PNTR GovStrucPtr;

static int LIBCALLBACK SortVnpByGvspTermFirst (VoidPtr ptr1, VoidPtr ptr2)

{
  int           compare;
  GovStrucPtr   gsp1, gsp2;
  ValNodePtr    vnp1, vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  gsp1 = (GovStrucPtr) vnp1->data.ptrvalue;
  gsp2 = (GovStrucPtr) vnp2->data.ptrvalue;
  if (gsp1 == NULL || gsp2 == NULL) return 0;

  compare = StringICmp (gsp1->term, gsp2->term);
  if (compare > 0) {
    return 1;
  } else if (compare < 0) {
    return -1;
  }

  compare = StringICmp (gsp1->goid, gsp2->goid);
  if (compare > 0) {
    return 1;
  } else if (compare < 0) {
    return -1;
  }

  compare = StringICmp (gsp1->evidence, gsp2->evidence);
  if (compare > 0) {
    return 1;
  } else if (compare < 0) {
    return -1;
  }

  if (gsp1->pmid == 0) return 1;
  if (gsp2->pmid == 0) return -1;
  if (gsp1->pmid > gsp2->pmid) {
    return 1;
  } else if (gsp1->pmid < gsp2->pmid) {
    return -1;
  }

  return 0;
}

static int LIBCALLBACK SortVnpByGvspGoIDFirst (VoidPtr ptr1, VoidPtr ptr2)

{
  int           compare;
  GovStrucPtr   gsp1, gsp2;
  ValNodePtr    vnp1, vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  gsp1 = (GovStrucPtr) vnp1->data.ptrvalue;
  gsp2 = (GovStrucPtr) vnp2->data.ptrvalue;
  if (gsp1 == NULL || gsp2 == NULL) return 0;

  compare = StringICmp (gsp1->goid, gsp2->goid);
  if (compare > 0) {
    return 1;
  } else if (compare < 0) {
    return -1;
  }

  compare = StringICmp (gsp1->term, gsp2->term);
  if (compare > 0) {
    return 1;
  } else if (compare < 0) {
    return -1;
  }

  compare = StringICmp (gsp1->evidence, gsp2->evidence);
  if (compare > 0) {
    return 1;
  } else if (compare < 0) {
    return -1;
  }

  if (gsp1->pmid == 0) return 1;
  if (gsp2->pmid == 0) return -1;
  if (gsp1->pmid > gsp2->pmid) {
    return 1;
  } else if (gsp1->pmid < gsp2->pmid) {
    return -1;
  }

  return 0;
}


static ValNodePtr ValidateGoTermQualifier (
  ValidStructPtr vsp,
  UserFieldPtr   field_list
)

{
  UserFieldPtr term, ufp;
  CharPtr      textstr, evidence, goid;
  Char         gid[255];
  Int4         pmid, j;
  ObjectIdPtr  oip;
  ValNodePtr   head = NULL, vnp;
  GovStrucPtr  gsp, lastgsp;

  for (term = field_list; term != NULL; term = term->next) {
    if (term->choice != 11 || term->data.ptrvalue == NULL) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadGeneOntologyFormat, "Bad GO term format");
    } else {
      textstr = NULL;
      evidence = NULL;
      goid = NULL;
      pmid = 0;
      for (ufp = (UserFieldPtr) term->data.ptrvalue; ufp != NULL; ufp = ufp->next) {
        oip = ufp->label;
        if (oip == NULL) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadGeneOntologyFormat, "No label on GO term qualifier field");
          continue;
        }
        for (j = 0; goFieldType [j] != NULL; j++) {
          if (StringICmp (oip->str, goFieldType [j]) == 0) break;
        }
        if (goFieldType [j] == NULL) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadGeneOntologyFormat, "Unrecognized label on GO term qualifier field %s", oip->str == NULL ? "[blank]" : oip->str);
          continue;
        }
        switch (j) {
          case 1 :
            if (ufp->choice == 1) {
              textstr = (CharPtr) ufp->data.ptrvalue;
            } else {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadGeneOntologyFormat, "Bad data format for GO term qualifier term");
            }
            break;
          case 2 :
            if (ufp->choice == 1) {
              goid = (CharPtr) ufp->data.ptrvalue;
            } else if (ufp->choice == 2) {
              sprintf (gid, "%ld", (long) (Int4) ufp->data.intvalue);
              goid = (CharPtr) gid;
            } else {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadGeneOntologyFormat, "Bad data format for GO term qualifier GO ID");
            }
            break;
          case 3 :
            if (ufp->choice == 2) {
              pmid = (Int4) ufp->data.intvalue;
            } else {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadGeneOntologyFormat, "Bad data format for GO term qualifier PMID");
            }
            break;
          case 5 :
            if (ufp->choice == 1) {
              evidence = (CharPtr) ufp->data.ptrvalue;
            } else {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadGeneOntologyFormat, "Bad data format for GO term qualifier evidence");
            }
            break;
          default :
            break;
        }
      }
      if (StringHasNoText (goid)) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_GeneOntologyTermMissingGOID, "GO term does not have GO identifier");
      }

      gsp = (GovStrucPtr) MemNew (sizeof (GovStruc));
      if (gsp != NULL) {
        gsp->term = StringSave (textstr);
        gsp->goid = StringSave (goid);
        gsp->evidence = StringSave (evidence);
        gsp->pmid = pmid;
        ValNodeAddPointer (&head, 0, (Pointer) gsp);
      }
    }
  }

  if (head == NULL || head->next == NULL) {
    return head;
  }
  head = ValNodeSort (head, SortVnpByGvspTermFirst);
  lastgsp = head->data.ptrvalue;
  for (vnp = head->next; vnp != NULL; vnp = vnp->next) {
    gsp = vnp->data.ptrvalue;
    if (StringICmp (gsp->term, lastgsp->term) == 0 || StringICmp (gsp->goid, lastgsp->goid) == 0) {
      if (gsp->pmid == lastgsp->pmid && StringICmp (gsp->evidence, lastgsp->evidence) == 0) {
        ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_DuplicateGeneOntologyTerm, "Duplicate GO term on feature");
      }
    }
    lastgsp = gsp;
  }
  return head;
}


static void ValidateGoTermUserObject (ValidStructPtr vsp, UserObjectPtr uop)
{
  ObjectIdPtr oip;
  UserFieldPtr ufp;
  ValNodePtr   term_list = NULL, vnp;
  GovStrucPtr  gsp, lastgsp;
  Int4 j;

  if (uop == NULL || vsp == NULL) return;
  oip = uop->type;
  if (oip == NULL) return;
  if (StringCmp (oip->str, "GeneOntology") == 0) {
    for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
      if (ufp->choice != 11 || ufp->label == NULL) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadGeneOntologyFormat, "Bad data format for GO term");
      } else {
        for (j = 0; goQualType [j] != NULL; j++) {
          if (StringICmp (ufp->label->str, goQualType [j]) == 0) break;
        }
        if (goQualType [j] == NULL) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadGeneOntologyFormat, "Unrecognized GO term label %s", ufp->label->str == NULL ? "[blank]" : ufp->label->str);
        } else {
          ValNodeLink (&term_list, ValidateGoTermQualifier (vsp, ufp->data.ptrvalue));
        }
      }
    }
    if (term_list != NULL) {
      term_list = ValNodeSort (term_list, SortVnpByGvspGoIDFirst);
      lastgsp = term_list->data.ptrvalue;
      for (vnp = term_list->next; vnp != NULL; vnp = vnp->next) {
        gsp = vnp->data.ptrvalue;
        if (gsp->goid != NULL
            && StringCmp (lastgsp->goid, gsp->goid) == 0
            && StringCmp (lastgsp->term, gsp->term) != 0) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_InconsistentGeneOntologyTermAndId, "Inconsistent GO terms for GO ID %s", gsp->goid);
        }
        lastgsp = gsp;
      }
      /* free term list */
      for (vnp = term_list; vnp != NULL; vnp = vnp->next) {
        gsp = vnp->data.ptrvalue;
        gsp->goid = MemFree (gsp->goid);
        gsp->term = MemFree (gsp->term);
        gsp->evidence = MemFree (gsp->evidence);
      }
    }
  }
}

static void LookForAccnLocs (SeqIdPtr sip, Pointer userdata)

{
  BoolPtr       bp;
  TextSeqIdPtr  tsip;

  if (sip == NULL || userdata == NULL) return;
  bp = (BoolPtr) userdata;

  switch (sip->choice) {
    case SEQID_GENBANK :
    case SEQID_EMBL :
    case SEQID_DDBJ :
    case SEQID_TPG :
    case SEQID_TPE :
    case SEQID_TPD :
    case SEQID_OTHER :
      tsip = (TextSeqIdPtr) sip->data.ptrvalue;
      if (tsip != NULL) {
        if (StringDoesHaveText (tsip->accession)) {
          *bp = TRUE;
        }
      }
      break;
    default :
      break;
  }
}

static Boolean LocationIsFar (SeqLocPtr location)

{
  BioseqPtr    bsp;
  DeltaSeqPtr  dsp;
  Boolean      is_far = FALSE;
  SeqLocPtr    loc;
  SeqEntryPtr  oldscope;
  SeqIdPtr     sip;
  SeqLocPtr    slp;

  if (location == NULL) return FALSE;

  oldscope = SeqEntrySetScope (NULL);

  slp = SeqLocFindNext (location, NULL);
  while (slp != NULL) {
    if (slp->choice != SEQLOC_NULL) {
      sip = SeqLocId (slp);
      bsp = BioseqFind (sip);
      if (bsp == NULL) {
        is_far = TRUE;
      } else if (bsp->repr == Seq_repr_delta && bsp->seq_ext_type == 4) {
        for (dsp = (DeltaSeqPtr) bsp->seq_ext;
             dsp != NULL && (! is_far);
             dsp = dsp->next) {
          if (dsp->choice != 1) continue;
          loc = (SeqLocPtr) dsp->data.ptrvalue;
          if (loc == NULL) continue;
          if (loc->choice == SEQLOC_NULL) continue;
          sip = SeqLocId (loc);
          bsp = BioseqFind (sip);
          if (bsp == NULL) {
            is_far = TRUE;
          }
        }
      } else if (bsp->repr == Seq_repr_seg && bsp->seq_ext_type == 1) {
        for (loc = (SeqLocPtr) bsp->seq_ext;
             loc != NULL && (! is_far);
             loc = loc->next) {
          if (loc == NULL) continue;
          if (loc->choice == SEQLOC_NULL) continue;
          sip = SeqLocId (loc);
          bsp = BioseqFind (sip);
          if (bsp == NULL) {
            is_far = TRUE;
          }
        }
      }
    }
    slp = SeqLocFindNext (location, slp);
  }

  SeqEntrySetScope (oldscope);

  return is_far;
}

static Boolean NoFetchFunctions (void)

{
  ObjMgrProcPtr  ompp = NULL;

  ompp = ObjMgrProcFindNext (NULL, OMPROC_FETCH, OBJ_SEQID, OBJ_BIOSEQ, NULL);

  return (Boolean) (ompp == NULL);
}

static CharPtr infMessage [] = {
  "unknown error",
  "empty inference string",
  "bad inference prefix",
  "bad inference body",
  "single inference field",
  "spaces in inference",
  "same species misused",
  "bad inference accession",
  "bad inference accession version",
  "accession.version not public",
  "bad accession type",
  NULL
};

static CharPtr rnaNameByType [] = {
  "unknown",
  "premsg",
  "mRNA",
  "tRNA",
  "rRNA",
  "snRNA",
  "scRNA",
  "snoRNA",
  "otherRNA",
  NULL
};

static Boolean ValStrandsMatch (Uint1 featstrand, Uint1 locstrand)

{
  if (featstrand == locstrand) return TRUE;
  if (locstrand == Seq_strand_unknown && featstrand != Seq_strand_minus) return TRUE;
  if (featstrand == Seq_strand_unknown && locstrand != Seq_strand_minus) return TRUE;
  if (featstrand == Seq_strand_both && locstrand != Seq_strand_minus) return TRUE;
  if (locstrand == Seq_strand_both) return TRUE;
  return FALSE;
}

static CharPtr badGeneSyn [] = {
  "alpha",
  "alternative",
  "beta",
  "cellular",
  "cytokine",
  "delta",
  "drosophila",
  "epsilon",
  "gamma",
  "HLA",
  "homolog",
  "mouse",
  "orf",
  "partial",
  "plasma",
  "precursor",
  "pseudogene",
  "putative",
  "rearranged",
  "small",
  "trna",
  "unknown",
  "unknown function",
  "unknown protein",
  "unnamed",
  NULL
};

static CharPtr badProtName [] = {
  "'hypothetical protein",
  "alpha",
  "alternative",
  "alternatively spliced",
  "bacteriophage hypothetical protein",
  "beta",
  "cellular",
  "cnserved hypothetical protein",
  "conesrved hypothetical protein",
  "conserevd hypothetical protein",
  "conserved archaeal protein",
  "conserved domain protein",
  "conserved hypohetical protein",
  "conserved hypotehtical protein",
  "conserved hypotheical protein",
  "conserved hypothertical protein",
  "conserved hypothetcial protein",
  "conserved hypothetical",
  "conserved hypothetical exported protein",
  "conserved hypothetical integral membrane protein",
  "conserved hypothetical membrane protein",
  "conserved hypothetical phage protein",
  "conserved hypothetical prophage protein",
  "conserved hypothetical protein",
  "conserved hypothetical protein - phage associated",
  "conserved hypothetical protein fragment 3",
  "conserved hypothetical protein, fragment",
  "conserved hypothetical protein, putative",
  "conserved hypothetical protein, truncated",
  "conserved hypothetical protein, truncation",
  "conserved hypothetical protein.",
  "conserved hypothetical protein; possible membrane protein",
  "conserved hypothetical protein; putative membrane protein",
  "conserved hypothetical proteins",
  "conserved hypothetical protien",
  "conserved hypothetical transmembrane protein",
  "conserved hypotheticcal protein",
  "conserved hypthetical protein",
  "conserved in bacteria",
  "conserved membrane protein",
  "conserved protein",
  "conserved protein of unknown function",
  "conserved protein of unknown function ; putative membrane protein",
  "conserved unknown protein",
  "conservedhypothetical protein",
  "conserverd hypothetical protein",
  "conservered hypothetical protein",
  "consrved hypothetical protein",
  "converved hypothetical protein",
  "cytokine",
  "delta",
  "drosophila",
  "duplicated hypothetical protein",
  "epsilon",
  "gamma",
  "HLA",
  "homeodomain",
  "homeodomain protein",
  "homolog",
  "hyopthetical protein",
  "hypotethical",
  "hypotheical protein",
  "hypothertical protein",
  "hypothetcical protein",
  "hypothetical",
  "hypothetical  protein",
  "hypothetical conserved protein",
  "hypothetical exported protein",
  "hypothetical novel protein",
  "hypothetical orf",
  "hypothetical phage protein",
  "hypothetical prophage protein",
  "hypothetical protein",
  "hypothetical protein (fragment)",
  "hypothetical protein (multi-domain)",
  "hypothetical protein (phage associated)",
  "hypothetical protein - phage associated",
  "hypothetical protein fragment",
  "hypothetical protein fragment 1",
  "hypothetical protein predicted by genemark",
  "hypothetical protein predicted by glimmer",
  "hypothetical protein predicted by glimmer/critica",
  "hypothetical protein, conserved",
  "hypothetical protein, phage associated",
  "hypothetical protein, truncated",
  "hypothetical protein-putative conserved hypothetical protein",
  "hypothetical protein.",
  "hypothetical proteins",
  "hypothetical protien",
  "hypothetical transmembrane protein",
  "hypothetoical protein",
  "hypothteical protein",
  "identified by sequence similarity; putative; ORF located~using Blastx/FrameD",
  "identified by sequence similarity; putative; ORF located~using Blastx/Glimmer/Genemark",
  "ion channel",
  "membrane protein, putative",
  "mouse",
  "narrowly conserved hypothetical protein",
  "novel protein",
  "orf",
  "orf, conserved hypothetical protein",
  "orf, hypothetical",
  "orf, hypothetical protein",
  "orf, hypothetical, fragment",
  "orf, partial conserved hypothetical protein",
  "orf; hypothetical protein",
  "orf; unknown function",
  "partial",
  "partial cds, hypothetical",
  "partially conserved hypothetical protein",
  "phage hypothetical protein",
  "phage-related conserved hypothetical protein",
  "phage-related protein",
  "plasma",
  "possible hypothetical protein",
  "precursor",
  "predicted coding region",
  "predicted protein",
  "predicted protein (pseudogene)",
  "predicted protein family",
  "product uncharacterised protein family",
  "protein family",
  "protein of unknown function",
  "pseudogene",
  "putative",
  "putative conserved protein",
  "putative exported protein",
  "putative hypothetical protein",
  "putative membrane protein",
  "putative orf; unknown function",
  "putative phage protein",
  "putative protein",
  "rearranged",
  "repeats containing protein",
  "reserved",
  "ribosomal protein",
  "similar to",
  "small",
  "small hypothetical protein",
  "transmembrane protein",
  "trna",
  "trp repeat",
  "trp-repeat protein",
  "truncated conserved hypothetical protein",
  "truncated hypothetical protein",
  "uncharacterized conserved membrane protein",
  "uncharacterized conserved protein",
  "uncharacterized conserved secreted protein",
  "uncharacterized protein",
  "uncharacterized protein conserved in archaea",
  "uncharacterized protein conserved in bacteria",
  "uniprot",
  "unique hypothetical",
  "unique hypothetical protein",
  "unknown",
  "unknown CDS",
  "unknown function",
  "unknown gene",
  "unknown protein",
  "unknown, conserved protein",
  "unknown, hypothetical",
  "unknown-related protein",
  "unknown; predicted coding region",
  "unnamed",
  "unnamed protein product",
  "very hypothetical protein",
  NULL
};

static Boolean NameInList (CharPtr name, CharPtr PNTR list, size_t numelements)

{
  Int2  L, R, mid;

  if (StringHasNoText (name) || list == NULL || numelements < 1) return FALSE;

  L = 0;
  R = numelements - 1; /* -1 because now NULL terminated */

  while (L < R) {
    mid = (L + R) / 2;
    if (StringICmp (list [mid], name) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  if (StringICmp (list [R], name) == 0) return TRUE;

  return FALSE;
}

static Boolean HasBadCharacter (CharPtr str)

{
  Char  ch;

  if (StringHasNoText (str)) return FALSE;

  ch = *str;
  while (ch != '\0') {
    if (ch == '?' || ch == '!' || ch == '~') return TRUE;
    str++;
    ch = *str;
  }

  return FALSE;
}

static Boolean EndsWithBadCharacter (CharPtr str)

{
  Char    ch;
  size_t  len;

  if (StringHasNoText (str)) return FALSE;

  len = StringLen (str);
  if (len < 1) return FALSE;

  ch = str [len - 1];
  if (ch == '_' || ch == '.' || ch == ',' || ch == ':' || ch == ';') return TRUE;

  return FALSE;
}

static Boolean EndsWithHyphen (CharPtr str)

{
  Char    ch;
  size_t  len;

  if (StringHasNoText (str)) return FALSE;

  len = StringLen (str);
  if (len < 1) return FALSE;

  ch = str [len - 1];
  if (ch == '-') return TRUE;

  return FALSE;
}


static Boolean CouldExtendPartial (SeqLocPtr slp, Boolean partial5)
{
  BioseqPtr bsp;
  Int4      pos;
  Uint1     strand;
  Char      str[4];
  Boolean   rval = FALSE;

  if (slp == NULL) {
    return FALSE;
  }

  bsp = BioseqFindFromSeqLoc (slp);
  if (bsp == NULL) {
    return FALSE;
  }
  strand = SeqLocStrand (slp);

  if ((strand != Seq_strand_minus && partial5) || (strand == Seq_strand_minus && !partial5)) {
    pos = SeqLocStart (slp);
    if (pos < 2) {
      rval = TRUE;
    } else if (bsp->repr == Seq_repr_delta) {
      /* wasn't close to the sequence end, but perhaps it is close to a gap */
      SeqPortStreamInt (bsp, pos - 3, pos - 1, Seq_strand_plus, EXPAND_GAPS_TO_DASHES, (Pointer) str, NULL);
      if (str[0] == '-' || str[1] == '-' || str[2] == '-') {
        rval = TRUE;
      }
    }
  } else {
    pos = SeqLocStop (slp);
    if (pos > bsp->length - 2) {
      rval = TRUE;
    } else {
      /* wasn't close to the sequence end, but perhaps it is close to a gap */
      SeqPortStreamInt (bsp, pos + 1, pos + 3, Seq_strand_plus, EXPAND_GAPS_TO_DASHES, (Pointer) str, NULL);
      if (str[0] == '-' || str[1] == '-' || str[2] == '-') {
        rval = TRUE;
      }
    }
  }

  return rval;
}

static Boolean LocationStrandsIncompatible (
  SeqLocPtr slp1,
  SeqLocPtr slp2
)

{
  Uint1  strand1, strand2;

  if (slp1 == NULL || slp2 == NULL) return FALSE;
  strand1 = SeqLocStrand (slp1);
  strand2 = SeqLocStrand (slp2);
  if (strand1 != strand2) {
    if ((strand1 == Seq_strand_unknown || strand1 == Seq_strand_plus) &&
        (strand2 == Seq_strand_unknown || strand2 == Seq_strand_plus)) return FALSE;
    /* if strands are mixed, need to check if interval is contained */
    if (strand1 == Seq_strand_other) {
      if (SeqLocCompareEx (slp1, slp2, TRUE) == SLC_B_IN_A) {
        return FALSE;
      }
    } else if (strand2 == Seq_strand_other) {
      if (SeqLocCompareEx (slp1, slp2, TRUE) == SLC_A_IN_B) {
        return FALSE;
      }
    }
    return TRUE;
  }
  return FALSE;
}


static CharPtr GetGeneXrefLabel (GeneRefPtr grp)
{
  SeqFeat sf;
  Char    buf[255];

  MemSet (&sf, 0, sizeof (SeqFeat));
  sf.data.choice = SEQFEAT_GENE;
  sf.data.value.ptrvalue = grp;

  FeatDefLabel (&sf, buf, sizeof (buf) - 1, OM_LABEL_CONTENT);
  return StringSave (buf);
}


static void TestForBracketsInProductName (CharPtr str, ValidStructPtr vsp)
{
  size_t  len;
  Boolean report_error = FALSE;
  CharPtr cp;

  if (StringHasNoText (str)) {
    return;
  }

  len = StringLen (str);
  if (len > 1) {
    if (str [len - 1] != ']') {
      /* doesn't end with bracket */
    } else if (len < 5) {
      /* too short to contain special text */
      report_error = TRUE;
    } else if ((cp = StringRChr (str, '[')) == NULL) {
      /* doesn't contain matched brackets */
      report_error = TRUE;
    } else if (StringNCmp (cp, "[NAD", 4) == 0) {
      /* contains special text */
    } else {
      report_error = TRUE;
    }
    if (report_error) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_ProteinNameEndsInBracket, "Protein name ends with bracket and may contain organism name");
    }
  }
}


static void ValidateRna (SeqFeatPtr sfp, ValidStructPtr vsp, GatherContextPtr gcp)
{
  RnaRefPtr rrp;
  Boolean   pseudo, ovgenepseudo = FALSE;
  Boolean   protidqual = FALSE, transidqual = FALSE;
  GBQualPtr gbq;
  tRNAPtr   trp;
  Boolean   badanticodon, anticodonqual, productqual;
  Int4      anticodonlen;
  SeqLocPtr slp;
  RNAGenPtr rgp;
  Int2      i;
  CharPtr   str;

    rrp = (RnaRefPtr) (sfp->data.value.ptrvalue);

    pseudo = sfp->pseudo;
    ovgenepseudo = FALSE;
    if (OverlappingGeneIsPseudo (sfp)) {
      pseudo = TRUE;
      ovgenepseudo = TRUE;
    }

    if (rrp->type == 2) {       /* mRNA */
      if (!pseudo) {
        MrnaTransCheck (vsp, sfp);      /* transcription check */
        SpliceCheck (vsp, sfp);
      }
      /* CheckForBothStrands (vsp, sfp); */
      CheckForBadGeneOverlap (vsp, sfp);
      CheckForCommonMRNAProduct (vsp, sfp);
      protidqual = FALSE;
      transidqual = FALSE;
      gbq = sfp->qual;
      while (gbq != NULL) {
        if (StringICmp (gbq->qual, "protein_id") == 0) {
          protidqual = TRUE;
        }
        if (StringICmp (gbq->qual, "transcript_id") == 0) {
          transidqual = TRUE;
        }
        gbq = gbq->next;
      }
      if (protidqual) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_WrongQualOnFeature, "protein_id should not be a gbqual on an mRNA feature");
      }
      if (transidqual) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_WrongQualOnFeature, "transcript_id should not be a gbqual on an mRNA feature");
      }
    }
    if (rrp->ext.choice == 2) { /* tRNA */
      trp = (tRNAPtr) (rrp->ext.value.ptrvalue);
      if (trp->anticodon != NULL) {
        badanticodon = FALSE;
        anticodonlen = 0;
        slp = SeqLocFindNext (trp->anticodon, NULL);
        while (slp != NULL) {
          anticodonlen += SeqLocLen (slp);
          i = SeqLocCompare (slp, sfp->location);
          if ((i != SLC_A_IN_B) && (i != SLC_A_EQ_B)) {
            badanticodon = TRUE;
          }
          slp = SeqLocFindNext (trp->anticodon, slp);
        }
        if (badanticodon) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_Range, "Anticodon location not in tRNA");
        }
        if (anticodonlen != 3) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_Range, "Anticodon is not 3 bases in length");
        }
        ValidateAnticodon (vsp, trp->anticodon);
      }
      CheckTrnaCodons (vsp, gcp, sfp, trp);
    }
    if (rrp->type == 3) {       /* tRNA */
      anticodonqual = FALSE;
      productqual = FALSE;
      gbq = sfp->qual;
      while (gbq != NULL) {
        if (StringICmp (gbq->qual, "anticodon") == 0) {
          anticodonqual = TRUE;
        } else if (StringICmp (gbq->qual, "product") == 0) {
          productqual = TRUE;
        }
        gbq = gbq->next;
      }
      if (anticodonqual) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "Unparsed anticodon qualifier in tRNA");
      }
      if (productqual) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "Unparsed product qualifier in tRNA");
      }
    }
    if (rrp->type == 3 && rrp->ext.choice == 1) { /* tRNA with string extension */
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "Unparsed product qualifier in tRNA");
    }
    if (rrp->type == 3 && rrp->ext.choice == 0) { /* tRNA with no extension */
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_MissingTrnaAA, "Missing encoded amino acid qualifier in tRNA");
    }
    if (rrp->type == 0) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_RNAtype0, "RNA type 0 (unknown) not supported");
    }
    if (rrp->type == 4 || rrp->type == 5 || rrp->type == 6 || rrp->type == 7) { /* rRNA, snRNA, scRNA, snoRNA */
      if (rrp->ext.choice != 1 || StringHasNoText ((CharPtr) rrp->ext.value.ptrvalue)) {
        if (! pseudo) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_InvalidQualifierValue, "%s has no name", rnaNameByType [(int) rrp->type]);
        }
      }
    }
    /*
    if (rrp->type == 255 && rrp->ext.choice == 1) {
      str = (CharPtr) rrp->ext.value.ptrvalue;
      if (StringICmp (str, "ncRNA") == 0) {
        for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
          if (StringICmp (gbq->qual, "ncRNA_class") != 0) continue;
          if (StringHasNoText (gbq->val)) continue;
          if (IsStringInNcRNAClassList (gbq->val)) continue;
          ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_InvalidQualifierValue, "Illegal ncRNA_class value '%s'", gbq->val);
        }
      }
    }
    */
    if (rrp->type == 2) {
      if (rrp->ext.choice == 1) {
        str = (CharPtr) rrp->ext.value.ptrvalue;
        if (StringDoesHaveText (str)) {
          if (HasBadCharacter (str)) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadInternalCharacter, "mRNA name contains undesired character");
          }
          if (EndsWithBadCharacter (str)) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadTrailingCharacter, "mRNA name ends with undesired character");
          }
          if (EndsWithHyphen (str)) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadTrailingHyphen, "mRNA name ends with hyphen");
          }
          if (StringHasSgml (vsp, str)) {
            ValidErr (vsp, SEV_WARNING, ERR_GENERIC_SgmlPresentInText, "mRNA name %s has SGML", str);
          }
        }
      }
    }
    if (rrp->type == 4) {
      if (rrp->ext.choice == 1) {
        str = (CharPtr) rrp->ext.value.ptrvalue;
        if (StringDoesHaveText (str)) {
          if (HasBadCharacter (str)) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadInternalCharacter, "rRNA name contains undesired character");
          }
          if (EndsWithBadCharacter (str)) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadTrailingCharacter, "rRNA name ends with undesired character");
          }
          if (EndsWithHyphen (str)) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadTrailingHyphen, "rRNA name ends with hyphen");
          }
          if (StringHasSgml (vsp, str)) {
            ValidErr (vsp, SEV_WARNING, ERR_GENERIC_SgmlPresentInText, "rRNA name %s has SGML", str);
          }
        }
      }
    }
    if (sfp->product != NULL) {
      CheckRnaProductType (vsp, gcp, sfp, rrp);
    }

    if (pseudo && sfp->product != NULL && StringISearch (sfp->except_text, "transcribed pseudogene") == NULL) {
      if (ovgenepseudo) {
        if (sfp->pseudo) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PseudoRnaHasProduct, "A pseudo RNA should not have a product");
        } else {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PseudoRnaViaGeneHasProduct, "An RNA overlapped by a pseudogene should not have a product");
        }
      } else {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PseudoRnaHasProduct, "A pseudo RNA should not have a product");
      }
    }

    if (rrp->ext.choice == 3
        && (rgp = (RNAGenPtr) rrp->ext.value.ptrvalue) != NULL
        && !StringHasNoText (rgp->_class)
        && rrp->type != 8) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_WrongQualOnFeature, "Only ncRNA should have ncRNA-class");
    }
}


NLM_EXTERN Boolean IsGeneXrefRedundant (SeqFeatPtr sfp)
{
  GeneRefPtr grp;
  SeqFeatPtr sfpx;
  GeneRefPtr grpx;
  Boolean    redundantgenexref = FALSE;
  CharPtr    syn1, syn2;
  DummySmfeData dsd;
  Int2          count;
  SeqMgrFeatContext fcontext;

  grp = SeqMgrGetGeneXref (sfp);
  if (grp == NULL) {
    return FALSE;
  }
  if (grp != NULL && SeqMgrGeneIsSuppressed (grp)) return FALSE;

  sfpx = SeqMgrGetOverlappingGene (sfp->location, &fcontext);
  if (sfpx == NULL || sfpx->data.choice != SEQFEAT_GENE)
    return FALSE;
  grpx = (GeneRefPtr) sfpx->data.value.ptrvalue;
  if (grpx == NULL)
    return FALSE;

  if (StringDoesHaveText (grp->locus_tag) && StringDoesHaveText (grp->locus_tag)) {
    if (StringICmp (grp->locus_tag, grpx->locus_tag) == 0) {
      redundantgenexref = TRUE;
    }
  } else if (StringDoesHaveText (grp->locus) && StringDoesHaveText (grp->locus)) {
    if (StringICmp (grp->locus, grpx->locus) == 0) {
      redundantgenexref = TRUE;
    }
  } else if (grp->syn != NULL && grpx->syn != NULL) {
    syn1 = (CharPtr) grp->syn->data.ptrvalue;
    syn2 = (CharPtr) grpx->syn->data.ptrvalue;
    if ((StringDoesHaveText (syn1)) && StringDoesHaveText (syn2)) {
      if (StringICmp (syn1, syn2) == 0) {
        redundantgenexref = TRUE;
      }
    }
  }
  if (redundantgenexref) {
    MemSet ((Pointer) &dsd, 0, sizeof (DummySmfeData));
    dsd.max = INT4_MAX;
    dsd.num_at_max = 0;
    dsd.num_trans_spliced = 0;
    dsd.equivalent_genes = FALSE;
    dsd.grp_at_max = NULL;
    count = SeqMgrGetAllOverlappingFeatures (sfp->location, FEATDEF_GENE, NULL, 0,
                                             LOCATION_SUBSET, (Pointer) &dsd, DummySMFEProc);
    if (dsd.num_at_max > 1) {
      redundantgenexref = FALSE;
    }
  }
  return redundantgenexref;
}


static void CheckCodingRegionAndProteinFeaturePartials (SeqFeatPtr sfp, ValidStructPtr vsp)
{
  BioseqPtr protbsp;
  SeqFeatPtr prot;
  SeqMgrFeatContext context;
  Boolean cds_partial5, cds_partial3, prot_partial5, prot_partial3, conflict = FALSE;
  SeqDescrPtr sdp;
  MolInfoPtr mip;
  Uint1 completeness;

  if (sfp == NULL || vsp == NULL) return;

  if (sfp->data.choice == SEQFEAT_CDREGION) {
    protbsp = BioseqFindFromSeqLoc (sfp->product);
    if (protbsp == NULL) return;
    prot = SeqMgrGetNextFeature (protbsp, NULL, 0, FEATDEF_PROT, &context);
    if (prot == NULL) return;
    CheckSeqLocForPartial (sfp->location, &cds_partial5, &cds_partial3);
    CheckSeqLocForPartial (prot->location, &prot_partial5, &prot_partial3);
    if ((cds_partial5 && !prot_partial5) || (!cds_partial5 && prot_partial5) ||
        (cds_partial3 && !prot_partial3) || (!cds_partial3 && prot_partial3)) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PartialsInconsistent, "Coding region and protein feature partials conflict");
    }
  } else if (sfp->data.choice == SEQFEAT_PROT) {
    protbsp = BioseqFindFromSeqLoc (sfp->location);
    if (protbsp == NULL) return;
    if (SeqMgrGetCDSgivenProduct (protbsp, &context) != NULL) return;
    prot = SeqMgrGetNextFeature (protbsp, NULL, 0, FEATDEF_PROT, &context);
    if (prot == NULL) return;
    sdp = GetNextDescriptorUnindexed (protbsp, Seq_descr_molinfo, NULL);
    if (sdp == NULL) return;
    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip == NULL) return;
    CheckSeqLocForPartial (prot->location, &prot_partial5, &prot_partial3);
    completeness = mip->completeness;
    if (completeness == 2 && ((! prot_partial5) && (! prot_partial3))) {
      conflict = TRUE;
    } else if (completeness == 3 && ((! prot_partial5) || prot_partial3)) {
      conflict = TRUE;
    } else if (completeness == 4 && (prot_partial5 || (! prot_partial3))) {
      conflict = TRUE;
    } else if (completeness == 5 && ((! prot_partial5) || (! prot_partial3))) {
      conflict = TRUE;
    } else if ((completeness < 2 || completeness > 5) && (prot_partial5 || prot_partial3)) {
      conflict = TRUE;
    }
    if (conflict) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PartialsInconsistent, "Molinfo completeness and protein feature partials conflict");
    }
  }
}


NLM_EXTERN void ValidateSeqFeat (GatherContextPtr gcp)
{
  Int2            type, i, j;
  static char    *parterrs[4] = {
    "Start does not include first/last residue of sequence",
    "Stop does not include first/last residue of sequence",
    "Internal partial intervals do not include first/last residue of sequence",
    "Improper use of partial (greater than or less than)"
  };
  Uint2           partials[2], errtype;
  Char            buf[80];
  CharPtr         tmp;
  ValidStructPtr  vsp;
  SeqFeatPtr      sfp;
  SeqFeatPtr      cds;
  CloneRefPtr     clrp;
  CdRegionPtr     crp;
  CodeBreakPtr    cbp, prevcbp;
  CharPtr         ctmp;
  GBQualPtr       gbq;
  Boolean         pseudo, excpt, conflict, codonqual,
                  protidqual,
                  transidqual, ovgenepseudo;
  ImpFeatPtr      ifp;
  GeneRefPtr      grp;
  SeqFeatPtr      gene;
  ProtRefPtr      prp;
  ValNodePtr      vnp;
  BioseqPtr       bsp, nbsp;
  BioseqContextPtr bcp = NULL;
  BioSourcePtr    biop, dbiop;
  OrgNamePtr      onp;
  OrgRefPtr       orp, dorp;
  SubSourcePtr    ssp;
  Boolean         transgenic;
  Int2            biopgencode;
  Int2            cdsgencode;
  Boolean         plastid;
  GeneticCodePtr  gc;
  PubdescPtr      pdp;
  /*
  DbtagPtr        db = NULL;
  Int4            id = -1;
  */
  SeqMgrDescContext context;
  GeneRefPtr      grpx;
  SeqFeatPtr      sfpx = NULL, sfpy = NULL, prt;
  SeqFeatPtr      operon;
  Boolean         redundantgenexref;
  SeqMgrFeatContext fcontext, gcontext;
  CharPtr         syn1, syn2, label = NULL, genexref_label;
  Uint2           oldEntityID;
  Uint4           oldItemID;
  SeqIdPtr        sip;
  TextSeqIdPtr    tsip;
  BioseqPtr       protBsp;
  ErrSev          sev;
  Boolean         multitoken;
  Char            ch;
  CharPtr         ptr;
  SeqLocPtr       slp;
  Int2            count;
  DummySmfeData   dsd;
  CharPtr         str;
  Boolean         isgap;
  Boolean         badseq;
  Boolean         is_seqloc_bond;
  SeqBondPtr      sbp;
  SeqFeatXrefPtr  xref, matchxref;
  SeqFeatPtr      matchsfp, origsfp;
  Boolean         hasxref;
  CharPtr         sfp_old_locus_tag;
  CharPtr         gene_old_locus_tag;
  Boolean         bypassGeneTest;
  Boolean         dicistronic = FALSE;
  Int2            inferenceCode;
  Boolean         hasInference = FALSE;
  Boolean         hasExperiment = FALSE;
  Boolean         accn_seqid;
  SeqDescrPtr     sdp;
  SeqMgrDescContext dcontext;
  MolInfoPtr      mip;
  Boolean         farFetchProd;
  Boolean         skip;
  Boolean         is_nc = FALSE;
  VariationRefPtr  vrfp;

  vsp = (ValidStructPtr) (gcp->userdata);
  sfp = (SeqFeatPtr) (gcp->thisitem);
  vsp->descr = NULL;
  vsp->sfp = sfp;
  type = (Int2) (sfp->data.choice);

  ValidateSeqLoc (vsp, sfp->location, "Location");

  ValidateSeqLoc (vsp, sfp->product, "Product");

  CheckForBothOrBothRev (vsp, sfp);

  if (vsp->feat_loc_has_gi) {
    accn_seqid = FALSE;
    VisitSeqIdsInSeqLoc (sfp->location, (Pointer) &accn_seqid, LookForAccnLocs);
    if (accn_seqid) {
      if (! vsp->is_smupd_in_sep) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FeatureRefersToAccession, "Feature location refers to accession");
      }
    }
  }

  if (vsp->feat_prod_has_gi) {
    accn_seqid = FALSE;
    VisitSeqIdsInSeqLoc (sfp->product, (Pointer) &accn_seqid, LookForAccnLocs);
    if (accn_seqid) {
      if (! vsp->is_smupd_in_sep) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FeatureRefersToAccession, "Feature product refers to accession");
      }
    }
  }

  farFetchProd = (Boolean) (vsp->farFetchCDSproducts || vsp->farFetchMRNAproducts);
  partials[0] = SeqLocPartialCheckEx (sfp->product, farFetchProd);
  partials[1] = SeqLocPartialCheck (sfp->location);

  CheckCodingRegionAndProteinFeaturePartials (sfp, vsp);

  if ((partials[0] != SLP_COMPLETE) || (partials[1] != SLP_COMPLETE) || (sfp->partial)) {       /* partialness */
    /* a feature on a partial sequence should be partial -- if often isn't */
    if ((!sfp->partial) && (partials[1] != SLP_COMPLETE) && (sfp->location->choice == SEQLOC_WHOLE)) {
      ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_PartialProblem, "On partial Bioseq, SeqFeat.partial should be TRUE");
    }
    /* a partial feature, with complete location, but partial product */
    else if ((sfp->partial) && (sfp->product != NULL) && (partials[1] == SLP_COMPLETE) && (sfp->product->choice == SEQLOC_WHOLE)
             && (partials[0] != SLP_COMPLETE)) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PartialProblem, "When SeqFeat.product is a partial Bioseq, SeqFeat.location should also be partial");
    }
    /* gene on segmented set is now 'order', should also be partial */
    else if (type == SEQFEAT_GENE && sfp->product == NULL && partials[1] == SLP_INTERNAL) {
      if (!sfp->partial) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PartialProblem, "Gene of 'order' with otherwise complete location should have partial flag set");
      }
    }
    /* inconsistent combination of partial/complete product,location,partial flag - part 1 */
    else if (((partials[0] == SLP_COMPLETE) && (sfp->product != NULL))) {
      sev = SEV_WARNING;
      bsp = GetBioseqGivenSeqLoc (sfp->product, gcp->entityID);
      /* if not local bioseq product, lower severity */
      if (bsp == NULL) {
        sev = SEV_INFO;
      }
      tmp = StringMove (buf, "Inconsistent: ");
      if (sfp->product != NULL) {
        tmp = StringMove (tmp, "Product= ");
        if (partials[0])
          tmp = StringMove (tmp, "partial, ");
        else
          tmp = StringMove (tmp, "complete, ");
      }
      tmp = StringMove (tmp, "Location= ");
      if (partials[1])
        tmp = StringMove (tmp, "partial, ");
      else
        tmp = StringMove (tmp, "complete, ");
      tmp = StringMove (tmp, "Feature.partial= ");
      if (sfp->partial)
        tmp = StringMove (tmp, "TRUE");
      else
        tmp = StringMove (tmp, "FALSE");
      if (bsp == NULL && LocationIsFar (sfp->product) && NoFetchFunctions ()) {
        vsp->far_fetch_failure = TRUE;
      } else {
        ValidErr (vsp, sev, ERR_SEQ_FEAT_PartialsInconsistent, buf);
      }
    /* inconsistent combination of partial/complete product,location,partial flag - part 2 */
    } else if ((partials[1] == SLP_COMPLETE) || (!sfp->partial)) {
      tmp = StringMove (buf, "Inconsistent: ");
      if (sfp->product != NULL) {
        tmp = StringMove (tmp, "Product= ");
        if (partials[0])
          tmp = StringMove (tmp, "partial, ");
        else
          tmp = StringMove (tmp, "complete, ");
      }
      tmp = StringMove (tmp, "Location= ");
      if (partials[1])
        tmp = StringMove (tmp, "partial, ");
      else
        tmp = StringMove (tmp, "complete, ");
      tmp = StringMove (tmp, "Feature.partial= ");
      if (sfp->partial)
        tmp = StringMove (tmp, "TRUE");
      else
        tmp = StringMove (tmp, "FALSE");
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PartialsInconsistent, buf);
    }
    /* 5' or 3' partial location giving unclassified partial product */
    else if (((partials [1] & SLP_START) != 0 || ((partials [1] & SLP_STOP) != 0)) && ((partials [0] & SLP_OTHER) != 0) && sfp->partial) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PartialProblem, "5' or 3' partial location should not have unclassified partial in product molinfo descriptor");
    }

    /* may have other error bits set as well */

    /* PartialProduct */
    errtype = SLP_NOSTART;
    for (j = 0; j < 4; j++) {
      bypassGeneTest = FALSE;
      if (partials[0] & errtype) {
        if (sfp->data.choice == SEQFEAT_CDREGION && sfp->excpt &&
                   StringStr (sfp->except_text, "rearrangement required for product") != NULL) {
        } else if (sfp->data.choice == SEQFEAT_CDREGION && j == 0) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PartialProblem,
                "PartialProduct: 5' partial is not at start AND is not at consensus splice site");
        } else if (sfp->data.choice == SEQFEAT_CDREGION && j == 1) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PartialProblem,
                "PartialProduct: 3' partial is not at stop AND is not at consensus splice site");
        } else {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PartialProblem,
            "PartialProduct: %s", parterrs[j]);
        }
      }
      errtype <<= 1;
    }

    /* PartialLocation */
    errtype = SLP_NOSTART;
    for (j = 0; j < 4; j++) {
      bypassGeneTest = FALSE;
      if (partials[1] & errtype) {
        if (j == 3) {
          if (LocationIsFar (sfp->location) && NoFetchFunctions ()) {
            vsp->far_fetch_failure = TRUE;
          } else if (sfp->data.choice == SEQFEAT_CDREGION && sfp->excpt &&
                     StringStr (sfp->except_text, "rearrangement required for product") != NULL) {
          } else {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PartialProblem,
              "PartialLocation: Improper use of partial (greater than or less than)");
          }
        } else if (j == 2) {
          if (LocationIsFar (sfp->location) && NoFetchFunctions ()) {
            vsp->far_fetch_failure = TRUE;
          } else {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PartialProblem,
              "PartialLocation: Internal partial intervals do not include first/last residue of sequence");
          }
        } else {
          if (IsCddFeat (sfp)) {
            /* suppresses  warning */
          } else if (sfp->data.choice == SEQFEAT_GENE && SameAsCDS (sfp, errtype, NULL)) {
          } else if (sfp->data.choice == SEQFEAT_GENE && SameAsMRNA (sfp, errtype)) {
          } else if (sfp->idx.subtype == FEATDEF_mRNA && SameAsCDS (sfp, errtype, &bypassGeneTest)) {
          } else if (sfp->idx.subtype == FEATDEF_mRNA && (! bypassGeneTest) && SameAsGene (sfp)) {
          } else if (sfp->idx.subtype == FEATDEF_exon && SameAsMRNA (sfp, errtype)) {
          } else if (LocationIsFar (sfp->location) && NoFetchFunctions ()) {
            vsp->far_fetch_failure = TRUE;
          } else if (sfp->data.choice == SEQFEAT_CDREGION && SameAsMRNA (sfp, errtype) &&
                     PartialAtSpliceSiteOrGap (vsp, sfp->location, errtype, &isgap, &badseq)) {
          } else if (PartialAtSpliceSiteOrGap (vsp, sfp->location, errtype, &isgap, &badseq)) {
            if (! isgap) {
              if (sfp->idx.subtype != FEATDEF_CDS || SplicingNotExpected (sfp)) {
                ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_PartialProblem,
                            "PartialLocation: %s (but is at consensus splice site)",
                            parterrs[j]);
              } else if (sfp->idx.subtype == FEATDEF_CDS) {
                bsp = BioseqFindFromSeqLoc (sfp->location);
                if (bsp != NULL) {
                  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
                  if (sdp != NULL) {
                    mip = (MolInfoPtr) sdp->data.ptrvalue;
                    if (mip != NULL) {
                      if (mip->biomol == MOLECULE_TYPE_MRNA) {
                        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PartialProblem,
                                  "PartialLocation: %s (but is at consensus splice site, but is on an mRNA that is already spliced)",
                                  parterrs[j]);
                      }
                    }
                  }
                }
              }
            }
          } else if (badseq) {
            ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_PartialProblem,
              "PartialLocation: %s (and is at bad sequence)",
              parterrs[j]);
          } else if (sfp->data.choice == SEQFEAT_CDREGION && sfp->excpt &&
                     StringStr (sfp->except_text, "rearrangement required for product") != NULL) {
          } else if (sfp->data.choice == SEQFEAT_CDREGION && j == 0) {
            if (PartialAtGapOrNs (vsp, sfp->location, errtype) && StringStr (sfp->comment, "coding region disrupted by sequencing gap") != NULL) {
            } else {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PartialProblem,
                    "PartialLocation: 5' partial is not at start AND is not at consensus splice site");
            }
          } else if (sfp->data.choice == SEQFEAT_CDREGION && j == 1) {
            if (PartialAtGapOrNs (vsp, sfp->location, errtype) && StringStr (sfp->comment, "coding region disrupted by sequencing gap") != NULL) {
            } else {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PartialProblem,
                    "PartialLocation: 3' partial is not at stop AND is not at consensus splice site");
            }
          } else if (j == 0) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PartialProblem,
              "PartialLocation: Start does not include first/last residue of sequence");
          } else if (j == 1) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PartialProblem,
              "PartialLocation: Stop does not include first/last residue of sequence");
          }
        }
      }
      errtype <<= 1;
    }

  }

  CheckForIllegalDbxref (vsp, gcp, sfp->dbxref);

  switch (type) {
  case 1:                      /* Gene-ref */
    grp = (GeneRefPtr) (sfp->data.value.ptrvalue);
    if (grp != NULL) {
      if (EmptyOrNullString (grp->locus) &&
          EmptyOrNullString (grp->allele) && EmptyOrNullString (grp->desc) &&
          EmptyOrNullString (grp->maploc) && EmptyOrNullString (grp->locus_tag) &&
          grp->db == NULL && grp->syn == NULL) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_GeneRefHasNoData, "There is a gene feature where all fields are empty");
      }
      if (StringDoesHaveText (grp->locus_tag)) {
        multitoken = FALSE;
        for (ptr = grp->locus_tag, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
          if (IS_WHITESP (ch)) {
            multitoken = TRUE;
          }
        }
        if (multitoken) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_LocusTagProblem, "Gene locus_tag '%s' should be a single word without any spaces", grp->locus_tag);
        }
        /* check for matching old_locus_tag */
        for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
          if (StringCmp (gbq->qual, "old_locus_tag") == 0 && StringCmp (grp->locus_tag, gbq->val) == 0) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_LocusTagProblem, "Gene locus_tag and old_locus_tag '%s' match", grp->locus_tag);
          }
        }
        if (StringDoesHaveText (grp->locus) && StringICmp (grp->locus, grp->locus_tag) == 0) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_LocusTagProblem, "Gene locus and locus_tag '%s' match", grp->locus);
        }
      }
      CheckForIllegalDbxref (vsp, gcp, grp->db);
      if (StringDoesHaveText (grp->allele)) {
        for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
          if (StringCmp (gbq->qual, "allele") == 0 && StringDoesHaveText (gbq->val)) {
            if (StringICmp (gbq->val, grp->allele) == 0) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_WrongQualOnFeature, "Redundant allele qualifier (%s) on gene", gbq->val);
            } else if (sfp->idx.subtype != FEATDEF_variation) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_WrongQualOnFeature, "Hidden allele qualifier (%s) on gene", gbq->val);
            }
          }
        }
      }
      /*
      for (vnp = grp->db; vnp != NULL; vnp = vnp->next) {
        id = -1;
        db = vnp->data.ptrvalue;
        if (db && db->db) {
          for (i = 0; i < DBNUM; i++) {
            if (StringCmp (db->db, dbtag[i]) == 0) {
              id = i;
              break;
            }
          }
          if (id == -1 || (type != SEQFEAT_CDREGION && id < 4)) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_IllegalDbXref, "Illegal db_xref type %s", db->db);
          }
        }
      }
      */
      if (grp->locus != NULL && sfp->comment != NULL) {
        if (StringCmp (grp->locus, sfp->comment) == 0) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_RedundantFields, "Comment has same value as gene locus");
        }
      }
      if (grp->locus != NULL) {
        if (HasBadCharacter (grp->locus)) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadInternalCharacter, "Gene locus contains undesired character");
        }
        if (EndsWithBadCharacter (grp->locus)) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadTrailingCharacter, "Gene locus ends with undesired character");
        }
        if (EndsWithHyphen (grp->locus)) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadTrailingHyphen, "Gene locus ends with hyphen");
        }
      }
      if (grp->locus_tag != NULL && sfp->comment != NULL) {
        if (StringCmp (grp->locus_tag, sfp->comment) == 0) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_RedundantFields, "Comment has same value as gene locus_tag");
        }
      }
      if (StringDoesHaveText (grp->locus_tag)) {
        for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
          if (StringCmp (gbq->qual, "old_locus_tag") == 0 && StringDoesHaveText (gbq->val)) {
            if (StringICmp (gbq->val, grp->locus_tag) == 0) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_RedundantFields, "old_locus_tag has same value as gene locus_tag");
            }
          }
        }
      }
      if (grp->syn != NULL && (vsp->is_refseq_in_sep /* || vsp->seqSubmitParent */)) {
        for (vnp = grp->syn; vnp != NULL; vnp = vnp->next) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (StringHasNoText (str)) continue;
          if (NameInList (str, badGeneSyn, sizeof (badGeneSyn) / sizeof (badGeneSyn [0]))) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UndesiredGeneSynonym, "Uninformative gene synonym '%s'", str);
          }
        }
      }
      if (grp->syn != NULL && (vsp->is_refseq_in_sep || vsp->seqSubmitParent)) {
        for (vnp = grp->syn; vnp != NULL; vnp = vnp->next) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (StringHasNoText (str)) continue;
          if (StringDoesHaveText (grp->locus) && StringCmp (grp->locus, str) == 0) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UndesiredGeneSynonym, "gene synonym has same value as gene locus");
          }
        }
      }
      if (grp->syn != NULL) {
        bsp = BioseqFindFromSeqLoc (sfp->location);
        for (vnp = grp->syn; vnp != NULL; vnp = vnp->next) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (StringHasNoText (str)) continue;
          sfpx = SeqMgrGetFeatureByLabel (bsp, str, SEQFEAT_GENE, 0, NULL);
          if (sfpx != NULL && sfpx != sfp) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_IdenticalGeneSymbolAndSynonym, "gene synonym has same value (%s) as locus of another gene feature", str);
          }
        }
      }
      if (StringDoesHaveText (grp->locus) && StringDoesHaveText (grp->desc) && StringCmp (grp->locus, grp->desc) == 0) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UndesiredGeneSynonym, "gene description has same value as gene locus");
      }
      if (StringHasNoText (grp->locus) && StringHasNoText (grp->desc) && grp->syn != NULL) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UndesiredGeneSynonym, "gene synonym without gene locus or description");
      }
      if (StringDoesHaveText (grp->desc) && StringStr (grp->desc,"..") != NULL) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_WrongQualOnFeature, "Possible location text (%s) on gene description", grp->desc);
      }
      /* - need to ignore if curated drosophila - add to vsp internal flags for efficiency?
      if (StringDoesHaveText (grp->locus)) {
        for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
          if (StringCmp (gbq->qual, "old_locus_tag") == 0 && StringDoesHaveText (gbq->val)) {
            if (StringICmp (gbq->val, grp->locus) == 0) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_RedundantFields, "old_locus_tag has same value as gene locus");
            }
          }
        }
      }
      */
      if (StringHasSgml (vsp, grp->locus)) {
        ValidErr (vsp, SEV_WARNING, ERR_GENERIC_SgmlPresentInText, "gene locus %s has SGML", grp->locus);
      }
      if (StringHasSgml (vsp, grp->locus_tag)) {
        ValidErr (vsp, SEV_WARNING, ERR_GENERIC_SgmlPresentInText, "gene locus_tag %s has SGML", grp->locus_tag);
      }
      if (StringHasSgml (vsp, grp->desc)) {
        ValidErr (vsp, SEV_WARNING, ERR_GENERIC_SgmlPresentInText, "gene description %s has SGML", grp->desc);
      }
      for (vnp = grp->syn; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasSgml (vsp, str)) {
          ValidErr (vsp, SEV_WARNING, ERR_GENERIC_SgmlPresentInText, "gene synonym %s has SGML", str);
        }
      }
      if (StringDoesHaveText (grp->locus)) {
        bsp = BioseqFindFromSeqLoc (sfp->location);
        sfpx = SeqMgrGetGeneByLocusTag (bsp, grp->locus, &fcontext);
        if (sfpx != NULL) {
          grpx = (GeneRefPtr) sfpx->data.value.ptrvalue;
          if (grpx != NULL) {
            if (grp == grpx) {
              /*
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_LocusCollidesWithLocusTag, "locus collides with locus_tag in same gene");
              */
            } else {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_LocusCollidesWithLocusTag, "locus collides with locus_tag in another gene");
            }
          }
        }
      }
    }
    break;
  case 2:                      /* Org-ref */
    break;
  case 3:                      /* Cdregion */
    pseudo = sfp->pseudo;       /* now also uses new feature pseudo flag */
    excpt = FALSE;
    conflict = FALSE;
    codonqual = FALSE;
    crp = (CdRegionPtr) (sfp->data.value.ptrvalue);
    if (crp != NULL) {
      conflict = crp->conflict;
    }
    protidqual = FALSE;
    transidqual = FALSE;
    ovgenepseudo = FALSE;
    gbq = sfp->qual;
    while (gbq != NULL) {
      if (StringICmp (gbq->qual, "pseudo") == 0) {
        pseudo = TRUE;
      }
      if (StringICmp (gbq->qual, "exception") == 0) {
        excpt = TRUE;
      }
      if (StringICmp (gbq->qual, "codon") == 0) {
        codonqual = TRUE;
      }
      if (StringICmp (gbq->qual, "protein_id") == 0) {
        protidqual = TRUE;
      }
      if (StringICmp (gbq->qual, "transcript_id") == 0) {
        transidqual = TRUE;
      }
      gbq = gbq->next;
    }
    if (protidqual) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_WrongQualOnFeature, "protein_id should not be a gbqual on a CDS feature");
    }
    if (transidqual) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_WrongQualOnFeature, "transcript_id should not be a gbqual on a CDS feature");
    }
    if (OverlappingGeneIsPseudo (sfp)) {
      pseudo = TRUE;
      ovgenepseudo = TRUE;
    }
    if ((!pseudo) && (!conflict)) {
      CdTransCheck (vsp, sfp);
      SpliceCheck (vsp, sfp);
    } else if (conflict) {
      CdConflictCheck (vsp, sfp);
    }
    CdsProductIdCheck (vsp, sfp);
    crp = (CdRegionPtr) (sfp->data.value.ptrvalue);
    if (crp != NULL) {
      if (crp->code_break != NULL && StringISearch (sfp->except_text, "RNA editing") != NULL) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_TranslExceptAndRnaEditing, "CDS has both RNA editing /exception and /transl_except qualifiers");
      }
      prevcbp = NULL;
      for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) {
        i = SeqLocCompare (cbp->loc, sfp->location);
        if ((i != SLC_A_IN_B) && (i != SLC_A_EQ_B)) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_Range, "Code-break location not in coding region");
        } else if (sfp->product != NULL) {
          slp = dnaLoc_to_aaLoc (sfp, cbp->loc, TRUE, NULL, TRUE);
          if (slp == NULL) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_Range, "Code-break location not in coding region - may be frame problem");
          }
          SeqLocFree (slp);
        }
        if (prevcbp != NULL) {
          i = SeqLocCompare (cbp->loc, prevcbp->loc);
          if (i == SLC_A_EQ_B) {
            ctmp = SeqLocPrint (cbp->loc);
            if (ctmp != NULL) {
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_DuplicateTranslExcept, "Multiple code-breaks at same location [%s]", ctmp);
            } else {
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_DuplicateTranslExcept, "Multiple code-breaks at same location");
            }
            MemFree (ctmp);
          }
        }
        prevcbp = cbp;
      }
      if (excpt && (!sfp->excpt)) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_ExceptInconsistent, "Exception flag should be set in coding region");
      }
      if (crp->orf && sfp->product != NULL) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_OrfCdsHasProduct, "An ORF coding region should not have a product");
      }
      if (pseudo && sfp->product != NULL) {
        if (ovgenepseudo) {
          if (sfp->pseudo) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_PseudoCdsHasProduct, "A pseudo coding region should not have a product");
          } else {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_PseudoCdsViaGeneHasProduct, "A coding region overlapped by a pseudogene should not have a product");
          }
        } else {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_PseudoCdsHasProduct, "A pseudo coding region should not have a product");
        }
      }
      if (pseudo && SeqMgrGetProtXref (sfp) != NULL) {
        if (NGorNT (vsp->sep, sfp->location, &is_nc) || IsEMBLAccn (vsp->sep, sfp->location)) {
          sev = SEV_WARNING;
        } else {
          sev = SEV_ERROR;
        }
        ValidErr (vsp, sev, ERR_SEQ_FEAT_PseudoCdsHasProtXref, "A pseudo coding region should not have a protein xref");
      }
      if (codonqual) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_CodonQualifierUsed, "Use the proper genetic code, if available, or set transl_excepts on specific codons");
      }
      biopgencode = 0;
      cdsgencode = 0;
      bsp = GetBioseqGivenSeqLoc (sfp->location, gcp->entityID);
      if (bsp != NULL) {
        vnp = NULL;
        if (vsp->useSeqMgrIndexes) {
          vnp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
        } else {
          bcp = BioseqContextNew (bsp);
          vnp = BioseqContextGetSeqDescr (bcp, Seq_descr_source, NULL, NULL);
        }
        if (vnp != NULL && vnp->data.ptrvalue != NULL) {
          plastid = FALSE;
          biop = (BioSourcePtr) vnp->data.ptrvalue;
          orp = biop->org;
          if (orp != NULL && orp->orgname != NULL) {
            onp = orp->orgname;
            if (biop->genome == GENOME_kinetoplast ||
                biop->genome == GENOME_mitochondrion ||
                biop->genome == GENOME_hydrogenosome) {
              biopgencode = onp->mgcode;
            } else if (biop->genome == GENOME_chloroplast ||
                       biop->genome == GENOME_chromoplast ||
                       biop->genome == GENOME_plastid ||
                       biop->genome == GENOME_cyanelle ||
                       biop->genome == GENOME_apicoplast ||
                       biop->genome == GENOME_leucoplast ||
                       biop->genome == GENOME_proplastid ||
                       biop->genome == GENOME_chromatophore) {
              if (onp->pgcode > 0) {
                biopgencode = onp->pgcode;
              } else {
                biopgencode = 11;
              }
              plastid = TRUE;
            } else {
              biopgencode = onp->gcode;
            }
            gc = crp->genetic_code;
            if (gc != NULL) {
              for (vnp = gc->data.ptrvalue; vnp != NULL; vnp = vnp->next) {
                if (vnp->choice == 2) {
                  cdsgencode = (Int2) vnp->data.intvalue;
                }
              }
            }
            if (biopgencode != cdsgencode) {
              if (! vsp->seqSubmitParent) { /* suppress when validator run from tbl2asn */
                if (plastid) {
                  ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_GenCodeMismatch,
                            "Genetic code conflict between CDS (code %d) and BioSource.genome biological context (%s) (uses code 11)", (int) cdsgencode, plastidtxt [biop->genome]);
                } else {
                  ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_GenCodeMismatch,
                            "Genetic code conflict between CDS (code %d) and BioSource (code %d)", (int) cdsgencode, (int) biopgencode);
                }
              }
            }
          }
        }
        if (!vsp->useSeqMgrIndexes) {
          BioseqContextFree (bcp);
        }
      }
    }
    /* CheckForBothStrands (vsp, sfp); */
    CheckForBadGeneOverlap (vsp, sfp);
    CheckForBadMRNAOverlap (vsp, sfp);
    CheckForCommonCDSProduct (vsp, sfp);
    CheckCDSPartial (vsp, sfp);
    if (StringDoesHaveText (sfp->comment)) {
      if (LookForECnumberPattern (sfp->comment)) {
        skip = FALSE;
        bsp = BioseqFindFromSeqLoc (sfp->product);
        if (bsp != NULL && ISA_aa (bsp->mol)) {
          prt = SeqMgrGetBestProteinFeature (bsp, NULL);
          if (prt != NULL && prt->data.choice == SEQFEAT_PROT) {
            prp = (ProtRefPtr) prt->data.value.ptrvalue;
            if (prp != NULL) {
              for (vnp = prp->ec; vnp != NULL; vnp = vnp->next) {
                str = (CharPtr) vnp->data.ptrvalue;
                if (StringHasNoText (str)) continue;
                if (StringStr (sfp->comment, str) != NULL) {
                  skip = TRUE;
                }
                skip = TRUE; /* now suppress even if EC numbers are different */
              }
            }
          }
        }
        if (! skip) {
          ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_EcNumberProblem, "Apparent EC number in CDS comment");
        }
      }
    }
    break;
  case 4:                      /* Prot-ref */
    prp = (ProtRefPtr) (sfp->data.value.ptrvalue);
    if (prp != NULL) {
      if (prp->processed != 3 && prp->processed != 4) {
        vnp = prp->name;
        if ((vnp == NULL || EmptyOrNullString ((CharPtr) vnp->data.ptrvalue)) &&
            EmptyOrNullString (prp->desc) && prp->ec == NULL && prp->activity == NULL && prp->db == NULL) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_ProtRefHasNoData, "There is a protein feature where all fields are empty");
        }
        if (vnp != NULL) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (StringDoesHaveText (str)) {
            TestForBracketsInProductName (str, vsp);
            if (StringNICmp (str, "hypothetical protein XP_", 24) == 0) {
              bsp = GetBioseqGivenSeqLoc (sfp->location, gcp->entityID);
              if (bsp != NULL) {
                for (sip = bsp->id; sip != NULL; sip = sip->next) {
                  if (sip->choice != SEQID_OTHER) continue;
                  tsip = (TextSeqIdPtr) sip->data.ptrvalue;
                  if (tsip == NULL) continue;
                  if (StringICmp (tsip->accession, str + 21) != 0) {
                    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_HpotheticalProteinMismatch, "Hypothetical protein reference does not match accession");
                  }
                }
              }
            }
            if (prp->ec != NULL) {
              if (StringCmp (str, "Hypothetical protein") == 0 ||
                  StringCmp (str, "hypothetical protein") == 0 ||
                  StringCmp (str, "Unknown protein") == 0 ||
                  StringCmp (str, "unknown protein") == 0) {
                ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadProteinName, "Unknown or hypothetical protein should not have EC number");
              }
            }
            if (LookForECnumberPattern (str)) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_EcNumberProblem, "Apparent EC number in protein title");
            }
            if (vsp->rubiscoTest && StringStr (str, "ribulose") != NULL && StringStr (str, "bisphosphate") != NULL) {
              if (StringStr (str, "methyltransferase") == NULL && StringStr (str, "activase") == NULL) {
                if (StringICmp (str, "ribulose-1,5-bisphosphate carboxylase/oxygenase") == 0) {
                  /* allow standard name without large or small subunit designation - later need kingdom test */
                } else if (StringICmp (str, "ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit") != 0 &&
                    StringICmp (str, "ribulose-1,5-bisphosphate carboxylase/oxygenase small subunit") != 0) {
                  ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_RubiscoProblem, "Nonstandard ribulose bisphosphate protein name");
                }
              }
            }
            if (StringHasPMID (str)) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_ProteinNameHasPMID, "Protein name has internal PMID");
            }
          }
          if (str != NULL && sfp->comment != NULL) {
            if (StringCmp (str, sfp->comment) == 0) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_RedundantFields, "Comment has same value as protein name");
            }
          }
          if (StringDoesHaveText (sfp->comment)) {
            if (LookForECnumberPattern (sfp->comment)) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_EcNumberProblem, "Apparent EC number in protein comment");
            }
          }
        }
      }
      CheckForIllegalDbxref (vsp, gcp, prp->db);
      /*
      for (vnp = prp->db; vnp != NULL; vnp = vnp->next) {
        id = -1;
        db = vnp->data.ptrvalue;
        if (db && db->db) {
          for (i = 0; i < DBNUM; i++) {
            if (StringCmp (db->db, dbtag[i]) == 0) {
              id = i;
              break;
            }
          }
          if (id == -1 || (type != SEQFEAT_CDREGION && id < 4)) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_IllegalDbXref, "Illegal db_xref type %s", db->db);
          }
        }
      }
      */
      if (prp->name == NULL && prp->processed != 3 && prp->processed != 4) {
        if (StringDoesHaveText (prp->desc)) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_NoNameForProtein, "Protein feature has description but no name");
        } else if (prp->activity != NULL) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_NoNameForProtein, "Protein feature has function but no name");
        } else if (prp->ec != NULL) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_NoNameForProtein, "Protein feature has EC number but no name");
        } else {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_NoNameForProtein, "Protein feature has no name");
        }
      }
      if (prp->desc != NULL && sfp->comment != NULL) {
        if (StringCmp (prp->desc, sfp->comment) == 0) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_RedundantFields, "Comment has same value as protein description");
        }
      }
      for (vnp = prp->ec; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringDoesHaveText (str)) {
          if (! ValidateECnumber (str)) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadEcNumberFormat, "%s is not in proper EC_number format", str);
          } else if (ECnumberNotInList (str)) {
            if (ECnumberWasDeleted (str)) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadEcNumberValue, "EC_number %s was deleted", str);
            } else if (ECnumberWasReplaced (str)) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadEcNumberValue, "EC_number %s was transferred and is no longer valid", str);
            } else {
              StringNCpy_0 (buf, str, sizeof (buf));
              ptr = StringChr (buf, 'n');
              if (ptr != NULL) {
                ch = ptr [1];
                if (IS_DIGIT (ch)) {
                  ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_BadEcNumberValue, "%s is not a legal preliminary value for qualifier EC_number", str);
                } else {
                  ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadEcNumberValue, "%s is not a legal value for qualifier EC_number", str);
                }
              } else {
                ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadEcNumberValue, "%s is not a legal value for qualifier EC_number", str);
              }
            }
          }
        } else {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_EcNumberProblem, "EC number should not be empty");
        }
      }
    }
    if (prp != NULL && prp->name != NULL && (vsp->is_refseq_in_sep /* || vsp->seqSubmitParent */)) {
      for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
          if (NameInList (str, badProtName, sizeof (badProtName) / sizeof (badProtName [0]))) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UndesiredProteinName, "Uninformative protein name '%s'", str);
          } else if (StringStr (str, "=") != NULL ||
                     StringStr (str, "~") != NULL ||
                     StringISearch (str, "uniprot") != NULL ||
                     StringISearch (str, "uniprotkb") != NULL ||
                     StringISearch (str, "pmid") != NULL ||
                     StringISearch (str, "dbxref") != NULL) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UndesiredProteinName, "Uninformative protein name '%s'", str);
          }
        }
      }
      if (prp != NULL && prp->name != NULL) {
        for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (StringHasNoText (str)) continue;
          if (HasBadCharacter (str)) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadInternalCharacter, "Protein name contains undesired character");
          }
          if (EndsWithBadCharacter (str)) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadTrailingCharacter, "Protein name ends with undesired character");
          }
          if (EndsWithHyphen (str)) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadTrailingHyphen, "Protein name ends with hyphen");
          }
        }
      }
      for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasSgml (vsp, str)) {
          ValidErr (vsp, SEV_WARNING, ERR_GENERIC_SgmlPresentInText, "protein name %s has SGML", str);
        }
      }
      if (StringHasSgml (vsp, prp->desc)) {
        ValidErr (vsp, SEV_WARNING, ERR_GENERIC_SgmlPresentInText, "protein description %s has SGML", prp->desc);
      }
    break;
  case 5:                      /* RNA-ref */
    ValidateRna(sfp, vsp, gcp);

    break;
  case 6:                      /* Pub */
    pdp = (PubdescPtr) sfp->data.value.ptrvalue;
    /*
       ValidatePubdesc (vsp, pdp);
     */
    break;
  case 7:                      /* Seq */
    break;
  case 8:                      /* Imp-feat */
    ifp = (ImpFeatPtr) sfp->data.value.ptrvalue;
    if (vsp->validateExons) {

      if (ifp != NULL && StringICmp (ifp->key, "exon") == 0 && (! sfp->pseudo)) {
        SpliceCheckEx (vsp, sfp, TRUE);
      }
    }
    if (ifp != NULL) {
      ValidateImpFeat (vsp, gcp, sfp, ifp);
    }
    break;
  case 9:                      /* Region */
    break;
  case 10:                     /* Comment */
    break;
  case 11:                     /* Bond */
    break;
  case 12:                     /* Site */
    break;
  case 13:                     /* Rsite-ref */
    break;
  case 14:                     /* User-object */
    break;
  case 15:                     /* TxInit */
    break;
  case 16:                     /* Numbering */
    break;
  case 17:                     /* Secondary Structure */
    break;
  case 18:                     /* NonStdRes */
    break;
  case 19:                     /* Heterogen */
    break;
  case 20:                     /* BioSource */
    biop = (BioSourcePtr) sfp->data.value.ptrvalue;
    if (biop != NULL && biop->is_focus) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_FocusOnBioSourceFeature, "Focus must be on BioSource descriptor, not BioSource feature.");
    }
    if (biop != NULL) {
      orp = biop->org;
      if (orp != NULL) {
        bsp = GetBioseqGivenSeqLoc (sfp->location, gcp->entityID);
        if (bsp != NULL) {
          vnp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
          if (vnp != NULL) {
            dbiop = (BioSourcePtr) vnp->data.ptrvalue;
            if (dbiop != NULL) {
              dorp = dbiop->org;
              if (dorp != NULL) {
                if (!StringHasNoText (orp->taxname)) {
                  if (StringICmp (orp->taxname, dorp->taxname) != 0) {
                    if (!dbiop->is_focus) {
                      transgenic = FALSE;
                      for (ssp = dbiop->subtype; ssp != NULL; ssp = ssp->next) {
                        if (ssp->subtype == SUBSRC_transgenic) {
                          transgenic = TRUE;
                        }
                      }
                      if (! transgenic) {
                        oldEntityID = gcp->entityID;
                        oldItemID = gcp->itemID;

                        gcp->entityID = context.entityID;
                        gcp->itemID = context.itemID;
                        gcp->thistype = OBJ_SEQDESC;

                        ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_BioSourceNeedsFocus,
                                  "BioSource descriptor must have focus or transgenic when BioSource feature with different taxname is present.");

                        gcp->entityID = oldEntityID;
                        gcp->itemID = oldItemID;
                        gcp->thistype = OBJ_SEQFEAT;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    /*
       ValidateBioSource (vsp, gcp, biop, sfp, NULL);
     */
    break;
  case 21:                     /* CloneRef */
    clrp = (CloneRefPtr) sfp->data.value.ptrvalue;
    if (clrp != NULL) {
    }
    break;
  case 22:                     /* VariationRef */
    vrfp = (VariationRefPtr) sfp->data.value.ptrvalue;
    if (vrfp != NULL) {
    }
    break;
  default:
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidType, "Invalid SeqFeat type [%d]", (int) (type));
    break;
  }
  if (type == SEQFEAT_HET) {
    /* heterogen can have mix of bonds with just "a" point specified */
    is_seqloc_bond = FALSE;
    slp = SeqLocFindNext (sfp->location, NULL);
    while (slp != NULL) {
      if (slp->choice == SEQLOC_BOND) {
        sbp = (SeqBondPtr) slp->data.ptrvalue;
        if (sbp != NULL) {
          if (sbp->a == NULL || sbp->b != NULL) {
            is_seqloc_bond = TRUE;
          }
        }
      }
      slp = SeqLocFindNext (sfp->location, slp);
    }
    if (is_seqloc_bond) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_ImproperBondLocation, "Bond location should only be on bond features");
    }
  } else if (type != SEQFEAT_BOND) {
    is_seqloc_bond = FALSE;
    slp = SeqLocFindNext (sfp->location, NULL);
    while (slp != NULL) {
      if (slp->choice == SEQLOC_BOND) {
        is_seqloc_bond = TRUE;
      }
      slp = SeqLocFindNext (sfp->location, slp);
    }
    if (is_seqloc_bond) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_ImproperBondLocation, "Bond location should only be on bond features");
    }
  }
  if (type != 8) {
    ValidateNonImpFeat (vsp, gcp, sfp);
  }
  if ((! sfp->excpt) && (! StringHasNoText (sfp->except_text))) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_ExceptInconsistent, "Exception text is present, but exception flag is not set");
  }
  if ((sfp->excpt) && (StringHasNoText (sfp->except_text))) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_ExceptInconsistent, "Exception flag is set, but exception text is empty");
  }
  if (! StringHasNoText (sfp->except_text)) {
    ValidateExceptText (vsp, gcp, sfp);
  }

  for (xref = sfp->xref; xref != NULL; xref = xref->next) {
    if (xref->id.choice == 0 && xref->data.choice == 0) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_SeqFeatXrefProblem, "SeqFeatXref with no id or data field");
    } else if (xref->id.choice != 0) {
      matchsfp = SeqMgrGetFeatureByFeatID (sfp->idx.entityID, NULL, NULL, xref, NULL);
      if (matchsfp == NULL) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_SeqFeatXrefFeatureMissing, "Cross-referenced feature cannot be found");
      } else {
        hasxref = FALSE;
        for (matchxref = matchsfp->xref; matchxref != NULL; matchxref = matchxref->next) {
          if (matchxref->id.choice != 0) {
            hasxref = TRUE;
            origsfp = SeqMgrGetFeatureByFeatID (matchsfp->idx.entityID, NULL, NULL, matchxref, NULL);
            if (origsfp != sfp) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_SeqFeatXrefNotReciprocal, "Cross-referenced feature does not link reciprocally");
            } else {
              if (sfp->idx.subtype == FEATDEF_CDS && matchsfp->idx.subtype == FEATDEF_mRNA) {
                /* okay */
              } else if (sfp->idx.subtype == FEATDEF_mRNA && matchsfp->idx.subtype == FEATDEF_CDS) {
                /* okay */
              } else {
                ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_SeqFeatXrefProblem, "Cross-references are not between CDS and mRNA pair");
              }
            }
          }
        }
        if (! hasxref) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_SeqFeatXrefProblem, "Cross-referenced feature does not have its own cross-reference");
        }
      }
    }
  }

  if (StringHasSgml (vsp, sfp->comment)) {
    ValidErr (vsp, SEV_WARNING, ERR_GENERIC_SgmlPresentInText, "feature comment %s has SGML", sfp->comment);
  }

  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    /* first check for anything other than replace */
    if (StringICmp (gbq->qual, "replace") != 0) {
      if (JustQuotes (gbq->val)) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_InvalidQualifierValue, "Qualifier other than replace has just quotation marks");
      }
    }
    /* now check specific gbual types */
    if (StringICmp (gbq->qual, "inference") == 0) {
      hasInference = TRUE;
      inferenceCode = ValidateInferenceQualifier (gbq->val, vsp->inferenceAccnCheck);
      if (inferenceCode != VALID_INFERENCE) {
        if (inferenceCode < VALID_INFERENCE || inferenceCode > BAD_ACCESSION_TYPE) {
          inferenceCode = VALID_INFERENCE;
        }
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_InvalidInferenceValue, "Inference qualifier problem - %s (%s)",
                  infMessage [(int) inferenceCode], (gbq->val != NULL)? gbq->val : "?");
      }
    } else if (StringICmp (gbq->qual, "experiment") == 0) {
      hasExperiment = TRUE;
    } else if (StringICmp (gbq->qual, "EC_number") == 0) {
      str = gbq->val;
      if (StringDoesHaveText (str)) {
        if (! ValidateECnumber (str)) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadEcNumberFormat, "%s is not in proper EC_number format", str);
        } else if (ECnumberNotInList (str)) {
          if (ECnumberWasDeleted (str)) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadEcNumberValue, "EC_number %s was deleted", str);
          } else if (ECnumberWasReplaced (str)) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadEcNumberValue, "EC_number %s was replaced", str);
          } else {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadEcNumberValue, "%s is not a legal value for qualifier EC_number", str);
          }
        }
      } else {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_EcNumberProblem, "EC number should not be empty");
      }
    } else if (StringICmp (gbq->qual, "old_locus_tag") == 0) {
      if (StringChr (gbq->val, ',') != NULL) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_LocusTagProblem,
                  "old_locus_tag has comma, may contain multiple values");
      }
      grp = SeqMgrGetGeneXref (sfp);
      if (grp == NULL) {
        if (sfp->data.choice == SEQFEAT_GENE) {
          grp = (GeneRefPtr) sfp->data.value.ptrvalue;
        } else {
          gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
          if (gene != NULL && ! gene->pseudo) {
            grp = (GeneRefPtr) gene->data.value.ptrvalue;
          }
        }
      }
      if (grp == NULL || SeqMgrGeneIsSuppressed (grp) || StringHasNoText (grp->locus_tag)) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_LocusTagProblem,
                  "old_locus_tag without inherited locus_tag");
      }
    }
    if (StringHasSgml (vsp, gbq->val)) {
      ValidErr (vsp, SEV_WARNING, ERR_GENERIC_SgmlPresentInText, "feature qualifier %s has SGML", gbq->val);
    }
  }
  if (sfp->exp_ev > 0 && (! hasInference) && (! hasExperiment) && (! vsp->feat_loc_has_gi)) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_InvalidInferenceValue,
              "Inference or experiment qualifier missing but obsolete experimental evidence qualifier set");
  }

  if (sfp->product != NULL) {
    sip = SeqLocId (sfp->product);
    if (sip != NULL) {
      switch (sip->choice) {
        case SEQID_LOCAL :
      break;
        case SEQID_GENBANK :
        case SEQID_EMBL :
        case SEQID_DDBJ :
        case SEQID_OTHER :
        case SEQID_TPG :
        case SEQID_TPE :
        case SEQID_TPD :
        case SEQID_GPIPE :
          tsip = (TextSeqIdPtr) sip->data.ptrvalue;
          if (tsip != NULL) {
            if (tsip->accession == NULL && (! StringHasNoText (tsip->name))) {
              if (ValidateAccn (tsip->name) == 0) {
                ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadProductSeqId,
                         "Feature product should not put an accession in the Textseq-id 'name' slot");
              } else {
                ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadProductSeqId,
                         "Feature product should not use Textseq-id 'name' slot");
              }
            }
          }
          break;
        default :
          break;
      }
    }
    bsp = BioseqFindFromSeqLoc (sfp->location);
    protBsp = BioseqFindFromSeqLoc (sfp->product);
    if (bsp != NULL && protBsp != NULL) {
      if (bsp == protBsp) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_SelfReferentialProduct, "Self-referential feature product");
      }
    }
    if (protBsp != NULL && protBsp->id != NULL) {
      for (sip = protBsp->id; sip != NULL; sip = sip->next) {
        switch (sip->choice) {
          case SEQID_GENBANK :
          case SEQID_EMBL :
          case SEQID_DDBJ :
          case SEQID_OTHER :
          case SEQID_TPG :
          case SEQID_TPE :
          case SEQID_TPD :
          case SEQID_GPIPE:
            tsip = (TextSeqIdPtr) sip->data.ptrvalue;
            if (tsip != NULL) {
              if (tsip->accession == NULL && (! StringHasNoText (tsip->name))) {
                if (ValidateAccn (tsip->name) == 0) {
                  ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadProductSeqId,
                            "Protein bioseq has Textseq-id 'name' that looks"
                            " like it is derived from a nucleotide accession");
                } else {
                  ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadProductSeqId,
                            "Protein bioseq has Textseq-id 'name' and no accession");
                }
              }
            }
            break;
          default :
            break;
        }
      }
    }

  }

  if (sfp->ext != NULL) {
    ValidateGoTermUserObject (vsp, sfp->ext);
  }

  if (type != SEQFEAT_GENE) {
    grp = SeqMgrGetGeneXref (sfp);
    if (grp == NULL) {
      sfpx = SeqMgrGetOverlappingGene (sfp->location, &fcontext);
      if (sfpx != NULL) {
        grp = (GeneRefPtr) sfpx->data.value.ptrvalue;
      }
    }
    if (grp != NULL && (! SeqMgrGeneIsSuppressed (grp))) {
      if (! StringHasNoText (grp->allele)) {
        for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
          if (StringCmp (gbq->qual, "allele") == 0 && StringDoesHaveText (gbq->val)) {
            if (StringICmp (gbq->val, grp->allele) == 0) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_InvalidQualifierValue, "Redundant allele qualifier (%s) on gene and feature", gbq->val);
            } else if (sfp->idx.subtype != FEATDEF_variation) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_InvalidQualifierValue, "Mismatched allele qualifier on gene (%s) and feature (%s)", grp->allele, gbq->val);
            }
          }
        }
      }
    }
    grp = SeqMgrGetGeneXref (sfp);
    if (grp != NULL && SeqMgrGeneIsSuppressed (grp)) return;

    if (grp == NULL) {
      sfpx = SeqMgrGetOverlappingGene (sfp->location, &fcontext);
      if (sfpx == NULL || sfpx->data.choice != SEQFEAT_GENE) return;
      sfp_old_locus_tag = NULL;
      gene_old_locus_tag = NULL;
      for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
        if (StringCmp (gbq->qual, "old_locus_tag") == 0 && StringDoesHaveText (gbq->val)) {
          sfp_old_locus_tag = gbq->val;
        }
      }
      for (gbq = sfpx->qual; gbq != NULL; gbq = gbq->next) {
        if (StringCmp (gbq->qual, "old_locus_tag") == 0 && StringDoesHaveText (gbq->val)) {
          gene_old_locus_tag = gbq->val;
        }
      }
      if (StringDoesHaveText (sfp_old_locus_tag) && StringDoesHaveText (gene_old_locus_tag)) {
        if (StringICmp (sfp_old_locus_tag, gene_old_locus_tag) != 0) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_OldLocusTagMismtach,
                    "Old locus tag on feature (%s) does not match that on gene (%s)",
                    sfp_old_locus_tag, gene_old_locus_tag);
        }
      }
      MemSet ((Pointer) &dsd, 0, sizeof (DummySmfeData));
      dsd.max = INT4_MAX;
      dsd.num_at_max = 0;
      dsd.num_trans_spliced = 0;
      dsd.equivalent_genes = FALSE;
      dsd.grp_at_max = NULL;
      count = SeqMgrGetAllOverlappingFeatures (sfp->location, FEATDEF_GENE, NULL, 0,
                                               LOCATION_SUBSET, (Pointer) &dsd, DummySMFEProc);
      if (dsd.num_at_max > 1 && sfp->idx.subtype != FEATDEF_repeat_region && sfp->idx.subtype != FEATDEF_mobile_element) {
        if (vsp->is_small_genome_set && dsd.num_at_max == dsd.num_trans_spliced) {
          /* suppress for trans-spliced genes on small genome set */
        } else if (dsd.equivalent_genes) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_GeneXrefNeeded,
                    "Feature overlapped by %d identical-length equivalent genes but has no cross-reference", (int) dsd.num_at_max);
        } else {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_MissingGeneXref,
                    "Feature overlapped by %d identical-length genes but has no cross-reference", (int) dsd.num_at_max);
        }
      }
      return;
    }

    if (StringDoesHaveText (grp->locus) /* && sfp->idx.subtype != FEATDEF_tRNA */) {
      bsp = BioseqFindFromSeqLoc (sfp->location);
      if (bsp != NULL) {
        sfpx = SeqMgrGetFeatureByLabel (bsp, grp->locus, SEQFEAT_GENE, 0, &fcontext);
        if (sfpx == NULL && ISA_aa (bsp->mol)) {
          cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
          if (cds != NULL) {
            nbsp = BioseqFindFromSeqLoc (cds->location);
            if (nbsp != NULL) {
              sfpx = SeqMgrGetFeatureByLabel (nbsp, grp->locus, SEQFEAT_GENE, 0, &fcontext);
            }
          }
        }
        if (sfpx != NULL) {
          sfpy = sfpx;
        }
        if (sfpx == NULL) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_GeneXrefWithoutGene,
                    "Feature has gene locus cross-reference but no equivalent gene feature exists");
        } else if (LocationStrandsIncompatible (sfp->location, sfpx->location)) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_GeneXrefStrandProblem,
                    "Gene cross-reference is not on expected strand");
        } else if (StringStr (sfpx->except_text, "dicistronic gene") != NULL) {
          dicistronic = TRUE;
        }
      }
    }
    if (StringDoesHaveText (grp->locus_tag)) {
      bsp = BioseqFindFromSeqLoc (sfp->location);
      if (bsp != NULL) {
        sfpx = SeqMgrGetGeneByLocusTag (bsp, grp->locus_tag, &fcontext);
        if (sfpx == NULL && ISA_aa (bsp->mol)) {
          cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
          if (cds != NULL) {
            nbsp = BioseqFindFromSeqLoc (cds->location);
            if (nbsp != NULL) {
              sfpx = SeqMgrGetFeatureByLabel (nbsp, grp->locus, SEQFEAT_GENE, 0, &fcontext);
            }
          }
        }
        if (sfpx != NULL) {
          sfpy = sfpx;
        }
        if (sfpx == NULL) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_GeneXrefWithoutGene,
                    "Feature has gene locus_tag cross-reference but no equivalent gene feature exists");
        } else if (LocationStrandsIncompatible (sfp->location, sfpx->location)) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_GeneXrefStrandProblem,
                    "Gene cross-reference is not on expected strand");
        } else if (StringStr (sfpx->except_text, "dicistronic gene") != NULL) {
          dicistronic = TRUE;
        }
        /* look for gene xrefs with locus_tag but no locus */
        if (StringHasNoText (grp->locus)
            && sfpx != NULL && sfpx->data.choice == SEQFEAT_GENE
            && sfpx->data.value.ptrvalue != NULL
            && StringDoesHaveText (((GeneRefPtr)sfpx->data.value.ptrvalue)->locus)) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_GeneXrefWithoutLocus,
                      "Feature has Gene Xref with locus_tag but no locus, gene with locus_tag and locus exists");
        }
      }
    }

    sfpx = NULL;
    if (SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL, 0, 0, sfp, &fcontext) == sfp) {
      if (fcontext.bad_order || fcontext.mixed_strand) {
        sfpx = SeqMgrGetOverlappingFeatureEx (sfp->location, FEATDEF_GENE, NULL, 0, NULL, LOCATION_SUBSET, &gcontext, TRUE);
      } else if (vsp->has_multi_int_genes) {
        sfpx = SeqMgrGetOverlappingFeatureEx (sfp->location, FEATDEF_GENE, NULL, 0, NULL, LOCATION_SUBSET, &gcontext, TRUE);
        if (sfpx == NULL && (vsp->has_seg_bioseqs || vsp->is_embl_ddbj_in_sep || vsp->is_old_gb_in_sep)) {
          sfpx = SeqMgrGetOverlappingGene (sfp->location, &gcontext);
        }
      } else {
        sfpx = SeqMgrGetOverlappingGene (sfp->location, &gcontext);
      }
    } else {
      sfpx = SeqMgrGetOverlappingGene (sfp->location, &gcontext);
    }
    if (sfpx == NULL || sfpx->data.choice != SEQFEAT_GENE)
      return;
    grpx = (GeneRefPtr) sfpx->data.value.ptrvalue;
    if (grpx == NULL)
      return;
    redundantgenexref = FALSE;
    label = gcontext.label;
    if (StringDoesHaveText (grp->locus_tag) && StringDoesHaveText (grp->locus_tag)) {
      if (StringICmp (grp->locus_tag, grpx->locus_tag) == 0) {
        redundantgenexref = TRUE;
        label = grp->locus_tag;
      }
    } else if (StringDoesHaveText (grp->locus) && StringDoesHaveText (grp->locus)) {
      if (StringICmp (grp->locus, grpx->locus) == 0) {
        redundantgenexref = TRUE;
        label = grp->locus;
      }
    } else if (grp->syn != NULL && grpx->syn != NULL) {
      syn1 = (CharPtr) grp->syn->data.ptrvalue;
      syn2 = (CharPtr) grpx->syn->data.ptrvalue;
      if ((StringDoesHaveText (syn1)) && StringDoesHaveText (syn2)) {
        if (StringICmp (syn1, syn2) == 0) {
          redundantgenexref = TRUE;
          label = syn1;
        }
      }
    }
    if (redundantgenexref) {
      MemSet ((Pointer) &dsd, 0, sizeof (DummySmfeData));
      dsd.max = INT4_MAX;
      dsd.num_at_max = 0;
      dsd.num_trans_spliced = 0;
      dsd.equivalent_genes = FALSE;
      dsd.grp_at_max = NULL;
      count = SeqMgrGetAllOverlappingFeatures (sfp->location, FEATDEF_GENE, NULL, 0,
                                               LOCATION_SUBSET, (Pointer) &dsd, DummySMFEProc);
      if (dsd.num_at_max > 1) {
        redundantgenexref = FALSE;
      }
    }
    if (redundantgenexref) {
      if (StringHasNoText (label)) {
        label = "?";
      }
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UnnecessaryGeneXref, "Unnecessary gene cross-reference %s", label);
    } else {
      if ((! dicistronic) && GPSorNTorNCorNGorNW (vsp->sep, sfp->location)) {
        /*
        SeqEntryToBioSource (vsp->sep, NULL, NULL, 0, &biop);
        */
        bsp = BioseqFindFromSeqLoc (sfp->location);
        BioseqToGeneticCode (bsp, NULL, NULL, NULL, NULL, 0, &biop);
        if (biop != NULL) {
          orp = biop->org;
          if (orp != NULL) {
            /* curated fly source still has duplicate features */
            if (StringICmp (orp->taxname, "Drosophila melanogaster") == 0) {
              if (StringHasNoText (label)) {
                label = "?";
              }
              if (sfpy != NULL && SeqLocAinB (sfp->location, sfpy->location) >= 0 &&
                  ValStrandsMatch (SeqLocStrand (sfp->location), SeqLocStrand (sfpy->location))) {
                /* cross-reference needed to disambiguate between multiple overlapping genes, ignore */
              } else {
                genexref_label = GetGeneXrefLabel (grp);
                ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_SuspiciousGeneXref, "Curated Drosophila record should not have gene cross-reference %s", genexref_label);
                genexref_label = MemFree (genexref_label);
              }
            }
          }
        }
      }
    }
  } else {
    grp = SeqMgrGetGeneXref (sfp);
    if (grp != NULL) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UnnecessaryGeneXref, "Gene feature has gene cross-reference");
    }
    operon = SeqMgrGetOverlappingOperon (sfp->location, &fcontext);
    if (operon != NULL) {
      if (SeqMgrGetDesiredFeature (sfp->idx.entityID, 0, 0, 0, sfp, &fcontext) == sfp) {
        if (! StringHasNoText (fcontext.label)) {
          for (gbq = operon->qual; gbq != NULL; gbq = gbq->next) {
            if (StringCmp (gbq->qual, "operon") == 0) {
              if (StringICmp (gbq->val, fcontext.label) == 0) {
                ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_InvalidQualifierValue, "Operon is same as gene - %s", gbq->val);
              }
            }
          }
        }
      }
    }
  }
}

/*****************************************************************************
*
*   MrnaTransCheck (sfp, vsp)
*
*****************************************************************************/

static CharPtr bypass_mrna_trans_check [] = {
  "RNA editing",
  "reasons given in citation",
  "artificial frameshift",
  "transcribed product replaced",
  "unclassified transcription discrepancy",
  "mismatches in transcription",
  "adjusted for low-quality genome",
  NULL
};

NLM_EXTERN void MrnaTransCheck (ValidStructPtr vsp, SeqFeatPtr sfp)
{
  BioseqPtr       bsp;
  Char            ch;
  Int4            counta, countnona;
  CharPtr         farstr = "";
  ErrSev          fetchsev;
  GatherContextPtr  gcp;
  Boolean         has_errors = FALSE, unclassified_except = FALSE,
                  mismatch_except = FALSE, other_than_mismatch = FALSE,
                  product_replaced = FALSE;
  Int2            i;
  Char            id [64];
  Boolean         is_refseq = FALSE;
  ErrSev          logsev;
  ErrSev          msgsev;
  Int4            mismatch, total;
  CharPtr         mrseq, pdseq;
  Int4            mlen, plen;
  CharPtr         ptr1, ptr2;
  Boolean         report_errors = TRUE;
  ErrSev          sev;
  SeqFeat         sf;
  SeqIdPtr        sip, sip2, sip3;
  Boolean         unlockProd = FALSE;
  ValNode         vn;
  SeqDescrPtr     sdp;
  MolInfoPtr      mip;
  TextSeqIdPtr    tsip;
  Boolean         rna_editing = FALSE;

  if (sfp == NULL)
    return;
  if (sfp->pseudo)
    return;
  if (sfp->product == NULL)
    return;

  if (sfp->excpt && (! vsp->ignoreExceptions) && (! StringHasNoText (sfp->except_text))) {
    for (i = 0; bypass_mrna_trans_check [i] != NULL; i++) {
      if (StringISearch (sfp->except_text,  bypass_mrna_trans_check [i]) != NULL) {
        report_errors = FALSE;  /* biological exception */
      }
    }
    if (StringISearch (sfp->except_text, "RNA editing") != NULL) {
      rna_editing = TRUE;
    }
    if (StringStr (sfp->except_text, "unclassified transcription discrepancy") != NULL) {
      unclassified_except = TRUE;
    }
    if (StringStr (sfp->except_text, "mismatches in transcription") != NULL) {
      mismatch_except = TRUE;
      report_errors = TRUE;
    }
    if (StringICmp (sfp->except_text, "transcribed product replaced") == 0) {
      product_replaced = TRUE;
    }
  }

  sip = SeqLocId (sfp->product);
  if (sip == NULL)
    return;

  msgsev = ErrSetMessageLevel (SEV_MAX);
  logsev = ErrSetLogLevel (SEV_MAX);

  mrseq = GetSequenceByFeature (sfp);

  ErrSetLogLevel (logsev);
  ErrSetMessageLevel (msgsev);

  if (mrseq == NULL) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_MrnaTransFail, "Unable to transcribe mRNA");
    return;
  }

  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp != NULL) {
    for (sip2 = bsp->id; sip2 != NULL; sip2 = sip2->next) {
      if (sip2->choice == SEQID_OTHER) {
        is_refseq = TRUE;
      }
    }
  }

  mismatch = 0;
  total = 0;

  sev = SEV_ERROR;
  gcp = vsp->gcp;
  if (gcp != NULL) {
    bsp = GetBioseqGivenSeqLoc (sfp->product, gcp->entityID);
    if (bsp == NULL) {
      /* if not local bioseq product, lower severity */
      sev = SEV_WARNING;
      if (is_refseq) {
        /* if refseq, restore higher severity */
        sev = SEV_ERROR;
      }
    }
    if (bsp == NULL && vsp->farFetchMRNAproducts) {
      if (sip != NULL && (sip->choice != SEQID_GI || sip->data.intvalue > 0)) {
        bsp = BioseqLockById (sip);
      }
      if (bsp != NULL) {
        unlockProd = TRUE;
        farstr = "(far) ";
        if (sfp->partial) {
          sdp = GetNextDescriptorUnindexed (bsp, Seq_descr_molinfo, NULL);
          if (sdp != NULL && sdp->choice == Seq_descr_molinfo) {
            mip = (MolInfoPtr) sdp->data.ptrvalue;
            if (mip != NULL) {
              if (mip->completeness < 2 || mip->completeness > 5) {
                for (sip3 = bsp->id; sip3 != NULL; sip3 = sip3->next) {
                  if (sip3->choice != SEQID_OTHER) continue;
                  tsip = (TextSeqIdPtr) sip3->data.ptrvalue;
                  if (tsip == NULL) continue;
                  if (StringNCmp (tsip->accession, "NM_", 3) == 0) {
                    /* if far NM_ record, return to lower severity */
                    sev = SEV_WARNING;
                  }
                }
              }
            }
          }
        }
      }
    }
    if (bsp == NULL && (! vsp->farFetchMRNAproducts)) {
      goto erret;
    }
    if (bsp == NULL && sfp->product != NULL && vsp->farFetchMRNAproducts) {
      SeqIdWrite (sip, id, PRINTID_FASTA_LONG, sizeof (id));
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_ProductFetchFailure, "Unable to fetch mRNA transcript '%s'", id);
      goto erret;
    }
  }
  if (is_refseq && unclassified_except) {
    /* if unclassified exception, drop back down to warning */
    sev = SEV_WARNING;
  }

  /* coerced feature on whole product for GetSequenceByFeature */

  MemSet ((Pointer) &sf, 0, sizeof (SeqFeat));
  MemSet ((Pointer) &vn, 0, sizeof (ValNode));
  sf.location = &vn;
  vn.choice = SEQLOC_WHOLE;
  vn.data.ptrvalue = sip;

  pdseq = GetSequenceByFeature (&sf);
  if (pdseq == NULL) {
    has_errors = TRUE;
    other_than_mismatch = TRUE;
    if (report_errors || unclassified_except) {
      fetchsev = SEV_ERROR;
      if (sip->choice != SEQID_GI) {
        fetchsev = SEV_WARNING;
      }
      ValidErr (vsp, fetchsev, ERR_SEQ_FEAT_MrnaTransFail, "Unable to fetch mRNA transcript");
    }
  }
  if (pdseq != NULL) {
    mlen = StringLen (mrseq);
    plen = StringLen (pdseq);
    if (mlen != plen) {
      if (mlen < plen) {
        ptr1 = pdseq + mlen;
        counta = 0;
        countnona = 0;
        ch = *ptr1;
        while (ch != '\0') {
          if (ch == 'A' || ch == 'a') {
            counta++;
          } else {
            countnona++;
          }
          ptr1++;
          ch = *ptr1;
        }
        if (counta < 19 * countnona) {
          has_errors = TRUE;
          other_than_mismatch = TRUE;
          if (report_errors || rna_editing) {
            ValidErr (vsp, sev, ERR_SEQ_FEAT_TranscriptLen, "Transcript length [%ld] less than %sproduct length [%ld], and tail < 95%s polyA", (long) mlen, farstr, (long) plen, "%");
          }
          plen = mlen; /* even if it fails polyA test, allow base-by-base comparison on common length */
        } else if (counta > 0 && countnona == 0) {
          has_errors = TRUE;
          other_than_mismatch = TRUE;
          if (report_errors || rna_editing) {
            ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_PolyATail, "Transcript length [%ld] less than %sproduct length [%ld], but tail is 100%s polyA", (long) mlen, farstr, (long) plen, "%");
          }
          plen = mlen; /* if it passes polyA test, allow base-by-base comparison on common length */
        } else {
          has_errors = TRUE;
          other_than_mismatch = TRUE;
          if (report_errors || rna_editing) {
            ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_PolyATail, "Transcript length [%ld] less than %sproduct length [%ld], but tail >= 95%s polyA", (long) mlen, farstr, (long) plen, "%");
          }
          plen = mlen; /* if it passes polyA test, allow base-by-base comparison on common length */
        }
      } else {
        has_errors = TRUE;
        other_than_mismatch = TRUE;
        if (report_errors || rna_editing) {
          ValidErr (vsp, sev, ERR_SEQ_FEAT_TranscriptLen, "Transcript length [%ld] greater than %sproduct length [%ld]", (long) mlen, farstr, (long) plen);
        }
      }
    }
    if (mlen == plen && mlen > 0 && StringICmp (mrseq, pdseq) != 0) {
      mismatch = 0;
      total = 0;
      ptr1 = mrseq;
      ptr2 = pdseq;
      while (total < mlen) {
        if (*ptr1 != *ptr2) {
          mismatch++;
        }
        ptr1++;
        ptr2++;
        total++;
      }
      if (mismatch > 0) {
        has_errors = TRUE;
        if (report_errors && (! mismatch_except)) {
          ValidErr (vsp, sev, ERR_SEQ_FEAT_TranscriptMismatches,
                    "There are %ld mismatches out of %ld bases between the transcript and %sproduct sequence", (long) mismatch, (long) total, farstr);
        }
      }
    }
    MemFree (pdseq);
  }

  if (! report_errors) {
    if (! has_errors) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UnnecessaryException, "mRNA has exception but passes transcription test");
    } else if (unclassified_except && (! other_than_mismatch)) {
      if (mismatch * 50 <= total) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_ErroneousException,
                  "mRNA has unclassified exception but only difference is %ld mismatches out of %ld bases",
                  (long) mismatch, (long) total);
      }
    } else if (product_replaced) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UnqualifiedException, "mRNA has unqualified transcribed product replaced exception");
    }
  }

erret:

  MemFree (mrseq);

  if (unlockProd) {
    BioseqUnlock (bsp);
  }

}

/*****************************************************************************
*
*   CdTransCheck(sfp)
*       Treatment of terminal 'X'
*          If either the protein or the translation end in 'X' (usually
*          due to partial last codon) it is ignored to minimize conflicts
*          between approaches to add the X or not in this case.
*
*****************************************************************************/
static CharPtr MapToNTCoords (SeqFeatPtr sfp, SeqIdPtr protID, Int4 pos)
{
  SeqLocPtr       nslp;
  SeqLocPtr       pslp;
  CharPtr         rsult;
  SeqPntPtr       spntp;

  rsult = NULL;
  if (sfp != NULL && protID != NULL && pos >= 0) {
    spntp = SeqPntNew ();
    pslp = ValNodeNew (NULL);
    pslp->choice = SEQLOC_PNT;
    pslp->data.ptrvalue = (Pointer) spntp;
    spntp->point = pos;
    spntp->id = SeqIdDup (protID);
    nslp = aaLoc_to_dnaLoc (sfp, pslp);
    if (nslp != NULL) {
      rsult = SeqLocPrint (nslp);
    }
    SeqLocFree (pslp);
    SeqLocFree (nslp);
  }
  return rsult;
}

static Boolean Loc_is_RefSeq (SeqLocPtr location)
{
  BioseqPtr       bsp;
  SeqIdPtr        sip;
  TextSeqIdPtr    tsip;

  if (location == NULL)
    return FALSE;
  sip = SeqLocId (location);
  if (sip == NULL)
    return FALSE;
  bsp = BioseqFind (sip);
  if (bsp == NULL)
    return FALSE;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_OTHER) {
      tsip = (TextSeqIdPtr) sip->data.ptrvalue;
      if (tsip != NULL) {
        if (StringNICmp (tsip->accession, "NM_", 3) == 0) {
          return TRUE;
        }
      }
    }
  }
  return FALSE;
}

static Boolean Loc_is_GEDL (SeqLocPtr location)
{
  BioseqPtr  bsp;
  SeqIdPtr   sip;

  if (location == NULL)
    return FALSE;
  sip = SeqLocId (location);
  if (sip == NULL)
    return FALSE;
  bsp = BioseqFind (sip);
  if (bsp == NULL)
    return FALSE;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GENBANK) return TRUE;
    if (sip->choice == SEQID_EMBL) return TRUE;
    if (sip->choice == SEQID_DDBJ) return TRUE;
    if (sip->choice == SEQID_LOCAL) return TRUE;
  }
  return FALSE;
}

static void CdConflictCheck (ValidStructPtr vsp, SeqFeatPtr sfp)

{
  ByteStorePtr  bs;
  BioseqPtr     bsp;
  CharPtr       str1, str2;

  if (sfp == NULL || vsp == NULL) return;

  bsp = BioseqFindFromSeqLoc (sfp->product);
  str1 = GetSequenceByBsp (bsp);
  bs = TransTableTranslateCdRegion (NULL, sfp, FALSE, FALSE, TRUE);
  str2 = (CharPtr) BSMerge (bs, NULL);
  BSFree (bs);

  if (str1 != NULL && str2 != NULL && StringCmp (str1, str2) == 0) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_BadConflictFlag, "Coding region conflict flag should not be set");
  } else {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_ConflictFlagSet, "Coding region conflict flag is set");
  }

  MemFree (str1);
  MemFree (str2);
}

static CharPtr bypass_cds_trans_check [] = {
  "RNA editing",
  "reasons given in citation",
  "artificial frameshift",
  "rearrangement required for product",
  "translated product replaced",
  "unclassified translation discrepancy",
  "mismatches in translation",
  "adjusted for low-quality genome",
  "annotated by transcript or proteomic data",
  /*
  "heterogeneous population sequenced",
  "low-quality sequence region",
  "artificial location",
  */
  NULL
};

static void ValidateTranslExcept (
  ValidStructPtr vsp,
  SeqFeatPtr sfp,
  ValNodePtr codebreakhead,
  Boolean farFetchProd,
  Uint1 frame,
  ValNodePtr genetic_code
)

{
  Boolean       alt_start = FALSE;
  CdRegion      cr;
  ByteStorePtr  newprot = NULL;
  CharPtr       protseq = NULL;
  Int4          prot2len, i;
  SeqFeat       sf;
  ValNodePtr    vnp;

  MemSet ((Pointer) &sf, 0, sizeof (SeqFeat));
  MemSet ((Pointer) &cr, 0, sizeof (CdRegion));
  sf.data.choice = SEQFEAT_CDREGION;
  sf.data.value.ptrvalue = (Pointer) &cr;
  sf.location = sfp->location;
  cr.frame = frame;
  cr.genetic_code = genetic_code;

  newprot = ProteinFromCdRegionExEx (&sf, TRUE, FALSE, &alt_start, farFetchProd);
  if (newprot == NULL) return;
  protseq = BSMerge (newprot, NULL);
  BSFree (newprot);
  if (protseq == NULL) return;
  prot2len = StringLen (protseq);
  for (vnp = codebreakhead; vnp != NULL; vnp = vnp->next) {
    i = vnp->data.intvalue;
    if (i >= 0 && i < prot2len) {
      if (protseq [i] == (Char) vnp->choice) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UnnecessaryTranslExcept,
                  "Unnecessary transl_except %c at position %ld",
                  (char) vnp->choice, (long) (i + 1));
      }
    } else if (i == prot2len) {
      if ((Char) vnp->choice != '*') {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UnnecessaryTranslExcept,
                  "Unexpected transl_except %c at position %ld just past end of protein",
                  (char) vnp->choice, (long) (i + 1));
      }
    }
  }
  MemFree (protseq);
}

typedef struct cdsmismatch {
  Int4 pos;
  Int2 cds_residue;
  Int2 prot_residue;
} CDSMismatchData, PNTR CDSMismatchPtr;

NLM_EXTERN void CdTransCheck (ValidStructPtr vsp, SeqFeatPtr sfp)

{
  ByteStorePtr    newprot = NULL;
  CharPtr         protseq = NULL;
  BioseqPtr       prot1seq = NULL, prot2seq = NULL;
  Int4            i, len = 0, x_count = 0,
                  nonx_count = 0, xcount1 = 0, xcount2 = 0;
  Int4            prot1len = 0; /* length of protein product sequence */
  Int4            prot2len = 0; /* length of translation of coding region */
  CdRegionPtr     crp;
  SeqIdPtr        protid = NULL;
  Int2            residue1, residue2, stop_count = 0, mismatch = 0, ragged = 0;
  CDSMismatchData mismatches[11];
  Boolean         got_stop = FALSE;
  /*
  SeqPortPtr      spp = NULL;
  */
  Uint2           part_loc = 0, part_prod = 0;
  Boolean         no_end = FALSE, no_beg = FALSE, show_stop = FALSE,
                  got_dash = FALSE, alt_start = FALSE, got_x = FALSE, done;
  GBQualPtr       gb;
  ValNodePtr      vnp, vnp2, code, codebreakhead = NULL;
  int             gccode = 0;
  Boolean         transl_except = FALSE, prot_ok = TRUE, is_nc = FALSE,
                  has_errors = FALSE, report_errors = TRUE,
                  unclassified_except = FALSE, mismatch_except = FALSE,
                  frameshift_except = FALSE, rearrange_except = FALSE,
                  other_than_mismatch = FALSE, product_replaced = FALSE,
                  mixed_population = FALSE, low_quality = FALSE,
                  artificial_location = FALSE;
  Boolean         partial5 = FALSE;
  Boolean         partial3 = FALSE;
  Boolean         rna_editing = FALSE;
  CharPtr         nuclocstr, farstr = "", loc2str;
  CodeBreakPtr    cbp;
  Int4            pos1, pos2, pos;
  SeqLocPtr       tmp;
  ErrSev          sev, trans_len_sev = SEV_ERROR;
  SeqEntryPtr     sep;
  Boolean         unlockProd = FALSE;
  StreamCache     sc;
  Boolean         isgap;
  Boolean         badseq;
  BioseqPtr       bsp;
  SeqIdPtr        sip, sip3;
  Char            id [64];
  Boolean         is_ged = FALSE;
  Boolean         is_refseq = FALSE;
  Boolean         has_gi = FALSE;
  Boolean         farFetchProd;
  SeqDescrPtr     sdp;
  MolInfoPtr      mip;
  TextSeqIdPtr    tsip;
  Boolean         annotated_by_transcript_or_proteomic = FALSE;

  if (sfp == NULL) return;

  crp = (CdRegionPtr) (sfp->data.value.ptrvalue);
  if (crp == NULL) return;

  for (gb = sfp->qual; gb != NULL; gb = gb->next) {     /* pseuogene */
    if (!StringICmp ("pseudo", gb->qual))
      return;
  }

  if (LocationIsFar (sfp->location) && NoFetchFunctions ()) {
    vsp->far_fetch_failure = TRUE;
    return;
  }

  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);

  if (sfp->excpt && (! vsp->ignoreExceptions) && (! StringHasNoText (sfp->except_text))) {
    for (i = 0; bypass_cds_trans_check [i] != NULL; i++) {
      if (StringISearch (sfp->except_text,  bypass_cds_trans_check [i]) != NULL) {
        report_errors = FALSE;  /* biological exception */
      }
    }
    if (StringStr (sfp->except_text, "unclassified translation discrepancy") != NULL) {
      unclassified_except = TRUE;
    }
    if (StringStr (sfp->except_text, "mismatches in translation") != NULL) {
      mismatch_except = TRUE;
      report_errors = TRUE;
    }
    if (StringStr (sfp->except_text, "artificial frameshift") != NULL) {
      frameshift_except = TRUE;
      report_errors = TRUE;
    }
    if (StringStr (sfp->except_text, "rearrangement required for product") != NULL) {
      rearrange_except = TRUE;
    }
    if (StringICmp (sfp->except_text, "translated product replaced") == 0) {
      product_replaced = TRUE;
    }
    if (StringICmp (sfp->except_text, "heterogeneous population sequenced") == 0) {
      mixed_population = TRUE;
    }
    if (StringICmp (sfp->except_text, "low-quality sequence region") == 0) {
      low_quality = TRUE;
    }
    if (StringICmp (sfp->except_text, "artificial location") == 0) {
      artificial_location = TRUE;
    }
    if (StringISearch (sfp->except_text, "RNA editing") != NULL) {
      rna_editing = TRUE;
    }
  }
  if (StringISearch (sfp->except_text, "annotated by transcript or proteomic data") != NULL) {
    annotated_by_transcript_or_proteomic = TRUE;
  }

  if (crp->code_break == NULL) {        /* check for unparsed transl_except */
    for (gb = sfp->qual; gb != NULL; gb = gb->next) {
      if (StringCmp (gb->qual, "transl_except") == 0) {
        transl_except = TRUE;
        break;
      }
    }
  } else {
    codebreakhead = MakeCodeBreakList (sfp->location, SeqLocLen (sfp->location), crp->code_break, crp->frame);
  }

  if (crp->genetic_code != NULL) {
    for (vnp = crp->genetic_code->data.ptrvalue; ((vnp != NULL) && (!gccode)); vnp = vnp->next) {
      switch (vnp->choice) {
      case 0:
        break;
      case 1:                  /* name */
        code = GeneticCodeFind (0, (CharPtr) (vnp->data.ptrvalue));
        if (code != NULL) {
          for (vnp2 = code->data.ptrvalue; ((vnp2 != NULL) && (!gccode)); vnp2 = vnp2->next) {
            if (vnp2->choice == 2)       /* id */
              gccode = (int) (vnp2->data.intvalue);
          }
        }
        break;
      case 2:                  /* id */
        gccode = (int) (vnp->data.intvalue);
        break;
      default:
        gccode = 255;
        break;
      }
    }
  }

  farFetchProd = (Boolean) (vsp->farFetchCDSproducts || vsp->farFetchMRNAproducts);
  newprot = ProteinFromCdRegionExEx (sfp, TRUE, FALSE, &alt_start, farFetchProd);   /* include stop codons, do not remove trailing X/B/Z */
  if (newprot == NULL) {
    has_errors = TRUE;
    other_than_mismatch = TRUE;
    if (report_errors || unclassified_except) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_CdTransFail, "Unable to translate");
    }
    prot_ok = FALSE;
    goto erret;
  }

  if (codebreakhead != NULL) {
    ValidateTranslExcept (vsp, sfp, codebreakhead, farFetchProd, crp->frame, crp->genetic_code);
  }

  protid = SeqLocId (sfp->product);
  if (protid != NULL) {
    prot1seq = BioseqFind (protid);
    if (prot1seq == NULL && vsp->farFetchCDSproducts) {
      if (protid != NULL && (protid->choice != SEQID_GI || protid->data.intvalue > 0)) {
        prot1seq = BioseqLockById (protid);
      }
      if (prot1seq != NULL) {
        unlockProd = TRUE;
        farstr = "(far) ";
        if (sfp->partial) {
          sdp = GetNextDescriptorUnindexed (prot1seq, Seq_descr_molinfo, NULL);
          if (sdp != NULL && sdp->choice == Seq_descr_molinfo) {
            mip = (MolInfoPtr) sdp->data.ptrvalue;
            if (mip != NULL) {
              if (mip->completeness < 2 || mip->completeness > 5) {
                for (sip3 = prot1seq->id; sip3 != NULL; sip3 = sip3->next) {
                  if (sip3->choice != SEQID_OTHER) continue;
                  tsip = (TextSeqIdPtr) sip3->data.ptrvalue;
                  if (tsip == NULL) continue;
                  if (StringNCmp (tsip->accession, "NP_", 3) == 0) {
                    /* if far NP_ record, return to lower severity */
                    trans_len_sev = SEV_WARNING;
                  }
                }
              }
            }
          }
        }
      }
    }
    if (prot1seq != NULL)
      prot1len = prot1seq->length;
  }

  if (alt_start && gccode == 1) {
    /* sev = SEV_WARNING; */
    sev = SEV_NONE; /* only enable for RefSeq, leave old code in for now */
    if (Loc_is_RefSeq (sfp->location)) {
      sev = /* SEV_ERROR */ SEV_NONE; /* now also disable for RefSeq */
    } else if (Loc_is_GEDL (sfp->location)) {
      sev = SEV_NONE;
    }
    if (sfp->excpt && StringDoesHaveText (sfp->except_text)) {
      if (StringStr (sfp->except_text, "alternative start codon") != NULL) {
        sev = SEV_NONE;
      }
    }
    if (sev > SEV_NONE) {
      has_errors = TRUE;
      other_than_mismatch = TRUE;
      if (report_errors) {
        ValidErr (vsp, sev, ERR_SEQ_FEAT_AltStartCodon, "Alternative start codon used");
      }
    }
  } else if (! alt_start) {
    if (sfp->excpt && StringDoesHaveText (sfp->except_text)) {
      if (StringStr (sfp->except_text, "alternative start codon") != NULL) {
        if (Loc_is_RefSeq (sfp->location)) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_AltStartCodon, "Unnecessary alternative start codon exception");
        }
      }
    }
  }

  part_loc = SeqLocPartialCheck (sfp->location);
  part_prod = SeqLocPartialCheckEx (sfp->product, farFetchProd);
  if ((part_loc & SLP_STOP) || (part_prod & SLP_STOP))
    no_end = TRUE;
  else {                        /* complete stop, so check for ragged end */

    len = SeqLocLen (sfp->location);
    if (crp->frame > 1)
      len -= (Int4) (crp->frame - 1);
    ragged = (Int2) (len % (Int4) (3));
    if (ragged) {
      len = SeqLocLen (sfp->location);
      cbp = crp->code_break;
      while (cbp != NULL) {
        pos1 = INT4_MAX;
        pos2 = -10;
        tmp = NULL;
        while ((tmp = SeqLocFindNext (cbp->loc, tmp)) != NULL) {
          pos = GetOffsetInLoc (tmp, sfp->location, SEQLOC_START);
          if (pos < pos1)
            pos1 = pos;
          pos = GetOffsetInLoc (tmp, sfp->location, SEQLOC_STOP);
          if (pos > pos2)
            pos2 = pos;
        }
        pos = pos2 - pos1;      /* codon length */
        if (pos >= 0 && pos <= 1 && pos2 == len - 1)
        {                       /*  a codon */
          /* allowing a partial codon at the end */
          ragged = 0;
        }

        cbp = cbp->next;
      }
    }
  }

  /* check for code break not on a codon */
  len = SeqLocLen (sfp->location);
  cbp = crp->code_break;
  while (cbp != NULL) {
    pos1 = INT4_MAX;
    pos2 = -10;
    tmp = NULL;
    while ((tmp = SeqLocFindNext (cbp->loc, tmp)) != NULL) {
      pos = GetOffsetInLoc (tmp, sfp->location, SEQLOC_START);
      if (pos < pos1)
        pos1 = pos;
      pos = GetOffsetInLoc (tmp, sfp->location, SEQLOC_STOP);
      if (pos > pos2)
        pos2 = pos;
    }
    pos = pos2 - pos1;          /* codon length */
    /* check for code break not on a codon */
    if (pos == 2 || (pos >= 0 && pos <= 1 && pos2 == len - 1)) {
      if (crp->frame == 2)
        pos = 1;
      else if (crp->frame == 3)
        pos = 2;
      else
        pos = 0;
      if ((pos1 % 3) != pos) {
        has_errors = TRUE;
        other_than_mismatch = TRUE;
        if (report_errors) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_TranslExceptPhase, "transl_except qual out of frame.");
        }
      }
    }


    cbp = cbp->next;
  }

  if (crp->frame > 1) {
    if (!(part_loc & SLP_START)) {
      sev = SEV_WARNING;
      if (Loc_is_RefSeq (sfp->location)) {
        sev = SEV_ERROR;
      }
      has_errors = TRUE;
      other_than_mismatch = TRUE;
      if (report_errors) {
        ValidErr (vsp, sev, ERR_SEQ_FEAT_SuspiciousFrame, "Suspicious CDS location - frame > 1 but not 5' partial");
      }
    } else if ((part_loc & SLP_NOSTART) && (!PartialAtSpliceSiteOrGap (vsp, sfp->location, SLP_NOSTART, &isgap, &badseq))) {
      if (PartialAtGapOrNs (vsp, sfp->location, SLP_NOSTART) && StringStr (sfp->comment, "coding region disrupted by sequencing gap") != NULL) {
        /* suppress */
      } else {
        sev = SEV_INFO;
        if (Loc_is_RefSeq (sfp->location)) {
          sev = SEV_ERROR;
        }
        has_errors = TRUE;
        other_than_mismatch = TRUE;
       if (report_errors) {
          ValidErr (vsp, sev, ERR_SEQ_FEAT_SuspiciousFrame, "Suspicious CDS location - frame > 1 and not at consensus splice site");
        }
      }
    }
  }

  if ((part_loc & SLP_START) || (part_prod & SLP_START))
    no_beg = TRUE;

  protseq = BSMerge (newprot, NULL);
  prot2len = StringLen (protseq);
  if (protseq != NULL) {
    len = prot2len;
    for (i = 0; i < len; i++) {
      residue1 = protseq [i];
      if (i == 0 && residue1 == '-') {
        got_dash = TRUE;
      }
      if (i == 0 && residue1 == 'X') {
        got_x = TRUE;
      }
      if (residue1 == '*') {
        if (i == (len - 1))
          got_stop = TRUE;
        else
          stop_count++;
      }
      if (residue1 == 'X') {
        x_count++;
      } else {
        nonx_count++;
      }
    }
    if (x_count > nonx_count) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_CDShasTooManyXs, "CDS translation consists of more than 50%s X residues", "%");
    }
  }

  if (annotated_by_transcript_or_proteomic) {
    if (1.2 * prot2len < prot1len) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_TransLen, "Protein product length [%ld] is more than 120%% of the %stranslation length [%ld]", prot1len, farstr, prot2len);
    }
  }

  /*
  prot2len = BSLen (newprot);
  len = prot2len;
  BSSeek (newprot, 0, SEEK_SET);
  for (i = 0; i < len; i++) {
    residue1 = BSGetByte (newprot);
    if ((i == 0) && (residue1 == '-'))
      got_dash = TRUE;
    if (residue1 == '*') {
      if (i == (len - 1))
        got_stop = TRUE;
      else
        stop_count++;
    }
  }
  */

  if (stop_count > 0) {
    if (got_dash) {
      has_errors = TRUE;
      other_than_mismatch = TRUE;
      sev = SEV_ERROR;
      if (unclassified_except) {
        sev = SEV_WARNING;
      }
      if (report_errors || unclassified_except) {
        ValidErr (vsp, sev, ERR_SEQ_FEAT_StartCodon,
                  "Illegal start codon (and %ld internal stops). Probably wrong genetic code [%d]", (long) stop_count, gccode);
        ValidErr (vsp, sev, ERR_SEQ_FEAT_InternalStop, "%ld internal stops (and illegal start codon). Genetic code [%d]", (long) stop_count, gccode);
      }
    } else if (got_x) {
      has_errors = TRUE;
      other_than_mismatch = TRUE;
      sev = SEV_ERROR;
      if (unclassified_except) {
        sev = SEV_WARNING;
      }
      if (report_errors || unclassified_except) {
        ValidErr (vsp, sev, ERR_SEQ_FEAT_StartCodon,
                    "Ambiguous start codon (and %ld internal stops). Possibly wrong genetic code [%d]", (long) stop_count, gccode);
        ValidErr (vsp, sev, ERR_SEQ_FEAT_InternalStop, "%ld internal stops (and ambiguous start codon). Genetic code [%d]", (long) stop_count, gccode);
      }
    } else {
      has_errors = TRUE;
      other_than_mismatch = TRUE;
      sev = SEV_ERROR;
      if (unclassified_except) {
        sev = SEV_WARNING;
      }
      if (report_errors || unclassified_except) {
        bsp = BioseqFindFromSeqLoc (sfp->location);
        if (bsp != NULL) {
          for (sip = bsp->id; sip != NULL; sip = sip->next) {
            switch (sip->choice) {
              case SEQID_GI :
                has_gi = TRUE;
                break;
              case SEQID_GENBANK :
              case SEQID_EMBL :
              case SEQID_DDBJ :
              case SEQID_TPG :
              case SEQID_TPE :
              case SEQID_TPD :
                is_ged = TRUE;
                break;
              case SEQID_OTHER :
                is_refseq = TRUE;
                break;
              default :
                break;
            }
          }
          if (has_gi && is_ged && (! is_refseq)) {
            sev = SEV_REJECT;
          }
        }
        ValidErr (vsp, sev, ERR_SEQ_FEAT_InternalStop, "%ld internal stops. Genetic code [%d]", (long) stop_count, gccode);
      }
    }
    prot_ok = FALSE;
    if (stop_count > 5)
      goto erret;
  } else if (got_dash) {
    has_errors = TRUE;
    other_than_mismatch = TRUE;
    if (report_errors) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_StartCodon, "Illegal start codon used. Wrong genetic code [%d] or protein should be partial", gccode);
    }
  } else if (got_x && (! partial5)) {
    has_errors = TRUE;
    other_than_mismatch = TRUE;
    if (report_errors) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_StartCodon, "Ambiguous start codon used. Wrong genetic code [%d] or protein should be partial", gccode);
    }
  }

  show_stop = TRUE;

  if (protid != NULL) {
    if (prot1seq == NULL && (! vsp->farFetchCDSproducts)) {
      goto erret;
    }
    if (prot1seq == NULL && sfp->product != NULL && vsp->farFetchCDSproducts) {
      SeqIdWrite (protid, id, PRINTID_FASTA_LONG, sizeof (id));
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_ProductFetchFailure, "Unable to fetch CDS product '%s'", id);
      goto erret;
    }
    if (prot1seq != NULL)
      prot1len = prot1seq->length;
  }

  if (prot1seq == NULL) {
    if (prot2len > 6) {
      if (! NGorNT (vsp->sep, sfp->location, &is_nc)) {
        sev = SEV_ERROR;
        if (DeltaOrFarSeg (vsp->sep, sfp->location)) {
          sev = SEV_WARNING;
        }
        if (is_nc) {
          sev = SEV_WARNING;
          sep = vsp->sep;
          if (sep != NULL && IS_Bioseq (sep)) {
            sev = SEV_NONE;
          }
        }
        if (sev != SEV_NONE) {
          has_errors = TRUE;
          other_than_mismatch = TRUE;
          if (report_errors) {
            ValidErr (vsp, sev, ERR_SEQ_FEAT_NoProtein, "No protein Bioseq given");
          }
        }
      }
    }
    goto erret;
  }

  len = prot2len;

  if ((got_stop) && (len == (prot1len + 1))) {  /* ok, got stop */
    len--;
  }

  if (! StreamCacheSetup (prot1seq, NULL, STREAM_EXPAND_GAPS, &sc)) {
    goto erret;
  }
  /*
  spp = SeqPortNew (prot1seq, 0, -1, 0, Seq_code_ncbieaa);
  if (spp == NULL)
    goto erret;
  */

  /* ignore terminal 'X' from partial last codon if present */

  done = FALSE;
  if ((!done) && (prot1len)) {
    /* prime the cache at a reasonable position near the end */
    if (prot1len > 4000) {
      StreamCacheSetPosition (&sc, prot1len - 2000);
    }
    residue1 = StreamCacheGetResidue (&sc);
  }
  while ((!done) && (prot1len)) {
    StreamCacheSetPosition (&sc, prot1len - 1);
    residue1 = StreamCacheGetResidue (&sc);
    /*
    SeqPortSeek (spp, (prot1len - 1), SEEK_SET);
    residue1 = SeqPortGetResidue (spp);
    */
    if (residue1 == 'X') {        /* remove terminal X */
      prot1len--;
      xcount1++;
    }
    else
      done = TRUE;
  }
  done = FALSE;
  while ((!done) && (len)) {
    /*
    BSSeek (newprot, (len - 1), SEEK_SET);
    residue2 = BSGetByte (newprot);
    */
    residue2 = protseq [len - 1];
    if (residue2 == 'X') {
      len--;
      xcount2++;
    }
    else
      done = TRUE;
  }

  if (xcount1 != xcount2) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_TerminalXDiscrepancy,
              "Terminal X count for CDS translation (%ld) and protein product sequence (%ld) are not equal",
              (long) xcount2, (long) xcount1);
  }

  if (len == prot1len) {        /* could be identical */
    StreamCacheSetPosition (&sc, 0);
    /*
    SeqPortSeek (spp, 0, SEEK_SET);
    BSSeek (newprot, 0, SEEK_SET);
    */
    for (i = 0; i < len; i++) {
      residue1 = protseq [i];
      residue2 = StreamCacheGetResidue (&sc);
      /*
      residue1 = BSGetByte (newprot);
      residue2 = SeqPortGetResidue (spp);
      */
      if (residue1 != residue2) {
        prot_ok = FALSE;
        if (residue2 == INVALID_RESIDUE)
          residue2 = '?';
        sev = SEV_ERROR;
        if (residue2 == 'X') {
          if (residue1 == 'B' || residue1 == 'Z' || residue1 == 'J') {
            sev = SEV_WARNING;
          }
        }
        if (i == 0) {
          if ((sfp->partial) && (!no_beg) && (!no_end)) { /* ok, it's partial */
            has_errors = TRUE;
            other_than_mismatch = TRUE;
            if (report_errors) {
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_PartialProblem, "Start of location should probably be partial");
            }
          } else if (residue1 == '-') {
            has_errors = TRUE;
            other_than_mismatch = TRUE;
            if (report_errors) {
              if (! got_dash) {
                ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_StartCodon, "Illegal start codon used. Wrong genetic code [%d] or protein should be partial", gccode);
              }
            }
          } else if (residue1 == 'X') {
            has_errors = TRUE;
            other_than_mismatch = TRUE;
            if (report_errors) {
              if (! got_x) {
                ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_StartCodon, "Ambiguous start codon used. Wrong genetic code [%d] or protein should be partial", gccode);
              }
            }
          } else {
            has_errors = TRUE;
            mismatches[mismatch].pos = i;
            mismatches[mismatch].cds_residue = residue1;
            mismatches[mismatch].prot_residue = residue2;
            mismatch++;
          }
        } else {
          has_errors = TRUE;
          if (mismatch >= 10) {
            mismatches[10].pos = i;
            mismatches[10].cds_residue = residue1;
            mismatches[10].prot_residue = residue2;
          } else {
            mismatches[mismatch].pos = i;
            mismatches[mismatch].cds_residue = residue1;
            mismatches[mismatch].prot_residue = residue2;
          }
          mismatch++;
        }
      }
    }

    if (report_errors && !mismatch_except) {
      if (mismatch > 10) {
        if (report_errors && !mismatch_except) {
          nuclocstr = MapToNTCoords (sfp, protid, mismatches[0].pos);
          loc2str = MapToNTCoords (sfp, protid, mismatches[10].pos);
          ValidErr (vsp, sev, ERR_SEQ_FEAT_MisMatchAA,
            "%d mismatches found.  First mismatch at %ld, residue in protein [%c] != translation [%c]%s%s.  Last mismatch at %ld, residue in protein [%c] != translation [%c]%s%s.  Genetic code [%d]",
            mismatch, 
            (long) (mismatches[0].pos + 1), mismatches[0].prot_residue, mismatches[0].cds_residue, 
            nuclocstr == NULL ? "" : " at ", nuclocstr == NULL ? "" : nuclocstr,
            (long) (mismatches[10].pos + 1), mismatches[10].prot_residue, mismatches[10].cds_residue, 
            loc2str == NULL ? "" : " at ", loc2str == NULL ? "" : loc2str,
            gccode);
          nuclocstr = MemFree (nuclocstr);
          loc2str = MemFree (loc2str);
        }
      } else {
        for (i = 0; i < mismatch; i++) {
          nuclocstr = MapToNTCoords (sfp, protid, mismatches[i].pos);
          ValidErr (vsp, sev, ERR_SEQ_FEAT_MisMatchAA,
                    "%sResidue %ld in protein [%c] != translation [%c]%s%s", farstr, 
                      (long) (mismatches[i].pos + 1), 
                      (char) mismatches[i].prot_residue, 
                      (char) mismatches[i].cds_residue,
                      nuclocstr == NULL ? "" : " at ",
                      nuclocstr == NULL ? "" : nuclocstr);
          nuclocstr = MemFree (nuclocstr);
        }
      }
    }

  } else {
    has_errors = TRUE;
    other_than_mismatch = TRUE;
    if (report_errors || (rna_editing && (prot1len < len - 1 || prot1len > len))) {
      ValidErr (vsp, rna_editing ? SEV_WARNING : trans_len_sev, ERR_SEQ_FEAT_TransLen, "Given protein length [%ld] does not match %stranslation length [%ld]", prot1len, farstr, len);
    }
  }

  if ((sfp->partial) && (!mismatch)) {
    if ((!no_beg) && (!no_end)) {       /* just didn't label */
      if (!got_stop) {
        has_errors = TRUE;
        other_than_mismatch = TRUE;
        if (report_errors) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_PartialProblem, "End of location should probably be partial");
        }
      } else {
        has_errors = TRUE;
        other_than_mismatch = TRUE;
        if (report_errors) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_PartialProblem, "This SeqFeat should not be partial");
        }
      }
      show_stop = FALSE;
    }
  }



erret:
  if (unlockProd) {
    BioseqUnlock (prot1seq);
  }

  if (show_stop) {
    if ((!got_stop) && (!no_end)) {
      has_errors = TRUE;
      other_than_mismatch = TRUE;
      if (report_errors) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_NoStop, "Missing stop codon");
      }
    } else if ((got_stop) && (no_end)) {
      has_errors = TRUE;
      other_than_mismatch = TRUE;
      if (report_errors) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_PartialProblem, "Got stop codon, but 3'end is labeled partial");
      }
    } else if ((got_stop) && (!no_end) && (ragged)) {
      has_errors = TRUE;
      other_than_mismatch = TRUE;
      sev = SEV_ERROR;
      if (unclassified_except) {
        sev = SEV_WARNING;
      }
      if (report_errors || unclassified_except) {
        ValidErr (vsp, sev, ERR_SEQ_FEAT_TransLen, "Coding region extends %d base(s) past stop codon", (int) ragged);
      }
    }
  }

  if (!prot_ok) {
    if (transl_except) {
      has_errors = TRUE;
      other_than_mismatch = TRUE;
      if (report_errors) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_TranslExcept, "Unparsed transl_except qual. Skipped");
      }
    }
  } else {
    if (transl_except) {
      has_errors = TRUE;
      other_than_mismatch = TRUE;
      if (report_errors) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_TranslExcept, "Unparsed transl_except qual (but protein is okay). Skipped");
      }
    }
  }

  if (prot2seq != NULL)
    BioseqFree (prot2seq);
  else
    BSFree (newprot);
  /*
  SeqPortFree (spp);
  */
  MemFree (protseq);
  ValNodeFree (codebreakhead);

  if (! report_errors) {
    if (! has_errors) {
      if ((! frameshift_except) && (! rearrange_except) && (! mixed_population) && (! low_quality) && (! artificial_location)) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UnnecessaryException, "CDS has exception but passes translation test");
      }
    } else if (unclassified_except && (! other_than_mismatch)) {
      if (mismatch * 50 <= len) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_ErroneousException,
                  "CDS has unclassified exception but only difference is %ld mismatches out of %ld residues",
                  (long) mismatch, (long) len);
      }
    } else if (product_replaced) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UnqualifiedException, "CDS has unqualified translated product replaced exception");
    }
  }
}


static void mRNAMatchesCompleteCDSEnd (SeqFeatPtr mrna, BoolPtr p5, BoolPtr p3)
{
  Boolean partial5, partial3;
  SeqFeatPtr cds;
  Uint2 strand;

  if (p5 != NULL) {
    *p5 = FALSE;
  } 
  if (p3 != NULL) {
    *p3 = FALSE;
  }

  cds = GetCDSformRNA (mrna);

  if (mrna == NULL || cds == NULL) {
    return;
  }

  strand = SeqLocStrand (mrna->location);

  CheckSeqLocForPartial (cds->location, &partial5, &partial3);
  if (p5 != NULL && !partial5) {
    if (strand == Seq_strand_minus) {
      if (SeqLocStop (cds->location) == SeqLocStop (mrna->location)) {
        *p5 = TRUE;
      }
    } else {
      if (SeqLocStart (cds->location) == SeqLocStart (mrna->location)) {
        *p5 = TRUE;
      }
    }
  }
   
  if (p3 != NULL && !partial3) {
    if (strand == Seq_strand_minus) {
      if (SeqLocStart (cds->location) == SeqLocStart (mrna->location)) {
        *p3 = TRUE;
      }
    } else {
      if (SeqLocStop (cds->location) == SeqLocStop (mrna->location)) {
        *p3 = TRUE;
      }
    }
  }
}


/*****************************************************************************
*
*   SpliceCheck(sfp)
*      checks for GT/AG rule at splice junctions
*
*****************************************************************************/
#define NOVALUE 0
#define HADGT 1
#define NOGT 2

static void SpliceCheckEx (ValidStructPtr vsp, SeqFeatPtr sfp, Boolean checkAll)
{
  SeqLocPtr       slp, nxt, head;
  Uint1           strand = Seq_strand_unknown;
  /*
  SeqPortPtr      spp = NULL;
  */
  SeqIdPtr        last_sip = NULL, sip;
  Int2            total, ctr;
  BioseqPtr       bsp = NULL;
  Int4            strt, stp, len = 0, donor, acceptor;
  Int2            residue1, residue2;
  Char            tbuf[40];
  Boolean         reportAsError, first, last, firstPartial, lastPartial, has_errors = FALSE,
                  report_errors = TRUE, checkExonDonor, checkExonAcceptor, pseudo;
  int             severity;
  Uint2           partialflag;
  SeqEntryPtr     sep;
  StreamCache     sc;
  SeqInt          sint;
  ValNode         vn;
  SeqMgrFeatContext  context;
  SeqFeatPtr      mrna, gene;
  GeneRefPtr      grp;
  Boolean         ignore_partial_mrna_5 = FALSE, ignore_partial_mrna_3 = FALSE;

  if (sfp == NULL)
    return;

  if (GetAppProperty ("NcbiSubutilValidation") != NULL)
    return;                     /* suppress if NCBISubValidate */

  /* suppress if organelle */
  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp != NULL && IsOrganelleBioseq(bsp)) {
    return;
  }

  /* specific biological exceptions suppress check */

  if (sfp->excpt) {
    if (StringISearch (sfp->except_text, "ribosomal slippage") != NULL||
        StringISearch (sfp->except_text, "artificial frameshift") != NULL ||
        StringISearch (sfp->except_text, "nonconsensus splice site") != NULL ||
        StringISearch (sfp->except_text, "adjusted for low-quality genome") != NULL ||
        StringISearch (sfp->except_text, "heterogeneous population sequenced") != NULL ||
        StringISearch (sfp->except_text, "low-quality sequence region") != NULL ||
        StringISearch (sfp->except_text, "artificial location") != NULL) {
      report_errors = FALSE;
    }
  }

  MemSet ((Pointer) &sint, 0, sizeof (SeqInt));
  MemSet ((Pointer) &vn, 0, sizeof (ValNode));

  head = sfp->location;
  if (head == NULL)
    return;

  if (LocationIsFar (sfp->location) && NoFetchFunctions ()) {
    vsp->far_fetch_failure = TRUE;
    return;
  }

  reportAsError = FALSE;
  if (GetAppProperty ("SpliceValidateAsError") != NULL) {
    reportAsError = TRUE;
  }

  slp = NULL;
  total = 0;
  while ((slp = SeqLocFindPart (head, slp, EQUIV_IS_ONE)) != NULL) {
    total++;
    if (slp->choice == SEQLOC_EQUIV)
      return;                   /* bail on this one */
    if (total == 1)
      strand = SeqLocStrand (slp);
    else {
      if (strand != SeqLocStrand (slp)) /* bail on mixed strand */
        return;
    }
  }

  if ((!checkAll) && total < 2)
    return;
  if (total < 1)
    return;

  slp = NULL;
  ctr = 0;

  first = TRUE;
  last = FALSE;
  firstPartial = FALSE;
  lastPartial = FALSE;

  if (sfp->idx.subtype == FEATDEF_mRNA) {
    mRNAMatchesCompleteCDSEnd (sfp, &ignore_partial_mrna_5, &ignore_partial_mrna_3);
  }


  /* genomic product set or NT_ contig always relaxes to SEV_WARNING */

  sep = vsp->sep;

  slp = SeqLocFindPart (head, slp, EQUIV_IS_ONE);
  while (slp != NULL) {
    nxt = SeqLocFindPart (head, slp, EQUIV_IS_ONE);
    last = (Boolean) (nxt == NULL);
    partialflag = SeqLocPartialCheck (slp);
    firstPartial = (Boolean) (first && (partialflag & SLP_START));
    lastPartial = (Boolean) (last && (partialflag & SLP_STOP));
    ctr++;
    sip = SeqLocId (slp);
    if (sip == NULL)
      break;

    bsp = BioseqFind (sip);

    if ((ctr == 1) || (!SeqIdMatch (sip, last_sip))) {
      /* spp = SeqPortFree (spp); */
      bsp = NULL;
      if (sip != NULL && (sip->choice != SEQID_GI || sip->data.intvalue > 0)) {
        bsp = BioseqLockById (sip);
      }
      if (bsp == NULL)
        break;
      len = bsp->length;
      if (strand != Seq_strand_minus) {
        if (! StreamCacheSetup (bsp, NULL, EXPAND_GAPS_TO_DASHES, &sc)) {
          BioseqUnlock (bsp);
          break;
        }
      } else {
        sint.from = 0;
        sint.to = len - 1;
        sint.strand = strand;
        sint.id = sip;
        vn.choice = SEQLOC_INT;
        vn.data.ptrvalue = (Pointer) &sint;
        vn.next = NULL;
        if (! StreamCacheSetup (NULL, &vn, EXPAND_GAPS_TO_DASHES, &sc)) {
          BioseqUnlock (bsp);
          break;
        }
      }
      /* spp = SeqPortNew (bsp, 0, -1, strand, Seq_code_ncbi4na); */
      BioseqUnlock (bsp);
      /*
      if (spp == NULL)
        break;
      */
      last_sip = sip;
    }

    acceptor = SeqLocStart (slp);
    donor = SeqLocStop (slp);

    if (acceptor < 0 || acceptor >= len || donor < 0 || donor >= len) {
      /*
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_Range,
                "Unable to check splice consensus because feature outside range of sequence");
      */
      return;
    }

    if (strand != Seq_strand_minus) {
      strt = acceptor;
      stp = donor;
    } else {
      strt = donor;
      donor = acceptor;
      acceptor = strt;
      stp = len - donor - 1;    /* orient to reverse complement seqport */
      strt = len - acceptor - 1;
    }

    checkExonDonor = FALSE;
    checkExonAcceptor = FALSE;
    if (checkAll) {
      pseudo = FALSE;
      grp = SeqMgrGetGeneXref (sfp);
      if (grp == NULL) {
        gene = SeqMgrGetOverlappingGene (sfp->location, &context);
        if (gene != NULL) {
          pseudo = gene->pseudo;
        }
      }
      if (! pseudo) {
        checkExonDonor = TRUE;
        checkExonAcceptor = TRUE;
        mrna = SeqMgrGetOverlappingmRNA (sfp->location, &context);
        if (mrna != NULL /* && (! mrna->partial) */ ) {
          if (strand != Seq_strand_minus) {
            if (donor == SeqLocStop (mrna->location) && (! context.partialR)) {
              checkExonDonor = FALSE;
            }
            if (acceptor == SeqLocStart (mrna->location) && (! context.partialL)) {
              checkExonAcceptor = FALSE;
            }
          } else {
            if (donor == SeqLocStart (mrna->location) && (! context.partialR)) {
              checkExonDonor = FALSE;
            }
            if (acceptor == SeqLocStop (mrna->location) && (! context.partialL)) {
              checkExonAcceptor = FALSE;
            }
          }
        }
      }
    }

    if (((checkExonDonor && (!lastPartial)) 
         || ctr < total 
         || (ctr == total && lastPartial && (sfp->idx.subtype != FEATDEF_mRNA || !ignore_partial_mrna_3) && sfp->idx.subtype != FEATDEF_exon)) 
        && (stp < (len - 2))) 
    {   /* check donor on all but last exon and on sequence */
      tbuf[0] = '\0';
      StreamCacheSetPosition (&sc, stp + 1);
      residue1 = StreamCacheGetResidue (&sc);
      residue2 = StreamCacheGetResidue (&sc);
      /*
      SeqPortSeek (spp, (stp + 1), SEEK_SET);
      residue1 = SeqPortGetResidue (spp);
      residue2 = SeqPortGetResidue (spp);
      */
      if (residue1 == '-' && residue2 == '-') {
        /* ignore gap, and suppress UnnecessaryException message */
        has_errors = TRUE;
      } else if (IS_residue (residue1) && IS_residue (residue2)) {
        if (residue1 != 'G' || residue2 != 'T') {        /* not T */
          if (residue1 == 'G' && residue2 == 'C') {       /* GC minor splice site */
            tbuf[0] = '\0';
            if (bsp == NULL) {
              StringCpy (tbuf, "?");
            } else if (vsp->suppressContext || vsp->convertGiToAccn) {
              WorstBioseqLabel (bsp, tbuf, 39, OM_LABEL_CONTENT);
            } else {
              BioseqLabel (bsp, tbuf, 39, OM_LABEL_CONTENT);
            }
            tbuf[39] = '\0';
            if (RareConsensusNotExpected (sfp)) {
              has_errors = TRUE;
              if (report_errors) {
                ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_RareSpliceConsensusDonor,
                          "Rare splice donor consensus (GC) found instead of (GT) after exon ending at position %ld of %s", (long) (donor + 1), tbuf);
              }
            }
          } else {
            if (checkExonDonor) {
              severity = SEV_WARNING;
            } else if (reportAsError) {
              severity = SEV_ERROR;
            } else {
              severity = SEV_WARNING;
            }
            tbuf[0] = '\0';
            if (bsp == NULL) {
              StringCpy (tbuf, "?");
            } else if (vsp->suppressContext || vsp->convertGiToAccn) {
              WorstBioseqLabel (bsp, tbuf, 39, OM_LABEL_CONTENT);
            } else {
              BioseqLabel (bsp, tbuf, 39, OM_LABEL_CONTENT);
            }
            tbuf[39] = '\0';
            has_errors = TRUE;
            if (report_errors) {
              ValidErr (vsp, severity, ERR_SEQ_FEAT_NotSpliceConsensusDonor,
                        "Splice donor consensus (GT) not found after exon ending at position %ld of %s", (long) (donor + 1), tbuf);
            }
          }
        }
      } else {
        has_errors = TRUE;
        if (report_errors) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_NotSpliceConsensusDonor,
                    "Bad sequence at splice donor after exon ending at position %ld of %s", (long) (donor + 1), tbuf);
        }
      }
    }

    if (((checkExonAcceptor && (!firstPartial))
         || ctr != 1 
         || (ctr == 1 && firstPartial && (sfp->idx.subtype != FEATDEF_mRNA || !ignore_partial_mrna_5) && sfp->idx.subtype != FEATDEF_exon)) 
        && (strt > 1)) 
    {
      StreamCacheSetPosition (&sc, strt - 2);
      residue1 = StreamCacheGetResidue (&sc);
      residue2 = StreamCacheGetResidue (&sc);
      /*
      SeqPortSeek (spp, (strt - 2), SEEK_SET);
      residue1 = SeqPortGetResidue (spp);
      residue2 = SeqPortGetResidue (spp);
      */
      if (residue1 == '-' && residue2 == '-') {
        /* ignore gap, and suppress UnnecessaryException message */
        has_errors = TRUE;
      } else if (IS_residue (residue1) && IS_residue (residue2)) {
        if (residue1 != 'A' || residue2 != 'G') {
          if (checkExonAcceptor) {
            severity = SEV_WARNING;
          } else if (reportAsError) {
            severity = SEV_ERROR;
          } else {
            severity = SEV_WARNING;
          }
          tbuf[0] = '\0';
          if (bsp == NULL) {
            StringCpy (tbuf, "?");
            SeqIdWrite (sip, tbuf, PRINTID_FASTA_SHORT, 39);
          } else if (vsp->suppressContext || vsp->convertGiToAccn) {
            WorstBioseqLabel (bsp, tbuf, 39, OM_LABEL_CONTENT);
          } else {
            BioseqLabel (bsp, tbuf, 39, OM_LABEL_CONTENT);
          }
          tbuf[39] = '\0';
          has_errors = TRUE;
          if (report_errors) {
            ValidErr (vsp, severity, ERR_SEQ_FEAT_NotSpliceConsensusAcceptor,
                      "Splice acceptor consensus (AG) not found before exon starting at position %ld of %s", (long) (acceptor + 1), tbuf);
          }
        }
      } else {
        has_errors = TRUE;
        if (report_errors) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_NotSpliceConsensusAcceptor,
                    "Bad sequence at splice acceptor before exon starting at position %ld of %s", (long) (acceptor + 1), tbuf);
        }
      }
    }

    first = FALSE;
    slp = nxt;
  }

  /* SeqPortFree (spp); */

  if (! report_errors) {
    if (! has_errors) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UnnecessaryException, "feature has exception but passes splice site test");
    }
  }
}

NLM_EXTERN void SpliceCheck (ValidStructPtr vsp, SeqFeatPtr sfp)
{
  SpliceCheckEx (vsp, sfp, FALSE);
}

/*****************************************************************************
*
*   CdsProductIdCheck (vsp, sfp)
*      code taken from asn2gnbk.c - release mode expects CDS product Bioseqs
*
*****************************************************************************/
static void CdsProductIdCheck (ValidStructPtr vsp, SeqFeatPtr sfp)

{
  SeqFeatPtr         gene;
  GeneRefPtr         grp;
  Boolean            juststop = FALSE;
  Boolean            okay = FALSE;
  SeqEntryPtr        oldscope;
  Boolean            partial5;
  Boolean            partial3;
  Boolean            pseudo = FALSE;
  SeqEntryPtr        sep;

   /* non-pseudo CDS must have /product */
   if (sfp->pseudo) {
     pseudo = TRUE;
   }
   grp = SeqMgrGetGeneXref (sfp);
   if (grp == NULL) {
     sep = GetTopSeqEntryForEntityID (sfp->idx.entityID);
     oldscope = SeqEntrySetScope (sep);
     gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
     SeqEntrySetScope (oldscope);
     if (gene != NULL) {
       grp = (GeneRefPtr) gene->data.value.ptrvalue;
       if (gene->pseudo) {
         pseudo = TRUE;
       }
     }
   }
   if (grp != NULL && grp->pseudo) {
     pseudo = TRUE;
   }
   if (sfp->location != NULL) {
     if (CheckSeqLocForPartial (sfp->location, &partial5, &partial3)) {
       if (partial5 && (! partial3)) {
         if (SeqLocLen (sfp->location) <= 5) {
           juststop = TRUE;
         }
       }
     }
   }
   if (pseudo || juststop) {
     okay = TRUE;
   } else if (sfp->product != NULL) {
     okay = TRUE;
   } else {
     if (sfp->excpt && (! StringHasNoText (sfp->except_text))) {
       if (StringStr (sfp->except_text, "rearrangement required for product") != NULL) {
         okay = TRUE;
       }
     }
   }
   if (! okay) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_MissingCDSproduct, "Expected CDS product absent");
   }
}

/*****************************************************************************
*
*   ValidateSeqLoc(vsp, slp, prefix)
*
*****************************************************************************/

static Int2 SeqLocMixCount (SeqLocPtr slp)

{
  Int2       count = 0;
  SeqLocPtr  loc;

  if (slp == NULL) return 0;

  while (slp != NULL) {
    if (slp->choice == SEQLOC_MIX) {
      count++;
      loc = (SeqLocPtr) slp->data.ptrvalue;
      count += SeqLocMixCount (loc);
    }
    slp = slp->next;
  }

  return count;
}

NLM_EXTERN void ValidateSeqLoc (ValidStructPtr vsp, SeqLocPtr slp, CharPtr prefix)
{
  SeqLocPtr       tmp, prev;
  Boolean         retval = TRUE, tmpval, mixed_strand = FALSE, unmarked_strand = FALSE,
                  ordered = TRUE, adjacent = FALSE, circular = FALSE, exception = FALSE,
                  bad = FALSE;
  CharPtr         ctmp;
  Uint1           strand2 = 0, strand1;
  ErrSev          sev, oldsev;
  SeqIntPtr       sip1, sip2, prevsip;
  SeqBondPtr      sbp;
  SeqPntPtr       spp;
  PackSeqPntPtr   pspp;
  SeqIdPtr        id1 = NULL, id2 = NULL;
  BioseqPtr       bsp;
  SeqFeatPtr      sfp = NULL;
  Int2            zeroGi = 0;
  Char            buf [32];
  SeqIdPtr        sip;

  if (slp == NULL)
    return;

  sfp = vsp->sfp;

  tmp = NULL;
  while ((tmp = SeqLocFindNext (slp, tmp)) != NULL) {
    sip = SeqLocId (tmp);
    if (sip != NULL && sip->choice == SEQID_GI && sip->data.intvalue <= 0) {
      zeroGi++;
    }
  }
  if (zeroGi > 0) {
    StringCpy (buf, "?");
    bsp = vsp->bsp;
    if (bsp != NULL) {
      SeqIdWrite (bsp->id, buf, PRINTID_FASTA_LONG, sizeof (buf) - 1);
    }
    if (zeroGi > 1) {
      ValidErr (vsp, SEV_REJECT, ERR_SEQ_FEAT_FeatureLocationIsGi0, "Feature has %d gi|0 locations on Bioseq %s",
                (int) zeroGi, buf);
    } else if (zeroGi > 0) {
      ValidErr (vsp, SEV_REJECT, ERR_SEQ_FEAT_FeatureLocationIsGi0, "Feature has %d gi|0 location on Bioseq %s",
                (int) zeroGi, buf);
    }
  }

  bsp = BioseqFindFromSeqLoc (slp);
  if (bsp != NULL && bsp->topology == 2) {
    circular = TRUE;
  }

  if (SeqLocMixCount (slp) > 1) {
      retval = FALSE;
      ctmp = SeqLocPrint (slp);
      if (ctmp != NULL && StringLen (ctmp) > 800) {
        StringCpy (ctmp + 797, "...");
      }
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_NestedSeqLocMix, "%s: SeqLoc [%s] has nested SEQLOC_MIX elements", prefix, ctmp);
      MemFree (ctmp);
  }

  tmp = NULL;
  prev = NULL;
  sip1 = NULL;
  prevsip = NULL;
  strand1 = Seq_strand_other;
  while ((tmp = SeqLocFindNext (slp, tmp)) != NULL) {
    tmpval = TRUE;
    switch (tmp->choice) {
    case SEQLOC_INT:
      sip1 = prevsip;
      sip2 = (SeqIntPtr) (tmp->data.ptrvalue);
      strand2 = sip2->strand;
      id2 = sip2->id;
      tmpval = SeqIntCheck (sip2);
      if ((tmpval) && (sip1 != NULL)) {
        if (SeqIdForSameBioseq (sip1->id, sip2->id)) {
          if (strand2 == Seq_strand_minus) {
            if (sip1->to < sip2->to && (! circular)) {
              ordered = FALSE;
            }
            if (sip2->to + 1 == sip1->from) {
              adjacent = TRUE;
            }
          } else {
            if (sip1->to > sip2->to && (! circular)) {
              ordered = FALSE;
            }
            if (sip1->to + 1 == sip2->from) {
              adjacent = TRUE;
            }
          }
        }
      }
      if (prevsip != NULL) {
        if (SeqIdForSameBioseq (prevsip->id, sip2->id)) {
          if (prevsip->strand == sip2->strand && prevsip->from == sip2->from && prevsip->to == sip2->to) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_DuplicateInterval, "Duplicate exons in location");
          }
        }
      }
      prevsip = sip2;
      break;
    case SEQLOC_PNT:
      spp = (SeqPntPtr) (tmp->data.ptrvalue);
      strand2 = spp->strand;
      id2 = spp->id;
      tmpval = SeqPntCheck (spp);
      prevsip = NULL;
      break;
    case SEQLOC_PACKED_PNT:
      pspp = (PackSeqPntPtr) (tmp->data.ptrvalue);
      strand2 = pspp->strand;
      id2 = pspp->id;
      tmpval = PackSeqPntCheck (pspp);
      prevsip = NULL;
      break;
    case SEQLOC_BOND:
      sbp = (SeqBondPtr) tmp->data.ptrvalue;
      if (sbp != NULL) {
        spp = (SeqPntPtr) sbp->a;
        if (spp != NULL) {
          tmpval = SeqPntCheck (spp);
        }
        /* if already failed, no need to check second point */
        if (tmpval) {
          spp = (SeqPntPtr) sbp->b;
          if (spp != NULL) {
            tmpval = SeqPntCheck (spp);
          }
        }
      }
    case SEQLOC_NULL:
      break;
    default:
      strand2 = Seq_strand_other;
      id2 = NULL;
      prevsip = NULL;
      break;
    }
    if (!tmpval) {
      retval = FALSE;
      ctmp = SeqLocPrint (tmp);
      if (ctmp != NULL && StringLen (ctmp) > 800) {
        StringCpy (ctmp + 797, "...");
      }
      ValidErr (vsp, SEV_REJECT, ERR_SEQ_FEAT_Range, "%s: SeqLoc [%s] out of range", prefix, ctmp);
      MemFree (ctmp);

    }

    if (tmp->choice != SEQLOC_NULL) {
      if ((strand1 != Seq_strand_other) && (strand2 != Seq_strand_other)) {
        if (SeqIdForSameBioseq (id1, id2)) {
          if (strand1 != strand2) {
            if (strand1 == Seq_strand_plus && strand2 == Seq_strand_unknown) {
              unmarked_strand = TRUE;
            } else if (strand1 == Seq_strand_unknown && strand2 == Seq_strand_plus) {
              unmarked_strand = TRUE;
            } else {
              mixed_strand = TRUE;
            }
          }
        }
      }

      strand1 = strand2;
      id1 = id2;
    }
  }

  if (sfp != NULL) {

    /* Publication intervals ordering does not matter */

    if (sfp->idx.subtype == FEATDEF_PUB) {
      ordered = TRUE;
      adjacent = FALSE;
    }

    /* ignore ordering of heterogen bonds */

    if (sfp->data.choice == SEQFEAT_HET) {
      ordered = TRUE;
      adjacent = FALSE;
    }

    /* misc_recomb intervals SHOULD be in reverse order */

    if (sfp->idx.subtype == FEATDEF_misc_recomb) {
      ordered = TRUE;
    }

    /* primer_bind intervals MAY be in on opposite strands */

    if (sfp->idx.subtype == FEATDEF_primer_bind) {
      mixed_strand = FALSE;
      unmarked_strand = FALSE;
      ordered = TRUE;
    }

    if (sfp->excpt) {
      exception = TRUE;
    }
  }

  if (adjacent) {
    ctmp = SeqLocPrint (slp);
    if (exception) {
      sev = SEV_WARNING;
    } else {
      sev = SEV_ERROR;
    }
    if (ctmp != NULL && StringLen (ctmp) > 800) {
      StringCpy (ctmp + 797, "...");
    }
    ValidErr (vsp, sev, ERR_SEQ_FEAT_AbuttingIntervals, "%s: Adjacent intervals in SeqLoc [%s]", prefix, ctmp);
    MemFree (ctmp);
  }

  if (exception) {
    /* trans splicing exception turns off both mixed_strand and out_of_order messages */
    if (StringISearch (sfp->except_text, "trans-splicing") != NULL) {
      return;
    }
  }

  if (mixed_strand || unmarked_strand || (!ordered)) {
    ctmp = SeqLocPrint (slp);
    if (ctmp != NULL && StringLen (ctmp) > 800) {
      StringCpy (ctmp + 797, "...");
    }
    sev = SEV_ERROR;
    if (vsp->is_small_genome_set) {
      sev = SEV_WARNING;
    }
    if (mixed_strand) {
      if (vsp->is_small_genome_set) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_MixedStrand, "%s: Mixed strands in SeqLoc [%s] in small genome set - set trans-splicing exception if appropriate", prefix, ctmp);
      } else {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_MixedStrand, "%s: Mixed strands in SeqLoc [%s]", prefix, ctmp);
      }
    } else if (unmarked_strand) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_MixedStrand, "%s: Mixed plus and unknown strands in SeqLoc [%s]", prefix, ctmp);
    }
    if (!ordered)
      ValidErr (vsp, sev, ERR_SEQ_FEAT_SeqLocOrder, "%s: Intervals out of order in SeqLoc [%s]", prefix, ctmp);
    MemFree (ctmp);
    return;
  }

  if (sfp != NULL) {

    /* ignore special case features here as well */

    if (sfp->idx.subtype == FEATDEF_PUB ||
        sfp->data.choice == SEQFEAT_HET ||
        sfp->idx.subtype == FEATDEF_misc_recomb ||
        sfp->idx.subtype == FEATDEF_primer_bind)
      return;
  }

  /* newer check for intervals out of order on segmented bioseq */

  if (bsp == NULL || bsp->repr != Seq_repr_seg) return;

  oldsev = ErrSetMessageLevel (SEV_ERROR);
  bad = SeqLocBadSortOrder (bsp, slp);
  ErrSetMessageLevel (oldsev);
  if (bad) {
    ctmp = SeqLocPrint (slp);
    if (ctmp != NULL && StringLen (ctmp) > 800) {
      StringCpy (ctmp + 797, "...");
    }
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_SeqLocOrder, "%s: Intervals out of order in SeqLoc [%s]", prefix, ctmp);
    MemFree (ctmp);
  }

  /* newer check for mixed strand on segmented bioseq */

  oldsev = ErrSetMessageLevel (SEV_ERROR);
  bad = SeqLocMixedStrands (bsp, slp);
  ErrSetMessageLevel (oldsev);
  if (bad) {
    ctmp = SeqLocPrint (slp);
    if (ctmp != NULL && StringLen (ctmp) > 800) {
      StringCpy (ctmp + 797, "...");
    }
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_MixedStrand, "%s: Mixed strands in SeqLoc [%s]", prefix, ctmp);
    MemFree (ctmp);
  }
}

/*****************************************************************************
*
*   SeqGraph validation section
*
*****************************************************************************/

typedef struct gphgetdata
{
  ValNodePtr      vnp;
  BioseqPtr       bsp;
}
GphGetData     , PNTR GphGetPtr;

typedef struct grphitem
{
  SeqGraphPtr     sgp;
  Int4            left;
  Int4            right;
  Int2            index;
}
GrphItem       , PNTR GrphItemPtr;

static void GetGraphsProc (SeqGraphPtr sgp, Pointer userdata)
{
  GphGetPtr       ggp;
  GrphItemPtr     gip;

  ggp = (GphGetPtr) userdata;
  if (ggp == NULL || sgp == NULL) return;
  /* only phrap or gap4 currently allowed */
  if (StringICmp (sgp->title, "Phrap Quality") == 0 || StringICmp (sgp->title, "Phred Quality") == 0 || StringICmp (sgp->title, "Gap4") == 0) {
    /* data type must be bytes */
    if (sgp->flags[2] == 3) {
      if (SeqIdIn (SeqLocId (sgp->loc), ggp->bsp->id)) {
        gip = (GrphItemPtr) MemNew (sizeof (GrphItem));
        if (gip == NULL) return;
        gip->sgp = sgp;
        gip->left = GetOffsetInBioseq (sgp->loc, ggp->bsp, SEQLOC_LEFT_END);
        gip->right = GetOffsetInBioseq (sgp->loc, ggp->bsp, SEQLOC_RIGHT_END);
        ValNodeAddPointer (&(ggp->vnp), 0, (Pointer) gip);
      }
    }
  }
  return;
}

static int LIBCALLBACK SortSeqGraphProc (VoidPtr ptr1, VoidPtr ptr2)
{
  GrphItemPtr     gip1, gip2;
  ValNodePtr      vnp1, vnp2;

  if (ptr1 == NULL || ptr2 == NULL)
    return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL)
    return 0;
  gip1 = (GrphItemPtr) vnp1->data.ptrvalue;
  gip2 = (GrphItemPtr) vnp2->data.ptrvalue;
  if (gip1 == NULL || gip2 == NULL)
    return 0;
  if (gip1->left > gip2->left) {
    return 1;
  } else if (gip1->left < gip2->left) {
    return -1;
  } else if (gip1->right > gip2->right) {
    return -1;
  } else if (gip2->right < gip2->right) {
    return 1;
  }
  return 0;
}

/* gets valnode list of sorted graphs in GrphItem structures */

static ValNodePtr GetSeqGraphsOnBioseq (Uint2 entityID, BioseqPtr bsp)
{
  GphGetData   ggd;
  GrphItemPtr  gip;
  Int2         index;
  ValNodePtr   vnp;

  ggd.vnp = NULL;
  ggd.bsp = bsp;
  VisitGraphsOnBsp (bsp, (Pointer) &ggd, GetGraphsProc);
  for (vnp = ggd.vnp, index = 1; vnp != NULL; vnp = vnp->next, index++) {
    gip = (GrphItemPtr) vnp->data.ptrvalue;
    if (gip != NULL) {
      gip->index = index;
    }
  }
  ggd.vnp = ValNodeSort (ggd.vnp, SortSeqGraphProc);
  return ggd.vnp;
}

static Boolean NextLitLength (DeltaSeqPtr next, Int4Ptr lenp)

{
  SeqLitPtr  slp;

  if (lenp == NULL) return FALSE;
  *lenp = 0;
  if (next == NULL || next->choice != 2) return FALSE;
  slp = (SeqLitPtr) next->data.ptrvalue;
  if (slp == NULL || slp->seq_data == NULL) return FALSE;
  *lenp = slp->length;
  return TRUE;
}

static void ValidateGraphsOnBioseq (GatherContextPtr gcp)
{
  Byte            scores [400];
  ByteStorePtr    bs;
  BioseqPtr       bsp;
  Int2            k, val, index, scount;
  Int4            curroffset = 0, gphlen = 0, seqlen = 0, slplen,
                  bslen, min = INT4_MAX, max = INT4_MIN, j, lastloc = -1,
                  numBases, NsWithScore, GapsWithScore, ACGTsWithoutScore,
                  ambigWithoutScore, valsBelowMin, valsAboveMax,
                  firstN, firstACGT, firstAmbig, pos, litlen, nxtlen;
  FloatHi         pct;
  DeltaSeqPtr     dsp, next;
  Uint2           entityID, olditemtype = 0, numdsp = 0, numsgp = 0;
  Uint4           firstsgitemid = 0;
  Uint4           olditemid = 0;
  GrphItemPtr     gip;
  ValNodePtr      head, vnp;
  Boolean         outOfOrder = FALSE, fa2htgsBug = FALSE, overlaps = FALSE;
  Uint1           residue;
  SeqGraphPtr     sgp;
  SeqIntPtr       sintp;
  SeqLocPtr       slocp;
  SeqLitPtr       slp;
  StreamCache     sc;
  ValidStructPtr  vsp;
  Boolean         single_report_mode = TRUE;
  CharPtr         ctmp;

  vsp = (ValidStructPtr) gcp->userdata;
  bsp = (BioseqPtr) gcp->thisitem;
  if (vsp == NULL || bsp == NULL)
    return;
  if (!ISA_na (bsp->mol))
    return;

  vsp->bsp = bsp;
  vsp->descr = NULL;
  vsp->sfp = NULL;
  vsp->bssp = (BioseqSetPtr) gcp->parentitem;

  if (SeqMgrGetParentOfPart (bsp, NULL) != NULL)
    return;

  entityID = ObjMgrGetEntityIDForPointer (bsp);
  head = GetSeqGraphsOnBioseq (entityID, bsp);
  if (head == NULL)
    return;

  olditemid = gcp->itemID;
  olditemtype = gcp->thistype;
  gcp->thistype = OBJ_SEQGRAPH;

  for (vnp = head, index = 1; vnp != NULL; vnp = vnp->next, index++) {
    gip = (GrphItemPtr) vnp->data.ptrvalue;
    if (gip == NULL)
      continue;

    sgp = gip->sgp;
    if (sgp == NULL)
      continue;
    gcp->itemID = sgp->idx.itemID;
    if (firstsgitemid == 0) {
      firstsgitemid = sgp->idx.itemID;
    }

    if (gip->index != index) {
      outOfOrder = TRUE;
      if (gip->index == 129 && index == 2) {
        fa2htgsBug = TRUE;
      }
    }
    if (gip->left <= lastloc) {
      overlaps = TRUE;
    }
    lastloc = gip->right;
    min = MIN ((Int4) min, (Int4) sgp->min.intvalue);
    max = MAX ((Int4) max, (Int4) sgp->max.intvalue);

    if (sgp->min.intvalue < 0 || sgp->min.intvalue > 100) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_GRAPH_GraphMin, "Graph min (%ld) out of range", (long) sgp->min.intvalue);
    }

    if (sgp->max.intvalue > 100) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_GRAPH_GraphMax, "Graph max (%ld) out of range", (long) sgp->max.intvalue);
    }
    if (sgp->max.intvalue <= 0) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_GRAPH_GraphMax, "Graph max (%ld) out of range", (long) sgp->max.intvalue);
    }

    gphlen += sgp->numval;
    bs = (ByteStorePtr) sgp->values;
    if (bs != NULL) {
      bslen = BSLen (bs);
      if (sgp->numval != bslen) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_GRAPH_GraphByteLen, "SeqGraph (%ld) and ByteStore (%ld) length mismatch", (long) sgp->numval, (long) bslen);
      }
    }
  }
  if (outOfOrder) {
    gcp->itemID = firstsgitemid;
    if (fa2htgsBug) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_GRAPH_GraphOutOfOrder, "Graph components are out of order - probably caused by old fa2htgs bug");
    } else {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_GRAPH_GraphOutOfOrder, "Graph components are out of order - may be a software bug");
    }
  }
  if (overlaps) {
    gcp->itemID = firstsgitemid;
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_GRAPH_GraphOverlap, "Graph components overlap, with multiple scores for a single base");
  }

  if (bsp->repr == Seq_repr_raw) {
    seqlen = bsp->length;
  } else if (bsp->repr == Seq_repr_delta) {
    for (dsp = (DeltaSeqPtr) (bsp->seq_ext); dsp != NULL; dsp = dsp->next) {
      switch (dsp->choice) {
      case 1:
        slocp = (SeqLocPtr) dsp->data.ptrvalue;
        if (slocp == NULL)
          break;
        if (slocp->choice != SEQLOC_NULL) {
          seqlen += SeqLocLen (slocp);
        }
        break;
      case 2:
        slp = (SeqLitPtr) dsp->data.ptrvalue;
        if (slp == NULL || slp->seq_data == NULL)
          break;
        seqlen += slp->length;
        break;
      default:
        break;
      }
    }
  }

  if (seqlen != gphlen && bsp->length != gphlen) {
    gcp->itemID = firstsgitemid;
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_GRAPH_GraphBioseqLen, "SeqGraph (%ld) and Bioseq (%ld) length mismatch", (long) gphlen, (long) seqlen);
  }

  if (bsp->repr == Seq_repr_delta) {
    if (head != NULL && head->next != NULL) {
      for (dsp = (DeltaSeqPtr) (bsp->seq_ext), vnp = head; dsp != NULL && vnp != NULL; dsp = next) {
        next = dsp->next;
        gip = (GrphItemPtr) vnp->data.ptrvalue;
        if (gip == NULL)
          continue;
        sgp = gip->sgp;
        if (sgp == NULL)
          continue;
        switch (dsp->choice) {
        case 1:
          slocp = (SeqLocPtr) dsp->data.ptrvalue;
          if (slocp != NULL && slocp->choice != SEQLOC_NULL) {
            slplen = SeqLocLen (slocp);
            curroffset += slplen;
            if (sgp->numval != slplen) {
              gcp->itemID = sgp->idx.itemID;
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_GRAPH_GraphSeqLocLen, "SeqGraph (%ld) and SeqLoc (%ld) length mismatch", (long) sgp->numval, (long) slplen);
            }
            numdsp++;
            if (vnp != NULL) {
              vnp = vnp->next;
              numsgp++;
            }
          }
          break;
        case 2:
          slp = (SeqLitPtr) dsp->data.ptrvalue;
          litlen = 0;
          if (slp != NULL) {
            litlen = slp->length;
          }
          if (slp != NULL && slp->seq_data != NULL) {
            while (NextLitLength (next, &nxtlen)) {
              litlen += nxtlen;
              next = next->next;
            }
            if (sgp->numval != litlen) {
              gcp->itemID = sgp->idx.itemID;
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_GRAPH_GraphSeqLitLen, "SeqGraph (%ld) and SeqLit (%ld) length mismatch",
                        (long) sgp->numval, (long) litlen);
            }
            slocp = sgp->loc;
            if (slocp != NULL && slocp->choice == SEQLOC_INT) {
              sintp = (SeqIntPtr) slocp->data.ptrvalue;
              if (sintp != NULL) {
                if (sintp->from != curroffset) {
                  gcp->itemID = sgp->idx.itemID;
                  ValidErr (vsp, SEV_ERROR, ERR_SEQ_GRAPH_GraphStartPhase, "SeqGraph (%ld) and SeqLit (%ld) start do not coincide",
                            (long) sintp->from, (long) curroffset);
                }
                if (sintp->to != litlen + curroffset - 1) {
                  gcp->itemID = sgp->idx.itemID;
                  ValidErr (vsp, SEV_ERROR, ERR_SEQ_GRAPH_GraphStopPhase, "SeqGraph (%ld) and SeqLit (%ld) stop do not coincide",
                            (long) sintp->to, (long) (litlen + curroffset - 1));
                }
              }
            }
            numdsp++;
            if (vnp != NULL) {
              vnp = vnp->next;
              numsgp++;
            }
          }
          if (slp != NULL) {
            curroffset += litlen;
          }
          break;
        default:
          break;
        }
      }
      for (dsp = (DeltaSeqPtr) (bsp->seq_ext), numdsp = 0; dsp != NULL; dsp = next) {
        next = dsp->next;
        switch (dsp->choice) {
        case 1:
          slocp = (SeqLocPtr) dsp->data.ptrvalue;
          if (slocp != NULL && slocp->choice != SEQLOC_NULL) {
            numdsp++;
          }
          break;
        case 2:
          slp = (SeqLitPtr) dsp->data.ptrvalue;
          if (slp != NULL && slp->seq_data != NULL) {
            while (NextLitLength (next, &nxtlen)) {
              next = next->next;
            }
            numdsp++;
          }
          break;
        default:
          break;
        }
      }
      for (vnp = head, numsgp = 0; vnp != NULL; vnp = vnp->next, numsgp++)
        continue;
      if (numdsp != numsgp) {
        gcp->itemID = firstsgitemid;
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_GRAPH_GraphDiffNumber, "Different number of SeqGraph (%d) and SeqLit (%d) components", (int) numsgp, (int) numdsp);
      }
    }
  }

  numBases = 0;
  NsWithScore = 0;
  GapsWithScore = 0;
  ACGTsWithoutScore = 0;
  ambigWithoutScore = 0;
  valsBelowMin = 0;
  valsAboveMax = 0;
  firstN = -1;
  firstACGT = -1;
  firstAmbig = -1;

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    gip = (GrphItemPtr) vnp->data.ptrvalue;
    if (gip == NULL)
      continue;
    sgp = gip->sgp;
    if (sgp == NULL)
      continue;

    if (sgp->loc == NULL ||
        SeqLocStart (sgp->loc) < 0 ||
        SeqLocStop (sgp->loc) >= bsp->length ||
        (sgp->loc->choice != SEQLOC_INT && sgp->loc->choice != SEQLOC_WHOLE) ||
        SeqLocStrand (sgp->loc) == Seq_strand_minus) {
      gcp->itemID = sgp->idx.itemID;
      ctmp = SeqLocPrint (sgp->loc);
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_GRAPH_GraphLocInvalid, "SeqGraph location (%s) is invalid", ctmp == NULL ? "Unknown" : ctmp);
      ctmp = MemFree (ctmp);
      continue;
    }

    if (! StreamCacheSetup (NULL, sgp->loc, EXPAND_GAPS_TO_DASHES, &sc)) continue;
    slplen = SeqLocLen (sgp->loc);

    bs = (ByteStorePtr) sgp->values;
    BSSeek (bs, 0, SEEK_SET);
    j = 0;
    val = 0;

    scount = (Int2) BSRead (bs, scores, sizeof (scores));
    k = 0;

    if (! single_report_mode) {
      numBases = 0;
      NsWithScore = 0;
      GapsWithScore = 0;
      ACGTsWithoutScore = 0;
      ambigWithoutScore = 0;
      valsBelowMin = 0;
      valsAboveMax = 0;
      firstN = -1;
      firstACGT = -1;
      firstAmbig = -1;
    }

    pos = gip->left;

    while ((residue = StreamCacheGetResidue (&sc)) != '\0' && j < sgp->numval) {
      if (IS_residue (residue)) {
        numBases++;
        /* val = (Int2) BSGetByte (bs); */
        if (k >= scount) {
          if (scount > 0) {
            scount = (Int2) BSRead (bs, scores, sizeof (scores));
          }
          k = 0;
        }
        if (scount > 0) {
          val = (Int2) scores [k];
          k++;
        } else {
          val = 0;
        }
        if (val < sgp->min.intvalue || val < 0) {
          valsBelowMin++;
        }
        if (val > sgp->max.intvalue || val > 100) {
          valsAboveMax++;
        }
        j++;
        switch (residue) {
        case '-': /* 0 */
          if (val > 0) {
            GapsWithScore++;
          }
          break;
        case 'A': /* 1, 2, 4, 8 */
        case 'C':
        case 'G':
        case 'T':
          if (val == 0) {
            ACGTsWithoutScore++;
            if (firstACGT == -1) {
              firstACGT = pos;
            }
          }
          break;
        case 'N': /* 15 */
          if (val > 0) {
            NsWithScore++;
            if (firstN == -1) {
              firstN = pos;
            }
          }
          break;
        default:
          if (val == 0) {
            ambigWithoutScore++;
            if (firstAmbig == -1) {
              firstAmbig = pos;
            }
          }
          break;
        }
      }
      pos++;
    }

    if (! single_report_mode) {
      gcp->itemID = sgp->idx.itemID;
      if (ACGTsWithoutScore > 0) {
        if (ACGTsWithoutScore * 10 >= numBases) {
          pct = (FloatHi) (ACGTsWithoutScore) * 100.0 / (FloatHi) numBases;
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_GRAPH_GraphACGTScoreMany, "%ld ACGT bases (%3.2f%s) have zero score value - first one at position %ld",
                    (long) ACGTsWithoutScore, (double) pct, "%", (long) (firstACGT + 1));
        } else {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_GRAPH_GraphACGTScore, "%ld ACGT bases have zero score value - first one at position %ld",
                    (long) ACGTsWithoutScore, (long) (firstACGT + 1));
        }
      }
      if (NsWithScore > 0) {
        if (NsWithScore * 10 >= numBases) {
          pct = (FloatHi) (NsWithScore) * 100.0 / (FloatHi) numBases;
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_GRAPH_GraphNScoreMany, "%ld N bases (%3.2f%s) have positive score value - first one at position %ld",
                    (long) NsWithScore, (double) pct, "%", (long) (firstN + 1));
        } else {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_GRAPH_GraphNScore, "%ld N bases have positive score value - first one at position %ld",
                    (long) NsWithScore, (long) (firstN + 1));
        }
      }
      if (GapsWithScore > 0) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_GRAPH_GraphGapScore, "%ld gap bases have positive score value", (long) GapsWithScore);
      }
      if (valsBelowMin > 0) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_GRAPH_GraphBelow, "%ld quality scores have values below the reported minimum or 0", (long) valsBelowMin);
      }
      if (valsAboveMax > 0) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_GRAPH_GraphAbove, "%ld quality scores have values above the reported maximum or 100", (long) valsAboveMax);
      }
    }
  }

  gcp->itemID = olditemid;
  gcp->thistype = olditemtype;

  if (single_report_mode) {
    if (ACGTsWithoutScore > 0) {
      if (ACGTsWithoutScore * 10 >= numBases) {
        pct = (FloatHi) (ACGTsWithoutScore) * 100.0 / (FloatHi) numBases;
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_GRAPH_GraphACGTScoreMany, "%ld ACGT bases (%3.2f%s) have zero score value - first one at position %ld",
                  (long) ACGTsWithoutScore, (double) pct, "%", (long) (firstACGT + 1));
      } else {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_GRAPH_GraphACGTScore, "%ld ACGT bases have zero score value - first one at position %ld",
                  (long) ACGTsWithoutScore, (long) (firstACGT + 1));
      }
    }
    if (NsWithScore > 0) {
      if (NsWithScore * 10 >= numBases) {
        pct = (FloatHi) (NsWithScore) * 100.0 / (FloatHi) numBases;
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_GRAPH_GraphNScoreMany, "%ld N bases (%3.2f%s) have positive score value - first one at position %ld",
                  (long) NsWithScore, (double) pct, "%", (long) (firstN + 1));
      } else {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_GRAPH_GraphNScore, "%ld N bases have positive score value - first one at position %ld",
                  (long) NsWithScore, (long) (firstN + 1));
      }
    }
    if (GapsWithScore > 0) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_GRAPH_GraphGapScore, "%ld gap bases have positive score value", (long) GapsWithScore);
    }
    if (valsBelowMin > 0) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_GRAPH_GraphBelow, "%ld quality scores have values below the reported minimum or 0", (long) valsBelowMin);
    }
    if (valsAboveMax > 0) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_GRAPH_GraphAbove, "%ld quality scores have values above the reported maximum or 100", (long) valsAboveMax);
    }
  }

  ValNodeFreeData (head);
}

/*****************************************************************************
*
*   PatchBadSequence(bsp)
*
*****************************************************************************/
NLM_EXTERN Boolean PatchBadSequence (BioseqPtr bsp)
{
  ByteStorePtr    newseq;
  SeqPortPtr      spp;
  Boolean         is_na;
  Uint1           seqcode;
  Int2            repchar, residue;
  Int4            i, len;

  if (bsp == NULL)
    return FALSE;
  if (!((bsp->repr == Seq_repr_raw) || (bsp->repr == Seq_repr_const)))
    return FALSE;

  is_na = ISA_na (bsp->mol);
  if (is_na) {
    seqcode = Seq_code_iupacna;
    repchar = (Int2) 'N';       /* N */
  } else {
    seqcode = Seq_code_iupacaa;
    repchar = (Int2) 'X';
  }

  spp = SeqPortNew (bsp, 0, -1, 0, seqcode);
  if (spp == NULL)
    return FALSE;

  len = bsp->length;
  newseq = BSNew (len);
  if (newseq == NULL) {
    SeqPortFree (spp);
    return FALSE;
  }

  for (i = 0; i < len; i++) {
    residue = SeqPortGetResidue (spp);
    if (residue == INVALID_RESIDUE) {
      residue = repchar;
    }
    BSPutByte (newseq, residue);
  }

  SeqPortFree (spp);
  SeqDataFree (bsp->seq_data, bsp->seq_data_type);
  bsp->seq_data = (SeqDataPtr) newseq;
  bsp->seq_data_type = seqcode;

  BioseqRawPack (bsp);

  return TRUE;
}

static void FindABioseq (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
  BioseqPtr PNTR  bp;
  BioseqPtr       bsp;

  bp = (BioseqPtr PNTR) data;
  if (*bp != NULL)              /* already got one */
    return;

  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) (sep->data.ptrvalue);
    *bp = bsp;
  }
  return;
}

NLM_EXTERN CharPtr FindIDForEntry (SeqEntryPtr sep, CharPtr buf)
{
  BioseqPtr       bsp = NULL;

  if ((sep == NULL) || (buf == NULL))
    return NULL;

  *buf = '\0';
  SeqEntryExplore (sep, (Pointer) (&bsp), FindABioseq);

  if (bsp == NULL)
    return NULL;

  SeqIdPrint (bsp->id, buf, PRINTID_FASTA_LONG);
  return buf;
}

static CharPtr TrimSpacesOnEitherSide (CharPtr str)
{
  Uchar           ch;
  CharPtr         dst;
  CharPtr         ptr;

  if (str != NULL && str[0] != '\0') {
    dst = str;
    ptr = str;
    ch = *ptr;
    while (ch != '\0' && ch <= ' ') {
      ptr++;
      ch = *ptr;
    }
    while (ch != '\0') {
      *dst = ch;
      dst++;
      ptr++;
      ch = *ptr;
    }
    *dst = '\0';
    dst = NULL;
    ptr = str;
    ch = *ptr;
    while (ch != '\0') {
      if (ch != ' ') {
        dst = NULL;
      } else if (dst == NULL) {
        dst = ptr;
      }
      ptr++;
      ch = *ptr;
    }
    if (dst != NULL) {
      *dst = '\0';
    }
  }
  return str;
}

static void CopyLetters (CharPtr dest, CharPtr source, size_t maxsize)
{
  Char            ch;
  CharPtr         tmp;

  if (dest == NULL || maxsize < 1)
    return;
  *dest = '\0';
  if (source == NULL)
    return;
  maxsize--;
  tmp = dest;
  ch = *source;
  while (maxsize > 1 && ch != '\0') {
    if (ch != '.') {
      *dest = ch;
      dest++;
      maxsize--;
    }
    source++;
    ch = *source;
  }
  *dest = '\0';
  TrimSpacesOnEitherSide (tmp);
}

static void LookForEtAl (ValidStructPtr vsp, ValNodePtr tmp)
{
  AuthorPtr       ap;
  AuthListPtr     authors = NULL;
  CitArtPtr       cap;
  CitBookPtr      cbp;
  CitGenPtr       cgp;
  CitSubPtr       csp;
  Char            first[64];
  Char            initials[16];
  Char            last[64];
  ValNodePtr      names;
  NameStdPtr      nsp;
  PersonIdPtr     pid;

  if (vsp == NULL || tmp == NULL)
    return;
  switch (tmp->choice) {
  case PUB_Article:
    cap = (CitArtPtr) (tmp->data.ptrvalue);
    authors = cap->authors;
    break;
  case PUB_Man:
  case PUB_Book:
  case PUB_Proc:
    cbp = (CitBookPtr) (tmp->data.ptrvalue);
    authors = cbp->authors;
    break;
  case PUB_Gen:
    cgp = (CitGenPtr) (tmp->data.ptrvalue);
    authors = cgp->authors;
    break;
  case PUB_Sub:
    csp = (CitSubPtr) (tmp->data.ptrvalue);
    authors = csp->authors;
    break;
  default:
    break;
  }
  if (authors == NULL || authors->choice != 1)
    return;
  for (names = authors->names; names != NULL; names = names->next) {
    ap = names->data.ptrvalue;
    if (ap != NULL) {
      pid = ap->name;
      if (pid != NULL && pid->choice == 2) {
        nsp = pid->data;
        if (nsp != NULL && nsp->names[0] != NULL) {
          CopyLetters (last, nsp->names[0], sizeof (last));
          CopyLetters (first, nsp->names[1], sizeof (first));
          CopyLetters (initials, nsp->names[4], sizeof (initials));
          if ((StringICmp (last, "et al") == 0) || (StringCmp (initials, "al") == 0 && StringCmp (last, "et") == 0 && first[0] == '\0')) {
            if (names->next == NULL) {
              ValidErr (vsp, SEV_WARNING, ERR_GENERIC_AuthorListHasEtAl, "Author list ends in et al.");
            } else {
              ValidErr (vsp, SEV_WARNING, ERR_GENERIC_AuthorListHasEtAl, "Author list contains et al.");
            }
          }
        }
      }
    }
  }
}

static void SpellCheckPub (ValidStructPtr vsp, ValNodePtr tmp)
{
  CitArtPtr       cap;
  CitBookPtr      cbp;
  CitGenPtr       cgp;
  ValNodePtr      titles = NULL;

  if ((vsp == NULL) || (tmp == NULL))
    return;

  switch (tmp->choice) {
  case PUB_Article:
    cap = (CitArtPtr) (tmp->data.ptrvalue);
    titles = cap->title;
    break;
  case PUB_Man:
  case PUB_Book:
  case PUB_Proc:
    cbp = (CitBookPtr) (tmp->data.ptrvalue);
    titles = cbp->title;
    break;
  case PUB_Gen:
    cgp = (CitGenPtr) (tmp->data.ptrvalue);
    if (cgp->cit != NULL)
      SpellCheckString (vsp, cgp->cit);
    if (cgp->title != NULL)
      SpellCheckString (vsp, cgp->title);
    break;
  default:
    break;
  }

  if (titles != NULL) {
    for (; titles != NULL; titles = titles->next) {
      if (titles->choice == Cit_title_name)
        SpellCheckString (vsp, (CharPtr) (titles->data.ptrvalue));
    }
  }

  return;
}

static void SpellCheckSeqDescr (GatherContextPtr gcp)
{
  PubdescPtr      pdp;
  ValNodePtr      tmp, vnp;
  ValidStructPtr  vsp;

  vsp = (ValidStructPtr) (gcp->userdata);
  if (vsp == NULL)
    return;

  vnp = (ValNodePtr) (gcp->thisitem);
  if (vnp == NULL)
    return;

  vsp->descr = vnp;
  vsp->sfp = NULL;

  if (vnp->choice == Seq_descr_pub) {
    pdp = (PubdescPtr) (vnp->data.ptrvalue);
    for (tmp = pdp->pub; tmp != NULL; tmp = tmp->next) {
      LookForEtAl (vsp, tmp);
    }
  }

  if (vsp->spellfunc == NULL)
    return;

  switch (vnp->choice) {
  case Seq_descr_title:
  case Seq_descr_region:
  case Seq_descr_comment:
    SpellCheckString (vsp, (CharPtr) (vnp->data.ptrvalue));
    break;
  case Seq_descr_pub:
    pdp = (PubdescPtr) (vnp->data.ptrvalue);
    for (tmp = pdp->pub; tmp != NULL; tmp = tmp->next) {
      SpellCheckPub (vsp, tmp);
    }
    SpellCheckString (vsp, pdp->comment);
    break;
  default:
    break;
  }
  return;
}

NLM_EXTERN void SpellCheckSeqFeat (GatherContextPtr gcp)
{
  PubdescPtr      pdp;
  SeqFeatPtr      sfp;
  ProtRefPtr      prp;
  ValidStructPtr  vsp;
  ValNodePtr      vnp;

  vsp = (ValidStructPtr) (gcp->userdata);
  if (vsp == NULL)
    return;

  sfp = (SeqFeatPtr) (gcp->thisitem);
  if (sfp == NULL)
    return;

  vsp->descr = NULL;
  vsp->sfp = sfp;

  if (sfp->data.choice == SEQFEAT_PUB) {
    pdp = (PubdescPtr) (sfp->data.value.ptrvalue);
    for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
      LookForEtAl (vsp, vnp);
    }
  }

  if (vsp->spellfunc == NULL)
    return;

  SpellCheckString (vsp, sfp->comment);

  switch (sfp->data.choice) {
  case 1:                      /* Gene-ref */
    break;
  case 2:                      /* Org-ref */
    break;
  case 3:                      /* Cdregion */
    break;
  case 4:                      /* Prot-ref */
    prp = (ProtRefPtr) (sfp->data.value.ptrvalue);
    for (vnp = prp->name; vnp != NULL; vnp = vnp->next)
      SpellCheckString (vsp, (CharPtr) (vnp->data.ptrvalue));
    SpellCheckString (vsp, prp->desc);
    break;
  case 5:                      /* RNA-ref */
    break;
  case 6:                      /* Pub */
    pdp = (PubdescPtr) (sfp->data.value.ptrvalue);
    for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
      SpellCheckPub (vsp, vnp);
    }
    SpellCheckString (vsp, pdp->comment);
    break;
  case 7:                      /* Seq */
    break;
  case 8:                      /* Imp-feat */
    break;
  case 9:                      /* Region */
    SpellCheckString (vsp, (CharPtr) (sfp->data.value.ptrvalue));
    break;
  case 10:                     /* Comment */
    break;
  case 11:                     /* Bond */
    break;
  case 12:                     /* Site */
    break;
  case 13:                     /* Rsite-ref */
    break;
  case 14:                     /* User-object */
    break;
  case 15:                     /* TxInit */
    break;
  case 16:                     /* Numbering */
    break;
  case 17:                     /* Secondary Structure */
    break;
  case 18:                     /* NonStdRes */
    break;
  case 19:                     /* Heterogen */
    break;
  case 20:                     /* BioSource */
    break;
  default:
    break;
  }

  return;
}

NLM_EXTERN void SpellCheckString (ValidStructPtr vsp, CharPtr str)
{
  if ((vsp == NULL) || (str == NULL))
    return;

  if (vsp->spellfunc == NULL)
    return;

  (*(vsp->spellfunc)) ((char *) str, (vsp->spellcallback));

  return;
}

NLM_EXTERN void SpellCallBack (char *str)
{
  ErrSev          sev;

  sev = SEV_ERROR;
  if (globalvsp != NULL && globalvsp->justwarnonspell) {
    sev = SEV_WARNING;
  }
  ValidErr (globalvsp, sev, ERR_GENERIC_Spell, "[ %s ]", (CharPtr) str);
  return;
}
