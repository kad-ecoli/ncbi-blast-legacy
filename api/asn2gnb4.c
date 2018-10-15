/*   asn2gnb4.c
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
* File Name:  asn2gnb4.c
*
* Author:  Karl Sirotkin, Tom Madden, Tatiana Tatusov, Jonathan Kans,
*          Mati Shomrat
*
* Version Creation Date:   10/21/98
*
* $Revision: 1.258 $
*
* File Description:  New GenBank flatfile generator - work in progress
*
* Modifications:
* --------------------------------------------------------------------------
* ==========================================================================
*/

#include <ncbi.h>
#include <objall.h>
#include <objsset.h>
#include <objsub.h>
#include <objfdef.h>
#include <objpubme.h>
#include <seqport.h>
#include <sequtil.h>
#include <sqnutils.h>
#include <subutil.h>
#include <tofasta.h>
#include <explore.h>
#include <gbfeat.h>
#include <gbftdef.h>
#include <edutil.h>
#include <alignmgr2.h>
#include <asn2gnbi.h>
#include <valid.h>

#ifdef WIN_MAC
#if __profile__
#include <Profiler.h>
#endif
#endif

static CharPtr link_muid = "http://www.ncbi.nlm.nih.gov/pubmed/";

static CharPtr link_go = "http://amigo.geneontology.org/cgi-bin/amigo/go.cgi?view=details&depth=1&query=GO:";

static CharPtr link_go_ref = "http://www.geneontology.org/cgi-bin/references.cgi#GO_REF:";

static CharPtr link_code = "http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?";

static CharPtr link_featn = "http://www.ncbi.nlm.nih.gov/nuccore/";
static CharPtr link_featp = "http://www.ncbi.nlm.nih.gov/protein/";

static CharPtr link_seqn = "http://www.ncbi.nlm.nih.gov/nuccore/";
static CharPtr link_seqp = "http://www.ncbi.nlm.nih.gov/protein/";

/*
static CharPtr ec_link = "http://www.expasy.org/cgi-bin/nicezyme.pl?";

static CharPtr ec_ambig = "http://www.chem.qmw.ac.uk/iubmb/enzyme/";
*/

static CharPtr ec_link = "http://www.expasy.org/enzyme/";

/* ordering arrays for qualifiers and note components */

static FtQualType feat_qual_order [] = {
  FTQUAL_partial,
  FTQUAL_gene,

  FTQUAL_locus_tag,
  FTQUAL_old_locus_tag,

  FTQUAL_gene_syn_refseq,
  FTQUAL_gene_syn,

  FTQUAL_gene_allele,

  FTQUAL_operon,

  FTQUAL_ncRNA_class,
  FTQUAL_ncRNA_other,

  FTQUAL_product,

  FTQUAL_prot_EC_number,
  FTQUAL_prot_activity,

  FTQUAL_standard_name,
  FTQUAL_coded_by,
  FTQUAL_derived_from,

  FTQUAL_prot_name,
  FTQUAL_region_name,
  FTQUAL_bond_type,
  FTQUAL_site_type,
  FTQUAL_sec_str_type,
  FTQUAL_heterogen,

  FTQUAL_tag_peptide,
  FTQUAL_tag_peptide_str,

  FTQUAL_evidence,
  FTQUAL_experiment,
  FTQUAL_experiment_string,
  FTQUAL_inference,
  FTQUAL_inference_string,
  FTQUAL_inference_good,
  FTQUAL_exception,
  FTQUAL_ribosomal_slippage,
  FTQUAL_trans_splicing,
  FTQUAL_artificial_location,
  FTQUAL_artificial_location_str,

  FTQUAL_note, 
  FTQUAL_citation,

  FTQUAL_number,

  FTQUAL_pseudo,
  FTQUAL_selenocysteine,
  FTQUAL_pyrrolysine,

  FTQUAL_codon_start,

  FTQUAL_anticodon,
  FTQUAL_trna_codons,
  FTQUAL_bound_moiety,
  FTQUAL_clone,
  FTQUAL_compare,
  FTQUAL_direction,
  FTQUAL_function,
  FTQUAL_frequency,
  FTQUAL_EC_number,
  FTQUAL_gene_map,
  FTQUAL_gene_cyt_map,
  FTQUAL_gene_gen_map,
  FTQUAL_gene_rad_map,
  FTQUAL_estimated_length,
  FTQUAL_gap_type,
  FTQUAL_linkage_evidence,
  FTQUAL_allele,
  FTQUAL_map,
  FTQUAL_mod_base,
  FTQUAL_PCR_conditions,
  FTQUAL_phenotype,
  FTQUAL_rpt_family,
  FTQUAL_rpt_type,
  FTQUAL_rpt_unit,
  FTQUAL_rpt_unit_range,
  FTQUAL_rpt_unit_seq,
  FTQUAL_satellite,
  FTQUAL_mobile_element,
  FTQUAL_mobile_element_type,
  FTQUAL_usedin,

  FTQUAL_illegal_qual,

  FTQUAL_replace,
  FTQUAL_delta_item,
  FTQUAL_variation_set,

  FTQUAL_transl_except,
  FTQUAL_transl_table,
  FTQUAL_codon,
  FTQUAL_organism,
  FTQUAL_label,
  FTQUAL_cds_product,
  FTQUAL_extra_products,
  FTQUAL_UniProtKB_evidence,
  FTQUAL_protein_id,
  FTQUAL_transcript_id,
  FTQUAL_db_xref, 
  FTQUAL_gene_xref,
  FTQUAL_variation_id,
  FTQUAL_mol_wt,
  FTQUAL_translation,
  FTQUAL_transcription,
  FTQUAL_peptide,
  (FtQualType) 0
};

/*
prot_names after seqfeat_note - gi|4210642|emb|AJ005084.1|HBVAJ5084
prot_conflict after prot_desc - gi|61183|emb|V01135.1|PIVM02
figure after prot_desc - gi|400553|gb|S64006.1|
seqfeat_note after prot_desc - gi|431713|gb|L20354.1|STVPATPOLB
  but prot_desc after seqfeat_note - AF252556.1
prot_names after figure - gi|234022|gb|S56149.1|S56149
seqfeat_note after prot_conflict after figure - gi|234046|gb|S51392.1|S51392
prot_method after prot_comment (descriptor) after prot_note after prot_desc
region after seqfeat_note - gi|6554164|gb|AF043644.3|AF043644
prot_desc after prot_names - gi|6581069|gb|AF202541.1|AF202541 - cannot do !!!
gene_syn after gene_desc - gi|3386543|gb|AF079528.1|AF079528
pseudo after note - gi|6598562|gb|AC006419.3|AC006419
*/

static FtQualType feat_note_order [] = {
  FTQUAL_transcript_id_note, /* !!! remove October 15, 2003 !!! */
  FTQUAL_gene_desc,
  FTQUAL_trna_codons_note,
  FTQUAL_encodes,
  FTQUAL_prot_desc,
  FTQUAL_prot_note,
  FTQUAL_prot_comment,
  FTQUAL_prot_method,
  FTQUAL_ncRNA_note,
  FTQUAL_figure,
  FTQUAL_maploc,
  FTQUAL_prot_conflict,
  FTQUAL_prot_missing,
  FTQUAL_seqfeat_note,
  FTQUAL_seqannot_note,
  FTQUAL_region,
  FTQUAL_selenocysteine_note,
  FTQUAL_pyrrolysine_note,
  FTQUAL_prot_names,
  FTQUAL_bond,
  FTQUAL_site,
  /*
  FTQUAL_rrna_its,
  */
  FTQUAL_xtra_prod_quals,
  FTQUAL_inference_bad,
  FTQUAL_modelev,
  FTQUAL_cdd_definition,
  /* GO terms appear as own qualifiers in RefSeq records, Sequin or Dump mode */
  FTQUAL_go_component,
  FTQUAL_go_function,
  FTQUAL_go_process,
  /* RefSeq-specific qualifiers have same display policy as GO terms */
  FTQUAL_nomenclature,
  FTQUAL_gene_nomen,
  FTQUAL_exception_note,
  (FtQualType) 0
};

typedef struct featurqual {
  CharPtr   name;
  QualType  qualclass;
} FeaturQual, PNTR FeaturQualPtr;

static FeaturQual asn2gnbk_featur_quals [ASN2GNBK_TOTAL_FEATUR]  =  {
  { "",                    Qual_class_ignore         },
  { "allele",              Qual_class_quote          },
  { "anticodon",           Qual_class_anti_codon     },
  { "artificial_location", Qual_class_boolean        },
  { "artificial_location", Qual_class_string         },
  { "bond",                Qual_class_bond           },
  { "bond_type",           Qual_class_bond           },
  { "bound_moiety",        Qual_class_quote          },
  { "cdd_definition",      Qual_class_string         },
  { "product",             Qual_class_string         },
  { "citation",            Qual_class_pubset         },
  { "clone",               Qual_class_quote          },
  { "coded_by",            Qual_class_seq_loc        },
  { "compare",             Qual_class_compare        },
  { "codon",               Qual_class_codon          },
  { "codon_start",         Qual_class_int            },
  { "cons_splice",         Qual_class_consplice      },
  { "db_xref",             Qual_class_db_xref        },
  { "delta_item",          Qual_class_delta_item     },
  { "derived_from",        Qual_class_seq_loc        },
  { "direction",           Qual_class_L_R_B          },
  { "EC_number",           Qual_class_EC_quote       },
  { "encodes",             Qual_class_encodes        },
  { "estimated_length",    Qual_class_number         },
  { "evidence",            Qual_class_evidence       },
  { "exception",           Qual_class_exception      },
  { "exception_note",      Qual_class_exception      },
  { "experiment",          Qual_class_experiment     },
  { "experiment",          Qual_class_string         },
  { "product",             Qual_class_valnode        },
  { "figure",              Qual_class_string         },
  { "frequency",           Qual_class_quote          },
  { "function",            Qual_class_quote          },
  { "gap_type",            Qual_class_quote          },
  { "gene",                Qual_class_sgml           },
  { "gene_desc",           Qual_class_string         },
  { "allele",              Qual_class_string         },
  { "map",                 Qual_class_string         },
  { "cyt_map",             Qual_class_map            },
  { "gen_map",             Qual_class_map            },
  { "rad_map",             Qual_class_map            },
  { "gene_synonym",        Qual_class_sep_gene_syn   },
  { "gene_synonym",        Qual_class_gene_syn       },
  { "gene_note",           Qual_class_string         },
  { "db_xref",             Qual_class_db_xref        },
  { "GO_component",        Qual_class_go             },
  { "GO_function",         Qual_class_go             },
  { "GO_process",          Qual_class_go             },
  { "heterogen",           Qual_class_string         },
  { "illegal",             Qual_class_illegal        },
  { "inference",           Qual_class_quote          },
  { "inference",           Qual_class_string         },
  { "inference",           Qual_class_valnode        },
  { "inference",           Qual_class_valnode        },
  { "insertion_seq",       Qual_class_quote          },
  { "label",               Qual_class_label          },
  { "linkage_evidence",    Qual_class_quote          },
  { "locus_tag",           Qual_class_locus_tag      },
  { "map",                 Qual_class_quote          },
  { "maploc",              Qual_class_string         },
  { "mobile_element",      Qual_class_mobile_element },
  { "mobile_element_type", Qual_class_mobile_element },
  { "mod_base",            Qual_class_noquote        },
  { "model_evidence",      Qual_class_model_ev       },
  { "calculated_mol_wt",   Qual_class_mol_wt         },
  { "ncRNA_class",         Qual_class_quote          },
  { "ncRNA_note",          Qual_class_string         },
  { "ncRNA_class",         Qual_class_string         },
  { "nomenclature",        Qual_class_nomenclature   },
  { "nomenclature",        Qual_class_gene_nomen     },
  { "note",                Qual_class_note           },
  { "number",              Qual_class_number         },
  { "old_locus_tag",       Qual_class_paren          },
  { "operon",              Qual_class_quote          },
  { "organism",            Qual_class_quote          },
  { "partial",             Qual_class_boolean        },
  { "PCR_conditions",      Qual_class_quote          },
  { "peptide",             Qual_class_peptide        },
  { "phenotype",           Qual_class_quote          },
  { "product",             Qual_class_product        },
  { "product",             Qual_class_quote          },
  { "function",            Qual_class_valnode        },
  { "prot_comment",        Qual_class_string         },
  { "EC_number",           Qual_class_EC_valnode     },
  { "prot_note",           Qual_class_string         },
  { "prot_method",         Qual_class_method         },
  { "prot_conflict",       Qual_class_string         },
  { "prot_desc",           Qual_class_string         },
  { "prot_missing",        Qual_class_string         },
  { "name",                Qual_class_tilde          },
  { "prot_names",          Qual_class_protnames      },
  { "protein_id",          Qual_class_prt_id         },
  { "pseudo",              Qual_class_boolean        },
  { "pyrrolysine",         Qual_class_boolean        },
  { "pyrrolysine",         Qual_class_string         },
  { "region",              Qual_class_region         },
  { "region_name",         Qual_class_string         },
  { "replace",             Qual_class_replace        },
  { "ribosomal_slippage",  Qual_class_boolean        },
  { "rpt_family",          Qual_class_quote          },
  { "rpt_type",            Qual_class_rpt            },
  { "rpt_unit",            Qual_class_rpt_unit       },
  { "rpt_unit_range",      Qual_class_rpt_unit       },
  { "rpt_unit_seq",        Qual_class_rpt_unit       },
  { "rrna_its",            Qual_class_its            },
  { "satellite",           Qual_class_quote          },
  { "sec_str_type",        Qual_class_sec_str        },
  { "selenocysteine",      Qual_class_boolean        },
  { "selenocysteine",      Qual_class_string         },
  { "seqannot_note",       Qual_class_string         },
  { "seqfeat_note",        Qual_class_string         },
  { "site",                Qual_class_site           },
  { "site_type",           Qual_class_site           },
  { "standard_name",       Qual_class_quote          },
  { "tag_peptide",         Qual_class_noquote        },
  { "tag_peptide",         Qual_class_tag_peptide    },
  { "transcription",       Qual_class_transcription  },
  { "transcript_id",       Qual_class_nuc_id         },
  { "tscpt_id_note",       Qual_class_nuc_id         },
  { "transl_except",       Qual_class_code_break     },
  { "transl_table",        Qual_class_int            },
  { "translation",         Qual_class_translation    },
  { "transposon",          Qual_class_quote          },
  { "trans_splicing",      Qual_class_boolean        },
  { "trna_aa",             Qual_class_ignore         },
  { "codon_recognized",    Qual_class_trna_codons    },
  { "trna_codons",         Qual_class_trna_codons    },
  { "UniProtKB_evidence",  Qual_class_quote          },
  { "usedin",              Qual_class_usedin         },
  { "db_xref",             Qual_class_variation_id   },
  { "variation_set",       Qual_class_variation_set  },
  { "xtra_products",       Qual_class_xtraprds       }
};


typedef struct qualfeatur {
  CharPtr     name;
  FtQualType  featurclass;
} QualFeatur, PNTR QualFeaturPtr;

#define NUM_GB_QUALS 44

static QualFeatur qualToFeature [NUM_GB_QUALS] = {
  { "allele",              FTQUAL_allele              },
  { "bound_moiety",        FTQUAL_bound_moiety        },
  { "clone",               FTQUAL_clone               },
  { "codon",               FTQUAL_codon               },
  { "compare",             FTQUAL_compare             },
  { "cons_splice",         FTQUAL_cons_splice         },
  { "cyt_map",             FTQUAL_gene_cyt_map        },
  { "direction",           FTQUAL_direction           },
  { "EC_number",           FTQUAL_EC_number           },
  { "estimated_length",    FTQUAL_estimated_length    },
  { "experiment",          FTQUAL_experiment          },
  { "frequency",           FTQUAL_frequency           },
  { "function",            FTQUAL_function            },
  { "gap_type",            FTQUAL_gap_type            },
  { "gen_map",             FTQUAL_gene_gen_map        },
  { "inference",           FTQUAL_inference           },
  { "insertion_seq",       FTQUAL_insertion_seq       },
  { "label",               FTQUAL_label               },
  { "linkage_evidence",    FTQUAL_linkage_evidence    },
  { "map",                 FTQUAL_map                 },
  { "mobile_element",      FTQUAL_mobile_element      },
  { "mobile_element_type", FTQUAL_mobile_element_type },
  { "mod_base",            FTQUAL_mod_base            },
  { "ncRNA_class",         FTQUAL_ncRNA_class         },
  { "number",              FTQUAL_number              },
  { "old_locus_tag",       FTQUAL_old_locus_tag       },
  { "operon",              FTQUAL_operon              },
  { "organism",            FTQUAL_organism            },
  { "PCR_conditions",      FTQUAL_PCR_conditions      },
  { "phenotype",           FTQUAL_phenotype           },
  { "product",             FTQUAL_product_quals       },
  { "rad_map",             FTQUAL_gene_rad_map        },
  { "replace",             FTQUAL_replace             },
  { "rpt_family",          FTQUAL_rpt_family          },
  { "rpt_type",            FTQUAL_rpt_type            },
  { "rpt_unit",            FTQUAL_rpt_unit            },
  { "rpt_unit_range",      FTQUAL_rpt_unit_range      },
  { "rpt_unit_seq",        FTQUAL_rpt_unit_seq        },
  { "satellite",           FTQUAL_satellite           },
  { "standard_name",       FTQUAL_standard_name       },
  { "tag_peptide",         FTQUAL_tag_peptide         },
  { "transposon",          FTQUAL_transposon          },
  { "UniProtKB_evidence",  FTQUAL_UniProtKB_evidence  },
  { "usedin",              FTQUAL_usedin              }
};

static Int2 GbqualToFeaturIndex (
  CharPtr qualname
)

{
  Int2  L, R, mid;

  if (qualname == NULL || *qualname == '\0') return 0;

  L = 0;
  R = NUM_GB_QUALS - 1;

  while (L < R) {
    mid = (L + R) / 2;
    if (StringICmp (qualToFeature [mid].name, qualname) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  if (StringICmp (qualToFeature [R].name, qualname) == 0) {
    return qualToFeature [R].featurclass;
  }

  return 0;
}

#define NUM_ILLEGAL_QUALS 14

static FeaturQual illegalGbqualList [NUM_ILLEGAL_QUALS] = {
  { "anticodon",      Qual_class_noquote },
  { "citation",       Qual_class_noquote },
  { "codon_start",    Qual_class_noquote },
  { "db_xref",        Qual_class_quote   },
  { "evidence",       Qual_class_noquote },
  { "exception",      Qual_class_quote   },
  { "gene",           Qual_class_quote   },
  { "note",           Qual_class_quote   },
  { "protein_id",     Qual_class_quote   },
  { "pseudo",         Qual_class_noquote },
  { "transcript_id",  Qual_class_quote   },
  { "transl_except",  Qual_class_noquote },
  { "transl_table",   Qual_class_noquote },
  { "translation",    Qual_class_quote   }
};

static Int2 IllegalGbqualToClass (
  CharPtr qualname
)

{
  Int2  L, R, mid;

  if (qualname == NULL || *qualname == '\0') return 0;

  L = 0;
  R = NUM_ILLEGAL_QUALS - 1;

  while (L < R) {
    mid = (L + R) / 2;
    if (StringICmp (illegalGbqualList [mid].name, qualname) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  if (StringICmp (illegalGbqualList [R].name, qualname) == 0) {
    return illegalGbqualList [R].qualclass;
  }

  return 0;
}

static CharPtr trnaList [] = {
  "tRNA-Gap",
  "tRNA-Ala",
  "tRNA-Asx",
  "tRNA-Cys",
  "tRNA-Asp",
  "tRNA-Glu",
  "tRNA-Phe",
  "tRNA-Gly",
  "tRNA-His",
  "tRNA-Ile",
  "tRNA-Xle",
  "tRNA-Lys",
  "tRNA-Leu",
  "tRNA-Met",
  "tRNA-Asn",
  "tRNA-Pyl",
  "tRNA-Pro",
  "tRNA-Gln",
  "tRNA-Arg",
  "tRNA-Ser",
  "tRNA-Thr",
  "tRNA-Sec",
  "tRNA-Val",
  "tRNA-Trp",
  "tRNA-OTHER",
  "tRNA-Tyr",
  "tRNA-Glx",
  "tRNA-TERM",
  NULL
};

static CharPtr evidenceText [] = {
  NULL, "experimental", "not_experimental"
};

NLM_EXTERN CharPtr secStrText [] = {
  NULL, "helix", "sheet", "turn"
};

static CharPtr  oops = "?";

static CharPtr SeqCodeNameGet (
  SeqCodeTablePtr table,
  Uint1 residue
)

{
  Uint1  index;

  if (table != NULL) {
    index = residue - table->start_at;
    if ( /*index >= 0 && */ index < table->num) {
      return (table->names) [index];
    }
  }

  return oops;
}

NLM_EXTERN CharPtr Get3LetterSymbol (
  IntAsn2gbJobPtr  ajp,
  Uint1 seq_code,
  SeqCodeTablePtr table,
  Uint1 residue
)

{
  Uint1            code = Seq_code_ncbieaa;
  Int2             index;
  Uint1            new_residue;
  CharPtr          ptr;
  CharPtr          retval = NULL;
  SeqMapTablePtr   smtp;
  SeqCodeTablePtr  table_3aa;

  if (residue == 42) { /* stop codon in NCBIeaa */
    retval = "TERM";
    return retval;
  }

  if (ajp != NULL && ajp->flags.iupacaaOnly) {
    code = Seq_code_iupacaa;
  } else {
    code = Seq_code_ncbieaa;
  }

  if (seq_code != code) {
    /* if code and seq_code are identical, then smtp is NULL?? */
    smtp = SeqMapTableFind (code, seq_code);
    new_residue = SeqMapTableConvert (smtp, residue);
  } else {
    new_residue = residue;
  }

  /* The following looks for non-symbols (255) and "Undetermined" (88) */
  if ((int) new_residue == 255 || (int) new_residue == 88) {
    retval = "OTHER";
    return retval;
  } else {
    if (table == NULL) {
      table = SeqCodeTableFind (Seq_code_ncbieaa);
      if (table == NULL) {
        retval = "OTHER";
        return retval;
      }
    }
    ptr = SeqCodeNameGet (table, residue);
    table_3aa = SeqCodeTableFind (Seq_code_iupacaa3);
    if (ptr != NULL && table_3aa != NULL) {
      for (index=0; index < (int) table_3aa->num; index++) {
        if (StringCmp(ptr, (table_3aa->names) [index]) == 0) {
          retval = (table_3aa->symbols) [index];
          return retval;
        }
      }
    }
  }

  retval = "OTHER";
  return retval;
}

static Boolean MatchCit (
  ValNodePtr ppr,
  RefBlockPtr rbp
)

{
  Char        buf [121];
  size_t      len;
  Int4        uid;
  ValNodePtr  vnp;

  if (ppr == NULL || rbp == NULL) return FALSE;
  switch (ppr->choice) {
    case PUB_Muid :
      uid = ppr->data.intvalue;
      if (rbp->muid == uid) return TRUE;
      break;
    case PUB_PMid :
      uid = ppr->data.intvalue;
      if (rbp->pmid == uid) return TRUE;
      break;
    case PUB_Equiv :
      for (vnp = (ValNodePtr) ppr->data.ptrvalue; vnp != NULL; vnp = vnp->next) {
        if (MatchCit (vnp, rbp)) return TRUE;
      }
      break;
    default :
      PubLabelUnique (ppr, buf, sizeof (buf) - 1, OM_LABEL_CONTENT, TRUE);
      len = StringLen (buf);
      if (len > 0 && buf [len - 1] == '>') {
        buf [len - 1] = '\0';
        len--;
      }
      len = MIN (len, StringLen (rbp->uniquestr));
      if (StringNICmp (rbp->uniquestr, buf, len) == 0) return TRUE;
      break;
  }
  return FALSE;
}

NLM_EXTERN Int2 MatchRef (
  ValNodePtr ppr,
  RefBlockPtr PNTR rbpp,
  Int2 numReferences
)

{
  Int2         j;
  RefBlockPtr  rbp;

  if (ppr == NULL || rbpp == NULL) return 0;

  for (j = 0; j < numReferences; j++) {
    rbp = rbpp [j];
    if (rbp == NULL) continue;
    if (MatchCit (ppr, rbp)) return rbp->serial;
  }
  return 0;
}

/* taken from asn2ff4.c */

static Boolean LookForFuzz (SeqLocPtr head)
{
  Boolean retval=FALSE;
  IntFuzzPtr ifp;
  PackSeqPntPtr pspp;
  SeqIntPtr sip;
  SeqLocPtr slp;
  SeqPntPtr spp;

  if (head == NULL)
    return retval;

  slp=NULL;
  while ((slp = SeqLocFindNext(head, slp)) != NULL)
  {
    switch (slp->choice)
    {
      case SEQLOC_INT:
        sip = (SeqIntPtr)(slp->data.ptrvalue);
        ifp = sip->if_from;
        if (ifp != NULL)
        {
          if (ifp->choice == 4)
          {
            if (ifp->a != 0)
              retval=TRUE;
          }
          else
            retval = TRUE;
        }
        ifp = sip->if_to;
        if (ifp != NULL)
        {
          if (ifp->choice == 4)
          {
            if (ifp->a != 0)
              retval=TRUE;
          }
          else
            retval = TRUE;
        }
        break;
      case SEQLOC_PNT:
        spp = (SeqPntPtr)(slp->data.ptrvalue);
        ifp = spp->fuzz;
        if (ifp != NULL)
        {
          if (ifp->choice == 4)
          {
            if (ifp->a != 0)
              retval=TRUE;
          }
          else
            retval = TRUE;
        }
        break;
      case SEQLOC_PACKED_PNT:
        pspp = (PackSeqPntPtr)(slp->data.ptrvalue);
        ifp = pspp->fuzz;
        if (ifp != NULL)
        {
          if (ifp->choice == 4)
          {
            if (ifp->a != 0)
              retval=TRUE;
          }
          else
            retval = TRUE;
        }
        break;
      default:
        break;
    }
    if (retval == TRUE)
      break;
  }
  return retval;
}

NLM_EXTERN CharPtr bondList [] = {
  NULL,
  "disulfide",
  "thiolester",
  "xlink",
  "thioether",
  "other"
};

NLM_EXTERN CharPtr siteList [] = {
  NULL,
  "active",
  "binding",
  "cleavage",
  "inhibit",
  "modified",
  "glycosylation",
  "myristoylation",
  "mutagenized",
  "metal-binding",
  "phosphorylation",
  "acetylation",
  "amidation",
  "methylation",
  "hydroxylation",
  "sulfatation",
  "oxidative-deamination",
  "pyrrolidone-carboxylic-acid",
  "gamma-carboxyglutamic-acid",
  "blocked",
  "lipid-binding",
  "np-binding",
  "DNA binding",
  "signal-peptide",
  "transit-peptide",
  "transmembrane-region",
  "nitrosylation",
  "other"
};

static CharPtr siteFFList [] = {
  NULL,
  "active",
  "binding",
  "cleavage",
  "inhibition",
  "modified",
  "glycosylation",
  "myristoylation",
  "mutagenized",
  "metal-binding",
  "phosphorylation",
  "acetylation",
  "amidation",
  "methylation",
  "hydroxylation",
  "sulfatation",
  "oxidative-deamination",
  "pyrrolidone-carboxylic-acid",
  "gamma-carboxyglutamic-acid",
  "blocked",
  "lipid-binding",
  "np-binding",
  "DNA binding",
  "signal peptide",
  "transit peptide",
  "transmembrane region",
  "nitrosylation",
  "other"
};

static CharPtr conflict_msg =
"Protein sequence is in conflict with the conceptual translation";

/*
static CharPtr no_protein_msg =
"Coding region translates with internal stops";
*/

/**/
/*  s_DisplayQVP () -- Displays the strings in a QVP structure.   */
/*                     This is a debugging function only.         */
/**/

#ifdef DISPLAY_STRINGS
static void s_DisplayQVP(QualValPtr qvp, Uint1Ptr notetbl)
{
  Int2 j;
  Int2 jdx;

  fprintf(stderr,"\n");
  for (j = 0, jdx = notetbl [j]; jdx != 0; j++, jdx = notetbl [j])
    {
      if (((int) qvp[jdx].str != 0x1000000) &&
      ((int) qvp[jdx].str != 0x1) &&
      ((int) qvp[jdx].str != 0xb) &&
      (qvp[jdx].str != NULL))
    fprintf(stderr, "%d\t%-25s %s\n", j, asn2gnbk_featur_quals[jdx].name,
        qvp[jdx].str);
      else
    fprintf(stderr, "%d\t%-25s %s\n", j, asn2gnbk_featur_quals[jdx].name,
        "NULL");
    }
}
#endif

/*
static Boolean NotInGeneSyn (
  CharPtr str,
  ValNodePtr gene_syn)

{
  CharPtr     syn;
  ValNodePtr  vnp;

  for (vnp = gene_syn; vnp != NULL; vnp = vnp->next) {
    syn = (CharPtr) vnp->data.ptrvalue;
    if (! StringHasNoText (syn)) {
      if (StringICmp (str, syn) == 0) return FALSE;
    }
  }
  return TRUE;
}
*/

typedef struct valqualstruc {
  Uint2       featdef;
  FtQualType  ftqual;
} ValQual, PNTR ValQualPtr;

/*
   WARNING - This list MUST be kept sorted in FEATDEF order as the primary
   key, and within a FEATDEF group sorted by FTQUAL as the secondary key
*/

static ValQual legalGbqualList [] = {

  { FEATDEF_GENE , FTQUAL_allele },
  { FEATDEF_GENE , FTQUAL_function },
  { FEATDEF_GENE , FTQUAL_label },
  { FEATDEF_GENE , FTQUAL_map },
  { FEATDEF_GENE , FTQUAL_old_locus_tag },
  { FEATDEF_GENE , FTQUAL_operon },
  { FEATDEF_GENE , FTQUAL_phenotype },
  { FEATDEF_GENE , FTQUAL_product },
  { FEATDEF_GENE , FTQUAL_standard_name },

  { FEATDEF_CDS , FTQUAL_allele },
  { FEATDEF_CDS , FTQUAL_codon },
  { FEATDEF_CDS , FTQUAL_label },
  { FEATDEF_CDS , FTQUAL_map },
  { FEATDEF_CDS , FTQUAL_number },
  { FEATDEF_CDS , FTQUAL_old_locus_tag },
  { FEATDEF_CDS , FTQUAL_operon },
  { FEATDEF_CDS , FTQUAL_standard_name },

  { FEATDEF_PROT , FTQUAL_product },
  { FEATDEF_PROT , FTQUAL_UniProtKB_evidence },

  { FEATDEF_preRNA , FTQUAL_allele },
  { FEATDEF_preRNA , FTQUAL_function },
  { FEATDEF_preRNA , FTQUAL_label },
  { FEATDEF_preRNA , FTQUAL_map },
  { FEATDEF_preRNA , FTQUAL_old_locus_tag },
  { FEATDEF_preRNA , FTQUAL_operon },
  { FEATDEF_preRNA , FTQUAL_product },
  { FEATDEF_preRNA , FTQUAL_standard_name },

  { FEATDEF_mRNA , FTQUAL_allele },
  { FEATDEF_mRNA , FTQUAL_function },
  { FEATDEF_mRNA , FTQUAL_label },
  { FEATDEF_mRNA , FTQUAL_map },
  { FEATDEF_mRNA , FTQUAL_old_locus_tag },
  { FEATDEF_mRNA , FTQUAL_operon },
  { FEATDEF_mRNA , FTQUAL_product },
  { FEATDEF_mRNA , FTQUAL_standard_name },

  { FEATDEF_tRNA , FTQUAL_allele },
  { FEATDEF_tRNA , FTQUAL_function },
  { FEATDEF_tRNA , FTQUAL_label },
  { FEATDEF_tRNA , FTQUAL_map },
  { FEATDEF_tRNA , FTQUAL_old_locus_tag },
  { FEATDEF_tRNA , FTQUAL_operon },
  { FEATDEF_tRNA , FTQUAL_product },
  { FEATDEF_tRNA , FTQUAL_standard_name },

  { FEATDEF_rRNA , FTQUAL_allele },
  { FEATDEF_rRNA , FTQUAL_function },
  { FEATDEF_rRNA , FTQUAL_label },
  { FEATDEF_rRNA , FTQUAL_map },
  { FEATDEF_rRNA , FTQUAL_old_locus_tag },
  { FEATDEF_rRNA , FTQUAL_operon },
  { FEATDEF_rRNA , FTQUAL_product },
  { FEATDEF_rRNA , FTQUAL_standard_name },

  { FEATDEF_snRNA , FTQUAL_allele },
  { FEATDEF_snRNA , FTQUAL_function },
  { FEATDEF_snRNA , FTQUAL_label },
  { FEATDEF_snRNA , FTQUAL_map },
  { FEATDEF_snRNA , FTQUAL_old_locus_tag },
  { FEATDEF_snRNA , FTQUAL_product },
  { FEATDEF_snRNA , FTQUAL_standard_name },

  { FEATDEF_scRNA , FTQUAL_allele },
  { FEATDEF_scRNA , FTQUAL_function },
  { FEATDEF_scRNA , FTQUAL_label },
  { FEATDEF_scRNA , FTQUAL_map },
  { FEATDEF_scRNA , FTQUAL_old_locus_tag },
  { FEATDEF_scRNA , FTQUAL_product },
  { FEATDEF_scRNA , FTQUAL_standard_name },

  { FEATDEF_otherRNA , FTQUAL_allele },
  { FEATDEF_otherRNA , FTQUAL_function },
  { FEATDEF_otherRNA , FTQUAL_label },
  { FEATDEF_otherRNA , FTQUAL_map },
  { FEATDEF_otherRNA , FTQUAL_old_locus_tag },
  { FEATDEF_otherRNA , FTQUAL_operon },
  { FEATDEF_otherRNA , FTQUAL_product },
  { FEATDEF_otherRNA , FTQUAL_standard_name },

  { FEATDEF_attenuator , FTQUAL_allele },
  { FEATDEF_attenuator , FTQUAL_label },
  { FEATDEF_attenuator , FTQUAL_map },
  { FEATDEF_attenuator , FTQUAL_old_locus_tag },
  { FEATDEF_attenuator , FTQUAL_operon },
  { FEATDEF_attenuator , FTQUAL_phenotype },

  { FEATDEF_C_region , FTQUAL_allele },
  { FEATDEF_C_region , FTQUAL_label },
  { FEATDEF_C_region , FTQUAL_map },
  { FEATDEF_C_region , FTQUAL_old_locus_tag },
  { FEATDEF_C_region , FTQUAL_product },
  { FEATDEF_C_region , FTQUAL_standard_name },

  { FEATDEF_CAAT_signal , FTQUAL_allele },
  { FEATDEF_CAAT_signal , FTQUAL_label },
  { FEATDEF_CAAT_signal , FTQUAL_map },
  { FEATDEF_CAAT_signal , FTQUAL_old_locus_tag },

  { FEATDEF_Imp_CDS , FTQUAL_codon },
  { FEATDEF_Imp_CDS , FTQUAL_EC_number },
  { FEATDEF_Imp_CDS , FTQUAL_function },
  { FEATDEF_Imp_CDS , FTQUAL_label },
  { FEATDEF_Imp_CDS , FTQUAL_map },
  { FEATDEF_Imp_CDS , FTQUAL_number },
  { FEATDEF_Imp_CDS , FTQUAL_old_locus_tag },
  { FEATDEF_Imp_CDS , FTQUAL_operon },
  { FEATDEF_Imp_CDS , FTQUAL_product },
  { FEATDEF_Imp_CDS , FTQUAL_standard_name },

  { FEATDEF_conflict , FTQUAL_allele },
  { FEATDEF_conflict , FTQUAL_compare },
  { FEATDEF_conflict , FTQUAL_label },
  { FEATDEF_conflict , FTQUAL_map },
  { FEATDEF_conflict , FTQUAL_old_locus_tag },
  { FEATDEF_conflict , FTQUAL_replace },

  { FEATDEF_D_loop , FTQUAL_allele },
  { FEATDEF_D_loop , FTQUAL_label },
  { FEATDEF_D_loop , FTQUAL_map },
  { FEATDEF_D_loop , FTQUAL_old_locus_tag },

  { FEATDEF_D_segment , FTQUAL_allele },
  { FEATDEF_D_segment , FTQUAL_label },
  { FEATDEF_D_segment , FTQUAL_map },
  { FEATDEF_D_segment , FTQUAL_old_locus_tag },
  { FEATDEF_D_segment , FTQUAL_product },
  { FEATDEF_D_segment , FTQUAL_standard_name },

  { FEATDEF_enhancer , FTQUAL_allele },
  { FEATDEF_enhancer , FTQUAL_bound_moiety },
  { FEATDEF_enhancer , FTQUAL_label },
  { FEATDEF_enhancer , FTQUAL_map },
  { FEATDEF_enhancer , FTQUAL_old_locus_tag },
  { FEATDEF_enhancer , FTQUAL_standard_name },

  { FEATDEF_exon , FTQUAL_allele },
  { FEATDEF_exon , FTQUAL_EC_number },
  { FEATDEF_exon , FTQUAL_function },
  { FEATDEF_exon , FTQUAL_label },
  { FEATDEF_exon , FTQUAL_map },
  { FEATDEF_exon , FTQUAL_number },
  { FEATDEF_exon , FTQUAL_old_locus_tag },
  { FEATDEF_exon , FTQUAL_product },
  { FEATDEF_exon , FTQUAL_standard_name },

  { FEATDEF_GC_signal , FTQUAL_allele },
  { FEATDEF_GC_signal , FTQUAL_label },
  { FEATDEF_GC_signal , FTQUAL_map },
  { FEATDEF_GC_signal , FTQUAL_old_locus_tag },

  { FEATDEF_iDNA , FTQUAL_allele },
  { FEATDEF_iDNA , FTQUAL_function },
  { FEATDEF_iDNA , FTQUAL_label },
  { FEATDEF_iDNA , FTQUAL_map },
  { FEATDEF_iDNA , FTQUAL_number },
  { FEATDEF_iDNA , FTQUAL_old_locus_tag },
  { FEATDEF_iDNA , FTQUAL_standard_name },

  { FEATDEF_intron , FTQUAL_allele },
  { FEATDEF_intron , FTQUAL_cons_splice },
  { FEATDEF_intron , FTQUAL_function },
  { FEATDEF_intron , FTQUAL_label },
  { FEATDEF_intron , FTQUAL_map },
  { FEATDEF_intron , FTQUAL_number },
  { FEATDEF_intron , FTQUAL_old_locus_tag },
  { FEATDEF_intron , FTQUAL_standard_name },

  { FEATDEF_J_segment , FTQUAL_allele },
  { FEATDEF_J_segment , FTQUAL_label },
  { FEATDEF_J_segment , FTQUAL_map },
  { FEATDEF_J_segment , FTQUAL_old_locus_tag },
  { FEATDEF_J_segment , FTQUAL_product },
  { FEATDEF_J_segment , FTQUAL_standard_name },

  { FEATDEF_LTR , FTQUAL_allele },
  { FEATDEF_LTR , FTQUAL_function },
  { FEATDEF_LTR , FTQUAL_label },
  { FEATDEF_LTR , FTQUAL_map },
  { FEATDEF_LTR , FTQUAL_old_locus_tag },
  { FEATDEF_LTR , FTQUAL_standard_name },

  { FEATDEF_mat_peptide , FTQUAL_allele },
  { FEATDEF_mat_peptide , FTQUAL_EC_number },
  { FEATDEF_mat_peptide , FTQUAL_function },
  { FEATDEF_mat_peptide , FTQUAL_label },
  { FEATDEF_mat_peptide , FTQUAL_map },
  { FEATDEF_mat_peptide , FTQUAL_old_locus_tag },
  { FEATDEF_mat_peptide , FTQUAL_product },
  { FEATDEF_mat_peptide , FTQUAL_standard_name },

  { FEATDEF_misc_binding , FTQUAL_allele },
  { FEATDEF_misc_binding , FTQUAL_bound_moiety },
  { FEATDEF_misc_binding , FTQUAL_function },
  { FEATDEF_misc_binding , FTQUAL_label },
  { FEATDEF_misc_binding , FTQUAL_map },
  { FEATDEF_misc_binding , FTQUAL_old_locus_tag },

  { FEATDEF_misc_difference , FTQUAL_allele },
  { FEATDEF_misc_difference , FTQUAL_clone },
  { FEATDEF_misc_difference , FTQUAL_compare },
  { FEATDEF_misc_difference , FTQUAL_label },
  { FEATDEF_misc_difference , FTQUAL_map },
  { FEATDEF_misc_difference , FTQUAL_old_locus_tag },
  { FEATDEF_misc_difference , FTQUAL_phenotype },
  { FEATDEF_misc_difference , FTQUAL_replace },
  { FEATDEF_misc_difference , FTQUAL_standard_name },

  { FEATDEF_misc_feature , FTQUAL_allele },
  { FEATDEF_misc_feature , FTQUAL_function },
  { FEATDEF_misc_feature , FTQUAL_label },
  { FEATDEF_misc_feature , FTQUAL_map },
  { FEATDEF_misc_feature , FTQUAL_number },
  { FEATDEF_misc_feature , FTQUAL_old_locus_tag },
  { FEATDEF_misc_feature , FTQUAL_phenotype },
  { FEATDEF_misc_feature , FTQUAL_product },
  { FEATDEF_misc_feature , FTQUAL_standard_name },

  { FEATDEF_misc_recomb , FTQUAL_allele },
  { FEATDEF_misc_recomb , FTQUAL_label },
  { FEATDEF_misc_recomb , FTQUAL_map },
  { FEATDEF_misc_recomb , FTQUAL_old_locus_tag },
  { FEATDEF_misc_recomb , FTQUAL_standard_name },

  { FEATDEF_misc_signal , FTQUAL_allele },
  { FEATDEF_misc_signal , FTQUAL_function },
  { FEATDEF_misc_signal , FTQUAL_label },
  { FEATDEF_misc_signal , FTQUAL_map },
  { FEATDEF_misc_signal , FTQUAL_old_locus_tag },
  { FEATDEF_misc_signal , FTQUAL_operon },
  { FEATDEF_misc_signal , FTQUAL_phenotype },
  { FEATDEF_misc_signal , FTQUAL_standard_name },

  { FEATDEF_misc_structure , FTQUAL_allele },
  { FEATDEF_misc_structure , FTQUAL_function },
  { FEATDEF_misc_structure , FTQUAL_label },
  { FEATDEF_misc_structure , FTQUAL_map },
  { FEATDEF_misc_structure , FTQUAL_old_locus_tag },
  { FEATDEF_misc_structure , FTQUAL_standard_name },

  { FEATDEF_modified_base , FTQUAL_allele },
  { FEATDEF_modified_base , FTQUAL_frequency },
  { FEATDEF_modified_base , FTQUAL_label },
  { FEATDEF_modified_base , FTQUAL_map },
  { FEATDEF_modified_base , FTQUAL_mod_base },
  { FEATDEF_modified_base , FTQUAL_old_locus_tag },

  { FEATDEF_N_region , FTQUAL_allele },
  { FEATDEF_N_region , FTQUAL_label },
  { FEATDEF_N_region , FTQUAL_map },
  { FEATDEF_N_region , FTQUAL_old_locus_tag },
  { FEATDEF_N_region , FTQUAL_product },
  { FEATDEF_N_region , FTQUAL_standard_name },

  { FEATDEF_old_sequence , FTQUAL_allele },
  { FEATDEF_old_sequence , FTQUAL_compare },
  { FEATDEF_old_sequence , FTQUAL_label },
  { FEATDEF_old_sequence , FTQUAL_map },
  { FEATDEF_old_sequence , FTQUAL_old_locus_tag },
  { FEATDEF_old_sequence , FTQUAL_replace },

  { FEATDEF_polyA_signal , FTQUAL_allele },
  { FEATDEF_polyA_signal , FTQUAL_label },
  { FEATDEF_polyA_signal , FTQUAL_map },
  { FEATDEF_polyA_signal , FTQUAL_old_locus_tag },

  { FEATDEF_polyA_site , FTQUAL_allele },
  { FEATDEF_polyA_site , FTQUAL_label },
  { FEATDEF_polyA_site , FTQUAL_map },
  { FEATDEF_polyA_site , FTQUAL_old_locus_tag },

  { FEATDEF_prim_transcript , FTQUAL_allele },
  { FEATDEF_prim_transcript , FTQUAL_function },
  { FEATDEF_prim_transcript , FTQUAL_label },
  { FEATDEF_prim_transcript , FTQUAL_map },
  { FEATDEF_prim_transcript , FTQUAL_old_locus_tag },
  { FEATDEF_prim_transcript , FTQUAL_operon },
  { FEATDEF_prim_transcript , FTQUAL_standard_name },

  { FEATDEF_primer_bind , FTQUAL_allele },
  { FEATDEF_primer_bind , FTQUAL_label },
  { FEATDEF_primer_bind , FTQUAL_map },
  { FEATDEF_primer_bind , FTQUAL_old_locus_tag },
  { FEATDEF_primer_bind , FTQUAL_PCR_conditions },
  { FEATDEF_primer_bind , FTQUAL_standard_name },

  { FEATDEF_promoter , FTQUAL_allele },
  { FEATDEF_promoter , FTQUAL_bound_moiety },
  { FEATDEF_promoter , FTQUAL_function },
  { FEATDEF_promoter , FTQUAL_label },
  { FEATDEF_promoter , FTQUAL_map },
  { FEATDEF_promoter , FTQUAL_old_locus_tag },
  { FEATDEF_promoter , FTQUAL_operon },
  { FEATDEF_promoter , FTQUAL_phenotype },
  { FEATDEF_promoter , FTQUAL_standard_name },

  { FEATDEF_protein_bind , FTQUAL_allele },
  { FEATDEF_protein_bind , FTQUAL_bound_moiety },
  { FEATDEF_protein_bind , FTQUAL_function },
  { FEATDEF_protein_bind , FTQUAL_label },
  { FEATDEF_protein_bind , FTQUAL_map },
  { FEATDEF_protein_bind , FTQUAL_old_locus_tag },
  { FEATDEF_protein_bind , FTQUAL_operon },
  { FEATDEF_protein_bind , FTQUAL_standard_name },

  { FEATDEF_RBS , FTQUAL_allele },
  { FEATDEF_RBS , FTQUAL_label },
  { FEATDEF_RBS , FTQUAL_map },
  { FEATDEF_RBS , FTQUAL_old_locus_tag },
  { FEATDEF_RBS , FTQUAL_standard_name },

  { FEATDEF_repeat_region , FTQUAL_allele },
  { FEATDEF_repeat_region , FTQUAL_function },
  { FEATDEF_repeat_region , FTQUAL_label },
  { FEATDEF_repeat_region , FTQUAL_map },
  { FEATDEF_repeat_region , FTQUAL_mobile_element },
  { FEATDEF_repeat_region , FTQUAL_old_locus_tag },
  { FEATDEF_repeat_region , FTQUAL_rpt_family },
  { FEATDEF_repeat_region , FTQUAL_rpt_type },
  { FEATDEF_repeat_region , FTQUAL_rpt_unit },
  { FEATDEF_repeat_region , FTQUAL_rpt_unit_range },
  { FEATDEF_repeat_region , FTQUAL_rpt_unit_seq },
  { FEATDEF_repeat_region , FTQUAL_satellite },
  { FEATDEF_repeat_region , FTQUAL_standard_name },

  { FEATDEF_repeat_unit , FTQUAL_allele },
  { FEATDEF_repeat_unit , FTQUAL_function },
  { FEATDEF_repeat_unit , FTQUAL_label },
  { FEATDEF_repeat_unit , FTQUAL_map },
  { FEATDEF_repeat_unit , FTQUAL_old_locus_tag },
  { FEATDEF_repeat_unit , FTQUAL_rpt_family },
  { FEATDEF_repeat_unit , FTQUAL_rpt_type },
  { FEATDEF_repeat_unit , FTQUAL_rpt_unit },
  { FEATDEF_repeat_unit , FTQUAL_rpt_unit_range },
  { FEATDEF_repeat_unit , FTQUAL_rpt_unit_seq },

  { FEATDEF_rep_origin , FTQUAL_allele },
  { FEATDEF_rep_origin , FTQUAL_direction },
  { FEATDEF_rep_origin , FTQUAL_label },
  { FEATDEF_rep_origin , FTQUAL_map },
  { FEATDEF_rep_origin , FTQUAL_old_locus_tag },
  { FEATDEF_rep_origin , FTQUAL_standard_name },

  { FEATDEF_S_region , FTQUAL_allele },
  { FEATDEF_S_region , FTQUAL_label },
  { FEATDEF_S_region , FTQUAL_map },
  { FEATDEF_S_region , FTQUAL_old_locus_tag },
  { FEATDEF_S_region , FTQUAL_product },
  { FEATDEF_S_region , FTQUAL_standard_name },

  { FEATDEF_satellite , FTQUAL_allele },
  { FEATDEF_satellite , FTQUAL_label },
  { FEATDEF_satellite , FTQUAL_map },
  { FEATDEF_satellite , FTQUAL_old_locus_tag },
  { FEATDEF_satellite , FTQUAL_rpt_family },
  { FEATDEF_satellite , FTQUAL_rpt_type },
  { FEATDEF_satellite , FTQUAL_rpt_unit },
  { FEATDEF_satellite , FTQUAL_rpt_unit_range },
  { FEATDEF_satellite , FTQUAL_rpt_unit_seq },
  { FEATDEF_satellite , FTQUAL_standard_name },

  { FEATDEF_sig_peptide , FTQUAL_allele },
  { FEATDEF_sig_peptide , FTQUAL_function },
  { FEATDEF_sig_peptide , FTQUAL_label },
  { FEATDEF_sig_peptide , FTQUAL_map },
  { FEATDEF_sig_peptide , FTQUAL_old_locus_tag },
  { FEATDEF_sig_peptide , FTQUAL_product },
  { FEATDEF_sig_peptide , FTQUAL_standard_name },

  { FEATDEF_stem_loop , FTQUAL_allele },
  { FEATDEF_stem_loop , FTQUAL_function },
  { FEATDEF_stem_loop , FTQUAL_label },
  { FEATDEF_stem_loop , FTQUAL_map },
  { FEATDEF_stem_loop , FTQUAL_old_locus_tag },
  { FEATDEF_stem_loop , FTQUAL_operon },
  { FEATDEF_stem_loop , FTQUAL_standard_name },

  { FEATDEF_STS , FTQUAL_allele },
  { FEATDEF_STS , FTQUAL_label },
  { FEATDEF_STS , FTQUAL_map },
  { FEATDEF_STS , FTQUAL_old_locus_tag },
  { FEATDEF_STS , FTQUAL_standard_name },

  { FEATDEF_TATA_signal , FTQUAL_allele },
  { FEATDEF_TATA_signal , FTQUAL_label },
  { FEATDEF_TATA_signal , FTQUAL_map },
  { FEATDEF_TATA_signal , FTQUAL_old_locus_tag },

  { FEATDEF_terminator , FTQUAL_allele },
  { FEATDEF_terminator , FTQUAL_label },
  { FEATDEF_terminator , FTQUAL_map },
  { FEATDEF_terminator , FTQUAL_old_locus_tag },
  { FEATDEF_terminator , FTQUAL_operon },
  { FEATDEF_terminator , FTQUAL_standard_name },

  { FEATDEF_transit_peptide , FTQUAL_allele },
  { FEATDEF_transit_peptide , FTQUAL_function },
  { FEATDEF_transit_peptide , FTQUAL_label },
  { FEATDEF_transit_peptide , FTQUAL_map },
  { FEATDEF_transit_peptide , FTQUAL_old_locus_tag },
  { FEATDEF_transit_peptide , FTQUAL_product },
  { FEATDEF_transit_peptide , FTQUAL_standard_name },

  { FEATDEF_unsure , FTQUAL_allele },
  { FEATDEF_unsure , FTQUAL_compare },
  { FEATDEF_unsure , FTQUAL_label },
  { FEATDEF_unsure , FTQUAL_map },
  { FEATDEF_unsure , FTQUAL_old_locus_tag },
  { FEATDEF_unsure , FTQUAL_replace },

  { FEATDEF_V_region , FTQUAL_allele },
  { FEATDEF_V_region , FTQUAL_label },
  { FEATDEF_V_region , FTQUAL_map },
  { FEATDEF_V_region , FTQUAL_old_locus_tag },
  { FEATDEF_V_region , FTQUAL_product },
  { FEATDEF_V_region , FTQUAL_standard_name },

  { FEATDEF_V_segment , FTQUAL_allele },
  { FEATDEF_V_segment , FTQUAL_label },
  { FEATDEF_V_segment , FTQUAL_map },
  { FEATDEF_V_segment , FTQUAL_old_locus_tag },
  { FEATDEF_V_segment , FTQUAL_product },
  { FEATDEF_V_segment , FTQUAL_standard_name },

  { FEATDEF_variation , FTQUAL_allele },
  { FEATDEF_variation , FTQUAL_compare },
  { FEATDEF_variation , FTQUAL_frequency },
  { FEATDEF_variation , FTQUAL_label },
  { FEATDEF_variation , FTQUAL_map },
  { FEATDEF_variation , FTQUAL_old_locus_tag },
  { FEATDEF_variation , FTQUAL_phenotype },
  { FEATDEF_variation , FTQUAL_product },
  { FEATDEF_variation , FTQUAL_replace },
  { FEATDEF_variation , FTQUAL_standard_name },

  { FEATDEF_3clip , FTQUAL_allele },
  { FEATDEF_3clip , FTQUAL_function },
  { FEATDEF_3clip , FTQUAL_label },
  { FEATDEF_3clip , FTQUAL_map },
  { FEATDEF_3clip , FTQUAL_old_locus_tag },
  { FEATDEF_3clip , FTQUAL_standard_name },

  { FEATDEF_3UTR , FTQUAL_allele },
  { FEATDEF_3UTR , FTQUAL_function },
  { FEATDEF_3UTR , FTQUAL_label },
  { FEATDEF_3UTR , FTQUAL_map },
  { FEATDEF_3UTR , FTQUAL_old_locus_tag },
  { FEATDEF_3UTR , FTQUAL_standard_name },

  { FEATDEF_5clip , FTQUAL_allele },
  { FEATDEF_5clip , FTQUAL_function },
  { FEATDEF_5clip , FTQUAL_label },
  { FEATDEF_5clip , FTQUAL_map },
  { FEATDEF_5clip , FTQUAL_old_locus_tag },
  { FEATDEF_5clip , FTQUAL_standard_name },

  { FEATDEF_5UTR , FTQUAL_allele },
  { FEATDEF_5UTR , FTQUAL_function },
  { FEATDEF_5UTR , FTQUAL_label },
  { FEATDEF_5UTR , FTQUAL_map },
  { FEATDEF_5UTR , FTQUAL_old_locus_tag },
  { FEATDEF_5UTR , FTQUAL_standard_name },

  { FEATDEF_10_signal , FTQUAL_allele },
  { FEATDEF_10_signal , FTQUAL_label },
  { FEATDEF_10_signal , FTQUAL_map },
  { FEATDEF_10_signal , FTQUAL_old_locus_tag },
  { FEATDEF_10_signal , FTQUAL_operon },
  { FEATDEF_10_signal , FTQUAL_standard_name },

  { FEATDEF_35_signal , FTQUAL_allele },
  { FEATDEF_35_signal , FTQUAL_label },
  { FEATDEF_35_signal , FTQUAL_map },
  { FEATDEF_35_signal , FTQUAL_old_locus_tag },
  { FEATDEF_35_signal , FTQUAL_operon },
  { FEATDEF_35_signal , FTQUAL_standard_name },

  { FEATDEF_REGION , FTQUAL_function },
  { FEATDEF_REGION , FTQUAL_label },
  { FEATDEF_REGION , FTQUAL_map },
  { FEATDEF_REGION , FTQUAL_number },
  { FEATDEF_REGION , FTQUAL_old_locus_tag },
  { FEATDEF_REGION , FTQUAL_phenotype },
  { FEATDEF_REGION , FTQUAL_product },
  { FEATDEF_REGION , FTQUAL_standard_name },

  { FEATDEF_preprotein , FTQUAL_allele },
  { FEATDEF_preprotein , FTQUAL_label },
  { FEATDEF_preprotein , FTQUAL_map },
  { FEATDEF_preprotein , FTQUAL_old_locus_tag },
  { FEATDEF_preprotein , FTQUAL_product },
  { FEATDEF_preprotein , FTQUAL_standard_name },

  { FEATDEF_mat_peptide_aa , FTQUAL_allele },
  { FEATDEF_mat_peptide_aa , FTQUAL_label },
  { FEATDEF_mat_peptide_aa , FTQUAL_map },
  { FEATDEF_mat_peptide_aa , FTQUAL_old_locus_tag },
  { FEATDEF_mat_peptide_aa , FTQUAL_product },
  { FEATDEF_mat_peptide_aa , FTQUAL_standard_name },

  { FEATDEF_sig_peptide_aa , FTQUAL_allele },
  { FEATDEF_sig_peptide_aa , FTQUAL_label },
  { FEATDEF_sig_peptide_aa , FTQUAL_map },
  { FEATDEF_sig_peptide_aa , FTQUAL_old_locus_tag },
  { FEATDEF_sig_peptide_aa , FTQUAL_product },
  { FEATDEF_sig_peptide_aa , FTQUAL_standard_name },

  { FEATDEF_transit_peptide_aa , FTQUAL_allele },
  { FEATDEF_transit_peptide_aa , FTQUAL_label },
  { FEATDEF_transit_peptide_aa , FTQUAL_map },
  { FEATDEF_transit_peptide_aa , FTQUAL_old_locus_tag },
  { FEATDEF_transit_peptide_aa , FTQUAL_product },
  { FEATDEF_transit_peptide_aa , FTQUAL_standard_name },

  { FEATDEF_snoRNA , FTQUAL_allele },
  { FEATDEF_snoRNA , FTQUAL_function },
  { FEATDEF_snoRNA , FTQUAL_label },
  { FEATDEF_snoRNA , FTQUAL_map },
  { FEATDEF_snoRNA , FTQUAL_old_locus_tag },
  { FEATDEF_snoRNA , FTQUAL_product },
  { FEATDEF_snoRNA , FTQUAL_standard_name },

  { FEATDEF_gap , FTQUAL_estimated_length },
  { FEATDEF_gap , FTQUAL_gap_type },
  { FEATDEF_gap , FTQUAL_linkage_evidence },
  { FEATDEF_gap , FTQUAL_map },

  { FEATDEF_operon , FTQUAL_allele },
  { FEATDEF_operon , FTQUAL_function },
  { FEATDEF_operon , FTQUAL_label },
  { FEATDEF_operon , FTQUAL_map },
  { FEATDEF_operon , FTQUAL_operon },
  { FEATDEF_operon , FTQUAL_phenotype },
  { FEATDEF_operon , FTQUAL_standard_name },

  { FEATDEF_oriT , FTQUAL_allele },
  { FEATDEF_oriT , FTQUAL_bound_moiety },
  { FEATDEF_oriT , FTQUAL_direction },
  { FEATDEF_oriT , FTQUAL_label },
  { FEATDEF_oriT , FTQUAL_map },
  { FEATDEF_oriT , FTQUAL_old_locus_tag },
  { FEATDEF_oriT , FTQUAL_rpt_type },
  { FEATDEF_oriT , FTQUAL_rpt_unit },
  { FEATDEF_oriT , FTQUAL_rpt_unit_range },
  { FEATDEF_oriT , FTQUAL_rpt_unit_seq },
  { FEATDEF_oriT , FTQUAL_standard_name },

  { FEATDEF_ncRNA , FTQUAL_allele },
  { FEATDEF_ncRNA , FTQUAL_function },
  { FEATDEF_ncRNA , FTQUAL_label },
  { FEATDEF_ncRNA , FTQUAL_map },
  { FEATDEF_ncRNA , FTQUAL_ncRNA_class },
  { FEATDEF_ncRNA , FTQUAL_old_locus_tag },
  { FEATDEF_ncRNA , FTQUAL_operon },
  { FEATDEF_ncRNA , FTQUAL_product },
  { FEATDEF_ncRNA , FTQUAL_standard_name },

  { FEATDEF_tmRNA , FTQUAL_allele },
  { FEATDEF_tmRNA , FTQUAL_function },
  { FEATDEF_tmRNA , FTQUAL_label },
  { FEATDEF_tmRNA , FTQUAL_map },
  { FEATDEF_tmRNA , FTQUAL_old_locus_tag },
  { FEATDEF_tmRNA , FTQUAL_operon },
  { FEATDEF_tmRNA , FTQUAL_product },
  { FEATDEF_tmRNA , FTQUAL_standard_name },
  { FEATDEF_tmRNA , FTQUAL_tag_peptide },

  { FEATDEF_VARIATIONREF , FTQUAL_allele },
  { FEATDEF_VARIATIONREF , FTQUAL_compare },
  { FEATDEF_VARIATIONREF , FTQUAL_frequency },
  { FEATDEF_VARIATIONREF , FTQUAL_label },
  { FEATDEF_VARIATIONREF , FTQUAL_map },
  { FEATDEF_VARIATIONREF , FTQUAL_old_locus_tag },
  { FEATDEF_VARIATIONREF , FTQUAL_phenotype },
  { FEATDEF_VARIATIONREF , FTQUAL_product },
  { FEATDEF_VARIATIONREF , FTQUAL_replace },
  { FEATDEF_VARIATIONREF , FTQUAL_standard_name },

  { FEATDEF_mobile_element , FTQUAL_allele },
  { FEATDEF_mobile_element , FTQUAL_function },
  { FEATDEF_mobile_element , FTQUAL_label },
  { FEATDEF_mobile_element , FTQUAL_map },
  { FEATDEF_mobile_element , FTQUAL_mobile_element_type },
  { FEATDEF_mobile_element , FTQUAL_old_locus_tag },
  { FEATDEF_mobile_element , FTQUAL_rpt_family },
  { FEATDEF_mobile_element , FTQUAL_rpt_type },
  { FEATDEF_mobile_element , FTQUAL_standard_name }
};

/* comparison of ValQual's -- first compare featdef then ftqual */

/* macro did not work properly on linux machine, so using function instead */
/* #define COMPARE_VALQUAL(av,aq,bv,bq) ( ((av)-(bv)) ? ((av)-(bv)) : ((aq)-(bq)) ) */

static Int2 CompareValQual (Uint2 av, FtQualType aq, Uint2 bv, FtQualType bq)

{
  if (av == bv) return (aq - bq);
  return (av - bv);
}

/* Returns TRUE if {featureKey, qualKey} exists in legalGbqualList */

static Boolean AllowedValQual (Uint2 featureKey, FtQualType qualKey, Boolean forGbRelease)

{
  Int2 L, R, mid;

  if (qualKey == FTQUAL_experiment || qualKey == FTQUAL_inference) return TRUE;

  L = 0;
  R = sizeof (legalGbqualList) / sizeof (ValQual) - 1;
  while (L < R) {
    mid = (L + R) / 2;
    if (CompareValQual (legalGbqualList [mid].featdef,
       legalGbqualList [mid].ftqual,
       featureKey, qualKey) < 0)
      L = mid + 1;
    else
      R = mid;
  }
  if (CompareValQual (legalGbqualList [R].featdef,
      legalGbqualList [R].ftqual, featureKey, qualKey) == 0) {
    return TRUE;
  }

  return FALSE;
}


static CharPtr validRptString [] = {
  "tandem", "inverted", "flanking", "terminal", "direct", "dispersed", "other", NULL
};

static CharPtr validLRBString [] = {
  "LEFT", "RIGHT", "BOTH", NULL
};

static CharPtr validConsSpliceString [] = {
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

static Boolean StringInStringList (CharPtr testString, CharPtr PNTR stringList) {
  Int2 i;
  i = 0;
  while (stringList [i] != NULL) {
    if (StringICmp (testString, stringList [i]) == 0)
      return 1;
    i++;
  }
  return 0;
}

static CharPtr validMobileElementString [] = {
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

static Boolean ValidateMobileElement (CharPtr testString)

{
  Boolean  found;
  Int2     i;
  size_t   len;
  CharPtr  ptr, str;

  found = FALSE;
  str = NULL;
  for (i = 0; validMobileElementString [i] != NULL; i++) {
    ptr = validMobileElementString [i];
    len = StringLen (ptr);
    if (StringNICmp (testString, ptr, len) == 0) {
      found = TRUE;
      str = testString + len;
      break;
    }
  }
  if (found) {
    if (StringDoesHaveText (str) && (str [0] != ':' || str [1] == '\0')) {
      return FALSE;
    } else if (StringNICmp (testString, "other", 5) == 0) {
      if (str [0] != ':' || str [1] == '\0') {
        return FALSE;
      }
    }
  }
  return found;
}

/*
Functions now public and prototyped in sequtil.h
Return values are:
 0: no problem - Accession is in proper format
-1: Accession did not start with a letter (or two letters)
-2: Accession did not contain five numbers (or six numbers after 2 letters)
-3: the original Accession number to be validated was NULL
-4: the original Accession number is too long (>16)
-5: missing version number (required by ValidateAccnDotVer)
-6: bad version number (required by ValidateAccnDotVer)
*/

static Int2 ValidateAccnInternal (
  CharPtr accession,
  CharPtr PNTR strptr
)

{
  Char     ch;
  Int2     numAlpha = 0;
  Int2     numDigits = 0;
  Int2     numUndersc = 0;
  CharPtr  str;

  if (accession == NULL || accession [0] == '\0') return -3;

  if (StringLen (accession) >= 16) return -4;

  if (accession [0] < 'A' || accession [0] > 'Z') return -1;

  str = accession;
  if (StringNCmp (str, "NZ_", 3) == 0) {
    str += 3;
  }
  ch = *str;
  while (IS_ALPHA (ch)) {
    numAlpha++;
    str++;
    ch = *str;
  }
  while (ch == '_') {
    numUndersc++;
    str++;
    ch = *str;
  }
  while (IS_DIGIT (ch)) {
    numDigits++;
    str++;
    ch = *str;
  }
  if (ch != '\0' && ch != ' ' && ch != '.') return -2;

  if (numUndersc > 1) return -2;

  if (strptr != NULL) {
    /* pass back current position for version check */
    *strptr = str;
  }

  if (numUndersc == 0) {
    if (numAlpha == 1 && numDigits == 5) return 0;
    if (numAlpha == 2 && numDigits == 6) return 0;
    if (numAlpha == 3 && numDigits == 5) return 0;
    if (numAlpha == 4 && numDigits == 8) return 0;
    if (numAlpha == 4 && numDigits == 9) return 0;
    if (numAlpha == 5 && numDigits == 7) return 0;
  } else if (numUndersc == 1) {
    if (numAlpha != 2 || (numDigits != 6 && numDigits != 8 && numDigits != 9)) return -2;
    if (accession [0] == 'N' || accession [0] == 'X' || accession [0] == 'Z') {
      if (accession [1] == 'M' ||
          accession [1] == 'C' ||
          accession [1] == 'T' ||
          accession [1] == 'P' ||
          accession [1] == 'G' ||
          accession [1] == 'R' ||
          accession [1] == 'S' ||
          accession [1] == 'W' ||
          accession [1] == 'Z') {
        return 0;
      }
    }
    if (accession [0] == 'A' || accession [0] == 'Y') {
      if (accession [1] == 'P') return 0;
    }
  }

  return -2;
}

NLM_EXTERN Int2 ValidateAccn (
  CharPtr accession
)

{
  return ValidateAccnInternal (accession, NULL);
}

NLM_EXTERN Int2 ValidateAccnDotVer (
  CharPtr accession
)

{
  Char     ch;
  Int2     numVersion = 0;
  Int2     rsult;
  CharPtr  str = NULL;

  rsult = ValidateAccnInternal (accession, &str);
  if (rsult != 0) return rsult;

  if (str == NULL) return -5;
  ch = *str;
  if (ch != '.') return -5;
  str++;
  ch = *str;
  while (IS_DIGIT (ch)) {
    numVersion++;
    str++;
    ch = *str;
  }
  if (numVersion < 1) return -5;
  if (ch != '\0' && ch != ' ') return -6;

  return 0;
}

NLM_EXTERN Int2 ValidateSeqID (
  SeqIdPtr sip
)

{
  Char  buf [41];

  if (sip == NULL) return -3;
  SeqIdWrite (sip, buf, PRINTID_TEXTID_ACC_VER, sizeof (buf) - 1);
  return ValidateAccn (buf);
}

static Boolean ValidateCompareQual (CharPtr accession, Boolean is_ged)

{
  if (ValidateAccnDotVer (accession) != 0) return FALSE;
  if (StringChr (accession, '_') == NULL) return TRUE;
  if (is_ged) return FALSE;
  return TRUE;
}

static CharPtr mrnaevtext1 = "Derived by automated computational analysis";
static CharPtr mrnaevtext2 = "using gene prediction method:";
static CharPtr mrnaevtext3 = "Supporting evidence includes similarity to:";

static void GetStrFormRNAEvidence (
  UserObjectPtr uop,
  Pointer userdata
)

{
  Int2          ce = 0, cm = 0, cp = 0, ne = 0, nm = 0, np = 0;
  Boolean       has_counts = FALSE;
  size_t        len;
  CharPtr       method = NULL, prefix = NULL;
  ObjectIdPtr   oip;
  CharPtr       str = NULL;
  CharPtr PNTR  strp;
  Char          tmp [20];
  UserFieldPtr  u, ufp, uu;

  if (uop == NULL) return;
  oip = uop->type;
  if (oip == NULL) return;
  if (StringCmp (oip->str, "ModelEvidence") != 0) return;
  strp = (CharPtr PNTR) userdata;

  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    oip = ufp->label;
    if (oip == NULL || ufp->data.ptrvalue == NULL) continue;
    if (StringCmp (oip->str, "Method") == 0) {
      method = StringSaveNoNull ((CharPtr) ufp->data.ptrvalue);
    } else if (StringCmp (oip->str, "mRNA") == 0) {
      for (u = (UserFieldPtr) ufp->data.ptrvalue; u != NULL; u = u->next) {
        if (u->data.ptrvalue == NULL) continue;
        for (uu = (UserFieldPtr) u->data.ptrvalue; uu != NULL; uu = uu->next) {
          oip = uu->label;
          if (oip == NULL) continue;
          if (StringCmp (oip->str, "accession") == 0) {
            nm++;
          }
        }
      }
    } else if (StringCmp (oip->str, "EST") == 0) {
      for (u = (UserFieldPtr) ufp->data.ptrvalue; u != NULL; u = u->next) {
        if (u->data.ptrvalue == NULL) continue;
        for (uu = (UserFieldPtr) u->data.ptrvalue; uu != NULL; uu = uu->next) {
          oip = uu->label;
          if (oip == NULL) continue;
          if (StringCmp (oip->str, "accession") == 0) {
            ne++;
          }
        }
      }
    } else if (StringCmp (oip->str, "Protein") == 0) {
      for (u = (UserFieldPtr) ufp->data.ptrvalue; u != NULL; u = u->next) {
        if (u->data.ptrvalue == NULL) continue;
        for (uu = (UserFieldPtr) u->data.ptrvalue; uu != NULL; uu = uu->next) {
          oip = uu->label;
          if (oip == NULL) continue;
          if (StringCmp (oip->str, "accession") == 0) {
            np++;
          }
        }
      }
    } else if (StringCmp (oip->str, "Counts") == 0) {
      has_counts = TRUE;
      for (u = (UserFieldPtr) ufp->data.ptrvalue; u != NULL; u = u->next) {
        if (u->data.ptrvalue == NULL) continue;
        if (u->choice != 2) continue;
        oip = u->label;
        if (oip == NULL) continue;
        if (StringCmp (oip->str, "mRNA") == 0) {
          cm = (Int2) u->data.intvalue;
        } else if (StringCmp (oip->str, "EST") == 0) {
          ce = (Int2) u->data.intvalue;
        } else if (StringCmp (oip->str, "Protein") == 0) {
          cp = (Int2) u->data.intvalue;
        }
      }
    }
  }

  if (has_counts) {
    nm = cm;
    ne = ce;
    np = cp;
  }

  len = StringLen (mrnaevtext1) + StringLen (mrnaevtext2) + StringLen (mrnaevtext3) + StringLen (method) + 80;
  str = (CharPtr) MemNew (len);
  if (str == NULL) return;

  if (method != NULL) {
    sprintf (str, "%s %s %s.", mrnaevtext1, mrnaevtext2, method);
  } else {
    sprintf (str, "%s.", mrnaevtext1);
  }
  if (nm > 0 || ne > 0 || np > 0) {
    StringCat (str, " ");
    StringCat (str, mrnaevtext3);
  }
  prefix = " ";
  if (nm > 0) {
    StringCat (str, prefix);
    if (nm > 1) {
      sprintf (tmp, "%d mRNAs", (int) nm);
    } else {
      sprintf (tmp, "%d mRNA", (int) nm);
    }
    StringCat (str, tmp);
    prefix = ", ";
  }
  if (ne > 0) {
    StringCat (str, prefix);
    if (ne > 1) {
      sprintf (tmp, "%d ESTs", (int) ne);
    } else {
      sprintf (tmp, "%d EST", (int) ne);
    }
    StringCat (str, tmp);
    prefix = ", ";
  }
  if (np > 0) {
    StringCat (str, prefix);
    if (np > 1) {
      sprintf (tmp, "%d Proteins", (int) np);
    } else {
      sprintf (tmp, "%d Protein", (int) np);
    }
    StringCat (str, tmp);
    prefix = ", ";
  }

  *strp = str;
}

static Boolean ValidateRptUnit (
  CharPtr buf
)

{
#if 0
  CharPtr  str;
  Char     tmp [255];

  StringNCpy_0 (tmp, buf, sizeof (tmp));
  TrimSpacesAroundString (tmp);

  str = tmp;
  /* first check for sequence letters with optional semicolons */
  while (IS_ALPHA (*str) || *str == ';') str++;
  if (*str == '\0') return TRUE;
  /* next check for letters, digits, commas, parentheses, dashes, and underscores */
  str = tmp;
  while (IS_ALPHANUM (*str) || *str == '(' || *str == ')' || *str == ',' || *str == ';' || *str == '-' || *str == '_') str++;
  if (*str == '\0') return TRUE;
  /* now check for officially legal styles */
  str = tmp;
  while (IS_ALPHANUM (*str)) str++;
  if (*str != '\0') { /* wasn't pure alphanumeric; now check for xxx..yyy */
    str = buf;
    while (IS_DIGIT (*str)) str++; /* xxx */
    if (*str == '\0' /* must be something after the xxx */
      || StringLen (str) < 3  /* need at least 2 '.'s and a digit*/
      || str[0] != '.' || str[1] != '.') return FALSE;
    str+=2;
    while (IS_DIGIT (*str)) str++;
    if (*str != '\0') return FALSE;  /* mustn't be anything after the yyy */
  }
#endif
  return TRUE;
}


NLM_EXTERN CharPtr goFieldType [] = {
  "", "text string", "go id", "pubmed id", "go ref", "evidence", NULL
};

typedef struct gostruc {
  CharPtr  term;
  CharPtr  goid;
  CharPtr  evidence;
  Int4     pmid;
  CharPtr  goref;
} GoStruc, PNTR GoStrucPtr;

static int LIBCALLBACK SortVnpByGsp (VoidPtr ptr1, VoidPtr ptr2)

{
  int           compare;
  GoStrucPtr    gsp1, gsp2;
  ValNodePtr    vnp1, vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  gsp1 = (GoStrucPtr) vnp1->data.ptrvalue;
  gsp2 = (GoStrucPtr) vnp2->data.ptrvalue;
  if (gsp1 == NULL || gsp2 == NULL) return 0;

  compare = StringICmp (gsp1->term, gsp2->term);
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

static CharPtr GetCombinedGOtext (
  UserFieldPtr entryhead,
  IntAsn2gbJobPtr ajp
)

{
  UserFieldPtr   entry, topufp, ufp;
  CharPtr        evidence, goid, goref, last = NULL,
                 str, textstr, prefix;
  StringItemPtr  ffstring;
  Char           gid [32], tmp [32];
  GoStrucPtr     gsp;
  ValNodePtr     head = NULL, vnp;
  Boolean        is_www;
  Int2           j;
  ObjectIdPtr    oip;
  Int4           pmid;

  if (entryhead == NULL || ajp == NULL) return NULL;
  is_www = GetWWW (ajp);

  for (entry = entryhead; entry != NULL; entry = entry->next) {
    if (entry == NULL || entry->choice != 11) break;
    topufp = (UserFieldPtr)  entry->data.ptrvalue;
    if (topufp == NULL) continue;

    textstr = NULL;
    evidence = NULL;
    goid = NULL;
    goref = NULL;
    pmid = 0;
    for (ufp = topufp; ufp != NULL; ufp = ufp->next) {
      oip = ufp->label;
      if (oip == NULL) continue;
      for (j = 0; goFieldType [j] != NULL; j++) {
        if (StringICmp (oip->str, goFieldType [j]) == 0) break;
      }
      if (goFieldType [j] == NULL) continue;
      switch (j) {
        case 1 :
          if (ufp->choice == 1) {
            textstr = (CharPtr) ufp->data.ptrvalue;
          }
          break;
        case 2 :
          if (ufp->choice == 1) {
            goid = (CharPtr) ufp->data.ptrvalue;
          } else if (ufp->choice == 2) {
            sprintf (gid, "%ld", (long) (Int4) ufp->data.intvalue);
            goid = (CharPtr) gid;
          }
          break;
        case 3 :
          if (ufp->choice == 2) {
            pmid = (Int4) ufp->data.intvalue;
          }
          break;
        case 4 :
          if (ufp->choice == 1) {
            goref = (CharPtr) ufp->data.ptrvalue;
          }
          break;
        case 5 :
          if (ufp->choice == 1) {
            evidence = (CharPtr) ufp->data.ptrvalue;
          }
          break;
        default :
          break;
      }
    }

    if (StringDoesHaveText (textstr)) {
      gsp = (GoStrucPtr) MemNew (sizeof (GoStruc));
      if (gsp != NULL) {
        gsp->term = StringSave (textstr);
        gsp->goid = StringSave (goid);
        gsp->evidence = StringSave (evidence);
        gsp->pmid = pmid;
        gsp->goref = StringSave (goref);
        ValNodeAddPointer (&head, 0, (Pointer) gsp);
      }
    }
  }

  if (head == NULL) return NULL;
  head = ValNodeSort (head, SortVnpByGsp);

  if (is_www) {
    ffstring = FFGetString (ajp);
    if (ffstring != NULL) {

      last = NULL;
      prefix = NULL;
      for (vnp = head; vnp != NULL; vnp = vnp->next) {
        gsp = (GoStrucPtr) vnp->data.ptrvalue;
        if (gsp == NULL) continue;
        if (StringICmp (gsp->term, last) != 0) {
          if (prefix != NULL) {
            FFAddOneString (ffstring, prefix, FALSE, FALSE, TILDE_IGNORE);
          }
          if (StringDoesHaveText (gsp->goid)) {
            FFAddOneString (ffstring, "GO:", FALSE, TRUE, TILDE_IGNORE);
            FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
            FF_Add_NCBI_Base_URL (ffstring, link_go);
            FFAddOneString (ffstring, gsp->goid, FALSE, FALSE, TILDE_IGNORE);
            FFAddOneString (ffstring, "\">", FALSE, FALSE, TILDE_IGNORE);
            FFAddOneString (ffstring, gsp->goid, FALSE, TRUE, TILDE_IGNORE);
            FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
          }
          if (StringDoesHaveText (gsp->term)) {
            FFAddOneString (ffstring, " - ", FALSE, TRUE, TILDE_IGNORE);
            FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
            FF_Add_NCBI_Base_URL (ffstring, link_go);
            FFAddOneString (ffstring, gsp->goid, FALSE, FALSE, TILDE_IGNORE);
            FFAddOneString (ffstring, "\">", FALSE, FALSE, TILDE_IGNORE);
            FFAddOneString (ffstring, gsp->term, FALSE, TRUE, TILDE_IGNORE);
            FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
          }
        }
        if (StringDoesHaveText (gsp->evidence)) {
          FFAddOneString (ffstring, " [Evidence ", FALSE, TRUE, TILDE_IGNORE);
          FFAddOneString (ffstring, gsp->evidence, FALSE, TRUE, TILDE_IGNORE);
          FFAddOneString (ffstring, "]", FALSE, TRUE, TILDE_IGNORE);
        }
        if (gsp->pmid > 0) {
          sprintf (tmp, "%ld", (long) gsp->pmid);
          FFAddOneString (ffstring, " [PMID <a href=\"", FALSE, FALSE, TILDE_IGNORE);
          FF_Add_NCBI_Base_URL (ffstring, link_muid);
          FFAddOneString (ffstring, tmp, FALSE, FALSE, TILDE_IGNORE);
          FFAddOneString (ffstring, "\">", FALSE, FALSE, TILDE_IGNORE);
          FFAddOneString (ffstring, tmp, FALSE, FALSE, TILDE_IGNORE);
          FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
          FFAddOneString (ffstring, "]", FALSE, TRUE, TILDE_IGNORE);
        } else if (StringDoesHaveText (gsp->goref)) {
          FFAddOneString (ffstring, " [GO Ref <a href=\"", FALSE, FALSE, TILDE_IGNORE);
          FF_Add_NCBI_Base_URL (ffstring, link_go_ref);
          FFAddOneString (ffstring, gsp->goref, FALSE, FALSE, TILDE_IGNORE);
          FFAddOneString (ffstring, "\">", FALSE, FALSE, TILDE_IGNORE);
          FFAddOneString (ffstring, gsp->goref, FALSE, FALSE, TILDE_IGNORE);
          FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
          FFAddOneString (ffstring, "]", FALSE, TRUE, TILDE_IGNORE);
        }
        prefix = "; ";
        last = gsp->term;
      }

      str = FFToCharPtr (ffstring);
      TrimSpacesAroundString (str);

      FFRecycleString (ajp, ffstring);

      for (vnp = head; vnp != NULL; vnp = vnp->next) {
        gsp = (GoStrucPtr) vnp->data.ptrvalue;
        if (gsp == NULL) continue;
        gsp->term = MemFree (gsp->term);
        gsp->goid = MemFree (gsp->goid);
        gsp->goref = MemFree (gsp->goref);
        gsp->evidence = MemFree (gsp->evidence);
      }
      ValNodeFreeData (head);

      return str;
    }
  }

  /* not is_www */

  ffstring = FFGetString (ajp);
  if (ffstring != NULL) {
    last = NULL;
    prefix = NULL;
    for (vnp = head; vnp != NULL; vnp = vnp->next) {
      gsp = (GoStrucPtr) vnp->data.ptrvalue;
      if (gsp == NULL) continue;
      if (StringICmp (gsp->term, last) != 0) {
        if (prefix != NULL) {
          FFAddOneString (ffstring, prefix, FALSE, FALSE, TILDE_IGNORE);
        }
        if (StringDoesHaveText (gsp->goid)) {
          FFAddOneString (ffstring, "GO:", FALSE, TRUE, TILDE_IGNORE);
          FFAddOneString (ffstring, gsp->goid, FALSE, TRUE, TILDE_IGNORE);
        }
        if (StringDoesHaveText (gsp->term)) {
          FFAddOneString (ffstring, " - ", FALSE, TRUE, TILDE_IGNORE);
          FFAddOneString (ffstring, gsp->term, FALSE, TRUE, TILDE_IGNORE);
        }
      }
      if (StringDoesHaveText (gsp->evidence)) {
        FFAddOneString (ffstring, " [Evidence ", FALSE, TRUE, TILDE_IGNORE);
        FFAddOneString (ffstring, gsp->evidence, FALSE, TRUE, TILDE_IGNORE);
        FFAddOneString (ffstring, "]", FALSE, TRUE, TILDE_IGNORE);
      }
      if (gsp->pmid > 0) {
        sprintf (tmp, "%ld", (long) gsp->pmid);
        FFAddOneString (ffstring, " [PMID ", FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, tmp, FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, "]", FALSE, TRUE, TILDE_IGNORE);
      } else if (StringDoesHaveText (gsp->goref)) {
        FFAddOneString (ffstring, " [GO Ref", FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, gsp->goref, FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, "]", FALSE, TRUE, TILDE_IGNORE);
      }
      prefix = "; ";
      last = gsp->term;
    }
  }

  str = FFToCharPtr (ffstring);
  TrimSpacesAroundString (str);

  FFRecycleString (ajp, ffstring);

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    gsp = (GoStrucPtr) vnp->data.ptrvalue;
    if (gsp == NULL) continue;
    gsp->term = MemFree (gsp->term);
    gsp->goid = MemFree (gsp->goid);
    gsp->goref = MemFree (gsp->goref);
    gsp->evidence = MemFree (gsp->evidence);
  }
  ValNodeFreeData (head);

  return str;
}

static CharPtr GetGOtext (
  UserFieldPtr topufp,
  IntAsn2gbJobPtr ajp,
  Boolean abbreviate
)

{
  CharPtr        evidence = NULL;
  StringItemPtr  ffstring;
  Char           gid [32];
  CharPtr        goid = NULL;
  CharPtr        goref = NULL;
  Boolean        is_www;
  Int2           j;
  ObjectIdPtr    oip;
  Int4           pmid = 0;
  CharPtr        str;
  CharPtr        textstr = NULL;
  Char           tmp [32];
  UserFieldPtr   ufp;

  if (topufp == NULL || ajp == NULL) return NULL;
  is_www = GetWWW (ajp);

  for (ufp = topufp; ufp != NULL; ufp = ufp->next) {
    oip = ufp->label;
    if (oip == NULL) continue;
    for (j = 0; goFieldType [j] != NULL; j++) {
      if (StringICmp (oip->str, goFieldType [j]) == 0) break;
    }
    if (goFieldType [j] == NULL) continue;
    switch (j) {
      case 1 :
        if (ufp->choice == 1) {
          textstr = (CharPtr) ufp->data.ptrvalue;
        }
        break;
      case 2 :
        if (ufp->choice == 1) {
          goid = (CharPtr) ufp->data.ptrvalue;
        } else if (ufp->choice == 2) {
          sprintf (gid, "%ld", (long) (Int4) ufp->data.intvalue);
          goid = (CharPtr) gid;
        }
        break;
      case 3 :
        if (ufp->choice == 2) {
          pmid = (Int4) ufp->data.intvalue;
        }
        break;
      case 4 :
        if (ufp->choice == 1) {
          goref = (CharPtr) ufp->data.ptrvalue;
        }
        break;
      case 5 :
        if (ufp->choice == 1) {
          evidence = (CharPtr) ufp->data.ptrvalue;
        }
        break;
      default :
        break;
    }
  }
  /* if (StringHasNoText (textstr)) return NULL; */

  if (is_www) {
    ffstring = FFGetString (ajp);
    if (ffstring != NULL) {
      if (StringDoesHaveText (goid)) {
        FFAddOneString (ffstring, "GO:", FALSE, TRUE, TILDE_IGNORE);
        FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
        FF_Add_NCBI_Base_URL (ffstring, link_go);
        FFAddOneString (ffstring, goid, FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, "\">", FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, goid, FALSE, TRUE, TILDE_IGNORE);
        FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
      }

      if (StringDoesHaveText (textstr)) {
        FFAddOneString (ffstring, " - ", FALSE, TRUE, TILDE_IGNORE);
        FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
        FF_Add_NCBI_Base_URL (ffstring, link_go);
        FFAddOneString (ffstring, goid, FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, "\">", FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, textstr, FALSE, TRUE, TILDE_IGNORE);
        FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
      }

      if (StringDoesHaveText (evidence)) {
        FFAddOneString (ffstring, " [Evidence ", FALSE, TRUE, TILDE_IGNORE);
        FFAddOneString (ffstring, evidence, FALSE, TRUE, TILDE_IGNORE);
        FFAddOneString (ffstring, "]", FALSE, TRUE, TILDE_IGNORE);
      }

      if (pmid != 0) {
        sprintf (tmp, "%ld", (long) pmid);
        FFAddOneString (ffstring, " [PMID ", FALSE, TRUE, TILDE_IGNORE);
        FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
        FF_Add_NCBI_Base_URL (ffstring, link_muid);
        FFAddOneString (ffstring, tmp, FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, "\">", FALSE, FALSE, TILDE_IGNORE); 
        FFAddOneString (ffstring, tmp, FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, "]", FALSE, TRUE, TILDE_IGNORE);
      } else if (StringDoesHaveText (goref)) {
        FFAddOneString (ffstring, " [GO Ref ", FALSE, TRUE, TILDE_IGNORE);
        FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
        FF_Add_NCBI_Base_URL (ffstring, link_go_ref);
        FFAddOneString (ffstring, goref, FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, "\">", FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, goref, FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, "]", FALSE, TRUE, TILDE_IGNORE);
      }

      str = FFToCharPtr (ffstring);
      TrimSpacesAroundString (str);

      FFRecycleString (ajp, ffstring);

      return str;
    }
  }

  /* not is_www */

  str = (CharPtr) MemNew (StringLen (goid) + StringLen (textstr) +
                          StringLen (evidence) + StringLen (goref) + 100);
  if (str == NULL) return NULL;

  if (StringDoesHaveText (goid)) {
    StringCat (str, "GO:");
    StringCat (str, goid);
  }

  if (StringDoesHaveText (textstr)) {
    StringCat (str, " - ");
    StringCat (str, textstr);
  }

  if (StringDoesHaveText (evidence)) {
    StringCat (str, " [Evidence ");
    StringCat (str, evidence);
    StringCat (str, "]");
  }

  if (pmid != 0) {
    sprintf (tmp, "%ld", (long) pmid);
    StringCat (str, " [PMID ");
    StringCat (str, tmp);
    StringCat (str, "]");
  } else if (StringDoesHaveText (goref)) {
    StringCat (str, " [GO Ref ");
    StringCat (str, goref);
    StringCat (str, "]");
  }

  TrimSpacesAroundString (str);

  return str;
}

static void GetNomenclatureText (
  UserObjectPtr uop,
  Pointer userdata
)

{
  CharPtr       ds = NULL, me = NULL, nm = NULL, sy = NULL;
  size_t        len;
  ObjectIdPtr   oip;
  CharPtr       str = NULL;
  CharPtr PNTR  strp;
  UserFieldPtr  ufp;

  if (uop == NULL) return;
  oip = uop->type;
  if (oip == NULL) return;
  if (StringCmp (oip->str, "OfficialNomenclature") != 0) return;
  strp = (CharPtr PNTR) userdata;

  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    oip = ufp->label;
    if (oip == NULL || oip->str == NULL) continue;
    if (StringICmp (oip->str, "Symbol") == 0) {
      if (ufp->choice == 1) {
        str = (CharPtr) ufp->data.ptrvalue;
        if (StringDoesHaveText (str)) {
          sy = str;
        }
      }
    } else if (StringICmp (oip->str, "Name") == 0) {
      if (ufp->choice == 1) {
        str = (CharPtr) ufp->data.ptrvalue;
        if (StringDoesHaveText (str)) {
          nm = str;
        }
      }
    } else if (StringICmp (oip->str, "DataSource") == 0) {
      if (ufp->choice == 1) {
        str = (CharPtr) ufp->data.ptrvalue;
        if (StringDoesHaveText (str)) {
          ds = str;
        }
      }
    } else if (StringICmp (oip->str, "Status") == 0) {
      if (ufp->choice == 1) {
        str = (CharPtr) ufp->data.ptrvalue;
        if (StringDoesHaveText (str)) {
          me = str;
        }
      }
    }
  }
  if (me == NULL) {
    me = "Unclassified";
  }

  if (StringHasNoText (sy)) return;

  len = StringLen (ds) + StringLen (me) + StringLen (nm) + StringLen (sy) + 80;
  str = (CharPtr) MemNew (len);
  if (str == NULL) return;

  StringCpy (str, me);
  StringCat (str, " Symbol: ");
  StringCat (str, sy);

  if (StringDoesHaveText (nm)) {
    StringCat (str, " | Name: ");
    StringCat (str, nm);
  }

  if (StringDoesHaveText (ds)) {
    StringCat (str, " | Provided by: ");
    StringCat (str, ds);
  }

  *strp = str;
}

static CharPtr GetNomenclature (
  GeneNomenclaturePtr gnp
)

{
  Char         buf [32];
  CharPtr      db = NULL, ds = NULL, me = NULL, nm = NULL, sy = NULL, str = NULL;
  DbtagPtr     dbt;
  size_t       len;
  ObjectIdPtr  oip;

  if (gnp == NULL) return NULL;

  if (StringDoesHaveText (gnp->symbol)) {
    sy = gnp->symbol;
  }
  if (StringHasNoText (sy)) return NULL;

  if (gnp->status == 1) {
    me = "Official";
  } else if (gnp->status == 2) {
    me = "Interim";
  }
  if (me == NULL) {
    me = "Unclassified";
  }

  if (StringDoesHaveText (gnp->name)) {
    nm = gnp->name;
  }

  dbt = gnp->source;
  if (dbt != NULL) {
    if (StringDoesHaveText (dbt->db)) {
      db = dbt->db;
    }
    oip = dbt->tag;
    if (oip != NULL) {
      if (StringDoesHaveText (oip->str)) {
        ds = oip->str;
      } else {
        sprintf (buf, "%ld", (long) oip->id);
        ds = buf;
      }
    }
  }

  len = StringLen (db) + StringLen (ds) + StringLen (me) + StringLen (nm) + StringLen (sy) + 80;
  str = (CharPtr) MemNew (sizeof (Char) * len);
  if (str == NULL) return NULL;

  StringCpy (str, me);
  StringCat (str, " Symbol: ");
  StringCat (str, sy);

  if (StringDoesHaveText (nm)) {
    StringCat (str, " | Name: ");
    StringCat (str, nm);
  }

  if (StringDoesHaveText (db) && StringDoesHaveText (ds)) {
    StringCat (str, " | Provided by: ");
    StringCat (str, db);
    StringCat (str, ":");
    StringCat (str, ds);
  }

  return str;
} 

static Boolean DbxrefAlreadyInGeneXref (
  DbtagPtr dbt,
  ValNodePtr dbxref
)

{
  DbtagPtr    gdbt;
  ValNodePtr  vnp;

  if (dbt == NULL) return FALSE;

  for (vnp = dbxref; vnp != NULL; vnp = vnp->next) {
    gdbt = (DbtagPtr) vnp->data.ptrvalue;
    if (gdbt == NULL) continue;
    if (DbtagMatch (dbt, gdbt)) return TRUE;
  }

  return FALSE;
}

static void FF_www_nuc_or_prot_id (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  CharPtr seqid,
  Int4 gi,
  Boolean is_na
)
{
  Char  buf [32];

  if ( GetWWW(ajp) ) {
    FFAddOneString(ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
    if (is_na) {
      FF_Add_NCBI_Base_URL (ffstring, link_seqn);
    } else {
      FF_Add_NCBI_Base_URL (ffstring, link_seqp);
    }
    if (gi > 0) {
      sprintf (buf, "%ld", (long) gi);
      FFAddOneString(ffstring, buf, FALSE, FALSE, TILDE_IGNORE);
    } else {
      FFAddOneString(ffstring, seqid, FALSE, FALSE, TILDE_IGNORE);
    }
    FFAddOneString(ffstring, "\">", FALSE, FALSE, TILDE_IGNORE); 
    FFAddOneString(ffstring, seqid, FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString(ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
  } else {
    FFAddOneString(ffstring, seqid, FALSE, FALSE, TILDE_IGNORE);
  }
}

static void FF_www_gcode (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  CharPtr gcode
)
{

  if ( GetWWW(ajp) ) {
    FFAddOneString(ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
    FF_Add_NCBI_Base_URL (ffstring, link_code);
    FFAddOneString(ffstring, "mode=c#SG", FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString(ffstring, gcode, FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString(ffstring, "\">", FALSE, FALSE, TILDE_IGNORE); 
    FFAddOneString(ffstring, gcode, FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString(ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
  } else {
    FFAddOneString(ffstring, gcode, FALSE, FALSE, TILDE_IGNORE);
  }
}

static void FF_AddECnumber (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  CharPtr str
)
{
  if (StringHasNoText (str)) return;
  if ( GetWWW(ajp) ) {
    /*
    if (StringChr (str, '-') != NULL) {
      FFAddOneString(ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
      FFAddOneString(ffstring, ec_ambig, FALSE, FALSE, TILDE_IGNORE);
      FFAddOneString(ffstring, "\">", FALSE, FALSE, TILDE_IGNORE);
      FFAddOneString(ffstring, str, FALSE, FALSE, TILDE_IGNORE);
      FFAddOneString(ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
    } else {
      FFAddOneString(ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
      FFAddOneString(ffstring, ec_link, FALSE, FALSE, TILDE_IGNORE);
      FFAddOneString(ffstring, str, FALSE, FALSE, TILDE_IGNORE);
      FFAddOneString(ffstring, "\">", FALSE, FALSE, TILDE_IGNORE); 
      FFAddOneString(ffstring, str, FALSE, FALSE, TILDE_IGNORE);
      FFAddOneString(ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
    }
    */
    FFAddOneString(ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString(ffstring, ec_link, FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString(ffstring, str, FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString(ffstring, "\">", FALSE, FALSE, TILDE_IGNORE); 
    FFAddOneString(ffstring, str, FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString(ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
  } else {
    FFAddOneString(ffstring, str, FALSE, FALSE, TILDE_IGNORE);
  }
}


/* FormatFeatureblockQuals should not be called directly,
   except from FormatFeatureBlock.  It performs no input
   validation.  (perhaps it should?) */

static void LIBCALLBACK SaveGBSeqTranslation (
  CharPtr sequence,
  Pointer userdata
)

{
  CharPtr       tmp;
  CharPtr PNTR  tmpp;

  tmpp = (CharPtr PNTR) userdata;
  tmp = *tmpp;

  tmp = StringMove (tmp, sequence);

  *tmpp = tmp;
}

static int LIBCALLBACK SortVnpByInt (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;

  if (vnp1->data.intvalue > vnp2->data.intvalue) {
    return 1;
  } else if (vnp1->data.intvalue < vnp2->data.intvalue) {
    return -1;
  }

  return 0;
}

static FloatHi MolWtForProtFeat (
  BioseqPtr bsp,
  SeqFeatPtr sfp,
  IntPrtBlockPtr ipp
)

{
  size_t      len;
  FloatHi     mol_wt = 0.0;
  ProtRefPtr  prp;
  CharPtr     str;

  if (bsp == NULL || sfp == NULL || ipp == NULL) return 0.0;
  if (sfp->partial) return 0.0;
  prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  if (prp == NULL) return 0.0;

  if (prp->processed >= 2) {
    return MolWtForLoc (sfp->location);
  }

  if (! ipp->is_whole_loc) {
    return MolWtForLoc (sfp->location);
  }

  if (ipp->suppress_mol_wt) return 0.0;

  if (ipp->sig_pept_trim_len > 0) {
    str = GetSequenceByFeature (sfp);
    if (str == NULL) return 0.0;
    len = StringLen (str);
    if (len > ipp->sig_pept_trim_len) {
      mol_wt = MolWtForStr (str + ipp->sig_pept_trim_len);
    } else {
      mol_wt = MolWtForStr (str);
    }
    MemFree (str);
    return mol_wt;
  }

  if (ipp->trim_initial_met) {
    str = GetSequenceByFeature (sfp);
    if (str == NULL) return 0.0;
    if (StringLen (str) > 1 && *str == 'M') {
      mol_wt = MolWtForStr (str + 1);
    } else {
      mol_wt = MolWtForStr (str);
    }
    MemFree (str);
    return mol_wt;
  }

  return MolWtForLoc (sfp->location);
}

/*
static void ChangeOandJtoX (CharPtr str)

{
  Char  ch;

  if (str == NULL) return;
  ch = *str;
  while (ch != '\0') {
    if (ch == 'O' || ch == 'J') {
      *str = 'X';
    } else if (ch == 'o' || ch == 'j') {
      *str = 'x';
    }
    str++;
    ch = *str;
  }
}
*/

static Boolean OnlyOneRealGeneral (SeqIdPtr sip)

{
  DbtagPtr  dbt;
  Int2      numGenerals = 0;

  while (sip != NULL) {
    if (sip->choice != SEQID_GENERAL) return FALSE;
    dbt = (DbtagPtr) sip->data.ptrvalue;
    if (dbt == NULL) return FALSE;
    if (!IsSkippableDbtag(dbt) &&
        StringICmp (dbt->db, "SMART") != 0) {
      numGenerals++;
    }
    sip = sip->next;
  }
  if (numGenerals == 1) return TRUE;
  return FALSE;
}

static void AddExperimentWithPMIDlinks(
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  CharPtr str
)

{
  Char     ch;
  Boolean  had_pmid;
  CharPtr  pmid;
  CharPtr  prefix = "PMID:";
  CharPtr  ptr;

  if (! GetWWW (ajp)) {
    FFAddOneString (ffstring, str, FALSE, TRUE, TILDE_IGNORE);
    return;
  }

  if (CommentHasSuspiciousHtml (ajp, str)) {
    FFAddOneString (ffstring, str, FALSE, TRUE, TILDE_IGNORE);
    return;
  }

  while (StringDoesHaveText (str)) {
    ptr = StringStr (str, prefix);
    if (ptr == NULL) {
      FFAddOneString (ffstring, str, FALSE, TRUE, TILDE_IGNORE);
      return;
    }
    *ptr = '\0';
    FFAddOneString (ffstring, str, FALSE, TRUE, TILDE_IGNORE);
    ptr += StringLen (prefix);
    pmid = ptr;
    ch = *ptr;
    while (ch == ' ') {
      ptr++;
      pmid = ptr;
      ch = *ptr;
    }
    while (IS_DIGIT (ch)) {
      ptr++;
      ch = *ptr;
    }
    *ptr = '\0';

    had_pmid = FALSE;
    if (StringDoesHaveText (pmid)) {
      FFAddOneString (ffstring, prefix, FALSE, TRUE, TILDE_IGNORE);
      FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
      FF_Add_NCBI_Base_URL (ffstring, link_muid);
      FFAddTextToString (ffstring, NULL, pmid, "\">", FALSE, FALSE, TILDE_IGNORE);
      FFAddOneString (ffstring, pmid, FALSE, FALSE, TILDE_IGNORE);
      FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
      had_pmid = TRUE;
    }

    *ptr = ch;
    str = ptr;

    prefix = "PMID:";
    ptr = str;
    ch = *ptr;
    if (had_pmid) {
      if (ch == ',') {
        ptr++;
        ch = *ptr;
        while (ch == ' ') {
          ptr++;
          ch = *ptr;
        }
        if (IS_DIGIT (ch)) {
          prefix = ",";
        }
      }
    }
  }
}

static void FormatFeatureBlockQuals (
  StringItemPtr    ffstring,
  IntAsn2gbJobPtr  ajp,
  Asn2gbSectPtr    asp,
  BioseqPtr        bsp,
  Uint1            featdeftype,
  ValNodePtr       gene_syn,
  CharPtr          lasttype,
  SeqLocPtr        location,
  BioseqPtr        prod,
  CharPtr          protein_pid_g,
  QualValPtr       qvp,
  Int4             left,
  Int4             right,
  Uint1            strand,
  SeqFeatPtr       sfp,
  BioseqPtr        target,
  IntFeatBlockPtr  ifp,
  Boolean          is_other,
  Boolean          is_journalscan,
  Boolean          is_gps,
  Boolean          is_ged
)

{
  Boolean              add_period;
  /*
  CharPtr              ascii;
  Int2                 ascii_len;
  */
  Boolean              at_end = FALSE;
  ByteStorePtr         bs;
  Char                 buf [80];
  Choice               cbaa;
  CodeBreakPtr         cbp;
  Char                 ch;
  Uint1                choice;
  ValNodePtr           citlist;
  Int4                 gi;
  Boolean              hadProtDesc = FALSE;
  DbtagPtr             dbt;
  DeltaItemPtr         dip;
  UserFieldPtr         entry;
  Int4                 exp_ev;
  GBQualPtr            gbq;
  GeneNomenclaturePtr  gnp;
  Int2                 i;
  FtQualType           idx;
  IntPrtBlockPtr       ipp;
  Boolean              isTRNA;
  Boolean              is_bc;
  Boolean              is_rf;
  Boolean              is_sc;
  Int2                 j;
  FtQualType           jdx;
  Int2                 k;
  Int2                 k_lower;
  Int2                 k_upper;
  Int4                 len;
  Boolean              link_is_na;
  FloatHi              molwt;
  SeqLocPtr            newloc;
  CharPtr              notestr;
  Char                 numbuf [32];
  Int2                 numcodons;
  Int2                 numsyns;
  ObjectIdPtr          oip;
  Boolean              okay = FALSE;
  Boolean              only_digits;
  BioseqPtr            pbsp;
  Int4                 pmid;
  Char                 pmidbuf [32];
  ValNodePtr           pmidlist;
  ValNodePtr           ppr;
  CharPtr              prefix;
  CharPtr              protein_seq = NULL;
  size_t               prtlen;
  CharPtr              ptr;
  RefBlockPtr          rbp;
  CharPtr              region;
  Uint1                residue;
  SeqCodeTablePtr      sctp;
  Int4                 sec_str;
  ValNodePtr           seq_seq;
  Uint1                seqcode;
  Char                 seqid [50];
  SeqIdPtr             sip;
  SeqLitPtr            slitp;
  SeqLocPtr            slp;
  Boolean              split;
  CharPtr              start;
  CharPtr              str;
  Boolean              suppress_period;
  CharPtr              tmp;
  tRNAPtr              trna;
  UserFieldPtr         ufp;
  UserObjectPtr        uop;
  ValNodePtr           vnp, vnp2, vnp3;
  VariationInstPtr     vip;
  VariationRefPtr      vrp;
  StringItemPtr        unique;
  Boolean              indexerVersion;

  unique = FFGetString(ajp);
  if ( unique == NULL ) return;

  indexerVersion = (Boolean) (GetAppProperty ("InternalNcbiSequin") != NULL);

  for (i = 0, idx = feat_qual_order [i]; idx != (FtQualType) 0; i++, idx = feat_qual_order [i]) {

    link_is_na = FALSE;

    lasttype = NULL;
    switch (asn2gnbk_featur_quals [idx].qualclass) {

      case Qual_class_ignore :
        break;

      case Qual_class_string :
        if (! StringHasNoText (qvp [idx].str)) {
          FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=",
            FALSE, TRUE, TILDE_TO_SPACES);
          FFAddTextToString(ffstring, "\"", qvp [idx].str, "\"",
            FALSE, TRUE, TILDE_TO_SPACES);
          FFAddOneChar(ffstring, '\n', FALSE);
        }
        break;

      case Qual_class_locus_tag :
        if (! StringHasNoText (qvp [idx].str)) {
          FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=",
            FALSE, TRUE, TILDE_TO_SPACES);
          FFAddTextToString(ffstring, "\"", qvp [idx].str, "\"",
            FALSE, TRUE, TILDE_TO_SPACES);
          FFAddOneChar(ffstring, '\n', FALSE);
        }
        break;

      case Qual_class_tilde :
        if (! StringHasNoText (qvp [idx].str)) {
          FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=",
            FALSE, TRUE, TILDE_EXPAND);
          FFAddTextToString(ffstring, "\"", qvp [idx].str, "\"",
            FALSE, TRUE, TILDE_EXPAND);
          FFAddOneChar(ffstring, '\n', FALSE);
        }
        break;

      case Qual_class_exception :
        if (! StringHasNoText (qvp [idx].str)) {
          FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=",
            FALSE, TRUE, TILDE_TO_SPACES);
          FFAddTextToString(ffstring, "\"", qvp [idx].str, "\"",
            FALSE, TRUE, TILDE_TO_SPACES);
          FFAddOneChar(ffstring, '\n', FALSE);
        }
        break;

      case Qual_class_product :
        if (StringHasNoText (qvp [idx].str) ||
            (ajp->flags.dropIllegalQuals &&
             (! AllowedValQual (featdeftype, idx, ajp->flags.forGbRelease)))) break;
          FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=",
            FALSE, TRUE, TILDE_TO_SPACES);
          FFAddTextToString(ffstring, "\"", qvp [idx].str, "\"",
            FALSE, TRUE, TILDE_TO_SPACES);
          FFAddOneChar(ffstring, '\n', FALSE);
        break;

      case Qual_class_sgml :
        if (! StringHasNoText (qvp [idx].str)) {
          /*
          if (is_journalscan) {
            ascii_len = Sgml2AsciiLen (qvp [idx].str);
            start = ascii = MemNew ((size_t) (10 + ascii_len));
            if (start != NULL) {
              ascii = Sgml2Ascii (qvp [idx].str, ascii, ascii_len + 1);

              FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=",
                FALSE, TRUE, TILDE_TO_SPACES);
              FFAddTextToString(ffstring, "\"", start, "\"",
                FALSE, TRUE, TILDE_TO_SPACES);
              FFAddOneChar(ffstring, '\n', FALSE);

              MemFree (start);
            }
          } else {
              FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=",
                FALSE, TRUE, TILDE_TO_SPACES);
              FFAddTextToString(ffstring, "\"", qvp[idx].str, "\"",
                FALSE, TRUE, TILDE_TO_SPACES);
              FFAddOneChar(ffstring, '\n', FALSE);
          }
          */
          FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=",
                            FALSE, TRUE, TILDE_TO_SPACES);
          FFAddTextToString(ffstring, "\"", qvp[idx].str, "\"",
                            FALSE, TRUE, TILDE_TO_SPACES);
          FFAddOneChar(ffstring, '\n', FALSE);
        }
        break;

      case Qual_class_boolean :
        if (qvp [idx].ble) {
          FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "\n",
                            FALSE, TRUE, TILDE_IGNORE);
        }
        break;

      case Qual_class_int :
        if (qvp [idx].num > 0) {
          if (idx == FTQUAL_transl_table) {
            sprintf (numbuf, "%ld", (long) qvp [idx].num);
            FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=",
              FALSE, TRUE, TILDE_IGNORE);
            FF_www_gcode (ajp, ffstring, numbuf);
          } else {
            sprintf (numbuf, "%ld", (long) qvp [idx].num);
            FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=",
              FALSE, TRUE, TILDE_IGNORE);
            FFAddTextToString(ffstring, NULL, numbuf, NULL,
              FALSE, TRUE, TILDE_IGNORE);
          }
          FFAddOneChar(ffstring, '\n', FALSE);
        }
        break;

      case Qual_class_evidence :
        exp_ev = qvp [idx].num;
        if (exp_ev > 0 && exp_ev <= 2) {
          FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=",
              FALSE, TRUE, TILDE_IGNORE);
          FFAddOneString(ffstring, evidenceText [exp_ev], FALSE, TRUE, TILDE_IGNORE);
          FFAddOneChar(ffstring, '\n', FALSE);
        }
        break;

      case Qual_class_valnode :
        for (vnp = qvp[idx].vnp; vnp != NULL; vnp = vnp->next) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (str != NULL) {
              FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=",
                  FALSE, TRUE, TILDE_TO_SPACES);
              FFAddTextToString(ffstring, "\"", str, "\"",
                  FALSE, TRUE, TILDE_TO_SPACES);
              FFAddOneChar(ffstring, '\n', FALSE);
          }
        }
        break;

      case Qual_class_sep_gene_syn :
        for (vnp = qvp[idx].vnp; vnp != NULL; vnp = vnp->next) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (StringHasNoText (str)) continue;
            FFAddTextToString (ffstring, "/", asn2gnbk_featur_quals [idx].name, "=",
                               FALSE, TRUE, TILDE_TO_SPACES);
            FFAddTextToString (ffstring, "\"", str, "\"",
                               FALSE, TRUE, TILDE_TO_SPACES);
            FFAddOneChar (ffstring, '\n', FALSE);
        }
        break;

      case Qual_class_gene_syn :
        numsyns = 0;
        for (vnp = qvp [idx].vnp; vnp != NULL; vnp = vnp->next) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (! StringHasNoText (str)) {
            numsyns++;
          }
        }
        if (numsyns > 0) {
          FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=\"",
                FALSE, TRUE, TILDE_TO_SPACES);
          prefix = NULL;
          for (vnp = qvp [idx].vnp; vnp != NULL; vnp = vnp->next) {
            str = (CharPtr) vnp->data.ptrvalue;
            if (! StringHasNoText (str)) {
              FFAddTextToString (ffstring, prefix, str, NULL, FALSE, FALSE, TILDE_IGNORE);
              prefix = "; ";
            }
          }
          FFAddOneChar(ffstring, '\"', FALSE);
          FFAddOneChar(ffstring, '\n', FALSE);
        }
        break;

      case Qual_class_map :
        gbq = qvp [idx].gbq;
        if (gbq == NULL || (ajp->flags.dropIllegalQuals &&
            (! AllowedValQual (featdeftype, idx, ajp->flags.forGbRelease)))) break;
        if (lasttype == NULL) {
          lasttype = gbq->qual;
        }
        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          if (! StringHasNoText (gbq->val)) {
            FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=\"",
                  FALSE, TRUE, TILDE_IGNORE);
            if (!StringIsJustQuotes (gbq->val)) {
              FFAddOneString(ffstring, gbq->val, FALSE, TRUE, TILDE_IGNORE);
            }
            FFAddOneChar(ffstring, '\"', FALSE);
            FFAddOneChar(ffstring, '\n', FALSE);
          }
          gbq = gbq->next;
        }
        break;

      case Qual_class_EC_valnode :
        for (vnp = qvp [idx].vnp; vnp != NULL; vnp = vnp->next) {
          str = (CharPtr) vnp->data.ptrvalue;
          okay = TRUE;

          if (str == NULL) continue;

          if (ajp->flags.dropIllegalQuals) {
            if (! ValidateECnumber (str)) {
              okay = FALSE;
            }
          }
          if (!okay) continue;

          FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=",
                FALSE, TRUE, TILDE_TO_SPACES);
          FFAddOneChar(ffstring, '\"', FALSE);
          FF_AddECnumber(ajp, ffstring, str);
          FFAddOneChar(ffstring, '\"', FALSE);
          FFAddOneChar(ffstring, '\n', FALSE);
        }
        break;

      case Qual_class_EC_quote :
        gbq = qvp [idx].gbq;
        if (gbq == NULL || (ajp->flags.dropIllegalQuals &&
            (! AllowedValQual (featdeftype, idx, ajp->flags.forGbRelease)))) break;
        if (lasttype == NULL) {
          lasttype = gbq->qual;
        }
        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          okay = TRUE;
          if (gbq->val == NULL) {
            okay = FALSE;
          }

          if (ajp->flags.dropIllegalQuals && okay) {
            if (! ValidateECnumber (gbq->val)) {
              okay = FALSE;
            }
          }

          if (StringHasNoText (gbq->val)) {
            okay = FALSE;
          }

          if (okay) {
            FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=",
                  FALSE, TRUE, TILDE_TO_SPACES);
            FFAddOneChar(ffstring, '\"', FALSE);
            if (!StringIsJustQuotes (gbq->val)) {
              FF_AddECnumber (ajp, ffstring, gbq->val);
            }
            FFAddOneChar(ffstring, '\"', FALSE);
            FFAddOneChar(ffstring, '\n', FALSE);
          }
          gbq = gbq->next;
        }
        break;

      case Qual_class_quote :
        gbq = qvp [idx].gbq;
        if (gbq == NULL || (ajp->flags.dropIllegalQuals &&
            (! AllowedValQual (featdeftype, idx, ajp->flags.forGbRelease)))) break;
        if (lasttype == NULL) {
          lasttype = gbq->qual;
        }
        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          if (! StringHasNoText (gbq->val)) {
            FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=\"",
                  FALSE, TRUE, TILDE_IGNORE);
            if (!StringIsJustQuotes (gbq->val)) {
              FFAddOneString(ffstring, gbq->val, FALSE, TRUE, TILDE_IGNORE);
            }
            FFAddOneChar(ffstring, '\"', FALSE);
            FFAddOneChar(ffstring, '\n', FALSE);
          }
          gbq = gbq->next;
        }
        break;

      case Qual_class_experiment :
        gbq = qvp [idx].gbq;
        if (gbq == NULL || (ajp->flags.dropIllegalQuals &&
            (! AllowedValQual (featdeftype, idx, ajp->flags.forGbRelease)))) break;
        if (lasttype == NULL) {
          lasttype = gbq->qual;
        }
        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          if (! StringHasNoText (gbq->val)) {
            FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=\"",
                  FALSE, TRUE, TILDE_IGNORE);
            if (!StringIsJustQuotes (gbq->val)) {
              AddExperimentWithPMIDlinks(ajp, ffstring, gbq->val);
            }
            FFAddOneChar(ffstring, '\"', FALSE);
            FFAddOneChar(ffstring, '\n', FALSE);
          }
          gbq = gbq->next;
        }
        break;

      case Qual_class_noquote :
        gbq = qvp [idx].gbq;
        if (gbq == NULL || (ajp->flags.dropIllegalQuals &&
            (! AllowedValQual (featdeftype, idx, ajp->flags.forGbRelease)))) break;
        if (lasttype == NULL) {
          lasttype = gbq->qual;
        }
        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          if (! StringHasNoText (gbq->val)) {
            FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=",
                  FALSE, TRUE, TILDE_IGNORE);
            FFAddOneString(ffstring, gbq->val, FALSE, TRUE, TILDE_TO_SPACES);
            FFAddOneChar(ffstring, '\n', FALSE);
          }
          gbq = gbq->next;
        }
        break;

      case Qual_class_label :
        gbq = qvp [idx].gbq;
        if (gbq == NULL || (ajp->flags.dropIllegalQuals &&
            (! AllowedValQual (featdeftype, idx, ajp->flags.forGbRelease)))) break;
        if (lasttype == NULL) {
          lasttype = gbq->qual;
        }
        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          if (! StringHasNoText (gbq->val)) {
            if (ajp->flags.checkQualSyntax) { /* single token, not just numeric */
              str = gbq->val;
              ch = *str;
              only_digits = TRUE;
              while (ch != '\0') {
                if (IS_WHITESP (ch)) break; /* only single token allowed */
                if (! IS_DIGIT (ch)) {
                  only_digits = FALSE;
                }
                str++;
                ch = *str;
              }
              if (only_digits) break; /* must not be just numeric */
            }
            FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                  FALSE, TRUE, TILDE_IGNORE);
            FFAddOneString(ffstring, gbq->val, FALSE, TRUE, TILDE_TO_SPACES);
            FFAddOneChar(ffstring, '\n', FALSE);
          }
          gbq = gbq->next;
        }
        break;

      case Qual_class_mobile_element :
        gbq = qvp [idx].gbq;
        if (gbq == NULL || (ajp->flags.dropIllegalQuals &&
            (! AllowedValQual (featdeftype, idx, ajp->flags.forGbRelease)))) break;
        if (lasttype == NULL) {
          lasttype = gbq->qual;
        }
        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          if (! StringHasNoText (gbq->val)) {
            str = gbq->val;
            if ((! ajp->flags.checkQualSyntax) || (ValidateMobileElement (str))) {

              /* mobile_element enabled as of 12/15/2006
              if (ajp->mode == RELEASE_MODE || ajp->mode == ENTREZ_MODE) {
                if (StringNICmp (str, "transposon", 10) == 0) {
                  str += 10;
                  if (*str == ':') {
                    str++;
                  }
                  FFAddTextToString(ffstring, "/", "transposon", "=\"",
                        FALSE, TRUE, TILDE_IGNORE);
                  if (!StringIsJustQuotes (str)) {
                    FFAddOneString(ffstring, str, FALSE, TRUE, TILDE_IGNORE);
                  }
                  FFAddOneChar(ffstring, '\"', FALSE);
                  FFAddOneChar(ffstring, '\n', FALSE);
                } else if (StringNICmp (str, "insertion sequence", 18) == 0) {
                  str += 18;
                  if (*str == ':') {
                    str++;
                  }
                  FFAddTextToString(ffstring, "/", "insertion_seq", "=\"",
                        FALSE, TRUE, TILDE_IGNORE);
                  if (!StringIsJustQuotes (str)) {
                    FFAddOneString(ffstring, str, FALSE, TRUE, TILDE_IGNORE);
                  }
                  FFAddOneChar(ffstring, '\"', FALSE);
                  FFAddOneChar(ffstring, '\n', FALSE);
                } else {
                  FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=\"",
                        FALSE, TRUE, TILDE_IGNORE);
                  if (!StringIsJustQuotes (str)) {
                    FFAddOneString(ffstring, str, FALSE, TRUE, TILDE_IGNORE);
                  }
                  FFAddOneChar(ffstring, '\"', FALSE);
                  FFAddOneChar(ffstring, '\n', FALSE);
                }
              } else {
                FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=\"",
                      FALSE, TRUE, TILDE_IGNORE);
                if (!StringIsJustQuotes (str)) {
                  FFAddOneString(ffstring, str, FALSE, TRUE, TILDE_IGNORE);
                }
                FFAddOneChar(ffstring, '\"', FALSE);
                FFAddOneChar(ffstring, '\n', FALSE);
              }
              */

              FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=\"",
                    FALSE, TRUE, TILDE_IGNORE);
              if (!StringIsJustQuotes (str)) {
                FFAddOneString(ffstring, str, FALSE, TRUE, TILDE_IGNORE);
              }
              FFAddOneChar(ffstring, '\"', FALSE);
              FFAddOneChar(ffstring, '\n', FALSE);
            }
          }
          gbq = gbq->next;
        }
        break;

      case Qual_class_number :
        gbq = qvp [idx].gbq;
        if (gbq == NULL || (ajp->flags.dropIllegalQuals &&
            (! AllowedValQual (featdeftype, idx, ajp->flags.forGbRelease)))) break;

        if (ajp->flags.checkQualSyntax) {
          str = gbq->val;

          if ( StringHasNoText (str) )
            break;
          while (!IS_WHITESP (*str) && *str != '\0')
            str++;
          if (! StringHasNoText (str) )
            break;
        }

        if (lasttype == NULL) {
          lasttype = gbq->qual;
        }
        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          if (! StringHasNoText (gbq->val)) {
            FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                  FALSE, TRUE, TILDE_IGNORE);
            FFAddOneString(ffstring, gbq->val, FALSE, TRUE, TILDE_TO_SPACES);
            FFAddOneChar(ffstring, '\n', FALSE);
          }
          gbq = gbq->next;
        }
        break;

      case Qual_class_usedin :
        gbq = qvp [idx].gbq;
        if (gbq == NULL || (ajp->flags.dropIllegalQuals &&
            (! AllowedValQual (featdeftype, idx, ajp->flags.forGbRelease)))) break;
        if (lasttype == NULL) {
          lasttype = gbq->qual;
        }
        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          if (! StringHasNoText (gbq->val)) {
            tmp = StringSave (gbq->val);
            str = tmp;
            len = StringLen (str);
            if (len > 1 && *str == '(' && str [len - 1] == ')' &&
                StringChr (str, ',') != NULL) {
              str++;
              while (! StringHasNoText (str)) {
                ptr = StringChr (str, ',');
                if (ptr == NULL) {
                  ptr = StringChr (str, ')');
                }
                if (ptr != NULL) {
                  *ptr = '\0';
                  ptr++;
                }
                FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                  FALSE, TRUE, TILDE_IGNORE);
                FFAddOneString(ffstring, str, FALSE, TRUE, TILDE_TO_SPACES);
                FFAddOneChar(ffstring, '\n', FALSE);
                str = ptr;
              }
            } else {
              FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                  FALSE, TRUE, TILDE_IGNORE);
              FFAddOneString(ffstring, str, FALSE, TRUE, TILDE_TO_SPACES);
              FFAddOneChar(ffstring, '\n', FALSE);
            }
            MemFree (tmp);
          }
          gbq = gbq->next;
        }
        break;

      case Qual_class_paren :
        gbq = qvp [idx].gbq;
        if (gbq == NULL || (ajp->flags.dropIllegalQuals &&
            (! AllowedValQual (featdeftype, idx, ajp->flags.forGbRelease)))) break;
        if (lasttype == NULL) {
          lasttype = gbq->qual;
        }
        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          if (! StringHasNoText (gbq->val)) {
            tmp = StringSave (gbq->val);
            str = tmp;
            len = StringLen (str);
            if (len > 1 && *str == '(' && str [len - 1] == ')' &&
                StringChr (str, ',') != NULL) {
              str++;
              while (! StringHasNoText (str)) {
                ptr = StringChr (str, ',');
                if (ptr == NULL) {
                  ptr = StringChr (str, ')');
                }
                if (ptr != NULL) {
                  *ptr = '\0';
                  ptr++;
                }
                FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=\"",
                  FALSE, TRUE, TILDE_IGNORE);
                FFAddOneString(ffstring, str, FALSE, TRUE, TILDE_TO_SPACES);
                FFAddOneChar(ffstring, '\"', FALSE);
                FFAddOneChar(ffstring, '\n', FALSE);
                str = ptr;
              }
            } else {
              FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=\"",
                  FALSE, TRUE, TILDE_IGNORE);
              FFAddOneString(ffstring, str, FALSE, TRUE, TILDE_TO_SPACES);
              FFAddOneChar(ffstring, '\"', FALSE);
              FFAddOneChar(ffstring, '\n', FALSE);
            }
            MemFree (tmp);
          }
          gbq = gbq->next;
        }
        break;

      case Qual_class_rpt :
        gbq = qvp [idx].gbq;
        if (gbq == NULL || (ajp->flags.dropIllegalQuals &&
            (! AllowedValQual (featdeftype, idx, ajp->flags.forGbRelease)))) break;

        if (lasttype == NULL) {
          lasttype = gbq->qual;
        }
        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          if (! StringHasNoText (gbq->val)) {
            tmp = StringSave (gbq->val);
            str = tmp;
            len = StringLen (str);
            if (len > 1 && *str == '(' && str [len - 1] == ')' &&
                StringChr (str, ',') != NULL) {
              str++;
              while (! StringHasNoText (str)) {
                ptr = StringChr (str, ',');
                if (ptr == NULL) {
                  ptr = StringChr (str, ')');
                }
                if (ptr != NULL) {
                  *ptr = '\0';
                  ptr++;
                }
                if ((! ajp->flags.checkQualSyntax) || (StringInStringList (str, validRptString))) {
                  FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                                    FALSE, TRUE, TILDE_IGNORE);
                  FFAddOneString(ffstring, str, FALSE, TRUE, TILDE_TO_SPACES);
                  FFAddOneChar(ffstring, '\n', FALSE);
                }
                str = ptr;
              }
            } else {
              if ((! ajp->flags.checkQualSyntax) || (StringInStringList (str, validRptString))) {
                FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                                  FALSE, TRUE, TILDE_IGNORE);
                FFAddOneString(ffstring, str, FALSE, TRUE, TILDE_TO_SPACES);
                FFAddOneChar(ffstring, '\n', FALSE);
              }
            }
            MemFree (tmp);
          }
          gbq = gbq->next;
        }
        break;

      case Qual_class_rpt_unit :
        gbq = qvp [idx].gbq;
        if (gbq == NULL || (ajp->flags.dropIllegalQuals &&
            (! AllowedValQual (featdeftype, idx, ajp->flags.forGbRelease)))) break;
        if (lasttype == NULL) {
          lasttype = gbq->qual;
        }

        /* in release_mode, must be of the form 123..4567 or a single-token label,
           or (technically illegal but common) letters and semicolons - NO LONGER CHECKED */

        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          if (! StringHasNoText (gbq->val)) {
            tmp = StringSave (gbq->val);
            str = tmp;
            len = StringLen (str);
#if 0
            if (len > 1 && *str == '(' && str [len - 1] == ')' /* &&
                StringChr (str + 1, '(') == NULL /* && StringChr (str, ',') != NULL */) {
              str++;
              while (! StringHasNoText (str)) {
                ptr = StringChr (str, ',');
                if (ptr == NULL) {
                  ptr = StringRChr (str, ')');
                }
                if (ptr != NULL) {
                  *ptr = '\0';
                  ptr++;
                }
                if ((! ajp->flags.checkQualSyntax) || (ValidateRptUnit (str))) {
                  TrimSpacesAroundString (str);
                  if (idx == FTQUAL_rpt_unit_range) {
                    FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                                      FALSE, TRUE, TILDE_IGNORE);
                    FFAddOneString(ffstring, str, FALSE, TRUE, TILDE_TO_SPACES);
                    FFAddOneChar(ffstring, '\n', FALSE);
                  } else {
                    FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=\"",
                                      FALSE, TRUE, TILDE_IGNORE);
                    FFAddOneString(ffstring, str, FALSE, TRUE, TILDE_TO_SPACES);
                    FFAddOneChar(ffstring, '\"', FALSE);
                    FFAddOneChar(ffstring, '\n', FALSE);
                  }
                }
                str = ptr;
              }
            } else {
#endif
              if ((! ajp->flags.checkQualSyntax) || (ValidateRptUnit (str))) {
                TrimSpacesAroundString (str);
                if (idx == FTQUAL_rpt_unit_range) {
                  FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                                    FALSE, TRUE, TILDE_IGNORE);
                  FFAddOneString(ffstring, str, FALSE, TRUE, TILDE_TO_SPACES);
                  FFAddOneChar(ffstring, '\n', FALSE);
                } else {
                  FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=\"",
                                    FALSE, TRUE, TILDE_IGNORE);
                  FFAddOneString(ffstring, str, FALSE, TRUE, TILDE_TO_SPACES);
                  FFAddOneChar(ffstring, '\"', FALSE);
                  FFAddOneChar(ffstring, '\n', FALSE);
                }
              }
#if 0
            }
#endif
            MemFree (tmp);
          }
          gbq = gbq->next;
        }
        break;

      case Qual_class_compare :
        gbq = qvp [idx].gbq;
        if (gbq == NULL || (ajp->flags.dropIllegalQuals &&
            (! AllowedValQual (featdeftype, idx, ajp->flags.forGbRelease)))) break;

        if (lasttype == NULL) {
          lasttype = gbq->qual;
        }
        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          if (! StringHasNoText (gbq->val)) {
            tmp = StringSave (gbq->val);
            str = tmp;
            len = StringLen (str);
            if (len > 1 && *str == '(' && str [len - 1] == ')' &&
                StringChr (str, ',') != NULL) {
              str++;
              while (! StringHasNoText (str)) {
                ptr = StringChr (str, ',');
                if (ptr == NULL) {
                  ptr = StringChr (str, ')');
                }
                if (ptr != NULL) {
                  *ptr = '\0';
                  ptr++;
                }
                if ((! ajp->flags.checkQualSyntax) || ValidateCompareQual (str, is_ged)) {
                  FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                                    FALSE, TRUE, TILDE_IGNORE);
                  FFAddOneString(ffstring, str, FALSE, TRUE, TILDE_TO_SPACES);
                  FFAddOneChar(ffstring, '\n', FALSE);
                }
                str = ptr;
              }
            } else {
              if ((! ajp->flags.checkQualSyntax) || ValidateCompareQual (str, is_ged)) {
                FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                                  FALSE, TRUE, TILDE_IGNORE);
                FFAddOneString(ffstring, str, FALSE, TRUE, TILDE_TO_SPACES);
                FFAddOneChar(ffstring, '\n', FALSE);
              }
            }
            MemFree (tmp);
          }
          gbq = gbq->next;
        }
        break;

      case Qual_class_replace :
        gbq = qvp [idx].gbq;
        if (gbq == NULL || (ajp->flags.dropIllegalQuals &&
            (! AllowedValQual (featdeftype, idx, ajp->flags.forGbRelease)))) break;
        if (lasttype == NULL) {
          lasttype = gbq->qual;
        }
        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                            FALSE, TRUE, TILDE_IGNORE);
          FFAddOneChar(ffstring, '\"', FALSE);
          if (!StringHasNoText (gbq->val)) {
            FFAddOneString(ffstring, gbq->val, FALSE, TRUE, TILDE_TO_SPACES);
          }
          FFAddOneChar(ffstring, '\"', FALSE);
          FFAddOneChar(ffstring, '\n', FALSE);
          gbq = gbq->next;
        }
        break;

      case Qual_class_consplice :
        gbq = qvp [idx].gbq;
        if (gbq == NULL || (ajp->flags.dropIllegalQuals &&
            (! AllowedValQual (featdeftype, idx, ajp->flags.forGbRelease)))) break;

        if (ajp->flags.checkQualSyntax && (! StringInStringList (gbq->val, validConsSpliceString)) ) {
          break;
        }

        if (lasttype == NULL) {
          lasttype = gbq->qual;
        }
        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          if (! StringHasNoText (gbq->val)) {
            FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                              FALSE, TRUE, TILDE_IGNORE);
            FFAddOneString(ffstring, gbq->val, FALSE, TRUE, TILDE_TO_SPACES);
            FFAddOneChar(ffstring, '\n', FALSE);
          }
          gbq = gbq->next;
        }
        break;

      case Qual_class_site :
        if (! StringHasNoText (qvp [idx].str)) {
          str = qvp [idx].str;
          FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                            FALSE, TRUE, TILDE_IGNORE);
          FFAddTextToString(ffstring, "\"", str, "\"", FALSE, TRUE, TILDE_TO_SPACES);
          FFAddOneChar(ffstring, '\n', FALSE);
        }
        break;

      case Qual_class_bond :
        if (! StringHasNoText (qvp [idx].str)) {
          FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                            FALSE, TRUE, TILDE_IGNORE);
          FFAddTextToString(ffstring, "\"", qvp[idx].str, "\"", FALSE, TRUE, TILDE_TO_SPACES);
          FFAddOneChar(ffstring, '\n', FALSE);
        }
        break;

      case Qual_class_L_R_B :
        gbq = qvp [idx].gbq;
        if (gbq == NULL || (ajp->flags.dropIllegalQuals &&
            (! AllowedValQual (featdeftype, idx, ajp->flags.forGbRelease)))) break;

        if (ajp->flags.checkQualSyntax && (! StringInStringList (gbq->val, validLRBString)) ) {
          break;
        }

        if (lasttype == NULL) {
          lasttype = gbq->qual;
        }
        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          if (! StringHasNoText (gbq->val)) {
            FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                              FALSE, TRUE, TILDE_IGNORE);
            FFAddOneString(ffstring, gbq->val, FALSE, TRUE, TILDE_TO_SPACES);
            FFAddOneChar(ffstring, '\n', FALSE);
          }
          gbq = gbq->next;
        }
        break;

      case Qual_class_sec_str :
        sec_str = qvp [idx].num;
        if (sec_str > 0 && sec_str <= 3) {
          FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                            FALSE, TRUE, TILDE_IGNORE);
          FFAddTextToString(ffstring, "\"", secStrText[sec_str], "\"", 
                            FALSE, FALSE, TILDE_IGNORE);
          FFAddOneChar(ffstring, '\n', FALSE);
        }
        break;

      case Qual_class_seq_loc :
        slp = qvp [idx].slp;
        if (slp != NULL) {
          FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                            FALSE, TRUE, TILDE_IGNORE);
          str = FFFlatLoc (ajp, target, slp, /* ajp->masterStyle */ FALSE, FALSE);
          FFAddTextToString(ffstring, "\"", str, "\"", 
                            FALSE, TRUE, TILDE_TO_SPACES);
          FFAddOneChar(ffstring, '\n', FALSE);
          MemFree (str);
        }
        break;

      case Qual_class_code_break :
        cbp = qvp [idx].cbp;
        seqcode = 0;
        sctp = NULL;
        while (cbp != NULL) {
          cbaa = cbp->aa;
          switch (cbaa.choice) {
            case 1 :
              seqcode = Seq_code_ncbieaa;
              break;
            case 2 :
              seqcode = Seq_code_ncbi8aa;
              break;
            case 3 :
              seqcode = Seq_code_ncbistdaa;
              break;
            default :
              break;
          }
          if (seqcode != 0) {
            sctp = SeqCodeTableFind (seqcode);
            if (sctp != NULL) {
              slp = NULL;
              while ((slp = SeqLocFindNext (cbp->loc, slp)) != NULL) {
                str = NULL;
                if (ajp->ajp.slp != NULL) {
                  sip = SeqIdParse ("lcl|dummy");
                  split = FALSE;
                  newloc = SeqLocReMapEx (sip, ajp->ajp.slp, slp, 0, FALSE, ajp->masterStyle);
 
                  SeqIdFree (sip);
                  if (newloc != NULL) {
                    A2GBSeqLocReplaceID (newloc, ajp->ajp.slp);
                    str = FFFlatLoc (ajp, target, newloc, ajp->masterStyle, FALSE);
                    SeqLocFree (newloc);
                  }
                } else {
                  str = FFFlatLoc (ajp, target, slp, ajp->masterStyle, FALSE);
                }
                if (str != NULL) {
                  residue = cbaa.value.intvalue;
                  ptr = Get3LetterSymbol (ajp, seqcode, sctp, residue);
                  /* O and J no longer quarantined */
                  /*
                  if (ajp->mode == RELEASE_MODE || ajp->mode == ENTREZ_MODE) {
                    if (StringICmp (ptr, "Pyl") == 0 || StringICmp (ptr, "Xle") == 0) {
                      ptr = "OTHER";
                    }
                  }
                  */
                  if (ptr == NULL) {
                    ptr = "OTHER";
                  }
                  FFAddOneString(ffstring, "/transl_except=", FALSE, FALSE, TILDE_IGNORE);
                  FFAddTextToString(ffstring, "(pos:", str, ",", FALSE, FALSE, TILDE_IGNORE);
                  FFAddTextToString(ffstring, "aa:", ptr, ")", FALSE, FALSE, TILDE_IGNORE);
                  FFAddOneChar(ffstring, '\n', FALSE);
                }
                MemFree (str);
              }
            }
          }
          cbp = cbp->next;
        }
        break;

      case Qual_class_anti_codon :
        slp = qvp [FTQUAL_anticodon].slp;
        newloc = NULL;
        if (slp != NULL && ajp->ajp.slp != NULL) {
          sip = SeqIdParse ("lcl|dummy");
          split = FALSE;
          newloc = SeqLocReMapEx (sip, ajp->ajp.slp, slp, 0, FALSE, ajp->masterStyle);
          /*
          newloc = SeqLocCopyRegion (sip, slp, bsp, left, right, strand, &split);
          */
          SeqIdFree (sip);
          slp = newloc;
          if (newloc != NULL) {
            A2GBSeqLocReplaceID (newloc, ajp->ajp.slp);
          }
        }
        str = qvp [FTQUAL_trna_aa].str;
        if (slp != NULL && StringDoesHaveText (str)) {
          tmp = FFFlatLoc (ajp, target, slp, ajp->masterStyle, FALSE);
          if (tmp != NULL) {
            if (ajp->mode == RELEASE_MODE &&
                (StringStr (tmp, "join") != NULL ||
                 StringStr (tmp, "order") != NULL ||
                 StringStr (tmp, "complement") != NULL)) {
              /* !!! join in anticodon quarantined pending collab approval !!! */
            } else {
              FFAddTextToString (ffstring, "/anticodon=(pos:", tmp, ",",
                                 FALSE, FALSE, TILDE_IGNORE);
              FFAddTextToString(ffstring, "aa:", str, ")",
                                FALSE, FALSE, TILDE_IGNORE);
              FFAddOneChar(ffstring, '\n', FALSE);
            }
          }
          MemFree (tmp);
        }
        if (newloc != NULL) {
          SeqLocFree (newloc);
        }
        break;

      case Qual_class_trna_codons :
        trna = qvp [idx].trp;
        if (trna) {
          numcodons = ComposeCodonsRecognizedString (trna, numbuf, sizeof (numbuf));
          if (numcodons < 1 || StringHasNoText (numbuf)) {
          } else {
            FFAddTextToString(ffstring, "/", "codon_recognized", "=\"",
                              FALSE, TRUE, TILDE_IGNORE);
            FFAddOneString(ffstring, numbuf, FALSE, TRUE, TILDE_TO_SPACES);
            FFAddOneChar(ffstring, '\"', FALSE);
            FFAddOneChar(ffstring, '\n', FALSE);
          }
        }
        break;

      case Qual_class_codon :
        gbq = qvp [idx].gbq;
        if (gbq == NULL || (ajp->flags.dropIllegalQuals &&
            (! AllowedValQual (featdeftype, idx, ajp->flags.forGbRelease)))) break;
        if (lasttype == NULL) {
          lasttype = gbq->qual;
        }
        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          if (! StringHasNoText (gbq->val)) {
            FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                              FALSE, FALSE, TILDE_IGNORE);
            FFAddOneString(ffstring, gbq->val, FALSE, FALSE, TILDE_TO_SPACES);
            FFAddOneChar(ffstring, '\n', FALSE);
          }
          gbq = gbq->next;
        }
        break;

      case Qual_class_pubset :
        vnp = qvp [idx].vnp;
        if (vnp != NULL && asp != NULL && asp->referenceArray != NULL) {
          citlist = NULL;
          pmidlist = NULL;
          for (ppr = vnp->data.ptrvalue; ppr != NULL; ppr = ppr->next) {
            j = MatchRef (ppr, asp->referenceArray, asp->numReferences);
            if (j > 0) {
              ValNodeAddInt (&citlist, 0, (Int4) j);
            } else if (is_other && ppr->choice == PUB_PMid && ajp->mode != RELEASE_MODE) {
              pmid = ppr->data.intvalue;
              ValNodeAddInt (&pmidlist, 0, (Int4) pmid);
            }
          }
          citlist = ValNodeSort (citlist, SortVnpByInt);
          pmidlist = ValNodeSort (pmidlist, SortVnpByInt);
          for (vnp = citlist; vnp != NULL; vnp = vnp->next) {
            j = (Int2) vnp->data.intvalue;
            if (j > 0) {
              sprintf (numbuf, "%d", (int) j);
              FFAddOneString(ffstring, "/citation=[", FALSE, TRUE, TILDE_TO_SPACES);
              pmid = 0;
              if( GetWWW (ajp) && asp->numReferences > 0 ) {
                /* binary search for reference that matches serial number j */
                k_lower = 0;
                k_upper = (asp->numReferences - 1);
                while( k_lower <= k_upper ) {
                  k = (k_upper + k_lower) / 2;
                  rbp = asp->referenceArray [k];
                  if( rbp->serial == j ) {
                    pmid = rbp->pmid;
                    break;
                  } else if( rbp->serial < j ) {
                    k_lower = (k+1);
                  } else { /* rbp->serial > j */
                    k_upper = (k-1);
                  }
                }
              }
              if (pmid > 0 && GetWWW (ajp)) {
                FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
                FF_Add_NCBI_Base_URL (ffstring, link_muid);
                sprintf (pmidbuf, "%ld", (long) pmid);
                FFAddTextToString(ffstring, NULL, pmidbuf, "\">", FALSE, FALSE, TILDE_IGNORE);
                FFAddOneString(ffstring, numbuf, FALSE, FALSE, TILDE_IGNORE);
                FFAddOneString(ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
              } else {
                FFAddOneString(ffstring, numbuf, FALSE, FALSE, TILDE_IGNORE);
              }
              FFAddOneString(ffstring, "]", FALSE, FALSE, TILDE_IGNORE);
              FFAddOneChar(ffstring, '\n', FALSE);
            }
          }
          for (vnp = pmidlist; vnp != NULL; vnp = vnp->next) {
            pmid = (Int4) vnp->data.intvalue;
            if (pmid > 0) {
              sprintf (pmidbuf, "%ld", (long) pmid);
              FFAddOneString(ffstring, "/citation=[PUBMED ", FALSE, TRUE, TILDE_TO_SPACES);
              if (GetWWW (ajp)) {

                FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
                FF_Add_NCBI_Base_URL (ffstring, link_muid);
                FFAddTextToString(ffstring, NULL, pmidbuf, "\">", FALSE, FALSE, TILDE_IGNORE);
                FFAddOneString(ffstring, pmidbuf, FALSE, FALSE, TILDE_IGNORE);
                FFAddOneString(ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
              } else {
                FFAddOneString(ffstring, pmidbuf, FALSE, FALSE, TILDE_IGNORE);
              }
              FFAddOneString(ffstring, "]", FALSE, FALSE, TILDE_IGNORE);
              /*
              FFAddTextToString(ffstring, "/citation=[PUBMED ", pmidbuf, "]",
                                FALSE, TRUE, TILDE_TO_SPACES);
              */
              FFAddOneChar(ffstring, '\n', FALSE);
            }
          }
          citlist = ValNodeFree (citlist);
          pmidlist = ValNodeFree (pmidlist);
        }
        break;

      case Qual_class_db_xref :
        for (vnp = qvp [idx].vnp; vnp != NULL; vnp = vnp->next) {
          buf [0] = '\0';
          dbt = (DbtagPtr) vnp->data.ptrvalue;
          if (dbt != NULL && (! StringHasNoText (dbt->db))) {
            oip = dbt->tag;
            if (oip != NULL) {

              okay = TRUE;
              if (ajp->flags.dropBadDbxref) {
                /* if RELEASE_MODE, drop unknown dbtag */

                okay = FALSE;
                if (DbxrefIsValid (dbt->db, &is_rf, &is_sc, &is_bc, NULL)) {
                  if (is_bc) {
                    /* case counts, so suppress if bad case */
                  } else if (is_rf) {
                    if (is_gps || is_other) {
                      okay = TRUE;
                    }
                  } else if (is_sc) {
                    /* show, but warn in validator */
                    okay = TRUE;
                  } else {
                    okay = TRUE;
                  }
                }

                /*okay = FALSE;
                for (j = 0; legalDbXrefs [j] != NULL; j++) {
                  if (StringCmp (dbt->db, legalDbXrefs [j]) == 0) {
                    okay = TRUE;
                  }
                }
                if (! okay) {
                  if (is_gps || is_other) {
                    for (j = 0; legalRefSeqDbXrefs [j] != NULL; j++) {
                      if (StringCmp (dbt->db, legalRefSeqDbXrefs [j]) == 0) {
                        okay = TRUE;
                      }
                    }
                  }
                }
                */
              }

              if (StringICmp (dbt->db, "taxon") == 0 ||
                  StringCmp (dbt->db, "PID") == 0 ||
                  StringCmp (dbt->db, "GI") == 0) {
                okay = FALSE;
              }
              if (okay && idx == FTQUAL_db_xref && qvp [FTQUAL_gene_xref].vnp != NULL) {
                if (DbxrefAlreadyInGeneXref (dbt, qvp [FTQUAL_gene_xref].vnp)) {
                  okay = FALSE;
                }
              }

              if (okay) {
                if (! StringHasNoText (oip->str)) {
                  if (StringLen (oip->str) < 80) {
                    sprintf (buf, "%s", oip->str);
                  }
                } else {
                  sprintf (buf, "%ld", (long) oip->id);
                }
              }
            }
          }
          if (! StringHasNoText (buf)) {
            if (StringICmp (buf, protein_pid_g) != 0) {
              /* already sorted and uniqued by BasicSeqEntryCleanup, per feature */
              if (dbt != NULL) {
                if (StringICmp (dbt->db, "LocusID") == 0 || StringICmp (dbt->db, "InterimID") == 0) {
                  if (FFStringSearch (ffstring, dbt->db, 0) >= 0) {
                    okay = FALSE;
                  }
                }
              }
              if (okay) {
                FFAddOneString(ffstring, "/db_xref=\"", FALSE, FALSE, TILDE_IGNORE);
                FF_www_db_xref(ajp, ffstring, dbt->db, buf, bsp);
                FFAddOneString(ffstring, "\"\n", FALSE, FALSE, TILDE_IGNORE);
              }
            }
          }
        }
        break;

      case Qual_class_variation_id :
        dbt = qvp [idx].dbt;
        if (dbt != NULL) {
          buf [0] = '\0';
          if (StringICmp (dbt->db, "dbSNP") == 0) {
            oip = dbt->tag;
            if (oip != NULL && StringDoesHaveText (oip->str)) {
              str = oip->str;
              if (StringNICmp (str, "rs", 2) == 0) {
                FFAddOneString(ffstring, "/db_xref=\"", FALSE, FALSE, TILDE_IGNORE);
                FF_www_db_xref(ajp, ffstring, dbt->db, str + 2, bsp);
                FFAddOneString(ffstring, "\"\n", FALSE, FALSE, TILDE_IGNORE);
              }
            }
          }
        }
        break;

      case Qual_class_delta_item :
        for (vnp = qvp [idx].vnp; vnp != NULL; vnp = vnp->next) {
          dip = (DeltaItemPtr) vnp->data.ptrvalue;
          if (dip == NULL) continue;
          seq_seq = dip->Seq_seq;
          if (seq_seq != NULL && seq_seq->choice == Seq_seq_literal) {
            slitp = (SeqLitPtr) seq_seq->data.ptrvalue;
            if (slitp != NULL) {
              if (slitp->length > 0 && slitp->seq_data_type != Seq_code_gap && slitp->seq_data != NULL) {
                str = (CharPtr) MemNew ((size_t) (slitp->length + 6));
                if (str != NULL) {
                  SeqPortStreamLit (slitp, 0, (Pointer) str, NULL);
                  FFAddOneString(ffstring, "/replace=\"", FALSE, FALSE, TILDE_IGNORE);
                  if (StringDoesHaveText (str)) {
                    ptr = str;
                    ch = *ptr;
                    while (ch != '\0') {
                      if (IS_UPPER (ch)) {
                        ch = TO_LOWER (ch);
                        *ptr = ch;
                      }
                      ptr++;
                      ch = *ptr;
                    }
                    FFAddOneString(ffstring, str, FALSE, FALSE, TILDE_IGNORE);
                  }
                  FFAddOneString(ffstring, "\"\n", FALSE, FALSE, TILDE_IGNORE);
                }
                MemFree (str);
              }
            }
          }
        }
        break;

      case Qual_class_variation_set:
        for (vnp = qvp [idx].vnp; vnp != NULL; vnp = vnp->next) {
          vrp = (VariationRefPtr) vnp->data.ptrvalue;
          if (vrp == NULL) continue;
          vnp2 = vrp->data;
          if (vnp2 == NULL) continue;
          if (vnp2->choice != VarRefData_instance) continue;
          vip = (VariationInstPtr) vnp2->data.ptrvalue;
          if (vip == NULL) continue;
          for (vnp3 = vip->delta; vnp3 != NULL; vnp3 = vnp3->next) {
            dip = (DeltaItemPtr) vnp3->data.ptrvalue;
            if (dip == NULL) continue;
            seq_seq = dip->Seq_seq;
            if (seq_seq != NULL && seq_seq->choice == Seq_seq_literal) {
              slitp = (SeqLitPtr) seq_seq->data.ptrvalue;
              if (slitp != NULL) {
                if (slitp->length > 0 && slitp->seq_data_type != Seq_code_gap && slitp->seq_data != NULL) {
                  str = (CharPtr) MemNew ((size_t) (slitp->length + 6));
                  if (str != NULL) {
                    SeqPortStreamLit (slitp, 0, (Pointer) str, NULL);
                    FFAddOneString(ffstring, "/replace=\"", FALSE, FALSE, TILDE_IGNORE);
                    if (StringDoesHaveText (str)) {
                      ptr = str;
                      ch = *ptr;
                      while (ch != '\0') {
                        if (IS_UPPER (ch)) {
                          ch = TO_LOWER (ch);
                          *ptr = ch;
                        }
                        ptr++;
                        ch = *ptr;
                      }
                      FFAddOneString(ffstring, str, FALSE, FALSE, TILDE_IGNORE);
                    }
                    FFAddOneString(ffstring, "\"\n", FALSE, FALSE, TILDE_IGNORE);
                  }
                  MemFree (str);
                }
              }
            }
          }
        }
        break;

      case Qual_class_nuc_id :
        link_is_na = TRUE;
        /* fall through */
      case Qual_class_prt_id :
        sip = qvp [idx].sip;
        if (sip != NULL) {
          /* should always be found above for protein_id or transcript_id
          prod = BioseqFind (sip);
          */
          if (prod != NULL) {
            gi = 0;
            for (sip = prod->id; sip != NULL; sip = sip->next) {
              if (sip->choice == SEQID_GI) {
                gi = sip->data.intvalue;
              }
            }
            choice = 0;
            for (sip = prod->id; sip != NULL; sip = sip->next) {
              if (sip->choice == SEQID_GENBANK ||
                  sip->choice == SEQID_EMBL ||
                  sip->choice == SEQID_DDBJ ||
                  sip->choice == SEQID_OTHER ||
                  sip->choice == SEQID_TPG ||
                  sip->choice == SEQID_TPE ||
                  sip->choice == SEQID_TPD ||
                  sip->choice == SEQID_GPIPE) {
                choice = sip->choice;
                if (SeqIdWrite (sip, seqid, PRINTID_TEXTID_ACC_VER, sizeof (seqid)) != NULL) {
                  FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=\"",
                                    FALSE, FALSE, TILDE_IGNORE);
                  FF_www_nuc_or_prot_id (ajp, ffstring, seqid, gi, link_is_na);
                  FFAddOneString(ffstring, "\"\n", FALSE, FALSE, TILDE_IGNORE);
                }
              } else if (sip->choice == SEQID_GI) {
                if (choice == 0) {
                  sprintf (seqid, "%ld", (long) sip->data.intvalue);
                  FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=\"",
                                    FALSE, FALSE, TILDE_IGNORE);
                  FF_www_nuc_or_prot_id (ajp, ffstring, seqid, gi, link_is_na);
                  FFAddOneString(ffstring, "\"\n", FALSE, FALSE, TILDE_IGNORE);
                }
                sprintf (seqid, "%ld", (long) sip->data.intvalue);
                FFAddOneString(ffstring, "/db_xref=\"", FALSE, FALSE, TILDE_IGNORE);
                FF_www_db_xref(ajp, ffstring, "GI", seqid, bsp);
                FFAddOneString(ffstring, "\"\n", FALSE, FALSE, TILDE_IGNORE);
              } else if (sip->choice == SEQID_GENERAL) {
                dbt = (DbtagPtr) sip->data.ptrvalue;
                if (dbt != NULL && StringCmp (dbt->db, "PID") == 0) {
                  /*
                  oip = dbt->tag;
                  if (oip != NULL) {
                    if (! StringHasNoText (oip->str)) {
                      sprintf (seqid, "PID:%s", oip->str);
                      NewContLine ();
                      gb_AddString ("/db_xref=\"", seqid, "\"", FALSE, TRUE, TILDE_TO_SPACES);
                    }
                  }
                  */
                } else if (dbt != NULL) {
                  pbsp = BioseqFind (sip);
                  if (pbsp != NULL && pbsp->id != NULL && /* pbsp->id->next == NULL && */ OnlyOneRealGeneral (pbsp->id)) {
                    dbt = (DbtagPtr) sip->data.ptrvalue;
                    if (dbt != NULL &&
                        !IsSkippableDbtag(dbt)) {
                      if (SeqIdWrite (sip, seqid, PRINTID_REPORT, sizeof (seqid)) != NULL) {
                        FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=\"",
                                          FALSE, FALSE, TILDE_IGNORE);
                        FF_www_nuc_or_prot_id (ajp, ffstring, seqid, gi, link_is_na);
                        FFAddOneString(ffstring, "\"\n", FALSE, FALSE, TILDE_IGNORE);
                      }
                    }
                  }
                }
              }
            }
          } else {
            if (sip->choice == SEQID_GI) {
              gi = sip->data.intvalue;
              if (GetAccnVerFromServer (gi, seqid)) {
#ifdef OS_UNIX
  if (getenv ("ASN2GB_PSF_DEBUG") != NULL) {
    printf ("GetAccnVerFromServer returned %s\n", seqid);
  }
#endif
                if ((! ajp->flags.dropIllegalQuals) || ValidateAccn (seqid) == 0) {
                  FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=\"",
                                    FALSE, FALSE, TILDE_IGNORE);
                  FF_www_nuc_or_prot_id (ajp, ffstring, seqid, gi, link_is_na);
                  FFAddOneString(ffstring, "\"\n", FALSE, FALSE, TILDE_IGNORE);
                } else {
                  ajp->relModeError = TRUE;
                }
              } else {
                sip = GetSeqIdForGI (gi);
                if (sip != NULL && SeqIdWrite (sip, seqid, PRINTID_TEXTID_ACC_VER, sizeof (seqid)) != NULL) {
#ifdef OS_UNIX
  if (getenv ("ASN2GB_PSF_DEBUG") != NULL) {
    printf ("GetSeqIdForGI returned %s\n", seqid);
  }
#endif
                  if ((! ajp->flags.dropIllegalQuals) || ValidateAccn (seqid) == 0) {
                    FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=\"",
                                    FALSE, FALSE, TILDE_IGNORE);
                    FF_www_nuc_or_prot_id (ajp, ffstring, seqid, gi, link_is_na);
                    FFAddOneString(ffstring, "\"\n", FALSE, FALSE, TILDE_IGNORE);
                  } else {
                    ajp->relModeError = TRUE;
                  }
                } else if (! ajp->flags.dropIllegalQuals) {
                  sprintf (seqid, "%ld", (long) gi);
                  FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=\"",
                                    FALSE, FALSE, TILDE_IGNORE);
                  FF_www_nuc_or_prot_id (ajp, ffstring, seqid, gi, link_is_na);
                  FFAddOneString(ffstring, "\"\n", FALSE, FALSE, TILDE_IGNORE);
                } else {
                  ajp->relModeError = TRUE;
                }
              }

              sprintf (seqid, "%ld", (long) gi);
              FFAddOneString(ffstring, "/db_xref=\"", FALSE, FALSE, TILDE_IGNORE);
              FF_www_db_xref(ajp, ffstring, "GI", seqid, bsp);
              FFAddOneString(ffstring, "\"\n", FALSE, FALSE, TILDE_IGNORE);
            } else if (SeqIdWrite (sip, seqid, PRINTID_TEXTID_ACC_VER, sizeof (seqid)) != NULL) {
              gi = GetGIForSeqId (sip);
              if ((! ajp->flags.dropIllegalQuals) || ValidateAccn (seqid) == 0) {
                FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=\"",
                                  FALSE, FALSE, TILDE_IGNORE);
                FF_www_nuc_or_prot_id (ajp, ffstring, seqid, gi, link_is_na);
                FFAddOneString(ffstring, "\"\n", FALSE, FALSE, TILDE_IGNORE);
              } else {
                ajp->relModeError = TRUE;
              }

              if (gi > 0) {
                sprintf (seqid, "%ld", (long) gi);
                FFAddOneString(ffstring, "/db_xref=\"", FALSE, FALSE, TILDE_IGNORE);
                FF_www_db_xref(ajp, ffstring, "GI", seqid, bsp);
                FFAddOneString(ffstring, "\"\n", FALSE, FALSE, TILDE_IGNORE);
              }
            }
          }
        }
        break;

      case Qual_class_mol_wt :
        if (qvp [idx].ble) {
          if (ifp != NULL && ifp->isPrt) {
            ipp = (IntPrtBlockPtr) ifp;
            molwt = MolWtForProtFeat (bsp, sfp, ipp);
            if (molwt > 0.01) {
              sprintf (buf, "%ld", (long) (molwt + 0.5));
              TrimSpacesAroundString (buf);
              FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                                FALSE, FALSE, TILDE_IGNORE);
              FFAddOneString(ffstring, buf, FALSE, FALSE, TILDE_TO_SPACES);
              FFAddOneChar(ffstring, '\n', FALSE);
            }
          }
        }
        break;

      case Qual_class_translation :
        if (qvp [idx].ble && (! ajp->hideTranslation)) {
          if ((prod == NULL && ajp->transIfNoProd) || ajp->alwaysTranslCds) {
            bs = ProteinFromCdRegionEx (sfp, TRUE, FALSE);
            if (bs != NULL) {
              str = BSMerge (bs, NULL);
              bs = BSFree (bs);
              if (str != NULL) {
                ptr = str;
                ch = *ptr;
                while (ch != '\0') {
                  *ptr = TO_UPPER (ch);
                  ptr++;
                  ch = *ptr;
                }
                prtlen = StringLen (str);
                if (prtlen > 1) {
                   if (str [prtlen - 1] == '*') {
                     str [prtlen - 1] = '\0';
                   }
                }
                if (! StringHasNoText (str)) {
                  /*
                  if (ajp->mode == RELEASE_MODE || ajp->mode == ENTREZ_MODE) {
                    ChangeOandJtoX (str);
                  }
                  */
                  FFAddTextToString(ffstring, "/translation=\"", str, "\"", 
                                    FALSE, TRUE, TILDE_TO_SPACES);
                  FFAddOneChar(ffstring, '\n', FALSE);
                }
                MemFree (str);
              }
            } else {
              ajp->relModeError = TRUE;
            }
          } else if (prod != NULL) {
            len = SeqLocLen (sfp->product);
            if (len > 0) {
              if (SeqLocStart (location) == 0 || (bsp != NULL && SeqLocStop (location) == bsp->length - 1)) {
                at_end = TRUE;
              }
              str = (CharPtr) MemNew ((size_t) (len + 1) * sizeof (Char));
              protein_seq = str;
              /*
              if (ajp->flags.iupacaaOnly) {
                code = Seq_code_iupacaa;
              } else {
                code = Seq_code_ncbieaa;
              }
              */
              SeqPortStreamLoc (sfp->product, STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL, (Pointer) &protein_seq, SaveGBSeqTranslation);
              if (! StringHasNoText (str)) {
                /*
                if (ajp->mode == RELEASE_MODE || ajp->mode == ENTREZ_MODE) {
                  ChangeOandJtoX (str);
                }
                */
                FFAddTextToString(ffstring, "/translation=\"", str, "\"", 
                                  FALSE, TRUE, TILDE_TO_SPACES);
                FFAddOneChar(ffstring, '\n', FALSE);
              }
              MemFree (str);
            } else {
              ajp->relModeError = TRUE;
            }
          }
        }
        break;
      
      case Qual_class_transcription :
        if (qvp [idx].ble && ajp->showTranscript) {
          if ((prod == NULL && ajp->transIfNoProd) || ajp->alwaysTranslCds) {
            str = GetSequenceByFeature (sfp);
            if (str != NULL) {
              ptr = str;
              ch = *ptr;
              while (ch != '\0') {
                *ptr = TO_UPPER (ch);
                ptr++;
                ch = *ptr;
              }
              if (! StringHasNoText (str)) {
                FFAddTextToString(ffstring, "/transcription=\"", str, "\"", 
                                  FALSE, TRUE, TILDE_TO_SPACES);
                FFAddOneChar(ffstring, '\n', FALSE);
              }
              MemFree (str);
            }
          } else if (prod != NULL) {
            len = SeqLocLen (sfp->product);
            if (len > 0) {
              str = (CharPtr) MemNew ((size_t) (len + 2) * sizeof (Char));
              SeqPortStreamLoc (sfp->product, STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL, (Pointer) str, NULL);
              if (! StringHasNoText (str)) {
                FFAddTextToString(ffstring, "/transcription=\"", str, "\"", 
                                  FALSE, TRUE, TILDE_TO_SPACES);
                FFAddOneChar(ffstring, '\n', FALSE);
              }
              MemFree (str);
            }
          } else {
            str = GetSequenceByFeature (sfp);
            if (str != NULL) {
              ptr = str;
              ch = *ptr;
              while (ch != '\0') {
                *ptr = TO_UPPER (ch);
                ptr++;
                ch = *ptr;
              }
              if (! StringHasNoText (str)) {
                FFAddTextToString(ffstring, "/transcription=\"", str, "\"", 
                                  FALSE, TRUE, TILDE_TO_SPACES);
                FFAddOneChar(ffstring, '\n', FALSE);
              }
              MemFree (str);
            }
          }
        }
        break;
      
      case Qual_class_peptide :
        if (qvp [idx].ble) {
          if (ajp->showPeptide) {
            str = GetSequenceByFeature (sfp);
            if (str != NULL) {
              ptr = str;
              ch = *ptr;
              while (ch != '\0') {
                *ptr = TO_UPPER (ch);
                ptr++;
                ch = *ptr;
              }
              if (! StringHasNoText (str)) {
                FFAddTextToString(ffstring, "/peptide=\"", str, "\"", 
                                  FALSE, TRUE, TILDE_TO_SPACES);
                FFAddOneChar(ffstring, '\n', FALSE);
              }
              MemFree (str);
            }
          }
        }
        break;
      
       case Qual_class_tag_peptide :
        if (! StringHasNoText (qvp [idx].str)) {
          FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=",
            FALSE, TRUE, TILDE_TO_SPACES);
          FFAddTextToString(ffstring, NULL, qvp [idx].str, NULL,
            FALSE, TRUE, TILDE_TO_SPACES);
          FFAddOneChar(ffstring, '\n', FALSE);
        }
        break;

     case Qual_class_illegal :
        for (vnp = qvp [idx].vnp; vnp != NULL; vnp = vnp->next) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (str != NULL) {
            if (ajp->mode == SEQUIN_MODE) {
              if (StringNICmp (str, "/orig_protein_id=", 17) == 0) continue;
              if (StringNICmp (str, "/orig_transcript_id=", 20) == 0) continue;
            }
            FFAddOneString(ffstring, str, FALSE, FALSE, TILDE_TO_SPACES);
            FFAddNewLine(ffstring);
          }
        }
        break;

      case Qual_class_note :
        if (! ajp->flags.goQualsToNote) {

          /* in GenBank sequin_mode and dump_mode, and in RefSeq, GO terms show up as separate /qualifiers */

          for (j = 0, jdx = feat_note_order [j]; jdx != 0; j++, jdx = feat_note_order [j]) {

            link_is_na = FALSE;

            switch (asn2gnbk_featur_quals [jdx].qualclass) {

              case Qual_class_go :
                if (qvp [jdx].ufp != NULL) {
                  if (ajp->mode == ENTREZ_MODE) {
                    str = GetCombinedGOtext (qvp [jdx].ufp, ajp);
                    if (StringDoesHaveText (str)) {
                      FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[jdx].name, "=",
                                        FALSE, TRUE, TILDE_IGNORE);
                      FFAddTextToString(ffstring, "\"", str, "\"", 
                                        FALSE, FALSE, TILDE_IGNORE);
                      FFAddOneChar(ffstring, '\n', FALSE);
                    }
                    MemFree (str);
                  } else {
                    for (entry = qvp [jdx].ufp; entry != NULL; entry = entry->next) {
                      if (entry == NULL || entry->choice != 11) break;
                      ufp = (UserFieldPtr)  entry->data.ptrvalue;
                      str = GetGOtext (ufp, ajp, (Boolean) (ajp->mode == ENTREZ_MODE));
                      if (! StringHasNoText (str)) {
                        FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[jdx].name, "=",
                                          FALSE, TRUE, TILDE_IGNORE);
                        FFAddTextToString(ffstring, "\"", str, "\"", 
                                          FALSE, FALSE, TILDE_IGNORE);
                        FFAddOneChar(ffstring, '\n', FALSE);
                      }
                      MemFree (str);
                    }
                  }
                }
                break;

              default :
                break;
            }
          }
        }

        if (! ajp->flags.refSeqQualsToNote) {

          /* in entrez_mode, sequin_mode and dump_mode in RefSeq, RefSeq-specific qualifiers show up as separate /qualifiers */

          for (j = 0, jdx = feat_note_order [j]; jdx != 0; j++, jdx = feat_note_order [j]) {
            switch (asn2gnbk_featur_quals [jdx].qualclass) {

              case Qual_class_nomenclature :
                uop = qvp [jdx].uop;
                if (uop != NULL) {
                  str = NULL;
                  VisitUserObjectsInUop (sfp->ext, (Pointer) &str, GetNomenclatureText);
                  if (! StringHasNoText (str)) {
                    FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[jdx].name, "=",
                                      FALSE, TRUE, TILDE_IGNORE);
                    FFAddTextToString(ffstring, "\"", str, "\"", 
                                      FALSE, FALSE, TILDE_IGNORE);
                    FFAddOneChar(ffstring, '\n', FALSE);
                    prefix = "; ";
                    add_period = FALSE;
                  }
                  MemFree (str);
                }
                break;

              case Qual_class_gene_nomen :
                gnp = qvp [jdx].gnp;
                if (gnp != NULL) {
                  str = GetNomenclature (gnp);
                  if (! StringHasNoText (str)) {
                    FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[jdx].name, "=",
                                      FALSE, TRUE, TILDE_IGNORE);
                    FFAddTextToString(ffstring, "\"", str, "\"", 
                                      FALSE, FALSE, TILDE_IGNORE);
                    FFAddOneChar(ffstring, '\n', FALSE);
                    prefix = "; ";
                    add_period = FALSE;
                  }
                  MemFree (str);
                }
                break;

              default :
                break;
            }
          }
        }

        /*head = NULL;*/
        notestr = NULL;
        prefix = NULL;
        add_period = FALSE;
        suppress_period = FALSE;
        lasttype = NULL;
        isTRNA = FALSE;
        

#ifdef DISPLAY_STRINGS
        s_DisplayQVP (qvp, feat_note_order);
#endif
        for (j = 0, jdx = feat_note_order [j]; jdx != 0; j++, jdx = feat_note_order [j]) {
          switch (asn2gnbk_featur_quals [jdx].qualclass) {

            case Qual_class_string :
              if (! StringHasNoText (qvp [jdx].str)) {
                if (jdx == FTQUAL_figure) {
                  if (!IsEllipsis (qvp [jdx].str))
                    s_RemovePeriodFromEnd (qvp [jdx].str);
                  sprintf (buf, "This sequence comes from %s", qvp [jdx].str);
                  FFAddString_NoRedund (unique, prefix, buf, NULL, TRUE);
                  add_period = FALSE;
                } else if (jdx == FTQUAL_maploc) {
                  if (!IsEllipsis (qvp [jdx].str))
                    s_RemovePeriodFromEnd (qvp [jdx].str);
                  sprintf (buf, "Map location %s", qvp [jdx].str);
                  FFAddString_NoRedund (unique, prefix, buf, NULL, TRUE);
                  add_period = FALSE;
                } else if (jdx == FTQUAL_seqannot_note) {
                  str = StringSave (qvp [jdx].str);
                  TrimSpacesAndJunkFromEnds (str, TRUE);
                  if (! IsEllipsis (str))
                    add_period = s_RemovePeriodFromEnd (str);
                  /*  NOTE -- The following function call cleans up some strings
                      (i.e., U34661 & U31565) but should be commented back
                      in only if the problem can't be fixed upstream of here

                  s_StringCleanup(str);

                  */
                  FFAddString_NoRedund (unique, prefix, str, NULL, TRUE);
                  MemFree (str);
                  if (hadProtDesc) {
                    suppress_period = TRUE;
                  }
                } else if (jdx == FTQUAL_seqfeat_note) {
                  str = StringSave (qvp [jdx].str);
                  if (indexerVersion) {
                    TrimSpacesAroundString (str);
                  } else {
                    TrimSpacesAndJunkFromEnds (str, TRUE);
                  }
                  if (! IsEllipsis (str))
                    add_period = s_RemovePeriodFromEnd (str);
                  /*  NOTE -- The following function call cleans up some strings
                      (i.e., U34661 & U31565) but should be commented back
                      in only if the problem can't be fixed upstream of here

                  s_StringCleanup(str);

                  */
                  FFAddString_NoRedund (unique, prefix, str, NULL, TRUE);
                  MemFree (str);
                  if (hadProtDesc) {
                    suppress_period = TRUE;
                  }
                } else if (jdx == FTQUAL_prot_note) {
                  str = StringSave (qvp [jdx].str);
                  TrimSpacesAndJunkFromEnds (str, TRUE);
                  if (! IsEllipsis (str))
                    s_RemovePeriodFromEnd (str);
                  FFAddString_NoRedund (unique, prefix, str, NULL, TRUE);
                  MemFree (str);
                  add_period = FALSE;
                } else if (jdx == FTQUAL_prot_desc) {
                  str = StringSave (qvp [jdx].str);
                  TrimSpacesAndJunkFromEnds (str, TRUE);
                  if (! IsEllipsis (str))
                    add_period = s_RemovePeriodFromEnd (str);
                  FFAddString_NoRedund (unique, prefix, str, NULL, TRUE);
                  MemFree (str);
                  hadProtDesc = TRUE; /* gi|347886|gb|M96268.1|ECOUBIA */
                } else {
                  if (! IsEllipsis (qvp [jdx].str)) {
                    s_RemovePeriodFromEnd (qvp [jdx].str);
                  }
                  FFAddString_NoRedund (unique, prefix, qvp [jdx].str, NULL, TRUE);
                  add_period = FALSE;
                }
                prefix = "; ";
              }
              break;

            case Qual_class_exception :
              if (! StringHasNoText (qvp [jdx].str)) {
                if (! IsEllipsis (qvp [jdx].str)) {
                  s_RemovePeriodFromEnd (qvp [jdx].str);
                }
                if (StringCmp (prefix, "; ") == 0) {
                  prefix = "~";
                }
                FFAddString_NoRedund (unique, prefix, qvp [jdx].str, NULL, TRUE);
                add_period = FALSE;
                prefix = "; ";
              }
              break;

            case Qual_class_encodes :
              if (! StringHasNoText (qvp [jdx].str)) {
                if (! IsEllipsis (qvp [jdx].str)) {
                  s_RemovePeriodFromEnd (qvp [jdx].str);
                }
                FFAddTextToString (unique, prefix, "encodes ", NULL, FALSE, FALSE, TILDE_IGNORE);
                FFAddString_NoRedund (unique, NULL, qvp [jdx].str, NULL, TRUE);
                prefix = "; ";
                add_period = FALSE;
              }
              break;

            case Qual_class_locus_tag :
              if (! StringHasNoText (qvp [jdx].str)) {
                if (! IsEllipsis (qvp [jdx].str)) {
                  s_RemovePeriodFromEnd (qvp [jdx].str);
                }
                FFAddTextToString (unique, prefix, "locus_tag: ", NULL, FALSE, FALSE, TILDE_IGNORE);
                FFAddString_NoRedund (unique, NULL, qvp [jdx].str, NULL, TRUE);
                prefix = "; ";
                add_period = FALSE;
              }
              break;

            case Qual_class_go :
              if (ajp->flags.goQualsToNote && qvp [jdx].ufp != NULL) {
                for (entry = qvp [jdx].ufp; entry != NULL; entry = entry->next) {
                  if (entry == NULL || entry->choice != 11) break;
                  ufp = (UserFieldPtr)  entry->data.ptrvalue;
                  str = GetGOtext (ufp, ajp, (Boolean) (ajp->mode == ENTREZ_MODE));
                  if (! StringHasNoText (str)) {
                    if (StringCmp (prefix, "; ") == 0) {
                      prefix = ";\n";
                    }
                    FFAddTextToString (unique, prefix, asn2gnbk_featur_quals[jdx].name, ": ", FALSE, FALSE, TILDE_IGNORE);
                    FFAddTextToString(unique, NULL, str, NULL, FALSE, FALSE, TILDE_IGNORE);
                  }
                  MemFree (str);
                  prefix = "; ";
                  add_period = FALSE;
                }
              }
              break;

            case Qual_class_nomenclature :
              if (ajp->flags.refSeqQualsToNote) {
                uop = qvp [jdx].uop;
                if (uop != NULL) {
                  str = NULL;
                  VisitUserObjectsInUop (sfp->ext, (Pointer) &str, GetNomenclatureText);
                  if (! StringHasNoText (str)) {
                    if (StringCmp (prefix, "; ") == 0) {
                      prefix = ";\n";
                    }
                    FFAddTextToString (unique, prefix, asn2gnbk_featur_quals[jdx].name, ": ", FALSE, FALSE, TILDE_IGNORE);
                    FFAddTextToString(unique, NULL, str, NULL, FALSE, TRUE, TILDE_IGNORE);
                    prefix = "; ";
                    add_period = FALSE;
                  }
                  MemFree (str);
                  prefix = "; ";
                  add_period = FALSE;
                }
              }
              break;

            case Qual_class_gene_nomen :
              if (ajp->flags.refSeqQualsToNote) {
                gnp = qvp [jdx].gnp;
                if (gnp != NULL) {
                  str = GetNomenclature (gnp);
                  if (! StringHasNoText (str)) {
                    FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[jdx].name, "=",
                                      FALSE, TRUE, TILDE_IGNORE);
                    FFAddTextToString(ffstring, "\"", str, "\"", 
                                      FALSE, FALSE, TILDE_IGNORE);
                    FFAddOneChar(ffstring, '\n', FALSE);
                    prefix = "; ";
                    add_period = FALSE;
                  }
                  MemFree (str);
                }
              }
              break;

            case Qual_class_method :
              if (! StringHasNoText (qvp [jdx].str)) {
                if ( FFEmpty(unique) ) {
                  prefix = "Method: ";
                } else {
                  prefix = "; Method: ";
                }
                FFAddString_NoRedund (unique, prefix, qvp [jdx].str, NULL, TRUE);
                prefix = "; ";
                add_period = FALSE;
              }
              break;

            case Qual_class_valnode :
              for (vnp = qvp [jdx].vnp; vnp != NULL; vnp = vnp->next) {
                str = (CharPtr) vnp->data.ptrvalue;
                if (! StringHasNoText (str)) {
                  FFAddString_NoRedund (unique, prefix, str, NULL, TRUE);
                  prefix = "; ";
                  add_period = FALSE;
                }
              }
              break;

            /*
            case Qual_class_gene_syn :
              numsyns = 0;
              for (vnp = qvp [jdx].vnp; vnp != NULL; vnp = vnp->next) {
                str = (CharPtr) vnp->data.ptrvalue;
                if (! StringHasNoText (str)) {
                  numsyns++;
                }
              }
              if (numsyns > 0) {
                if (numsyns > 1) {
                  FFAddTextToString (unique, prefix, "synonyms: ", NULL, FALSE, FALSE, TILDE_IGNORE);
                } else {
                  FFAddTextToString (unique, prefix, "synonym: ", NULL, FALSE, FALSE, TILDE_IGNORE);
                }
                prefix = NULL;
                for (vnp = qvp [jdx].vnp; vnp != NULL; vnp = vnp->next) {
                  str = (CharPtr) vnp->data.ptrvalue;
                  if (! StringHasNoText (str)) {
                    FFAddTextToString (unique, prefix, str, NULL, FALSE, TRUE, TILDE_IGNORE);
                    prefix = ", ";
                  }
                }
                prefix = "; ";
                add_period = FALSE;
              }
              break;
            */

            case Qual_class_region :
#ifdef ASN2GNBK_STRIP_NOTE_PERIODS
              FFAddTextToString(unique, prefix, qvp [jdx].str, NULL, FALSE, TRUE, TILDE_IGNORE);
#else
              region = NULL;
              if (! StringHasNoText (qvp [jdx].str)) {
                if ( FFEmpty(unique) ) {
                    prefix = "Region: ";
                } else {
                  prefix = "; Region: ";
                }
                region = MemNew(StringLen(prefix) + StringLen(qvp [jdx].str) + 1);
                if ( region != NULL ) {
                    sprintf(region, "%s%s", prefix, (qvp [jdx].str));
                    FFAddString_NoRedund(unique, NULL, region, NULL, TRUE);
                    region = MemFree(region);
                } else {
                    FFAddTextToString(unique, prefix, qvp [jdx].str, NULL, FALSE, TRUE, TILDE_IGNORE);
                }
                prefix = "; ";
                add_period = FALSE;
              }
#endif
              break;

            case Qual_class_site :
              if (! StringHasNoText (qvp [jdx].str)) {
                str = qvp [jdx].str;
                if (StringCmp (str, "signal peptide") == 0 ||
                    StringCmp (str, "transit peptide") == 0 ||
                    StringCmp (str, "transmembrane region") == 0) {
                  FFAddString_NoRedund (unique, prefix, str, NULL, TRUE);
                } else {
                  FFAddString_NoRedund (unique, prefix, str, " site", TRUE);
                }
                add_period = FALSE;
                prefix = "\n";
              }
              break;

            case Qual_class_bond :
              if (! StringHasNoText (qvp [jdx].str)) {
                FFAddString_NoRedund (unique, prefix, qvp [jdx].str, " bond", TRUE);
                add_period = FALSE;
                prefix = "\n";
              }
              break;

            case Qual_class_protnames :
              /* process gene sgml for check against subsequent protein names */
              start = NULL;
              if (! StringHasNoText (qvp [FTQUAL_gene].str)) {
                /*
                if (is_journalscan) {
                  ascii_len = Sgml2AsciiLen (qvp [FTQUAL_gene].str);
                  start = ascii = MemNew ((size_t) (10 + ascii_len));
                  if (start != NULL) {
                    ascii = Sgml2Ascii (qvp [FTQUAL_gene].str, ascii, ascii_len + 1);
                  }
                } else {
                  start = StringSaveNoNull (qvp [FTQUAL_gene].str);
                }
                */
                start = StringSaveNoNull (qvp [FTQUAL_gene].str);
              }
              for (vnp = qvp [jdx].vnp; vnp != NULL; vnp = vnp->next) {
                str = (CharPtr) vnp->data.ptrvalue;
                if (! StringHasNoText (str)) {
                  /* case sensitive - gi|4973426|gb|AF148501.1|AF148501 */
                  /* check with and without sgml conversion */
                  if (StringCmp (start, str) != 0 &&
                      StringCmp (qvp [FTQUAL_gene].str, str) != 0) {
                    if (! StringStr (qvp [FTQUAL_prot_desc].str, str)) {
                      /* if (NotInGeneSyn (str, gene_syn)) { */
                        FFAddString_NoRedund (unique, prefix, str, NULL, TRUE);
                        prefix = "; ";
                        add_period = FALSE;
                      /* } */
                    }
                  }
                }
              }
              MemFree (start);
              break;

            case Qual_class_xtraprds :
              gbq = qvp [jdx].gbq;
              if (lasttype == NULL && gbq != NULL) {
                lasttype = gbq->qual;
              }
              while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
                if (! StringHasNoText (gbq->val)) {
                  if (StringCmp(gbq->val,qvp[FTQUAL_gene].str) != 0 &&
                      StringCmp(gbq->val,qvp[FTQUAL_product].str) != 0) {
                    if (!isTRNA || !StringStr (gbq->val, "RNA")) {
                      FFAddString_NoRedund (unique, prefix, gbq->val, NULL, TRUE);
                      prefix = "; ";
                      add_period = FALSE;
                    }
                  }
                }
                gbq = gbq->next;
              }
              break;

            case Qual_class_its :
              str = qvp [jdx].str;
              if (! StringHasNoText (str)) {
                if (sfp->comment == NULL || StringStr (sfp->comment, str) == NULL) {
                  FFAddString_NoRedund (unique, prefix, str, NULL, TRUE);
                  prefix = "; ";
                  add_period = FALSE;
                }
              }
              break;

            case Qual_class_trna_codons :
              trna = qvp [jdx].trp;
              if (trna) {
                numcodons = ComposeCodonsRecognizedString (trna, numbuf, sizeof (numbuf));
                if (numcodons < 1 || StringHasNoText (numbuf)) {
                } else if (numcodons == 1) {
                  isTRNA = TRUE;
                  sprintf (buf, "codon recognized: %s", numbuf);
                  if (StringStr (qvp [FTQUAL_seqfeat_note].str, buf) == NULL) {
                    FFAddString_NoRedund (unique, prefix, "codon recognized: ", numbuf, TRUE);
                    prefix = "; ";
                  }
                } else {
                  isTRNA = TRUE;
                  FFAddString_NoRedund (unique, prefix, "codons recognized: ", numbuf, TRUE);
                  prefix = "; ";
                  add_period = FALSE;
                }
              }
              break;

            case Qual_class_model_ev :
              uop = qvp [jdx].uop;
              if (uop != NULL) {
                str = NULL;
                VisitUserObjectsInUop (sfp->ext, (Pointer) &str, GetStrFormRNAEvidence);
                if (! StringHasNoText (str)) {
                  FFAddString_NoRedund (unique, prefix, str, NULL, TRUE);
                  prefix = "; ";
                  add_period = FALSE;
                }
              }
              break;

            case Qual_class_nuc_id :
              link_is_na = TRUE;
              /* fall through */
            case Qual_class_prt_id :
              sip = qvp [jdx].sip;
              if (sip != NULL) {
                /* should always be found above for protein_id or transcript_id
                prod = BioseqFind (sip);
                */
                if (prod != NULL) {
                  choice = 0;
                  for (sip = prod->id; sip != NULL; sip = sip->next) {
                    if (sip->choice == SEQID_GENBANK ||
                        sip->choice == SEQID_EMBL ||
                        sip->choice == SEQID_DDBJ ||
                        sip->choice == SEQID_OTHER ||
                        sip->choice == SEQID_TPG ||
                        sip->choice == SEQID_TPE ||
                        sip->choice == SEQID_TPD ||
                        sip->choice == SEQID_GPIPE) {
                      choice = sip->choice;
                      if (SeqIdWrite (sip, seqid, PRINTID_TEXTID_ACC_VER, sizeof (seqid)) != NULL) {
                        FFAddTextToString(unique, prefix, "transcript found in: ", seqid,
                                          FALSE, TRUE, TILDE_IGNORE);
                        prefix = "; ";
                      }
                    } else if (sip->choice == SEQID_GI) {
                      if (choice == 0) {
                        sprintf (seqid, "%ld", (long) sip->data.intvalue);
                        FFAddTextToString(unique, prefix, "transcript found in: ", seqid,
                                          FALSE, TRUE, TILDE_IGNORE);
                        prefix = "; ";
                      }
                    }
                  }
                } else {
                  if (sip->choice == SEQID_GI) {
                    gi = sip->data.intvalue;
                    if (GetAccnVerFromServer (gi, seqid)) {
                      if ((! ajp->flags.dropIllegalQuals) || ValidateAccn (seqid) == 0) {
                        FFAddTextToString(unique, prefix, "transcript found in: ", seqid,
                                          FALSE, TRUE, TILDE_IGNORE);
                        prefix = "; ";
                      }
                    } else {
                      sip = GetSeqIdForGI (gi);
                      if (sip != NULL && SeqIdWrite (sip, seqid, PRINTID_TEXTID_ACC_VER, sizeof (seqid)) != NULL) {
                        if ((! ajp->flags.dropIllegalQuals) || ValidateAccn (seqid) == 0) {
                          FFAddTextToString(unique, prefix, "transcript found in: ", seqid,
                                          FALSE, TRUE, TILDE_IGNORE);
                          prefix = "; ";
                        }
                      } else if (! ajp->flags.dropIllegalQuals) {
                        sprintf (seqid, "%ld", (long) gi);
                        FFAddTextToString(unique, prefix, "transcript found in: ", seqid,
                                          FALSE, TRUE, TILDE_IGNORE);
                        prefix = "; ";
                      }
                    }
                  } else if (SeqIdWrite (sip, seqid, PRINTID_TEXTID_ACC_VER, sizeof (seqid)) != NULL) {
                    if ((! ajp->flags.dropIllegalQuals) || ValidateAccn (seqid) == 0) {
                      FFAddTextToString(unique, prefix, "transcript found in: ", seqid,
                                          FALSE, TRUE, TILDE_IGNORE);
                      prefix = "; ";
                    }
                  }
                }
                add_period = FALSE;
              }
              break;
            default :
              break;
          }
        }

        if ( !FFEmpty(unique) ) {
          notestr = FFToCharPtr(unique);
          TrimSpacesAroundString (notestr);
          if (add_period) {
            if (! suppress_period) {
              s_AddPeriodToEnd (notestr);
            }
          }

#ifdef ASN2GNBK_STRIP_NOTE_PERIODS
          if (! IsEllipsis (notestr))
            s_RemovePeriodFromEnd (notestr);
#endif
    
          FFAddOneString(ffstring, "/note=\"", FALSE, FALSE, TILDE_IGNORE);
          FFAddOneString(ffstring, notestr, FALSE, FALSE, TILDE_SEMICOLON);
          FFAddOneString(ffstring, "\"\n", FALSE, FALSE, TILDE_IGNORE);

          MemFree (notestr);
          /*ValNodeFreeData (head);*/
        }
        break;

      default:
        break;

    }
  }
  FFRecycleString(ajp, unique);
}


NLM_EXTERN void FF_asn2gb_www_featkey (
  StringItemPtr ffstring,
  CharPtr key,
  SeqFeatPtr sfp,
  SeqLocPtr slp,
  Int4 from,
  Int4 to,
  Uint1 strand,
  Uint4 itemID
)

{
  BioseqPtr    bsp;
  Char         buf [16];
  Int4         featID = 0;
  Int4         ffrom = 0;
  Int4         fto = 0;
  Int4         gi = 0;
  Char         gi_buf[16];
  Boolean      is_aa = FALSE;
  ObjectIdPtr  oip;
  CharPtr      prefix = "?";
  SeqIntPtr    sintp;
  SeqIdPtr     sip;

  bsp = BioseqFindFromSeqLoc (slp);
  if (bsp != NULL) {
    is_aa = ISA_aa (bsp->mol);
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_GI) {
        gi = (Int4) sip->data.intvalue;
      }
    }
  } else {
    if (sfp != NULL && sfp->id.choice == 3) {
      oip = (ObjectIdPtr) sfp->id.value.ptrvalue;
      if (oip != NULL && oip->str == NULL) {
        featID = oip->id;
      }
    }
    sip = SeqLocId (slp);
    if (sip != NULL && sip->choice == SEQID_GI) {
      gi = (Int4) sip->data.intvalue;
    }
  }
  if (slp->choice == SEQLOC_INT) {
    sintp = (SeqIntPtr) slp->data.ptrvalue;
    if (sintp != NULL) {
      ffrom = sintp->from + 1;
      fto = sintp->to + 1;
      sip = sintp->id;
      if (sip->choice == SEQID_GI) {
        gi = (Int4) sip->data.intvalue;
      }
    }
  }

  sprintf (gi_buf, "%ld", (long)gi);

  if (gi > 0) {
    FFAddOneString(ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
    if (is_aa) {
      FF_Add_NCBI_Base_URL (ffstring, link_featp);
    } else {
      FF_Add_NCBI_Base_URL (ffstring, link_featn);
    }
    /* FFAddOneString(ffstring, "val=", FALSE, FALSE, TILDE_IGNORE); */
    FFAddOneString(ffstring, gi_buf, FALSE, FALSE, TILDE_IGNORE);
    if (featID > 0) {
      sprintf (buf, "%ld", (long) featID);
      FFAddOneString(ffstring, "?featID=", FALSE, FALSE, TILDE_IGNORE);
      FFAddOneString(ffstring, buf, FALSE, FALSE, TILDE_IGNORE);
      prefix = "&";
    } else if (ffrom > 0 && fto > 0) {
      sprintf (buf, "%ld", (long) ffrom);
      FFAddOneString(ffstring, "?from=", FALSE, FALSE, TILDE_IGNORE);
      FFAddOneString(ffstring, buf, FALSE, FALSE, TILDE_IGNORE);
      sprintf (buf, "%ld", (long) fto);
      FFAddOneString(ffstring, "&to=", FALSE, FALSE, TILDE_IGNORE);
      FFAddOneString(ffstring, buf, FALSE, FALSE, TILDE_IGNORE);
      prefix = "&";
    } else if (itemID > 0) {
      sprintf (buf, "%ld", (long) itemID);
      FFAddOneString(ffstring, "?itemid=", FALSE, FALSE, TILDE_IGNORE);
      FFAddOneString(ffstring, buf, FALSE, FALSE, TILDE_IGNORE);
      prefix = "&";
    }
    /*
    if ( is_aa ) {
      FFAddOneString(ffstring, prefix, FALSE, FALSE, TILDE_IGNORE);
      FFAddOneString(ffstring, "report=gpwithparts", FALSE, FALSE, TILDE_IGNORE);
    } else {
      FFAddOneString(ffstring, prefix, FALSE, FALSE, TILDE_IGNORE);
      FFAddOneString(ffstring, "report=gbwithparts", FALSE, FALSE, TILDE_IGNORE);
    }
    */
    FFAddOneString(ffstring, "\">", FALSE, FALSE, TILDE_IGNORE);
  }

  FFAddOneString(ffstring, key, FALSE, FALSE, TILDE_IGNORE);

  if (gi > 0) {
    FFAddOneString(ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
  }
}


NLM_EXTERN SeqIdPtr SeqLocIdForProduct (
  SeqLocPtr product
)

{
  SeqIdPtr   sip;
  SeqLocPtr  slp;

  /* in case product is a SEQLOC_EQUIV */

  if (product == NULL) return NULL;
  sip = SeqLocId (product);
  if (sip != NULL) return sip;
  slp = SeqLocFindNext (product, NULL);
  while (slp != NULL) {
    sip = SeqLocId (slp);
    if (sip != NULL) return sip;
    slp = SeqLocFindNext (product, slp);
  }
  return NULL;
}

NLM_EXTERN CharPtr goQualType [] = {
  "", "Process", "Component", "Function", NULL
};

static void RecordGoFieldsInQVP (
  UserFieldPtr ufp,
  Pointer userdata
)

{
  UserFieldPtr  entry;
  Int2          i;
  ObjectIdPtr   oip;
  QualValPtr    qvp;

  qvp = (QualValPtr) userdata;

  if (ufp == NULL || ufp->choice != 11) return;
  oip = ufp->label;
  if (oip == NULL) return;
  for (i = 0; goQualType [i] != NULL; i++) {
    if (StringICmp (oip->str, goQualType [i]) == 0) break;
  }
  if (goQualType [i] == NULL) return;

  entry = ufp->data.ptrvalue;
  if (entry == NULL || entry->choice != 11) return;

  /* ufp = (UserFieldPtr)  entry->data.ptrvalue; */
  switch (i) {
    case 1 :
      qvp [FTQUAL_go_process].ufp = entry;
      break;
    case 2 :
      qvp [FTQUAL_go_component].ufp = entry;
      break;
    case 3 :
      qvp [FTQUAL_go_function].ufp = entry;
      break;
    default :
      break;
  }
}

static void RecordUserObjectsInQVP (
  UserObjectPtr uop,
  Pointer userdata
)

{
  ObjectIdPtr  oip;
  QualValPtr   qvp;

  if (uop == NULL || userdata == NULL) return;
  qvp = (QualValPtr) userdata;
  oip = uop->type;
  if (oip == NULL) return;
  if (StringCmp (oip->str, "ModelEvidence") == 0) {
    qvp [FTQUAL_modelev].uop = uop;
  } else if (StringCmp (oip->str, "GeneOntology") == 0) {
    VisitUserFieldsInUop (uop, (Pointer) qvp, RecordGoFieldsInQVP);
  } else if (StringCmp (oip->str, "OfficialNomenclature") == 0) {
    qvp [FTQUAL_nomenclature].uop = uop;
  }
}

NLM_EXTERN void AddIntervalsToGbfeat (
  GBFeaturePtr gbfeat,
  SeqLocPtr location,
  BioseqPtr target
)

{
  Char           accn [41];
  SeqLocPtr      copy = NULL;
  Int4           from;
  IntFuzzPtr     fuzz;
  GBIntervalPtr  gbint;
  Int4           gi;
  Boolean        interbp;
  Boolean        iscomp;
  GBIntervalPtr  last = NULL;
  Int4           point;
  SeqIntPtr      sint;
  SeqIdPtr       sip;
  SeqLocPtr      slp;
  SeqPntPtr      spp;
  Int4           to;
  Int4           swap;

  if (gbfeat == NULL || location == NULL) return;
  if (target != NULL) {
    copy = SeqLocMerge (target, location, NULL, FALSE, TRUE, FALSE);
    location = copy;
  }

  slp = SeqLocFindNext (location, NULL);
  while (slp != NULL) {
    from = 0;
    to = 0;
    point = 0;
    iscomp = FALSE;
    interbp = FALSE;
    sip = NULL;
    switch (slp->choice) {
      case SEQLOC_WHOLE :
        sip = (SeqIdPtr) slp->data.ptrvalue;
        if (sip != NULL) {
          from = 1;
          to = SeqLocLen (slp);
          if (to < 0) {
            sip = NULL;
          }
        }
        break;
      case SEQLOC_INT :
        sint = (SeqIntPtr) slp->data.ptrvalue;
        if (sint != NULL) {
          from = sint->from + 1;
          to = sint->to + 1;
          sip = sint->id;
          if (sint->strand == Seq_strand_minus && from < to) {
            swap = from;
            from = to;
            to = swap;
          }
          if (sint->strand == Seq_strand_minus) {
            iscomp = TRUE;
          }
        }
        break;
      case SEQLOC_PNT :
        spp = (SeqPntPtr) slp->data.ptrvalue;
        if (spp != NULL) {
          point = spp->point + 1;
          sip = spp->id;
          if (spp->strand == Seq_strand_minus) {
            iscomp = TRUE;
          }
          fuzz = spp->fuzz;
          if (fuzz != NULL) {
            if (fuzz->choice == 4) {
              if (fuzz->a == 3) { /* space to right */
                from = point;
                to = point + 1;
                point = 0;
                interbp = TRUE;
              } else if (fuzz->a == 4 && point > 1) { /* space to left */
                from = point - 1;
                to = point;
                point = 0;
                interbp = TRUE;
              }
            }
          }
        }
        break;
      default :
        break;
    }
    if (sip != NULL) {
      accn [0] = '\0';
      if (sip->choice == SEQID_GI) {
        gi = sip->data.intvalue;
        if (! GetAccnVerFromServer (gi, accn)) {
          accn [0] = '\0';
        }
        if (StringHasNoText (accn)) {
          sip = GetSeqIdForGI (gi);
          SeqIdWrite (sip, accn, PRINTID_TEXTID_ACC_VER, sizeof (accn));
          SeqIdFree (sip);
        }
      } else {
        SeqIdWrite (sip, accn, PRINTID_TEXTID_ACC_VER, sizeof (accn));
      }
      if (! StringHasNoText (accn)) {
        gbint = GBIntervalNew ();
        if (gbint != NULL) {
          gbint->from = from;
          gbint->to = to;
          gbint->point = point;
          gbint->iscomp = iscomp;
          gbint->interbp = interbp;
          gbint->accession = StringSave (accn);
          if (gbfeat->intervals == NULL) {
            gbfeat->intervals = gbint;
          } else if (last != NULL) {
            last->next = gbint;
          }
          last = gbint;
        }
      }
    }
    slp = SeqLocFindNext (location, slp);
  }

  SeqLocFree (copy);
}

static CharPtr validExceptionString [] = {
  "RNA editing",
  "reasons given in citation",
  "rearrangement required for product",
  "annotated by transcript or proteomic data",
  NULL
};

static CharPtr validRefSeqExceptionString [] = {
  "alternative processing",
  "artificial frameshift",
  "nonconsensus splice site",
  "modified codon recognition",
  "alternative start codon",
  "dicistronic gene",
  "unclassified transcription discrepancy",
  "unclassified translation discrepancy",
  "mismatches in transcription",
  "mismatches in translation",
  "adjusted for low-quality genome",
  "transcribed product replaced",
  "translated product replaced",
  "transcribed pseudogene",
  /*
  "heterogeneous population sequenced",
  "low-quality sequence region",
  */
  "unextendable partial coding region",
  NULL
};

/* ribosomal slippage, trans-splicing, and artificial location now are separate qualifiers */

static void ParseException (
  CharPtr original,
  CharPtr PNTR exception_string,
  CharPtr PNTR exception_note,
  Boolean isRefSeq,
  Boolean isRelaxed,
  Uint1 subtype,
  BoolPtr riboSlipP,
  BoolPtr transSpliceP,
  BoolPtr artLocP,
  BoolPtr hetPopP,
  BoolPtr lowQualP
)

{
  ValNodePtr  excpt = NULL, note = NULL, vnp;
  Boolean     first, found;
  Int2        i;
  size_t      len;
  CharPtr     ptr, str, tmp;

  *exception_string = NULL;
  *exception_note = NULL;
  *riboSlipP = FALSE;
  *transSpliceP = FALSE;
  *artLocP = FALSE;
  *hetPopP = FALSE;
  *lowQualP = FALSE;

  if (StringHasNoText (original)) return;

  str = StringSave (original);
  if (str == NULL) return;

  tmp = str;
  while (! StringHasNoText (tmp)) {
    ptr = StringChr (tmp, ',');
    if (ptr != NULL) {
      *ptr = '\0';
      ptr++;
    }
    TrimSpacesAroundString (tmp);
    if (! StringHasNoText (tmp)) {
      found = FALSE;
      for (i = 0; validExceptionString [i] != NULL; i++) {
        if (StringICmp (tmp, validExceptionString [i]) == 0) {
          if (isRefSeq || isRelaxed || subtype == FEATDEF_CDS) {
            ValNodeCopyStr (&excpt, 0, tmp);
          } else {
            ValNodeCopyStr (&note, 0, tmp);
          }
          found = TRUE;
          break;
        }
      }
      if (! found) {
        for (i = 0; validRefSeqExceptionString [i] != NULL; i++) {
          if (StringICmp (tmp, validRefSeqExceptionString [i]) == 0) {
            if (isRefSeq || isRelaxed) {
              ValNodeCopyStr (&excpt, 0, tmp);
            } else {
              ValNodeCopyStr (&note, 0, tmp);
            }
            found = TRUE;
            break;
          }
        }
      }
      if (! found) {
        if (StringICmp (tmp, "ribosomal slippage") == 0) {
          if (subtype == FEATDEF_CDS) {
            *riboSlipP = TRUE;
          } else {
            ValNodeCopyStr (&note, 0, tmp);
          }
          found = TRUE;
        } else if (StringICmp (tmp, "trans-splicing") == 0) {
          if (subtype == FEATDEF_GENE ||
              subtype == FEATDEF_CDS ||
              subtype == FEATDEF_mRNA ||
              subtype == FEATDEF_tRNA ||
              subtype == FEATDEF_preRNA ||
              subtype == FEATDEF_otherRNA ||
              subtype == FEATDEF_3clip ||
              subtype == FEATDEF_3UTR ||
              subtype == FEATDEF_5clip ||
              subtype == FEATDEF_5UTR) {
            *transSpliceP = TRUE;
          } else {
            ValNodeCopyStr (&note, 0, tmp);
          }
          found = TRUE;
        } else if (StringICmp (tmp, "artificial location") == 0) {
          if (subtype == FEATDEF_CDS ||
              subtype == FEATDEF_mRNA) {
            *artLocP = TRUE;
          } else {
            ValNodeCopyStr (&note, 0, tmp);
          }
          found = TRUE;
        } else if (StringICmp (tmp, "heterogeneous population sequenced") == 0) {
          if (subtype == FEATDEF_CDS ||
              subtype == FEATDEF_mRNA) {
            *hetPopP = TRUE;
          } else {
            ValNodeCopyStr (&note, 0, tmp);
          }
          found = TRUE;
        } else if (StringICmp (tmp, "low-quality sequence region") == 0) {
          if (subtype == FEATDEF_CDS ||
              subtype == FEATDEF_mRNA) {
            *lowQualP = TRUE;
          } else {
            ValNodeCopyStr (&note, 0, tmp);
          }
          found = TRUE;
        }
      }
      if (! found) {
        if (isRelaxed) {
          ValNodeCopyStr (&excpt, 0, tmp);
        } else {
          ValNodeCopyStr (&note, 0, tmp);
        }
      }
    }
    tmp = ptr;
  }

  if (excpt != NULL) {
    for (vnp = excpt, len = 0; vnp != NULL; vnp = vnp->next) {
      tmp = (CharPtr) vnp->data.ptrvalue;
      len += StringLen (tmp) + 3;
    }
    ptr = (CharPtr) MemNew (len + 2);
    if (ptr != NULL) {
      for (vnp = excpt, first = TRUE; vnp != NULL; vnp = vnp->next) {
        if (! first) {
          StringCat (ptr, ", ");
        }
        tmp = (CharPtr) vnp->data.ptrvalue;
        StringCat (ptr, tmp);
        first = FALSE;
      }
    }
    *exception_string = ptr;
  }

  if (note != NULL) {
    for (vnp = note, len = 0; vnp != NULL; vnp = vnp->next) {
      tmp = (CharPtr) vnp->data.ptrvalue;
      len += StringLen (tmp) + 3;
    }
    ptr = (CharPtr) MemNew (len + 2);
    if (ptr != NULL) {
      for (vnp = note, first = TRUE; vnp != NULL; vnp = vnp->next) {
        if (! first) {
          StringCat (ptr, ", ");
        }
        tmp = (CharPtr) vnp->data.ptrvalue;
        StringCat (ptr, tmp);
        first = FALSE;
      }
    }
    *exception_note = ptr;
  }

  ValNodeFreeData (excpt);
  ValNodeFreeData (note);
  MemFree (str);
}

static CharPtr legalCategoryPrefixes [] = {
  "",
  "COORDINATES: ",
  "DESCRIPTION: ",
  "EXISTENCE: ",
  NULL
};

static CharPtr legalInferencePrefixes [] = {
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

static void ParseInference (
  GBQualPtr quals,
  ValNodePtr PNTR good_inferenceP,
  ValNodePtr PNTR bad_inferenceP
)

{
  Int2        best, j;
  ValNodePtr  good = NULL, bad = NULL;
  GBQualPtr   gbq;
  size_t      len;
  CharPtr     skip, val;

  *good_inferenceP = NULL;
  *bad_inferenceP = NULL;

  if (quals == NULL) return;

  for (gbq = quals; gbq != NULL; gbq = gbq->next) {
    if (StringICmp (gbq->qual, "inference") != 0) continue;
    if (StringHasNoText (gbq->val)) continue;
    val = gbq->val;
    skip = NULL;
    for (j = 0; legalCategoryPrefixes [j] != NULL; j++) {
      len = StringLen (legalCategoryPrefixes [j]);
      if (StringNICmp (val, legalCategoryPrefixes [j], len) != 0) continue;
      skip = val + len;
    }
    if (skip != NULL) {
      val = skip;
    }
    best = -1;
    for (j = 0; legalInferencePrefixes [j] != NULL; j++) {
      len = StringLen (legalInferencePrefixes [j]);
      if (StringNICmp (val, legalInferencePrefixes [j], len) != 0) continue;
      best = j;
    }
    if (best >= 0 && legalInferencePrefixes [best] != NULL) {
      ValNodeCopyStr (&good, 0, gbq->val);
    } else {
      ValNodeCopyStr (&bad, 0, gbq->val);
    }
  }

  *good_inferenceP = good;
  *bad_inferenceP = bad;
}

typedef struct geneprot {
  SeqFeatPtr  gene;
  SeqFeatPtr  cds;
  Boolean     failed;
} GeneProtData, PNTR GeneProtPtr;

static void CheckGeneOnIsolatedProtein (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  GeneProtPtr  gpp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_GENE) return;
  gpp = (GeneProtPtr) userdata;
  if (gpp == NULL) return;

  if (SeqLocAinB (gpp->cds->location, sfp->location) < 0) return;
  if (gpp->gene != NULL) {
    gpp->failed = TRUE;
  } else {
    gpp->gene = sfp;
  }
}

static SeqFeatPtr FindGeneOnIsolatedProtein (
  SeqEntryPtr sep,
  SeqFeatPtr cds
)

{
  GeneProtData  gpd;

  if (sep == NULL || cds == NULL) return NULL;

  MemSet ((Pointer) &gpd, 0, sizeof (GeneProtData));
  gpd.cds = cds;
  VisitFeaturesInSep (sep, (Pointer) &gpd, CheckGeneOnIsolatedProtein);

  if (gpd.failed) return NULL;

  return gpd.gene;
}

static SeqFeatPtr GetOverlappingGeneInEntity (
  Uint2 entityID,
  SeqMgrFeatContextPtr fcontext,
  SeqMgrFeatContextPtr gcontext,
  SeqLocPtr locforgene,
  Boolean is_ed,
  Boolean is_oldgb,
  IntAsn2gbJobPtr ajp
)

{
  SeqFeatPtr   gene = NULL;
  SeqEntryPtr  sep, oldscope;
  SeqInt       sint;
  SeqIntPtr    sintp;
  SeqPntPtr    spp;
  SeqPnt       spt;
  ValNode      vn;

  sep = GetTopSeqEntryForEntityID (entityID);
  oldscope = SeqEntrySetScope (sep);
  if (fcontext->featdeftype == FEATDEF_variation && locforgene != NULL) {
    /* first check same strand for variation */
    gene = SeqMgrGetOverlappingGene (locforgene, gcontext);
    if (gene == NULL) {
      /* special case variation feature - copy location but set strand both */
      if (locforgene->choice == SEQLOC_INT && locforgene->data.ptrvalue != NULL) {
        sintp = (SeqIntPtr) locforgene->data.ptrvalue;
        MemSet ((Pointer) &sint, 0, sizeof (SeqInt));
        MemSet ((Pointer) &vn, 0, sizeof (ValNode));
        sint.from = sintp->from;
        sint.to = sintp->to;
        sint.id = sintp->id;
        sint.if_from = sintp->if_from;
        sint.if_to = sintp->if_to;
        sint.strand = Seq_strand_both;
        vn.choice = SEQLOC_INT;
        vn.data.ptrvalue = (Pointer) &sint;
        gene = SeqMgrGetOverlappingGene (&vn, gcontext);

      } else if (locforgene->choice == SEQLOC_PNT && locforgene->data.ptrvalue != NULL) {
        spp = (SeqPntPtr) locforgene->data.ptrvalue;
        MemSet ((Pointer) &spt, 0, sizeof (SeqPnt));
        MemSet ((Pointer) &vn, 0, sizeof (ValNode));
        spt.point = spp->point;
        spt.id = spp->id;
        spt.fuzz = spp->fuzz;
        spt.strand = Seq_strand_both;
        vn.choice = SEQLOC_PNT;
        vn.data.ptrvalue = (Pointer) &spt;
        gene = SeqMgrGetOverlappingGene (&vn, gcontext);

      /*
      } else {
        gene = SeqMgrGetOverlappingGene (locforgene, gcontext);
      */
      }
    }
  } else {
    if (fcontext->bad_order || fcontext->mixed_strand) {
      gene = SeqMgrGetOverlappingFeatureEx (locforgene, FEATDEF_GENE, NULL, 0, NULL, LOCATION_SUBSET, gcontext, TRUE);
    } else if (ajp->multiIntervalGenes) {
      gene = SeqMgrGetOverlappingFeatureEx (locforgene, FEATDEF_GENE, NULL, 0, NULL, LOCATION_SUBSET, gcontext, TRUE);
      if (gene == NULL && (ajp->segmentedBioseqs || is_ed || is_oldgb)) {
        gene = SeqMgrGetOverlappingGene (locforgene, gcontext);
      }
    } else {
      gene = SeqMgrGetOverlappingGene (locforgene, gcontext);
    }
  }
  SeqEntrySetScope (oldscope);
  return gene;
}

static Boolean LocStrandsMatch (SeqLocPtr loc1, SeqLocPtr loc2)

{
  Uint1  featstrand;
  Uint1  locstrand;

  if (loc1 == NULL || loc2 == NULL) return FALSE;
  featstrand = SeqLocStrand (loc1);
  locstrand = SeqLocStrand (loc2);
  if (featstrand == locstrand) return TRUE;
  if (locstrand == Seq_strand_unknown && featstrand != Seq_strand_minus) return TRUE;
  if (featstrand == Seq_strand_unknown && locstrand != Seq_strand_minus) return TRUE;
  if (featstrand == Seq_strand_both && locstrand != Seq_strand_minus) return TRUE;
  if (locstrand == Seq_strand_both) return TRUE;
  return FALSE;
}

/*
static CharPtr SeqLoc2Str (
  SeqLocPtr slp
)

{
  AsnIoBSPtr    aibp;
  ByteStorePtr  bs;
  Char          ch;
  CharPtr       ptr;
  CharPtr       str;

  if (slp == NULL) return NULL;

  bs = BSNew (1000);
  if (bs == NULL) return NULL;
  aibp = AsnIoBSOpen ("w", bs);
  if (aibp == NULL) return NULL;

  SeqLocAsnWrite (slp, aibp->aip, NULL);

  AsnIoFlush (aibp->aip);
  AsnIoBSClose (aibp);

  str = BSMerge (bs, NULL);
  BSFree (bs);

  if (str == NULL) return NULL;

  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    if (ch == '\n' || ch == '\r' || ch == '\t') {
      *ptr = ' ';
    }
    ptr++;
    ch = *ptr;
  }

  TrimSpacesAndSemicolons (str);
  Asn2gnbkCompressSpaces (str);

  return str;
}
*/

static CharPtr AddJsPush (
  BioseqPtr target,
  SeqLocPtr location
)

{
  ValNodePtr  head = NULL, tail = NULL;
  IntFuzzPtr  ifp;
  SeqLocPtr   slp;
  SeqPntPtr   spp;
  Int4        start, stop;
  Char        str [64];
  CharPtr     tmp;

  if (target == NULL || location == NULL) return NULL;

  if (location->choice == SEQLOC_PNT) {
    spp = (SeqPntPtr) location->data.ptrvalue;
    if (spp != NULL) {
      ifp = spp->fuzz;
      if (ifp != NULL && ifp->choice == 4 && ifp->a == 3) {
        sprintf (str, "[[%ld, %ld]]", (long) (spp->point + 1), (long) (spp->point + 2));
        return StringSave (str);
      }
    }
  }

  slp = SeqLocFindNext (location, NULL);
  if (slp == NULL) return NULL;

  start = GetOffsetInBioseq (slp, target, SEQLOC_START) + 1;
  stop = GetOffsetInBioseq (slp, target, SEQLOC_STOP) + 1;
  if (start > 0 && stop > 0) {
    sprintf (str, "[%ld, %ld]", (long) start, (long) stop);
    ValNodeCopyStrEx (&head, &tail, 0, str);
  }

  while ((slp = SeqLocFindNext (location, slp)) != NULL) {
    start = GetOffsetInBioseq (slp, target, SEQLOC_START) + 1;
    stop = GetOffsetInBioseq (slp, target, SEQLOC_STOP) + 1;
    if (start > 0 && stop > 0) {
      sprintf (str, "[%ld, %ld]", (long) start, (long) stop);
      ValNodeCopyStrEx (&head, &tail, 0, str);
    }
  }

  tmp = ValNodeMergeStrsExEx (head, ",", "[", "]");
  ValNodeFreeData (head);

  return tmp;
}

NLM_EXTERN CharPtr AddJsInterval (
  IntAsn2gbSectPtr iasp,
  CharPtr pfx,
  BioseqPtr target,
  Uint1 featdeftype,
  SeqLocPtr location
)

{
  Char        buf [512];
  ValNodePtr  head = NULL, tail = NULL;
  CharPtr     ivls;
  CharPtr     key = NULL;
  CharPtr     tmp;

  if (iasp == NULL || target == NULL || location == NULL) return NULL;
  if (featdeftype >= FEATDEF_MAX) return NULL;

  if (StringICmp (iasp->feat_key [featdeftype], "misc_feature") == 0) {
    featdeftype = FEATDEF_misc_feature;
    if (iasp->feat_key [featdeftype] == NULL) {
      iasp->feat_key [featdeftype] = StringSave ("misc_feature");
    }
  }

  if (StringICmp (iasp->feat_key [featdeftype], "source") == 0) {
    featdeftype = FEATDEF_source;
    if (iasp->feat_key [featdeftype] == NULL) {
      iasp->feat_key [featdeftype] = StringSave ("source");
    }
  }

  key = iasp->feat_key [featdeftype];
  if (StringHasNoText (key)) return NULL;

  if (StringDoesHaveText (pfx)) {
    ValNodeCopyStrEx (&head, &tail, 0, pfx);
  }

  ValNodeCopyStrEx (&head, &tail, 0, "<script type=\"text/javascript\">");
  if (! iasp->feat_js_prefix_added) {
    sprintf (buf, "if (typeof(oData) == \"undefined\") oData = []; oData.push ({gi:%s,acc:\"%s\",features: {}});",
             iasp->gi, iasp->acc);
    ValNodeCopyStrEx (&head, &tail, 0, buf);

    iasp->feat_js_prefix_added = TRUE;
  }

  ValNodeCopyStrEx (&head, &tail, 0, "if (!oData[oData.length - 1].features[\"");
  ValNodeCopyStrEx (&head, &tail, 0, key);
  ValNodeCopyStrEx (&head, &tail, 0, "\"]) oData[oData.length - 1].features[\"");
  ValNodeCopyStrEx (&head, &tail, 0, key);
  ValNodeCopyStrEx (&head, &tail, 0, "\"] = [];");
  ValNodeCopyStrEx (&head, &tail, 0, "oData[oData.length - 1].features[\"");
  ValNodeCopyStrEx (&head, &tail, 0, key);
  ValNodeCopyStrEx (&head, &tail, 0, "\"].push(");

  ivls = AddJsPush (target, location);
  ValNodeCopyStrEx (&head, &tail, 0, ivls);

  ValNodeCopyStrEx (&head, &tail, 0, ");</script>");

  tmp = ValNodeMergeStrs (head);
  ValNodeFreeData (head);
  return tmp;
}

static CharPtr FormatFeatureBlockEx (
  Asn2gbFormatPtr afp,
  IntAsn2gbJobPtr ajp,
  Asn2gbSectPtr asp,
  BioseqPtr bsp,
  BioseqPtr target,
  SeqFeatPtr sfp,
  SeqMgrFeatContextPtr fcontext,
  QualValPtr qvp,
  FmtType format,
  IntFeatBlockPtr ifp,
  Boolean isProt,
  Boolean doKey
)

{
  Uint1              aa;
  AnnotDescrPtr      adp;
  Boolean            annotDescCommentToComment;
  Boolean            artLoc = FALSE;
  ValNodePtr         bad_inference = NULL;
  Int2               bondidx;
  BioseqPtr          bspx = NULL;
  BioseqSetPtr       bssp;
  Choice             cbaa;
  CodeBreakPtr       cbp;
  BioseqPtr          cdna;
  SeqFeatPtr         cds = NULL;
  Char               ch;
  Uint1              code = Seq_code_ncbieaa;
  CdRegionPtr        crp;
  Int4               currGi = 0;
  SeqMgrDescContext  dcontext;
  Boolean            encode_prefix = FALSE;
  CharPtr            exception_note = NULL;
  CharPtr            exception_string = NULL;
  Char               fbuf [32];
  Uint1              featdeftype;
  CharPtr            featid = NULL;
  ObjectIdPtr        fid = NULL;
  Uint1              from;
  GBQualPtr          gbq;
  GBFeaturePtr       gbfeat = NULL;
  GBSeqPtr           gbseq;
  SeqMgrFeatContext  gcontext;
  ValNodePtr         gcp;
  SeqFeatPtr         gene = NULL;
  ValNodePtr         gene_syn = NULL;
  ValNodePtr         good_inference = NULL;
  GeneRefPtr         grp = NULL;
  Boolean            hetPop = FALSE;
  IntAsn2gbSectPtr   iasp;
  IntCdsBlockPtr     icp;
  Uint2              idx;
  ValNodePtr         illegal = NULL;
  ImpFeatPtr         imp = NULL;
  IndxPtr            index;
  Boolean            is_ed = FALSE;
  Boolean            is_ged = FALSE;
  Boolean            is_gps = FALSE;
  Boolean            is_journalscan = FALSE;
  Boolean            is_oldgb = FALSE;
  Boolean            is_other = FALSE;
  Boolean            is_misc_rna = FALSE;
  Boolean            isGap = FALSE;
  Uint4              itemID;
  CharPtr            its_prod = NULL;
  CharPtr            js = NULL;
  CharPtr            key = NULL;
  CharPtr            lasttype = NULL;
  Int4               left = -1;
  SeqLocPtr          loc = NULL;
  SeqLocPtr          location = NULL;
  SeqLocPtr          locforgene = NULL;
  SeqLocPtr          locformatpep = NULL;
  Boolean            lowQual = FALSE;
  SeqMgrFeatContext  mcontext;
  MolInfoPtr         mip;
  SeqFeatPtr         mrna;
  SeqLocPtr          newloc;
  Boolean            noLeft;
  Boolean            noRight;
  SeqLocPtr          nslp = NULL;
  SeqMgrFeatContext  ocontext;
  ObjectIdPtr        oip;
  SeqEntryPtr        oldscope;
  SeqFeatPtr         operon = NULL;
  Uint2              partial;
  SeqMgrFeatContext  pcontext;
  Char               pfx [128], sfx [128];
  BioseqPtr          prd;
  CharPtr            precursor_comment = NULL;
  BioseqPtr          prod = NULL;
  SeqFeatPtr         prot;
  Boolean            protein = FALSE;
  Char               protein_pid_g [32];
  ProtRefPtr         prp;
  ProtRefPtr         prpxref;
  Boolean            pseudo = FALSE;
  CharPtr            ptr;
  Uint2              pEID;
  Int2               qualclass;
  Uint1              residue;
  RNAGenPtr          rgp;
  Boolean            riboSlippage = FALSE;
  Int4               right = -1;
  RNAQualSetPtr      rqsp;
  RnaRefPtr          rrp;
  SeqAnnotPtr        sap;
  SeqCodeTablePtr    sctp;
  SeqDescrPtr        sdp;
  SeqEntryPtr        sep;
  Uint1              seqcode;
  Uint1              seqfeattype;
  SeqIdPtr           sip;
  Int2               siteidx;
  SeqMapTablePtr     smtp;
  Boolean            split;
  CharPtr            str;
  Uint1              strand = Seq_strand_unknown;
  Boolean            suppressed = FALSE;
  CharPtr            tmp;
  Boolean            transSplice = FALSE;
  tRNAPtr            trna;
  TextSeqIdPtr       tsip;
  UserFieldPtr       ufp;
  BioseqPtr          unlockme = NULL;
  UserObjectPtr      uop;
  VariationInstPtr   vip;
  VariationRefPtr    vrp;
  VarRefDataSetPtr   vsp;
  ValNodePtr         vnp;
  SeqLocPtr          xslp = NULL;
  StringItemPtr      ffstring;
  /*
  CharPtr            firstloc = NULL;
  CharPtr            secondloc = NULL;
  CharPtr            thirdloc = NULL;
  */

  if (ajp == NULL || fcontext == NULL || qvp == NULL || ifp == NULL) return NULL;

  ffstring = FFGetString(ajp);
  if ( ffstring == NULL ) return NULL;

  if (ajp->index && asp != NULL) {
    index = &asp->index;
  } else {
    index = NULL;
  }

  if (ajp->gbseq && asp != NULL) {
    gbseq = &asp->gbseq;
  } else {
    gbseq = NULL;
  }

  pfx [0] = '\0';
  sfx [0] = '\0';

  protein_pid_g [0] = '\0';

  itemID = fcontext->itemID;

  featdeftype = fcontext->featdeftype;

  if (featdeftype < FEATDEF_GENE || featdeftype >= FEATDEF_MAX) {
    featdeftype = FEATDEF_BAD;
  }
  if (featdeftype == 0) {
    featdeftype = sfp->idx.subtype;
  }
  
  seqfeattype = fcontext->seqfeattype;
  if (seqfeattype == 0) {
    seqfeattype = sfp->data.choice;
  }


  if (doKey) {
    /* may need to map location between aa and dna */
  
    if (ifp->mapToNuc) {
  
      /* map mat_peptide, etc., to nucleotide coordinates */
  
      sip = SeqLocId (sfp->location);
      prd = BioseqFind (sip);
      cds = SeqMgrGetCDSgivenProduct (prd, NULL);
      CheckSeqLocForPartial (sfp->location, &noLeft, &noRight);
      location = aaFeatLoc_to_dnaFeatLoc (cds, sfp->location);
      SetSeqLocPartial (location, noLeft, noRight);
      /*
      locforgene = location;
      */
      if (cds != NULL) {
        grp = SeqMgrGetGeneXref (cds); /* mat_peptide first obeys any CDS gene xref */
        locformatpep = location;       /* mat_peptide next gets exact match for /gene */
        locforgene = cds->location;    /* mat_peptide last gets parent CDS /gene */
      }
      loc = location;
  
    } else if (ifp->mapToProt) {
  
      /* map CDS to protein product coordinates */
  
      sip = SeqLocIdForProduct (sfp->product);
      prd = BioseqFind (sip);
      cds = SeqMgrGetCDSgivenProduct (prd, NULL);
      location = dnaLoc_to_aaLoc (cds, sfp->location, TRUE, NULL, FALSE);
      SetSeqLocPartial (location, FALSE, FALSE);
      locforgene = sfp->location;
      loc = location;
  
    } else if (ifp->mapToGen) {
  
      /* map CDS from cDNA to genomic Bioseq */
  
      cdna = BioseqFindFromSeqLoc (sfp->location);
      mrna = SeqMgrGetRNAgivenProduct (cdna, &mcontext);
      CheckSeqLocForPartial (sfp->location, &noLeft, &noRight);
      location = productLoc_to_locationLoc (mrna, sfp->location);
      SetSeqLocPartial (location, noLeft, noRight);
      locforgene = location;
      loc = location;
  
    } else if (ifp->mapToMrna) {
  
      /* map gene from genomic to cDNA Bioseq */
  
      sep = SeqMgrGetSeqEntryForData (bsp);
      location = CreateWholeInterval (sep);
      SetSeqLocPartial (location, FALSE, FALSE);
      locforgene = location;
      loc = location;
  
    } else if (ifp->mapToPep) {
  
      /* map protein processing from precursor to subpeptide Bioseq */
  
      sep = SeqMgrGetSeqEntryForData (bsp);
      location = CreateWholeInterval (sep);
      SetSeqLocPartial (location, FALSE, FALSE);
      locforgene = location;
      loc = location;
  
    } else {
  
      /* no aa-dna or dna-aa mapping, just use location */
  
      location = sfp->location;
      locforgene = sfp->location;
    }
    if (location == NULL) return NULL;

    if (loc != NULL) {
      NormalizeNullsBetween (loc);
    }

    sep = GetTopSeqEntryForEntityID (ajp->ajp.entityID);
    if (sep != NULL && IS_Bioseq_set (sep)) {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      if (bssp != NULL && bssp->_class == BioseqseqSet_class_gen_prod_set) {
        is_gps = TRUE;
      }
    }

    if (bsp != NULL) {
      for (sip = bsp->id; sip != NULL; sip = sip->next) {
        switch (sip->choice) {
          case SEQID_OTHER :
            is_other = TRUE;
            break;
          case SEQID_GIBBSQ :
          case SEQID_GIBBMT :
          case SEQID_GIIM :
            is_journalscan = TRUE;
            break;
          case SEQID_GENBANK :
          case SEQID_TPG :
            is_ged = TRUE;
            tsip = (TextSeqIdPtr) sip->data.ptrvalue;
            if (tsip != NULL) {
              if (StringLen (tsip->accession) == 6) {
                is_oldgb = TRUE;
              }
            }
            break;
          case SEQID_EMBL :
          case SEQID_TPE :
            is_ged = TRUE;
            is_ed = TRUE;
            break;
          case SEQID_DDBJ :
          case SEQID_TPD :
            is_ged = TRUE;
            is_ed = TRUE;
            break;
          default :
          break;
        }
      }
    }
 
    if (ajp->refseqConventions) {
      is_other = TRUE;
    }

    key = FindKeyFromFeatDefType (featdeftype, TRUE);
    isGap = (Boolean) (featdeftype == FEATDEF_gap);
  
    if (format == GENPEPT_FMT && isProt) {
      if (featdeftype == FEATDEF_REGION) {
        key = "Region";
      } else if (featdeftype == FEATDEF_BOND) {
        key = "Bond";
      } else if (featdeftype == FEATDEF_SITE) {
        key = "Site";
      }
      if (ifp->mapToPep) {
        if (featdeftype >= FEATDEF_preprotein && featdeftype <= FEATDEF_transit_peptide_aa) {
          key = "Precursor";
          itemID = 0;
        }
      }
    }
    if (! isProt) {
      if (featdeftype == FEATDEF_preprotein) {
        if (! is_other) {
          key = "misc_feature";
          encode_prefix = TRUE;
        }
      }
    }
    if (featdeftype == FEATDEF_CLONEREF) {
      if (ajp->mode == RELEASE_MODE || ajp->mode == ENTREZ_MODE) {
        key = "misc_feature";
      }
    }
    /*
    if (featdeftype == FEATDEF_VARIATIONREF) {
      if (ajp->mode == RELEASE_MODE || ajp->mode == ENTREZ_MODE) {
        key = "misc_feature";
      }
    }
    */
  
    /* deal with unmappable impfeats */
  
    if (featdeftype == FEATDEF_BAD && seqfeattype == SEQFEAT_IMP) {
      imp = (ImpFeatPtr) sfp->data.value.ptrvalue;
      if (imp != NULL) {
        key = imp->key;
      }
    }

    /* prior to BSEC conversion, map for release and web, allow old feature to be seen in Sequin */
    if (featdeftype == FEATDEF_repeat_unit && (ajp->mode == RELEASE_MODE || ajp->mode == ENTREZ_MODE)) {
      key = "repeat_region";
    }

    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_GI) {
        currGi = (Int4) sip->data.intvalue;
      }
    }

    iasp = (IntAsn2gbSectPtr) asp;

    if (iasp != NULL && GetWWW (ajp) && ajp->mode == ENTREZ_MODE && featdeftype < FEATDEF_MAX) {
      if (iasp->feat_key [featdeftype] == NULL) {
        iasp->feat_key [featdeftype] = StringSave (key);
      }
    }

    if( isGap ) {
        /* look through quals for any that might indicate we
           should be an assembly_gap */
        for( gbq = sfp->qual; gbq != NULL; gbq = gbq->next ) {
            if( 0 == StrCmp(gbq->qual, "gap_type") || 
                0 == StrCmp(gbq->qual, "linkage_evidence") ) 
            {
                key = "assembly_gap";
                break;
            }
        }
    }


    if (afp != NULL && GetWWW (ajp) && ajp->mode == ENTREZ_MODE && ajp->seqspans &&
        (ajp->format == GENBANK_FMT || ajp->format == GENPEPT_FMT)) {
      sprintf (pfx, "<span id=\"feature_%ld_%s_%ld\" class=\"feature\">", (long) currGi, key, (long) ifp->feat_count);
    }

    FFStartPrint(ffstring, format, 5, 21, NULL, 0, 5, 21, "FT", /* ifp->firstfeat */ FALSE);

    if (ajp->ajp.slp != NULL) {
      FFAddOneString(ffstring, key, FALSE, FALSE, TILDE_IGNORE);
    } else if ( 
        GetWWW(ajp) &&
        StringICmp (key, "gap") != 0 &&
        StringICmp (key, "assembly_gap") != 0 &&
        bsp != NULL /* && SeqMgrGetParentOfPart (bsp, NULL) == NULL */ ) {
      FF_asn2gb_www_featkey (ffstring, key, sfp, sfp->location, fcontext->left + 1, fcontext->right + 1,
                             fcontext->strand, itemID);
    } else {
      FFAddOneString(ffstring, key, FALSE, FALSE, TILDE_IGNORE);
    }
    FFAddNChar(ffstring, ' ', 21 - 5 - StringLen(key), FALSE);
  
    if (gbseq != NULL) {
      gbfeat = GBFeatureNew ();
      if (gbfeat != NULL) {
        gbfeat->key = StringSave (key);
      }
    }
  
    if (imp == NULL || StringHasNoText (imp->loc)) {
  
      if (ajp->ajp.slp != NULL) {
        sip = SeqIdParse ("lcl|dummy");
        left = GetOffsetInBioseq (ajp->ajp.slp, bsp, SEQLOC_LEFT_END);
        right = GetOffsetInBioseq (ajp->ajp.slp, bsp, SEQLOC_RIGHT_END);
        strand = SeqLocStrand (ajp->ajp.slp);
        split = FALSE;
        newloc = SeqLocReMapEx (sip, ajp->ajp.slp, location, 0, FALSE, ajp->masterStyle);
        /*
        newloc = SeqLocCopyRegion (sip, location, bsp, left, right, strand, &split);
        */
        SeqIdFree (sip);
        if (newloc == NULL) {
          return NULL;
        }
        /*
        firstloc = SeqLoc2Str (newloc);
        */
        A2GBSeqLocReplaceID (newloc, ajp->ajp.slp);
        /*
        secondloc = SeqLoc2Str (newloc);
        */
        str = FFFlatLoc (ajp, target, newloc, ajp->masterStyle, isGap);
        if (iasp != NULL && GetWWW (ajp) && ajp->mode == ENTREZ_MODE && ajp->seqspans && featdeftype < FEATDEF_MAX) {
          js = AddJsInterval (iasp, pfx, target, featdeftype, newloc);
        }
        SeqLocFree (newloc);
        /*
        thirdloc = SeqLoc2Str (ajp->ajp.slp);
        if (StringCmp (str, "?") != 0) {
          firstloc = MemFree (firstloc);
          secondloc = MemFree (secondloc);
          thirdloc = MemFree (thirdloc);
        }
        */
      } else {
        str = FFFlatLoc (ajp, target, location, ajp->masterStyle, isGap);
        if (iasp != NULL && GetWWW (ajp) && ajp->mode == ENTREZ_MODE && ajp->seqspans && featdeftype < FEATDEF_MAX) {
          js = AddJsInterval (iasp, pfx, target, featdeftype, location);
        }
        /*
        if (StringCmp (str, "?") == 0) {
          firstloc = SeqLoc2Str (location);
          SeqIdWrite (target->id, buf, PRINTID_FASTA_LONG, sizeof (buf));
          secondloc = StringSave (buf);
          thirdloc = NULL;
        }
        */
      }
    } else {
      str = StringSave (imp->loc);
    }
    if ( GetWWW(ajp) ) {
      FF_www_featloc (ffstring, str);
    } else {
      FFAddOneString(ffstring, str, FALSE, FALSE, TILDE_IGNORE);
    }
  
    if (gbseq != NULL) {
      if (gbfeat != NULL) {
        gbfeat->location = StringSave (str);
        if (StringDoesHaveText (str)) {
          if (StringStr (str, "join") != NULL) {
            gbfeat->operator__ = StringSave ("join");
          } else if (StringStr (str, "order") != NULL) {
            gbfeat->operator__ = StringSave ("order");
          }
        }
        gbfeat->partial5 = fcontext->partialL;
        gbfeat->partial3 = fcontext->partialR;
        if (ajp->masterStyle) {
          AddIntervalsToGbfeat (gbfeat, location, target);
        } else {
          AddIntervalsToGbfeat (gbfeat, location, NULL);
        }
      }
    }
  
    MemFree (str);

  } else {

    location = sfp->location;
    locforgene = sfp->location;
  }

  if (location != NULL) {
   ifp->left = GetOffsetInBioseq (location, bsp, SEQLOC_LEFT_END);
   ifp->right = GetOffsetInBioseq (location, bsp, SEQLOC_RIGHT_END);
  }

  /* populate qualifier table from feature fields */

  /*
  if (sfp->partial == TRUE)
    sfp->partial = FlatAnnotPartial(sfp, use_product);
  */

  if (sfp->partial) {
    partial = SeqLocPartialCheck (location);
    if (partial == SLP_COMPLETE /* || partial > SLP_OTHER */ ) {
      qvp [FTQUAL_partial].ble = TRUE;
    }
    if (LookForFuzz (location)) {
      qvp [FTQUAL_partial].ble = FALSE;
    }
    if (imp != NULL) {
      if (StringChr (imp->loc, '<') != NULL || StringChr (imp->loc, '>') != NULL) {
        qvp [FTQUAL_partial].ble = FALSE;
      }
    }

    /* hide unclassified /partial in RELEASE_MODE and ENTREZ_MODE */

    if (ajp->mode == RELEASE_MODE || ajp->mode == ENTREZ_MODE) {
      qvp [FTQUAL_partial].ble = FALSE;
    }
    /*
    if (ajp->flags.checkQualSyntax) {
      switch (featdeftype) {
      case FEATDEF_conflict:
      case FEATDEF_mutation:
      case FEATDEF_N_region:
      case FEATDEF_polyA_site:
        qvp [FTQUAL_partial].ble = FALSE;
        break;
      default:
        break;
      }
    }
    */
  }
  if (ifp->mapToProt) {
    qvp [FTQUAL_partial].ble = FALSE;
  }

  if (sfp->pseudo) {
    pseudo = TRUE;
  }

  if (seqfeattype == SEQFEAT_GENE) {
    grp = (GeneRefPtr) sfp->data.value.ptrvalue;
    if (grp != NULL) {
      if (! StringHasNoText (grp->locus)) {
        qvp [FTQUAL_gene].str = grp->locus;
        qvp [FTQUAL_locus_tag].str = grp->locus_tag;
        qvp [FTQUAL_gene_desc].str = grp->desc;
        qvp [FTQUAL_gene_syn].vnp = grp->syn;
      } else if (grp->locus_tag != NULL) {
        qvp [FTQUAL_locus_tag].str = grp->locus_tag;
        qvp [FTQUAL_gene_desc].str = grp->desc;
        qvp [FTQUAL_gene_syn].vnp = grp->syn;
      } else if (! StringHasNoText (grp->desc)) {
        qvp [FTQUAL_gene].str = grp->desc;
        qvp [FTQUAL_gene_syn].vnp = grp->syn;
      } else if (grp->syn != NULL) {
        vnp = grp->syn;
        qvp [FTQUAL_gene].str = (CharPtr) vnp->data.ptrvalue;
        vnp = vnp->next;
        qvp [FTQUAL_gene_syn].vnp = vnp;
      }
      qvp [FTQUAL_gene_map].str = grp->maploc;
      qvp [FTQUAL_gene_allele].str = grp->allele;
      qvp [FTQUAL_gene_xref].vnp = grp->db;
      if (grp->pseudo) {
        pseudo = TRUE;
      }
      qvp [FTQUAL_gene_nomen].gnp = grp->formal_name;
    }
    if (! ajp->flags.separateGeneSyns) {
      qvp [FTQUAL_gene_syn_refseq].vnp = qvp [FTQUAL_gene_syn].vnp;
      qvp [FTQUAL_gene_syn].vnp = NULL;
    }
    operon = SeqMgrGetOverlappingOperon (locforgene, &ocontext);
    if (operon != NULL) {
      for (gbq = operon->qual; gbq != NULL; gbq = gbq->next) {
        if (StringCmp (gbq->qual, "operon") == 0) {
          qvp [FTQUAL_operon].gbq = gbq;
        }
      }
      if (operon->pseudo) {
        pseudo = TRUE;
      }
    }

  } else if (featdeftype != FEATDEF_operon && featdeftype != FEATDEF_gap) {

    /* if mat_peptide, grp is already be set based on parent CDS, otherwise check current feature */

    if (grp == NULL) {
      grp = SeqMgrGetGeneXrefEx (sfp, &fid);
      if (fid != NULL) {
        if (StringDoesHaveText (fid->str)) {
          featid = fid->str;
        } else {
          sprintf (fbuf, "%ld", (long) fid->id);
          featid = fbuf;
        }
      }
    }

    /* if gene xref, then find referenced gene, take everything as if it overlapped */

    if (grp != NULL) {
      if (SeqMgrGeneIsSuppressed (grp)) {
        suppressed = TRUE;
      } else {
        if (grp->pseudo) {
          pseudo = TRUE;
        }
        bspx = BioseqFindFromSeqLoc (sfp->location);
        if (bspx != NULL) {
          if (featid != NULL) {
            gene = SeqMgrGetFeatureByFeatID (0, bspx, featid, NULL, &gcontext);
          } else if (StringDoesHaveText (grp->locus_tag)) {
            gene = SeqMgrGetGeneByLocusTag (bspx, grp->locus_tag, &gcontext);
          } else if (StringDoesHaveText (grp->locus)) {
            gene = SeqMgrGetFeatureByLabel (bspx, grp->locus, SEQFEAT_GENE, 0, &gcontext);
          }
          if (gene != NULL) {
            grp = (GeneRefPtr) gene->data.value.ptrvalue;
            if (gene->pseudo) {
              pseudo = TRUE;
            }
            if (grp != NULL && grp->db != NULL) {
              qvp [FTQUAL_gene_xref].vnp = grp->db;
            } else {
              qvp [FTQUAL_gene_xref].vnp = gene->dbxref;
            }
          }
        }
      }
    }

    if (! suppressed) {

      /* first look for gene that exactly matches mat_peptide DNA projection */

      if (gene == NULL && grp == NULL && locformatpep != NULL) {
        gene = GetOverlappingGeneInEntity (ajp->ajp.entityID, fcontext, &gcontext, locformatpep, is_ed, is_oldgb, ajp);
        if (gene == NULL && ajp->ajp.entityID != sfp->idx.entityID) {
          gene = GetOverlappingGeneInEntity (sfp->idx.entityID, fcontext, &gcontext, locformatpep, is_ed, is_oldgb, ajp);
        }

        if (gene != NULL) {
          xslp = AsnIoMemCopy (gene->location, (AsnReadFunc) SeqLocAsnRead, (AsnWriteFunc) SeqLocAsnWrite);
          if (xslp != NULL) {
            if (xslp->choice == SEQLOC_MIX) {
              nslp = ValNodeExtractList ((ValNodePtr PNTR) &(xslp->data.ptrvalue), SEQLOC_NULL);
            }
            if (SeqLocCompare (xslp, locformatpep) == SLC_A_EQ_B &&
                LocStrandsMatch (xslp, locformatpep)) {
              qvp [FTQUAL_gene_note].str = gene->comment;

              grp = (GeneRefPtr) gene->data.value.ptrvalue;
              if (gene->pseudo) {
                pseudo = TRUE;
              }
              if (grp != NULL && grp->db != NULL) {
                qvp [FTQUAL_gene_xref].vnp = grp->db;
              } else {
                qvp [FTQUAL_gene_xref].vnp = gene->dbxref;
              }
            } else {
              gene = NULL;
            }
          } else {
            gene = NULL;
          }
          SeqLocFree (xslp);
          SeqLocFree (nslp);
        }
      }

      /* otherwise, if not suppressed and no gene xref, get gene by overlap */

      if (gene == NULL && grp == NULL) {
        if (featdeftype != FEATDEF_primer_bind) {
          gene = GetOverlappingGeneInEntity (ajp->ajp.entityID, fcontext, &gcontext, locforgene, is_ed, is_oldgb, ajp);
          if (gene == NULL && ajp->ajp.entityID != sfp->idx.entityID) {
            gene = GetOverlappingGeneInEntity (sfp->idx.entityID, fcontext, &gcontext, locforgene, is_ed, is_oldgb, ajp);
          }
        }

        /* special case to get gene by overlap for coded_by cds on isolated protein bioseq */
        if (ifp->mapToProt && seqfeattype == SEQFEAT_CDREGION) {
          sep = GetTopSeqEntryForEntityID (ajp->ajp.entityID);
          if (sep != NULL && IS_Bioseq (sep)) {
            bspx = (BioseqPtr) sep->data.ptrvalue;
            if (bspx != NULL && ISA_aa (bspx->mol)) {
              gene = FindGeneOnIsolatedProtein (sep, sfp);
            }
          }
        }

        if (gene != NULL) {
          qvp [FTQUAL_gene_note].str = gene->comment;

          grp = (GeneRefPtr) gene->data.value.ptrvalue;
          if (gene->pseudo) {
            pseudo = TRUE;
          }
          if (grp != NULL && grp->db != NULL) {
            qvp [FTQUAL_gene_xref].vnp = grp->db;
          } else {
            qvp [FTQUAL_gene_xref].vnp = gene->dbxref;
          }
        }
      }

      if (grp != NULL && grp->pseudo) {
        pseudo = TRUE;
      }

      if (grp != NULL &&
          ((featdeftype != FEATDEF_repeat_region && featdeftype != FEATDEF_mobile_element) ||
           is_ed || gene == NULL)) {
        if (! StringHasNoText (grp->locus)) {
          qvp [FTQUAL_gene].str = grp->locus;
          qvp [FTQUAL_locus_tag].str = grp->locus_tag;
          qvp [FTQUAL_gene_syn].vnp = grp->syn;
          gene_syn = grp->syn;
        } else if (! StringHasNoText (grp->locus_tag)) {
          qvp [FTQUAL_locus_tag].str = grp->locus_tag;
          qvp [FTQUAL_gene_syn].vnp = grp->syn;
          gene_syn = grp->syn;
        } else if (! StringHasNoText (grp->desc)) {
          qvp [FTQUAL_gene].str = grp->desc;
          qvp [FTQUAL_gene_syn].vnp = grp->syn;
          gene_syn = grp->syn;
        } else if (grp->syn != NULL) {
          vnp = grp->syn;
          qvp [FTQUAL_gene].str = (CharPtr) vnp->data.ptrvalue;
          vnp = vnp->next;
          qvp [FTQUAL_gene_syn].vnp = vnp;
          gene_syn = vnp;
        }
      }
      if (! ajp->flags.separateGeneSyns) {
        qvp [FTQUAL_gene_syn_refseq].vnp = qvp [FTQUAL_gene_syn].vnp;
        qvp [FTQUAL_gene_syn].vnp = NULL;
      }
      if (grp != NULL &&
          featdeftype != FEATDEF_variation &&
          ((featdeftype != FEATDEF_repeat_region && featdeftype != FEATDEF_mobile_element) || is_ed)) {
        qvp [FTQUAL_gene_allele].str = grp->allele; /* now propagating /allele */
      }

      if (gene != NULL && ((featdeftype != FEATDEF_repeat_region && featdeftype != FEATDEF_mobile_element) || is_ed)) {
        /* now propagate old_locus_tag to almost any underlying feature */
        for (gbq = gene->qual; gbq != NULL; gbq = gbq->next) {
          if (StringHasNoText (gbq->val)) continue;
          idx = GbqualToFeaturIndex (gbq->qual);
          if (idx == FTQUAL_old_locus_tag) {
            qvp [FTQUAL_old_locus_tag].gbq = gbq;
            break; /* record first old_locus_tag gbqual to display all */
          }
        }
      }
      if (seqfeattype != SEQFEAT_CDREGION && seqfeattype != SEQFEAT_RNA) {
        qvp [FTQUAL_gene_xref].vnp = NULL;
      }

      if (featdeftype != FEATDEF_operon) {
        operon = SeqMgrGetOverlappingOperon (locforgene, &ocontext);
        if (operon != NULL) {
          for (gbq = operon->qual; gbq != NULL; gbq = gbq->next) {
            if (StringCmp (gbq->qual, "operon") == 0) {
              qvp [FTQUAL_operon].gbq = gbq;
            }
          }
          if (operon->pseudo) {
            pseudo = TRUE;
          }
        }
      }
    }

    /* specific fields set here */

    switch (seqfeattype) {
      case SEQFEAT_CDREGION :
        if (! ifp->mapToProt) {
          crp = (CdRegionPtr) sfp->data.value.ptrvalue;
          if (crp != NULL) {

            qvp [FTQUAL_codon_start].num = crp->frame;
            if (qvp [FTQUAL_codon_start].num == 0) {
              qvp [FTQUAL_codon_start].num = 1;
            }
            qvp [FTQUAL_transl_except].cbp = crp->code_break;
            for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) {
              seqcode = 0;
              sctp = NULL;
              cbaa = cbp->aa;
              switch (cbaa.choice) {
                case 1 :
                  seqcode = Seq_code_ncbieaa;
                  break;
                case 2 :
                  seqcode = Seq_code_ncbi8aa;
                  break;
                case 3 :
                  seqcode = Seq_code_ncbistdaa;
                  break;
                default :
                  break;
              }
              if (seqcode != 0) {
                sctp = SeqCodeTableFind (seqcode);
                if (sctp != NULL) {
                  residue = cbaa.value.intvalue;
                  if (residue != 42) {
                    if (seqcode != Seq_code_ncbieaa) {
                      smtp = SeqMapTableFind (seqcode, Seq_code_ncbieaa);
                      residue = SeqMapTableConvert (smtp, residue);
                    }
                    /*
                    if (residue == 'U') {
                      if (ajp->flags.selenocysteineToNote) {
                        qvp [FTQUAL_selenocysteine_note].str = "selenocysteine";
                      } else {
                        qvp [FTQUAL_selenocysteine].ble = TRUE;
                      }
                    } else if (residue == 'O') {
                      if (ajp->flags.pyrrolysineToNote) {
                        qvp [FTQUAL_pyrrolysine_note].str = "pyrrolysine";
                      } else {
                        qvp [FTQUAL_pyrrolysine].ble = TRUE;
                      }
                    }
                    */
                  }
                }
              }
            }

            gcp = crp->genetic_code;
            if (gcp != NULL) {
              for (vnp = gcp->data.ptrvalue; vnp != NULL; vnp = vnp->next) {
                if (vnp->choice == 2 && vnp->data.intvalue != 0) {
                  qvp [FTQUAL_transl_table].num = vnp->data.intvalue;
                }
              }

              /* suppress table 1, but always show it in GBSeqXML */

              if (qvp [FTQUAL_transl_table].num == 1 && ajp->gbseq == NULL) {
                qvp [FTQUAL_transl_table].num = 0;
              }
            }

            if (sfp->product != NULL && SeqLocLen (sfp->product) != 0) {
              protein = TRUE;
            }
            if (crp->conflict && (protein || (! sfp->excpt))) {
              if (protein) {
                qvp [FTQUAL_prot_conflict].str = conflict_msg;
              } else {
                /*
                qvp [FTQUAL_prot_missing].str = no_protein_msg;
                */
              }
            }
          }

          sip = SeqLocIdForProduct (sfp->product);
          qvp [FTQUAL_protein_id].sip = sip;

          sep = GetTopSeqEntryForEntityID (ajp->ajp.entityID);

          if (! ajp->alwaysTranslCds) {

            /* by default only show /translation if product bioseq is within entity */

            oldscope = SeqEntrySetScope (sep);
            prod = BioseqFind (sip);
            SeqEntrySetScope (oldscope);

            if (prod == NULL && ajp->showFarTransl) {

              /* but flag can override and force far /translation */

              prod = BioseqLockById (sip);
              unlockme = prod;
            }
          }

          prp = NULL;

          if (prod != NULL) {
            for (sip = prod->id; sip != NULL; sip = sip->next) {
              if (sip->choice == SEQID_GI) {
                sprintf (protein_pid_g, "PID:g%ld", (long) sip->data.intvalue);
              }
            }
            sdp = SeqMgrGetNextDescriptor (prod, NULL, Seq_descr_comment, &dcontext);
            if (sdp != NULL && dcontext.level == 0) {
              if (! StringHasNoText ((CharPtr) sdp->data.ptrvalue)) {
                qvp [FTQUAL_prot_comment].str = (CharPtr) sdp->data.ptrvalue;
              }
            }
            sdp = SeqMgrGetNextDescriptor (prod, NULL, Seq_descr_molinfo, &dcontext);
            if (sdp != NULL && dcontext.level == 0) {
              mip = (MolInfoPtr) sdp->data.ptrvalue;
              if (mip != NULL && mip->tech > 1 &&
                  mip->tech != MI_TECH_concept_trans &&
                  mip->tech != MI_TECH_concept_trans_a) {
                str = StringForSeqTech (mip->tech);
                if (! StringHasNoText (str)) {
                  qvp [FTQUAL_prot_method].str = str;
                }
              }
            }
            pEID = ObjMgrGetEntityIDForPointer (prod);
            if (pEID != 0 && pEID != ajp->ajp.entityID &&
                SeqMgrFeaturesAreIndexed (pEID) == 0) {
              /* index far record so SeqMgrGetBestProteinFeature can work */
              SeqMgrIndexFeatures (pEID, NULL);
            }
            prot = SeqMgrGetBestProteinFeature (prod, &pcontext);
            if (prot != NULL) {
              prp = (ProtRefPtr) prot->data.value.ptrvalue;
              if (prp != NULL && prp->processed < 2) {
                qvp [FTQUAL_prot_note].str = prot->comment;
              }
            }
          }

          /* protein xref overrides names, but should not prevent /protein_id, etc. */

          prpxref = SeqMgrGetProtXref (sfp);
          if (prpxref != NULL) {
            prp = prpxref;
          }
          if (prp != NULL) {
            vnp = prp->name;
            if (vnp != NULL && (! StringHasNoText ((CharPtr) vnp->data.ptrvalue))) {
              qvp [FTQUAL_cds_product].str = (CharPtr) vnp->data.ptrvalue;
              vnp = vnp->next;
              if (ajp->flags.extraProductsToNote) {
                qvp [FTQUAL_prot_names].vnp = vnp;
              } else {
                qvp [FTQUAL_extra_products].vnp = vnp;
              }
            }
            qvp [FTQUAL_prot_desc].str = prp->desc;
            qvp [FTQUAL_prot_activity].vnp = prp->activity;
            qvp [FTQUAL_prot_EC_number].vnp = prp->ec;
          }

          if (! pseudo) {
            if (prod != NULL || ajp->transIfNoProd || ajp->alwaysTranslCds) {
              if (doKey) {
                if (! ajp->hideTranslation) {
                  qvp [FTQUAL_translation].ble = TRUE;
                }
              }
            }
          }

          if (ifp->isCDS) {
            icp = (IntCdsBlockPtr) ifp;
            qvp [FTQUAL_figure].str = icp->fig;
            qvp [FTQUAL_maploc].str = icp->maploc;
          }
        } else {
          qvp [FTQUAL_coded_by].slp = sfp->location;

          /* show /evidence on coded_by CDS */

          qvp [FTQUAL_evidence].num = sfp->exp_ev;

          crp = (CdRegionPtr) sfp->data.value.ptrvalue;
          if (crp != NULL) {
            if (crp->frame > 1) {
              qvp [FTQUAL_codon_start].num = crp->frame;
            }
            gcp = crp->genetic_code;
            if (gcp != NULL) {
              for (vnp = gcp->data.ptrvalue; vnp != NULL; vnp = vnp->next) {
                if (vnp->choice == 2 && vnp->data.intvalue != 0) {
                  qvp [FTQUAL_transl_table].num = vnp->data.intvalue;
                }
              }

              /* suppress table 1 */

              if (qvp [FTQUAL_transl_table].num == 1) {
                qvp [FTQUAL_transl_table].num = 0;
              }
            }
            for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) {
              seqcode = 0;
              sctp = NULL;
              cbaa = cbp->aa;
              switch (cbaa.choice) {
                case 1 :
                  seqcode = Seq_code_ncbieaa;
                  break;
                case 2 :
                  seqcode = Seq_code_ncbi8aa;
                  break;
                case 3 :
                  seqcode = Seq_code_ncbistdaa;
                  break;
                default :
                  break;
              }
              if (seqcode != 0) {
                sctp = SeqCodeTableFind (seqcode);
                if (sctp != NULL) {
                  residue = cbaa.value.intvalue;
                  if (residue != 42) {
                    if (seqcode != Seq_code_ncbieaa) {
                      smtp = SeqMapTableFind (seqcode, Seq_code_ncbieaa);
                      residue = SeqMapTableConvert (smtp, residue);
                    }
                    /*
                    if (residue == 'U') {
                      if (ajp->flags.selenocysteineToNote) {
                        qvp [FTQUAL_selenocysteine_note].str = "selenocysteine";
                      } else {
                        qvp [FTQUAL_selenocysteine].ble = TRUE;
                      }
                    } else if (residue == 'O') {
                      if (ajp->flags.pyrrolysineToNote) {
                        qvp [FTQUAL_pyrrolysine_note].str = "pyrrolysine";
                      } else {
                        qvp [FTQUAL_pyrrolysine].ble = TRUE;
                      }
                    }
                    */
                  }
                }
              }
            }
          }
        }
        break;
      case SEQFEAT_PROT :
        if (! ifp->mapToPep) {
          prp = (ProtRefPtr) sfp->data.value.ptrvalue;
          if (prp != NULL) {
            vnp = prp->name;
            if (vnp != NULL && (! StringHasNoText ((CharPtr) vnp->data.ptrvalue))) {
              qvp [FTQUAL_product].str = (CharPtr) vnp->data.ptrvalue;
              vnp = vnp->next;
              qvp [FTQUAL_prot_names].vnp = vnp;
            }
            if (format != GENPEPT_FMT) {
              qvp [FTQUAL_prot_desc].str = prp->desc;
            } else {
              qvp [FTQUAL_prot_name].str = prp->desc;
            }
            if (format != GENPEPT_FMT || prp->processed != 2) {
              qvp [FTQUAL_prot_activity].vnp = prp->activity;
            }
            qvp [FTQUAL_prot_EC_number].vnp = prp->ec;
          }
          sip = SeqLocIdForProduct (sfp->product);
          if (sip != NULL) {
            /* for RefSeq records or GenBank not release_mode */
            if (is_other || (! ajp->flags.forGbRelease)) {
              qvp [FTQUAL_protein_id].sip = sip;
            }
            prod = BioseqFind (sip);
          }
        } else {
          qvp [FTQUAL_derived_from].slp = sfp->location;
          sip = SeqLocIdForProduct (sfp->product);
          if (sip != NULL) {
            prod = BioseqFind (sip);
            if (prod != NULL) {
              prot = SeqMgrGetBestProteinFeature (prod, NULL);
              if (prot != NULL) {
                precursor_comment = prot->comment;
              }
            }
          }
        }
        prp = (ProtRefPtr) sfp->data.value.ptrvalue;
        if (prp != NULL) {
          if (! pseudo) {
            if (ajp->showPeptide) {
              if (prp->processed == 2 || prp->processed == 3 || prp->processed == 4) {
                qvp [FTQUAL_peptide].ble = TRUE;
              }
            }
            if (format == GENPEPT_FMT && isProt && is_other) {
              /* enable calculated_mol_wt qualifier for RefSeq proteins */
              qvp [FTQUAL_mol_wt].ble = TRUE;
            }
          }
          if (prp->processed == 3 || prp->processed == 4) {
            if (! is_other) {
              /* Only RefSeq allows product on signal or transit peptide */
              qvp [FTQUAL_product].str = NULL;
            }
          }
          if (prp->processed == 1 && encode_prefix && (! is_other)) {
            qvp [FTQUAL_encodes].str = qvp [FTQUAL_product].str;
            qvp [FTQUAL_product].str = NULL;
          }
        }
        break;
      case SEQFEAT_RNA :
        rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
        if (rrp != NULL) {
          if (rrp->pseudo) {
            pseudo = TRUE;
          }
          sip = SeqLocIdForProduct (sfp->product);
          if (sip != NULL) {
            /* for RefSeq records or GenBank not release_mode or entrez_mode */
            if (is_other || (ajp->mode == SEQUIN_MODE || ajp->mode == DUMP_MODE)) {
              qvp [FTQUAL_transcript_id].sip = sip;
            } else {
              /* otherwise now goes in note */
              qvp [FTQUAL_transcript_id_note].sip = sip; /* !!! remove October 15, 2003 !!! */
            }

            if (! ajp->alwaysTranslCds) {

              /* by default only show /transcription if product bioseq is within entity */

              oldscope = SeqEntrySetScope (sep);
              prod = BioseqFind (sip);
              SeqEntrySetScope (oldscope);

              if (prod == NULL && ajp->showFarTransl) {

                /* but flag can override and force far /transcription */

                prod = BioseqLockById (sip);
                unlockme = prod;
              }
            }
          }
          if (rrp->type == 2) {
            if (! pseudo) {
              if (ajp->showTranscript) {
                qvp [FTQUAL_transcription].ble = TRUE;
              }
            }
          }
          if (rrp->type == 3) {
            if (rrp->ext.choice == 1) {
              /* amino acid could not be parsed into structured form */
              if (! ajp->flags.dropIllegalQuals) {
                str = (CharPtr) rrp->ext.value.ptrvalue;
                qvp [FTQUAL_product].str = str;
              } else {
                qvp [FTQUAL_product].str = "tRNA-OTHER";
              }
            } else if (rrp->ext.choice == 2) {
              trna = (tRNAPtr) rrp->ext.value.ptrvalue;
              if (trna != NULL) {
                aa = 0;
                if (trna->aatype == 2) {
                  aa = trna->aa;
                } else {
                  from = 0;
                  switch (trna->aatype) {
                    case 0 :
                      from = 0;
                      break;
                    case 1 :
                      from = Seq_code_iupacaa;
                      break;
                    case 2 :
                      from = Seq_code_ncbieaa;
                      break;
                    case 3 :
                      from = Seq_code_ncbi8aa;
                      break;
                    case 4 :
                      from = Seq_code_ncbistdaa;
                      break;
                    default:
                      break;
                  }
                  if (ajp->flags.iupacaaOnly) {
                    code = Seq_code_iupacaa;
                  } else {
                    code = Seq_code_ncbieaa;
                  }
                  smtp = SeqMapTableFind (code, from);
                  if (smtp != NULL) {
                    aa = SeqMapTableConvert (smtp, trna->aa);
                    if (aa == 255 && from == Seq_code_iupacaa) {
                      if (trna->aa == 'U') {
                        aa = 'U';
                      } else if (trna->aa == 'O') {
                        aa = 'O';
                      }
                    }
                  }
                }
                if (ajp->flags.iupacaaOnly) {
                  smtp = SeqMapTableFind (Seq_code_iupacaa, Seq_code_ncbieaa);
                  if (smtp != NULL) {
                    aa = SeqMapTableConvert (smtp, trna->aa);
                  } else {
                    aa = 'X';
                  }
                }
                if (aa > 0 && aa != 255) {
                  /*
                  if (aa == 'U') {
                     if (ajp->flags.selenocysteineToNote) {
                      qvp [FTQUAL_selenocysteine_note].str = "selenocysteine";
                    } else {
                      qvp [FTQUAL_selenocysteine].ble = TRUE;
                    }
                  } else if (aa == 'O') {
                    if (ajp->flags.pyrrolysineToNote) {
                      qvp [FTQUAL_pyrrolysine_note].str = "pyrrolysine";
                    } else {
                      qvp [FTQUAL_pyrrolysine].ble = TRUE;
                    }
                  }
                  */
                  if (ajp->mode == RELEASE_MODE || ajp->mode == ENTREZ_MODE) {
                    /* O and J no longer quarantined */
                    /*
                    if (aa == 79 || aa == 74) {
                      aa = 88;
                    }
                    */
                  }
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
                    idx = 25;
                  }
                  if (idx > 0 && idx < 28) {
                    str = trnaList [idx];
                    qvp [FTQUAL_product].str = str;
                    if (StringNICmp (str, "tRNA-", 5) == 0) {
                      qvp [FTQUAL_trna_aa].str = str + 5;
                    }
                  }
                }
                qvp [FTQUAL_anticodon].slp = trna->anticodon;
                if (ajp->flags.codonRecognizedToNote) {
                  qvp [FTQUAL_trna_codons_note].trp = trna;
                } else {
                  qvp [FTQUAL_trna_codons].trp = trna;
                }
              }
            }
          } else {
            if (rrp->ext.choice == 1) {
              if (rrp->type == 255) {
                str = (CharPtr) rrp->ext.value.ptrvalue;
                if (StringCmp (str, "ncRNA") == 0 ||
                    StringCmp (str, "tmRNA") == 0) {
                   /* pick up product from gbqual */
                 } else if (StringICmp (str, "misc_RNA") == 0) {
                   is_misc_rna = TRUE;
                   /* pick up product from gbqual */
                 } else {
                  str = (CharPtr) rrp->ext.value.ptrvalue;
                  qvp [FTQUAL_product].str = str;
                 }
              } else {
                str = (CharPtr) rrp->ext.value.ptrvalue;
                qvp [FTQUAL_product].str = str;
              }
            }
          }
          if (rrp->type == 10) {
            is_misc_rna = TRUE;
          }
          if (rrp->ext.choice == 3) {
            rgp = (RNAGenPtr) rrp->ext.value.ptrvalue;
            if (rgp != NULL) {
              if (StringDoesHaveText (rgp->product)) {
                qvp [FTQUAL_product].str = rgp->product;
              }
              if (StringDoesHaveText (rgp->_class)) {
                if (IsStringInNcRNAClassList (rgp->_class)) {
                  qvp [FTQUAL_ncRNA_other].str = rgp->_class;
                } else {
                  qvp [FTQUAL_ncRNA_other].str = "other";
                  qvp [FTQUAL_ncRNA_note].str = rgp->_class;
                }
              }
              for (rqsp = rgp->quals; rqsp != NULL; rqsp = rqsp->next) {
                if (StringICmp (rqsp->qual, "tag_peptide") == 0) {
                  if (StringDoesHaveText (rqsp->val)) {
                    qvp [FTQUAL_tag_peptide_str].str = rqsp->val;
                  }
                }
              }
            }
          }
        }
        if (rrp != NULL && rrp->ext.choice == 3) {
          rgp = (RNAGenPtr) rrp->ext.value.ptrvalue;
          if (rgp != NULL) {
            if (StringDoesHaveText (rgp->product)) {
              str = rgp->product;
              if (StringNICmp (str, "internal transcribed spacer ", 28) == 0) {
                str += 28;
                if (IS_DIGIT (*str) && str [1] == '\0') {
                  its_prod = str;
                }
              }
            }
          }
        }
        for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
          if (StringCmp (gbq->qual, "product") != 0) continue;
          if (StringDoesHaveText (gbq->val)) {
            str = gbq->val;
            if (StringNICmp (str, "internal transcribed spacer ", 28) == 0) {
              str += 28;
              if (IS_DIGIT (*str) && str [1] == '\0') {
                its_prod = str;
              }
            }
          }
        }
        if (is_misc_rna && StringDoesHaveText (its_prod)) {
          if (StringCmp (its_prod, "1") == 0) {
            qvp [FTQUAL_rrna_its].str = "ITS1";
          } else if (StringCmp (its_prod, "2") == 0) {
            qvp [FTQUAL_rrna_its].str = "ITS2";
          } else if (StringCmp (its_prod, "3") == 0) {
            qvp [FTQUAL_rrna_its].str = "ITS3";
          }
        }
        break;
      case SEQFEAT_REGION :
        if (format == GENPEPT_FMT && featdeftype == FEATDEF_REGION && isProt) {
          qvp [FTQUAL_region_name].str = (CharPtr) sfp->data.value.ptrvalue;
        } else {
          qvp [FTQUAL_region].str = (CharPtr) sfp->data.value.ptrvalue;
        }
        if (sfp->ext != NULL) {
          uop = FindUopByTag (sfp->ext, "cddScoreData");
          if (uop != NULL) {
            for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
              if (ufp->choice != 1) continue;
              oip = ufp->label;
              if (oip == NULL) continue;
              if (StringICmp (oip->str, "definition") == 0) {
                str = (CharPtr) ufp->data.ptrvalue;
                if (StringDoesHaveText (str)) {
                  if (StringICmp (str, (CharPtr) sfp->data.value.ptrvalue) != 0) {
                    qvp [FTQUAL_cdd_definition].str = str;
                  }
                }
              }
            }
          }
        }
        break;
      case SEQFEAT_COMMENT :
        break;
      case SEQFEAT_BOND :
        bondidx = (Int2) sfp->data.value.intvalue;
        if (bondidx == 255) {
          bondidx = 5;
        }
        if (bondidx > 0 && bondidx < 6) {
          if (format == GENPEPT_FMT && isProt) {
            qvp [FTQUAL_bond_type].str = bondList [bondidx];
          } else {
            qvp [FTQUAL_bond].str = bondList [bondidx];
          }
        }
        break;
      case SEQFEAT_SITE :
        siteidx = (Int2) sfp->data.value.intvalue;
        if (siteidx == 255) {
          siteidx = 27;
        }
        if (siteidx > 0 && siteidx < 28) {
          if (format == GENPEPT_FMT && isProt) {
            qvp [FTQUAL_site_type].str = siteFFList [siteidx];
          } else {
            qvp [FTQUAL_site].str = siteFFList [siteidx];
          }
        }
        break;
      case SEQFEAT_PSEC_STR :
        qvp [FTQUAL_sec_str_type].num = sfp->data.value.intvalue;
        break;
      case SEQFEAT_HET :
        qvp [FTQUAL_heterogen].str = (CharPtr) sfp->data.value.ptrvalue;
        break;
      case SEQFEAT_VARIATIONREF :
        vrp = (VariationRefPtr) sfp->data.value.ptrvalue;
        if (vrp != NULL) {
          qvp [FTQUAL_variation_id].dbt = vrp->id;
          vnp = vrp->data;
          if (vnp != NULL) {
            if (vnp->choice == VarRefData_instance) {
              vip = (VariationInstPtr) vnp->data.ptrvalue;
              if (vip != NULL) {
                qvp [FTQUAL_delta_item].vnp = vip->delta;
              }
            } else if (vnp->choice == VarRefData_set) {
              vsp = (VarRefDataSetPtr) vnp->data.ptrvalue;
              if (vsp != NULL) {
                qvp [FTQUAL_variation_set].vnp = vsp->variations;
              }
            }
          }
        }
        break;
      default :
        break;
    }
  }

  /* common fields set here */

  if (ajp->mode == DUMP_MODE && qvp [FTQUAL_gene_syn_refseq].vnp != NULL) {
    qvp [FTQUAL_gene_syn].vnp = qvp [FTQUAL_gene_syn_refseq].vnp;
    qvp [FTQUAL_gene_syn_refseq].vnp = NULL;
  }

  VisitUserObjectsInUop (sfp->ext, (Pointer) qvp, RecordUserObjectsInQVP);

  if (ajp->hideGoTerms) {
    qvp [FTQUAL_go_process].ufp = NULL;
    qvp [FTQUAL_go_component].ufp = NULL;
    qvp [FTQUAL_go_function].ufp = NULL;
  }

  if (featdeftype == FEATDEF_repeat_region || featdeftype == FEATDEF_mobile_element) {
    pseudo = FALSE;
  }

  qvp [FTQUAL_pseudo].ble = pseudo;

  qvp [FTQUAL_seqfeat_note].str = sfp->comment;

  sap = fcontext->sap;
  if (sap != NULL) {
    annotDescCommentToComment = FALSE;
    for (adp = sap->desc; adp != NULL; adp = adp->next) {
      if (adp->choice == Annot_descr_comment) {
        if (StringDoesHaveText ((CharPtr) adp->data.ptrvalue)) {
          qvp [FTQUAL_seqannot_note].str = (CharPtr) adp->data.ptrvalue;
        }
      } else if (adp->choice == Annot_descr_user) {
        uop = (UserObjectPtr) adp->data.ptrvalue;
        if (uop == NULL) continue;
        oip = uop->type;
        if (oip == NULL) continue;
        if (StringCmp (oip->str, "AnnotDescCommentPolicy") == 0) {
          for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
            oip = ufp->label;
            if (oip == NULL || ufp->data.ptrvalue == NULL) continue;
            if (StringCmp (oip->str, "Policy") == 0) {
              if (StringICmp ((CharPtr) ufp->data.ptrvalue, "ShowInComment") == 0) {
                annotDescCommentToComment = TRUE;
              }
            }
          }
        }
      }
    }
    if (annotDescCommentToComment) {
      qvp [FTQUAL_seqannot_note].str = NULL;
    }
  }

  /* if RELEASE_MODE, check list of features that can have /pseudo */

  if (ajp->flags.dropIllegalQuals && pseudo  &&
      (seqfeattype == SEQFEAT_RNA || seqfeattype == SEQFEAT_IMP) ) {
    switch (featdeftype) {

    case  FEATDEF_allele:
    case  FEATDEF_attenuator:
    case  FEATDEF_CAAT_signal:
    case  FEATDEF_conflict:
    case  FEATDEF_D_loop:
    case  FEATDEF_enhancer:
    case  FEATDEF_GC_signal:
    case  FEATDEF_iDNA:
    case  FEATDEF_LTR:
    case  FEATDEF_misc_binding:
    case  FEATDEF_misc_difference:
    case  FEATDEF_misc_recomb:
    case  FEATDEF_misc_signal:
    case  FEATDEF_misc_structure:
    case  FEATDEF_modified_base:
    case  FEATDEF_mobile_element:
    case  FEATDEF_mutation:
    case  FEATDEF_old_sequence:
    case  FEATDEF_polyA_signal:
    case  FEATDEF_polyA_site:
    case  FEATDEF_precursor_RNA:
    case  FEATDEF_prim_transcript:
    case  FEATDEF_primer_bind:
    case  FEATDEF_protein_bind:
    case  FEATDEF_RBS:
    case  FEATDEF_repeat_region:
    case  FEATDEF_repeat_unit:
    case  FEATDEF_rep_origin:
    case  FEATDEF_satellite:
    case  FEATDEF_stem_loop:
    case  FEATDEF_STS:
    case  FEATDEF_TATA_signal:
    case  FEATDEF_terminator:
    case  FEATDEF_unsure:
    case  FEATDEF_variation:
    case  FEATDEF_3clip:
    case  FEATDEF_3UTR:
    case  FEATDEF_5clip:
    case  FEATDEF_5UTR:
    case  FEATDEF_10_signal:
    case  FEATDEF_35_signal:
      qvp [FTQUAL_pseudo].ble = FALSE;
        break;
    default:
        break;
    }
  }

  /*
  if (format != GENPEPT_FMT) {
    qvp [FTQUAL_evidence].num = sfp->exp_ev;
  }
  */
  qvp [FTQUAL_evidence].num = sfp->exp_ev;

  if (sfp->excpt && StringDoesHaveText (sfp->except_text)) {
    ParseException (sfp->except_text,
                    &exception_string,
                    &exception_note,
                    is_other,
                    (Boolean) (! ajp->flags.dropIllegalQuals),
                    sfp->idx.subtype,
                    &riboSlippage,
                    &transSplice,
                    &artLoc,
                    &hetPop,
                    &lowQual);

    qvp [FTQUAL_exception].str = exception_string;
    qvp [FTQUAL_exception_note].str = exception_note;
    qvp [FTQUAL_ribosomal_slippage].ble = riboSlippage;
    qvp [FTQUAL_trans_splicing].ble = transSplice;
    qvp [FTQUAL_artificial_location].ble = artLoc;
    if (hetPop) {
      qvp [FTQUAL_artificial_location_str].str = "heterogeneous population sequenced";
    }
    if (lowQual) {
      qvp [FTQUAL_artificial_location_str].str = "low-quality sequence region";
    }

    /*
    if (StringHasNoText (qvp [FTQUAL_exception].str)) {
      qvp [FTQUAL_exception].str = NULL;
    }
    if (StringHasNoText (qvp [FTQUAL_exception_note].str)) {
      qvp [FTQUAL_exception_note].str = NULL;
    }
    */
  }

  qvp [FTQUAL_db_xref].vnp = sfp->dbxref;
  qvp [FTQUAL_citation].vnp = sfp->cit;

  /* /product same as sfp->comment will suppress /note */

  if (! StringHasNoText (qvp [FTQUAL_product].str) &&
      StringICmp (sfp->comment, qvp [FTQUAL_product].str) == 0) {
    qvp [FTQUAL_seqfeat_note].str = NULL;
  }
  /* case sensitive AJ011317.1 */
  if (! StringHasNoText (qvp [FTQUAL_cds_product].str) &&
      StringCmp (sfp->comment, qvp [FTQUAL_cds_product].str) == 0) {
    qvp [FTQUAL_seqfeat_note].str = NULL;
  }

  /* /gene same as sfp->comment will suppress /note */
  /* case sensitive -gi|6572973|gb|AF195052.1|AF195052 */

  if (! StringHasNoText (qvp [FTQUAL_gene].str) &&
      StringCmp (sfp->comment, qvp [FTQUAL_gene].str) == 0) {
    qvp [FTQUAL_seqfeat_note].str = NULL;
  }

  /* gene /note same as sfp->comment will suppress /note - U92435.1 says do not do this */

  /*
  if (! StringHasNoText (qvp [FTQUAL_gene_note].str) &&
      StringICmp (sfp->comment, qvp [FTQUAL_gene_note].str) == 0) {
    qvp [FTQUAL_seqfeat_note].str = NULL;
  }
  */

  /* if site sfp->comment contains site name, suppress site in /note */

  if (! StringHasNoText (qvp [FTQUAL_site].str) &&
      StringStr (sfp->comment, qvp [FTQUAL_site].str) != NULL) {
    qvp [FTQUAL_site].str = NULL;
  }

  /* /EC_number same as sfp->comment will suppress /note */

  for (vnp = qvp [FTQUAL_prot_EC_number].vnp; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if ((! StringHasNoText (str)) &&
        StringICmp (sfp->comment, str) == 0) {
      qvp [FTQUAL_seqfeat_note].str = NULL;
    }
  }

  /* mat_peptide precursor same as sfp->comment will suppress /note in GenPept */

  if (precursor_comment != NULL && StringCmp (precursor_comment, sfp->comment) == 0) {
    qvp [FTQUAL_seqfeat_note].str = NULL;
  }


  /* now go through gbqual list */

  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    idx = GbqualToFeaturIndex (gbq->qual);
    if (idx > 0 && idx < ASN2GNBK_TOTAL_FEATUR) {
      if (qvp [idx].gbq == NULL) {
        if (idx == FTQUAL_product_quals) {
          if (qvp [FTQUAL_product].str == NULL) {
            qvp [FTQUAL_product].str = gbq->val;
          } else if (qvp [FTQUAL_xtra_prod_quals].gbq == NULL) {
            /* chain will include remaining product gbquals */
            qvp [FTQUAL_xtra_prod_quals].gbq = gbq;
          }
        } else {
          qvp [idx].gbq = gbq;
        }
      }

    } else if (idx == 0) {

      qualclass = IllegalGbqualToClass (gbq->qual);
      if (qualclass == 0) {
        qualclass = Qual_class_quote;
      }
      tmp = StringSave (gbq->val);
      if (tmp != NULL) {
        str = MemNew (sizeof (Char) * (StringLen (gbq->qual) + StringLen (tmp) + 10));
        if (str != NULL) {
          if (qualclass == Qual_class_quote) {
            if (StringIsJustQuotes (tmp)) {
              sprintf (str, "/%s", gbq->qual);
            } else {
              ptr = tmp;
              ch = *ptr;
              while (ch != '\0') {
                if (ch == '"') {
                  *ptr = '\'';
                }
                ptr++;
                ch = *ptr;
              }
              sprintf (str, "/%s=\"%s\"", gbq->qual, tmp);
            }
            ValNodeCopyStr (&illegal, 0, str);
          } else if (qualclass == Qual_class_noquote || qualclass == Qual_class_label) {
            if (StringIsJustQuotes (tmp)) {
              sprintf (str, "/%s", gbq->qual);
            } else {
              sprintf (str, "/%s=%s", gbq->qual, tmp);
            }
            ValNodeCopyStr (&illegal, 0, str);
          }
          MemFree (str);
        }
        MemFree (tmp);
      }
    }
  }

  /* experiment and inference qualifiers supercede evidence qualifier */

  if (qvp [FTQUAL_evidence].num > 0) {
    if (qvp [FTQUAL_experiment].gbq != NULL || qvp [FTQUAL_inference].gbq != NULL) {
      qvp [FTQUAL_evidence].num = 0;
    } else if (qvp [FTQUAL_evidence].num == 1) {
      qvp [FTQUAL_experiment_string].str = "experimental evidence, no additional details recorded";
      qvp [FTQUAL_evidence].num = 0;
    } else if (qvp [FTQUAL_evidence].num == 2) {
      qvp [FTQUAL_inference_string].str = "non-experimental evidence, no additional details recorded";
      qvp [FTQUAL_evidence].num = 0;
    }
  }

  if (qvp [FTQUAL_ncRNA_class].gbq != NULL) {
    gbq = qvp [FTQUAL_ncRNA_class].gbq;
    if (StringDoesHaveText (gbq->val)) {
      if (! IsStringInNcRNAClassList (gbq->val)) {
        qvp [FTQUAL_ncRNA_other].str = "other";
        qvp [FTQUAL_ncRNA_note].str = gbq->val;
        qvp [FTQUAL_ncRNA_class].gbq = NULL;
      }
    }
  }

  if (ajp->mode != DUMP_MODE) {
    ParseInference (qvp [FTQUAL_inference].gbq, &good_inference, &bad_inference);
    qvp [FTQUAL_inference_good].vnp = good_inference;
    qvp [FTQUAL_inference_bad].vnp = bad_inference;
    qvp [FTQUAL_inference].gbq = NULL;
  }

  /* optionally suppress evidence, inference and experiment qualifiers */

  if (ajp->hideEvidence) {
    qvp [FTQUAL_inference_good].vnp = NULL;
    qvp [FTQUAL_inference_bad].vnp = NULL;
    qvp [FTQUAL_inference].gbq = NULL;
    qvp [FTQUAL_experiment].gbq = NULL;
    qvp [FTQUAL_evidence].num = 0;
    qvp [FTQUAL_experiment_string].gbq = NULL;
    qvp [FTQUAL_inference_string].gbq = NULL;
  }

  /* special handling for cyt_map, gen_map, rad_map */

  if (ajp->flags.hideSpecificGeneMaps) {
    if (qvp [FTQUAL_gene_map].str == NULL) {
      gbq = qvp [FTQUAL_gene_cyt_map].gbq;
      if (gbq != NULL) {
        qvp [FTQUAL_gene_map].str = gbq->val;
      } else {
        gbq = qvp [FTQUAL_gene_gen_map].gbq;
        if (gbq != NULL) {
          qvp [FTQUAL_gene_map].str = gbq->val;
        } else {
          gbq = qvp [FTQUAL_gene_rad_map].gbq;
          if (gbq != NULL) {
            qvp [FTQUAL_gene_map].str = gbq->val;
          }
        }
      }
    }
    qvp [FTQUAL_gene_cyt_map].gbq = NULL;
    qvp [FTQUAL_gene_gen_map].gbq = NULL;
    qvp [FTQUAL_gene_rad_map].gbq = NULL;
  }

  /* illegal qualifiers are copied and formatted in valnode chain */

  if (! ajp->flags.dropIllegalQuals) {
    qvp [FTQUAL_illegal_qual].vnp = illegal;
  }

  /* remove protein description that equals the gene name, case sensitive */

  if (StringCmp (qvp [FTQUAL_gene].str, qvp [FTQUAL_prot_desc].str) == 0) {
    qvp [FTQUAL_prot_desc].str = NULL;
  }

  /* remove protein description that equals the cds product, case sensitive */

  if (StringCmp (qvp [FTQUAL_cds_product].str, qvp [FTQUAL_prot_desc].str) == 0) {
    qvp [FTQUAL_prot_desc].str = NULL;
  }

  /* remove comment contained in prot_desc - gi|4530123|gb|AF071539.1|AF071539 */

  if (StringStr (qvp [FTQUAL_prot_desc].str, qvp [FTQUAL_seqfeat_note].str) != NULL) {
    qvp [FTQUAL_seqfeat_note].str = NULL;
  }

  /* remove protein description that equals the standard name */

  if (qvp [FTQUAL_standard_name].gbq != NULL && qvp [FTQUAL_prot_desc].str != NULL) {
    gbq = qvp [FTQUAL_standard_name].gbq;
    lasttype = gbq->qual;
    while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
      if (StringICmp (gbq->val, qvp [FTQUAL_prot_desc].str) == 0) {
        qvp [FTQUAL_prot_desc].str = NULL;
      }
      gbq = gbq->next;
    }
  }

  /* remove protein description that equals a gene synonym - case insensitive AF109216.1 */

  for (vnp = gene_syn; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if ((! StringHasNoText (str)) &&
        StringCmp (str, qvp [FTQUAL_prot_desc].str) == 0) {
      /* NC_001823 leave in prot_desc if no cds_product */
      if (qvp [FTQUAL_cds_product].str != NULL) {
        qvp [FTQUAL_prot_desc].str = NULL;
      }
    }
  }

  /* remove comment that equals a gene synonym */

  if (format != GENPEPT_FMT && (! ifp->mapToProt)) {
    for (vnp = gene_syn; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if ((! StringHasNoText (str)) &&
          StringICmp (str, qvp [FTQUAL_seqfeat_note].str) == 0) {
        qvp [FTQUAL_seqfeat_note].str = NULL;
      }
    }
  }

  /* remove protein comment descriptor that equals the protein note */

  if (StringCmp (qvp [FTQUAL_prot_note].str, qvp [FTQUAL_prot_comment].str) == 0) {
    qvp [FTQUAL_prot_comment].str = NULL;
  }

  /* suppress cds comment if a subset of protein note - AF002218.1 */

  if (StringStr (qvp [FTQUAL_prot_note].str, qvp [FTQUAL_seqfeat_note].str) != NULL) {
    qvp [FTQUAL_seqfeat_note].str = NULL;
  }

  /* suppress selenocysteine note if already in comment */

  if (StringStr (sfp->comment, "selenocysteine") != NULL) {
    qvp [FTQUAL_selenocysteine_note].str = NULL;
  }

  /* suppress pyrrolysine note if already in comment */

  if (StringStr (sfp->comment, "pyrrolysine") != NULL) {
    qvp [FTQUAL_pyrrolysine_note].str = NULL;
  }

  /* if /allele inherited from gene, suppress allele gbqual on feature */

  if (qvp [FTQUAL_gene_allele].str != NULL) {
    qvp [FTQUAL_allele].gbq = NULL;
  }

  /* now print qualifiers from table */

#ifdef DISPLAY_STRINGS
  s_DisplayQVP(qvp, feat_note_order);
#endif

  /* Strip duplicate notes */

  if ((StringCmp(qvp[FTQUAL_product].str,
         qvp[FTQUAL_seqfeat_note].str) == 0)) {
    qvp[FTQUAL_seqfeat_note].str = NULL;
  }

  if ((qvp[FTQUAL_standard_name].gbq != NULL) &&
      (qvp[FTQUAL_standard_name].gbq->val != NULL)) {
    if ((StringCmp(qvp[FTQUAL_seqfeat_note].str,
           qvp[FTQUAL_standard_name].gbq->val) == 0)) {
      qvp[FTQUAL_seqfeat_note].str = NULL;
    }
  }

  /* Display strings for debugging purposes */

#ifdef DISPLAY_STRINGS
  s_DisplayQVP(qvp, feat_qual_order);
#endif

  /*
  qvp[FTQUAL_loc_debug_str1].str = firstloc;
  qvp[FTQUAL_loc_debug_str2].str = secondloc;
  qvp[FTQUAL_loc_debug_str3].str = thirdloc;
  */

  /* optionally populate indexes for NCBI internal database */

  if (index != NULL) {
    if (! StringHasNoText (qvp [FTQUAL_gene].str)) {
      ValNodeCopyStrToHead (&(index->genes), 0, qvp [FTQUAL_gene].str);
    }
  }

  if (doKey) {
    FFAddOneChar(ffstring, '\n', FALSE);
  }

  /* Build the flat file */
  FormatFeatureBlockQuals (ffstring, ajp, asp, bsp, featdeftype, gene_syn,
                           lasttype, location, prod,
                           protein_pid_g, qvp,
                           left, right, strand,
                           sfp, target, ifp, is_other,
                           is_journalscan, is_gps, is_ged);

  /* ??? and then deal with the various note types separately (not in order table) ??? */

  /* free aa-dna or dna-aa mapped location */

  SeqLocFree (loc);

  /*
  MemFree (firstloc);
  MemFree (secondloc);
  MemFree (thirdloc);
  */

  ValNodeFreeData (illegal);
  MemFree (exception_string);
  MemFree (exception_note);
  ValNodeFreeData (good_inference);
  ValNodeFreeData (bad_inference);

  BioseqUnlock (unlockme);

  if (afp != NULL && GetWWW (ajp) && ajp->mode == ENTREZ_MODE && ajp->seqspans &&
    (ajp->format == GENBANK_FMT || ajp->format == GENPEPT_FMT)) {
    sprintf (sfx, "</span>");
  }

  str = NULL;

  if (js != NULL) {
    str = FFEndPrintEx (ajp, ffstring, format, 21, 21, 21, 21, "FT", js, sfx);
  } else {
    str = FFEndPrintEx (ajp, ffstring, format, 21, 21, 21, 21, "FT", pfx, sfx);
  }

  MemFree (js);

  /* optionally populate gbseq for XML-ized GenBank format */

  if (gbseq != NULL) {
    if (gbfeat != NULL) {
      AddFeatureToGbseq (gbseq, gbfeat, str, sfp);
    }
  }

  FFRecycleString(ajp, ffstring);
  return str;
}

NLM_EXTERN CharPtr FormatFeatureBlock (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp
)

{
  IntAsn2gbJobPtr    ajp;
  Asn2gbSectPtr      asp;
  BioseqPtr          bsp;
  SeqMgrFeatContext  fcontext;
  FmtType            format;
  ValNodePtr         head;
  QualValPtr         qvp;
  SeqFeatPtr         sfp;
  CharPtr            str;
  BioseqPtr          target;

  if (afp == NULL || bbp == NULL) return NULL;
  asp = afp->asp;
  if (asp == NULL) return NULL;
  target = asp->target;
  bsp = asp->bsp;
  if (target == NULL || bsp == NULL) return NULL;
  ajp = afp->ajp;
  if (ajp == NULL) return NULL;
  qvp = afp->qvp;
  if (qvp == NULL) return NULL;
  format = afp->format;

  /* all features in this list are known to be valid for the designated mode */

  sfp = SeqMgrGetDesiredFeature (bbp->entityID, NULL, bbp->itemID, 0, NULL, &fcontext);
  if (sfp == NULL) return NULL;

  /* five-column feature table uses special code for formatting */

  if (ajp->format == FTABLE_FMT) {
    head = NULL;
    PrintFtableLocAndQuals (ajp, &head, target, sfp, &fcontext);
    str = MergeFFValNodeStrs (head);
    ValNodeFreeData (head);
    return str;
  }

  /* otherwise do regular flatfile formatting */

  return FormatFeatureBlockEx (afp, ajp, asp, bsp, target, sfp, &fcontext, qvp,
                               format, (IntFeatBlockPtr) bbp, ISA_aa (bsp->mol), TRUE);
}

const CharPtr feature_table_header_format = ">Feature %s\n";

NLM_EXTERN CharPtr FormatFeatHeaderBlock (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp
)

{
  IntAsn2gbJobPtr  ajp;
  Asn2gbSectPtr    asp;
  BioseqPtr        bsp;
  Char             ch;
  Boolean          has_space;
  Char             id [128];
  ObjectIdPtr      oip;
  CharPtr          ptr;
  SeqIdPtr         sip;
  SeqIdPtr         sip2;
  CharPtr          str = NULL;
  BioseqPtr        target;
  CharPtr          tmp = NULL;

  if (afp == NULL || bbp == NULL) return NULL;
  ajp = afp->ajp;
  if (ajp == NULL) return NULL;
  asp = afp->asp;
  if (asp == NULL) return NULL;
  target = asp->target;
  bsp = asp->bsp;
  if (target == NULL || bsp == NULL) return NULL;

  /* five-column feature table uses special code for formatting */

  if (ajp->format == FTABLE_FMT) {
    sip = SeqIdFindBest (target->id, 0);
    if (sip == NULL) return NULL;
    id [0] = '\0';

    if (sip->choice == SEQID_GI) {
      sip2 = GetSeqIdForGI (sip->data.intvalue);
      if (sip2 != NULL) {
        sip = sip2;
      }
    }
    if (sip->choice == SEQID_LOCAL) {
      oip = (ObjectIdPtr) sip->data.ptrvalue;
      if (oip != NULL && StringDoesHaveText (oip->str)) {
        has_space = FALSE;
        ptr = oip->str;
        ch = *ptr;
        while (ch != '\0') {
          if (IS_WHITESP (ch)) {
            has_space = TRUE;
          }
          ptr++;
          ch = *ptr;
        }
        if (has_space) {
          sprintf (id, "lcl|%c%s%c", (char) '"', oip->str, (char) '"');
        }
      }
    }

    if (id [0] == '\0') {
      SeqIdWrite (sip, id, PRINTID_FASTA_LONG, sizeof (id) - 1);
    }
    if (! StringHasNoText (id)) {
      tmp = (CharPtr) MemNew ((StringLen(feature_table_header_format) + StringLen(id)) * sizeof(Char));
      sprintf (tmp, ">Feature %s\n", id);
      str = tmp;
    }
    return str;
  }

  /* otherwise do regular flatfile formatting */

  return StringSaveNoNull (bbp->string);
}


/* stand alone function to produce qualifiers in genbank style */

static void StripLeadingSpaces (
  CharPtr str
)

{
  Uchar    ch;
  CharPtr  dst;
  CharPtr  ptr;


  if (str == NULL || str [0] == '\0') return;

  dst = str;
  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    while (ch == ' ') {
      ptr++;
      ch = *ptr;
    }
    while (ch != '\n' && ch != '\r') {
      *dst = ch;
      dst++;
      ptr++;
      ch = *ptr;
    }
    *dst = ch;
    dst++;
    ptr++;
    ch = *ptr;
  }
  *dst = '\0';
}

NLM_EXTERN void DoImmediateRemoteFeatureFormat (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp,
  SeqFeatPtr sfp
)

{
  IntAsn2gbJobPtr    ajp;
  Asn2gbSectPtr      asp;
  BlockType          blocktype;
  BioseqPtr          bsp;
  SeqMgrFeatContext  fcontext;
  size_t             max;
  SeqEntryPtr        oldscope;
  QualValPtr         qvp = NULL;
  SeqEntryPtr        sep;
  SeqLocPtr          slp;
  CharPtr            str = NULL;
  BioseqPtr          target;

  if (afp == NULL || bbp == NULL || sfp == NULL) return;
  asp = afp->asp;
  if (asp == NULL) return;
  target = asp->target;
  bsp = asp->bsp;
  if (target == NULL || bsp == NULL) return;
  ajp = afp->ajp;
  if (ajp == NULL) return;

  blocktype = bbp->blocktype;
  if (blocktype < LOCUS_BLOCK || blocktype > SLASH_BLOCK) return;

  max = (size_t) (MAX (ASN2GNBK_TOTAL_SOURCE, ASN2GNBK_TOTAL_FEATUR));
  qvp = MemNew (sizeof (QualVal) * (max + 5));
  if (qvp == NULL) return;

  MemSet ((Pointer) &fcontext, 0, sizeof (SeqMgrFeatContext));
  fcontext.itemID = 0;
  fcontext.featdeftype = sfp->idx.subtype;
  fcontext.seqfeattype = sfp->data.choice;
  slp = sfp->location;
  fcontext.left = GetOffsetInBioseq (slp, bsp, SEQLOC_LEFT_END);
  fcontext.right = GetOffsetInBioseq (slp, bsp, SEQLOC_RIGHT_END);
  fcontext.strand = SeqLocStrand (slp);

  sep = GetTopSeqEntryForEntityID (bbp->entityID);

  oldscope = SeqEntrySetScope (sep);

  if (ajp->format != FTABLE_FMT) {
    str = FormatFeatureBlockEx (afp, ajp, asp, bsp, target, sfp, &fcontext, qvp,
                                ajp->format, (IntFeatBlockPtr) bbp, ISA_aa (bsp->mol), TRUE);
  }

  SeqEntrySetScope (oldscope);

  if (str != NULL) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "%s", str);
    }
    if (afp->ffwrite != NULL) {
      afp->ffwrite (str, afp->userdata, blocktype, sfp->idx.entityID, OBJ_SEQFEAT,
                    sfp->idx.itemID, fcontext.left, fcontext.right);
    }
  } else {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "?\n");
    }
    if (afp->ffwrite != NULL) {
      afp->ffwrite ("?\n", afp->userdata, blocktype, sfp->idx.entityID, OBJ_SEQFEAT,
                    sfp->idx.itemID, fcontext.left, fcontext.right);
    }
  }

  MemFree (str);
  MemFree (qvp
  );
}

NLM_EXTERN CharPtr FormatFeatureQuals (
  SeqFeatPtr sfp
)

{
  IntAsn2gbJob       ajb;
  IntAsn2gbJobPtr    ajp;
  BioseqPtr          bsp;
  SeqMgrFeatContext  fcontext;
  IntCdsBlock        ifb;
  IntFeatBlockPtr    ifp;
  size_t             max;
  QualValPtr         qvp;
  CharPtr            str;

  if (sfp == NULL) return NULL;
  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) return NULL;

  if (SeqMgrGetDesiredFeature (0, bsp, 0, 0, sfp, &fcontext) != sfp) return NULL;

  MemSet ((Pointer) &ajb, 0, sizeof (IntAsn2gbJob));
  ajp = &ajb;
  MemSet ((Pointer) &ifb, 0, sizeof (IntCdsBlock));
  ifp = (IntFeatBlockPtr) &ifb;

  max = (size_t) (MAX (ASN2GNBK_TOTAL_SOURCE, ASN2GNBK_TOTAL_FEATUR));
  qvp = MemNew (sizeof (QualVal) * (max + 5));
  if (qvp == NULL) return NULL;

  str = FormatFeatureBlockEx (NULL, ajp, NULL, NULL, NULL, sfp, &fcontext, qvp,
                              GENBANK_FMT, ifp, FALSE, FALSE);

  MemFree (qvp);
  StripLeadingSpaces (str);
  return str;
}

