/*   diffshift.c
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
* File Name:  diffshift.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   11/26/10
*
* $Revision: 1.1 $
*
* File Description:
*
*  Processes diff output, shifting spaces to minimize offsets per group
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
*
* ==========================================================================
*/

#include <ncbi.h>

typedef struct diffblock {
  ValNodePtr  head;
  ValNodePtr  tail;
} DiffBlock, PNTR DiffBlockPtr;

static void RecordDiffBlock (
  DiffBlockPtr dbp,
  CharPtr str
)

{
  ValNodePtr  vnp;

  if (dbp == NULL || StringHasNoText (str)) return;

  vnp = ValNodeCopyStr (&(dbp->tail), 0, str);
  if (dbp->head == NULL) {
    dbp->head = vnp;
  }
  dbp->tail = vnp;
}

static void WriteDiffBlock (
  DiffBlockPtr dbp,
  FILE *ofp
)

{
  Char        ch;
  Int2        idx;
  Int2        margin = INT2_MAX;
  CharPtr     ptr;
  Int2        spaces;
  CharPtr     str;
  ValNodePtr  vnp;

  if (dbp == NULL || dbp->head == NULL || ofp == NULL) return;

  for (vnp = dbp->head; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    ch = str [0];
    if (ch == '<' || ch == '>') {
      ptr = str + 1;
      ch = *ptr;
      spaces = 0;
      while (ch == ' ') {
        spaces++;
        ptr++;
        ch = *ptr;
      }
      if (spaces < margin) {
        margin = spaces;
      }
    }
  }

  if (margin > 80) {
    margin = 80;
  }

  for (vnp = dbp->head; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    ch = str [0];
    if (ch == '<' || ch == '>') {
      ptr = str + 1;
      ch = *ptr;
      idx = 0;
      while (idx < margin && ch == ' ') {
        idx++;
        ptr++;
        ch = *ptr;
      }
      fprintf (ofp, "%c %s\n", str [0], ptr);
    } else if (ch == '-') {
      fprintf (ofp, "---\n");
    } else if (ch == '=') {
      fprintf (ofp, "===\n");
    }
  }

  fprintf (ofp, "\n");
  fflush (ofp);
}

static void ResetDiffBlock (
  DiffBlockPtr dbp
)

{
  if (dbp == NULL) return;

  dbp->head = ValNodeFreeData (dbp->head);
  dbp->tail = NULL;
}

static void ProcessDiff (
  FILE *ifp,
  FILE *ofp
)

{
  Char       ch;
  DiffBlock  db;
  FileCache  fc;
  Char       line [1024];
  CharPtr    str;

  db.head = NULL;
  db.tail = NULL;

  if (! FileCacheSetup (&fc, ifp)) return;

  str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
  while (str != NULL) {
    ch = line [0];
    if (ch == '<' || ch == '>') {
      RecordDiffBlock (&db, line);
    } else if (ch == '-') {
      RecordDiffBlock (&db, "---");
    } else if (IS_DIGIT (ch)) {
      WriteDiffBlock (&db, ofp);
      ResetDiffBlock (&db);
      RecordDiffBlock (&db, "===");
    } else if (StringHasNoText (str)) {
      WriteDiffBlock (&db, ofp);
      ResetDiffBlock (&db);
      fprintf (ofp, "\n");
    } else {
      fprintf (ofp, "%s\n", line);
    }
    str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
  }

  WriteDiffBlock (&db, ofp);
  ResetDiffBlock (&db);

  fprintf (ofp, "//\n\n");
}

#define i_argInputFile  0
#define o_argOutputFile 1

Args myargs [] = {
  {"Input File", "stdin", NULL, NULL,
    FALSE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  {"Output File", "stdout", NULL, NULL,
    FALSE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
};

Int2 Main (void)

{
  FILE     *ifp, *ofp;
  CharPtr  infile, outfile;

  /* standard setup */

  ErrSetFatalLevel (SEV_MAX);
  ErrClearOptFlags (EO_SHOW_USERSTR);
  ErrPathReset ();

  if (! GetArgs ("diffshift", sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  infile = (CharPtr) myargs [i_argInputFile].strvalue;
  outfile = (CharPtr) myargs [o_argOutputFile].strvalue;

  ifp = FileOpen (infile, "r");
  if (ifp == NULL) {
    Message (MSG_FATAL, "Unable to open input file");
    return 1;
  }

  ofp = FileOpen (outfile, "w");
  if (ofp == NULL) {
    Message (MSG_FATAL, "Unable to open output file");
    return 1;
  }

  ProcessDiff (ifp, ofp);

  FileClose (ofp);
  FileClose (ifp);

  return 0;
}

