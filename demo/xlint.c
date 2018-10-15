/*   xlint.c
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
* File Name:  xlint.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   4/22/11
*
* $Revision: 1.2 $
*
* File Description:
*
*  Lint for XML
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
#include <urlquery.h>

static void WriteStringCallback (XmlObjPtr xop, XmlObjPtr parent, Int2 level, Pointer userdata)

{
  FILE  *ofp;

  if (xop == NULL || userdata == NULL) return;
  ofp = (FILE*) userdata;

  if (StringHasNoText (xop->contents)) return;

  fprintf (ofp, "%s\n", xop->contents);
}

#define i_argInputFile    0
#define o_argOutputFile   1
#define s_argAltSelf      2
#define t_argUseTabs      3
#define n_argNodeFilter   4
#define p_argParentFilter 5
#define a_argAttTagFilter 6
#define v_argAttValFilter 7

Args myargs [] = {
  {"Input File", "stdin", NULL, NULL,
    FALSE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  {"Output File", "stdout", NULL, NULL,
    FALSE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
  {"Self-Closing Tag Alternative Style", "F", NULL, NULL,
    TRUE, 's', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Indent with Tabs instead of Spaces", "F", NULL, NULL,
    TRUE, 't', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Node Filter", NULL, NULL, NULL,
    TRUE, 'n', ARG_STRING, 0.0, 0, NULL},
  {"Parent Filter", NULL, NULL, NULL,
    TRUE, 'p', ARG_STRING, 0.0, 0, NULL},
  {"Attribute Tag Filter", NULL, NULL, NULL,
    TRUE, 'a', ARG_STRING, 0.0, 0, NULL},
  {"Attribute Value Filter", NULL, NULL, NULL,
    TRUE, 'v', ARG_STRING, 0.0, 0, NULL},
};

Int2 Main (void)

{
  Boolean     altSelf, useTabs;
  FileCache   fc;
  FILE        *ifp, *ofp;
  CharPtr     infile, outfile, str;
  Char        line [4096];
  CharPtr     nodeFilt, argFilt, attFilt, valFilt;
  ValNodePtr  head = NULL, tail = NULL;
  XmlObjPtr   xop;

  ErrSetFatalLevel (SEV_MAX);
  ErrClearOptFlags (EO_SHOW_USERSTR);
  ErrPathReset ();

  if (! GetArgs ("xlint", sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  infile = (CharPtr) myargs [i_argInputFile].strvalue;
  outfile = (CharPtr) myargs [o_argOutputFile].strvalue;
  altSelf = (Boolean) myargs [s_argAltSelf].intvalue;
  useTabs = (Boolean) myargs [t_argUseTabs].intvalue;

  nodeFilt = (CharPtr) myargs [n_argNodeFilter].strvalue;
  argFilt = (CharPtr) myargs [p_argParentFilter].strvalue;
  attFilt = (CharPtr) myargs [a_argAttTagFilter].strvalue;
  valFilt = (CharPtr) myargs [v_argAttValFilter].strvalue;

  if (StringHasNoText (infile)) return 1;
  if (StringHasNoText (outfile)) return 1;

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

  FileCacheSetup (&fc, ifp);
  str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
  while (str != NULL) {
    if (StringDoesHaveText (str)) {
      ValNodeCopyStrEx (&head, &tail, 0, str);
    }
    str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
  }

  if (head == NULL) {
    Message (MSG_FATAL, "Unable to read input file");
    return 1;
  }

  str = ValNodeMergeStrs (head);
  ValNodeFreeData (head);

  if (str == NULL) {
    Message (MSG_FATAL, "Unable to merge valnodes");
    return 1;
  }

  xop = ParseXmlString (str);
  if (xop == NULL) {
    Message (MSG_FATAL, "Unable to ParseXmlString");
    return 1;
  }

  if (StringDoesHaveText (nodeFilt) ||
      StringDoesHaveText (argFilt) ||
      StringDoesHaveText (attFilt) ||
      StringDoesHaveText (valFilt)) {
    VisitXmlNodes (xop, (Pointer) ofp, WriteStringCallback, nodeFilt, argFilt, attFilt, valFilt, 0);
  } else {
    WriteXmlObjectEx (xop, ofp, useTabs, altSelf);
  }

  FreeXmlObject (xop);

  MemFree (str);

  FileClose (ofp);

  return 0;
}

