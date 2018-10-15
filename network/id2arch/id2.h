/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "id2.h20";
static AsnValxNode avnx[63] = {
    {20,"set-value" ,1,0.0,&avnx[1] } ,
    {20,"get-value" ,2,0.0,&avnx[2] } ,
    {20,"force-value" ,3,0.0,&avnx[3] } ,
    {20,"use-package" ,4,0.0,NULL } ,
    {3,NULL,1,0.0,NULL } ,
    {20,"any" ,0,0.0,&avnx[6] } ,
    {20,"gi" ,1,0.0,&avnx[7] } ,
    {20,"text" ,2,0.0,&avnx[8] } ,
    {20,"general" ,4,0.0,&avnx[9] } ,
    {20,"all" ,127,0.0,NULL } ,
    {3,NULL,0,0.0,NULL } ,
    {20,"main" ,0,0.0,&avnx[12] } ,
    {20,"snp" ,1,0.0,&avnx[13] } ,
    {20,"snp-graph" ,4,0.0,&avnx[14] } ,
    {20,"cdd" ,8,0.0,&avnx[15] } ,
    {20,"mgc" ,16,0.0,&avnx[16] } ,
    {20,"hprd" ,32,0.0,&avnx[17] } ,
    {20,"sts" ,64,0.0,&avnx[18] } ,
    {20,"trna" ,128,0.0,&avnx[19] } ,
    {20,"exon" ,512,0.0,NULL } ,
    {3,NULL,0,0.0,NULL } ,
    {3,NULL,1,0.0,NULL } ,
    {3,NULL,1,0.0,NULL } ,
    {3,NULL,0,0.0,NULL } ,
    {3,NULL,0,0.0,NULL } ,
    {3,NULL,0,0.0,NULL } ,
    {20,"none" ,0,0.0,&avnx[27] } ,
    {20,"seq-map" ,1,0.0,&avnx[28] } ,
    {20,"seq-data" ,2,0.0,NULL } ,
    {3,NULL,0,0.0,NULL } ,
    {20,"warning" ,1,0.0,&avnx[31] } ,
    {20,"failed-command" ,2,0.0,&avnx[32] } ,
    {20,"failed-connection" ,3,0.0,&avnx[33] } ,
    {20,"failed-server" ,4,0.0,&avnx[34] } ,
    {20,"no-data" ,5,0.0,&avnx[35] } ,
    {20,"restricted-data" ,6,0.0,&avnx[36] } ,
    {20,"unsupported-command" ,7,0.0,&avnx[37] } ,
    {20,"invalid-arguments" ,8,0.0,NULL } ,
    {3,NULL,0,0.0,NULL } ,
    {20,"live" ,0,0.0,&avnx[40] } ,
    {20,"suppressed-temp" ,1,0.0,&avnx[41] } ,
    {20,"suppressed" ,2,0.0,&avnx[42] } ,
    {20,"dead" ,3,0.0,&avnx[43] } ,
    {20,"protected" ,4,0.0,&avnx[44] } ,
    {20,"withdrawn" ,5,0.0,NULL } ,
    {20,"seq-entry" ,0,0.0,&avnx[46] } ,
    {20,"seq-annot" ,1,0.0,&avnx[47] } ,
    {20,"id2s-split-info" ,2,0.0,&avnx[48] } ,
    {20,"id2s-chunk" ,3,0.0,NULL } ,
    {3,NULL,0,0.0,NULL } ,
    {20,"asn-binary" ,0,0.0,&avnx[51] } ,
    {20,"asn-text" ,1,0.0,&avnx[52] } ,
    {20,"xml" ,2,0.0,NULL } ,
    {3,NULL,0,0.0,NULL } ,
    {20,"none" ,0,0.0,&avnx[55] } ,
    {20,"gzip" ,1,0.0,&avnx[56] } ,
    {20,"nlmzip" ,2,0.0,&avnx[57] } ,
    {20,"bzip2" ,3,0.0,NULL } ,
    {3,NULL,0,0.0,NULL } ,
    {3,NULL,0,0.0,NULL } ,
    {20,"reply" ,0,0.0,&avnx[61] } ,
    {20,"last-octet-string" ,1,0.0,&avnx[62] } ,
    {20,"nothing" ,2,0.0,NULL } };

static AsnType atx[148] = {
  {401, "Seq-id" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[1]} ,
  {402, "Seq-loc" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[2]} ,
  {403, "ID2S-Chunk-Id" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[3]} ,
  {404, "ID2S-Seq-annot-Info" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[4]} ,
  {405, "ID2-Request-Packet" ,1,0,0,0,0,0,0,0,NULL,&atx[17],&atx[5],0,&atx[6]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[6],NULL,0,NULL} ,
  {406, "ID2-Request" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[7],0,&atx[10]} ,
  {0, "serial-number" ,128,0,0,1,0,0,0,0,NULL,&atx[8],NULL,0,&atx[9]} ,
  {302, "INTEGER" ,0,2,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "params" ,128,1,0,1,0,0,0,0,NULL,&atx[10],NULL,0,&atx[21]} ,
  {407, "ID2-Params" ,1,0,0,0,0,0,0,0,NULL,&atx[17],&atx[11],0,&atx[25]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[12],NULL,0,NULL} ,
  {431, "ID2-Param" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[13],0,NULL} ,
  {0, "name" ,128,0,0,0,0,0,0,0,NULL,&atx[14],NULL,0,&atx[15]} ,
  {323, "VisibleString" ,0,26,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "value" ,128,1,0,1,0,0,0,0,NULL,&atx[17],&atx[16],0,&atx[18]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[14],NULL,0,NULL} ,
  {312, "SEQUENCE OF" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "type" ,128,2,0,0,1,0,0,0,&avnx[4],&atx[19],&avnx[0],0,NULL} ,
  {310, "ENUMERATED" ,0,10,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {311, "SEQUENCE" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "request" ,128,2,0,0,0,0,0,0,NULL,&atx[35],&atx[22],0,NULL} ,
  {0, "init" ,128,0,0,0,0,0,0,0,NULL,&atx[23],NULL,0,&atx[24]} ,
  {305, "NULL" ,0,5,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "get-packages" ,128,1,0,0,0,0,0,0,NULL,&atx[25],NULL,0,&atx[29]} ,
  {408, "ID2-Request-Get-Packages" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[26],0,&atx[30]} ,
  {0, "names" ,128,0,0,1,0,0,0,0,NULL,&atx[17],&atx[27],0,&atx[28]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[14],NULL,0,NULL} ,
  {0, "no-contents" ,128,1,0,1,0,0,0,0,NULL,&atx[23],NULL,0,NULL} ,
  {0, "get-seq-id" ,128,2,0,0,0,0,0,0,NULL,&atx[30],NULL,0,&atx[37]} ,
  {409, "ID2-Request-Get-Seq-id" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[31],0,&atx[38]} ,
  {0, "seq-id" ,128,0,0,0,0,0,0,0,NULL,&atx[32],NULL,0,&atx[36]} ,
  {414, "ID2-Seq-id" ,1,0,0,0,0,0,0,0,NULL,&atx[35],&atx[33],0,&atx[47]} ,
  {0, "string" ,128,0,0,0,0,0,0,0,NULL,&atx[14],NULL,0,&atx[34]} ,
  {0, "seq-id" ,128,1,0,0,0,0,0,0,NULL,&atx[0],NULL,0,NULL} ,
  {315, "CHOICE" ,0,-1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "seq-id-type" ,128,1,0,0,1,0,0,0,&avnx[10],&atx[8],&avnx[5],0,NULL} ,
  {0, "get-blob-id" ,128,3,0,0,0,0,0,0,NULL,&atx[38],NULL,0,&atx[43]} ,
  {410, "ID2-Request-Get-Blob-Id" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[39],0,&atx[44]} ,
  {0, "seq-id" ,128,0,0,0,0,0,0,0,NULL,&atx[30],NULL,0,&atx[40]} ,
  {0, "sources" ,128,1,0,1,0,0,0,0,NULL,&atx[17],&atx[41],0,&atx[42]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[14],NULL,0,NULL} ,
  {0, "external" ,128,2,0,1,0,0,0,0,NULL,&atx[23],NULL,0,NULL} ,
  {0, "get-blob-info" ,128,4,0,0,0,0,0,0,NULL,&atx[44],NULL,0,&atx[66]} ,
  {411, "ID2-Request-Get-Blob-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[45],0,&atx[67]} ,
  {0, "blob-id" ,128,0,0,0,0,0,0,0,NULL,&atx[35],&atx[46],0,&atx[56]} ,
  {0, "blob-id" ,128,0,0,0,0,0,0,0,NULL,&atx[47],NULL,0,&atx[52]} ,
  {415, "ID2-Blob-Id" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[48],0,&atx[58]} ,
  {0, "sat" ,128,0,0,0,0,0,0,0,NULL,&atx[8],NULL,0,&atx[49]} ,
  {0, "sub-sat" ,128,1,0,0,1,0,0,0,&avnx[20],&atx[8],&avnx[11],0,&atx[50]} ,
  {0, "sat-key" ,128,2,0,0,0,0,0,0,NULL,&atx[8],NULL,0,&atx[51]} ,
  {0, "version" ,128,3,0,1,0,0,0,0,NULL,&atx[8],NULL,0,NULL} ,
  {0, "resolve" ,128,1,0,0,0,0,0,0,NULL,&atx[20],&atx[53],0,NULL} ,
  {0, "request" ,128,0,0,0,0,0,0,0,NULL,&atx[38],NULL,0,&atx[54]} ,
  {0, "exclude-blobs" ,128,1,0,1,0,0,0,0,NULL,&atx[17],&atx[55],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[47],NULL,0,NULL} ,
  {0, "get-seq-ids" ,128,1,0,1,0,0,0,0,NULL,&atx[23],NULL,0,&atx[57]} ,
  {0, "get-data" ,128,2,0,1,0,0,0,0,NULL,&atx[58],NULL,0,NULL} ,
  {416, "ID2-Get-Blob-Details" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[59],0,&atx[77]} ,
  {0, "location" ,128,0,0,1,0,0,0,0,NULL,&atx[1],NULL,0,&atx[60]} ,
  {0, "seq-class-level" ,128,1,0,0,1,0,0,0,&avnx[21],&atx[8],NULL,0,&atx[61]} ,
  {0, "descr-level" ,128,2,0,0,1,0,0,0,&avnx[22],&atx[8],NULL,0,&atx[62]} ,
  {0, "descr-type-mask" ,128,3,0,0,1,0,0,0,&avnx[23],&atx[8],NULL,0,&atx[63]} ,
  {0, "annot-type-mask" ,128,4,0,0,1,0,0,0,&avnx[24],&atx[8],NULL,0,&atx[64]} ,
  {0, "feat-type-mask" ,128,5,0,0,1,0,0,0,&avnx[25],&atx[8],NULL,0,&atx[65]} ,
  {0, "sequence-level" ,128,6,0,0,1,0,0,0,&avnx[29],&atx[19],&avnx[26],0,NULL} ,
  {0, "reget-blob" ,128,5,0,0,0,0,0,0,NULL,&atx[67],NULL,0,&atx[71]} ,
  {412, "ID2-Request-ReGet-Blob" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[68],0,&atx[72]} ,
  {0, "blob-id" ,128,0,0,0,0,0,0,0,NULL,&atx[47],NULL,0,&atx[69]} ,
  {0, "split-version" ,128,1,0,0,0,0,0,0,NULL,&atx[8],NULL,0,&atx[70]} ,
  {0, "offset" ,128,2,0,0,0,0,0,0,NULL,&atx[8],NULL,0,NULL} ,
  {0, "get-chunks" ,128,6,0,0,0,0,0,0,NULL,&atx[72],NULL,0,NULL} ,
  {413, "ID2S-Request-Get-Chunks" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[73],0,&atx[32]} ,
  {0, "blob-id" ,128,0,0,0,0,0,0,0,NULL,&atx[47],NULL,0,&atx[74]} ,
  {0, "chunks" ,128,1,0,0,0,0,0,0,NULL,&atx[17],&atx[75],0,&atx[76]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "split-version" ,128,2,0,1,0,0,0,0,NULL,&atx[8],NULL,0,NULL} ,
  {417, "ID2-Reply" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[78],0,&atx[82]} ,
  {0, "serial-number" ,128,0,0,1,0,0,0,0,NULL,&atx[8],NULL,0,&atx[79]} ,
  {0, "params" ,128,1,0,1,0,0,0,0,NULL,&atx[10],NULL,0,&atx[80]} ,
  {0, "error" ,128,2,0,1,0,0,0,0,NULL,&atx[17],&atx[81],0,&atx[86]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[82],NULL,0,NULL} ,
  {418, "ID2-Error" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[83],0,&atx[91]} ,
  {0, "severity" ,128,0,0,0,0,0,0,0,NULL,&atx[19],&avnx[30],0,&atx[84]} ,
  {0, "retry-delay" ,128,1,0,1,0,0,0,0,NULL,&atx[8],NULL,0,&atx[85]} ,
  {0, "message" ,128,2,0,1,0,0,0,0,NULL,&atx[14],NULL,0,NULL} ,
  {0, "end-of-reply" ,128,3,0,1,0,0,0,0,NULL,&atx[23],NULL,0,&atx[87]} ,
  {0, "reply" ,128,4,0,0,0,0,0,0,NULL,&atx[35],&atx[88],0,&atx[142]} ,
  {0, "init" ,128,0,0,0,0,0,0,0,NULL,&atx[23],NULL,0,&atx[89]} ,
  {0, "empty" ,128,1,0,0,0,0,0,0,NULL,&atx[23],NULL,0,&atx[90]} ,
  {0, "get-package" ,128,2,0,0,0,0,0,0,NULL,&atx[91],NULL,0,&atx[94]} ,
  {419, "ID2-Reply-Get-Package" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[92],0,&atx[95]} ,
  {0, "name" ,128,0,0,0,0,0,0,0,NULL,&atx[14],NULL,0,&atx[93]} ,
  {0, "params" ,128,1,0,1,0,0,0,0,NULL,&atx[10],NULL,0,NULL} ,
  {0, "get-seq-id" ,128,3,0,0,0,0,0,0,NULL,&atx[95],NULL,0,&atx[100]} ,
  {420, "ID2-Reply-Get-Seq-id" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[96],0,&atx[101]} ,
  {0, "request" ,128,0,0,0,0,0,0,0,NULL,&atx[30],NULL,0,&atx[97]} ,
  {0, "seq-id" ,128,1,0,1,0,0,0,0,NULL,&atx[17],&atx[98],0,&atx[99]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[0],NULL,0,NULL} ,
  {0, "end-of-reply" ,128,2,0,1,0,0,0,0,NULL,&atx[23],NULL,0,NULL} ,
  {0, "get-blob-id" ,128,4,0,0,0,0,0,0,NULL,&atx[101],NULL,0,&atx[110]} ,
  {421, "ID2-Reply-Get-Blob-Id" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[102],0,&atx[111]} ,
  {0, "seq-id" ,128,0,0,0,0,0,0,0,NULL,&atx[0],NULL,0,&atx[103]} ,
  {0, "blob-id" ,128,1,0,1,0,0,0,0,NULL,&atx[47],NULL,0,&atx[104]} ,
  {0, "split-version" ,128,2,0,0,1,0,0,0,&avnx[38],&atx[8],NULL,0,&atx[105]} ,
  {0, "annot-info" ,128,3,0,1,0,0,0,0,NULL,&atx[17],&atx[106],0,&atx[107]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[3],NULL,0,NULL} ,
  {0, "end-of-reply" ,128,4,0,1,0,0,0,0,NULL,&atx[23],NULL,0,&atx[108]} ,
  {0, "blob-state" ,128,5,0,1,0,0,0,0,NULL,&atx[109],NULL,0,NULL} ,
  {427, "ID2-Blob-State" ,1,0,0,0,0,0,0,0,NULL,&atx[19],&avnx[39],0,&atx[114]} ,
  {0, "get-blob-seq-ids" ,128,5,0,0,0,0,0,0,NULL,&atx[111],NULL,0,&atx[121]} ,
  {422, "ID2-Reply-Get-Blob-Seq-ids" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[112],0,&atx[122]} ,
  {0, "blob-id" ,128,0,0,0,0,0,0,0,NULL,&atx[47],NULL,0,&atx[113]} ,
  {0, "ids" ,128,1,0,1,0,0,0,0,NULL,&atx[114],NULL,0,NULL} ,
  {428, "ID2-Reply-Data" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[115],0,&atx[143]} ,
  {0, "data-type" ,128,0,0,0,1,0,0,0,&avnx[49],&atx[8],&avnx[45],0,&atx[116]} ,
  {0, "data-format" ,128,1,0,0,1,0,0,0,&avnx[53],&atx[8],&avnx[50],0,&atx[117]} ,
  {0, "data-compression" ,128,2,0,0,1,0,0,0,&avnx[58],&atx[8],&avnx[54],0,&atx[118]} ,
  {0, "data" ,128,3,0,0,0,0,0,0,NULL,&atx[17],&atx[119],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[120],NULL,0,NULL} ,
  {304, "OCTET STRING" ,0,4,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "get-blob" ,128,6,0,0,0,0,0,0,NULL,&atx[122],NULL,0,&atx[126]} ,
  {423, "ID2-Reply-Get-Blob" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[123],0,&atx[127]} ,
  {0, "blob-id" ,128,0,0,0,0,0,0,0,NULL,&atx[47],NULL,0,&atx[124]} ,
  {0, "split-version" ,128,1,0,0,1,0,0,0,&avnx[59],&atx[8],NULL,0,&atx[125]} ,
  {0, "data" ,128,2,0,1,0,0,0,0,NULL,&atx[114],NULL,0,NULL} ,
  {0, "reget-blob" ,128,7,0,0,0,0,0,0,NULL,&atx[127],NULL,0,&atx[132]} ,
  {424, "ID2-Reply-ReGet-Blob" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[128],0,&atx[133]} ,
  {0, "blob-id" ,128,0,0,0,0,0,0,0,NULL,&atx[47],NULL,0,&atx[129]} ,
  {0, "split-version" ,128,1,0,0,0,0,0,0,NULL,&atx[8],NULL,0,&atx[130]} ,
  {0, "offset" ,128,2,0,0,0,0,0,0,NULL,&atx[8],NULL,0,&atx[131]} ,
  {0, "data" ,128,3,0,1,0,0,0,0,NULL,&atx[114],NULL,0,NULL} ,
  {0, "get-split-info" ,128,8,0,0,0,0,0,0,NULL,&atx[133],NULL,0,&atx[137]} ,
  {425, "ID2S-Reply-Get-Split-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[134],0,&atx[138]} ,
  {0, "blob-id" ,128,0,0,0,0,0,0,0,NULL,&atx[47],NULL,0,&atx[135]} ,
  {0, "split-version" ,128,1,0,0,0,0,0,0,NULL,&atx[8],NULL,0,&atx[136]} ,
  {0, "data" ,128,2,0,1,0,0,0,0,NULL,&atx[114],NULL,0,NULL} ,
  {0, "get-chunk" ,128,9,0,0,0,0,0,0,NULL,&atx[138],NULL,0,NULL} ,
  {426, "ID2S-Reply-Get-Chunk" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[139],0,&atx[109]} ,
  {0, "blob-id" ,128,0,0,0,0,0,0,0,NULL,&atx[47],NULL,0,&atx[140]} ,
  {0, "chunk-id" ,128,1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[141]} ,
  {0, "data" ,128,2,0,1,0,0,0,0,NULL,&atx[114],NULL,0,NULL} ,
  {0, "discard" ,128,5,0,1,0,0,0,0,NULL,&atx[19],&avnx[60],0,NULL} ,
  {429, "ID2-Blob-Seq-ids" ,1,0,0,0,0,0,0,0,NULL,&atx[17],&atx[144],0,&atx[145]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[145],NULL,0,NULL} ,
  {430, "ID2-Blob-Seq-id" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[146],0,&atx[12]} ,
  {0, "seq-id" ,128,0,0,0,0,0,0,0,NULL,&atx[0],NULL,0,&atx[147]} ,
  {0, "replaced" ,128,1,0,1,0,0,0,0,NULL,&atx[23],NULL,0,NULL} };

static AsnModule ampx[1] = {
  { "NCBI-ID2Access" , "id2.h20",&atx[0],NULL,NULL,0,0} };

static AsnValxNodePtr avn = avnx;
static AsnTypePtr at = atx;
static AsnModulePtr amp = ampx;



/**************************************************
*
*    Defines for Module NCBI-ID2Access
*
**************************************************/

#define ID2_REQUEST_PACKET &at[4]
#define ID2_REQUEST_PACKET_E &at[5]

#define ID2_REQUEST &at[6]
#define ID2_REQUEST_serial_number &at[7]
#define ID2_REQUEST_params &at[9]
#define ID2_REQUEST_request &at[21]
#define ID2_REQUEST_request_init &at[22]
#define REQUEST_request_get_packages &at[24]
#define ID2_REQUEST_request_get_seq_id &at[29]
#define ID2_REQUEST_request_get_blob_id &at[37]
#define REQUEST_request_get_blob_info &at[43]
#define ID2_REQUEST_request_reget_blob &at[66]
#define ID2_REQUEST_request_get_chunks &at[71]

#define ID2_PARAMS &at[10]
#define ID2_PARAMS_E &at[11]

#define ID2_REQUEST_GET_PACKAGES &at[25]
#define ID2_REQUEST_GET_PACKAGES_names &at[26]
#define REQUEST_GET_PACKAGES_names_E &at[27]
#define GET_PACKAGES_no_contents &at[28]

#define ID2_REQUEST_GET_SEQ_ID &at[30]
#define ID2_REQUEST_GET_SEQ_ID_seq_id &at[31]
#define REQUEST_GET_SEQ_ID_seq_id_type &at[36]

#define ID2_REQUEST_GET_BLOB_ID &at[38]
#define ID2_REQUEST_GET_BLOB_ID_seq_id &at[39]
#define ID2_REQUEST_GET_BLOB_ID_sources &at[40]
#define REQUEST_GET_BLOB_ID_sources_E &at[41]
#define REQUEST_GET_BLOB_ID_external &at[42]

#define ID2_REQUEST_GET_BLOB_INFO &at[44]
#define REQUEST_GET_BLOB_INFO_blob_id &at[45]
#define GET_BLOB_INFO_blob_id_blob_id &at[46]
#define GET_BLOB_INFO_blob_id_resolve &at[52]
#define INFO_blob_id_resolve_request &at[53]
#define blob_id_resolve_exclude_blobs &at[54]
#define id_resolve_exclude_blobs_E &at[55]
#define GET_BLOB_INFO_get_seq_ids &at[56]
#define REQUEST_GET_BLOB_INFO_get_data &at[57]

#define ID2_REQUEST_REGET_BLOB &at[67]
#define ID2_REQUEST_REGET_BLOB_blob_id &at[68]
#define REGET_BLOB_split_version &at[69]
#define ID2_REQUEST_REGET_BLOB_offset &at[70]

#define ID2S_REQUEST_GET_CHUNKS &at[72]
#define ID2S_REQUEST_GET_CHUNKS_blob_id &at[73]
#define ID2S_REQUEST_GET_CHUNKS_chunks &at[74]
#define REQUEST_GET_CHUNKS_chunks_E &at[75]
#define GET_CHUNKS_split_version &at[76]

#define ID2_SEQ_ID &at[32]
#define ID2_SEQ_ID_string &at[33]
#define ID2_SEQ_ID_seq_id &at[34]

#define ID2_BLOB_ID &at[47]
#define ID2_BLOB_ID_sat &at[48]
#define ID2_BLOB_ID_sub_sat &at[49]
#define ID2_BLOB_ID_sat_key &at[50]
#define ID2_BLOB_ID_version &at[51]

#define ID2_GET_BLOB_DETAILS &at[58]
#define ID2_GET_BLOB_DETAILS_location &at[59]
#define BLOB_DETAILS_seq_class_level &at[60]
#define GET_BLOB_DETAILS_descr_level &at[61]
#define BLOB_DETAILS_descr_type_mask &at[62]
#define BLOB_DETAILS_annot_type_mask &at[63]
#define BLOB_DETAILS_feat_type_mask &at[64]
#define BLOB_DETAILS_sequence_level &at[65]

#define ID2_REPLY &at[77]
#define ID2_REPLY_serial_number &at[78]
#define ID2_REPLY_params &at[79]
#define ID2_REPLY_error &at[80]
#define ID2_REPLY_error_E &at[81]
#define ID2_REPLY_end_of_reply &at[86]
#define ID2_REPLY_reply &at[87]
#define ID2_REPLY_reply_init &at[88]
#define ID2_REPLY_reply_empty &at[89]
#define ID2_REPLY_reply_get_package &at[90]
#define ID2_REPLY_reply_get_seq_id &at[94]
#define ID2_REPLY_reply_get_blob_id &at[100]
#define REPLY_reply_get_blob_seq_ids &at[110]
#define ID2_REPLY_reply_get_blob &at[121]
#define ID2_REPLY_reply_reget_blob &at[126]
#define ID2_REPLY_reply_get_split_info &at[132]
#define ID2_REPLY_reply_get_chunk &at[137]
#define ID2_REPLY_discard &at[142]

#define ID2_ERROR &at[82]
#define ID2_ERROR_severity &at[83]
#define ID2_ERROR_retry_delay &at[84]
#define ID2_ERROR_message &at[85]

#define ID2_REPLY_GET_PACKAGE &at[91]
#define ID2_REPLY_GET_PACKAGE_name &at[92]
#define ID2_REPLY_GET_PACKAGE_params &at[93]

#define ID2_REPLY_GET_SEQ_ID &at[95]
#define ID2_REPLY_GET_SEQ_ID_request &at[96]
#define ID2_REPLY_GET_SEQ_ID_seq_id &at[97]
#define ID2_REPLY_GET_SEQ_ID_seq_id_E &at[98]
#define REPLY_GET_SEQ_ID_end_of_reply &at[99]

#define ID2_REPLY_GET_BLOB_ID &at[101]
#define ID2_REPLY_GET_BLOB_ID_seq_id &at[102]
#define ID2_REPLY_GET_BLOB_ID_blob_id &at[103]
#define GET_BLOB_ID_split_version &at[104]
#define REPLY_GET_BLOB_ID_annot_info &at[105]
#define REPLY_GET_BLOB_ID_annot_info_E &at[106]
#define REPLY_GET_BLOB_ID_end_of_reply &at[107]
#define REPLY_GET_BLOB_ID_blob_state &at[108]

#define ID2_REPLY_GET_BLOB_SEQ_IDS &at[111]
#define REPLY_GET_BLOB_SEQ_IDS_blob_id &at[112]
#define ID2_REPLY_GET_BLOB_SEQ_IDS_ids &at[113]

#define ID2_REPLY_GET_BLOB &at[122]
#define ID2_REPLY_GET_BLOB_blob_id &at[123]
#define REPLY_GET_BLOB_split_version &at[124]
#define ID2_REPLY_GET_BLOB_data &at[125]

#define ID2_REPLY_REGET_BLOB &at[127]
#define ID2_REPLY_REGET_BLOB_blob_id &at[128]
#define REPLY_REGET_BLOB_split_version &at[129]
#define ID2_REPLY_REGET_BLOB_offset &at[130]
#define ID2_REPLY_REGET_BLOB_data &at[131]

#define ID2S_REPLY_GET_SPLIT_INFO &at[133]
#define REPLY_GET_SPLIT_INFO_blob_id &at[134]
#define GET_SPLIT_INFO_split_version &at[135]
#define ID2S_REPLY_GET_SPLIT_INFO_data &at[136]

#define ID2S_REPLY_GET_CHUNK &at[138]
#define ID2S_REPLY_GET_CHUNK_blob_id &at[139]
#define ID2S_REPLY_GET_CHUNK_chunk_id &at[140]
#define ID2S_REPLY_GET_CHUNK_data &at[141]

#define ID2_BLOB_STATE &at[109]

#define ID2_REPLY_DATA &at[114]
#define ID2_REPLY_DATA_data_type &at[115]
#define ID2_REPLY_DATA_data_format &at[116]
#define ID2_REPLY_DATA_data_compression &at[117]
#define ID2_REPLY_DATA_data &at[118]
#define ID2_REPLY_DATA_data_E &at[119]

#define ID2_BLOB_SEQ_IDS &at[143]
#define ID2_BLOB_SEQ_IDS_E &at[144]

#define ID2_BLOB_SEQ_ID &at[145]
#define ID2_BLOB_SEQ_ID_seq_id &at[146]
#define ID2_BLOB_SEQ_ID_replaced &at[147]

#define ID2_PARAM &at[12]
#define ID2_PARAM_name &at[13]
#define ID2_PARAM_value &at[15]
#define ID2_PARAM_value_E &at[16]
#define ID2_PARAM_type &at[18]
