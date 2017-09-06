!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve a coupled Monodomain equation Finite Elasticity problem using OpenCMISS calls.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s): Thomas Heidlauf
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!> \example MultiPhysics/BioelectricFiniteElasticity/GudunovMonodomainElasticity1D3DMesh/EntireTA/src/EntireTAExample.f90
!! Example program to solve a coupled Finite Elasticity Monodomain equation using OpenCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/MultiPhysics/BioelectricFiniteElasticity/GudunovMonodomainElasticity1D3DMesh/EntireTA/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/MultiPhysics/BioelectricFiniteElasticity/GudunovMonodomainElasticity1D3DMesh/EntireTA/build-gnu'>Linux GNU Build</a>
!!
!<

!> Main program
PROGRAM EntireTAEXAMPLE

  USE OpenCMISS
  USE OpenCMISS_Iron
#ifndef NOMPIMOD
  USE MPI
#endif

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

  !--------------------------------------------------------------------------------------------------------------------------------
  !Test program parameters
  REAL(CMISSRP), PARAMETER :: tol=1.0E-8_CMISSRP

  REAL :: elapsed(2)     ! For receiving user and system time
!  REAL :: total 

  REAL(CMISSRP), PARAMETER :: P_max=7.5_CMISSRP ! N/cm^2  !with new CellML file
!  REAL(CMISSRP), PARAMETER :: P_max=0.0_CMISSRP ! test (wrong)

  LOGICAL, PARAMETER :: biophysical=.FALSE.
  LOGICAL, PARAMETER :: independent_field_auto_create=.FALSE.

  INTEGER(CMISSIntg), PARAMETER :: NumberOfElementsFE=39
  INTEGER(CMISSIntg), PARAMETER :: NumberOfNodesFE=543

!  INTEGER(CMISSIntg), PARAMETER :: NumberOfElementsM=22784!7424
!  INTEGER(CMISSIntg), PARAMETER :: NumberOfNodesM=23040!7680
!  INTEGER(CMISSIntg), PARAMETER :: NumberOfElementsM=58!22784!58!29
!  INTEGER(CMISSIntg), PARAMETER :: NumberOfNodesM=60!23040!60!30
  INTEGER(CMISSIntg) :: NumberOfElementsM,NumberOfElementsMmuscle,NumberOfElementsMskin
  INTEGER(CMISSIntg) :: NumberOfNodesM,NumberOfNodesMmuscle,NumberOfNodesMskin
  INTEGER(CMISSIntg) :: NumberOfNodesPerFibre
  INTEGER(CMISSIntg) :: NumberOfElementsPerElasticityElement
  INTEGER(CMISSIntg), PARAMETER :: NumberOfInSeriesFibres=1

!  !STA
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element1= &
!   & [6,21,103,31,33,155,8,22,104,46,47,182,48,49,184,70,71,230,13,23,125,32,34,156,15,24,127]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element2= &
!   & [13,23,125,32,34,156,15,24,127,54,55,202,56,57,204,74,75,240,9,25,121,35,37,169,11,26,123]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element3= &
!   & [9,25,121,35,37,169,11,26,123,50,51,198,52,53,200,72,73,238,17,27,129,36,38,170,19,28,131]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element4= &
!   & [17,27,129,36,38,170,19,28,131,58,59,206,60,61,208,76,77,242,120,194,101,166,196,165,117,236,102]

!  !DTA
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element5= &
!   & [103,137,105,155,157,173,104,138,106,182,183,214,184,185,215,230,231,246,125,139,126,156,158,174,127,140,128]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element6= &
!   & [109,141,113,159,161,175,110,142,114,186,187,216,188,189,217,232,233,247,133,143,134,160,162,176,135,144,136]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element7= &
!   & [111,145,115,163,164,177,112,146,116,190,191,218,192,193,219,234,235,248,103,137,105,155,157,173,104,138,106]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element8= &
!   & [101,147,107,165,167,178,102,148,108,194,195,220,196,197,221,236,237,249,120,149,119,166,168,179,117,150,118]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element9= &
!   & [121,151,122,169,171,180,123,152,124,198,199,222,200,201,223,238,239,250,129,153,130,170,172,181,131,154,132]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element10= &
!   & [125,139,126,156,158,174,127,140,128,202,203,224,204,205,225,240,241,251,121,151,122,169,171,180,123,152,124]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element11= &
!   & [129,153,130,170,172,181,131,154,132,206,207,226,208,209,227,242,243,252,101,147,107,165,167,178,102,148,108]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element12= &
!   & [133,143,134,160,162,176,135,144,136,210,211,228,212,213,229,244,245,253,111,145,115,163,164,177,112,146,116]

!  !SKIN
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element13= &
!   & [1017,1075,111,1110,1112,163,1068,1076,112,1190,1191,190,1192,1193,192,1238, &
!   &  1239,234,1020,1073,103,1109,1111,155,1070,1074,104]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element14= &
!   & [1051,1077,1001,1113,1115,1165,1033,1078,1035,1194,1195,1288,1196,1197,1289, &
!   &  1334,1335,1350,1030,1079,1026,1114,1116,1166,1034,1080,1016]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element15= &
!   & [1065,1081,1015,1117,1118,1167,1066,1082,1027,1198,1199,1290,1200,1201,1291, &
!   &  1336,1337,1351,1051,1077,1001,1113,1115,1165,1033,1078,1035]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element16= &
!   & [1059,1083,1009,1119,1120,1168,1031,1084,1029,1202,1203,1292,1204,1205,1293, &
!   &  1338,1339,1352,1065,1081,1015,1117,1118,1167,1066,1082,1027]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element17= &
!   & [1062,1085,1012,1121,1122,1169,1019,1086,1055,1206,1207,1294,1208,1209,1295, &
!   &  1340,1341,1353,1059,1083,1009,1119,1120,1168,1031,1084,1029]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element18= &
!   & [1058,1088,11,1124,1126,1171,1059,1083,1009,1210,1211,72, &
!   &  1212,1213,1297,1202,1203,1292,1064,1087,19,1123,1125,1170,1065,1081,1015]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element19= &
!   & [1020,1073,103,1109,1111,155,1070,1074,104,1214,1215,21,1216,1217,33, &
!   &  1230,1231,22,1054,1089,6,1127,1128,31,1056,1090,8]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element20= &
!   & [1061,1091,15,1129,1130,1173,1062,1085,1012,1218,1219,74, &
!   &  1220,1221,1301,1206,1207,1294,1058,1088,11,1124,1126,1171,1059,1083,1009]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element21= &
!   & [1056,1090,8,1131,1132,1174,1018,1092,1028,1222,1223,70, &
!   &  1224,1225,1303,1226,1227,1304,1061,1091,15,1129,1130,1173,1062,1085,1012]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element22= &
!   & [1018,1092,1028,1133,1134,1175,1032,1093,1025,1226,1227,1304, &
!   &  1228,1229,1305,1342,1343,1354,1062,1085,1012,1121,1122,1169,1019,1086,1055]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element23= &
!   & [1070,1074,104,1135,1136,1176,1053,1094,1024,1230,1231,22, &
!   &  1232,1233,1307,1234,1235,1308,1056,1090,8,1131,1132,1174,1018,1092,1028]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element24= &
!   & [1053,1094,1024,1137,1138,1177,1037,1095,1041,1234,1235,1308, &
!   &  1236,1237,1309,1344,1345,1355,1018,1092,1028,1133,1134,1175,1032,1093,1025]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element25= &
!   & [1068,1076,112,1139,1140,1178,1036,1096,1039,1238,1239,234, &
!   &  1240,1241,1311,1242,1243,1312,1070,1074,104,1135,1136,1176,1053,1094,1024]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element26= &
!   & [1036,1096,1039,1141,1142,1179,1040,1097,1038,1242,1243,1312, &
!   &  1244,1245,1313,1346,1347,1356,1053,1094,1024,1137,1138,1177,1037,1095,1041]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element27= &
!   & [1246,1247,133,1248,1249,160,1250,1251,135,1358,1360,210, &
!   &  1364,1366,212,1370,1372,244,1017,1075,111,1110,1112,163,1068,1076,112]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element28= &
!   & [1250,1251,135,1252,1253,1317,1254,1255,1318,1370,1372,244, &
!   &  1376,1378,1380,1382,1384,1386,1068,1076,112,1139,1140,1178,1036,1096,1039]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element29= &
!   & [1254,1255,1318,1256,1257,1319,1348,1349,1357,1382,1384,1386, &
!   &  1388,1390,1392,1394,1396,1398,1036,1096,1039,1141,1142,1179,1040,1097,1038]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element30= &
!   & [1054,1089,6,1127,1128,31,1056,1090,8,1258,1259,46,1260,1261,48, &
!   &  1222,1223,70,1060,1102,13,1149,1150,32,1061,1091,15]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element31= &
!   & [1060,1102,13,1149,1150,32,1061,1091,15,1262,1263,54,1264,1265,56, &
!   &  1218,1219,74,1057,1103,9,1151,1152,35,1058,1088,11]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element32= &
!   & [1057,1103,9,1151,1152,35,1058,1088,11,1266,1267,50,1268,1269,52, &
!   &  1210,1211,72,1063,1104,17,1153,1154,36,1064,1087,19]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element33= &
!   & [1063,1104,17,1153,1154,36,1064,1087,19,1270,1271,58,1272,1273,60, &
!   &  1278,1279,76,1050,1105,120,1155,1156,166,1052,1106,117]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element34= &
!   & [1050,1105,120,1155,1156,166,1052,1106,117,1274,1275,149,1276,1277,168, &
!   &  1282,1283,150,1072,1107,119,1157,1158,179,1071,1108,118]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element35= &
!   & [1064,1087,19,1123,1125,1170,1065,1081,1015,1278,1279,76,1280,1281,1331, &
!   &  1198,1199,1290,1052,1106,117,1159,1160,1188,1051,1077,1001]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element36= &
!   & [1052,1106,117,1159,1160,1188,1051,1077,1001,1282,1283,150,1284,1285,1333, &
!   &  1194,1195,1288,1071,1108,118,1161,1162,1189,1030,1079,1026]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element37= &
!   & [1043,1098,109,1143,1144,159,1049,1099,110,1359,1361,186, &
!   &  1365,1367,188,1371,1373,232,1246,1247,133,1248,1249,160,1250,1251,135]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element38= &
!   & [1049,1099,110,1145,1146,1181,1044,1100,1046,1371,1373,232, &
!   &  1377,1379,1381,1383,1385,1387,1250,1251,135,1252,1253,1317,1254,1255,1318]
!  INTEGER(CMISSIntg), DIMENSION(27), PARAMETER :: Element39= &
!   & [1044,1100,1046,1147,1148,1182,1047,1101,1045,1383,1385,1387, &
!   &  1389,1391,1393,1395,1397,1399,1254,1255,1318,1256,1257,1319,1348,1349,1357]

  INTEGER(CMISSIntg), DIMENSION(NumberOfElementsFE,27), PARAMETER :: AllElementNodes=transpose(reshape( &
   & [6,21,103,31,33,155,8,22,104,46,47,182,48,49,184,70,71,230,13,23,125,32,34,156,15,24,127, &
   &  13,23,125,32,34,156,15,24,127,54,55,202,56,57,204,74,75,240,9,25,121,35,37,169,11,26,123, &
   &  9,25,121,35,37,169,11,26,123,50,51,198,52,53,200,72,73,238,17,27,129,36,38,170,19,28,131, &
   &  17,27,129,36,38,170,19,28,131,58,59,206,60,61,208,76,77,242,120,194,101,166,196,165,117,236,102, &
   &  103,137,105,155,157,173,104,138,106,182,183,214,184,185,215,230,231,246,125,139,126,156,158,174,127,140,128, &
   &  109,141,113,159,161,175,110,142,114,186,187,216,188,189,217,232,233,247,133,143,134,160,162,176,135,144,136, &
   &  111,145,115,163,164,177,112,146,116,190,191,218,192,193,219,234,235,248,103,137,105,155,157,173,104,138,106, &
   &  101,147,107,165,167,178,102,148,108,194,195,220,196,197,221,236,237,249,120,149,119,166,168,179,117,150,118, &
   &  121,151,122,169,171,180,123,152,124,198,199,222,200,201,223,238,239,250,129,153,130,170,172,181,131,154,132, &
   &  125,139,126,156,158,174,127,140,128,202,203,224,204,205,225,240,241,251,121,151,122,169,171,180,123,152,124, &
   &  129,153,130,170,172,181,131,154,132,206,207,226,208,209,227,242,243,252,101,147,107,165,167,178,102,148,108, &
   &  133,143,134,160,162,176,135,144,136,210,211,228,212,213,229,244,245,253,111,145,115,163,164,177,112,146,116, &
   &  1017,1075,111,1110,1112,163,1068,1076,112,1190,1191,190,1192,1193,192,1238, &
   &  1239,234,1020,1073,103,1109,1111,155,1070,1074,104, &
   &  1051,1077,1001,1113,1115,1165,1033,1078,1035,1194,1195,1288,1196,1197,1289, &
   &  1334,1335,1350,1030,1079,1026,1114,1116,1166,1034,1080,1016, &
   &  1065,1081,1015,1117,1118,1167,1066,1082,1027,1198,1199,1290,1200,1201,1291, &
   &  1336,1337,1351,1051,1077,1001,1113,1115,1165,1033,1078,1035, &
   &  1059,1083,1009,1119,1120,1168,1031,1084,1029,1202,1203,1292,1204,1205,1293, &
   &  1338,1339,1352,1065,1081,1015,1117,1118,1167,1066,1082,1027, &
   &  1062,1085,1012,1121,1122,1169,1019,1086,1055,1206,1207,1294,1208,1209,1295, &
   &  1340,1341,1353,1059,1083,1009,1119,1120,1168,1031,1084,1029, &
   &  1058,1088,11,1124,1126,1171,1059,1083,1009,1210,1211,72, &
   &  1212,1213,1297,1202,1203,1292,1064,1087,19,1123,1125,1170,1065,1081,1015, &
   &  1020,1073,103,1109,1111,155,1070,1074,104,1214,1215,21,1216,1217,33, &
   &  1230,1231,22,1054,1089,6,1127,1128,31,1056,1090,8, &
   &  1061,1091,15,1129,1130,1173,1062,1085,1012,1218,1219,74, &
   &  1220,1221,1301,1206,1207,1294,1058,1088,11,1124,1126,1171,1059,1083,1009, &
   &  1056,1090,8,1131,1132,1174,1018,1092,1028,1222,1223,70, &
   &  1224,1225,1303,1226,1227,1304,1061,1091,15,1129,1130,1173,1062,1085,1012, &
   &  1018,1092,1028,1133,1134,1175,1032,1093,1025,1226,1227,1304, &
   &  1228,1229,1305,1342,1343,1354,1062,1085,1012,1121,1122,1169,1019,1086,1055, &
   &  1070,1074,104,1135,1136,1176,1053,1094,1024,1230,1231,22, &
   &  1232,1233,1307,1234,1235,1308,1056,1090,8,1131,1132,1174,1018,1092,1028, &
   &  1053,1094,1024,1137,1138,1177,1037,1095,1041,1234,1235,1308, &
   &  1236,1237,1309,1344,1345,1355,1018,1092,1028,1133,1134,1175,1032,1093,1025, &
   &  1068,1076,112,1139,1140,1178,1036,1096,1039,1238,1239,234, &
   &  1240,1241,1311,1242,1243,1312,1070,1074,104,1135,1136,1176,1053,1094,1024, &
   &  1036,1096,1039,1141,1142,1179,1040,1097,1038,1242,1243,1312, &
   &  1244,1245,1313,1346,1347,1356,1053,1094,1024,1137,1138,1177,1037,1095,1041, &
   &  1246,1247,133,1248,1249,160,1250,1251,135,1358,1360,210, &
   &  1364,1366,212,1370,1372,244,1017,1075,111,1110,1112,163,1068,1076,112, &
   &  1250,1251,135,1252,1253,1317,1254,1255,1318,1370,1372,244, &
   &  1376,1378,1380,1382,1384,1386,1068,1076,112,1139,1140,1178,1036,1096,1039, &
   &  1254,1255,1318,1256,1257,1319,1348,1349,1357,1382,1384,1386, &
   &  1388,1390,1392,1394,1396,1398,1036,1096,1039,1141,1142,1179,1040,1097,1038, &
   &  1054,1089,6,1127,1128,31,1056,1090,8,1258,1259,46,1260,1261,48, &
   &  1222,1223,70,1060,1102,13,1149,1150,32,1061,1091,15, &
   &  1060,1102,13,1149,1150,32,1061,1091,15,1262,1263,54,1264,1265,56, &
   &  1218,1219,74,1057,1103,9,1151,1152,35,1058,1088,11, &
   &  1057,1103,9,1151,1152,35,1058,1088,11,1266,1267,50,1268,1269,52, &
   &  1210,1211,72,1063,1104,17,1153,1154,36,1064,1087,19, &
   &  1063,1104,17,1153,1154,36,1064,1087,19,1270,1271,58,1272,1273,60, &
   &  1278,1279,76,1050,1105,120,1155,1156,166,1052,1106,117, &
   &  1050,1105,120,1155,1156,166,1052,1106,117,1274,1275,149,1276,1277,168, &
   &  1282,1283,150,1072,1107,119,1157,1158,179,1071,1108,118, &
   &  1064,1087,19,1123,1125,1170,1065,1081,1015,1278,1279,76,1280,1281,1331, &
   &  1198,1199,1290,1052,1106,117,1159,1160,1188,1051,1077,1001, &
   &  1052,1106,117,1159,1160,1188,1051,1077,1001,1282,1283,150,1284,1285,1333, &
   &  1194,1195,1288,1071,1108,118,1161,1162,1189,1030,1079,1026, &
   &  1043,1098,109,1143,1144,159,1049,1099,110,1359,1361,186, &
   &  1365,1367,188,1371,1373,232,1246,1247,133,1248,1249,160,1250,1251,135, &
   &  1049,1099,110,1145,1146,1181,1044,1100,1046,1371,1373,232, &
   &  1377,1379,1381,1383,1385,1387,1250,1251,135,1252,1253,1317,1254,1255,1318, &
   &  1044,1100,1046,1147,1148,1182,1047,1101,1045,1383,1385,1387, &
   &  1389,1391,1393,1395,1397,1399,1254,1255,1318,1256,1257,1319,1348,1349,1357], &
   & [27,NumberOfElementsFE]))

  INTEGER(CMISSIntg), DIMENSION(NumberOfNodesFE), PARAMETER :: AllNodes= &
   & [6,8,9,11,13,15,17,19,21,22,23,24,25,26,27,28,31,32,33,34,35,36,37,38,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,70,71, &
   & 72,73,74,75,76,77,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126, &
   & 127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157, &
   & 158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188, &
   & 189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219, &
   & 220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250, &
   & 251,252,253,1001,1009,1012,1015,1016,1017,1018,1019,1020,1024,1025,1026,1027,1028,1029,1030,1031,1032,1033,1034,1035,1036, &
   & 1037,1038,1039,1040,1041,1043,1044,1045,1046,1047,1049,1050,1051,1052,1053,1054,1055,1056,1057,1058,1059,1060,1061,1062, &
   & 1063,1064,1065,1066,1068,1070,1071,1072,1073,1074,1075,1076,1077,1078,1079,1080,1081,1082,1083,1084,1085,1086,1087,1088, &
   & 1089,1090,1091,1092,1093,1094,1095,1096,1097,1098,1099,1100,1101,1102,1103,1104,1105,1106,1107,1108,1109,1110,1111,1112, &
   & 1113,1114,1115,1116,1117,1118,1119,1120,1121,1122,1123,1124,1125,1126,1127,1128,1129,1130,1131,1132,1133,1134,1135,1136,1137,&
   & 1138,1139,1140,1141,1142,1143,1144,1145,1146,1147,1148,1149,1150,1151,1152,1153,1154,1155,1156,1157,1158,1159,1160,1161,1162,&
   & 1165,1166,1167,1168,1169,1170,1171,1173,1174,1175,1176,1177,1178,1179,1181,1182,1188,1189,1190,1191,1192,1193,1194,1195,1196,&
   & 1197,1198,1199,1200,1201,1202,1203,1204,1205,1206,1207,1208,1209,1210,1211,1212,1213,1214,1215,1216,1217,1218,1219,1220,1221,&
   & 1222,1223,1224,1225,1226,1227,1228,1229,1230,1231,1232,1233,1234,1235,1236,1237,1238,1239,1240,1241,1242,1243,1244,1245,1246,&
   & 1247,1248,1249,1250,1251,1252,1253,1254,1255,1256,1257,1258,1259,1260,1261,1262,1263,1264,1265,1266,1267,1268,1269,1270,1271,&
   & 1272,1273,1274,1275,1276,1277,1278,1279,1280,1281,1282,1283,1284,1285,1288,1289,1290,1291,1292,1293,1294,1295,1297,1301,1303,&
   & 1304,1305,1307,1308,1309,1311,1312,1313,1317,1318,1319,1331,1333,1334,1335,1336,1337,1338,1339,1340,1341,1342,1343,1344,1345,&
   & 1346,1347,1348,1349,1350,1351,1352,1353,1354,1355,1356,1357,1358,1359,1360,1361,1364,1365,1366,1367,1370,1371,1372,1373,1376,&
   & 1377,1378,1379,1380,1381,1382,1383,1384,1385,1386,1387,1388,1389,1390,1391,1392,1393,1394,1395,1396,1397,1398,1399]

!   & [6,8,9,11,13,15,17,19,21,22,23,24,25,26,27,28,31,32,33,34,35,36,37,38,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,70,71, &
!   & 72,73,74,75,76,77,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,&
!   & 128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158, &
!   & 159,160,161,162,163,164,165,166,167,168,176,169,170,171,172,173,174,175,177,178,179,180,181,182,183,184,185,186,187,188,189, &
!   & 190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220, &
!   & 221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251, &
!   & 252,253]
  
  ! gaussian distribution with mean 0 and std 2
  INTEGER(CMISSIntg), DIMENSION(2700), PARAMETER :: IZ_offset= &
   & [-1,-1,-2,0,-3,-1,2,2,-4,1,-1,1,0,-1,0,1,1,-1,-3,-1, & 
   & 2,1,1,-1,-2,4,0,0,-1,0,-1,-2,1,3,0,-1,0,2,4,-3, & 
   & 0,-1,1,0,-1,-4,-2,-1,-1,2,-1,1,2,0,0,-1,2,3,1,1, & 
   & 0,2,-1,-1,1,-1,2,-3,3,-1,0,-1,1,-3,-2,3,-2,-1,0,-2, & 
   & 0,2,0,3,-3,-1,-2,1,-3,-4,-1,-4,-2,-1,-3,1,1,-1,-3,-1, & 
   & 0,-4,3,1,-4,4,0,-3,-2,1,-3,-1,0,-1,0,-1,-1,0,-4,-4, & 
   & -2,0,0,-1,3,2,0,7,-1,2,0,-3,4,-2,0,4,1,1,0,0, & 
   & 4,-1,0,-1,-1,0,1,0,-2,3,-2,-1,0,-3,-2,-3,-5,4,-1,2, & 
   & 0,-2,0,1,-2,2,2,0,1,-2,1,0,3,-2,4,0,1,3,2,-2, & 
   & -2,-2,2,0,-1,0,2,0,1,-5,-1,-2,-1,-2,-1,-1,1,3,-3,4, & 
   & 3,0,3,0,-2,2,4,-2,-2,2,0,-3,2,0,-2,-1,-1,2,-2,1, & 
   & 2,-2,1,-3,1,-2,2,1,-1,-3,3,0,3,0,-4,3,-1,2,-1,1, & 
   & 2,2,1,0,0,1,-1,1,-2,-4,0,0,0,4,2,4,0,2,-1,0, & 
   & -3,0,2,-2,0,-1,-2,2,0,1,1,2,-5,-2,0,1,2,-1,-1,-1, & 
   & 1,0,2,-2,1,1,-3,1,3,-2,-1,-2,2,-1,2,-3,-4,2,1,1, & 
   & 2,0,0,-2,0,-2,-2,-1,0,3,1,1,1,2,1,0,-1,1,-3,0, & 
   & 2,0,1,0,1,-1,3,1,0,-2,0,-3,-1,-2,0,0,-1,3,3,3, & 
   & 1,2,-1,0,-2,0,2,-2,0,-3,1,2,2,-2,-1,0,-1,-1,-5,1, & 
   & 0,1,1,-6,-1,2,2,-1,1,0,2,4,-2,4,0,-2,-2,0,0,1, & 
   & 2,1,-2,-1,-2,-1,1,1,0,-2,-1,1,0,3,-2,1,2,-2,0,-2, & 
   & -2,-3,-2,1,0,1,2,1,2,2,5,2,-3,1,1,2,2,-1,0,1, & 
   & 0,4,-3,2,2,-1,1,0,-2,0,3,-1,1,-2,-2,-2,-3,-2,-2,-1, & 
   & -3,-1,-2,-2,-1,2,0,1,3,1,0,2,4,-1,2,-2,-2,-1,0,-4, & 
   & -3,0,1,2,1,0,1,5,4,0,3,-1,1,0,2,-2,1,0,1,-2, & 
   & -1,-3,-1,-2,-1,0,2,-2,1,2,2,-3,1,0,3,2,-1,1,1,-5, & 
   & -1,-1,-2,-4,2,-1,2,-2,2,1,1,-2,0,-1,-2,-1,-3,0,0,-2, & 
   & 1,-2,0,1,0,-5,-2,-1,-1,-3,0,0,3,0,1,0,-4,-3,-2,-1, & 
   & 1,-2,-1,1,-1,-2,-1,-1,1,0,1,-1,4,0,-2,-3,0,-1,-3,-1, & 
   & -3,1,4,1,-1,3,1,-4,0,0,1,0,2,-3,5,3,0,1,-2,-2, & 
   & 0,2,2,4,-1,0,-5,-1,0,0,-1,-3,-1,4,-2,0,-2,2,2,0, & 
   & 0,3,1,0,0,3,1,1,1,0,-2,1,-3,-1,2,-1,1,0,1,1, & 
   & -1,-4,2,-1,1,2,-2,3,1,0,-1,2,1,1,5,1,2,0,4,0, & 
   & -2,0,-3,3,-2,1,2,1,1,-2,0,2,-2,-2,-1,2,-1,0,-3,-2, & 
   & -1,-3,2,2,-3,0,-2,0,0,3,1,-1,2,2,0,4,0,-1,1,-3, & 
   & -1,-1,3,1,1,2,2,0,-3,-3,0,3,2,2,4,-1,-1,2,-1,4, & 
   & 2,-1,-1,1,1,2,1,-2,0,-3,1,1,0,-2,3,2,-1,-3,2,1, & 
   & 0,1,1,-1,3,1,-2,-1,3,-1,-2,2,1,4,-1,1,0,1,-3,0, & 
   & 2,0,-1,0,-4,2,1,1,1,1,3,-2,0,2,0,1,-3,1,2,-4, & 
   & 0,-3,-2,1,2,0,-1,4,0,1,-1,1,-2,0,-2,0,1,4,-1,1, & 
   & -1,0,6,0,-1,2,2,-3,-1,-1,-1,0,-1,1,1,-2,1,-2,2,0, & 
   & 1,-2,2,-1,0,-1,1,4,2,0,-1,0,-2,0,2,0,2,-1,1,1, & 
   & -3,-1,2,-2,0,-2,-2,-2,2,0,1,-1,0,2,-3,-2,2,3,-2,1, & 
   & -2,1,0,-2,1,3,0,-2,-2,-2,1,1,-2,2,-2,3,-4,-1,2,0, & 
   & 0,-1,4,-1,1,2,-1,3,1,-3,-1,0,1,1,1,-3,-3,3,2,1, & 
   & 3,4,-2,0,1,0,0,0,2,-5,0,-1,-1,0,1,-1,0,-3,-2,-1, & 
   & 0,-2,-2,1,0,3,1,-2,-4,0,0,1,2,0,0,-1,2,0,-1,-1, & 
   & 4,1,1,1,1,1,0,3,0,-2,-1,2,2,2,0,0,-4,2,0,0, & 
   & 2,-1,1,-1,-1,1,0,-3,-1,6,1,0,-1,-4,0,-6,-2,-3,0,2, & 
   & 1,-1,-3,-1,-2,1,-1,3,0,-3,0,1,-1,1,1,1,-1,0,0,0, & 
   & -3,0,-2,0,-2,3,2,-1,4,2,2,-4,-3,-1,-1,-3,1,-1,0,-1, & 
   & 0,3,-4,-2,-1,1,-2,-1,1,-2,2,1,1,0,-2,-3,0,1,1,0, & 
   & 3,0,-1,0,-1,1,3,2,0,-3,-1,1,0,3,1,-3,1,1,4,0, & 
   & 1,-3,-1,-3,0,-2,-3,0,-2,3,4,0,1,-4,3,-1,1,0,1,0, & 
   & 0,-2,0,-3,0,0,-1,1,2,-3,-3,-1,1,1,-1,1,4,-2,-2,-4, & 
   & 0,2,-1,2,0,0,2,1,1,3,-2,-1,1,-1,-1,1,-1,1,0,-1, & 
   & 1,1,-4,-1,-2,0,-1,0,3,3,2,-3,3,2,1,5,-1,4,-2,-1, & 
   & -2,1,1,-1,-3,-5,-2,-3,1,2,-4,-3,-1,2,3,3,-1,-1,-2,1, & 
   & 2,-3,0,-1,1,2,2,1,-1,-1,-2,2,-2,-4,-2,3,1,-1,0,1, & 
   & 1,0,5,2,-5,-2,0,1,-1,0,-1,-3,4,0,5,-2,-1,-3,-2,2, & 
   & 2,-2,0,-1,2,-2,3,-3,1,-3,-1,-2,3,0,1,-3,-2,1,0,1, & 
   & 2,-2,-2,1,0,0,-2,-2,1,3,-1,2,3,2,3,1,-3,-2,0,-3, & 
   & -1,3,-2,0,1,4,-1,1,-2,-2,1,1,-1,1,-1,1,-1,1,1,3, & 
   & -2,0,-3,0,2,-2,0,2,3,0,3,0,-2,0,-3,-2,-1,3,0,4, & 
   & 1,1,0,-2,0,-1,-3,1,0,1,-1,1,-1,1,1,0,-2,0,2,-4, & 
   & 0,0,1,-1,1,-1,1,2,0,-2,0,1,-1,3,-1,-1,3,1,1,3, & 
   & -1,0,2,-2,0,-2,1,4,0,-1,3,-2,-1,0,1,-2,-2,1,0,-3, & 
   & -3,-3,0,1,1,-2,-4,1,0,0,-1,-1,4,0,1,-4,-3,-2,0,1, & 
   & 5,-4,-3,-1,-2,-2,2,-1,1,0,-3,2,-1,-1,0,1,2,0,-3,2, & 
   & -1,4,1,1,-1,-2,-1,-3,-2,2,-1,0,3,-1,0,1,0,1,-1,0, & 
   & 3,4,0,-2,0,-1,-1,-4,-3,0,-1,3,-4,-2,1,1,-1,-2,-5,-1, & 
   & 1,0,-3,2,1,0,0,0,0,2,-1,3,-2,-1,-2,-1,-3,2,0,0, & 
   & 0,-1,-1,2,0,1,1,2,-1,-1,1,-2,2,-2,2,3,-3,2,1,0, & 
   & 0,0,2,-1,-1,6,3,2,2,0,-3,0,2,1,3,2,1,-1,-1,1, & 
   & -1,-1,3,0,-3,-4,0,-1,-1,0,0,-2,0,-4,-2,-1,1,1,3,3, & 
   & 1,-2,1,3,2,0,-3,-3,-1,1,-1,-3,3,-1,2,1,-3,0,1,-1, & 
   & -1,1,-2,-1,3,-1,1,-1,-3,-1,-1,0,-1,2,-4,1,-1,2,0,3, & 
   & -1,-1,1,2,-2,-5,0,0,3,-1,0,1,-1,-2,0,-4,1,-2,-3,0, & 
   & 1,0,1,2,1,0,1,-1,3,-1,-2,-1,-1,-3,-3,3,0,-2,2,-2, & 
   & 3,-1,-2,3,-1,0,1,1,-1,-4,-2,1,2,0,2,1,3,3,-1,0, & 
   & 1,3,0,-1,0,-3,-1,-2,0,0,2,0,-2,1,-2,0,5,-1,0,2, & 
   & -2,3,-1,2,-1,3,2,-1,4,-1,2,0,0,4,-3,-3,2,-4,2,3, & 
   & -3,0,-1,2,-1,2,-2,2,1,-2,-2,0,-4,2,3,1,-2,0,4,3, & 
   & 1,2,3,0,0,-1,0,0,-2,-2,1,-1,-2,-3,-3,1,-1,1,2,-1, & 
   & 1,6,-1,1,-1,2,-1,-1,1,0,-1,-4,-3,-1,2,1,1,2,0,1, & 
   & 0,1,0,2,-3,5,-1,-2,-2,0,3,-2,1,-2,2,-1,-2,4,1,3, & 
   & 2,0,3,0,1,-3,-2,2,-1,-1,1,-5,-1,-2,1,0,-3,0,0,1, & 
   & -1,-3,0,1,1,1,-2,-2,-1,1,1,1,2,1,0,0,4,3,-1,-4, & 
   & -1,0,-4,3,3,0,-1,2,1,0,-1,1,4,-3,3,-4,-2,1,2,-2, & 
   & 0,-2,-3,2,-1,0,4,1,0,1,0,-2,1,0,0,2,0,-1,-2,1, & 
   & 2,0,0,0,0,0,3,2,-1,1,-1,-1,0,0,-1,-3,2,-1,0,1, & 
   & 0,0,1,1,-1,1,1,1,-3,2,-2,4,-1,0,2,-3,1,1,1,-1, & 
   & 1,-2,1,-1,0,-2,1,-1,-2,-4,1,3,-1,3,0,-2,0,-1,1,-2, & 
   & 0,1,0,1,-1,-3,1,-2,2,-1,2,3,-1,-1,2,1,0,0,-1,2, & 
   & -1,2,-1,3,1,-6,0,3,-1,-1,1,3,1,0,1,1,0,1,0,1, & 
   & 1,-3,-1,-2,0,-5,0,2,1,3,0,-1,-4,5,1,2,1,-1,-1,-5, & 
   & -2,-1,-1,1,1,-1,-2,2,0,1,1,0,1,0,0,-2,0,2,-1,1, & 
   & 1,7,-1,0,-1,-3,0,1,-1,0,0,-1,-2,-1,-2,1,0,1,-1,-4, & 
   & 0,1,2,-4,1,0,0,0,4,-2,0,1,-4,1,3,-1,-1,-3,1,2, & 
   & -3,-3,-1,0,2,1,0,-2,0,1,-1,0,1,3,3,-4,2,-1,-3,-1, & 
   & -1,1,0,-3,-1,-4,3,-1,1,1,-1,-3,-3,2,0,2,-1,-6,-2,2, & 
   & 0,4,1,-4,5,2,-3,0,0,-2,-4,-2,1,-1,-2,0,5,0,2,-1, & 
   & -3,-2,1,0,0,2,4,2,-1,1,2,0,3,2,1,-3,-2,-2,2,0, & 
   & 4,0,0,0,3,2,0,-1,-1,-1,1,2,3,2,-2,2,0,-1,0,0, & 
   & 0,0,0,-1,2,3,1,0,3,-1,-1,0,2,0,1,0,0,-4,3,2, & 
   & -2,-5,3,0,2,1,1,2,1,-3,1,-1,0,0,1,-1,-1,-2,2,0, & 
   & 1,2,-2,-2,0,1,-2,-2,1,1,1,1,2,-3,0,1,-3,0,2,-3, & 
   & 0,-1,1,-2,-2,3,1,3,2,1,0,-3,0,0,1,-2,0,3,-1,0, & 
   & 0,-2,1,0,-5,4,3,2,2,-1,-2,-1,1,4,-2,4,-2,2,0,-4, & 
   & 5,-2,-1,-1,-1,0,1,4,3,2,-1,2,0,-1,1,-2,-5,-1,-2,2, & 
   & -1,3,-1,1,1,1,-2,0,0,-1,-1,-1,-1,-2,-2,-1,-1,2,-1,3, & 
   & -3,4,0,-4,1,-4,1,5,0,1,1,-1,2,-1,2,-2,0,-1,-2,0, & 
   & -3,-3,-3,-2,2,0,2,2,-1,-3,0,-1,2,-1,3,1,1,-1,-2,0, & 
   & 0,-3,0,-2,1,2,-4,-1,-1,-5,-5,-3,2,3,-1,-1,0,2,2,5, & 
   & 1,-2,5,-2,-1,2,0,0,0,-3,-3,-3,2,-2,-3,1,-1,-3,1,1, & 
   & 1,-1,-1,-1,0,1,3,-1,-1,0,2,-1,2,-1,-2,4,2,1,0,-1, & 
   & 2,2,-4,-2,-1,3,1,3,-1,1,-1,-4,1,0,-2,-1,-2,1,-2,1, & 
   & -1,0,-2,-2,-1,-2,-1,-1,-3,-3,-6,1,-2,0,-1,2,1,3,0,1, & 
   & -1,-2,-3,-3,4,-3,1,0,1,0,-1,1,-1,0,3,1,-1,-1,3,-4, & 
   & 1,3,0,0,2,5,2,-1,0,-1,0,-5,-2,1,1,-3,1,0,2,0, & 
   & -2,-2,-2,0,2,-1,-1,-1,-4,-3,-1,1,-2,-2,-2,2,1,0,2,2, & 
   & -1,-3,-3,2,1,-3,-1,0,0,0,-2,2,-1,-2,-4,-3,6,-1,3,-1, & 
   & 2,-1,-1,1,0,1,-1,1,-1,2,1,1,-2,-3,-1,0,-3,1,3,-1, & 
   & 0,-2,0,-2,2,2,0,-1,3,-4,1,1,0,3,-2,0,3,-2,3,-4, & 
   & 0,3,1,0,1,-2,-3,-4,1,-4,1,-4,1,1,0,-1,3,2,-2,-1, & 
   & 1,-1,2,-3,1,1,1,0,-3,3,0,-3,-2,-3,-2,1,-1,0,2,-1, & 
   & 2,0,-1,1,0,-1,3,0,0,-1,-3,3,-2,0,0,2,-1,-2,-1,0, & 
   & 0,-3,2,-1,5,0,-1,0,0,1,5,-4,-2,-1,1,-1,1,-3,2,1, & 
   & 1,-2,-1,0,-2,2,-1,-2,1,-2,2,-3,-1,4,-2,-1,1,1,6,3, & 
   & -1,1,3,0,-5,-1,-3,1,0,2,-2,2,-2,3,0,2,-3,2,-1,-3, & 
   & 3,1,0,0,-1,-3,-4,-3,-1,1,0,0,-1,1,2,-1,-3,1,-1,-4, & 
   & 0,1,1,1,0,0,4,0,-1,3,2,3,2,1,3,0,-1,-2,3,5, & 
   & 2,-1,-2,-2,0,-3,-3,-3,1,-4,1,1,-2,-2,-1,-1,2,0,-2,1, & 
   & 2,-3,2,-2,-4,0,2,2,-1,1,0,0,2,0,2,1,-1,1,0,-2, & 
   & 0,0,-1,0,-2,4,-3,-3,2,-2,-1,-1,2,0,1,-4,-4,0,-4,4, & 
   & 1,3,-4,1,1,2,-1,-3,-1,-1,-1,0,1,4,-2,-1,4,5,-1,-4]
   

  INTEGER(CMISSIntg) :: FT_1,FT_2,FT_3,FT_4,FT_5,FT_6,FT_7,FT_8,FT_9,FT_10

    
  !the nodes in each element for quadratic basis functions
  INTEGER(CMISSIntg), DIMENSION(NumberOfElementsFE,27) :: AllQuadraticElements
  !the nodes in each element for linear basis functions
!  INTEGER(CMISSIntg), DIMENSION(NumberOfElementsFE,8) :: AllLinearElements

  !the entries in MyElemen1,... to be used in linearElement1,...
  INTEGER(CMISSIntg), DIMENSION(8), PARAMETER :: ENTRIES = [1,3,7,9,19,21,25,27]


  real(CMISSRP) :: posX,posY,posZ

  ! this should be 1 !!!!!
  real(CMISSRP), parameter :: scalefactor = 10.0_CMISSRP

  integer(CMISSIntg) :: stat
  character(len=256) :: pathname,filename,filename2

  character(len=8), parameter :: element_str = "Element:"
  character(len=6), parameter :: nodes_str = "Nodes:"
  character(len=4), parameter :: node_str = "Node"
  character(len=5), parameter :: nodep_str = "Node:"
  character(len=132) :: blub,blub2,blub3,blub4,blub5,blub6

  INTEGER(CMISSIntg) :: mu_nr,Ftype,fibre_nr,NearestGP,InElement

  real(CMISSRP), dimension(1000,10) :: FIRING_TIMES
  real(CMISSRP) :: time
  
  logical :: less_info,fast_twitch

  !--------------------------------------------------------------------------------------------------------------------------------

  !all times in [ms]
  REAL(CMISSRP), PARAMETER :: STIM_STOP=0.1_CMISSRP
!  REAL(CMISSRP), PARAMETER :: PERIODD=5.00_CMISSRP
  REAL(CMISSRP), PARAMETER :: PERIODD=1.00_CMISSRP
  REAL(CMISSRP), PARAMETER :: TIME_STOP=1000.0_CMISSRP !10.0_CMISSRP !30.0_CMISSRP
  
  REAL(CMISSRP), PARAMETER :: ODE_TIME_STEP=0.0001_CMISSRP
  REAL(CMISSRP), PARAMETER :: PDE_TIME_STEP=0.0005_CMISSRP
  REAL(CMISSRP), PARAMETER :: ELASTICITY_TIME_STEP=0.10000000001_CMISSRP!0.5_CMISSRP!0.05_CMISSRP!0.8_CMISSRP
!tomo keep ELASTICITY_TIME_STEP and STIM_STOP at the same values

  INTEGER(CMISSIntg), PARAMETER :: OUTPUT_FREQUENCY=2
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !stimulation current in [uA/cm^2]
  REAL(CMISSRP) :: STIM_VALUE,STIM_VALUE2

  !condctivity in [mS/cm]
  REAL(CMISSRP), PARAMETER :: CONDUCTIVITY=3.828_CMISSRP

  !surface area to volume ratio in [cm^-1]
  REAL(CMISSRP), PARAMETER :: Am=500.0_CMISSRP

  !membrane capacitance in [uF/cm^2]
  REAL(CMISSRP), PARAMETER :: Cm_fast=1.0_CMISSRP 
  REAL(CMISSRP), PARAMETER :: Cm_slow=0.58_CMISSRP

  !--------------------------------------------------------------------------------------------------------------------------------
  !maximum contraction velocity in [cm/ms]
  REAL(CMISSRP), PARAMETER :: Vmax=-0.02_CMISSRP ! 0.02 cm/ms = 0.2 m/s, rat GM    ! this is how it should be!!!
!  REAL(CMISSRP), PARAMETER :: Vmax=-0.2_CMISSRP !this value is 10 times too large, which stabilizes the problem !TODO
  !all lengths in [cm]
  REAL(CMISSRP), PARAMETER :: LENGTH= 2.5_CMISSRP ! X-direction
!  REAL(CMISSRP), PARAMETER :: WIDTH= 0.05_CMISSRP ! Y-direction
!  REAL(CMISSRP), PARAMETER :: HEIGHT=0.05_CMISSRP ! Z-direction

  !--------------------------------------------------------------------------------------------------------------------------------
  !CAUTION - what are the units???
!  REAL(CMISSRP), PARAMETER, DIMENSION(4) :: MAT_FE= &
!    & [0.0000000000635201_CMISSRP,0.3626712895523322_CMISSRP,0.0000027562837093_CMISSRP,43.372873938671383_CMISSRP] ![N/cm^2]

  !Material-Parameters C=[mu_1, mu_2, mu_3, alpha_1, alpha_2, alpha_3, mu_0, XB_stiffness]
  REAL(CMISSRP), PARAMETER, DIMENSION(8) :: C = &
    & [0.0085_CMISSRP*5.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP, & ! [N/cm^2 = 10^4 J/m^3]
    &  11.0_CMISSRP,1.0_CMISSRP,6.0_CMISSRP, &
    &  1.0_CMISSRP,2.2e-9_CMISSRP]      

  REAL(CMISSRP) :: INIT_PRESSURE

  !--------------------------------------------------------------------------------------------------------------------------------
  
  INTEGER(CMISSIntg), PARAMETER, DIMENSION(21) :: TOP_NODES = &
   & [66,67,127,206,213,217,221,254,255,262,263,290,291,297,299,340,341,344,345,347,363]
  INTEGER(CMISSIntg), PARAMETER, DIMENSION(21) :: BOTTOM_NODES = &
   & [57,58,107,229,230,231,232,232,234,281,282,283,284,326,327,328,329,330,331,360,361]
  INTEGER(CMISSIntg), PARAMETER, DIMENSION(52) :: SIDE_NODES = &
   & [484,485,500,220,222,261,486,487,501,214,251,265,488,489,502,216,218,267,490,491,503, &
   &  209,240,269,492,493,504,212,219,276,494,495,505,224,228,278,496,497,506,225,227,280, &
   &  538,540,542,233,498,499,507,539,541,543]


  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumberFE=1
  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumberM=2
  INTEGER(CMISSIntg), PARAMETER :: NumberOfSpatialCoordinates=3

  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumberFE=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumberM=2

  INTEGER(CMISSIntg), PARAMETER :: QuadraticBasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: LinearBasisUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: LinearBasisUserNumberM=3
  INTEGER(CMISSIntg), PARAMETER :: NumberOfXiCoordinates=NumberOfSpatialCoordinates
  INTEGER(CMISSIntg), PARAMETER :: NumberOfGaussPoints=3

  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshDimensionsFE=3

  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumberFE=1
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumberM=2
  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshComponentsFE=2
  INTEGER(CMISSIntg), PARAMETER :: QuadraticMeshComponentNumber=1
  INTEGER(CMISSIntg), PARAMETER :: LinearMeshComponentNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MonodomainMeshComponentNumber=1

  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumberFE=1
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumberM=2

  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumberFE=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumberM=2

  INTEGER(CMISSIntg), PARAMETER :: FieldFibreUserNumber=3

  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumberM=4
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumberFE=5

  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumberM=6
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumberFE=7

  INTEGER(CMISSIntg), PARAMETER :: FieldIndependentUserNumberFE=8
  INTEGER(CMISSIntg), PARAMETER :: FieldIndependentUserNumberM=9

  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumberM=10
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumberFE=11

  INTEGER(CMISSIntg), PARAMETER :: CellMLUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: CellMLModelsFieldUserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: CellMLStateFieldUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: CellMLIntermediateFieldUserNumber=14
  INTEGER(CMISSIntg), PARAMETER :: CellMLParametersFieldUserNumber=15

  INTEGER(CMISSIntg), PARAMETER :: EquationsSetsUserNumberM=1
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetsUserNumberFE=2

  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

  INTEGER(CMISSIntg), PARAMETER :: SolverDAEIndex=1
  INTEGER(CMISSIntg), PARAMETER :: SolverParabolicIndex=2
  INTEGER(CMISSIntg), PARAMETER :: SolverFEIndex=1

  INTEGER(CMISSIntg), PARAMETER :: ControlLoopMonodomainNumber=1
  INTEGER(CMISSIntg), PARAMETER :: ControlLoopElasticityNumber=2
  
  !Program types
  
  !Program variables

  INTEGER(CMISSIntg) :: EquationsSetIndexM,EquationsSetIndexFE
  INTEGER(CMISSIntg) :: CellMLIndex
  INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  
  INTEGER(CMISSIntg) :: NodeNumber,NodeDomain,ElementDomain,node_idx,NumberOfDomains,domain_idx!,ComponentNumber
  INTEGER(CMISSIntg) :: NumberOfNodesInXi1,NumberOfNodesInXi2,NumberOfNodesInXi3,NumberOfNodesInXi1skin
  INTEGER(CMISSIntg) :: i,j,k,m,my_node_idx,elem_idx,node1,node2,elem_idx2,val

  LOGICAL :: EXPORT_FIELD

  INTEGER(CMISSIntg) :: shortenModelIndex,shortenModelIndex2
  INTEGER(CMISSIntg) :: stimcomponent

  REAL(CMISSRP) :: VALUE

  INTEGER(CMISSIntg) :: Err


  !CMISS variables

  TYPE(cmfe_BasisType) :: QuadraticBasis,LinearBasis,LinearBasisM
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditionsM,BoundaryConditionsFE
  TYPE(cmfe_CellMLType) :: CellML
  TYPE(cmfe_CellMLEquationsType) :: CellMLEquations
  TYPE(cmfe_ControlLoopType) :: ControlLoopMain
  TYPE(cmfe_ControlLoopType) :: ControlLoopM,ControlLoopFE
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystemFE,CoordinateSystemM,WorldCoordinateSystem
  TYPE(cmfe_DecompositionType) :: DecompositionFE,DecompositionM
  TYPE(cmfe_EquationsType) :: EquationsM,EquationsFE
  TYPE(cmfe_EquationsSetType) :: EquationsSetM,EquationsSetFE
  TYPE(cmfe_FieldType) :: EquationsSetFieldM,EquationsSetFieldFE
  TYPE(cmfe_FieldType) :: GeometricFieldM,GeometricFieldFE
  TYPE(cmfe_FieldType) :: DependentFieldM,DependentFieldFE
  TYPE(cmfe_FieldType) :: IndependentFieldFE,IndependentFieldM
  TYPE(cmfe_FieldType) :: MaterialFieldM,MaterialFieldFE
  TYPE(cmfe_FieldType) :: FibreField
  TYPE(cmfe_FieldType) :: CellMLModelsField,CellMLStateField,CellMLIntermediateField,CellMLParametersField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_MeshType) :: MeshFE,MeshM
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_RegionType) :: RegionFE,RegionM,WorldRegion
  TYPE(cmfe_SolverType) :: SolverDAE,SolverParabolic
  TYPE(cmfe_SolverType) :: SolverFE,LinearSolverFE
  TYPE(cmfe_SolverEquationsType) :: SolverEquationsM,SolverEquationsFE
  TYPE(cmfe_NodesType) :: Nodes
  TYPE(cmfe_MeshElementsType) :: QuadraticElements
  TYPE(cmfe_MeshElementsType) :: LinearElements
  TYPE(cmfe_MeshElementsType) :: ElementsM


#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables

#ifdef WIN32
  !Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  !Set the window parameters
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  !If attempt fails set with system estimated values
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif


!##################################################################################################################################
  less_info = .true.!.false.!
  if(less_info) then
    NumberOfNodesInXi1=50!500!240
    NumberOfNodesInXi1skin=10
    NumberOfNodesInXi2=2!15
    NumberOfNodesInXi3=2!15
  else
    NumberOfNodesInXi1=50!500
    NumberOfNodesInXi2=4
    NumberOfNodesInXi3=4
  endif
!##################################################################################################################################

!  NumberOfNodesPerFibre=(NumberOfNodesInXi1-1)*NumberGlobalXElements/NumberOfInSeriesFibres+1
!  NumberOfNodesM=NumberOfNodesPerFibre*NumberOfInSeriesFibres*NumberGlobalYElements*NumberGlobalZElements*NumberOfNodesInXi2* &
!    & NumberOfNodesInXi3
!  NumberOfElementsM=(NumberOfNodesPerFibre-1)*NumberOfInSeriesFibres*NumberGlobalYElements*NumberGlobalZElements* &
!    & NumberOfNodesInXi2*NumberOfNodesInXi3
  NumberOfNodesPerFibre=NumberOfNodesInXi1
  NumberOfElementsPerElasticityElement=(NumberOfNodesInXi1-1)*NumberOfNodesInXi2*NumberOfNodesInXi3
!  NumberOfNodesM=NumberOfNodesInXi1*NumberOfNodesInXi2*NumberOfNodesInXi3*NumberOfElementsFE
!  NumberOfElementsM=(NumberOfNodesInXi1-1)*NumberOfNodesInXi2*NumberOfNodesInXi3*NumberOfElementsFE
  NumberOfNodesMmuscle=NumberOfNodesInXi1*NumberOfNodesInXi2*NumberOfNodesInXi3*12
  NumberOfNodesMskin=NumberOfNodesInXi1skin*NumberOfNodesInXi2*NumberOfNodesInXi3*4
  NumberOfNodesM=NumberOfNodesMmuscle+NumberOfNodesMskin

  NumberOfElementsMmuscle=(NumberOfNodesInXi1-1)*NumberOfNodesInXi2*NumberOfNodesInXi3*12
  NumberOfElementsMskin=(NumberOfNodesInXi1skin-1)*NumberOfNodesInXi2*NumberOfNodesInXi3*4
  NumberOfElementsM=NumberOfElementsMmuscle+NumberOfElementsMskin

!##################################################################################################################################
  !pathname="/home/heidlauf/OpenCMISS/OpenCMISS_cmake_new/examples/TA_w_skin/Fortran/src/"
  pathname="input/"
!!!  fast_twitch=.true.
!!!  if(fast_twitch) then
!!!    filename=trim(pathname)//"fast_2014_03_25_no_Fl_no_Fv.xml" !FAST currently valid (new) CellML file
!!!!    filename2= &
!!!!   &"/home/heidlauf/OpenCMISS/opencmiss/examples/MultiPhysics/BioelectricFiniteElasticity/cellModelFiles/fast_stim_2012_07_23.xml"
!!!!   &"/home/heidlauf/OpenCMISS/opencmiss/examples/MultiPhysics/BioelectricFiniteElasticity/cellModelFiles/shorten_mod_2011_07_04.xml"
!!!    STIM_VALUE=1200.0_CMISSRP
!!!!    STIM_VALUE=8000.0_CMISSRP
!!!  else !slow twitch
!!!    filename=trim(pathname)//"slow_2014_01_10.xml" !SLOW currently valid (new) CellML file
!!!!  &"/home/heidlauf/OpenCMISS/opencmiss/examples/MultiPhysics/BioelectricFiniteElasticity/cellModelFiles/slow_2012_07_23.xml_0.401"
!!!!  &"/home/heidlauf/OpenCMISS/opencmiss/examples/MultiPhysics/BioelectricFiniteElasticity/cellModelFiles/slow_twitch_2012_01_27.xml"
!!!    STIM_VALUE=2000.0_CMISSRP
!!!  endif

  IF(biophysical) THEN
    !biophysical half-sarcomere model
    filename=trim(pathname)//"slow_2014_10_01_no_Fl_no_Fv.xml" !SLOW
    filename2=trim(pathname)//"fast_2014_03_25_no_Fl_no_Fv.xml" !FAST
    STIM_VALUE=2000.0_CMISSRP
    STIM_VALUE2=1200.0_CMISSRP
  ELSE
    !phenomenological half-sarcomere model
    filename=trim(pathname)//"Aliev_Panfilov_Razumova_2016_08_22.cellml"
    STIM_VALUE=90.0_CMISSRP
    STIM_VALUE2=90.0_CMISSRP
  END IF


!##################################################################################################################################




  !================================================================================================================================
  !  G E N E R A L   F E A T U R E S
  !================================================================================================================================

  !--------------------------------------------------------------------------------------------------------------------------------
  !Intialise OpenCMISS
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)

  !Trap errors
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)
  
!  CALL cmfe_DiagnosticsSetOn(CMFE_FROM_DIAG_TYPE,[1,2,3,4,5],"dignostics",["FIELD_MAPPINGS_CALCULATE"],err)
  
  !Get the computational nodes information
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  CALL cmfe_OutputSetOn("EntireTA",Err)
  
  NumberOfDomains=NumberOfComputationalNodes
  
  !Broadcast the number of elements in the X & Y directions and the number of partitions to the other computational nodes
!  CALL MPI_BCAST(NumberGlobalXElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
!  CALL MPI_BCAST(NumberGlobalYElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
!  CALL MPI_BCAST(NumberGlobalZElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NumberOfDomains,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystemFE,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumberFE,CoordinateSystemFE,Err)
  !Set the coordinate system to be 3D
  CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystemFE,3,Err)
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystemFE,Err)


  ! CREATE A 1D COORDINATE SYSTEM
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystemM,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumberM,CoordinateSystemM,Err)
  !Set the coordinate system to be 1D
  CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystemM,1,Err)
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystemM,Err)


  !--------------------------------------------------------------------------------------------------------------------------------
  !Start the creation of the region
  CALL cmfe_Region_Initialise(RegionFE,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumberFE,WorldRegion,RegionFE,Err)
  CALL cmfe_Region_CoordinateSystemSet(RegionFE,CoordinateSystemFE,Err)
  CALL cmfe_Region_LabelSet(RegionFE,"Region3D",Err)
  CALL cmfe_Region_CreateFinish(RegionFE,Err)


  ! CREATE A SECOND REGION
  CALL cmfe_Region_Initialise(RegionM,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumberM,WorldRegion,RegionM,Err)
  CALL cmfe_Region_CoordinateSystemSet(RegionM,CoordinateSystemM,Err)
!  CALL cmfe_Region_CoordinateSystemSet(RegionM,CoordinateSystemFE,Err)
  CALL cmfe_Region_LabelSet(RegionM,"Region1D",Err)
  CALL cmfe_Region_CreateFinish(RegionM,Err)


  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the bases
  !Define basis functions - tri-Quadratic Lagrange 
  CALL cmfe_Basis_Initialise(QuadraticBasis,Err)
  CALL cmfe_Basis_CreateStart(QuadraticBasisUserNumber,QuadraticBasis,Err)
  CALL cmfe_Basis_TypeSet(QuadraticBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(QuadraticBasis,NumberOfXiCoordinates,Err)
  CALL cmfe_Basis_InterpolationXiSet(QuadraticBasis,[CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
   & CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(QuadraticBasis, &
   & [NumberOfGaussPoints,NumberOfGaussPoints,NumberOfGaussPoints],Err)
  CALL cmfe_Basis_CreateFinish(QuadraticBasis,Err)

  !Define basis functions - tri-Linear Lagrange
  CALL cmfe_Basis_Initialise(LinearBasis,Err)
  CALL cmfe_Basis_CreateStart(LinearBasisUserNumber,LinearBasis,Err)
  CALL cmfe_Basis_TypeSet(LinearBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(LinearBasis,NumberOfXiCoordinates,Err)
  CALL cmfe_Basis_InterpolationXiSet(LinearBasis,[CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
   & CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION,CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(LinearBasis, &
   & [NumberOfGaussPoints,NumberOfGaussPoints,NumberOfGaussPoints],Err)
  CALL cmfe_Basis_CreateFinish(LinearBasis,Err)


  ! CREATE A SECOND LINEAR BASIS FOR THE 1D GRID
  !Define basis functions - tri-Linear Lagrange
  CALL cmfe_Basis_Initialise(LinearBasisM,Err)
  CALL cmfe_Basis_CreateStart(LinearBasisUserNumberM,LinearBasisM,Err)
  CALL cmfe_Basis_TypeSet(LinearBasisM,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(LinearBasisM,1,Err)
  CALL cmfe_Basis_InterpolationXiSet(LinearBasisM,[CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(LinearBasisM,[2],Err)
  CALL cmfe_Basis_CreateFinish(LinearBasisM,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a mesh with three-dimensional elements
  CALL cmfe_Mesh_Initialise(MeshFE,Err)
  CALL cmfe_Mesh_CreateStart(MeshUserNumberFE,RegionFE,NumberOfMeshDimensionsFE,MeshFE,Err)
  CALL cmfe_Mesh_NumberOfComponentsSet(MeshFE,NumberOfMeshComponentsFE,Err) 
  CALL cmfe_Mesh_NumberOfElementsSet(MeshFE,NumberOfElementsFE,Err)  
  !Define nodes for the mesh
  CALL cmfe_Nodes_Initialise(Nodes,Err)
  CALL cmfe_Nodes_CreateStart(RegionFE,NumberOfNodesFE,Nodes,Err)
  CALL cmfe_Nodes_CreateFinish(Nodes,Err)


  ! README: 
  !
  ! Arbitrary user numbers in OpenCMISS does not really work. It is easier to leave the user numbers 
  ! identical to the global numbers (default).
  ! A mapping from the original CMISS numbering to the OpenCMISS consecutive numbering is therefore required:
  !
  ! AllNodes contains all node numbers that occur in the original CMISS mesh. These numbers are not consecutive. 
  !     This is a problem in OpenCMISS.
  ! Element1,... contain the node numbers corresponding to each Finite Element in the original CMISS numbering.
  ! 
  ! quadraticElement1,... contain the node numbers corresponding to each Finite Element in the new OpenCMISS numbering (1-81).
  DO i=1,27
    DO k=1,NumberOfNodesFE
      DO j=1,NumberOfElementsFE
        IF(AllElementNodes(j,i)==AllNodes(k)) AllQuadraticElements(j,i)=k
      ENDDO
    ENDDO
  ENDDO

  CALL cmfe_MeshElements_Initialise(QuadraticElements,Err)
  CALL cmfe_MeshElements_CreateStart(MeshFE,QuadraticMeshComponentNumber,QuadraticBasis,QuadraticElements,Err)
  DO j=1,NumberOfElementsFE
    CALL cmfe_MeshElements_NodesSet(QuadraticElements,j,AllQuadraticElements(j,:),Err)
  ENDDO
  CALL cmfe_MeshElements_CreateFinish(QuadraticElements,Err)

  !for the linear elements, we need only specific entries from the elements above
  CALL cmfe_MeshElements_Initialise(LinearElements,Err)
  CALL cmfe_MeshElements_CreateStart(MeshFE,LinearMeshComponentNumber,LinearBasis,LinearElements,Err)
  DO j=1,NumberOfElementsFE
    CALL cmfe_MeshElements_NodesSet(LinearElements,j,AllQuadraticElements(j,ENTRIES),Err)
  ENDDO
  CALL cmfe_MeshElements_CreateFinish(LinearElements,Err)

  CALL cmfe_Mesh_CreateFinish(MeshFE,Err) 
  


  ! CREATE A SECOND MESH
  !Create a mesh in the region
  CALL cmfe_Mesh_Initialise(MeshM,Err)
  CALL cmfe_Mesh_CreateStart(MeshUserNumberM,RegionM,1,MeshM,Err) ! 1=NumberOfMeshDimensionsM
  CALL cmfe_Mesh_NumberOfComponentsSet(MeshM,1,Err) ! 1=NumberOfComponentsM
  CALL cmfe_Mesh_NumberOfElementsSet(MeshM,NumberOfElementsM,Err)

  CALL cmfe_MeshElements_Initialise(ElementsM,Err)
  CALL cmfe_MeshElements_CreateStart(MeshM,MonodomainMeshComponentNumber,LinearBasisM,ElementsM,Err)

  !Define nodes for the mesh
  CALL cmfe_Nodes_Initialise(Nodes,Err)
  CALL cmfe_Nodes_CreateStart(RegionM,NumberOfNodesM,Nodes,Err)
  CALL cmfe_Nodes_CreateFinish(Nodes,Err)

    
  elem_idx=0
!  DO node1=1,NumberOfNodesMmuscle
  DO node1=1,NumberOfNodesMmuscle-1
    IF(mod(node1,NumberOfNodesInXi1)==0) CYCLE
    elem_idx=elem_idx+1
    CALL cmfe_MeshElements_NodesSet(ElementsM,elem_idx,[node1,node1+1],Err)
!    WRITE(*,*) elem_idx,node1,node1+1
  ENDDO    
  DO node1=NumberOfNodesMmuscle+1,NumberOfNodesM-1
!    IF(mod(node1-NumberOfNodesMmuscle,NumberOfNodesInXi1skin)==0) CYCLE
    IF(mod(node1,NumberOfNodesInXi1skin)==0) CYCLE
    elem_idx=elem_idx+1
    CALL cmfe_MeshElements_NodesSet(ElementsM,elem_idx,[node1,node1+1],Err)
!    WRITE(*,*) elem_idx,node1,node1+1
  ENDDO    
  write(*,*) "Finished setting up 1D elements"


  CALL cmfe_MeshElements_CreateFinish(ElementsM,Err)
  CALL cmfe_Mesh_CreateFinish(MeshM,Err) 



  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a decomposition
!  CALL cmfe_Decomposition_Initialise(DecompositionFE,Err)
!  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumberFE,MeshFE,DecompositionFE,Err)
!  CALL cmfe_Decomposition_TypeSet(DecompositionFE,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
!  CALL cmfe_Decomposition_NumberOfDomainsSet(DecompositionFE,NumberOfDomains,Err)
!  CALL cmfe_Decomposition_CalculateFacesSet(DecompositionFE,.TRUE.,Err)
!  CALL cmfe_Decomposition_CreateFinish(DecompositionFE,Err)

  CALL cmfe_Decomposition_Initialise(DecompositionFE,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumberFE,MeshFE,DecompositionFE,Err)
  IF(NumberOfDomains>1) THEN

    if(NumberOfDomains==2) then
    
      CALL cmfe_Decomposition_TypeSet(DecompositionFE,CMFE_DECOMPOSITION_USER_DEFINED_TYPE,Err)

      !set FE elements 2,3,4 and 31,32,33 to domain 0
      DO elem_idx=1,4
        CALL cmfe_Decomposition_ElementDomainSet(DecompositionFE,elem_idx,0,Err)
      ENDDO !elem_idx
      DO elem_idx=31,33
        CALL cmfe_Decomposition_ElementDomainSet(DecompositionFE,elem_idx,0,Err)
      ENDDO !elem_idx

      !set all other FE elements to domain 1
      DO elem_idx=5,31
        CALL cmfe_Decomposition_ElementDomainSet(DecompositionFE,elem_idx,1,Err)
      ENDDO !elem_idx
      DO elem_idx=34,39
        CALL cmfe_Decomposition_ElementDomainSet(DecompositionFE,elem_idx,1,Err)
      ENDDO !elem_idx
      CALL cmfe_Decomposition_NumberOfDomainsSet(DecompositionFE,NumberOfDomains,Err)
    
    else
      CALL cmfe_Decomposition_TypeSet(DecompositionFE,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
      CALL cmfe_Decomposition_NumberOfDomainsSet(DecompositionFE,NumberOfDomains,Err)
      CALL cmfe_Decomposition_CalculateFacesSet(DecompositionFE,.TRUE.,Err)
    endif
    
  ELSE
    CALL cmfe_Decomposition_TypeSet(DecompositionFE,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
    CALL cmfe_Decomposition_NumberOfDomainsSet(DecompositionFE,NumberOfDomains,Err)
  ENDIF
  CALL cmfe_Decomposition_CreateFinish(DecompositionFE,Err)






  ! CREATE A SECOND DECOMPOSITION
  CALL cmfe_Decomposition_Initialise(DecompositionM,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumberM,MeshM,DecompositionM,Err)
  IF(NumberOfDomains>1) THEN
!    CALL cmfe_Decomposition_NumberOfDomainsSet(DecompositionM,NumberOfDomains,Err)
    CALL cmfe_Decomposition_TypeSet(DecompositionM,CMFE_DECOMPOSITION_USER_DEFINED_TYPE,Err)
    DO elem_idx=1,NumberOfElementsFE
      CALL cmfe_Decomposition_ElementDomainGet(DecompositionFE,elem_idx,domain_idx,Err)
      DO elem_idx2=(elem_idx-1)*NumberOfElementsPerElasticityElement+1,elem_idx*NumberOfElementsPerElasticityElement
        CALL cmfe_Decomposition_ElementDomainSet(DecompositionM,elem_idx2,domain_idx,Err)
      ENDDO
    ENDDO !elem_idx
    CALL cmfe_Decomposition_NumberOfDomainsSet(DecompositionM,NumberOfDomains,Err)
  ELSE
    CALL cmfe_Decomposition_TypeSet(DecompositionM,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
    CALL cmfe_Decomposition_NumberOfDomainsSet(DecompositionM,NumberOfDomains,Err)
  ENDIF
  CALL cmfe_Decomposition_CreateFinish(DecompositionM,Err)


  !================================================================================================================================
  !  F I N I T E   E L A S T C I T Y
  !================================================================================================================================
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a geometric field for finite elasticity - quadratic interpolation
  CALL cmfe_Field_Initialise(GeometricFieldFE,Err)
  CALL cmfe_Field_CreateStart(FieldGeometryUserNumberFE,RegionFE,GeometricFieldFE,Err)
  CALL cmfe_Field_TypeSet(GeometricFieldFE,CMFE_FIELD_GEOMETRIC_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(GeometricFieldFE,DecompositionFE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(GeometricFieldFE,1,Err)
  CALL cmfe_Field_NumberOfComponentsSet(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,3,Err)  
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_VariableLabelSet(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
  CALL cmfe_Field_CreateFinish(GeometricFieldFE,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Update the geometric field parameters for the finite elasticity geometric field 
  open(unit=2,file="input/FinalModell.ipnode",iostat=stat)
  do
    read(2,*,iostat=stat) blub,blub2,blub3,blub4,node_idx
    if(stat<0) exit !end of file
    !search for the lines that start with "Node" and store the node number
    if(blub==node_str) then
      !read the Xj(1), Xj(2), and Xj(3) position of the node
      read(2,*,iostat=stat) blub,blub2,blub3,blub4,blub5,blub6,posX
      if(stat/=0) write(*,*) "ERROR!!!"
      read(2,*,iostat=stat) blub,blub2,blub3,blub4,blub5,blub6,posY
      if(stat/=0) write(*,*) "ERROR!!!"
      read(2,*,iostat=stat) blub,blub2,blub3,blub4,blub5,blub6,posZ
      if(stat/=0) write(*,*) "ERROR!!!"
      !find the index of this node in AllNodes
      do i=1,size(AllNodes)
        if(node_idx==AllNodes(i)) then
          my_node_idx=i
          exit
        endif
      enddo
      !print the node number and its position
      !write(*,*) my_node_idx,posX,posY,posZ
      CALL cmfe_Decomposition_NodeDomainGet(DecompositionFE,my_node_idx,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL cmfe_Field_ParameterSetUpdateNode(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
         & my_node_idx,1,posX/scalefactor,Err)
        CALL cmfe_Field_ParameterSetUpdateNode(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
         & my_node_idx,2,posY/scalefactor,Err)
        CALL cmfe_Field_ParameterSetUpdateNode(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
         & my_node_idx,3,posZ/scalefactor,Err)
      ENDIF
    endif
  enddo
  close(unit=2)
  write(*,*) "Finished reading file: input/FinalModell.ipnode"

  CALL cmfe_Field_ParameterSetUpdateStart(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a fibre field and attach it to the geometric field - quadratic interpolation
  CALL cmfe_Field_Initialise(FibreField,Err)
  CALL cmfe_Field_CreateStart(FieldFibreUserNumber,RegionFE,FibreField,Err)
  CALL cmfe_Field_TypeSet(FibreField,CMFE_FIELD_FIBRE_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(FibreField,DecompositionFE,Err)
  CALL cmfe_Field_GeometricFieldSet(FibreField,GeometricFieldFE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(FibreField,1,Err)
  CALL cmfe_Field_NumberOfComponentsSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,3,Err)  
  CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err) 
  CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_VariableLabelSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,"Fibre",Err)
  CALL cmfe_Field_CreateFinish(FibreField,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a material field for Finite Elasticity and attach it to the geometric field - quadratic interpolation
  CALL cmfe_Field_Initialise(MaterialFieldFE,Err)
  CALL cmfe_Field_CreateStart(FieldMaterialUserNumberFE,RegionFE,MaterialFieldFE,Err)
  CALL cmfe_Field_TypeSet(MaterialFieldFE,CMFE_FIELD_MATERIAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(MaterialFieldFE,DecompositionFE,Err)
  CALL cmfe_Field_GeometricFieldSet(MaterialFieldFE,GeometricFieldFE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(MaterialFieldFE,2,Err)
  CALL cmfe_Field_VariableTypesSet(MaterialFieldFE,[CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_V_VARIABLE_TYPE],Err)
  CALL cmfe_Field_NumberOfComponentsSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,8,Err)
  CALL cmfe_Field_NumberOfComponentsSet(MaterialFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,1,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,2,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,3,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,4,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,5,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,6,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,7,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,8,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,1,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_VariableLabelSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,"MaterialFE",Err)
  CALL cmfe_Field_VariableLabelSet(MaterialFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,"Gravity",Err)
  CALL cmfe_Field_CreateFinish(MaterialFieldFE,Err)

  !Set Material-Parameters [mu(1) mu(2) mu(3) alpha(1) alpha(2) alpha(3) mu_0]
  CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,C(1),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,C(2),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,C(3),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4,C(4),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,5,C(5),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,6,C(6),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,7,C(7),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,8,C(8),Err)
!  CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,0.0_CMISSRP,Err)

  !change the material in the fat/skin to passive
  do elem_idx=13,NumberOfElementsFE
    CALL cmfe_Decomposition_ElementDomainGet(DecompositionFE,elem_idx,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetUpdateElement(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & elem_idx,8,0.0_CMISSRP,Err)
    ENDIF
  enddo

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the dependent field for FE with 2 variables and * components 
  CALL cmfe_Field_Initialise(DependentFieldFE,Err)
  CALL cmfe_Field_CreateStart(FieldDependentUserNumberFE,RegionFE,DependentFieldFE,Err)
  CALL cmfe_Field_TypeSet(DependentFieldFE,CMFE_FIELD_GEOMETRIC_GENERAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(DependentFieldFE,DecompositionFE,Err)
  CALL cmfe_Field_GeometricFieldSet(DependentFieldFE,GeometricFieldFE,Err)
  CALL cmfe_Field_DependentTypeSet(DependentFieldFE,CMFE_FIELD_DEPENDENT_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(DependentFieldFE,2,Err)
  CALL cmfe_Field_VariableTypesSet(DependentFieldFE,[CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE],Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,4,Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentFieldFE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,4,LinearMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4,LinearMeshComponentNumber,Err)
!  CALL cmfe_Field_ScalingTypeSet(DependentFieldFE,CMFE_FIELD_UNIT_SCALING,Err)
  CALL cmfe_Field_VariableLabelSet(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,"DependentFE",Err)
  CALL cmfe_Field_VariableLabelSet(DependentFieldFE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,"Reaction_Force",Err)
  CALL cmfe_Field_CreateFinish(DependentFieldFE,Err)


  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the independent field for the active stress in the mechanics mesh
!  independent_field_auto_create = .TRUE.
  CALL cmfe_Field_Initialise(IndependentFieldFE,Err)
  IF(.NOT. independent_field_auto_create) THEN
    CALL cmfe_Field_CreateStart(FieldIndependentUserNumberFE,RegionFE,IndependentFieldFE,Err)
    CALL cmfe_Field_TypeSet(IndependentFieldFE,CMFE_FIELD_GENERAL_TYPE,Err)
    CALL cmfe_Field_MeshDecompositionSet(IndependentFieldFE,DecompositionFE,Err)
    CALL cmfe_Field_GeometricFieldSet(IndependentFieldFE,GeometricFieldFE,Err)
    CALL cmfe_Field_DependentTypeSet(IndependentFieldFE,CMFE_FIELD_INDEPENDENT_TYPE,Err)
    CALL cmfe_Field_NumberOfVariablesSet(IndependentFieldFE,2,Err)
    CALL cmfe_Field_VariableTypesSet(IndependentFieldFE,[CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_V_VARIABLE_TYPE],Err)
    CALL cmfe_Field_NumberOfComponentsSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,5,Err)
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1, &
     & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) !A_1
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,2, &
     & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) !A_2
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,3, &
     & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) !x_1
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,4, &
     & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) !x_2
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,5, &
     & CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) !lambda_a
    CALL cmfe_Field_ComponentMeshComponentSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
    CALL cmfe_Field_DataTypeSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_DP_TYPE,Err)
    CALL cmfe_Field_VariableLabelSet(IndependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,"XB_state_variables_FE",Err)

    CALL cmfe_Field_NumberOfComponentsSet(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,4,Err)
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,1, &
     & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err) !nodes in Xi1
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,2, &
     & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err) !nodes in Xi2
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,3, &
     & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err) !nodes in Xi3
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,4, &
     & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err) !do fibres begin in this element?
    CALL cmfe_Field_DataTypeSet(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_INTG_TYPE,Err)
    CALL cmfe_Field_VariableLabelSet(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,"subgrid_info",Err)

    CALL cmfe_Field_CreateFinish(IndependentFieldFE,Err)
  ENDIF

  !UPDATE THE INDEPENDENT FIELD IndependentFieldFE
  !second variable of IndependentFieldFE
  !  components:
  !    1) number of nodes in Xi(1) direction per element
  !    2) number of nodes in Xi(2) direction per element
  !    3) number of nodes in Xi(3) direction per element
  !    4) beginning of fibres in this FE element??? 1=yes, 0=no
  !
  !initialise as if the fibres would start in any element, and adjust below
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
   & NumberOfNodesInXi1,Err)
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2, &
   & NumberOfNodesInXi2,Err)
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3, &
   & NumberOfNodesInXi3,Err)
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4, &
   & 1,Err)


  !fewer nodes per fibre in the skin
  DO elem_idx=30,33
    CALL cmfe_Decomposition_ElementDomainGet(DecompositionFE,elem_idx,ElementDomain,Err)
    IF(ElementDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetUpdateElement(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
       & elem_idx,1,NumberOfNodesInXi1skin,Err)
    ENDIF
  ENDDO

  !no nodes in the rest of the skin
  DO elem_idx=13,29
    CALL cmfe_Decomposition_ElementDomainGet(DecompositionFE,elem_idx,ElementDomain,Err)
    IF(ElementDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetUpdateElement(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
       & elem_idx,1,0,Err)
      CALL cmfe_Field_ParameterSetUpdateElement(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
       & elem_idx,2,0,Err)
      CALL cmfe_Field_ParameterSetUpdateElement(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
       & elem_idx,3,0,Err)
    ENDIF
  ENDDO
  DO elem_idx=34,39
    CALL cmfe_Decomposition_ElementDomainGet(DecompositionFE,elem_idx,ElementDomain,Err)
    IF(ElementDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetUpdateElement(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
       & elem_idx,1,0,Err)
      CALL cmfe_Field_ParameterSetUpdateElement(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
       & elem_idx,2,0,Err)
      CALL cmfe_Field_ParameterSetUpdateElement(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
       & elem_idx,3,0,Err)
    ENDIF
  ENDDO


!  !fibres are starting in elements 1,4,7,10,...
!  DO elem_idx=1,NumberOfElementsFE,NumberGlobalXElements/NumberOfInSeriesFibres
!    CALL cmfe_Decomposition_ElementDomainGet(DecompositionFE,elem_idx,ElementDomain,Err)
!    IF(ElementDomain==ComputationalNodeNumber) THEN
!      !fibres begin in this element
!      CALL cmfe_Field_ParameterSetUpdateElement(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
!       & elem_idx,4,1,Err) 
!      CALL cmfe_Field_ParameterSetUpdateElement(IndependentFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
!       & elem_idx,1,NumberOfNodesInXi1,Err)
!    ENDIF
!  ENDDO

  !================================================================================================================================
  !  M O N O D O M A I N
  !================================================================================================================================
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a geometric field for monodomain - quadratic interpolation
  CALL cmfe_Field_Initialise(GeometricFieldM,Err)
  CALL cmfe_Field_CreateStart(FieldGeometryUserNumberM,RegionM,GeometricFieldM,Err)
  CALL cmfe_Field_TypeSet(GeometricFieldM,CMFE_FIELD_GEOMETRIC_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(GeometricFieldM,DecompositionM,Err)
  CALL cmfe_Field_TypeSet(GeometricFieldM,CMFE_FIELD_GEOMETRIC_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(GeometricFieldM,1,Err)
  CALL cmfe_Field_NumberOfComponentsSet(GeometricFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,Err)
  CALL cmfe_Field_VariableLabelSet(GeometricFieldM,CMFE_FIELD_U_VARIABLE_TYPE,"GeometryM",Err)
  CALL cmfe_Field_CreateFinish(GeometricFieldM,Err)


  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a materials field for monodomain and attach it to the geometric field - constant interpolation
  CALL cmfe_Field_Initialise(MaterialFieldM,Err)
  CALL cmfe_Field_CreateStart(FieldMaterialUserNumberM,RegionM,MaterialFieldM,Err)
  CALL cmfe_Field_TypeSet(MaterialFieldM,CMFE_FIELD_MATERIAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(MaterialFieldM,DecompositionM,Err)
  CALL cmfe_Field_GeometricFieldSet(MaterialFieldM,GeometricFieldM,Err)
  CALL cmfe_Field_NumberOfVariablesSet(MaterialFieldM,1,Err)
  CALL cmfe_Field_NumberOfComponentsSet(MaterialFieldM,CMFE_FIELD_U_VARIABLE_TYPE,3,Err) !FieldMaterialNumberOfComponentsM
  CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
!  CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldM,CMFE_FIELD_U_VARIABLE_TYPE,2,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldM,CMFE_FIELD_U_VARIABLE_TYPE,2,CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldM,CMFE_FIELD_U_VARIABLE_TYPE,3,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_VariableLabelSet(MaterialFieldM,CMFE_FIELD_U_VARIABLE_TYPE,"MaterialM",Err)
  CALL cmfe_Field_CreateFinish(MaterialFieldM,Err)

  IF(biophysical) THEN
    !Set Am
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,Am,Err)
    !Set Cm
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,Cm_slow,Err)
    !Set Conductivity
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3, &
     & CONDUCTIVITY,Err)
  ELSE
    !Set Am
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
     & 1.0_CMISSRP,Err)
    !Set Cm
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2, &
     & 1.0_CMISSRP,Err)
    !Set Conductivity
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3, &
     & 0.5_CMISSRP,Err)
  ENDIF
  

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the dependent field for monodomain with 2 variables and 1 components 
  CALL cmfe_Field_Initialise(DependentFieldM,Err)
  CALL cmfe_Field_CreateStart(FieldDependentUserNumberM,RegionM,DependentFieldM,Err)
  CALL cmfe_Field_TypeSet(DependentFieldM,CMFE_FIELD_GENERAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(DependentFieldM,DecompositionM,Err)
  CALL cmfe_Field_GeometricFieldSet(DependentFieldM,GeometricFieldM,Err)
  CALL cmfe_Field_DependentTypeSet(DependentFieldM,CMFE_FIELD_DEPENDENT_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(DependentFieldM,3,Err)
  CALL cmfe_Field_VariableTypesSet(DependentFieldM,[CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE, &
   & CMFE_FIELD_V_VARIABLE_TYPE],Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentFieldM,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,3,Err)
  !additional v_variable_type with 3 components for the 3D position of the monodomain mesh nodes
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,MonodomainMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldM,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,MonodomainMeshComponentNumber,Err)
!  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,1,MonodomainMeshComponentNumber,Err)
  CALL cmfe_Field_DimensionSet(DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL cmfe_Field_DimensionSet(DependentFieldM,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,CMFE_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL cmfe_Field_VariableLabelSet(DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,"Vm",Err)
  CALL cmfe_Field_VariableLabelSet(DependentFieldM,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,"dVm/dt",Err)
  CALL cmfe_Field_VariableLabelSet(DependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,"GeometryM3D",Err)
  CALL cmfe_Field_CreateFinish(DependentFieldM,Err)


  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the independent field for the active stress in the electrics mesh
!  independent_field_auto_create = .TRUE.
  CALL cmfe_Field_Initialise(IndependentFieldM,Err)
  IF(.NOT. independent_field_auto_create) THEN
    CALL cmfe_Field_CreateStart(FieldIndependentUserNumberM,RegionM,IndependentFieldM,Err)
    CALL cmfe_Field_TypeSet(IndependentFieldM,CMFE_FIELD_GENERAL_TYPE,Err)
    CALL cmfe_Field_MeshDecompositionSet(IndependentFieldM,DecompositionM,Err)
    CALL cmfe_Field_GeometricFieldSet(IndependentFieldM,GeometricFieldM,Err)
    CALL cmfe_Field_DependentTypeSet(IndependentFieldM,CMFE_FIELD_INDEPENDENT_TYPE,Err)
    CALL cmfe_Field_NumberOfVariablesSet(IndependentFieldM,4,Err)
    CALL cmfe_Field_VariableTypesSet(IndependentFieldM,[CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_V_VARIABLE_TYPE, &
     & CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_U2_VARIABLE_TYPE],Err)
    
    !first variable:   CMFE_FIELD_U_VARIABLE_TYPE
    CALL cmfe_Field_DataTypeSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_DP_TYPE,Err)
    CALL cmfe_Field_NumberOfComponentsSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,4,Err)
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err) !A_1
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,2, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err) !A_2
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,3, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err) !x_1
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,4, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err) !x_2
    CALL cmfe_Field_VariableLabelSet(IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,"XB_state_variables_M",Err)

    !second variable:   CMFE_FIELD_V_VARIABLE_TYPE -- 1) motor unit number   2) fibre type   3) fibre number   4) nearest Gauss point   5) in element number (LOCAL NODE NUMBERING!!!)
    CALL cmfe_Field_DataTypeSet(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_INTG_TYPE,Err)
    CALL cmfe_Field_NumberOfComponentsSet(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,5,Err)
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,1, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,2, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,3, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,4, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,5, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_VariableLabelSet(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,"fibre_info",Err)

    !third variable:   FIELD_U1_VARIABLE_TYPE -- 1) half-sarcomere length   2) inital half-sarcomere length   3) initial node distance
    CALL cmfe_Field_DataTypeSet(IndependentFieldM,CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_DP_TYPE,Err)
    CALL cmfe_Field_NumberOfComponentsSet(IndependentFieldM,CMFE_FIELD_U1_VARIABLE_TYPE,3,Err)
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U1_VARIABLE_TYPE,1, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U1_VARIABLE_TYPE,2, &
     & CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U1_VARIABLE_TYPE,3, &
     & CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL cmfe_Field_VariableLabelSet(IndependentFieldM,CMFE_FIELD_U1_VARIABLE_TYPE,"half-sarcomere_length",Err)

    !fourth variable:   FIELD_U2_VARIABLE_TYPE -- 1) old node distance   2) maximum contraction velocity   3) relative contraction velocity   4) velocity before 1 time step   5) velocity before 2 time step   6) velocity before 3 time steps
    CALL cmfe_Field_DataTypeSet(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_DP_TYPE,Err)
    CALL cmfe_Field_NumberOfComponentsSet(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,6,Err)
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,1, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,2, &
     & CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,3, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,4, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,5, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,6, &
     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_VariableLabelSet(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,"contraction_velocity",Err)

    CALL cmfe_Field_CreateFinish(IndependentFieldM,Err)
  ENDIF
  
  !initialise all motor unit numbers to 11, which is not stimulated (e.g. in the fat/skin)
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,11,Err)

  !UPDATE THE INDEPENDENT FIELD IndependentFieldM
  !first variable
  !  components: 
  !    1) active stress
  !
  !second variable
  !  components: 
  !    1) motor unit number
  !    2) fibre type
  !    3) fibre number
  !    4) nearest Gauss point
  !    5) in element number (LOCAL NODE NUMBERING!!!)
  !
  !subgrid info
!  open(unit=4,file="input/EntireTA.exgrid",iostat=stat)
!  read(4,*,iostat=stat) blub
!  k=0
!  m=1
!  do while(k<NumberOfNodesInXi2*NumberOfNodesInXi3*12)
!    read(4,*,iostat=stat) node_idx,mu_nr,Ftype,fibre_nr,posX,posY,posZ,NearestGP,InElement
!    if(stat<0) exit !end of file

!    do i=1,NumberOfNodesInXi1
!!      CALL cmfe_Decomposition_NodeDomainGet(DecompositionM,node_idx,1,NodeDomain,Err)
!      CALL cmfe_Decomposition_NodeDomainGet(DecompositionM,m,1,NodeDomain,Err)
!      IF(NodeDomain==ComputationalNodeNumber) THEN
!!        if((mu_nr.LE.0).OR.(mu_nr.GT.11).OR.(node_idx.GT.(NumberOfNodesMmuscle))) mu_nr=11 !tomo TODO
!!        CALL cmfe_Field_ParameterSetUpdateNode(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
!!         & node_idx,1,mu_nr,Err)
!!        CALL cmfe_Field_ParameterSetUpdateNode(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
!!         & node_idx,2,Ftype,Err)
!        if((mu_nr.LE.0).OR.(mu_nr.GT.11).OR.(m.GT.(NumberOfNodesMmuscle))) mu_nr=11 !tomo TODO
!        CALL cmfe_Field_ParameterSetUpdateNode(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
!         & m,1,mu_nr,Err)
!        CALL cmfe_Field_ParameterSetUpdateNode(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
!         & m,2,Ftype,Err)
!      ENDIF
!      m=m+1
!    enddo

!!    k=k+1
!    k=k+5
!  enddo
!  close(unit=4)

!  write(*,*) "Finished reading file: input/EntireTA.exgrid"

  open(unit=4,file="input/MU_fibres.txt",iostat=stat)
  k=0
  m=1
  Ftype=1;
  do while(k<NumberOfNodesInXi2*NumberOfNodesInXi3*12)
    read(4,*,iostat=stat) mu_nr
    if(stat<0) exit !end of file

    CALL cmfe_Decomposition_NodeDomainGet(DecompositionM,m,1,NodeDomain,Err) !this only works if the fibers are not cut in the decomposition
    do i=1,NumberOfNodesInXi1
!      CALL cmfe_Decomposition_NodeDomainGet(DecompositionM,node_idx,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        if((mu_nr.LE.0).OR.(mu_nr.GT.11).OR.(m.GT.NumberOfNodesMmuscle)) mu_nr=11 !tomo TODO
        CALL cmfe_Field_ParameterSetUpdateNode(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
         & m,1,mu_nr,Err)
        CALL cmfe_Field_ParameterSetUpdateNode(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
         & m,2,Ftype,Err)
      ENDIF
      m=m+1
    enddo

    k=k+1
  enddo
  close(unit=4)
  write(*,*) "Finished reading file: input/MU_fibres.txt"

  !init the fibre number, the nearest Gauss point info and the inElem info to 0
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,0,Err)
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4,0,Err)
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,5,0,Err) !(LOCAL NODE NUMBERING!!!)

  !third variable:
  !  components:
  !    1) sarcomere half length
  !    2) initial sarcomere half length
  !    3) initial node distance
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2, &
   & 1.0_CMISSRP,Err)
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3, &
   & LENGTH/(NumberOfInSeriesFibres*NumberOfNodesPerFibre-1),Err)
  !fourth variable:
  !  components:
  !    1) old node distance
  !    2) maximum contraction velocity
  !    3) relative contraction velocity
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
   & LENGTH/(NumberOfInSeriesFibres*NumberOfNodesPerFibre-1),Err)
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2, &
   & Vmax/(NumberOfInSeriesFibres*NumberOfNodesPerFibre-1),Err)
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3, &
   & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4, &
   & Vmax/(NumberOfInSeriesFibres*NumberOfNodesPerFibre-1),Err)
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,5, &
   & Vmax/(NumberOfInSeriesFibres*NumberOfNodesPerFibre-1),Err)
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,6, &
   & Vmax/(NumberOfInSeriesFibres*NumberOfNodesPerFibre-1),Err)


  !================================================================================================================================
  !  EQUATIONS SET
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the equations_set for Finite Elasticity
  CALL cmfe_Field_Initialise(EquationsSetFieldFE,Err)
  CALL cmfe_EquationsSet_Initialise(EquationsSetFE,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationsSetsUserNumberFE,RegionFE,FibreField,[CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
!   & CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMFE_EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE], & 
   & CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMFE_EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE], & 
   & EquationsSetFieldUserNumberFE,EquationsSetFieldFE,EquationsSetFE,Err)
  CALL cmfe_EquationsSet_CreateFinish(EquationsSetFE,Err)

  !Create the equations set dependent field variables for Finite Elasticity
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSetFE,FieldDependentUserNumberFE,DependentFieldFE,Err)
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSetFE,Err)

  !Create the equations set independent field variable for the active stress component for Finite Elasticity
  CALL cmfe_EquationsSet_IndependentCreateStart(EquationsSetFE,FieldIndependentUserNumberFE,IndependentFieldFE,Err)
  CALL cmfe_EquationsSet_IndependentCreateFinish(EquationsSetFE,Err)

  !Create the equations set materials field variables for Finite Elasticity
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSetFE,FieldMaterialUserNumberFE,MaterialFieldFE,Err)  
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSetFE,Err)


  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the equations_set for monodomain
  CALL cmfe_Field_Initialise(EquationsSetFieldM,Err)
  CALL cmfe_EquationsSet_Initialise(EquationsSetM,Err)
  !Set the equations set to be a Monodomain equations set
  !> \todo solve the monodomain problem on the fibre field rather than on the geometric field: GeometricField <--> FibreField
  CALL cmfe_EquationsSet_CreateStart(EquationsSetsUserNumberM,RegionM,GeometricFieldM,[CMFE_EQUATIONS_SET_BIOELECTRICS_CLASS, &
!   & CMFE_EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE,CMFE_EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE], &
   & CMFE_EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE,CMFE_EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE], &
   & EquationsSetFieldUserNumberM,EquationsSetFieldM,EquationsSetM,Err)
  CALL cmfe_EquationsSet_CreateFinish(EquationsSetM,Err)

  !Create the equations set dependent field variables for monodomain
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSetM,FieldDependentUserNumberM,DependentFieldM,Err)
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSetM,Err)

  !Create the equations set independent field variable for the active stress component for Finite Elasticity
  CALL cmfe_EquationsSet_IndependentCreateStart(EquationsSetM,FieldIndependentUserNumberM,IndependentFieldM,Err)
  CALL cmfe_EquationsSet_IndependentCreateFinish(EquationsSetM,Err)

  !Create the equations set materials field variables for monodomain
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSetM,FieldMaterialUserNumberM,MaterialFieldM,Err)
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSetM,Err)


  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the equations set equations for monodomain
  CALL cmfe_Equations_Initialise(EquationsM,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSetM,EquationsM,Err)
  CALL cmfe_Equations_SparsityTypeSet(EquationsM,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  CALL cmfe_Equations_OutputTypeSet(EquationsM,CMFE_EQUATIONS_NO_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_TIMING_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_MATRIX_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSetM,Err)

  !Create the equations set equations for Finite Elasticity
  CALL cmfe_Equations_Initialise(EquationsFE,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSetFE,EquationsFE,Err)
  CALL cmfe_Equations_SparsityTypeSet(EquationsFE,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  CALL cmfe_Equations_OutputTypeSet(EquationsFE,CMFE_EQUATIONS_NO_OUTPUT,Err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSetFE,Err)


  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the CellML environment
  CALL cmfe_CellML_Initialise(CellML,Err)
  CALL cmfe_CellML_CreateStart(CellMLUserNumber,RegionM,CellML,Err)
  !Import the Shorten et al. 2007 model from a file
  ! Now we have imported all the models we are able to specify which variables from the model we want:
  !,- to set from this side

  IF(biophysical) THEN
    !biophysical half-sarcomere model
    CALL cmfe_CellML_ModelImport(CellML,filename,shortenModelIndex,Err)
    CALL cmfe_CellML_ModelImport(CellML,filename2,shortenModelIndex2,Err)
    CALL cmfe_CellML_VariableSetAsKnown(CellML,shortenModelIndex,"wal_environment/I_HH",Err)
    CALL cmfe_CellML_VariableSetAsKnown(CellML,shortenModelIndex2,"wal_environment/I_HH",Err)
    !,- to get from the CellML side
  !  CALL cmfe_CellML_VariableSetAsWanted(CellML,shortenModelIndex,"wal_environment/I_T",Err)
  !  CALL cmfe_CellML_VariableSetAsWanted(CellML,shortenModelIndex,"wal_environment/I_ionic_s",Err)
  !  CALL cmfe_CellML_VariableSetAsWanted(CellML,shortenModelIndex,"wal_environment/I_ionic_t",Err)
    !
    !NOTE: If an INTERMEDIATE (or ALGEBRAIC) variable should be used in a mapping, it has to be set as known or wanted first!
    !,  --> set "razumova/stress" as wanted!
    !,  --> no need to set "wal_environment/vS" since all STATE variables are automatically set as wanted! 
    CALL cmfe_CellML_VariableSetAsWanted(CellML,shortenModelIndex,"razumova/stress",Err)
    CALL cmfe_CellML_VariableSetAsWanted(CellML,shortenModelIndex2,"razumova/stress",Err)
    !,- and override constant parameters without needing to set up fields
    !> \todo Need to allow parameter values to be overridden for the case when user has non-spatially varying parameter value.
    !Finish the CellML environment
  ELSE
    !phenomenological half-sarcomere model
    CALL cmfe_CellML_ModelImport(CellML,filename,shortenModelIndex,Err)
    CALL cmfe_CellML_VariableSetAsKnown(CellML,shortenModelIndex,"Aliev_Panfilov/I_HH",Err)
    CALL cmfe_CellML_VariableSetAsKnown(CellML,shortenModelIndex,"Razumova/l_hs",Err)
    CALL cmfe_CellML_VariableSetAsKnown(CellML,shortenModelIndex,"Razumova/velo",Err)
  END IF
  
  CALL cmfe_CellML_CreateFinish(CellML,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the CellML <--> OpenCMISS field maps
  CALL cmfe_CellML_FieldMapsCreateStart(CellML,Err)
  IF(biophysical) THEN
    !Map the transmembrane voltage V_m
    CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
     & shortenModelIndex,"wal_environment/vS",CMFE_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
     & shortenModelIndex2,"wal_environment/vS",CMFE_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,shortenModelIndex,"wal_environment/vS",CMFE_FIELD_VALUES_SET_TYPE, &
     & DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,shortenModelIndex2,"wal_environment/vS",CMFE_FIELD_VALUES_SET_TYPE, &
     & DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,shortenModelIndex,"razumova/stress",CMFE_FIELD_VALUES_SET_TYPE, &
     & IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,shortenModelIndex2,"razumova/stress",CMFE_FIELD_VALUES_SET_TYPE, &
     & IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)
  ELSE
    !Map the half-sarcomere length L_S
    CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,IndependentFieldM,CMFE_FIELD_U1_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
     & shortenModelIndex,"Razumova/l_hs",CMFE_FIELD_VALUES_SET_TYPE,Err)
    !Map the sarcomere relative contraction velocity
    CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,3,CMFE_FIELD_VALUES_SET_TYPE, &
     & shortenModelIndex,"Razumova/velo",CMFE_FIELD_VALUES_SET_TYPE,Err)
    !Map the transmembrane voltage V_m
    CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
     & shortenModelIndex,"Aliev_Panfilov/V_m",CMFE_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,shortenModelIndex,"Aliev_Panfilov/V_m",CMFE_FIELD_VALUES_SET_TYPE, &
     & DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)
    !Map the active stress
    CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,shortenModelIndex,"Razumova/A_1",CMFE_FIELD_VALUES_SET_TYPE, &
     & IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,shortenModelIndex,"Razumova/A_2",CMFE_FIELD_VALUES_SET_TYPE, &
     & IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,2,CMFE_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,shortenModelIndex,"Razumova/x_1",CMFE_FIELD_VALUES_SET_TYPE, &
     & IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,3,CMFE_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,shortenModelIndex,"Razumova/x_2",CMFE_FIELD_VALUES_SET_TYPE, &
     & IndependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,4,CMFE_FIELD_VALUES_SET_TYPE,Err)
  END IF
  CALL cmfe_CellML_FieldMapsCreateFinish(CellML,Err)
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Initialise dependent field for monodomain
  !> \todo - get V_m initialial value.
  
  IF(biophysical) THEN
    !biophysical half-sarcomere model
    CALL cmfe_Field_ComponentValuesInitialise(DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
     & -79.974_CMISSRP,Err)
  ELSE
    !phenomenological half-sarcomere model
    CALL cmfe_Field_ComponentValuesInitialise(DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
     & 0.0_CMISSRP,Err)
  END IF
  
  !Initialise dependent field for Finite Elasticity from undeformed geometry and set hydrostatic pressure
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE, &
   & CMFE_FIELD_VALUES_SET_TYPE,1,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,Err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE, & 
   & CMFE_FIELD_VALUES_SET_TYPE,2,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,Err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE, &
   & CMFE_FIELD_VALUES_SET_TYPE,3,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,Err)

!  INIT_PRESSURE=-2.0_CMISSRP*MAT_FE(2)-MAT_FE(1)
  INIT_PRESSURE=0.0_CMISSRP
  CALL cmfe_Field_ComponentValuesInitialise(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4, & 
   & INIT_PRESSURE,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the CellML models field
  CALL cmfe_Field_Initialise(CellMLModelsField,Err)
  CALL cmfe_CellML_ModelsFieldCreateStart(CellML,CellMLModelsFieldUserNumber,CellMLModelsField,Err)
  CALL cmfe_CellML_ModelsFieldCreateFinish(CellML,Err)

  !Set up the models field
  IF(biophysical) THEN
    CALL cmfe_Field_ComponentValuesInitialise(CellMLModelsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
     & 0,Err)
    do NodeNumber=1,NumberOfNodesMmuscle
      CALL cmfe_Decomposition_NodeDomainGet(DecompositionM,NodeNumber,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL cmfe_Field_ParameterSetGetNode(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
          & NodeNumber,1,mu_nr,Err)
        if(mu_nr.LE.6) then
          CALL cmfe_Field_ParameterSetUpdateNode(CellMLModelsField,CMFE_FIELD_U_VARIABLE_TYPE, &
            & CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,1,shortenModelIndex,Err)
        else
          CALL cmfe_Field_ParameterSetUpdateNode(CellMLModelsField,CMFE_FIELD_U_VARIABLE_TYPE, &
            & CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,1,shortenModelIndex2,Err)
          CALL cmfe_Field_ParameterSetUpdateNode(MaterialFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
            & NodeNumber,2,Cm_fast,Err)
        endif
      ENDIF
    enddo
  ELSE
    CALL cmfe_Field_ComponentValuesInitialise(CellMLModelsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
     & 0,Err)
    do NodeNumber=1,NumberOfNodesMmuscle
      CALL cmfe_Decomposition_NodeDomainGet(DecompositionM,NodeNumber,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL cmfe_Field_ParameterSetGetNode(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
          & NodeNumber,1,mu_nr,Err)
        CALL cmfe_Field_ParameterSetUpdateNode(CellMLModelsField,CMFE_FIELD_U_VARIABLE_TYPE, &
          & CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,1,shortenModelIndex,Err)
      ENDIF
    enddo
  ENDIF

!  DO NodeNumber=NumberOfNodesPerFibre/2,NumberOfNodesM,NumberOfNodesPerFibre
!    CALL cmfe_Decomposition_NodeDomainGet(DecompositionM,NodeNumber,1,NodeDomain,Err)
!    IF(NodeDomain==ComputationalNodeNumber) THEN
!      CALL cmfe_Field_ParameterSetUpdateNode(CellMLModelsField,CMFE_FIELD_U_VARIABLE_TYPE, &
!        & CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,1,shortenModelIndex2,Err)
!    ENDIF
!  ENDDO

  CALL cmfe_Field_ParameterSetUpdateStart(CellMLModelsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(CellMLModelsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)


  !Create the CellML state field
  CALL cmfe_Field_Initialise(CellMLStateField,Err)
  CALL cmfe_CellML_StateFieldCreateStart(CellML,CellMLStateFieldUserNumber,CellMLStateField,Err)
  CALL cmfe_CellML_StateFieldCreateFinish(CellML,Err)

  IF(biophysical) THEN
    !Create the CellML intermediate field
    CALL cmfe_Field_Initialise(CellMLIntermediateField,Err)
    CALL cmfe_CellML_IntermediateFieldCreateStart(CellML,CellMLIntermediateFieldUserNumber,CellMLIntermediateField,Err)
    CALL cmfe_CellML_IntermediateFieldCreateFinish(CellML,Err)
  END IF
  
  !Create the CellML parameters field
  CALL cmfe_Field_Initialise(CellMLParametersField,Err)
  CALL cmfe_CellML_ParametersFieldCreateStart(CellML,CellMLParametersFieldUserNumber,CellMLParametersField,Err)
  CALL cmfe_CellML_ParametersFieldCreateFinish(CellML,Err)
  

  !--------------------------------------------------------------------------------------------------------------------------------
  !Define the problem
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_MULTI_PHYSICS_CLASS, &
!   & CMFE_PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE,CMFE_PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE],Problem,Err)
   & CMFE_PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE,CMFE_PROBLEM_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE],Problem,Err)
  CALL cmfe_Problem_CreateFinish(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)

  !set the main control loop (time loop type)
  CALL cmfe_ControlLoop_Initialise(ControlLoopMain,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoopMain,Err)
  CALL cmfe_ControlLoop_LabelSet(ControlLoopMain,'MAIN_TIME_LOOP',Err)
  !Loop in time for STIM_STOP with the Stimulus applied.
  CALL cmfe_ControlLoop_TimesSet(ControlLoopMain,0.0_CMISSRP,ELASTICITY_TIME_STEP,ELASTICITY_TIME_STEP,Err)
  CALL cmfe_ControlLoop_TimeOutputSet(ControlLoopMain,OUTPUT_FREQUENCY,Err)
  CALL cmfe_ControlLoop_OutputTypeSet(ControlLoopMain,CMFE_CONTROL_LOOP_TIMING_OUTPUT,Err)

  !set the monodomain loop (time loop type)
  CALL cmfe_ControlLoop_Initialise(ControlLoopM,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,[ControlLoopMonodomainNumber,CMFE_CONTROL_LOOP_NODE],ControlLoopM,Err)
  CALL cmfe_ControlLoop_LabelSet(ControlLoopM,'MONODOMAIN_TIME_LOOP',Err)
  CALL cmfe_ControlLoop_TimesSet(ControlLoopM,0.0_CMISSRP,ELASTICITY_TIME_STEP,PDE_TIME_STEP,Err)
  CALL cmfe_ControlLoop_OutputTypeSet(ControlLoopM,CMFE_CONTROL_LOOP_NO_OUTPUT,Err)

  !set the finite elasticity loop (simple type)
  CALL cmfe_ControlLoop_Initialise(ControlLoopFE,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,[ControlLoopElasticityNumber,CMFE_CONTROL_LOOP_NODE],ControlLoopFE,Err)
  CALL cmfe_ControlLoop_LabelSet(ControlLoopFE,'ELASTICITY_LOOP',Err)

  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the problem solvers
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)

  !Create the DAE solver
  CALL cmfe_Solver_Initialise(SolverDAE,Err)
  CALL cmfe_Problem_SolverGet(Problem,[ControlLoopMonodomainNumber,CMFE_CONTROL_LOOP_NODE], &
   & SolverDAEIndex,SolverDAE,Err)
  CALL cmfe_Solver_DAETimeStepSet(SolverDAE,ODE_TIME_STEP,Err)
  !> \todo - solve the CellML equations on the GPU for efficiency (later)
  !CALL cmfe_Solver_DAESolverTypeSet(SolverDAE,CMFE_SOLVER_DAE_EXTERNAL,Err) 
  CALL cmfe_Solver_OutputTypeSet(SolverDAE,CMFE_SOLVER_NO_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverDAE,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverDAE,CMFE_SOLVER_TIMING_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverDAE,CMFE_SOLVER_SOLVER_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverDAE,CMFE_SOLVER_MATRIX_OUTPUT,Err)

  !Create the parabolic solver
  CALL cmfe_Solver_Initialise(SolverParabolic,Err)
  CALL cmfe_Problem_SolverGet(Problem,[ControlLoopMonodomainNumber,CMFE_CONTROL_LOOP_NODE], &
   & SolverParabolicIndex,SolverParabolic,Err)
  CALL cmfe_Solver_DynamicSchemeSet(SolverParabolic,CMFE_SOLVER_DYNAMIC_BACKWARD_EULER_SCHEME,Err)
  CALL cmfe_Solver_OutputTypeSet(SolverParabolic,CMFE_SOLVER_NO_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverParabolic,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverParabolic,CMFE_SOLVER_TIMING_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverParabolic,CMFE_SOLVER_SOLVER_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverParabolic,CMFE_SOLVER_MATRIX_OUTPUT,Err)

  !Create the Finte Elasticity solver
  CALL cmfe_Solver_Initialise(SolverFE,Err)
  CALL cmfe_Solver_Initialise(LinearSolverFE,Err)
  CALL cmfe_Problem_SolverGet(Problem,[ControlLoopElasticityNumber,CMFE_CONTROL_LOOP_NODE], &
   & SolverFEIndex,SolverFE,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverFE,CMFE_SOLVER_NO_OUTPUT,Err)
  CALL cmfe_Solver_OutputTypeSet(SolverFE,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverFE,CMFE_SOLVER_TIMING_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverFE,CMFE_SOLVER_SOLVER_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverFE,CMFE_SOLVER_MATRIX_OUTPUT,Err)
!  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(SolverFE,CMFE_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(SolverFE,CMFE_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,Err)
  CALL cmfe_Solver_NewtonMaximumIterationsSet(SolverFE,100,Err)
  CALL cmfe_Solver_NewtonAbsoluteToleranceSet(SolverFE,1.E-6_CMISSRP,Err) !6
  CALL cmfe_Solver_NewtonSolutionToleranceSet(SolverFE,2.E-6_CMISSRP,Err)
  CALL cmfe_Solver_NewtonLinearSolverGet(SolverFE,LinearSolverFE,Err)
  CALL cmfe_Solver_LinearTypeSet(LinearSolverFE,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)

  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the problem solver CellML equations
  CALL cmfe_Problem_CellMLEquationsCreateStart(Problem,Err)

  CALL cmfe_Solver_Initialise(SolverDAE,Err)
  CALL cmfe_CellMLEquations_Initialise(CellMLEquations,Err)
  CALL cmfe_Problem_SolverGet(Problem,[ControlLoopMonodomainNumber,CMFE_CONTROL_LOOP_NODE],SolverDAEIndex,SolverDAE,Err)
  CALL cmfe_Solver_CellMLEquationsGet(SolverDAE,CellMLEquations,Err)
  CALL cmfe_CellMLEquations_CellMLAdd(CellMLEquations,CellML,CellMLIndex,Err)

  CALL cmfe_Problem_CellMLEquationsCreateFinish(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)

  !Create the problem solver parabolic equations (Monodomain)
  CALL cmfe_Solver_Initialise(SolverParabolic,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquationsM,Err)
  CALL cmfe_Problem_SolverGet(Problem,[ControlLoopMonodomainNumber,CMFE_CONTROL_LOOP_NODE], &
   & SolverParabolicIndex,SolverParabolic,Err)
  CALL cmfe_Solver_SolverEquationsGet(SolverParabolic,SolverEquationsM,Err)
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquationsM,CMFE_SOLVER_SPARSE_MATRICES,Err)
  !CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquationsM,CMFE_SOLVER_FULL_MATRICES,Err)  
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquationsM,EquationsSetM,EquationsSetIndexM,Err)

  !Create the problem solver Finite Elasticity equations
  CALL cmfe_Solver_Initialise(SolverFE,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquationsFE,Err)
  CALL cmfe_Problem_SolverGet(Problem,[ControlLoopElasticityNumber,CMFE_CONTROL_LOOP_NODE], &
   & SolverFEIndex,SolverFE,Err)
  CALL cmfe_Solver_SolverEquationsGet(SolverFE,SolverEquationsFE,Err)
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquationsFE,CMFE_SOLVER_SPARSE_MATRICES,Err)
  !CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquationsFE,CMFE_SOLVER_FULL_MATRICES,Err)  
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquationsFE,EquationsSetFE,EquationsSetIndexFE,Err)

  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !boundary conditions

  !Prescribe boundary conditions for monodomain
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditionsM,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquationsM,BoundaryConditionsM,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquationsM,Err)

  !Prescribe boundary conditions for Finite Elasticity (absolute nodal parameters)
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditionsFE,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquationsFE,BoundaryConditionsFE,Err)



    
  !Fix TOP_NODES in all directions
  do i=1,21
    NodeNumber=TOP_NODES(i)
    CALL cmfe_Decomposition_NodeDomainGet(DecompositionFE,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetGetNode(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber, &
       & 1,VALUE,Err)
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
       & CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
      CALL cmfe_Field_ParameterSetGetNode(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber, & 
       & 2,VALUE,Err)
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,2, &
       & CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
      CALL cmfe_Field_ParameterSetGetNode(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber, &
       & 3,VALUE,Err)
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,3, &
       & CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
    ENDIF
  enddo

  !Fix BOTTOM_NODES in all directions
  do i=1,21
    NodeNumber=BOTTOM_NODES(i)
    CALL cmfe_Decomposition_NodeDomainGet(DecompositionFE,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetGetNode(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber, &
       & 1,VALUE,Err)
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
       & CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
      CALL cmfe_Field_ParameterSetGetNode(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber, & 
       & 2,VALUE,Err)
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,2, &
       & CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
      CALL cmfe_Field_ParameterSetGetNode(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber, &
       & 3,VALUE,Err)
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,3, &
       & CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
    ENDIF
  enddo

  !Fix SIDE_NODES in all directions
  do i=1,52
    NodeNumber=SIDE_NODES(i)
    CALL cmfe_Decomposition_NodeDomainGet(DecompositionFE,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetGetNode(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber, &
       & 1,VALUE,Err)
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
       & CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
      CALL cmfe_Field_ParameterSetGetNode(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber, & 
       & 2,VALUE,Err)
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,2, &
       & CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
      CALL cmfe_Field_ParameterSetGetNode(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber, &
       & 3,VALUE,Err)
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,3, &
       & CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
    ENDIF
  enddo


  
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquationsFE,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Calculate the bioelectrics geometric field 
  CALL cmfe_ControlLoop_Initialise(ControlLoopM,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,[ControlLoopMonodomainNumber,CMFE_CONTROL_LOOP_NODE],ControlLoopM,Err)
  CALL cmfe_BioelectricsFiniteElasticity_UpdateGeometricField(ControlLoopM,.TRUE.,Err)

  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL cmfe_Fields_Initialise(Fields,Err)
    CALL cmfe_Fields_Create(RegionM,Fields,Err)
    CALL cmfe_Fields_NodesExport(Fields,"EntireTAExample_M","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields,"EntireTAExample_M","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields,Err)

    CALL cmfe_Fields_Initialise(Fields,Err)
    CALL cmfe_Fields_Create(RegionFE,Fields,Err)
    CALL cmfe_Fields_NodesExport(Fields,"EntireTAExample_FE","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields,"EntireTAExample_FE","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields,Err)
  ENDIF




  !--------------------------------------------------------------------------------------------------------------------------------
  !Solve the problem


  !Solve the problem -- bring to new length before applying the stimulus
  WRITE(*,'(A)') "Start solve before stimulation"
  CALL cmfe_Problem_Solve(Problem,Err)

  !reset the relative contraction velocity to 0
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldM,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3, &
   & 0.0_CMISSRP,Err)

  CALL cmfe_ControlLoop_Initialise(ControlLoopFE,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,[ControlLoopElasticityNumber,CMFE_CONTROL_LOOP_NODE],ControlLoopFE,Err)
  CALL cmfe_ControlLoop_MaximumIterationsSet(ControlLoopFE,1,Err)


! no change for BCs -- fix at this length!!!


  !--------------------------------------------------------------------------------------------------------------------------------
  !Read in the MU firing times

  k=1
  open(unit=5,file="input/MU_firing_times.txt",iostat=stat)
  do while(.true.)
    read(5,*,iostat=stat) FT_1,FT_2,FT_3,FT_4,FT_5,FT_6,FT_7,FT_8,FT_9,FT_10 
    if(stat<0) exit !end of file

    FIRING_TIMES(k,:)=[FT_1,FT_2,FT_3,FT_4,FT_5,FT_6,FT_7,FT_8,FT_9,FT_10]

    k=k+1
  enddo
  close(unit=5)
  write(*,*) "Finished reading file: input/MU_firing_times.txt"
  




  !--------------------------------------------------------------------------------------------------------------------------------
  !Solve the problem
  time = 0.0_CMISSRP
  m=1 
  k=1
  do while(time <= TIME_STOP)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Set the Stimulus for monodomain at the middle of the fibres
  IF(biophysical) THEN
    CALL cmfe_CellML_FieldComponentGet(CellML,shortenModelIndex,CMFE_CELLML_PARAMETERS_FIELD, &
     & "wal_environment/I_HH",stimcomponent,Err)
  ELSE
    CALL cmfe_CellML_FieldComponentGet(CellML,shortenModelIndex,CMFE_CELLML_PARAMETERS_FIELD, &
     & "Aliev_Panfilov/I_HH",stimcomponent,Err)
  ENDIF

  NodeNumber=NumberOfNodesInXi1/2
  !loop over all neuromuscular junctions (middle point of the fibres)
  DO WHILE(NodeNumber<NumberOfNodesMmuscle)
    CALL cmfe_Decomposition_NodeDomainGet(DecompositionM,NodeNumber+IZ_offset(m),1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetGetNode(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
       & NodeNumber+IZ_offset(m),1,mu_nr,Err)
      if((mu_nr.LE.0).OR.(mu_nr.GE.11)) then
        mu_nr=11
      else
        val=FIRING_TIMES(k,mu_nr)
        if(val==1) then
          if(mu_nr.LE.6) then
            CALL cmfe_Field_ParameterSetUpdateNode(CellMLParametersField, &
              & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber+IZ_offset(m),stimcomponent,STIM_VALUE,Err)
          else
            CALL cmfe_Field_ParameterSetUpdateNode(CellMLParametersField, &
              & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber+IZ_offset(m),stimcomponent,STIM_VALUE2,Err)
          endif
        endif
      endif
    ENDIF
    NodeNumber=NodeNumber+NumberOfNodesInXi1
    m=m+1
  ENDDO
  m=m-NumberOfNodesInXi2*NumberOfNodesInXi3*12

  !--------------------------------------------------------------------------------------------------------------------------------
  !Solve the problem for the stimulation time
  CALL cmfe_ControlLoop_TimesSet(ControlLoopMain,time,time+STIM_STOP,ELASTICITY_TIME_STEP,Err)
  CALL cmfe_Problem_Solve(Problem,Err)


  !--------------------------------------------------------------------------------------------------------------------------------
  !Now turn the stimulus off
  NodeNumber=NumberOfNodesInXi1/2
  DO WHILE(NodeNumber<NumberOfNodesMmuscle)
    CALL cmfe_Decomposition_NodeDomainGet(DecompositionM,NodeNumber+IZ_offset(m),1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) CALL cmfe_Field_ParameterSetUpdateNode(CellMLParametersField, & 
     & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber+IZ_offset(m),stimcomponent,0.0_CMISSRP,Err)
    NodeNumber=NodeNumber+NumberOfNodesInXi1
    m=m+1
  ENDDO
  m=m-NumberOfNodesInXi2*NumberOfNodesInXi3*12

  !--------------------------------------------------------------------------------------------------------------------------------
  !Solve the problem for the rest of the period
  CALL cmfe_ControlLoop_TimesSet(ControlLoopMain,time+STIM_STOP,time+PERIODD,ELASTICITY_TIME_STEP,Err)
  CALL cmfe_Problem_Solve(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  time = time + PERIODD
  k=k+1

  end do
  
  !--------------------------------------------------------------------------------------------------------------------------------
  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL cmfe_Fields_Initialise(Fields,Err)
    CALL cmfe_Fields_Create(RegionM,Fields,Err)
    CALL cmfe_Fields_NodesExport(Fields,"EntireTAExample_M","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields,"EntireTAExample_M","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields,Err)

    CALL cmfe_Fields_Initialise(Fields,Err)
    CALL cmfe_Fields_Create(RegionFE,Fields,Err)
    CALL cmfe_Fields_NodesExport(Fields,"EntireTAExample_FE","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields,"EntireTAExample_FE","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields,Err)
  ENDIF
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Finialise CMISS
  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
END PROGRAM EntireTAEXAMPLE
