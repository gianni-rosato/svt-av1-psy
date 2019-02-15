/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

/*
* Copyright (c) 2016, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at www.aomedia.org/license/software. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at www.aomedia.org/license/patent.
*/

#include <stdlib.h>

#include "EbDefinitions.h"
#include "EbUtility.h"
#include "EbSequenceControlSet.h"
#include "EbPictureControlSet.h"
#include "EbModeDecisionConfigurationProcess.h"
#include "EbRateControlResults.h"
#include "EbEncDecTasks.h"
#include "EbModeDecisionConfiguration.h"
#include "EbLambdaRateTables.h"
#include "EbReferenceObject.h"
#include "EbModeDecisionProcess.h"


// Shooting states
#define UNDER_SHOOTING                        0
#define OVER_SHOOTING                         1
#define TBD_SHOOTING                          2


// Set a cost to each search method (could be modified)
// EB30 @ Revision 12879
#define PRED_OPEN_LOOP_1_NFL_COST    97 // PRED_OPEN_LOOP_1_NFL_COST is ~03% faster than PRED_OPEN_LOOP_COST
#define PRED_OPEN_LOOP_COST         100 // Let's assume PRED_OPEN_LOOP_COST costs ~100 U
#define U_104                       104

#define LIGHT_OPEN_LOOP_COST        106 // L_MDC is ~06% slower than PRED_OPEN_LOOP_COST
#define U_108                       108
#define U_109                       109
#define OPEN_LOOP_COST              110 // F_MDC is ~10% slower than PRED_OPEN_LOOP_COST
#define U_111                       111
#define U_112                       112
#define U_113                       113
#define U_114                       114
#define U_115                       115
#define U_116                       116
#define U_117                       117
#define U_120                       120
#define U_121                       121
#define LIGHT_BDP_COST              123 // L_BDP is ~23% slower than PRED_OPEN_LOOP_COST
#define U_125                       125
#define U_127                       127
#define BDP_COST                    129 // F_BDP is ~29% slower than PRED_OPEN_LOOP_COST
#define U_130                       130
#define U_133                       133
#define U_134                       134
#define AVC_COST                    138 // L_BDP is ~38% slower than PRED_OPEN_LOOP_COST
#define U_145                       145
#define U_150                       150
#define FULL_SEARCH_COST            155 // FS    is ~55% slower than PRED_OPEN_LOOP_COST


// ADP SB score Manipulation
#define ADP_CLASS_SHIFT_DIST_0    50

#define ADP_DARK_LCU_TH           25
#define ADP_LIGHT_LCU_TH         225


#define ADP_NON_MOVING_INDEX_TH_0 15
#define ADP_NON_MOVING_INDEX_TH_1 30
#define LUMINOSITY_CHANGE_TH       3



#define VALID_SLOT_TH                        2


// Coefficient scaling and quantization with AV1 TX are tailored to
// the AV1 TX transforms.  Regardless of the bit-depth of the input,
// the transform stages scale the coefficient values up by a factor of
// 8 (3 bits) over the scale of the pixel values.  Thus, for 8-bit
// input, the coefficients have effectively 11 bits of scale depth
// (8+3), 10-bit input pixels result in 13-bit coefficient depth
// (10+3) and 12-bit pixels yield 15-bit (12+3) coefficient depth.
// All quantizers are built using this invariant of x8, 3-bit scaling,
// thus the Q3 suffix.

// A partial exception to this rule is large transforms; to avoid
// overflow, TX blocks with > 256 pels (>16x16) are scaled only
// 4-times unity (2 bits) over the pixel depth, and TX blocks with
// over 1024 pixels (>32x32) are scaled up only 2x unity (1 bit).
// This descaling is found via av1_tx_get_scale().  Thus, 16x32, 32x16
// and 32x32 transforms actually return Q2 coefficients, and 32x64,
// 64x32 and 64x64 transforms return Q1 coefficients.  However, the
// quantizers are de-scaled down on-the-fly by the same amount
// (av1_tx_get_scale()) during quantization, and as such the
// dequantized/decoded coefficients, even for large TX blocks, are always
// effectively Q3. Meanwhile, quantized/coded coefficients are Q0
// because Qn quantizers are applied to Qn tx coefficients.

// Note that encoder decision making (which uses the quantizer to
// generate several bespoke lamdas for RDO and other heuristics)
// expects quantizers to be larger for higher-bitdepth input.  In
// addition, the minimum allowable quantizer is 4; smaller values will
// underflow to 0 in the actual quantization routines.

static const int16_t dc_qlookup_Q3[QINDEX_RANGE] = {
    4, 8, 8, 9, 10, 11, 12, 12, 13, 14, 15, 16, 17, 18,
    19, 19, 20, 21, 22, 23, 24, 25, 26, 26, 27, 28, 29, 30,
    31, 32, 32, 33, 34, 35, 36, 37, 38, 38, 39, 40, 41, 42,
    43, 43, 44, 45, 46, 47, 48, 48, 49, 50, 51, 52, 53, 53,
    54, 55, 56, 57, 57, 58, 59, 60, 61, 62, 62, 63, 64, 65,
    66, 66, 67, 68, 69, 70, 70, 71, 72, 73, 74, 74, 75, 76,
    77, 78, 78, 79, 80, 81, 81, 82, 83, 84, 85, 85, 87, 88,
    90, 92, 93, 95, 96, 98, 99, 101, 102, 104, 105, 107, 108, 110,
    111, 113, 114, 116, 117, 118, 120, 121, 123, 125, 127, 129, 131, 134,
    136, 138, 140, 142, 144, 146, 148, 150, 152, 154, 156, 158, 161, 164,
    166, 169, 172, 174, 177, 180, 182, 185, 187, 190, 192, 195, 199, 202,
    205, 208, 211, 214, 217, 220, 223, 226, 230, 233, 237, 240, 243, 247,
    250, 253, 257, 261, 265, 269, 272, 276, 280, 284, 288, 292, 296, 300,
    304, 309, 313, 317, 322, 326, 330, 335, 340, 344, 349, 354, 359, 364,
    369, 374, 379, 384, 389, 395, 400, 406, 411, 417, 423, 429, 435, 441,
    447, 454, 461, 467, 475, 482, 489, 497, 505, 513, 522, 530, 539, 549,
    559, 569, 579, 590, 602, 614, 626, 640, 654, 668, 684, 700, 717, 736,
    755, 775, 796, 819, 843, 869, 896, 925, 955, 988, 1022, 1058, 1098, 1139,
    1184, 1232, 1282, 1336,
};

static const int16_t dc_qlookup_10_Q3[QINDEX_RANGE] = {
    4, 9, 10, 13, 15, 17, 20, 22, 25, 28, 31, 34, 37,
    40, 43, 47, 50, 53, 57, 60, 64, 68, 71, 75, 78, 82,
    86, 90, 93, 97, 101, 105, 109, 113, 116, 120, 124, 128, 132,
    136, 140, 143, 147, 151, 155, 159, 163, 166, 170, 174, 178, 182,
    185, 189, 193, 197, 200, 204, 208, 212, 215, 219, 223, 226, 230,
    233, 237, 241, 244, 248, 251, 255, 259, 262, 266, 269, 273, 276,
    280, 283, 287, 290, 293, 297, 300, 304, 307, 310, 314, 317, 321,
    324, 327, 331, 334, 337, 343, 350, 356, 362, 369, 375, 381, 387,
    394, 400, 406, 412, 418, 424, 430, 436, 442, 448, 454, 460, 466,
    472, 478, 484, 490, 499, 507, 516, 525, 533, 542, 550, 559, 567,
    576, 584, 592, 601, 609, 617, 625, 634, 644, 655, 666, 676, 687,
    698, 708, 718, 729, 739, 749, 759, 770, 782, 795, 807, 819, 831,
    844, 856, 868, 880, 891, 906, 920, 933, 947, 961, 975, 988, 1001,
    1015, 1030, 1045, 1061, 1076, 1090, 1105, 1120, 1137, 1153, 1170, 1186, 1202,
    1218, 1236, 1253, 1271, 1288, 1306, 1323, 1342, 1361, 1379, 1398, 1416, 1436,
    1456, 1476, 1496, 1516, 1537, 1559, 1580, 1601, 1624, 1647, 1670, 1692, 1717,
    1741, 1766, 1791, 1817, 1844, 1871, 1900, 1929, 1958, 1990, 2021, 2054, 2088,
    2123, 2159, 2197, 2236, 2276, 2319, 2363, 2410, 2458, 2508, 2561, 2616, 2675,
    2737, 2802, 2871, 2944, 3020, 3102, 3188, 3280, 3375, 3478, 3586, 3702, 3823,
    3953, 4089, 4236, 4394, 4559, 4737, 4929, 5130, 5347,
};

static const int16_t dc_qlookup_12_Q3[QINDEX_RANGE] = {
    4, 12, 18, 25, 33, 41, 50, 60, 70, 80, 91,
    103, 115, 127, 140, 153, 166, 180, 194, 208, 222, 237,
    251, 266, 281, 296, 312, 327, 343, 358, 374, 390, 405,
    421, 437, 453, 469, 484, 500, 516, 532, 548, 564, 580,
    596, 611, 627, 643, 659, 674, 690, 706, 721, 737, 752,
    768, 783, 798, 814, 829, 844, 859, 874, 889, 904, 919,
    934, 949, 964, 978, 993, 1008, 1022, 1037, 1051, 1065, 1080,
    1094, 1108, 1122, 1136, 1151, 1165, 1179, 1192, 1206, 1220, 1234,
    1248, 1261, 1275, 1288, 1302, 1315, 1329, 1342, 1368, 1393, 1419,
    1444, 1469, 1494, 1519, 1544, 1569, 1594, 1618, 1643, 1668, 1692,
    1717, 1741, 1765, 1789, 1814, 1838, 1862, 1885, 1909, 1933, 1957,
    1992, 2027, 2061, 2096, 2130, 2165, 2199, 2233, 2267, 2300, 2334,
    2367, 2400, 2434, 2467, 2499, 2532, 2575, 2618, 2661, 2704, 2746,
    2788, 2830, 2872, 2913, 2954, 2995, 3036, 3076, 3127, 3177, 3226,
    3275, 3324, 3373, 3421, 3469, 3517, 3565, 3621, 3677, 3733, 3788,
    3843, 3897, 3951, 4005, 4058, 4119, 4181, 4241, 4301, 4361, 4420,
    4479, 4546, 4612, 4677, 4742, 4807, 4871, 4942, 5013, 5083, 5153,
    5222, 5291, 5367, 5442, 5517, 5591, 5665, 5745, 5825, 5905, 5984,
    6063, 6149, 6234, 6319, 6404, 6495, 6587, 6678, 6769, 6867, 6966,
    7064, 7163, 7269, 7376, 7483, 7599, 7715, 7832, 7958, 8085, 8214,
    8352, 8492, 8635, 8788, 8945, 9104, 9275, 9450, 9639, 9832, 10031,
    10245, 10465, 10702, 10946, 11210, 11482, 11776, 12081, 12409, 12750, 13118,
    13501, 13913, 14343, 14807, 15290, 15812, 16356, 16943, 17575, 18237, 18949,
    19718, 20521, 21387,
};
static const int16_t ac_qlookup_Q3[QINDEX_RANGE] = {
    4, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
    20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
    33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45,
    46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58,
    59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71,
    72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84,
    85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97,
    98, 99, 100, 101, 102, 104, 106, 108, 110, 112, 114, 116, 118,
    120, 122, 124, 126, 128, 130, 132, 134, 136, 138, 140, 142, 144,
    146, 148, 150, 152, 155, 158, 161, 164, 167, 170, 173, 176, 179,
    182, 185, 188, 191, 194, 197, 200, 203, 207, 211, 215, 219, 223,
    227, 231, 235, 239, 243, 247, 251, 255, 260, 265, 270, 275, 280,
    285, 290, 295, 300, 305, 311, 317, 323, 329, 335, 341, 347, 353,
    359, 366, 373, 380, 387, 394, 401, 408, 416, 424, 432, 440, 448,
    456, 465, 474, 483, 492, 501, 510, 520, 530, 540, 550, 560, 571,
    582, 593, 604, 615, 627, 639, 651, 663, 676, 689, 702, 715, 729,
    743, 757, 771, 786, 801, 816, 832, 848, 864, 881, 898, 915, 933,
    951, 969, 988, 1007, 1026, 1046, 1066, 1087, 1108, 1129, 1151, 1173, 1196,
    1219, 1243, 1267, 1292, 1317, 1343, 1369, 1396, 1423, 1451, 1479, 1508, 1537,
    1567, 1597, 1628, 1660, 1692, 1725, 1759, 1793, 1828,
};

static const int16_t ac_qlookup_10_Q3[QINDEX_RANGE] = {
    4, 9, 11, 13, 16, 18, 21, 24, 27, 30, 33, 37, 40,
    44, 48, 51, 55, 59, 63, 67, 71, 75, 79, 83, 88, 92,
    96, 100, 105, 109, 114, 118, 122, 127, 131, 136, 140, 145, 149,
    154, 158, 163, 168, 172, 177, 181, 186, 190, 195, 199, 204, 208,
    213, 217, 222, 226, 231, 235, 240, 244, 249, 253, 258, 262, 267,
    271, 275, 280, 284, 289, 293, 297, 302, 306, 311, 315, 319, 324,
    328, 332, 337, 341, 345, 349, 354, 358, 362, 367, 371, 375, 379,
    384, 388, 392, 396, 401, 409, 417, 425, 433, 441, 449, 458, 466,
    474, 482, 490, 498, 506, 514, 523, 531, 539, 547, 555, 563, 571,
    579, 588, 596, 604, 616, 628, 640, 652, 664, 676, 688, 700, 713,
    725, 737, 749, 761, 773, 785, 797, 809, 825, 841, 857, 873, 889,
    905, 922, 938, 954, 970, 986, 1002, 1018, 1038, 1058, 1078, 1098, 1118,
    1138, 1158, 1178, 1198, 1218, 1242, 1266, 1290, 1314, 1338, 1362, 1386, 1411,
    1435, 1463, 1491, 1519, 1547, 1575, 1603, 1631, 1663, 1695, 1727, 1759, 1791,
    1823, 1859, 1895, 1931, 1967, 2003, 2039, 2079, 2119, 2159, 2199, 2239, 2283,
    2327, 2371, 2415, 2459, 2507, 2555, 2603, 2651, 2703, 2755, 2807, 2859, 2915,
    2971, 3027, 3083, 3143, 3203, 3263, 3327, 3391, 3455, 3523, 3591, 3659, 3731,
    3803, 3876, 3952, 4028, 4104, 4184, 4264, 4348, 4432, 4516, 4604, 4692, 4784,
    4876, 4972, 5068, 5168, 5268, 5372, 5476, 5584, 5692, 5804, 5916, 6032, 6148,
    6268, 6388, 6512, 6640, 6768, 6900, 7036, 7172, 7312,
};

static const int16_t ac_qlookup_12_Q3[QINDEX_RANGE] = {
    4, 13, 19, 27, 35, 44, 54, 64, 75, 87, 99,
    112, 126, 139, 154, 168, 183, 199, 214, 230, 247, 263,
    280, 297, 314, 331, 349, 366, 384, 402, 420, 438, 456,
    475, 493, 511, 530, 548, 567, 586, 604, 623, 642, 660,
    679, 698, 716, 735, 753, 772, 791, 809, 828, 846, 865,
    884, 902, 920, 939, 957, 976, 994, 1012, 1030, 1049, 1067,
    1085, 1103, 1121, 1139, 1157, 1175, 1193, 1211, 1229, 1246, 1264,
    1282, 1299, 1317, 1335, 1352, 1370, 1387, 1405, 1422, 1440, 1457,
    1474, 1491, 1509, 1526, 1543, 1560, 1577, 1595, 1627, 1660, 1693,
    1725, 1758, 1791, 1824, 1856, 1889, 1922, 1954, 1987, 2020, 2052,
    2085, 2118, 2150, 2183, 2216, 2248, 2281, 2313, 2346, 2378, 2411,
    2459, 2508, 2556, 2605, 2653, 2701, 2750, 2798, 2847, 2895, 2943,
    2992, 3040, 3088, 3137, 3185, 3234, 3298, 3362, 3426, 3491, 3555,
    3619, 3684, 3748, 3812, 3876, 3941, 4005, 4069, 4149, 4230, 4310,
    4390, 4470, 4550, 4631, 4711, 4791, 4871, 4967, 5064, 5160, 5256,
    5352, 5448, 5544, 5641, 5737, 5849, 5961, 6073, 6185, 6297, 6410,
    6522, 6650, 6778, 6906, 7034, 7162, 7290, 7435, 7579, 7723, 7867,
    8011, 8155, 8315, 8475, 8635, 8795, 8956, 9132, 9308, 9484, 9660,
    9836, 10028, 10220, 10412, 10604, 10812, 11020, 11228, 11437, 11661, 11885,
    12109, 12333, 12573, 12813, 13053, 13309, 13565, 13821, 14093, 14365, 14637,
    14925, 15213, 15502, 15806, 16110, 16414, 16734, 17054, 17390, 17726, 18062,
    18414, 18766, 19134, 19502, 19886, 20270, 20670, 21070, 21486, 21902, 22334,
    22766, 23214, 23662, 24126, 24590, 25070, 25551, 26047, 26559, 27071, 27599,
    28143, 28687, 29247,
};


int16_t av1_dc_quant_Q3(int32_t qindex, int32_t delta, aom_bit_depth_t bit_depth) {
    switch (bit_depth) {
    case AOM_BITS_8: return dc_qlookup_Q3[clamp(qindex + delta, 0, MAXQ)];
    case AOM_BITS_10: return dc_qlookup_10_Q3[clamp(qindex + delta, 0, MAXQ)];
    case AOM_BITS_12: return dc_qlookup_12_Q3[clamp(qindex + delta, 0, MAXQ)];
    default:
        assert(0 && "bit_depth should be AOM_BITS_8, AOM_BITS_10 or AOM_BITS_12");
        return -1;
    }
}
int16_t av1_ac_quant_Q3(int32_t qindex, int32_t delta, aom_bit_depth_t bit_depth) {
    switch (bit_depth) {
    case AOM_BITS_8: return ac_qlookup_Q3[clamp(qindex + delta, 0, MAXQ)];
    case AOM_BITS_10: return ac_qlookup_10_Q3[clamp(qindex + delta, 0, MAXQ)];
    case AOM_BITS_12: return ac_qlookup_12_Q3[clamp(qindex + delta, 0, MAXQ)];
    default:
        assert(0 && "bit_depth should be AOM_BITS_8, AOM_BITS_10 or AOM_BITS_12");
        return -1;
    }
}

static int32_t get_qzbin_factor(int32_t q, aom_bit_depth_t bit_depth) {
    const int32_t quant = av1_dc_quant_Q3(q, 0, bit_depth);
    switch (bit_depth) {
    case AOM_BITS_8: return q == 0 ? 64 : (quant < 148 ? 84 : 80);
    case AOM_BITS_10: return q == 0 ? 64 : (quant < 592 ? 84 : 80);
    case AOM_BITS_12: return q == 0 ? 64 : (quant < 2368 ? 84 : 80);
    default:
        assert(0 && "bit_depth should be AOM_BITS_8, AOM_BITS_10 or AOM_BITS_12");
        return -1;
    }
}

// In AV1 TX, the coefficients are always scaled up a factor of 8 (3
// bits), so QTX == Q3.

int16_t av1_dc_quant_QTX(int32_t qindex, int32_t delta, aom_bit_depth_t bit_depth) {
    return av1_dc_quant_Q3(qindex, delta, bit_depth);
}
int16_t av1_ac_quant_QTX(int32_t qindex, int32_t delta, aom_bit_depth_t bit_depth) {
    return av1_ac_quant_Q3(qindex, delta, bit_depth);
}

static void invert_quant(int16_t *quant, int16_t *shift, int32_t d) {
    uint32_t t;
    int32_t l, m;
    t = d;
    for (l = 0; t > 1; l++) t >>= 1;
    m = 1 + (1 << (16 + l)) / d;
    *quant = (int16_t)(m - (1 << 16));
    *shift = 1 << (16 - l);
}

static INLINE int32_t aom_get_qmlevel(int32_t qindex, int32_t first, int32_t last) {
    return first + (qindex * (last + 1 - first)) / QINDEX_RANGE;
}

void SetGlobalMotionField(
    PictureControlSet_t                    *picture_control_set_ptr)
{
    // Init Global Motion Vector
    uint8_t frameIndex;
    for (frameIndex = INTRA_FRAME; frameIndex <= ALTREF_FRAME; ++frameIndex) {
        picture_control_set_ptr->parent_pcs_ptr->global_motion[frameIndex].wmtype = IDENTITY;
        picture_control_set_ptr->parent_pcs_ptr->global_motion[frameIndex].alpha = 0;
        picture_control_set_ptr->parent_pcs_ptr->global_motion[frameIndex].beta = 0;
        picture_control_set_ptr->parent_pcs_ptr->global_motion[frameIndex].delta = 0;
        picture_control_set_ptr->parent_pcs_ptr->global_motion[frameIndex].gamma = 0;
        picture_control_set_ptr->parent_pcs_ptr->global_motion[frameIndex].invalid = 0;
        picture_control_set_ptr->parent_pcs_ptr->global_motion[frameIndex].wmmat[0] = 0;
        picture_control_set_ptr->parent_pcs_ptr->global_motion[frameIndex].wmmat[1] = 0;
        picture_control_set_ptr->parent_pcs_ptr->global_motion[frameIndex].wmmat[2] = (1 << WARPEDMODEL_PREC_BITS);
        picture_control_set_ptr->parent_pcs_ptr->global_motion[frameIndex].wmmat[3] = 0;
        picture_control_set_ptr->parent_pcs_ptr->global_motion[frameIndex].wmmat[4] = 0;
        picture_control_set_ptr->parent_pcs_ptr->global_motion[frameIndex].wmmat[5] = (1 << WARPEDMODEL_PREC_BITS);
        picture_control_set_ptr->parent_pcs_ptr->global_motion[frameIndex].wmmat[6] = 0;
        picture_control_set_ptr->parent_pcs_ptr->global_motion[frameIndex].wmmat[7] = 0;

    }

    //Update MV
    if (picture_control_set_ptr->parent_pcs_ptr->is_pan && picture_control_set_ptr->parent_pcs_ptr->is_tilt) {
        picture_control_set_ptr->parent_pcs_ptr->global_motion[LAST_FRAME].wmtype = TRANSLATION;
        picture_control_set_ptr->parent_pcs_ptr->global_motion[LAST_FRAME].wmmat[1] = ((picture_control_set_ptr->parent_pcs_ptr->panMvx + picture_control_set_ptr->parent_pcs_ptr->tiltMvx) / 2) << 1 << GM_TRANS_ONLY_PREC_DIFF;
        picture_control_set_ptr->parent_pcs_ptr->global_motion[LAST_FRAME].wmmat[0] = ((picture_control_set_ptr->parent_pcs_ptr->panMvy + picture_control_set_ptr->parent_pcs_ptr->tiltMvy) / 2) << 1 << GM_TRANS_ONLY_PREC_DIFF;
    }
    else if (picture_control_set_ptr->parent_pcs_ptr->is_pan) {

        picture_control_set_ptr->parent_pcs_ptr->global_motion[LAST_FRAME].wmtype = TRANSLATION;
        picture_control_set_ptr->parent_pcs_ptr->global_motion[LAST_FRAME].wmmat[1] = picture_control_set_ptr->parent_pcs_ptr->panMvx << 1 << GM_TRANS_ONLY_PREC_DIFF;
        picture_control_set_ptr->parent_pcs_ptr->global_motion[LAST_FRAME].wmmat[0] = picture_control_set_ptr->parent_pcs_ptr->panMvy << 1 << GM_TRANS_ONLY_PREC_DIFF;

    }
    else if (picture_control_set_ptr->parent_pcs_ptr->is_tilt) {

        picture_control_set_ptr->parent_pcs_ptr->global_motion[LAST_FRAME].wmtype = TRANSLATION;
        picture_control_set_ptr->parent_pcs_ptr->global_motion[LAST_FRAME].wmmat[1] = picture_control_set_ptr->parent_pcs_ptr->tiltMvx << 1 << GM_TRANS_ONLY_PREC_DIFF;
        picture_control_set_ptr->parent_pcs_ptr->global_motion[LAST_FRAME].wmmat[0] = picture_control_set_ptr->parent_pcs_ptr->tiltMvy << 1 << GM_TRANS_ONLY_PREC_DIFF;

    }

    picture_control_set_ptr->parent_pcs_ptr->global_motion[LAST_FRAME].wmmat[1] = (int32_t)clamp(picture_control_set_ptr->parent_pcs_ptr->global_motion[LAST_FRAME].wmmat[1], GM_TRANS_MIN*GM_TRANS_DECODE_FACTOR, GM_TRANS_MAX*GM_TRANS_DECODE_FACTOR);
    picture_control_set_ptr->parent_pcs_ptr->global_motion[LAST_FRAME].wmmat[0] = (int32_t)clamp(picture_control_set_ptr->parent_pcs_ptr->global_motion[LAST_FRAME].wmmat[0], GM_TRANS_MIN*GM_TRANS_DECODE_FACTOR, GM_TRANS_MAX*GM_TRANS_DECODE_FACTOR);

    picture_control_set_ptr->parent_pcs_ptr->global_motion[BWDREF_FRAME].wmtype = TRANSLATION;
    picture_control_set_ptr->parent_pcs_ptr->global_motion[BWDREF_FRAME].wmmat[1] = 0 - picture_control_set_ptr->parent_pcs_ptr->global_motion[LAST_FRAME].wmmat[1];
    picture_control_set_ptr->parent_pcs_ptr->global_motion[BWDREF_FRAME].wmmat[0] = 0 - picture_control_set_ptr->parent_pcs_ptr->global_motion[LAST_FRAME].wmmat[0];
    picture_control_set_ptr->parent_pcs_ptr->global_motion[BWDREF_FRAME].wmmat[1] = (int32_t)clamp(picture_control_set_ptr->parent_pcs_ptr->global_motion[BWDREF_FRAME].wmmat[1], GM_TRANS_MIN*GM_TRANS_DECODE_FACTOR, GM_TRANS_MAX*GM_TRANS_DECODE_FACTOR);
    picture_control_set_ptr->parent_pcs_ptr->global_motion[BWDREF_FRAME].wmmat[0] = (int32_t)clamp(picture_control_set_ptr->parent_pcs_ptr->global_motion[BWDREF_FRAME].wmmat[0], GM_TRANS_MIN*GM_TRANS_DECODE_FACTOR, GM_TRANS_MAX*GM_TRANS_DECODE_FACTOR);



    //convert_to_trans_prec(
    //    picture_control_set_ptr->parent_pcs_ptr->allow_high_precision_mv,
    //    picture_control_set_ptr->parent_pcs_ptr->global_motion[LAST_FRAME].wmmat[0]) *GM_TRANS_ONLY_DECODE_FACTOR;

    //convert_to_trans_prec(
    //    picture_control_set_ptr->parent_pcs_ptr->allow_high_precision_mv,
    //    picture_control_set_ptr->parent_pcs_ptr->global_motion[LAST_FRAME].wmmat[1]) *GM_TRANS_ONLY_DECODE_FACTOR;

}


void av1_set_quantizer(
    PictureParentControlSet_t                    *picture_control_set_ptr,
    int32_t q)
{
    // quantizer has to be reinitialized with av1_init_quantizer() if any
    // delta_q changes.
    picture_control_set_ptr->using_qmatrix = 0;
    picture_control_set_ptr->min_qmlevel = 5;
    picture_control_set_ptr->max_qmlevel = 9;

    picture_control_set_ptr->base_qindex = (uint16_t)AOMMAX(picture_control_set_ptr->delta_q_present_flag, q);
    picture_control_set_ptr->y_dc_delta_q = 0;
#if TUNE_CHROMA_OFFSET
    picture_control_set_ptr->u_dc_delta_q = (picture_control_set_ptr->slice_type == I_SLICE) ? -10 : -20;
    picture_control_set_ptr->u_ac_delta_q = (picture_control_set_ptr->slice_type == I_SLICE) ? -10 : -20;
    picture_control_set_ptr->v_dc_delta_q = (picture_control_set_ptr->slice_type == I_SLICE) ? -10 : -20;
    picture_control_set_ptr->v_ac_delta_q = (picture_control_set_ptr->slice_type == I_SLICE) ? -10 : -20;
#else
    picture_control_set_ptr->u_dc_delta_q = 0;
    picture_control_set_ptr->u_ac_delta_q = 0;
    picture_control_set_ptr->v_dc_delta_q = 0;
    picture_control_set_ptr->v_ac_delta_q = 0;
#endif
    picture_control_set_ptr->qm_y = aom_get_qmlevel(picture_control_set_ptr->base_qindex, picture_control_set_ptr->min_qmlevel, picture_control_set_ptr->max_qmlevel);
    picture_control_set_ptr->qm_u = aom_get_qmlevel(picture_control_set_ptr->base_qindex + picture_control_set_ptr->u_ac_delta_q,
        picture_control_set_ptr->min_qmlevel, picture_control_set_ptr->max_qmlevel);

    if (!picture_control_set_ptr->separate_uv_delta_q)
        picture_control_set_ptr->qm_v = picture_control_set_ptr->qm_u;
    else
        picture_control_set_ptr->qm_v = aom_get_qmlevel(picture_control_set_ptr->base_qindex + picture_control_set_ptr->v_ac_delta_q,
            picture_control_set_ptr->min_qmlevel, picture_control_set_ptr->max_qmlevel);
}

void av1_build_quantizer(
    aom_bit_depth_t bit_depth,
    int32_t y_dc_delta_q,
    int32_t u_dc_delta_q,
    int32_t u_ac_delta_q,
    int32_t v_dc_delta_q,
    int32_t v_ac_delta_q,
    Quants *const quants,
    Dequants *const deq)
{
    int32_t i, q, quant_Q3, quant_QTX;

    for (q = 0; q < QINDEX_RANGE; q++) {
        const int32_t qzbin_factor = get_qzbin_factor(q, bit_depth);
        const int32_t qrounding_factor = q == 0 ? 64 : 48;

        for (i = 0; i < 2; ++i) {
            int32_t qrounding_factor_fp = 64;
            // y quantizer setup with original coeff shift of Q3
            quant_Q3 = i == 0 ? av1_dc_quant_Q3(q, y_dc_delta_q, bit_depth)
                : av1_ac_quant_Q3(q, 0, bit_depth);
            // y quantizer with TX scale
            quant_QTX = i == 0 ? av1_dc_quant_QTX(q, y_dc_delta_q, bit_depth)
                : av1_ac_quant_QTX(q, 0, bit_depth);
            invert_quant(&quants->y_quant[q][i], &quants->y_quant_shift[q][i],
                quant_QTX);
            quants->y_quant_fp[q][i] = (int16_t)((1 << 16) / quant_QTX);
            quants->y_round_fp[q][i] = (int16_t)((qrounding_factor_fp * quant_QTX) >> 7);
            quants->y_zbin[q][i] = (int16_t)ROUND_POWER_OF_TWO(qzbin_factor * quant_QTX, 7);
            quants->y_round[q][i] = (int16_t)((qrounding_factor * quant_QTX) >> 7);
            deq->y_dequant_QTX[q][i] = (int16_t)quant_QTX;
            deq->y_dequant_Q3[q][i] = (int16_t)quant_Q3;

            // u quantizer setup with original coeff shift of Q3
            quant_Q3 = i == 0 ? av1_dc_quant_Q3(q, u_dc_delta_q, bit_depth)
                : av1_ac_quant_Q3(q, u_ac_delta_q, bit_depth);
            // u quantizer with TX scale
            quant_QTX = i == 0 ? av1_dc_quant_QTX(q, u_dc_delta_q, bit_depth)
                : av1_ac_quant_QTX(q, u_ac_delta_q, bit_depth);
            invert_quant(&quants->u_quant[q][i], &quants->u_quant_shift[q][i],
                quant_QTX);
            quants->u_quant_fp[q][i] = (int16_t)((1 << 16) / quant_QTX);
            quants->u_round_fp[q][i] = (int16_t)((qrounding_factor_fp * quant_QTX) >> 7);
            quants->u_zbin[q][i] = (int16_t)ROUND_POWER_OF_TWO(qzbin_factor * quant_QTX, 7);
            quants->u_round[q][i] = (int16_t)((qrounding_factor * quant_QTX) >> 7);
            deq->u_dequant_QTX[q][i] = (int16_t)quant_QTX;
            deq->u_dequant_Q3[q][i] = (int16_t)quant_Q3;

            // v quantizer setup with original coeff shift of Q3
            quant_Q3 = i == 0 ? av1_dc_quant_Q3(q, v_dc_delta_q, bit_depth)
                : av1_ac_quant_Q3(q, v_ac_delta_q, bit_depth);
            // v quantizer with TX scale
            quant_QTX = i == 0 ? av1_dc_quant_QTX(q, v_dc_delta_q, bit_depth)
                : av1_ac_quant_QTX(q, v_ac_delta_q, bit_depth);
            invert_quant(&quants->v_quant[q][i], &quants->v_quant_shift[q][i],
                quant_QTX);
            quants->v_quant_fp[q][i] = (int16_t)((1 << 16) / quant_QTX);
            quants->v_round_fp[q][i] = (int16_t)((qrounding_factor_fp * quant_QTX) >> 7);
            quants->v_zbin[q][i] = (int16_t)ROUND_POWER_OF_TWO(qzbin_factor * quant_QTX, 7);
            quants->v_round[q][i] = (int16_t)((qrounding_factor * quant_QTX) >> 7);
            deq->v_dequant_QTX[q][i] = (int16_t)quant_QTX;
            deq->v_dequant_Q3[q][i] = (int16_t)quant_Q3;
        }

        for (i = 2; i < 8; i++) {  // 8: SIMD width
            quants->y_quant[q][i] = quants->y_quant[q][1];
            quants->y_quant_fp[q][i] = quants->y_quant_fp[q][1];
            quants->y_round_fp[q][i] = quants->y_round_fp[q][1];
            quants->y_quant_shift[q][i] = quants->y_quant_shift[q][1];
            quants->y_zbin[q][i] = quants->y_zbin[q][1];
            quants->y_round[q][i] = quants->y_round[q][1];
            deq->y_dequant_QTX[q][i] = deq->y_dequant_QTX[q][1];
            deq->y_dequant_Q3[q][i] = deq->y_dequant_Q3[q][1];

            quants->u_quant[q][i] = quants->u_quant[q][1];
            quants->u_quant_fp[q][i] = quants->u_quant_fp[q][1];
            quants->u_round_fp[q][i] = quants->u_round_fp[q][1];
            quants->u_quant_shift[q][i] = quants->u_quant_shift[q][1];
            quants->u_zbin[q][i] = quants->u_zbin[q][1];
            quants->u_round[q][i] = quants->u_round[q][1];
            deq->u_dequant_QTX[q][i] = deq->u_dequant_QTX[q][1];
            deq->u_dequant_Q3[q][i] = deq->u_dequant_Q3[q][1];
            quants->v_quant[q][i] = quants->u_quant[q][1];
            quants->v_quant_fp[q][i] = quants->v_quant_fp[q][1];
            quants->v_round_fp[q][i] = quants->v_round_fp[q][1];
            quants->v_quant_shift[q][i] = quants->v_quant_shift[q][1];
            quants->v_zbin[q][i] = quants->v_zbin[q][1];
            quants->v_round[q][i] = quants->v_round[q][1];
            deq->v_dequant_QTX[q][i] = deq->v_dequant_QTX[q][1];
            deq->v_dequant_Q3[q][i] = deq->v_dequant_Q3[q][1];
        }
    }
}


void av1_qm_init(
    PictureParentControlSet_t                   *picture_control_set_ptr
)
{
    const uint8_t num_planes = 3;// MAX_MB_PLANE;// NM- No monochroma
    uint8_t q, c, t;
    int32_t current;
    for (q = 0; q < NUM_QM_LEVELS; ++q) {
        for (c = 0; c < num_planes; ++c) {
            current = 0;
            for (t = 0; t < TX_SIZES_ALL; ++t) {
                const int32_t size = tx_size_2d[t];
                const TxSize qm_tx_size = av1_get_adjusted_tx_size(t);
                if (q == NUM_QM_LEVELS - 1) {
                    picture_control_set_ptr->gqmatrix[q][c][t] = NULL;
                    picture_control_set_ptr->giqmatrix[q][c][t] = NULL;
                }
                else if (t != qm_tx_size) {  // Reuse matrices for 'qm_tx_size'
                    picture_control_set_ptr->gqmatrix[q][c][t] = picture_control_set_ptr->gqmatrix[q][c][qm_tx_size];
                    picture_control_set_ptr->giqmatrix[q][c][t] = picture_control_set_ptr->giqmatrix[q][c][qm_tx_size];
                }
                else {
                    assert(current + size <= QM_TOTAL_SIZE);
                    picture_control_set_ptr->gqmatrix[q][c][t] = &wt_matrix_ref[q][c >= 1][current];
                    picture_control_set_ptr->giqmatrix[q][c][t] = &iwt_matrix_ref[q][c >= 1][current];
                    current += size;
                }
            }
        }
    }
}



/******************************************************
* Compute picture and slice level chroma QP offsets
******************************************************/
void SetSliceAndPictureChromaQpOffsets(
    PictureControlSet_t                    *picture_control_set_ptr

)
{
    // This is a picture level chroma QP offset and is sent in the PPS
    picture_control_set_ptr->cb_qp_offset = 0;
    picture_control_set_ptr->cr_qp_offset = 0;


    //In order to have QP offsets for chroma at a slice level set slice_level_chroma_qp_flag flag in picture_control_set_ptr (can be done in the PCS Ctor)

    // The below are slice level chroma QP offsets and is sent for each slice when slice_level_chroma_qp_flag is set

    // IMPORTANT: Lambda tables assume that the cb and cr have the same QP offsets.
    // However the offsets for each component can be very different for ENC DEC and we are conformant.
    picture_control_set_ptr->slice_cb_qp_offset = 0;
    picture_control_set_ptr->slice_cr_qp_offset = 0;

    if (picture_control_set_ptr->parent_pcs_ptr->pic_noise_class >= PIC_NOISE_CLASS_6) {

        picture_control_set_ptr->slice_cb_qp_offset = 10;
        picture_control_set_ptr->slice_cr_qp_offset = 10;

    }
    else if (picture_control_set_ptr->parent_pcs_ptr->pic_noise_class >= PIC_NOISE_CLASS_4) {

        picture_control_set_ptr->slice_cb_qp_offset = 8;
        picture_control_set_ptr->slice_cr_qp_offset = 8;

    }
    else {
        if (picture_control_set_ptr->temporal_layer_index == 1) {
            picture_control_set_ptr->slice_cb_qp_offset = 2;
            picture_control_set_ptr->slice_cr_qp_offset = 2;
        }
        else {
            picture_control_set_ptr->slice_cb_qp_offset = 0;
            picture_control_set_ptr->slice_cr_qp_offset = 0;
        }
    }


}

#if FAST_SG
/******************************************************
* Set the reference sg ep for a given picture
******************************************************/
void set_reference_sg_ep(
    PictureControlSet_t                    *picture_control_set_ptr)
{

    Av1Common* cm = picture_control_set_ptr->parent_pcs_ptr->av1_cm;
    EbReferenceObject_t  * refObjL0, *refObjL1;
    memset(cm->sg_frame_ep_cnt, 0, SGRPROJ_PARAMS * sizeof(int32_t));
    cm->sg_frame_ep = 0;

    // NADER: set cm->sg_ref_frame_ep[0] = cm->sg_ref_frame_ep[1] = -1 to perform all iterations
    switch(picture_control_set_ptr->slice_type){
    case I_SLICE:
        cm->sg_ref_frame_ep[0] = cm->sg_ref_frame_ep[1] = -1;
        break;
    case B_SLICE:
        refObjL0 = (EbReferenceObject_t*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_0]->objectPtr;
        refObjL1 = (EbReferenceObject_t*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_1]->objectPtr;
        cm->sg_ref_frame_ep[0] = refObjL0->sg_frame_ep;
        cm->sg_ref_frame_ep[1] = refObjL1->sg_frame_ep;
        break;
    case P_SLICE:
        refObjL0 = (EbReferenceObject_t*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_0]->objectPtr;
        cm->sg_ref_frame_ep[0] = refObjL0->sg_frame_ep;
        break;
    default:
        printf("SG: Not supported picture type");
        break;
    }
}
#endif
#if FAST_CDEF
/******************************************************
* Set the reference cdef strength for a given picture
******************************************************/
void set_reference_cdef_strength(
    PictureControlSet_t                    *picture_control_set_ptr)
{
    EbReferenceObject_t  * refObjL0, *refObjL1;
    int32_t strength;
    // NADER: set picture_control_set_ptr->parent_pcs_ptr->use_ref_frame_cdef_strength 0 to test all strengths
    switch (picture_control_set_ptr->slice_type) {
    case I_SLICE:
        picture_control_set_ptr->parent_pcs_ptr->use_ref_frame_cdef_strength = 0;
        picture_control_set_ptr->parent_pcs_ptr->cdf_ref_frame_strenght = 0;
        break;
    case B_SLICE:
        refObjL0 = (EbReferenceObject_t*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_0]->objectPtr;
        refObjL1 = (EbReferenceObject_t*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_1]->objectPtr;
        strength = (refObjL0->cdef_frame_strength + refObjL1->cdef_frame_strength) / 2;
        picture_control_set_ptr->parent_pcs_ptr->use_ref_frame_cdef_strength = 1;
        picture_control_set_ptr->parent_pcs_ptr->cdf_ref_frame_strenght = strength;
        break;
    case P_SLICE:
        refObjL0 = (EbReferenceObject_t*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_0]->objectPtr;
        strength = refObjL0->cdef_frame_strength;
        picture_control_set_ptr->parent_pcs_ptr->use_ref_frame_cdef_strength = 1;
        picture_control_set_ptr->parent_pcs_ptr->cdf_ref_frame_strenght = strength;
        break;
    default:
        printf("CDEF: Not supported picture type");
        break;
    }
}
#endif
/******************************************************
* Compute Tc, and Beta offsets for a given picture
******************************************************/
void AdaptiveDlfParameterComputation(
    ModeDecisionConfigurationContext_t     *context_ptr,
    SequenceControlSet_t                   *sequence_control_set_ptr,
    PictureControlSet_t                    *picture_control_set_ptr,
    EbBool                                    scene_transition_flag)
{
    EbReferenceObject_t  * refObjL0, *refObjL1;
    uint8_t                  highIntra = 0;
    picture_control_set_ptr->high_intra_slection = highIntra;

    (void)sequence_control_set_ptr;
    (void)context_ptr;



    if (picture_control_set_ptr->slice_type == B_SLICE) {

        refObjL0 = (EbReferenceObject_t*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_0]->objectPtr;
        refObjL1 = (EbReferenceObject_t*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_1]->objectPtr;

        if ((refObjL0->intra_coded_area > INTRA_AREA_TH_CLASS_1[picture_control_set_ptr->parent_pcs_ptr->hierarchical_levels][refObjL0->tmpLayerIdx]) && (refObjL1->intra_coded_area > INTRA_AREA_TH_CLASS_1[picture_control_set_ptr->parent_pcs_ptr->hierarchical_levels][refObjL1->tmpLayerIdx]))

            highIntra = 1;
        else
            highIntra = 0;
    }

    if (picture_control_set_ptr->slice_type == B_SLICE) {

        if (!scene_transition_flag && !picture_control_set_ptr->parent_pcs_ptr->fade_out_from_black && !picture_control_set_ptr->parent_pcs_ptr->fade_in_to_black) {
            refObjL0 = (EbReferenceObject_t*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_0]->objectPtr;
            refObjL1 = (EbReferenceObject_t*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_1]->objectPtr;

            if (picture_control_set_ptr->temporal_layer_index == 0) {
                highIntra = (picture_control_set_ptr->parent_pcs_ptr->intra_coded_block_probability > INTRA_AREA_TH_CLASS_1[picture_control_set_ptr->parent_pcs_ptr->hierarchical_levels][picture_control_set_ptr->temporal_layer_index]) ? 1 : 0;
            }
            else {
                highIntra = (refObjL0->penalizeSkipflag || refObjL1->penalizeSkipflag) ? 1 : 0;
            }
        }
    }

    picture_control_set_ptr->high_intra_slection = highIntra;
}


/******************************************************
 * Mode Decision Configuration Context Constructor
 ******************************************************/
EbErrorType ModeDecisionConfigurationContextCtor(
    ModeDecisionConfigurationContext_t **context_dbl_ptr,
    EbFifo_t                            *rateControlInputFifoPtr,

    EbFifo_t                            *modeDecisionConfigurationOutputFifoPtr,
    uint16_t                                 sb_total_count)

{
    ModeDecisionConfigurationContext_t *context_ptr;

    EB_MALLOC(ModeDecisionConfigurationContext_t*, context_ptr, sizeof(ModeDecisionConfigurationContext_t), EB_N_PTR);

    *context_dbl_ptr = context_ptr;

    // Input/Output System Resource Manager FIFOs
    context_ptr->rateControlInputFifoPtr = rateControlInputFifoPtr;
    context_ptr->modeDecisionConfigurationOutputFifoPtr = modeDecisionConfigurationOutputFifoPtr;
    // Rate estimation
    EB_MALLOC(MdRateEstimationContext_t*, context_ptr->md_rate_estimation_ptr, sizeof(MdRateEstimationContext_t), EB_N_PTR);


    // Budgeting
    EB_MALLOC(uint32_t*, context_ptr->lcuScoreArray, sizeof(uint32_t) * sb_total_count, EB_N_PTR);
    EB_MALLOC(uint8_t *, context_ptr->lcuCostArray, sizeof(uint8_t) * sb_total_count, EB_N_PTR);


    return EB_ErrorNone;
}

/******************************************************
* Predict the SB partitionning
******************************************************/
void PerformEarlyLcuPartitionning(
    ModeDecisionConfigurationContext_t     *context_ptr,
    SequenceControlSet_t                   *sequence_control_set_ptr,
    PictureControlSet_t                    *picture_control_set_ptr) {

    LargestCodingUnit_t            *sb_ptr;
    uint32_t                         sb_index;
    uint32_t                         slice_type;

    // MD Conf Rate Estimation Array from encodeContext
    MdRateEstimationContext_t    *mdConfRateEstimationArray;

    // Lambda Assignement
    if (sequence_control_set_ptr->static_config.pred_structure == EB_PRED_RANDOM_ACCESS) {

        if (picture_control_set_ptr->temporal_layer_index == 0) {

            context_ptr->lambda = lambdaModeDecisionRaSad[context_ptr->qp];
        }
        else if (picture_control_set_ptr->temporal_layer_index < 3) {
            context_ptr->lambda = lambdaModeDecisionRaSadQpScalingL1[context_ptr->qp];
        }
        else {
            context_ptr->lambda = lambdaModeDecisionRaSadQpScalingL3[context_ptr->qp];
        }
    }
    else {
        if (picture_control_set_ptr->temporal_layer_index == 0) {
            context_ptr->lambda = lambdaModeDecisionLdSad[context_ptr->qp];
        }
        else {
            context_ptr->lambda = lambdaModeDecisionLdSadQpScaling[context_ptr->qp];
        }
    }

#if NEW_QPS
    context_ptr->qp_index = (uint8_t)picture_control_set_ptr->parent_pcs_ptr->base_qindex;
#else
    context_ptr->qp_index = quantizer_to_qindex[context_ptr->qp];
#endif
    uint32_t lambdaSse;
    uint32_t lambdaSad;
    (*av1_lambda_assignment_function_table[picture_control_set_ptr->parent_pcs_ptr->pred_structure])(
        &lambdaSad,
        &lambdaSse,
        &lambdaSad,
        &lambdaSse,
        (uint8_t)picture_control_set_ptr->parent_pcs_ptr->enhanced_picture_ptr->bit_depth,
        context_ptr->qp_index);
    context_ptr->lambda = (uint64_t)lambdaSad;

    // Slice Type
    slice_type =
        (picture_control_set_ptr->parent_pcs_ptr->idr_flag == EB_TRUE) ? I_SLICE :
        picture_control_set_ptr->slice_type;

    // Increment the MD Rate Estimation array pointer to point to the right address based on the QP and slice type
    mdConfRateEstimationArray = (MdRateEstimationContext_t*)sequence_control_set_ptr->encode_context_ptr->md_rate_estimation_array;
    mdConfRateEstimationArray += slice_type * TOTAL_NUMBER_OF_QP_VALUES + context_ptr->qp;

    // Reset MD rate Estimation table to initial values by copying from md_rate_estimation_array
    EB_MEMCPY(&(context_ptr->md_rate_estimation_ptr->splitFlagBits[0]), &(mdConfRateEstimationArray->splitFlagBits[0]), sizeof(MdRateEstimationContext_t));

    picture_control_set_ptr->parent_pcs_ptr->average_qp = (uint8_t)picture_control_set_ptr->parent_pcs_ptr->picture_qp;

    // SB Loop : Partitionnig Decision
    for (sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; ++sb_index) {

        sb_ptr = picture_control_set_ptr->sb_ptr_array[sb_index];
        {
            sb_ptr->qp = (uint8_t)picture_control_set_ptr->parent_pcs_ptr->picture_qp;
        }

        EarlyModeDecisionLcu(
            sequence_control_set_ptr,
            picture_control_set_ptr,
            sb_ptr,
            sb_index,
            context_ptr);

    } // End of SB Loop

}

void PerformEarlyLcuPartitionningLcu(
    ModeDecisionConfigurationContext_t     *context_ptr,
    SequenceControlSet_t                   *sequence_control_set_ptr,
    PictureControlSet_t                    *picture_control_set_ptr,
    uint32_t                                    sb_index) {

    LargestCodingUnit_t            *sb_ptr;

    // SB Loop : Partitionnig Decision
    sb_ptr = picture_control_set_ptr->sb_ptr_array[sb_index];
    {
        sb_ptr->qp = (uint8_t)picture_control_set_ptr->parent_pcs_ptr->picture_qp;
    }

    EarlyModeDecisionLcu(
        sequence_control_set_ptr,
        picture_control_set_ptr,
        sb_ptr,
        sb_index,
        context_ptr);
}


void Forward85CuToModeDecisionLCU(
    SequenceControlSet_t  *sequence_control_set_ptr,
    PictureControlSet_t   *picture_control_set_ptr,
    uint32_t                 sb_index) {

    const CodedUnitStats_t  *cuStatsPtr;
    EbBool split_flag;
    // SB Loop : Partitionnig Decision

    SbParams_t  *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];
    MdcLcuData_t *resultsPtr = &picture_control_set_ptr->mdc_sb_array[sb_index];
    uint32_t cuIndexInRaterScan;   uint16_t cuVar;

    resultsPtr->leaf_count = 0;
    uint8_t cu_index = 0;
    while (cu_index < CU_MAX_COUNT)
    {
        split_flag = EB_TRUE;
        cuStatsPtr = GetCodedUnitStats(cu_index);
        if (sb_params->raster_scan_cu_validity[MD_SCAN_TO_RASTER_SCAN[cu_index]])
        {
            switch (cuStatsPtr->depth) {

            case 0:

                resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = cu_index;
                resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = split_flag = EB_TRUE;

                break;
            case 1:
                cuIndexInRaterScan = MD_SCAN_TO_RASTER_SCAN[cu_index];
                cuVar = (picture_control_set_ptr->parent_pcs_ptr->variance[sb_index][cuIndexInRaterScan]);
#if ENCODER_MODE_CLEANUP
                if (picture_control_set_ptr->slice_type == I_SLICE && cuVar > 40)

#else
                if (picture_control_set_ptr->enc_mode <= ENC_M2 && picture_control_set_ptr->slice_type == I_SLICE && cuVar > 40)
#endif
                    split_flag = EB_TRUE;
                else {
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = cu_index;
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = split_flag = EB_TRUE;
                }

                break;

            case 2:

                resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = cu_index;
                resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = split_flag = EB_TRUE;

                break;
            case 3:

                resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = cu_index;
                resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = split_flag = EB_FALSE;

                break;

            default:
                resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = cu_index;
                resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = split_flag = EB_TRUE;
                break;
            }
        }

        cu_index += (split_flag == EB_TRUE) ? 1 : DepthOffset[cuStatsPtr->depth];

    } // End CU Loop
}

void Forward84CuToModeDecisionLCU(
    SequenceControlSet_t  *sequence_control_set_ptr,
    PictureControlSet_t   *picture_control_set_ptr,
    uint32_t                 sb_index) {

    const CodedUnitStats_t  *cuStatsPtr;
    EbBool split_flag;
    // SB Loop : Partitionnig Decision

    SbParams_t  *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];
    MdcLcuData_t *resultsPtr = &picture_control_set_ptr->mdc_sb_array[sb_index];
    uint32_t cuIndexInRaterScan;   uint16_t cuVar;

    resultsPtr->leaf_count = 0;
    uint8_t cu_index = 0;
    while (cu_index < CU_MAX_COUNT)
    {
        split_flag = EB_TRUE;
        cuStatsPtr = GetCodedUnitStats(cu_index);
        if (sb_params->raster_scan_cu_validity[MD_SCAN_TO_RASTER_SCAN[cu_index]])
        {
            switch (cuStatsPtr->depth) {

            case 0:

                split_flag = EB_TRUE;

                break;
            case 1:
                cuIndexInRaterScan = MD_SCAN_TO_RASTER_SCAN[cu_index];
                cuVar = (picture_control_set_ptr->parent_pcs_ptr->variance[sb_index][cuIndexInRaterScan]);
#if ENCODER_MODE_CLEANUP
                if (picture_control_set_ptr->slice_type == I_SLICE && cuVar > 40)

#else
                if (picture_control_set_ptr->enc_mode <= ENC_M2 && picture_control_set_ptr->slice_type == I_SLICE && cuVar > 40)
#endif
                    split_flag = EB_TRUE;
                else {
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = cu_index;
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = split_flag = EB_TRUE;
                }

                break;

            case 2:

                resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = cu_index;
                resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = split_flag = EB_TRUE;

                break;
            case 3:

                resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = cu_index;
                resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = split_flag = EB_FALSE;

                break;

            default:
                resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = cu_index;
                resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = split_flag = EB_TRUE;
                break;
            }
        }

        cu_index += (split_flag == EB_TRUE) ? 1 : DepthOffset[cuStatsPtr->depth];

    } // End CU Loop
}

void forward_all_blocks_to_md(
    SequenceControlSet_t                   *sequence_control_set_ptr,
    PictureControlSet_t                    *picture_control_set_ptr)
{

    uint32_t                   sb_index;
    EbBool                  split_flag;

    UNUSED(split_flag);

    for (sb_index = 0; sb_index < sequence_control_set_ptr->sb_tot_cnt; ++sb_index)
    {

        MdcLcuData_t *resultsPtr = &picture_control_set_ptr->mdc_sb_array[sb_index];

        resultsPtr->leaf_count = 0;

        uint32_t  blk_index = 0;



        while (blk_index < sequence_control_set_ptr->max_block_cnt)
        {
            split_flag = EB_TRUE;

            const BlockGeom * blk_geom = Get_blk_geom_mds(blk_index);

            //if the parentSq is inside inject this block
            uint8_t is_blk_allowed = picture_control_set_ptr->slice_type != I_SLICE ? 1 : (blk_geom->sq_size < 128) ? 1 : 0;

            if (sequence_control_set_ptr->sb_geom[sb_index].block_is_inside_md_scan[blk_index] && is_blk_allowed)

            {
                resultsPtr->leaf_data_array[resultsPtr->leaf_count].tot_d1_blocks =

                    blk_geom->sq_size == 128 ? 17 :
                    blk_geom->sq_size > 8 ? 25 :
                    blk_geom->sq_size == 8 ? 5 : 1;

                resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = 0;//valid only for square 85 world. will be removed.
                resultsPtr->leaf_data_array[resultsPtr->leaf_count].mds_idx = blk_index;
                if (blk_geom->sq_size > 4)
                {
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = EB_TRUE;
                    split_flag = EB_TRUE;
                }
                else {
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = EB_FALSE;
                    split_flag = EB_FALSE;
                }


            }

            blk_index++;

        }

    }

    picture_control_set_ptr->parent_pcs_ptr->average_qp = (uint8_t)picture_control_set_ptr->parent_pcs_ptr->picture_qp;

}

#if INJECT_ONLY_SQ
void forward_sq_blocks_to_md(
    SequenceControlSet_t                   *sequence_control_set_ptr,
    PictureControlSet_t                    *picture_control_set_ptr)
{

    uint32_t                   sb_index;
    EbBool                  split_flag;


    for (sb_index = 0; sb_index < sequence_control_set_ptr->sb_tot_cnt; ++sb_index)
    {

        MdcLcuData_t *resultsPtr = &picture_control_set_ptr->mdc_sb_array[sb_index];

        resultsPtr->leaf_count = 0;

        uint32_t  blk_index = picture_control_set_ptr->slice_type == I_SLICE && sequence_control_set_ptr->sb_size == BLOCK_128X128 ? 17 : 0;



        while (blk_index < sequence_control_set_ptr->max_block_cnt)
        {
            split_flag = EB_TRUE;

            const BlockGeom * blk_geom = Get_blk_geom_mds(blk_index);

            //if the parentSq is inside inject this block
            if (sequence_control_set_ptr->sb_geom[sb_index].block_is_inside_md_scan[blk_index])

            {

                //int32_t offset_d1 = ns_blk_offset[(int32_t)from_shape_to_part[blk_geom->shape]]; //cu_ptr->best_d1_blk; // TOCKECK
                //int32_t num_d1_block = ns_blk_num[(int32_t)from_shape_to_part[blk_geom->shape]]; // context_ptr->blk_geom->totns; // TOCKECK
                //
                //                                                  // for (int32_t d1_itr = blk_it; d1_itr < blk_it + num_d1_block; d1_itr++) {
                // for (int32_t d1_itr = (int32_t)blk_index ; d1_itr < (int32_t)blk_index +  num_d1_block ; d1_itr++) {

                resultsPtr->leaf_data_array[resultsPtr->leaf_count].tot_d1_blocks = 1;


                resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = 0;//valid only for square 85 world. will be removed.
                resultsPtr->leaf_data_array[resultsPtr->leaf_count].mds_idx = blk_index;

                if (blk_geom->sq_size > 4)
                {
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = EB_TRUE;
                    split_flag = EB_TRUE;
                }
                else {
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = EB_FALSE;
                    split_flag = EB_FALSE;
                }


            }

            blk_index += split_flag ? d1_depth_offset[sequence_control_set_ptr->sb_size == BLOCK_128X128][blk_geom->depth] : ns_depth_offset[sequence_control_set_ptr->sb_size == BLOCK_128X128][blk_geom->depth];

        }

    }

    picture_control_set_ptr->parent_pcs_ptr->average_qp = (uint8_t)picture_control_set_ptr->parent_pcs_ptr->picture_qp;

}
#endif
void Forward85CuToModeDecision(
    SequenceControlSet_t                   *sequence_control_set_ptr,
    PictureControlSet_t                    *picture_control_set_ptr) {

    const CodedUnitStats_t  *cuStatsPtr;
    uint32_t                   sb_index;
    EbBool split_flag;
    // SB Loop : Partitionnig Decision
    for (sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; ++sb_index) {

        SbParams_t  *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];
        MdcLcuData_t *resultsPtr = &picture_control_set_ptr->mdc_sb_array[sb_index];
#if !REMOVE_INTRA_CONST
        uint32_t cuIndexInRaterScan;   uint16_t cuVar;
#endif
        resultsPtr->leaf_count = 0;
        uint8_t cu_index = 0;
        while (cu_index < CU_MAX_COUNT)
        {
            split_flag = EB_TRUE;
            cuStatsPtr = GetCodedUnitStats(cu_index);
            if (sb_params->raster_scan_cu_validity[MD_SCAN_TO_RASTER_SCAN[cu_index]])
            {
                switch (cuStatsPtr->depth) {

                case 0:
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = cu_index;
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = split_flag = EB_TRUE;

                    break;
                case 1:
#if REMOVE_INTRA_CONST
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = cu_index;
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = split_flag = EB_TRUE;
#else
                    //OMK To revisit : add Varpart flag and move to MD
                    cuIndexInRaterScan = MD_SCAN_TO_RASTER_SCAN[cu_index];
                    cuVar = (picture_control_set_ptr->parent_pcs_ptr->variance[sb_index][cuIndexInRaterScan]);
#if ENCODER_MODE_CLEANUP
                    if ((picture_control_set_ptr->slice_type == I_SLICE && cuVar > 40) || (sequence_control_set_ptr->input_resolution < INPUT_SIZE_4K_RANGE&& picture_control_set_ptr->slice_type == I_SLICE && cuVar>40))

#else
                    if ((picture_control_set_ptr->enc_mode <= ENC_M2 && picture_control_set_ptr->slice_type == I_SLICE && cuVar > 40) || (sequence_control_set_ptr->input_resolution < INPUT_SIZE_4K_RANGE&& picture_control_set_ptr->slice_type == I_SLICE && cuVar>40))
#endif
                        split_flag = EB_TRUE;
                    else {
                        resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = cu_index;
                        resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = split_flag = EB_TRUE;
                    }
#endif
                    break;

                case 2:

                    resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = cu_index;
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = split_flag = EB_TRUE;

                    break;
                case 3:

                    resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = cu_index;
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = split_flag = EB_FALSE;

                    break;

                default:
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = cu_index;
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = split_flag = EB_TRUE;
                    break;
                }
            }

            cu_index += (split_flag == EB_TRUE) ? 1 : DepthOffset[cuStatsPtr->depth];

        } // End CU Loop
    }

    picture_control_set_ptr->parent_pcs_ptr->average_qp = (uint8_t)picture_control_set_ptr->parent_pcs_ptr->picture_qp;

}

void Forward84CuToModeDecision(
    SequenceControlSet_t                   *sequence_control_set_ptr,
    PictureControlSet_t                    *picture_control_set_ptr) {

    const CodedUnitStats_t  *cuStatsPtr;
    uint32_t                   sb_index;
    EbBool split_flag;
    // SB Loop : Partitionnig Decision
    for (sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; ++sb_index) {

        SbParams_t  *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];
        MdcLcuData_t *resultsPtr = &picture_control_set_ptr->mdc_sb_array[sb_index];
        uint32_t cuIndexInRaterScan;   uint16_t cuVar;

        resultsPtr->leaf_count = 0;
        uint8_t cu_index = 0;
        while (cu_index < CU_MAX_COUNT)
        {
            split_flag = EB_TRUE;
            cuStatsPtr = GetCodedUnitStats(cu_index);
            if (sb_params->raster_scan_cu_validity[MD_SCAN_TO_RASTER_SCAN[cu_index]])
            {
                switch (cuStatsPtr->depth) {

                case 0:
                    if (picture_control_set_ptr->slice_type == I_SLICE) {
                        split_flag = EB_TRUE;
                    }
                    else {
                        resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = cu_index;
                        resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = split_flag = EB_TRUE;
                    }

                    break;
                case 1:

                    //OMK To revisit : add Varpart flag and move to MD
                    cuIndexInRaterScan = MD_SCAN_TO_RASTER_SCAN[cu_index];
                    cuVar = (picture_control_set_ptr->parent_pcs_ptr->variance[sb_index][cuIndexInRaterScan]);
#if ENCODER_MODE_CLEANUP
                    if ((picture_control_set_ptr->slice_type == I_SLICE && cuVar > 40) || (sequence_control_set_ptr->input_resolution < INPUT_SIZE_4K_RANGE&& picture_control_set_ptr->slice_type == I_SLICE && cuVar>40))

#else
                    if ((picture_control_set_ptr->enc_mode <= ENC_M2 && picture_control_set_ptr->slice_type == I_SLICE && cuVar > 40) || (sequence_control_set_ptr->input_resolution < INPUT_SIZE_4K_RANGE&& picture_control_set_ptr->slice_type == I_SLICE && cuVar>40))
#endif
                        split_flag = EB_TRUE;
                    else {
                        resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = cu_index;
                        resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = split_flag = EB_TRUE;
                    }

                    break;

                case 2:

                    resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = cu_index;
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = split_flag = EB_TRUE;

                    break;
                case 3:



                    resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = cu_index;
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = split_flag = EB_FALSE;

                    break;

                default:
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = cu_index;
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = split_flag = EB_TRUE;
                    break;
                }
            }

            cu_index += (split_flag == EB_TRUE) ? 1 : DepthOffset[cuStatsPtr->depth];

        } // End CU Loop
    }

    picture_control_set_ptr->parent_pcs_ptr->average_qp = (uint8_t)picture_control_set_ptr->parent_pcs_ptr->picture_qp;

}

void PartitioningInitialization(
    SequenceControlSet_t                   *sequence_control_set_ptr,
    PictureControlSet_t                    *picture_control_set_ptr,
    ModeDecisionConfigurationContext_t     *context_ptr) {

    uint32_t                         slice_type;

    // HG: to collapse
    // MD Conf Rate Estimation Array from encodeContext
    MdRateEstimationContext_t    *mdConfRateEstimationArray;

    // Lambda Assignement
    if (sequence_control_set_ptr->static_config.pred_structure == EB_PRED_RANDOM_ACCESS) {

        if (picture_control_set_ptr->temporal_layer_index == 0) {

            context_ptr->lambda = lambdaModeDecisionRaSad[context_ptr->qp];
        }
        else if (picture_control_set_ptr->temporal_layer_index < 3) {
            context_ptr->lambda = lambdaModeDecisionRaSadQpScalingL1[context_ptr->qp];
        }
        else {
            context_ptr->lambda = lambdaModeDecisionRaSadQpScalingL3[context_ptr->qp];
        }
    }
    else {
        if (picture_control_set_ptr->temporal_layer_index == 0) {
            context_ptr->lambda = lambdaModeDecisionLdSad[context_ptr->qp];
        }
        else {
            context_ptr->lambda = lambdaModeDecisionLdSadQpScaling[context_ptr->qp];
        }
    }

    // Slice Type
    slice_type =
        (picture_control_set_ptr->parent_pcs_ptr->idr_flag == EB_TRUE) ? I_SLICE :
        picture_control_set_ptr->slice_type;

    // Increment the MD Rate Estimation array pointer to point to the right address based on the QP and slice type
    mdConfRateEstimationArray = (MdRateEstimationContext_t*)sequence_control_set_ptr->encode_context_ptr->md_rate_estimation_array;
    mdConfRateEstimationArray += slice_type * TOTAL_NUMBER_OF_QP_VALUES + context_ptr->qp;

    // Reset MD rate Estimation table to initial values by copying from md_rate_estimation_array
    EB_MEMCPY(&(context_ptr->md_rate_estimation_ptr->splitFlagBits[0]), &(mdConfRateEstimationArray->splitFlagBits[0]), sizeof(MdRateEstimationContext_t)/*sizeof(EB_BitFraction)* MAX_SIZE_OF_MD_RATE_ESTIMATION_CASES*/);

    picture_control_set_ptr->parent_pcs_ptr->average_qp = (uint8_t)picture_control_set_ptr->parent_pcs_ptr->picture_qp;
}


/******************************************************
* Derive MD parameters
******************************************************/
void SetMdSettings(
    SequenceControlSet_t                   *sequence_control_set_ptr,
    PictureControlSet_t                    *picture_control_set_ptr) {
    // Initialize the homogeneous area threshold
    // Set the MD Open Loop Flag
    // HG - to clean up the intra_md_open_loop_flag derivation
#if ENCODER_MODE_CLEANUP
    if (0)
#else
    if (picture_control_set_ptr->enc_mode >= ENC_M3)
#endif
        picture_control_set_ptr->intra_md_open_loop_flag = picture_control_set_ptr->slice_type == I_SLICE ? EB_FALSE : EB_TRUE;
    else
        picture_control_set_ptr->intra_md_open_loop_flag = picture_control_set_ptr->temporal_layer_index == 0 ? EB_FALSE : EB_TRUE;
#if ENCODER_MODE_CLEANUP
    if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE && sequence_control_set_ptr->input_resolution < INPUT_SIZE_4K_RANGE)
#else
    if (picture_control_set_ptr->enc_mode <= ENC_M1 && picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE && sequence_control_set_ptr->input_resolution < INPUT_SIZE_4K_RANGE)
#endif
        picture_control_set_ptr->intra_md_open_loop_flag = EB_FALSE;

    picture_control_set_ptr->constrained_intra_flag = EB_FALSE;

    if (sequence_control_set_ptr->static_config.constrained_intra == EB_TRUE &&
        picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_FALSE)
    {
        picture_control_set_ptr->constrained_intra_flag = EB_TRUE;
    }
#if ! ENCODER_MODE_CLEANUP
    picture_control_set_ptr->limit_intra = EB_FALSE;
    if (picture_control_set_ptr->enc_mode >= ENC_M5 && picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_FALSE)
        picture_control_set_ptr->limit_intra = EB_TRUE;
#endif
    picture_control_set_ptr->limit_intra = EB_FALSE;

    picture_control_set_ptr->intra_md_open_loop_flag = EB_FALSE;
}

/******************************************************
* Detect complex/non-flat/moving SB in a non-complex area (used to refine MDC depth control in Gold)
******************************************************/
void DetectComplexNonFlatMovingLcu(
    PictureControlSet_t       *picture_control_set_ptr,
    uint32_t                    picture_width_in_sb,
    uint32_t                    picture_height_in_sb) {


    LargestCodingUnit_t *sb_ptr;
    uint32_t               sb_index;
    uint32_t                 sb_x;
    uint32_t                 sb_y;

    if (picture_control_set_ptr->parent_pcs_ptr->non_moving_index_average >= 10 && picture_control_set_ptr->temporal_layer_index <= 2) {
        // Determine deltaQP and assign QP value for each leaf
        for (sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; ++sb_index) {

            sb_ptr = picture_control_set_ptr->sb_ptr_array[sb_index];

            sb_x = sb_index % picture_width_in_sb;
            sb_y = sb_index / picture_width_in_sb;
            EbBool condition = EB_FALSE;

            if (!picture_control_set_ptr->parent_pcs_ptr->similar_colocated_sb_array[sb_index] &&
                sb_ptr->picture_control_set_ptr->parent_pcs_ptr->edge_results_ptr[sb_index].edge_block_num > 0) {
                condition = EB_TRUE;
            }

            if (condition) {
                uint32_t  counter = 0;
                if (sb_x > 0 && sb_x < (uint32_t)(picture_width_in_sb - 1) && sb_y >  0 && sb_y < (uint32_t)(picture_height_in_sb - 1)) {
                    // Top
                    if (picture_control_set_ptr->parent_pcs_ptr->edge_results_ptr[sb_index - picture_width_in_sb].edge_block_num == 0)
                        counter++;
                    // Bottom
                    if (picture_control_set_ptr->parent_pcs_ptr->edge_results_ptr[sb_index + picture_width_in_sb].edge_block_num == 0)
                        counter++;
                    // Left
                    if (picture_control_set_ptr->parent_pcs_ptr->edge_results_ptr[sb_index - 1].edge_block_num == 0)
                        counter++;
                    // right
                    if (picture_control_set_ptr->parent_pcs_ptr->edge_results_ptr[sb_index + 1].edge_block_num == 0)
                        counter++;
                }
            }
        }
    }
}

EbAuraStatus AuraDetection64x64Gold(
    PictureControlSet_t           *picture_control_set_ptr,
    uint8_t                          picture_qp,
    uint32_t                         xLcuIndex,
    uint32_t                         yLcuIndex
)
{

    SequenceControlSet_t  *sequence_control_set_ptr = (SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->objectPtr;
    int32_t                 picture_width_in_sb = (sequence_control_set_ptr->luma_width + BLOCK_SIZE_64 - 1) >> LOG2_64_SIZE;
    int32_t                 picture_height_in_sb = (sequence_control_set_ptr->luma_height + BLOCK_SIZE_64 - 1) >> LOG2_64_SIZE;
    uint32_t                 sb_index = yLcuIndex * picture_width_in_sb + xLcuIndex;
    uint32_t                 currDist;
    uint32_t                 topDist, topLDist, topRDist;
    uint32_t                 localAvgDist, distThresh0, distThresh1;
    int32_t                 lcuOffset;

    EbAuraStatus            auraClass = AURA_STATUS_0;
    uint8_t                    auraClass1 = 0;
    uint8_t                    auraClass2 = 0;

    int32_t                 xMv0 = 0;
    int32_t                 yMv0 = 0;
    int32_t                 xMv1 = 0;
    int32_t                 yMv1 = 0;

    uint32_t                 leftDist, rightDist;


    distThresh0 = picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag || sequence_control_set_ptr->input_resolution == INPUT_SIZE_4K_RANGE ? 15 : 14;
    distThresh1 = 23;


    if (picture_qp > 38) {
        distThresh0 = distThresh0 << 2;
        distThresh1 = distThresh1 << 2;
    }

    if (xLcuIndex > 0 && xLcuIndex < (uint32_t)(picture_width_in_sb - 1) && yLcuIndex>0 && yLcuIndex < (uint32_t)(picture_height_in_sb - 1)) {

        uint32_t k;



        MeCuResults_t * mePuResult = &picture_control_set_ptr->parent_pcs_ptr->me_results[sb_index][0];

        //Curr Block

        for (k = 0; k < mePuResult->totalMeCandidateIndex; k++) {

            if (mePuResult->distortionDirection[k].direction == UNI_PRED_LIST_0) {
                // Get reference list 0 / reference index 0 MV
                xMv0 = mePuResult->xMvL0;
                yMv0 = mePuResult->yMvL0;
            }
            if (mePuResult->distortionDirection[k].direction == UNI_PRED_LIST_1) {
                // Get reference list  1 / reference index 0 MV
                xMv1 = mePuResult->xMvL1;
                yMv1 = mePuResult->yMvL1;
            }

        }
        currDist = mePuResult->distortionDirection[0].distortion;



        if ((currDist > 64 * 64) &&
            // Only mark a block as aura when it is moving (MV amplitude higher than X; X is temporal layer dependent)
            (abs(xMv0) > GLOBAL_MOTION_THRESHOLD[picture_control_set_ptr->parent_pcs_ptr->hierarchical_levels][picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index] ||
                abs(yMv0) > GLOBAL_MOTION_THRESHOLD[picture_control_set_ptr->parent_pcs_ptr->hierarchical_levels][picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index] ||
                abs(xMv1) > GLOBAL_MOTION_THRESHOLD[picture_control_set_ptr->parent_pcs_ptr->hierarchical_levels][picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index] ||
                abs(yMv1) > GLOBAL_MOTION_THRESHOLD[picture_control_set_ptr->parent_pcs_ptr->hierarchical_levels][picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index]))
        {




            //Top Distortion
            lcuOffset = -picture_width_in_sb;
            topDist = picture_control_set_ptr->parent_pcs_ptr->me_results[sb_index + lcuOffset]->distortionDirection[0].distortion;


            //TopLeft Distortion
            lcuOffset = -picture_width_in_sb - 1;
            topLDist = picture_control_set_ptr->parent_pcs_ptr->me_results[sb_index + lcuOffset]->distortionDirection[0].distortion;


            //TopRightDistortion
            lcuOffset = -picture_width_in_sb + 1;
            topRDist = picture_control_set_ptr->parent_pcs_ptr->me_results[sb_index + lcuOffset]->distortionDirection[0].distortion;


            topRDist = (xLcuIndex < (uint32_t)(picture_width_in_sb - 2)) ? topRDist : currDist;

            //left Distortion
            lcuOffset = -1;
            leftDist = picture_control_set_ptr->parent_pcs_ptr->me_results[sb_index + lcuOffset]->distortionDirection[0].distortion;



            //RightDistortion
            lcuOffset = 1;
            rightDist = picture_control_set_ptr->parent_pcs_ptr->me_results[sb_index + lcuOffset]->distortionDirection[0].distortion;



            rightDist = (xLcuIndex < (uint32_t)(picture_width_in_sb - 2)) ? topRDist : currDist;

            localAvgDist = MIN(MIN(MIN(topLDist, MIN(topRDist, topDist)), leftDist), rightDist);

            if (10 * currDist > distThresh0*localAvgDist) {
                if (10 * currDist > distThresh1*localAvgDist)
                    auraClass2++;
                else
                    auraClass1++;
            }
        }

    }

    auraClass = (auraClass2 > 0 || auraClass1 > 0) ? AURA_STATUS_1 : AURA_STATUS_0;

    return   auraClass;

}


/******************************************************
* Aura detection
******************************************************/
void AuraDetection(
    SequenceControlSet_t         *sequence_control_set_ptr,
    PictureControlSet_t            *picture_control_set_ptr,
    uint32_t                          picture_width_in_sb,
    uint32_t                          picture_height_in_sb)

{
    uint32_t                      sb_index;
    uint32_t                      sb_x;
    uint32_t                      sb_y;

    for (sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; ++sb_index) {

        SbParams_t    *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];
        LargestCodingUnit_t*        sb_ptr = picture_control_set_ptr->sb_ptr_array[sb_index];
        sb_x = sb_params->horizontal_index;
        sb_y = sb_params->vertical_index;
        // Aura status intialization
        sb_ptr->aura_status_iii = INVALID_AURA_STATUS;

        if (picture_control_set_ptr->slice_type == B_SLICE) {
            if ((sb_x > 0) && (sb_x < picture_width_in_sb - 1) && (sb_y < picture_height_in_sb - 1)) {
                sb_ptr->aura_status_iii = AuraDetection64x64Gold(
                    picture_control_set_ptr,
                    (uint8_t)picture_control_set_ptr->picture_qp,
                    sb_x,
                    sb_y);
            }
        }


    }
    return;
}

EbErrorType DeriveDefaultSegments(
    PictureControlSet_t                 *picture_control_set_ptr,
    ModeDecisionConfigurationContext_t  *context_ptr)
{
    EbErrorType                    return_error = EB_ErrorNone;

    // @ BASE MDC is not considered as similar to BDP_L in term of speed
    if (picture_control_set_ptr->temporal_layer_index == 0) {

        if (context_ptr->depthSensitivePictureFlag && context_ptr->budget >= (uint32_t)(picture_control_set_ptr->parent_pcs_ptr->sb_total_count * LIGHT_BDP_COST)) {

            if (context_ptr->budget > (uint32_t)(picture_control_set_ptr->parent_pcs_ptr->sb_total_count * BDP_COST)) {
                context_ptr->numberOfSegments = 2;
                context_ptr->scoreTh[0] = (int8_t)((1 * 100) / context_ptr->numberOfSegments);
                context_ptr->intervalCost[0] = context_ptr->costDepthMode[LCU_BDP_DEPTH_MODE - 1];
                context_ptr->intervalCost[1] = context_ptr->costDepthMode[LCU_FULL84_DEPTH_MODE - 1];
            }
            else {
                context_ptr->numberOfSegments = 2;
                context_ptr->scoreTh[0] = (int8_t)((1 * 100) / context_ptr->numberOfSegments);
                context_ptr->intervalCost[0] = context_ptr->costDepthMode[LCU_LIGHT_BDP_DEPTH_MODE - 1];
                context_ptr->intervalCost[1] = context_ptr->costDepthMode[LCU_BDP_DEPTH_MODE - 1];
            }

        }
        else {
            if (context_ptr->budget > (uint32_t)(picture_control_set_ptr->parent_pcs_ptr->sb_total_count * BDP_COST)) {

                context_ptr->numberOfSegments = 2;

                context_ptr->scoreTh[0] = (int8_t)((1 * 100) / context_ptr->numberOfSegments);

                context_ptr->intervalCost[0] = context_ptr->costDepthMode[LCU_BDP_DEPTH_MODE - 1];
                context_ptr->intervalCost[1] = context_ptr->costDepthMode[LCU_FULL84_DEPTH_MODE - 1];
            }
            else if (context_ptr->budget > (uint32_t)(picture_control_set_ptr->parent_pcs_ptr->sb_total_count * U_120)) {
                context_ptr->numberOfSegments = 4;

                context_ptr->scoreTh[0] = (int8_t)((1 * 100) / context_ptr->numberOfSegments);
                context_ptr->scoreTh[1] = (int8_t)((2 * 100) / context_ptr->numberOfSegments);
                context_ptr->scoreTh[2] = (int8_t)((3 * 100) / context_ptr->numberOfSegments);

                context_ptr->intervalCost[0] = context_ptr->costDepthMode[LCU_PRED_OPEN_LOOP_DEPTH_MODE - 1];
                context_ptr->intervalCost[1] = context_ptr->costDepthMode[LCU_LIGHT_OPEN_LOOP_DEPTH_MODE - 1];
                context_ptr->intervalCost[2] = context_ptr->costDepthMode[LCU_LIGHT_BDP_DEPTH_MODE - 1];
                context_ptr->intervalCost[3] = context_ptr->costDepthMode[LCU_BDP_DEPTH_MODE - 1];
            }
            else {
                context_ptr->numberOfSegments = 5;

                context_ptr->scoreTh[0] = (int8_t)((1 * 100) / context_ptr->numberOfSegments);
                context_ptr->scoreTh[1] = (int8_t)((2 * 100) / context_ptr->numberOfSegments);
                context_ptr->scoreTh[2] = (int8_t)((3 * 100) / context_ptr->numberOfSegments);
                context_ptr->scoreTh[3] = (int8_t)((4 * 100) / context_ptr->numberOfSegments);

                context_ptr->intervalCost[0] = context_ptr->costDepthMode[LCU_PRED_OPEN_LOOP_1_NFL_DEPTH_MODE - 1];
                context_ptr->intervalCost[1] = context_ptr->costDepthMode[LCU_PRED_OPEN_LOOP_DEPTH_MODE - 1];
                context_ptr->intervalCost[2] = context_ptr->costDepthMode[LCU_LIGHT_OPEN_LOOP_DEPTH_MODE - 1];
                context_ptr->intervalCost[3] = context_ptr->costDepthMode[LCU_LIGHT_BDP_DEPTH_MODE - 1];
                context_ptr->intervalCost[4] = context_ptr->costDepthMode[LCU_BDP_DEPTH_MODE - 1];
            }

        }
    }
    else {

        if (context_ptr->budget > (uint32_t)(picture_control_set_ptr->parent_pcs_ptr->sb_total_count * U_120)) {

            context_ptr->numberOfSegments = 6;

            context_ptr->scoreTh[0] = (int8_t)((1 * 100) / context_ptr->numberOfSegments);
            context_ptr->scoreTh[1] = (int8_t)((2 * 100) / context_ptr->numberOfSegments);
            context_ptr->scoreTh[2] = (int8_t)((3 * 100) / context_ptr->numberOfSegments);
            context_ptr->scoreTh[3] = (int8_t)((4 * 100) / context_ptr->numberOfSegments);
            context_ptr->scoreTh[4] = (int8_t)((5 * 100) / context_ptr->numberOfSegments);

            context_ptr->intervalCost[0] = context_ptr->costDepthMode[LCU_PRED_OPEN_LOOP_DEPTH_MODE - 1];
            context_ptr->intervalCost[1] = context_ptr->costDepthMode[LCU_LIGHT_OPEN_LOOP_DEPTH_MODE - 1];
            context_ptr->intervalCost[2] = context_ptr->costDepthMode[LCU_OPEN_LOOP_DEPTH_MODE - 1];
            context_ptr->intervalCost[3] = context_ptr->costDepthMode[LCU_LIGHT_BDP_DEPTH_MODE - 1];
            context_ptr->intervalCost[4] = context_ptr->costDepthMode[LCU_BDP_DEPTH_MODE - 1];
            context_ptr->intervalCost[5] = context_ptr->costDepthMode[LCU_FULL85_DEPTH_MODE - 1];
        }
        else if (context_ptr->budget > (uint32_t)(picture_control_set_ptr->parent_pcs_ptr->sb_total_count * U_115)) {

            context_ptr->numberOfSegments = 5;

            context_ptr->scoreTh[0] = (int8_t)((1 * 100) / context_ptr->numberOfSegments);
            context_ptr->scoreTh[1] = (int8_t)((2 * 100) / context_ptr->numberOfSegments);
            context_ptr->scoreTh[2] = (int8_t)((3 * 100) / context_ptr->numberOfSegments);
            context_ptr->scoreTh[3] = (int8_t)((4 * 100) / context_ptr->numberOfSegments);

            context_ptr->intervalCost[0] = context_ptr->costDepthMode[LCU_PRED_OPEN_LOOP_DEPTH_MODE - 1];
            context_ptr->intervalCost[1] = context_ptr->costDepthMode[LCU_LIGHT_OPEN_LOOP_DEPTH_MODE - 1];
            context_ptr->intervalCost[2] = context_ptr->costDepthMode[LCU_OPEN_LOOP_DEPTH_MODE - 1];
            context_ptr->intervalCost[3] = context_ptr->costDepthMode[LCU_LIGHT_BDP_DEPTH_MODE - 1];
            context_ptr->intervalCost[4] = context_ptr->costDepthMode[LCU_BDP_DEPTH_MODE - 1];

        }
        else if (context_ptr->budget > (uint32_t)(picture_control_set_ptr->parent_pcs_ptr->sb_total_count * OPEN_LOOP_COST)) {

            context_ptr->numberOfSegments = 4;

            context_ptr->scoreTh[0] = (int8_t)((1 * 100) / context_ptr->numberOfSegments);
            context_ptr->scoreTh[1] = (int8_t)((2 * 100) / context_ptr->numberOfSegments);
            context_ptr->scoreTh[2] = (int8_t)((3 * 100) / context_ptr->numberOfSegments);

            context_ptr->intervalCost[0] = context_ptr->costDepthMode[LCU_PRED_OPEN_LOOP_DEPTH_MODE - 1];
            context_ptr->intervalCost[1] = context_ptr->costDepthMode[LCU_LIGHT_OPEN_LOOP_DEPTH_MODE - 1];
            context_ptr->intervalCost[2] = context_ptr->costDepthMode[LCU_OPEN_LOOP_DEPTH_MODE - 1];
            context_ptr->intervalCost[3] = context_ptr->costDepthMode[LCU_LIGHT_BDP_DEPTH_MODE - 1];

        }
        else {

            context_ptr->numberOfSegments = 4;

            context_ptr->scoreTh[0] = (int8_t)((1 * 100) / context_ptr->numberOfSegments);
            context_ptr->scoreTh[1] = (int8_t)((2 * 100) / context_ptr->numberOfSegments);
            context_ptr->scoreTh[2] = (int8_t)((3 * 100) / context_ptr->numberOfSegments);

            context_ptr->intervalCost[0] = context_ptr->costDepthMode[LCU_PRED_OPEN_LOOP_1_NFL_DEPTH_MODE - 1];
            context_ptr->intervalCost[1] = context_ptr->costDepthMode[LCU_PRED_OPEN_LOOP_DEPTH_MODE - 1];
            context_ptr->intervalCost[2] = context_ptr->costDepthMode[LCU_LIGHT_OPEN_LOOP_DEPTH_MODE - 1];
            context_ptr->intervalCost[3] = context_ptr->costDepthMode[LCU_OPEN_LOOP_DEPTH_MODE - 1];
        }
    }

    return return_error;
}





/******************************************************
* Set the target budget
    Input   : cost per depth
    Output  : budget per picture
******************************************************/


void SetTargetBudget(
    SequenceControlSet_t                *sequence_control_set_ptr,
    PictureControlSet_t                 *picture_control_set_ptr,
    ModeDecisionConfigurationContext_t  *context_ptr)
{

    uint32_t budget;
#if ENCODER_MODE_CLEANUP
    if (1) {
#else
    if (picture_control_set_ptr->enc_mode == ENC_M1) {
#endif
        if (sequence_control_set_ptr->input_resolution <= INPUT_SIZE_1080i_RANGE) {
            if (picture_control_set_ptr->temporal_layer_index == 0)
                budget = picture_control_set_ptr->sb_total_count * FULL_SEARCH_COST;
            else if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
                budget = picture_control_set_ptr->sb_total_count * U_150;
            else
                budget = picture_control_set_ptr->sb_total_count * U_145;
        }
        else

            if (sequence_control_set_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE) {
                if (picture_control_set_ptr->temporal_layer_index == 0)
                    budget = picture_control_set_ptr->sb_total_count * FULL_SEARCH_COST;
                else if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
                    budget = picture_control_set_ptr->sb_total_count * AVC_COST;
                else
                    budget = picture_control_set_ptr->sb_total_count * U_134;

            }
            else {

                if (picture_control_set_ptr->temporal_layer_index == 0)
                    budget = picture_control_set_ptr->sb_total_count * FULL_SEARCH_COST;
                else if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
                    budget = picture_control_set_ptr->sb_total_count * U_125;
                else
                    budget = picture_control_set_ptr->sb_total_count * U_121;

            }

    }
#if !ENCODER_MODE_CLEANUP
    else if (picture_control_set_ptr->enc_mode == ENC_M2) {
        if (sequence_control_set_ptr->input_resolution <= INPUT_SIZE_1080i_RANGE) {
            if (picture_control_set_ptr->temporal_layer_index == 0)
                budget = picture_control_set_ptr->sb_total_count * FULL_SEARCH_COST;
            else if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
                budget = picture_control_set_ptr->sb_total_count * U_133;
            else
                budget = picture_control_set_ptr->sb_total_count * U_120;
        }
        else

            if (sequence_control_set_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE) {
                if (picture_control_set_ptr->temporal_layer_index == 0)
                    budget = picture_control_set_ptr->parent_pcs_ptr->sb_total_count * FULL_SEARCH_COST;
                else if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
                    budget = picture_control_set_ptr->sb_total_count * U_130;
                else
                    budget = picture_control_set_ptr->sb_total_count * U_120;
            }
            else {

                if (picture_control_set_ptr->temporal_layer_index == 0)
                    budget = picture_control_set_ptr->parent_pcs_ptr->sb_total_count * FULL_SEARCH_COST;
                else if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)

                    budget = picture_control_set_ptr->sb_total_count * U_125;
                else
                    budget = picture_control_set_ptr->sb_total_count * U_115;

            }

    }
    else if (picture_control_set_ptr->enc_mode == ENC_M3) {
        if (sequence_control_set_ptr->input_resolution <= INPUT_SIZE_1080i_RANGE) {
            if (picture_control_set_ptr->temporal_layer_index == 0)
                budget = picture_control_set_ptr->sb_total_count * FULL_SEARCH_COST;
            else if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
                budget = picture_control_set_ptr->sb_total_count * U_115;
            else
                budget = picture_control_set_ptr->sb_total_count * U_114;
        }
        else

            if (sequence_control_set_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE) {

                if (picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index == 0)
                    budget = picture_control_set_ptr->parent_pcs_ptr->sb_total_count * FULL_SEARCH_COST;
                else if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
                    budget = picture_control_set_ptr->parent_pcs_ptr->sb_total_count * U_117;
                else
                    budget = picture_control_set_ptr->parent_pcs_ptr->sb_total_count * U_114;

            }
            else {
                if (picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index == 0)
                    budget = picture_control_set_ptr->parent_pcs_ptr->sb_total_count * AVC_COST;
                else if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
                    budget = picture_control_set_ptr->parent_pcs_ptr->sb_total_count * LIGHT_BDP_COST;
                else
                    budget = picture_control_set_ptr->parent_pcs_ptr->sb_total_count * U_111;
            }
    }
    else if (picture_control_set_ptr->enc_mode == ENC_M4) {
        if (sequence_control_set_ptr->input_resolution <= INPUT_SIZE_1080i_RANGE) {
            if (picture_control_set_ptr->temporal_layer_index == 0)
                budget = picture_control_set_ptr->sb_total_count * BDP_COST;
            else if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
                budget = picture_control_set_ptr->sb_total_count * U_113;
            else
                budget = picture_control_set_ptr->sb_total_count * U_112;
        }
        else

            if (sequence_control_set_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE) {


                if (picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index == 0)
                    budget = picture_control_set_ptr->parent_pcs_ptr->sb_total_count * BDP_COST;
                else if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
                    budget = picture_control_set_ptr->sb_total_count * U_116;
                else
                    budget = picture_control_set_ptr->sb_total_count * OPEN_LOOP_COST;

            }
            else {
                if (picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index == 0)
                    budget = picture_control_set_ptr->parent_pcs_ptr->sb_total_count * U_127;
                else if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
                    budget = picture_control_set_ptr->parent_pcs_ptr->sb_total_count * U_114;
                else
                    budget = picture_control_set_ptr->parent_pcs_ptr->sb_total_count * OPEN_LOOP_COST;

            }
    }
    else if (picture_control_set_ptr->enc_mode == ENC_M5) {
        if (sequence_control_set_ptr->input_resolution <= INPUT_SIZE_1080i_RANGE) {
            if (picture_control_set_ptr->temporal_layer_index == 0)
                budget = picture_control_set_ptr->sb_total_count * BDP_COST;
            else if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
                budget = picture_control_set_ptr->sb_total_count * U_112;
            else
                budget = picture_control_set_ptr->sb_total_count * U_111;
        }
        else

            if (sequence_control_set_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE) {

                if (picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index == 0)
                    budget = picture_control_set_ptr->parent_pcs_ptr->sb_total_count * BDP_COST;
                else if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
                    budget = picture_control_set_ptr->sb_total_count * U_113;
                else
                    budget = picture_control_set_ptr->sb_total_count * OPEN_LOOP_COST;

            }
            else {
                if (picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index == 0)
                    budget = (context_ptr->depthSensitivePictureFlag) ?
                    picture_control_set_ptr->parent_pcs_ptr->sb_total_count * U_127 :
                    picture_control_set_ptr->parent_pcs_ptr->sb_total_count * U_125;
                else if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag)
                    budget = picture_control_set_ptr->parent_pcs_ptr->sb_total_count * OPEN_LOOP_COST;
                else
                    budget = picture_control_set_ptr->parent_pcs_ptr->sb_total_count * LIGHT_OPEN_LOOP_COST;


            }
    }
    else {
        if (sequence_control_set_ptr->input_resolution <= INPUT_SIZE_1080i_RANGE) {
            if (picture_control_set_ptr->temporal_layer_index == 0)
                budget = picture_control_set_ptr->sb_total_count * BDP_COST;
            else if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
                budget = picture_control_set_ptr->sb_total_count * OPEN_LOOP_COST;
            else
                budget = picture_control_set_ptr->sb_total_count * U_109;
        }
        else

            if (sequence_control_set_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE) {

                if (picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index == 0)
                    budget = picture_control_set_ptr->parent_pcs_ptr->sb_total_count * BDP_COST;
                else if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
                    budget = picture_control_set_ptr->sb_total_count * U_112;
                else
                    budget = picture_control_set_ptr->sb_total_count * U_108;

            }
            else {
                if (picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index == 0)
                    budget = (context_ptr->depthSensitivePictureFlag) ?
                    picture_control_set_ptr->parent_pcs_ptr->sb_total_count * U_127 :
                    picture_control_set_ptr->parent_pcs_ptr->sb_total_count * U_125;
                else if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag)
                    budget = (context_ptr->depthSensitivePictureFlag) ?
                    picture_control_set_ptr->parent_pcs_ptr->sb_total_count * OPEN_LOOP_COST :
                    picture_control_set_ptr->parent_pcs_ptr->sb_total_count * U_108;
                else
                    budget = picture_control_set_ptr->parent_pcs_ptr->sb_total_count * U_104;


            }
    }
#endif
    context_ptr->budget = budget;
}





/******************************************************
 * IsAvcPartitioningMode()
 * Returns TRUE for LCUs where only Depth2 & Depth3
 * (AVC Partitioning) are goind to be tested by MD
 * The SB is marked if Sharpe Edge or Potential Aura/Grass
 * or B-Logo or S-Logo or Potential Blockiness Area
 * Input: Sharpe Edge, Potential Aura/Grass, B-Logo, S-Logo, Potential Blockiness Area signals
 * Output: TRUE if one of the above is TRUE
 ******************************************************/
EbBool IsAvcPartitioningMode(
    SequenceControlSet_t  *sequence_control_set_ptr,
    PictureControlSet_t   *picture_control_set_ptr,
    LargestCodingUnit_t   *sb_ptr)
{

    uint32_t          sb_index = sb_ptr->index;
    SbParams_t    *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];
    EB_SLICE        slice_type = picture_control_set_ptr->slice_type;
    uint8_t           edge_block_num = picture_control_set_ptr->parent_pcs_ptr->edge_results_ptr[sb_index].edge_block_num;
    SbStat_t      *sb_stat_ptr = &(picture_control_set_ptr->parent_pcs_ptr->sb_stat_array[sb_index]);
    uint8_t           stationary_edge_over_time_flag = sb_stat_ptr->stationary_edge_over_time_flag;
    uint8_t           aura_status_iii = sb_ptr->aura_status_iii;

    // No refinment for sub-1080p
    if (sequence_control_set_ptr->input_resolution <= INPUT_SIZE_1080i_RANGE)
        return EB_FALSE;

    // Sharpe Edge
    if (picture_control_set_ptr->parent_pcs_ptr->high_dark_low_light_area_density_flag && picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index > 0 && picture_control_set_ptr->parent_pcs_ptr->sharp_edge_sb_flag[sb_index] && !picture_control_set_ptr->parent_pcs_ptr->similar_colocated_sb_array_ii[sb_index]) {
        return EB_TRUE;
    }

    // Potential Aura/Grass
    if ((slice_type != I_SLICE && picture_control_set_ptr->high_intra_slection == 0) && (sb_params->is_complete_sb)) {
        if (picture_control_set_ptr->scene_caracteristic_id == EB_FRAME_CARAC_0) {
            if (picture_control_set_ptr->parent_pcs_ptr->grass_percentage_in_picture > 60 && aura_status_iii) {
                return EB_TRUE;
            }
        }
    }

    // B-Logo
    if (picture_control_set_ptr->parent_pcs_ptr->logo_pic_flag && edge_block_num)
        return EB_TRUE;

    // S-Logo
    if (stationary_edge_over_time_flag > 0)
        return EB_TRUE;

    // Potential Blockiness Area
    if (picture_control_set_ptr->parent_pcs_ptr->complex_sb_array[sb_ptr->index] == SB_COMPLEXITY_STATUS_2)

        return EB_TRUE;

    return EB_FALSE;
}






/******************************************************
* Load the cost of the different partitioning method into a local array and derive sensitive picture flag
    Input   : the offline derived cost per search method, detection signals
    Output  : valid costDepthMode and valid sensitivePicture
******************************************************/
void ConfigureAdp(
    PictureControlSet_t                 *picture_control_set_ptr,
    ModeDecisionConfigurationContext_t  *context_ptr)
{

    context_ptr->costDepthMode[LCU_FULL85_DEPTH_MODE - 1] = FULL_SEARCH_COST;
    context_ptr->costDepthMode[LCU_FULL84_DEPTH_MODE - 1] = FULL_SEARCH_COST;
    context_ptr->costDepthMode[LCU_BDP_DEPTH_MODE - 1] = BDP_COST;
    context_ptr->costDepthMode[LCU_LIGHT_BDP_DEPTH_MODE - 1] = LIGHT_BDP_COST;
    context_ptr->costDepthMode[LCU_OPEN_LOOP_DEPTH_MODE - 1] = OPEN_LOOP_COST;
    context_ptr->costDepthMode[LCU_LIGHT_OPEN_LOOP_DEPTH_MODE - 1] = LIGHT_OPEN_LOOP_COST;
    context_ptr->costDepthMode[LCU_AVC_DEPTH_MODE - 1] = AVC_COST;
    context_ptr->costDepthMode[LCU_PRED_OPEN_LOOP_DEPTH_MODE - 1] = PRED_OPEN_LOOP_COST;
    context_ptr->costDepthMode[LCU_PRED_OPEN_LOOP_1_NFL_DEPTH_MODE - 1] = PRED_OPEN_LOOP_1_NFL_COST;


    // Initialize the score based TH
    context_ptr->scoreTh[0] = ~0;
    context_ptr->scoreTh[1] = ~0;
    context_ptr->scoreTh[2] = ~0;
    context_ptr->scoreTh[3] = ~0;
    context_ptr->scoreTh[4] = ~0;
    context_ptr->scoreTh[5] = ~0;
    context_ptr->scoreTh[6] = ~0;

    // Initialize the predicted budget
    context_ptr->predictedCost = (uint32_t)~0;

    // Initialize the predicted budget
    context_ptr->predictedCost = (uint32_t)~0;

    // Derive the sensitive picture flag
    context_ptr->depthSensitivePictureFlag = EB_FALSE;
    EbBool luminosityChange = EB_FALSE;
    // Derived for REF P & B & kept false otherwise (for temporal distance equal to 1 luminosity changes are easier to handle)
    if (picture_control_set_ptr->slice_type != I_SLICE) {
        if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag) {
            EbReferenceObject_t  * refObjL0, *refObjL1;
            refObjL0 = (EbReferenceObject_t*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_0]->objectPtr;
            refObjL1 = (picture_control_set_ptr->parent_pcs_ptr->slice_type == B_SLICE) ? (EbReferenceObject_t*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_1]->objectPtr : (EbReferenceObject_t*)EB_NULL;
            luminosityChange = ((ABS(picture_control_set_ptr->parent_pcs_ptr->average_intensity[0] - refObjL0->average_intensity) >= LUMINOSITY_CHANGE_TH) || (refObjL1 != EB_NULL && ABS(picture_control_set_ptr->parent_pcs_ptr->average_intensity[0] - refObjL1->average_intensity) >= LUMINOSITY_CHANGE_TH));
        }
    }
    context_ptr->depthSensitivePictureFlag = (picture_control_set_ptr->parent_pcs_ptr->non_moving_index_average != INVALID_ZZ_COST &&
        picture_control_set_ptr->parent_pcs_ptr->non_moving_index_average < ADP_NON_MOVING_INDEX_TH_1 && // could not seen by the eye if very active
        ((picture_control_set_ptr->parent_pcs_ptr->non_moving_index_average >= ADP_NON_MOVING_INDEX_TH_0 && picture_control_set_ptr->parent_pcs_ptr->pic_noise_class >= PIC_NOISE_CLASS_3) || // potential complex picture: moderate activity and high variance (noise or a lot of edge)
            picture_control_set_ptr->parent_pcs_ptr->high_dark_low_light_area_density_flag || // potential complex picture: light foreground and dark background (e.g. flash, light..)
            luminosityChange)) ? // potential complex picture: luminosity Change (e.g. fade, light..)
        EB_TRUE :
        EB_FALSE;



}


/******************************************************
* Assign a search method based on the allocated cost
    Input   : allocated budget per LCU
    Output  : search method per LCU
******************************************************/
void DeriveSearchMethod(
    PictureControlSet_t                 *picture_control_set_ptr,
    ModeDecisionConfigurationContext_t  *context_ptr)
{

    uint32_t sb_index;

    for (sb_index = 0; sb_index < picture_control_set_ptr->parent_pcs_ptr->sb_total_count; sb_index++) {

        if (context_ptr->lcuCostArray[sb_index] == context_ptr->costDepthMode[LCU_PRED_OPEN_LOOP_1_NFL_DEPTH_MODE - 1]) {
            picture_control_set_ptr->parent_pcs_ptr->sb_md_mode_array[sb_index] = LCU_PRED_OPEN_LOOP_1_NFL_DEPTH_MODE;
        }
        else
            if (context_ptr->lcuCostArray[sb_index] == context_ptr->costDepthMode[LCU_PRED_OPEN_LOOP_DEPTH_MODE - 1]) {
                picture_control_set_ptr->parent_pcs_ptr->sb_md_mode_array[sb_index] = LCU_PRED_OPEN_LOOP_DEPTH_MODE;
            }
            else if (context_ptr->lcuCostArray[sb_index] == context_ptr->costDepthMode[LCU_LIGHT_OPEN_LOOP_DEPTH_MODE - 1]) {
                picture_control_set_ptr->parent_pcs_ptr->sb_md_mode_array[sb_index] = LCU_LIGHT_OPEN_LOOP_DEPTH_MODE;
            }
            else if (context_ptr->lcuCostArray[sb_index] == context_ptr->costDepthMode[LCU_OPEN_LOOP_DEPTH_MODE - 1]) {
                picture_control_set_ptr->parent_pcs_ptr->sb_md_mode_array[sb_index] = LCU_OPEN_LOOP_DEPTH_MODE;
            }
            else if (context_ptr->lcuCostArray[sb_index] == context_ptr->costDepthMode[LCU_LIGHT_BDP_DEPTH_MODE - 1]) {
                picture_control_set_ptr->parent_pcs_ptr->sb_md_mode_array[sb_index] = LCU_LIGHT_BDP_DEPTH_MODE;
            }
            else if (context_ptr->lcuCostArray[sb_index] == context_ptr->costDepthMode[LCU_BDP_DEPTH_MODE - 1]) {
                picture_control_set_ptr->parent_pcs_ptr->sb_md_mode_array[sb_index] = LCU_BDP_DEPTH_MODE;
            }
            else if (context_ptr->lcuCostArray[sb_index] == context_ptr->costDepthMode[LCU_AVC_DEPTH_MODE - 1]) {
                picture_control_set_ptr->parent_pcs_ptr->sb_md_mode_array[sb_index] = LCU_AVC_DEPTH_MODE;
            }
            else if (picture_control_set_ptr->temporal_layer_index == 0) {
                picture_control_set_ptr->parent_pcs_ptr->sb_md_mode_array[sb_index] = LCU_FULL84_DEPTH_MODE;
            }
            else {
                picture_control_set_ptr->parent_pcs_ptr->sb_md_mode_array[sb_index] = LCU_FULL85_DEPTH_MODE;
            }

    }

}
/******************************************************
* Set SB budget
    Input   : SB score, detection signals, iteration
    Output  : predicted budget for the LCU
******************************************************/
void SetLcuBudget(
    SequenceControlSet_t                *sequence_control_set_ptr,
    PictureControlSet_t                 *picture_control_set_ptr,
    LargestCodingUnit_t                 *sb_ptr,
    ModeDecisionConfigurationContext_t  *context_ptr) {

    uint32_t      sb_index = sb_ptr->index;
    uint32_t      maxToMinScore, scoreToMin;

    if (context_ptr->performRefinement && IsAvcPartitioningMode(sequence_control_set_ptr, picture_control_set_ptr, sb_ptr)) {

        context_ptr->lcuCostArray[sb_index] = context_ptr->costDepthMode[LCU_AVC_DEPTH_MODE - 1];
        context_ptr->predictedCost += context_ptr->costDepthMode[LCU_AVC_DEPTH_MODE - 1];

    }
    else {
        context_ptr->lcuScoreArray[sb_index] = CLIP3(context_ptr->lcuMinScore, context_ptr->lcuMaxScore, context_ptr->lcuScoreArray[sb_index]);
        scoreToMin = context_ptr->lcuScoreArray[sb_index] - context_ptr->lcuMinScore;
        maxToMinScore = context_ptr->lcuMaxScore - context_ptr->lcuMinScore;

        if ((scoreToMin <= (maxToMinScore * context_ptr->scoreTh[0]) / 100 && context_ptr->scoreTh[0] != 0) || context_ptr->numberOfSegments == 1 || context_ptr->scoreTh[1] == 100) {
            context_ptr->lcuCostArray[sb_index] = context_ptr->intervalCost[0];
            context_ptr->predictedCost += context_ptr->intervalCost[0];
        }
        else if ((scoreToMin <= (maxToMinScore * context_ptr->scoreTh[1]) / 100 && context_ptr->scoreTh[1] != 0) || context_ptr->numberOfSegments == 2 || context_ptr->scoreTh[2] == 100) {
            context_ptr->lcuCostArray[sb_index] = context_ptr->intervalCost[1];
            context_ptr->predictedCost += context_ptr->intervalCost[1];
        }
        else if ((scoreToMin <= (maxToMinScore * context_ptr->scoreTh[2]) / 100 && context_ptr->scoreTh[2] != 0) || context_ptr->numberOfSegments == 3 || context_ptr->scoreTh[3] == 100) {
            context_ptr->lcuCostArray[sb_index] = context_ptr->intervalCost[2];
            context_ptr->predictedCost += context_ptr->intervalCost[2];
        }
        else if ((scoreToMin <= (maxToMinScore * context_ptr->scoreTh[3]) / 100 && context_ptr->scoreTh[3] != 0) || context_ptr->numberOfSegments == 4 || context_ptr->scoreTh[4] == 100) {
            context_ptr->lcuCostArray[sb_index] = context_ptr->intervalCost[3];
            context_ptr->predictedCost += context_ptr->intervalCost[3];
        }
        else if ((scoreToMin <= (maxToMinScore * context_ptr->scoreTh[4]) / 100 && context_ptr->scoreTh[4] != 0) || context_ptr->numberOfSegments == 5 || context_ptr->scoreTh[5] == 100) {
            context_ptr->lcuCostArray[sb_index] = context_ptr->intervalCost[4];
            context_ptr->predictedCost += context_ptr->intervalCost[4];
        }
        else if ((scoreToMin <= (maxToMinScore * context_ptr->scoreTh[5]) / 100 && context_ptr->scoreTh[5] != 0) || context_ptr->numberOfSegments == 6 || context_ptr->scoreTh[6] == 100) {
            context_ptr->lcuCostArray[sb_index] = context_ptr->intervalCost[5];
            context_ptr->predictedCost += context_ptr->intervalCost[5];
        }
        else {
            context_ptr->lcuCostArray[sb_index] = context_ptr->intervalCost[6];
            context_ptr->predictedCost += context_ptr->intervalCost[6];
        }
        // Switch to AVC mode if the SB cost is higher than the AVC cost and the the SB is marked + adjust the current picture cost accordingly
        if (IsAvcPartitioningMode(sequence_control_set_ptr, picture_control_set_ptr, sb_ptr) && context_ptr->lcuCostArray[sb_index] > context_ptr->costDepthMode[LCU_AVC_DEPTH_MODE - 1]) {
            context_ptr->predictedCost -= (context_ptr->lcuCostArray[sb_index] - context_ptr->costDepthMode[LCU_AVC_DEPTH_MODE - 1]);
            context_ptr->lcuCostArray[sb_index] = context_ptr->costDepthMode[LCU_AVC_DEPTH_MODE - 1];
        }


    }
}

/******************************************************
* Loop multiple times over the LCUs in order to derive the optimal budget per LCU
    Input   : budget per picture, ditortion, detection signals, iteration
    Output  : optimal budget for each LCU
******************************************************/
void  DeriveOptimalBudgetPerLcu(
    SequenceControlSet_t                *sequence_control_set_ptr,
    PictureControlSet_t                 *picture_control_set_ptr,
    ModeDecisionConfigurationContext_t  *context_ptr) {


    uint32_t sb_index;
    // Initialize the deviation between the picture predicted cost & the target budget to 100,
    uint32_t deviationToTarget = 1000;

    // Set the adjustment step to 1 (could be increased for faster convergence),
    int8_t  adjustementStep = 1;

    // Set the initial shooting state & the final shooting state to TBD
    uint32_t initialShooting = TBD_SHOOTING;
    uint32_t finalShooting = TBD_SHOOTING;

    uint8_t maxAdjustementIteration = 100;
    uint8_t adjustementIteration = 0;

    while (deviationToTarget != 0 && (initialShooting == finalShooting) && adjustementIteration <= maxAdjustementIteration) {

        if (context_ptr->predictedCost < context_ptr->budget) {
            initialShooting = UNDER_SHOOTING;
        }
        else {
            initialShooting = OVER_SHOOTING;
        }

        // reset running cost
        context_ptr->predictedCost = 0;

        for (sb_index = 0; sb_index < picture_control_set_ptr->parent_pcs_ptr->sb_total_count; sb_index++) {

            LargestCodingUnit_t* sb_ptr = picture_control_set_ptr->sb_ptr_array[sb_index];

            SetLcuBudget(
                sequence_control_set_ptr,
                picture_control_set_ptr,
                sb_ptr,
                context_ptr);
        }

        // Compute the deviation between the predicted budget & the target budget
        deviationToTarget = (ABS((int32_t)(context_ptr->predictedCost - context_ptr->budget)) * 1000) / context_ptr->budget;

        // Derive shooting status
        if (context_ptr->predictedCost < context_ptr->budget) {
            context_ptr->scoreTh[0] = MAX((context_ptr->scoreTh[0] - adjustementStep), 0);
            context_ptr->scoreTh[1] = MAX((context_ptr->scoreTh[1] - adjustementStep), 0);
            context_ptr->scoreTh[2] = MAX((context_ptr->scoreTh[2] - adjustementStep), 0);
            context_ptr->scoreTh[3] = MAX((context_ptr->scoreTh[3] - adjustementStep), 0);
            context_ptr->scoreTh[4] = MAX((context_ptr->scoreTh[4] - adjustementStep), 0);
            finalShooting = UNDER_SHOOTING;
        }
        else {
            context_ptr->scoreTh[0] = (context_ptr->scoreTh[0] == 0) ? 0 : MIN(context_ptr->scoreTh[0] + adjustementStep, 100);
            context_ptr->scoreTh[1] = (context_ptr->scoreTh[1] == 0) ? 0 : MIN(context_ptr->scoreTh[1] + adjustementStep, 100);
            context_ptr->scoreTh[2] = (context_ptr->scoreTh[2] == 0) ? 0 : MIN(context_ptr->scoreTh[2] + adjustementStep, 100);
            context_ptr->scoreTh[3] = (context_ptr->scoreTh[3] == 0) ? 0 : MIN(context_ptr->scoreTh[3] + adjustementStep, 100);
            context_ptr->scoreTh[4] = (context_ptr->scoreTh[4] == 0) ? 0 : MIN(context_ptr->scoreTh[4] + adjustementStep, 100);
            finalShooting = OVER_SHOOTING;
        }

        if (adjustementIteration == 0)
            initialShooting = finalShooting;

        adjustementIteration++;
    }
}


/******************************************************
* Compute the refinment cost
    Input   : budget per picture, and the cost of the refinment
    Output  : the refinment flag
******************************************************/
void ComputeRefinementCost(
    SequenceControlSet_t                *sequence_control_set_ptr,
    PictureControlSet_t                 *picture_control_set_ptr,
    ModeDecisionConfigurationContext_t  *context_ptr)
{

    uint32_t  sb_index;
    uint32_t  refinementCost = 0;

    for (sb_index = 0; sb_index < picture_control_set_ptr->parent_pcs_ptr->sb_total_count; sb_index++) {
        LargestCodingUnit_t* sb_ptr = picture_control_set_ptr->sb_ptr_array[sb_index];
        if (IsAvcPartitioningMode(sequence_control_set_ptr, picture_control_set_ptr, sb_ptr)) {
            refinementCost += context_ptr->costDepthMode[LCU_AVC_DEPTH_MODE - 1];
        }
        // assumes the fastest mode will be used otherwise
        else {
            refinementCost += context_ptr->intervalCost[0];
        }

    }

    // Shut the refinement if the refinement cost is higher than the allocated budget
    if (refinementCost > context_ptr->budget) {
        context_ptr->performRefinement = EB_FALSE;
    }
    else {
        context_ptr->performRefinement = EB_TRUE;
    }
}
/******************************************************
* Compute the score of each LCU
    Input   : distortion, detection signals
    Output  : SB score
******************************************************/
void DeriveLcuScore(
    SequenceControlSet_t               *sequence_control_set_ptr,
    PictureControlSet_t                *picture_control_set_ptr,
    ModeDecisionConfigurationContext_t *context_ptr)
{
    uint32_t  sb_index;
    uint32_t  lcuScore;
    uint32_t  distortion;

    context_ptr->lcuMinScore = ~0u;
    context_ptr->lcuMaxScore = 0u;

    for (sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; sb_index++) {

        SbParams_t sb_params = sequence_control_set_ptr->sb_params_array[sb_index];

        if (picture_control_set_ptr->slice_type == I_SLICE) {

            if (sb_params.raster_scan_cu_validity[RASTER_SCAN_CU_INDEX_64x64] == EB_FALSE) {

                uint8_t cu8x8Index;
                uint8_t validCu8x8Count = 0;
                distortion = 0;
                for (cu8x8Index = 0; cu8x8Index < 64; cu8x8Index++) {
                    if (sb_params.raster_scan_cu_validity[cu8x8Index]) {
                        distortion += picture_control_set_ptr->parent_pcs_ptr->ois_cu8_results[sb_index]->sorted_ois_candidate[cu8x8Index][0].distortion;
                        validCu8x8Count++;
                    }
                }
                distortion = CLIP3(picture_control_set_ptr->parent_pcs_ptr->intra_complexity_min_pre, picture_control_set_ptr->parent_pcs_ptr->intra_complexity_max_pre, (distortion / validCu8x8Count) * 64);
            }
            else {
                distortion = picture_control_set_ptr->parent_pcs_ptr->ois_cu32_cu16_results[sb_index]->sorted_ois_candidate[RASTER_SCAN_CU_INDEX_32x32_0][0].distortion +
                    picture_control_set_ptr->parent_pcs_ptr->ois_cu32_cu16_results[sb_index]->sorted_ois_candidate[RASTER_SCAN_CU_INDEX_32x32_1][0].distortion +
                    picture_control_set_ptr->parent_pcs_ptr->ois_cu32_cu16_results[sb_index]->sorted_ois_candidate[RASTER_SCAN_CU_INDEX_32x32_2][0].distortion +
                    picture_control_set_ptr->parent_pcs_ptr->ois_cu32_cu16_results[sb_index]->sorted_ois_candidate[RASTER_SCAN_CU_INDEX_32x32_3][0].distortion;
            }
            lcuScore = distortion;

        }
        else {
            if (sb_params.raster_scan_cu_validity[RASTER_SCAN_CU_INDEX_64x64] == EB_FALSE) {

                uint8_t cu8x8Index;
                uint8_t validCu8x8Count = 0;
                distortion = 0;
                for (cu8x8Index = RASTER_SCAN_CU_INDEX_8x8_0; cu8x8Index <= RASTER_SCAN_CU_INDEX_8x8_63; cu8x8Index++) {
                    if (sb_params.raster_scan_cu_validity[cu8x8Index]) {
                        distortion += picture_control_set_ptr->parent_pcs_ptr->me_results[sb_index][cu8x8Index].distortionDirection[0].distortion;
                        validCu8x8Count++;
                    }
                }
                distortion = CLIP3(picture_control_set_ptr->parent_pcs_ptr->inter_complexity_min_pre, picture_control_set_ptr->parent_pcs_ptr->inter_complexity_max_pre, (distortion / validCu8x8Count) * 64);
            }
            else {
                distortion = picture_control_set_ptr->parent_pcs_ptr->me_results[sb_index][RASTER_SCAN_CU_INDEX_64x64].distortionDirection[0].distortion;
            }

            lcuScore = distortion;

            // Use uncovered area
            if (picture_control_set_ptr->parent_pcs_ptr->failing_motion_sb_flag[sb_index]) {
                lcuScore = picture_control_set_ptr->parent_pcs_ptr->inter_complexity_max_pre;
            }
            else {

                // Use SB variance & activity
                if (picture_control_set_ptr->parent_pcs_ptr->non_moving_index_array[sb_index] == ADP_NON_MOVING_INDEX_TH_1 && picture_control_set_ptr->parent_pcs_ptr->variance[sb_index][RASTER_SCAN_CU_INDEX_64x64] > IS_COMPLEX_LCU_VARIANCE_TH)
                    lcuScore -= (((lcuScore - picture_control_set_ptr->parent_pcs_ptr->inter_complexity_min_pre) * ADP_CLASS_SHIFT_DIST_0) / 100);
                // Use SB luminosity
                if (picture_control_set_ptr->parent_pcs_ptr->yMean[sb_index][RASTER_SCAN_CU_INDEX_64x64] < ADP_DARK_LCU_TH || picture_control_set_ptr->parent_pcs_ptr->yMean[sb_index][RASTER_SCAN_CU_INDEX_64x64] > ADP_LIGHT_LCU_TH)
                    lcuScore -= (((lcuScore - picture_control_set_ptr->parent_pcs_ptr->inter_complexity_min_pre) * ADP_CLASS_SHIFT_DIST_0) / 100);
                else
                    lcuScore += (((picture_control_set_ptr->parent_pcs_ptr->inter_complexity_max_pre - lcuScore) * ADP_CLASS_SHIFT_DIST_0) / 100);
            }
        }

        context_ptr->lcuScoreArray[sb_index] = lcuScore;

        // Track MIN & MAX SB scores
        context_ptr->lcuMinScore = MIN(context_ptr->lcuScoreArray[sb_index], context_ptr->lcuMinScore);
        context_ptr->lcuMaxScore = MAX(context_ptr->lcuScoreArray[sb_index], context_ptr->lcuMaxScore);
    }



}

/******************************************************
* BudgetingOutlierRemovalLcu
    Input   : SB score histogram
    Output  : Adjusted min & max SB score (to be used to clip the SB score @ SetLcuBudget)
 Performs outlier removal by:
 1. dividing the total distance between the maximum lcuScore and the minimum lcuScore into NI intervals(NI = 10).
 For each lcuScore interval,
 2. computing the number of LCUs NV with lcuScore belonging to the subject lcuScore interval.
 3. marking the subject lcuScore interval as not valid if its NV represent less than validity threshold V_TH per - cent of the total number of the processed LCUs in the picture. (V_TH = 2)
 4. computing the distance MIN_D from 0 and the index of the first, in incremental way, lcuScore interval marked as valid in the prior step.
 5. computing the distance MAX_D from NI and the index of the first, in decreasing way, lcuScore interval marked as valid in the prior step.
 6. adjusting the minimum and maximum lcuScore as follows :
    MIN_SCORE = MIN_SCORE + MIN_D * I_Value.
    MAX_SCORE = MAX_SCORE - MAX_D * I_Value.
******************************************************/

void PerformOutlierRemoval(
    SequenceControlSet_t                *sequence_control_set_ptr,
    PictureParentControlSet_t           *picture_control_set_ptr,
    ModeDecisionConfigurationContext_t  *context_ptr)
{

    uint32_t maxInterval = 0;
    uint32_t subInterval = 0;
    uint32_t lcuScoreHistogram[10] = { 0 };
    uint32_t sb_index;
    uint32_t lcuScore;
    uint32_t processedlcuCount = 0;
    int32_t slot = 0;

    maxInterval = context_ptr->lcuMaxScore - context_ptr->lcuMinScore;
    // Consider 10 bins
    subInterval = maxInterval / 10;

    // Count # of LCUs at each bin
    for (sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; sb_index++) {

        SbParams_t sb_params = sequence_control_set_ptr->sb_params_array[sb_index];

        if (sb_params.raster_scan_cu_validity[RASTER_SCAN_CU_INDEX_64x64]) {

            processedlcuCount++;

            lcuScore = context_ptr->lcuScoreArray[sb_index] + context_ptr->lcuMinScore;
            if (lcuScore < (subInterval + context_ptr->lcuMinScore)) {
                lcuScoreHistogram[0]++;
            }
            else if (lcuScore < ((2 * subInterval) + context_ptr->lcuMinScore)) {
                lcuScoreHistogram[1]++;
            }
            else if (lcuScore < ((3 * subInterval) + context_ptr->lcuMinScore)) {
                lcuScoreHistogram[2]++;
            }
            else if (lcuScore < ((4 * subInterval) + context_ptr->lcuMinScore)) {
                lcuScoreHistogram[3]++;
            }
            else if (lcuScore < ((5 * subInterval) + context_ptr->lcuMinScore)) {
                lcuScoreHistogram[4]++;
            }
            else if (lcuScore < ((6 * subInterval) + context_ptr->lcuMinScore)) {
                lcuScoreHistogram[5]++;
            }
            else if (lcuScore < ((7 * subInterval) + context_ptr->lcuMinScore)) {
                lcuScoreHistogram[6]++;
            }
            else if (lcuScore < ((8 * subInterval) + context_ptr->lcuMinScore)) {
                lcuScoreHistogram[7]++;
            }
            else if (lcuScore < ((9 * subInterval) + context_ptr->lcuMinScore)) {
                lcuScoreHistogram[8]++;
            }
            else if (lcuScore < ((10 * subInterval) + context_ptr->lcuMinScore)) {
                lcuScoreHistogram[9]++;
            }
        }
    }

    // Zero-out the bin if percentage lower than VALID_SLOT_TH
    for (slot = 0; slot < 10; slot++) {
        if ((lcuScoreHistogram[slot] * 100 / processedlcuCount) < VALID_SLOT_TH) {
            lcuScoreHistogram[slot] = 0;
        }
    }

    // Ignore null bins
    for (slot = 0; slot < 10; slot++) {
        if (lcuScoreHistogram[slot]) {
            context_ptr->lcuMinScore = context_ptr->lcuMinScore + (slot * subInterval);
            break;
        }
    }

    for (slot = 9; slot >= 0; slot--) {
        if (lcuScoreHistogram[slot]) {
            context_ptr->lcuMaxScore = context_ptr->lcuMaxScore - ((9 - slot) * subInterval);
            break;
        }
    }
}
/******************************************************
* Assign a search method for each LCU
    Input   : SB score, detection signals
    Output  : search method for each LCU
******************************************************/
void DeriveLcuMdMode(
    SequenceControlSet_t                *sequence_control_set_ptr,
    PictureControlSet_t                 *picture_control_set_ptr,
    ModeDecisionConfigurationContext_t  *context_ptr) {

    // Configure ADP
    ConfigureAdp(
        picture_control_set_ptr,
        context_ptr);

    // Set the target budget
    SetTargetBudget(
        sequence_control_set_ptr,
        picture_control_set_ptr,
        context_ptr);

    // Set the percentage based thresholds
    DeriveDefaultSegments(
        picture_control_set_ptr,
        context_ptr);

    // Compute the cost of the refinements
    ComputeRefinementCost(
        sequence_control_set_ptr,
        picture_control_set_ptr,
        context_ptr);

    // Derive SB score
    DeriveLcuScore(
        sequence_control_set_ptr,
        picture_control_set_ptr,
        context_ptr);

    // Remove the outliers
    PerformOutlierRemoval(
        sequence_control_set_ptr,
        picture_control_set_ptr->parent_pcs_ptr,
        context_ptr);

    // Perform Budgetting
    DeriveOptimalBudgetPerLcu(
        sequence_control_set_ptr,
        picture_control_set_ptr,
        context_ptr);

    // Set the search method using the SB cost (mapping)
    DeriveSearchMethod(
        picture_control_set_ptr,
        context_ptr);

}
void forward_sq_non4_blocks_to_md(
    SequenceControlSet_t                   *sequence_control_set_ptr,
    PictureControlSet_t                    *picture_control_set_ptr)
{

    uint32_t                   sb_index;
    EbBool                  split_flag;


    for (sb_index = 0; sb_index < sequence_control_set_ptr->sb_tot_cnt; ++sb_index)
    {

        MdcLcuData_t *resultsPtr = &picture_control_set_ptr->mdc_sb_array[sb_index];

        resultsPtr->leaf_count = 0;

        uint32_t  blk_index = picture_control_set_ptr->slice_type == I_SLICE && sequence_control_set_ptr->sb_size == BLOCK_128X128 ? 17 : 0;



        while (blk_index < sequence_control_set_ptr->max_block_cnt)
        {
            split_flag = EB_TRUE;

            const BlockGeom * blk_geom = Get_blk_geom_mds(blk_index);


            //if the parentSq is inside inject this block
            if (sequence_control_set_ptr->sb_geom[sb_index].block_is_inside_md_scan[blk_index])

            {

                //int32_t offset_d1 = ns_blk_offset[(int32_t)from_shape_to_part[blk_geom->shape]]; //cu_ptr->best_d1_blk; // TOCKECK
                //int32_t num_d1_block = ns_blk_num[(int32_t)from_shape_to_part[blk_geom->shape]]; // context_ptr->blk_geom->totns; // TOCKECK
                //
                //                                                  // for (int32_t d1_itr = blk_it; d1_itr < blk_it + num_d1_block; d1_itr++) {
                // for (int32_t d1_itr = (int32_t)blk_index ; d1_itr < (int32_t)blk_index +  num_d1_block ; d1_itr++) {

                resultsPtr->leaf_data_array[resultsPtr->leaf_count].tot_d1_blocks = 1;


                resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = 0;//valid only for square 85 world. will be removed.
                resultsPtr->leaf_data_array[resultsPtr->leaf_count].mds_idx = blk_index;

                if (blk_geom->sq_size > 8)
                {
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = EB_TRUE;
                    split_flag = EB_TRUE;
                }
                else {
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = EB_FALSE;
                    split_flag = EB_FALSE;
                }


            }

            blk_index += split_flag ? d1_depth_offset[sequence_control_set_ptr->sb_size == BLOCK_128X128][blk_geom->depth] : ns_depth_offset[sequence_control_set_ptr->sb_size == BLOCK_128X128][blk_geom->depth];

        }

    }

    picture_control_set_ptr->parent_pcs_ptr->average_qp = (uint8_t)picture_control_set_ptr->parent_pcs_ptr->picture_qp;

}

void forward_all_c_blocks_to_md(
    SequenceControlSet_t   *sequence_control_set_ptr,
    PictureControlSet_t    *picture_control_set_ptr){

    uint32_t                sb_index;
    for (sb_index = 0; sb_index < sequence_control_set_ptr->sb_tot_cnt; ++sb_index){

        MdcLcuData_t *resultsPtr = &picture_control_set_ptr->mdc_sb_array[sb_index];
        resultsPtr->leaf_count = 0;
        uint32_t blk_index = 0;
        uint32_t tot_d1_blocks;

        while (blk_index < sequence_control_set_ptr->max_block_cnt)
        {
            tot_d1_blocks = 0;
            const BlockGeom * blk_geom = Get_blk_geom_mds(blk_index);

            //if the parentSq is inside inject this block
            uint8_t is_blk_allowed = picture_control_set_ptr->slice_type != I_SLICE ? 1 : (blk_geom->sq_size < 128) ? 1 : 0;

            if (sequence_control_set_ptr->sb_geom[sb_index].block_is_inside_md_scan[blk_index] && is_blk_allowed){

                tot_d1_blocks = resultsPtr->leaf_data_array[resultsPtr->leaf_count].tot_d1_blocks =

                    blk_geom->sq_size == 128 ? 17 :
                    blk_geom->sq_size > 16 ? 25 :
                    blk_geom->sq_size == 16 ? 17 :
                    blk_geom->sq_size == 8 ? 1 : 1;

                for (uint32_t idx = 0; idx < tot_d1_blocks; ++idx) {
                    blk_geom = Get_blk_geom_mds(blk_index);

                    //if the parentSq is inside inject this block
                    if (sequence_control_set_ptr->sb_geom[sb_index].block_is_inside_md_scan[blk_index]){

                        resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = 0;//valid only for square 85 world. will be removed.
                        resultsPtr->leaf_data_array[resultsPtr->leaf_count].mds_idx = blk_index;
                        if (blk_geom->sq_size > 4)
                        {
                            resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = EB_TRUE;
                        }
                        else {
                            resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = EB_FALSE;
                        }
                    }
                    blk_index++;
                }
            }
            blk_index += (d1_depth_offset[sequence_control_set_ptr->sb_size == BLOCK_128X128][blk_geom->depth] - tot_d1_blocks);

        }

    }

    picture_control_set_ptr->parent_pcs_ptr->average_qp = (uint8_t)picture_control_set_ptr->parent_pcs_ptr->picture_qp;
}
/******************************************************
 * Mode Decision Configuration Kernel
 ******************************************************/
void* ModeDecisionConfigurationKernel(void *input_ptr)
{
    // Context & SCS & PCS
    ModeDecisionConfigurationContext_t         *context_ptr = (ModeDecisionConfigurationContext_t*)input_ptr;
    PictureControlSet_t                        *picture_control_set_ptr;
    SequenceControlSet_t                       *sequence_control_set_ptr;

    // Input
    EbObjectWrapper_t                          *rateControlResultsWrapperPtr;
    RateControlResults_t                       *rateControlResultsPtr;

    // Output
    EbObjectWrapper_t                          *encDecTasksWrapperPtr;
    EncDecTasks_t                              *encDecTasksPtr;

    for (;;) {

        // Get RateControl Results
        EbGetFullObject(
            context_ptr->rateControlInputFifoPtr,
            &rateControlResultsWrapperPtr);

        rateControlResultsPtr = (RateControlResults_t*)rateControlResultsWrapperPtr->objectPtr;
        picture_control_set_ptr = (PictureControlSet_t*)rateControlResultsPtr->pictureControlSetWrapperPtr->objectPtr;
        sequence_control_set_ptr = (SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->objectPtr;
      
        context_ptr->qp = picture_control_set_ptr->picture_qp;

        picture_control_set_ptr->parent_pcs_ptr->average_qp = 0;

        picture_control_set_ptr->intra_coded_area = 0;

        picture_control_set_ptr->scene_caracteristic_id = EB_FRAME_CARAC_0;
#if ENCODER_MODE_CLEANUP
        EbPicnoiseClass picNoiseClassTH = PIC_NOISE_CLASS_1;
#else
        EbPicnoiseClass picNoiseClassTH = (picture_control_set_ptr->enc_mode <= ENC_M3) ? PIC_NOISE_CLASS_1 : PIC_NOISE_CLASS_3;
#endif

        picture_control_set_ptr->scene_caracteristic_id = (
            (!picture_control_set_ptr->parent_pcs_ptr->is_pan) &&
            (!picture_control_set_ptr->parent_pcs_ptr->is_tilt) &&
            (picture_control_set_ptr->parent_pcs_ptr->grass_percentage_in_picture > 0) &&
            (picture_control_set_ptr->parent_pcs_ptr->grass_percentage_in_picture <= 35) &&
            (picture_control_set_ptr->parent_pcs_ptr->pic_noise_class >= picNoiseClassTH) &&
            (picture_control_set_ptr->parent_pcs_ptr->pic_homogenous_over_time_sb_percentage < 50)) ? EB_FRAME_CARAC_1 : picture_control_set_ptr->scene_caracteristic_id;

        picture_control_set_ptr->scene_caracteristic_id = (
            (picture_control_set_ptr->parent_pcs_ptr->is_pan) &&
            (!picture_control_set_ptr->parent_pcs_ptr->is_tilt) &&
            (picture_control_set_ptr->parent_pcs_ptr->grass_percentage_in_picture > 35) &&
            (picture_control_set_ptr->parent_pcs_ptr->grass_percentage_in_picture <= 70) &&
            (picture_control_set_ptr->parent_pcs_ptr->pic_noise_class >= picNoiseClassTH) &&
            (picture_control_set_ptr->parent_pcs_ptr->pic_homogenous_over_time_sb_percentage < 50)) ? EB_FRAME_CARAC_2 : picture_control_set_ptr->scene_caracteristic_id;

        picture_control_set_ptr->adjust_min_qp_flag = (EbBool)((!picture_control_set_ptr->parent_pcs_ptr->is_pan) &&
            (!picture_control_set_ptr->parent_pcs_ptr->is_tilt) &&
            (picture_control_set_ptr->parent_pcs_ptr->grass_percentage_in_picture > 2) &&
            (picture_control_set_ptr->parent_pcs_ptr->grass_percentage_in_picture <= 35) &&
            (picture_control_set_ptr->parent_pcs_ptr->pic_homogenous_over_time_sb_percentage < 70) &&
            (picture_control_set_ptr->parent_pcs_ptr->zz_cost_average > 15) &&
            (picture_control_set_ptr->parent_pcs_ptr->pic_noise_class >= picNoiseClassTH));

        // Compute picture and slice level chroma QP offsets
        SetSliceAndPictureChromaQpOffsets( // HT done
            picture_control_set_ptr);

        // Compute Tc, and Beta offsets for a given picture
        AdaptiveDlfParameterComputation( // HT done
            context_ptr,
            sequence_control_set_ptr,
            picture_control_set_ptr,
            picture_control_set_ptr->slice_type == I_SLICE ? EB_FALSE : picture_control_set_ptr->parent_pcs_ptr->scene_transition_flag[REF_LIST_0]);

#if FAST_CDEF
        // Set reference cdef strength 
        set_reference_cdef_strength(
            picture_control_set_ptr);
#endif

#if FAST_SG
        // Set reference sg ep 
        set_reference_sg_ep(
            picture_control_set_ptr);
#endif
        SetGlobalMotionField(
            picture_control_set_ptr);

        av1_qm_init(
            picture_control_set_ptr->parent_pcs_ptr);

        Quants *const quants = &picture_control_set_ptr->parent_pcs_ptr->quants;
        Dequants *const dequants = &picture_control_set_ptr->parent_pcs_ptr->deq;

        av1_set_quantizer(
            picture_control_set_ptr->parent_pcs_ptr,
#if NEW_QPS
            picture_control_set_ptr->parent_pcs_ptr->base_qindex);
#else
            quantizer_to_qindex[picture_control_set_ptr->picture_qp]);
#endif
        av1_build_quantizer(
            (aom_bit_depth_t)sequence_control_set_ptr->static_config.encoder_bit_depth,
            picture_control_set_ptr->parent_pcs_ptr->y_dc_delta_q,
            picture_control_set_ptr->parent_pcs_ptr->u_dc_delta_q,
            picture_control_set_ptr->parent_pcs_ptr->u_ac_delta_q,
            picture_control_set_ptr->parent_pcs_ptr->v_dc_delta_q,
            picture_control_set_ptr->parent_pcs_ptr->v_ac_delta_q,
            quants,
            dequants);

#if MD_10BIT_FIX
        Quants *const quantsMd = &picture_control_set_ptr->parent_pcs_ptr->quantsMd;
        Dequants *const dequantsMd = &picture_control_set_ptr->parent_pcs_ptr->deqMd;
        av1_build_quantizer(
            (aom_bit_depth_t)8,
            picture_control_set_ptr->parent_pcs_ptr->y_dc_delta_q,
            picture_control_set_ptr->parent_pcs_ptr->u_dc_delta_q,
            picture_control_set_ptr->parent_pcs_ptr->u_ac_delta_q,
            picture_control_set_ptr->parent_pcs_ptr->v_dc_delta_q,
            picture_control_set_ptr->parent_pcs_ptr->v_ac_delta_q,
            quantsMd,
            dequantsMd);
#endif

        if (picture_control_set_ptr->parent_pcs_ptr->pic_depth_mode == PIC_SB_SWITCH_DEPTH_MODE) {

            DeriveLcuMdMode(
                sequence_control_set_ptr,
                picture_control_set_ptr,
                context_ptr);


            uint32_t sb_index;

            // Rate estimation/QP
            PartitioningInitialization(
                sequence_control_set_ptr,
                picture_control_set_ptr,
                context_ptr);

            // SB Loop : Partitionnig Decision
            for (sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; ++sb_index) {


                if (picture_control_set_ptr->parent_pcs_ptr->sb_md_mode_array[sb_index] == LCU_FULL85_DEPTH_MODE) {

                    Forward85CuToModeDecisionLCU(
                        sequence_control_set_ptr,
                        picture_control_set_ptr,
                        sb_index);
                }
                else if (picture_control_set_ptr->parent_pcs_ptr->sb_md_mode_array[sb_index] == LCU_FULL84_DEPTH_MODE) {
                    Forward84CuToModeDecisionLCU(
                        sequence_control_set_ptr,
                        picture_control_set_ptr,
                        sb_index);
                }
                else if (picture_control_set_ptr->parent_pcs_ptr->sb_md_mode_array[sb_index] == LCU_OPEN_LOOP_DEPTH_MODE || picture_control_set_ptr->parent_pcs_ptr->sb_md_mode_array[sb_index] == LCU_LIGHT_OPEN_LOOP_DEPTH_MODE || picture_control_set_ptr->parent_pcs_ptr->sb_md_mode_array[sb_index] == LCU_AVC_DEPTH_MODE || picture_control_set_ptr->parent_pcs_ptr->sb_md_mode_array[sb_index] == LCU_PRED_OPEN_LOOP_DEPTH_MODE || picture_control_set_ptr->parent_pcs_ptr->sb_md_mode_array[sb_index] == LCU_PRED_OPEN_LOOP_1_NFL_DEPTH_MODE) {

                    // Predict the SB partitionning
                    PerformEarlyLcuPartitionningLcu( // HT done
                        context_ptr,
                        sequence_control_set_ptr,
                        picture_control_set_ptr,
                        sb_index);
                }
            }
        }
        else  if (picture_control_set_ptr->parent_pcs_ptr->pic_depth_mode == PIC_ALL_DEPTH_MODE) {

            forward_all_blocks_to_md(
                sequence_control_set_ptr,
                picture_control_set_ptr);
        }
        else  if (picture_control_set_ptr->parent_pcs_ptr->pic_depth_mode == PIC_ALL_C_DEPTH_MODE) {

            forward_all_c_blocks_to_md(
                sequence_control_set_ptr,
                picture_control_set_ptr);
        }
        else  if (picture_control_set_ptr->parent_pcs_ptr->pic_depth_mode == PIC_SQ_DEPTH_MODE) {

            forward_sq_blocks_to_md(
                sequence_control_set_ptr,
                picture_control_set_ptr);
        }
        else  if (picture_control_set_ptr->parent_pcs_ptr->pic_depth_mode == PIC_SQ_NON4_DEPTH_MODE) {

            forward_sq_non4_blocks_to_md(
                sequence_control_set_ptr,
                picture_control_set_ptr);
        }
        else if (picture_control_set_ptr->parent_pcs_ptr->pic_depth_mode >= PIC_OPEN_LOOP_DEPTH_MODE) {
            // Predict the SB partitionning
            PerformEarlyLcuPartitionning( // HT done
                context_ptr,
                sequence_control_set_ptr,
                picture_control_set_ptr);
        }
        else {   // (picture_control_set_ptr->parent_pcs_ptr->mdMode == PICT_BDP_DEPTH_MODE || picture_control_set_ptr->parent_pcs_ptr->mdMode == PICT_LIGHT_BDP_DEPTH_MODE )
            picture_control_set_ptr->parent_pcs_ptr->average_qp = (uint8_t)picture_control_set_ptr->parent_pcs_ptr->picture_qp;
        }


        // Derive MD parameters
        SetMdSettings( // HT Done
            sequence_control_set_ptr,
            picture_control_set_ptr);

        // Post the results to the MD processes
        EbGetEmptyObject(
            context_ptr->modeDecisionConfigurationOutputFifoPtr,
            &encDecTasksWrapperPtr);

        encDecTasksPtr = (EncDecTasks_t*)encDecTasksWrapperPtr->objectPtr;
        encDecTasksPtr->pictureControlSetWrapperPtr = rateControlResultsPtr->pictureControlSetWrapperPtr;
        encDecTasksPtr->inputType = ENCDEC_TASKS_MDC_INPUT;

        // Post the Full Results Object
        EbPostFullObject(encDecTasksWrapperPtr);

        // Release Rate Control Results
        EbReleaseObject(rateControlResultsWrapperPtr);

    }

    return EB_NULL;
}
