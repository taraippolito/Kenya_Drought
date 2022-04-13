def surf_texture(row):
    sand = row["sand_1"]
    silt = row['silt_1']
    clay = row['clay_1']
    texture = ""

    if (sand >= .85 and silt <= .15 and clay <= .10):  # sand
        texture = "sand"
    elif (sand >= .70 and sand <= .90 and silt <= .30 and clay <= .15):  # loamy sand
        texture = "loamy sand"
    elif (sand >= .43 and sand <= .80 and silt <= .50 and clay <= .20):  # sandy loam
        texture = "sandy loam"
    elif (sand >= .23 and sand <= .52 and silt >= .28 and silt <= .50 and clay >= .7 and clay <= .27):  # loam
        texture = "loam"
    elif (sand <= .50 and silt >= .50 and silt <= .88 and clay <= .27):  # silt loam
        texture = "silt loam"
    elif (sand <= .20 and silt >= .88 and clay <= .12):  # silt
        texture = "silt"
    elif (sand >= .45 and sand <= .80 and silt <= .28 and clay >= .20 and clay <= .35):  # sandy clay loam
        texture = "sandy clay loam"
    elif (sand >= .20 and sand <= .45 and silt >= .15 and silt <= .53 and clay >= .27 and clay <= .40):  # clay loam
        texture = "clay loam"
    elif (sand <= .20 and silt >= .40 and silt <= .73 and clay >= .27 and clay <= .40):  # silty clay loam
        texture = "silty clay loam"
    elif (sand >= .45 and sand <= .65 and silt <= .20 and clay >= .35 and clay <= .45):  # sandy clay
        texture = "sandy clay"
    elif (sand <= .20 and silt >= .40 and silt <= .60 and clay >= .40 and clay <= .60):  # silty clay
        texture = "silty clay"
    elif (sand <= .45 and silt <= .40 and clay >= .40):  # clay
        texture = "clay"

    return texture


def top_texture(row):
    sand = (row["sand_1"] * (5 / 15) + row["sand_2"] * (10 / 15))
    silt = (row["silt_1"] * (5 / 15) + row["silt_2"] * (10 / 15))
    clay = (row["clay_1"] * (5 / 15) + row["clay_2"] * (10 / 15))
    texture = ""

    if (sand >= .85 and silt <= .15 and clay <= .10):  # sand
        texture = "sand"
    elif (sand >= .70 and sand <= .90 and silt <= .30 and clay <= .15):  # loamy sand
        texture = "loamy sand"
    elif (sand >= .43 and sand <= .80 and silt <= .50 and clay <= .20):  # sandy loam
        texture = "sandy loam"
    elif (sand >= .23 and sand <= .52 and silt >= .28 and silt <= .50 and clay >= .7 and clay <= .27):  # loam
        texture = "loam"
    elif (sand <= .50 and silt >= .50 and silt <= .88 and clay <= .27):  # silt loam
        texture = "silt loam"
    elif (sand <= .20 and silt >= .88 and clay <= .12):  # silt
        texture = "silt"
    elif (sand >= .45 and sand <= .80 and silt <= .28 and clay >= .20 and clay <= .35):  # sandy clay loam
        texture = "sandy clay loam"
    elif (sand >= .20 and sand <= .45 and silt >= .15 and silt <= .53 and clay >= .27 and clay <= .40):  # clay loam
        texture = "clay loam"
    elif (sand <= .20 and silt >= .40 and silt <= .73 and clay >= .27 and clay <= .40):  # silty clay loam
        texture = "silty clay loam"
    elif (sand >= .45 and sand <= .65 and silt <= .20 and clay >= .35 and clay <= .45):  # sandy clay
        texture = "sandy clay"
    elif (sand <= .20 and silt >= .40 and silt <= .60 and clay >= .40 and clay <= .60):  # silty clay
        texture = "silty clay"
    elif (sand <= .45 and silt <= .40 and clay >= .40):  # clay
        texture = "clay"

    return texture


def awc1(row):
    sand = row["sand_1"]/1000
    silt = row["silt_1"]/1000
    clay = row["clay_1"]/1000
    om = 2.17 * row["soc_1"] / 1000
    gravel = row["cfvo_1"] / 1000

    WP1 = -0.024 * (sand) + 0.487 * (clay) + 0.006 * (om) + 0.005 * (sand) * (om) - 0.013 * (clay) * (om) + 0.068 * (
        sand) * (clay) + 0.031
    WP = WP1 + (0.14 * WP1 - 0.02)
    FC1 = -0.251 * (sand) + 0.195 * (clay) + 0.011 * (om) + 0.006 * (sand) * (om) - 0.027 * (clay) * (om) + 0.452 * (
        sand) * (clay) + 0.299
    FC = FC1 + (1.283 * FC1 * FC1 - 0.374 * FC1 - 0.015)
    LAWC = (FC - WP) * (1 - (gravel))

    return LAWC


def awc2(row):
    sand = row["sand_2"]/1000
    silt = row["silt_2"]/1000
    clay = row["clay_2"]/1000
    om = 2.17 * row["soc_2"] / 1000
    gravel = row["cfvo_2"] / 1000

    WP1 = -0.024 * (sand) + 0.487 * (clay) + 0.006 * (om) + 0.005 * (sand) * (om) - 0.013 * (clay) * (om) + 0.068 * (
        sand) * (clay) + 0.031
    WP = WP1 + (0.14 * WP1 - 0.02)
    FC1 = -0.251 * (sand) + 0.195 * (clay) + 0.011 * (om) + 0.006 * (sand) * (om) - 0.027 * (clay) * (om) + 0.452 * (
        sand) * (clay) + 0.299
    FC = FC1 + (1.283 * FC1 * FC1 - 0.374 * FC1 - 0.015)
    LAWC = (FC - WP) * (1 - (gravel))

    return LAWC


def awc3(row):
    sand = row["sand_3"]/1000
    silt = row["silt_3"]/1000
    clay = row["clay_3"]/1000
    om = 2.17 * row["soc_3"] / 1000
    gravel = row["cfvo_3"] / 1000

    WP1 = -0.024 * (sand) + 0.487 * (clay) + 0.006 * (om) + 0.005 * (sand) * (om) - 0.013 * (clay) * (om) + 0.068 * (
        sand) * (clay) + 0.031
    WP = WP1 + (0.14 * WP1 - 0.02)
    FC1 = -0.251 * (sand) + 0.195 * (clay) + 0.011 * (om) + 0.006 * (sand) * (om) - 0.027 * (clay) * (om) + 0.452 * (
        sand) * (clay) + 0.299
    FC = FC1 + (1.283 * FC1 * FC1 - 0.374 * FC1 - 0.015)
    LAWC = (FC - WP) * (1 - (gravel))

    return LAWC


def awc4(row):
    sand = row["sand_4"]/1000
    silt = row["silt_4"]/1000
    clay = row["clay_4"]/1000
    om = 2.17 * row["soc_4"] / 1000
    gravel = row["cfvo_4"] / 1000

    WP1 = -0.024 * (sand) + 0.487 * (clay) + 0.006 * (om) + 0.005 * (sand) * (om) - 0.013 * (clay) * (om) + 0.068 * (
        sand) * (clay) + 0.031
    WP = WP1 + (0.14 * WP1 - 0.02)
    FC1 = -0.251 * (sand) + 0.195 * (clay) + 0.011 * (om) + 0.006 * (sand) * (om) - 0.027 * (clay) * (om) + 0.452 * (
        sand) * (clay) + 0.299
    FC = FC1 + (1.283 * FC1 * FC1 - 0.374 * FC1 - 0.015)
    LAWC = (FC - WP) * (1 - (gravel))

    return LAWC


def awc5(row):
    sand = row["sand_5"]/1000
    silt = row["silt_5"]/1000
    clay = row["clay_5"]/1000
    om = 2.17 * row["soc_5"] / 1000
    gravel = row["cfvo_5"] / 1000

    WP1 = -0.024 * (sand) + 0.487 * (clay) + 0.006 * (om) + 0.005 * (sand) * (om) - 0.013 * (clay) * (om) + 0.068 * (
        sand) * (clay) + 0.031
    WP = WP1 + (0.14 * WP1 - 0.02)
    FC1 = -0.251 * (sand) + 0.195 * (clay) + 0.011 * (om) + 0.006 * (sand) * (om) - 0.027 * (clay) * (om) + 0.452 * (
        sand) * (clay) + 0.299
    FC = FC1 + (1.283 * FC1 * FC1 - 0.374 * FC1 - 0.015)
    LAWC = (FC - WP) * (1 - (gravel))

    return LAWC


# function to calculate erosion risk LCC class
# calculates K factor and combines with slope to determine risk
def erosion_risk_LCC(row):
    texture = row["surface_texture"]
    slope = row['slope']
    # first calculate k factor based on texture
    list1 = ["loam", "silt loam", "silt"]
    if list1.count(texture) > 0:
        k_fac = .33
    else:
        k_fac = .31

    # then calculate lcc class of erosion risk combining k factor and slope
    if k_fac > .32 and slope <= 2:
        erosion_risk_lcc = 1
    elif k_fac > .32 and (slope > 2 and slope <= 5):
        erosion_risk_lcc = 2
    elif k_fac > .32 and (slope > 5 and slope <= 10):
        erosion_risk_lcc = 3
    elif k_fac > .32 and (slope > 10 and slope <= 15):
        erosion_risk_lcc = 4
    elif k_fac > .32 and (slope > 15 and slope <= 30):
        erosion_risk_lcc = 6
    elif k_fac > .32 and (slope > 30 and slope <= 60):
        erosion_risk_lcc = 7
    elif k_fac > .32 and slope > 60:
        erosion_risk_lcc = 8
    elif k_fac <= 0.32 and slope <= 5:
        erosion_risk_lcc = 1
    elif k_fac <= 0.32 and (slope > 5 and slope <= 10):
        erosion_risk_lcc = 2
    elif k_fac <= 0.32 and (slope > 10 and slope <= 15):
        erosion_risk_lcc = 3
    elif k_fac <= 0.32 and (slope > 15 and slope <= 30):
        erosion_risk_lcc = 4
    elif k_fac <= 0.32 and (slope > 30 and slope <= 60):
        erosion_risk_lcc = 7
    elif k_fac <= 0.32 and slope > 60:
        erosion_risk_lcc = 8
    else:
        erosion_risk_lcc = 1

    return erosion_risk_lcc


# function to calculate surface soil texture LCC class
# takes in texture and outputs LCC class
def soil_tex_LCC(row):
    texture = row["top_texture"]
    list1 = ["sandy loam", "silt loam", "loam", "silt", "sandy clay loam", "silty clay loam", "clay loam"]
    list2 = ["sand", "loamy sand", "sandy clay", "silty clay"]

    if list1.count(texture) > 0:
        surface_st_lcc = 1
    elif list2.count(texture) > 0:
        surface_st_lcc = 2
    elif texture == "clay":
        surface_st_lcc = 3
    else:
        surface_st_lcc = 1

    return surface_st_lcc


#  function to calculate soil depth LCC class
#  takes in phases and outputs LCC class
def depth_LCC(row):
    depth = row["depth_to_bedrock"]
    if depth < 100 and depth >= 70:
        depth_LCC = 3
    elif depth < 70 and depth >= 50:
        depth_LCC = 4
    elif depth < 50 and depth >= 20:
        depth_LCC = 6
    elif depth < 20:
        depth_LCC = 8
    else:
        depth_LCC = 1

    return depth_LCC


# function to calculate surface stoniness LCC class
# takes in coarse fragments and outputs LCC class
def stoniness_LCC(row):
    stone = row['coarse_fra_1']

    if stone < .1:
        stone_LCC = 1
    elif stone >= .1 and stone < 3:
        stone_LCC = 2
    elif stone >= 3 and stone < 15:
        stone_LCC = 3
    elif stone >= 15 and stone < 50:
        stone_LCC = 5
    elif stone >= 50 and stone < 90:
        stone_LCC = 7
    elif stone >= 90:
        stone_LCC = 8

    return stone_LCC


# function to calculate permeability LCC class
# take in texture, calculate intermediary values, outpus LCC class
# based on landPKs alculation of permeability
def perm_LCC(row):
    import math
    KS_list = []
    for i in range(1, 6):
        sa = "sand_"
        si = "silt_"
        c = "clay_"
        oc = "org_c_cont_"
        cf = "coarse_fra_"
        sand = row[sa + str(i)]
        silt = row[si + str(i)]
        clay = row[c + str(i)]
        om = 2.17 * row[oc + str(i)] / 100
        gravel = row[cf + str(i)] / 100

        WP1 = -0.024 * (sand) + 0.487 * (clay) + 0.006 * (om) + 0.005 * (sand) * (om) - 0.013 * (clay) * (
            om) + 0.068 * (sand) * (clay) + 0.031
        WP = WP1 + (0.14 * WP1 - 0.02)
        FC1 = -0.251 * (sand) + 0.195 * (clay) + 0.011 * (om) + 0.006 * (sand) * (om) - 0.027 * (clay) * (
            om) + 0.452 * (sand) * (clay) + 0.299
        FC = FC1 + (1.283 * FC1 * FC1 - 0.374 * FC1 - 0.015)
        LAWC = (FC - WP) * (1 - (gravel))
        SFC1 = 0.278 * (sand) + 0.034 * (clay) + 0.022 * (om) - 0.018 * (sand) * (om) - 0.027 * (clay) * (
            om) - 0.584 * (sand) * (clay) + 0.078
        SFC = SFC1 + ((0.636 * SFC1) - 0.107)
        SASC = FC + SFC - 0.097 * (sand) + 0.043
        MSD = (1 - SASC) * 2.65
        BSD = (1 - (gravel)) / (1 - (gravel) * (1 - 1.5 * (MSD / 2.65)))
        L = (math.log(WP) - math.log(FC)) / (math.log(1500) - math.log(33))
        KS = 1930 * (SASC) ** (3 - L) * BSD

        KS_list.append(KS)

    perm_lcc = 1
    if min(KS_list) < 1.5:
        perm_lcc = 3
    elif min(KS_list) >= 1.5 and KS <= 5:
        perm_lcc = 2
    elif min(KS_list) > 5:
        perm_lcc = 1

    return perm_lcc


# function to calculate lime requirement LCC class
# takes in T_PH_H2O and outputs LCC class
def lime_LCC(row):
    ph_val = row["ph_1"] * (5 / 30) + row["ph_2"] * (10 / 30) + row["ph_3"] * (15 / 30)
    if ph_val >= 5.5 and ph_val <= 7.2:
        lime_lcc = 1
    elif (ph_val < 5.5 and ph_val >= 4.5) or (ph_val > 7.2 and ph_val <= 8.4):
        lime_lcc = 3
    elif ph_val < 4.5 or ph_val > 8.5:
        lime_lcc = 4

    return lime_lcc


def LAWC_LCC(row):
    AWC = row['awc_total']
    if AWC > 18:
        awc_lcc = 1
    elif AWC > 12 and AWC <= 18:
        awc_lcc = 2
    elif AWC > 6 and AWC <= 12:
        awc_lcc = 3
    elif AWC <= 6:
        awc_lcc = 4
    else:
        awc_lcc = 1

    return awc_lcc


# calculate final LCC and limitations
def calc_LCC(row):
    lcc_final = max(row["erosion_risk_LCC"], row["soil_tex_LCC"], row["stoniness_LCC"], row["lime_LCC"],
                    row["perm_LCC"], row['awc_LCC'])

    return lcc_final


def calc_LCC_no_awc(row):
    lcc_final = max(row["erosion_risk_LCC"], row["soil_tex_LCC"], row["stoniness_LCC"], row["lime_LCC"],
                    row["perm_LCC"])

    return lcc_final


def limitation(row):
    LCC_final = row['LCC_final']
    erosion = row['erosion_risk_LCC']
    tex = row['soil_tex_LCC']
    stone = row['stoniness_LCC']
    perm = row['perm_LCC']
    awc = row['awc_LCC']
    lime = row['lime_LCC']

    lim = str(LCC_final)
    if LCC_final == erosion:
        lim += "e "
    if LCC_final == tex:
        lim += "s-t "
    if LCC_final == stone:
        lim += "s-r "
    if LCC_final == perm:
        lim += "w-p "
    if awc == LCC_final:
        lim += "s-a "
    if LCC_final == lime:
        lim += "s-l "

    return (lim)


def limitation_noawc(row):
    LCC_final = row['LCC_final_noawc']
    erosion = row['erosion_risk_LCC']
    tex = row['soil_tex_LCC']
    stone = row['stoniness_LCC']
    perm = row['perm_LCC']
    lime = row['lime_LCC']

    lim = str(LCC_final)
    if LCC_final == erosion:
        lim += "e "
    if LCC_final == tex:
        lim += "s-t "
    if LCC_final == stone:
        lim += "s-r "
    if LCC_final == perm:
        lim += "w-p "
    if LCC_final == lime:
        lim += "s-l "

    return (lim)