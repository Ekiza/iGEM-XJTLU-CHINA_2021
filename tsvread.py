import pandas as pd
import numpy as np
from scipy import stats
import openpyxl


def phage_ratio(p2, c2):
    moi = p2 / c2
    x = np.arange(1, 300, 1)
    pdf = stats.poisson.pmf(x, moi)
    pdf1 = pdf * c2
    pdfall = dict(zip(x, pdf1))
    for inf in x:
        tempnum = pdfall[inf]
        pdfall[inf] = round(tempnum)
        if pdfall[inf] <= .01 * c2:
            del pdfall[inf]
    above20 = []
    remainkey = list(dict.keys(pdfall))
    for inf1 in remainkey:
        if inf1 >= 20:
            above20.append(pdfall[inf1])
            del pdfall[inf1]
            pdfall.update({20: sum(above20)})
    if not pdfall:
        pdfall.update({20: c2})
    print(pdfall)
    return pdfall


def execute(phn, cn):
    dictpro = phage_ratio(phn, cn)
    infectionnum = dict.keys(dictpro)
    LuxR = []
    lysis = []
    for n in infectionnum:
        tempdata = express(n, dictpro[n])
        LuxR.append(tempdata[0])
        lysis.append(tempdata[1])
    print([phn, cn])
    return [sum(LuxR), max(lysis)]


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


def express(infnum, cellnum):
    holinpass = (infnum - 1) * (-124.65) + 12341
    # holinpass = 12341
    existing = np.array([0, 1, 2, 3, 4, 5])
    datanumber = find_nearest(existing, (infnum - 1))
    strfile = "real." + str(datanumber) + ".tsv"
    # if datanumber == 16 or datanumber == 11:
    #     train = pd.read_csv(strfile, sep='\t', header=0, usecols=[0, 1, 2])
    #     train1 = train[train['species'] == "LuxR"]
    #     train2 = train1['protein'].tail(1)
    #     train3 = train1['time'].tail(1)
    #     lytime = train3.iloc[0]
    #     protein = train2.iloc[0]
    #     return [protein, lytime]
    # else:
    train = pd.read_csv(strfile, sep='\t', header=0, usecols=[0, 1, 2])
    train1 = train[train['species'] == "holin"]
    train2 = train1[train1['protein'] >= holinpass].head(1)
    train3 = train1[train1['protein'] < holinpass].tail(1)
    train4 = train[train['species'] == "LuxR"]
    pro = train2.iloc[0].iat[2] - train3.iloc[0].iat[2]
    time = train2.iloc[0].iat[0] - train3.iloc[0].iat[0]
    slope = pro / time
    lysistime = (holinpass - train3.iloc[0].iat[2]) / slope + train3.iloc[0].iat[0]
    train5 = train4[train4['time'] >= lysistime].head(1)
    train6 = train4[train4['time'] < lysistime].tail(1)
    pro1 = train5.iloc[0].iat[2] - train6.iloc[0].iat[2]
    time1 = train5.iloc[0].iat[0] - train6.iloc[0].iat[0]
    slope1 = pro1 / time1
    output = slope1 * (lysistime - train6.iloc[0].iat[0]) + train6.iloc[0].iat[2]
    outputfinal = output * cellnum
    return [outputfinal, lysistime]


if __name__ == '__main__':
    timemap = np.zeros(shape=(100, 100))
    LuxRmap = np.zeros(shape=(100, 100))
    for p in range(200, 20200, 200):
        for c in range(20, 2020, 20):
            result = execute(p, c)
            # print(result)
            pm = (p / 200) - 1
            cm = (c / 20) - 1
            LuxRmap[int(pm)][int(cm)] = result[0]
            timemap[int(pm)][int(cm)] = result[1]
    LuxRdata = pd.DataFrame(LuxRmap)
    timedata = pd.DataFrame(timemap)
    writer1 = pd.ExcelWriter('LuxRmapf1.xlsx')
    writer2 = pd.ExcelWriter('timemapf1.xlsx')
    LuxRdata.to_excel(writer1, 'page_1', float_format='%5.f')
    timedata.to_excel(writer2, 'page_1', float_format='%5.f')
    np.save("LuxRmaptestf1.npy", LuxRmap)
    np.save("timemaptestf1.npy", timemap)
    writer1.save()
    writer2.save()
