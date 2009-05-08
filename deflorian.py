from xsec import datapoint

class Thing(object):
    pass

minus = Thing()

minus.mrst = [
datapoint(x=0.05, xlow=0.0, binwidth=0.1, y=0.1057E-01, sys=0.2928E-02),
datapoint(x=0.15, xlow=0.1, binwidth=0.1, y=0.1296E+00, sys=0.4335E-02),
datapoint(x=0.25, xlow=0.2, binwidth=0.1, y=0.7599E-01, sys=0.1943E-02),
datapoint(x=0.35, xlow=0.3, binwidth=0.1, y=0.3122E-01, sys=0.4656E-03),
datapoint(x=0.45, xlow=0.4, binwidth=0.1, y=0.1308E-01, sys=0.4491E-03),
datapoint(x=0.55, xlow=0.5, binwidth=0.1, y=0.5818E-02, sys=0.1392E-03),
datapoint(x=0.65, xlow=0.6, binwidth=0.1, y=0.2628E-02, sys=0.8203E-04),
datapoint(x=0.75, xlow=0.7, binwidth=0.1, y=0.8632E-03, sys=0.1523E-03),
datapoint(x=0.85, xlow=0.8, binwidth=0.1, y=0.3275E-03, sys=0.1294E-04),
datapoint(x=0.95, xlow=0.9, binwidth=0.1, y=0.7719E-04, sys=0.3246E-05)
]

minus.dssv = [
datapoint(x=0.05, xlow=0.0, binwidth=0.1, y= 0.4330E-04, sys=0.8099E-05),
datapoint(x=0.15, xlow=0.1, binwidth=0.1, y= 0.4508E-03, sys=0.9698E-05),
datapoint(x=0.25, xlow=0.2, binwidth=0.1, y= 0.1327E-03, sys=0.6698E-05),
datapoint(x=0.35, xlow=0.3, binwidth=0.1, y= 0.4078E-04, sys=0.2158E-05),
datapoint(x=0.45, xlow=0.4, binwidth=0.1, y= 0.1024E-04, sys=0.1669E-05),
datapoint(x=0.55, xlow=0.5, binwidth=0.1, y= 0.3412E-05, sys=0.4987E-06),
datapoint(x=0.65, xlow=0.6, binwidth=0.1, y=-0.1323E-05, sys=0.1332E-06),
datapoint(x=0.75, xlow=0.7, binwidth=0.1, y=-0.1919E-05, sys=0.7114E-07),
datapoint(x=0.85, xlow=0.8, binwidth=0.1, y=-0.1494E-05, sys=0.5914E-07),
datapoint(x=0.95, xlow=0.9, binwidth=0.1, y=-0.5313E-06, sys=0.1878E-07)
]

minus.std = [
datapoint(x=0.05, xlow=0.0, binwidth=0.1, y= 0.9511E-04, sys=0.8099E-05),
datapoint(x=0.15, xlow=0.1, binwidth=0.1, y= 0.2300E-02, sys=0.9698E-05),
datapoint(x=0.25, xlow=0.2, binwidth=0.1, y= 0.1038E-02, sys=0.6698E-05),
datapoint(x=0.35, xlow=0.3, binwidth=0.1, y= 0.3960E-03, sys=0.2158E-05),
datapoint(x=0.45, xlow=0.4, binwidth=0.1, y= 0.1501E-03, sys=0.1669E-05),
datapoint(x=0.55, xlow=0.5, binwidth=0.1, y= 0.6809E-04, sys=0.4987E-06),
datapoint(x=0.65, xlow=0.6, binwidth=0.1, y= 0.2304E-04, sys=0.1332E-06),
datapoint(x=0.75, xlow=0.7, binwidth=0.1, y= 0.6188E-05, sys=0.7114E-07),
datapoint(x=0.85, xlow=0.8, binwidth=0.1, y= 0.2943E-06, sys=0.5914E-07),
datapoint(x=0.95, xlow=0.9, binwidth=0.1, y=-0.3824E-06, sys=0.1878E-07)
]

minus.max = [
datapoint(x=0.05, xlow=0.0, binwidth=0.1, y= 0.7335E-03, sys=0.2076E-03),
datapoint(x=0.15, xlow=0.1, binwidth=0.1, y= 0.1134E-01, sys=0.3260E-03),
datapoint(x=0.25, xlow=0.2, binwidth=0.1, y= 0.5664E-02, sys=0.1576E-03),
datapoint(x=0.35, xlow=0.3, binwidth=0.1, y= 0.2351E-02, sys=0.5031E-04),
datapoint(x=0.45, xlow=0.4, binwidth=0.1, y= 0.8677E-03, sys=0.5013E-04),
datapoint(x=0.55, xlow=0.5, binwidth=0.1, y= 0.4260E-03, sys=0.1616E-04),
datapoint(x=0.65, xlow=0.6, binwidth=0.1, y= 0.1661E-03, sys=0.6103E-05),
datapoint(x=0.75, xlow=0.7, binwidth=0.1, y= 0.5497E-04, sys=0.2429E-05),
datapoint(x=0.85, xlow=0.8, binwidth=0.1, y= 0.1234E-04, sys=0.5243E-06),
datapoint(x=0.95, xlow=0.9, binwidth=0.1, y= 0.1443E-05, sys=0.1229E-06)
]

minus.gsc = [
datapoint(x=0.05, xlow=0.0, binwidth=0.1, y=-0.7625E-05, sys=0.1296E-05),
datapoint(x=0.15, xlow=0.1, binwidth=0.1, y=-0.1982E-03, sys=0.3064E-05),
datapoint(x=0.25, xlow=0.2, binwidth=0.1, y=-0.1129E-03, sys=0.8317E-06),
datapoint(x=0.35, xlow=0.3, binwidth=0.1, y=-0.5069E-04, sys=0.2994E-05),
datapoint(x=0.45, xlow=0.4, binwidth=0.1, y=-0.2539E-04, sys=0.6922E-06),
datapoint(x=0.55, xlow=0.5, binwidth=0.1, y=-0.1389E-04, sys=0.2118E-06),
datapoint(x=0.65, xlow=0.6, binwidth=0.1, y=-0.7587E-05, sys=0.1711E-06),
datapoint(x=0.75, xlow=0.7, binwidth=0.1, y=-0.3874E-05, sys=0.1466E-06),
datapoint(x=0.85, xlow=0.8, binwidth=0.1, y=-0.1869E-05, sys=0.7134E-07),
datapoint(x=0.95, xlow=0.9, binwidth=0.1, y=-0.5389E-06, sys=0.1785E-07)
]

plus = Thing()

plus.mrst = [
datapoint(x=0.05, xlow=0.0, binwidth=0.1, y=0.7250E-02, sys=0.2121E-02),
datapoint(x=0.15, xlow=0.1, binwidth=0.1, y=0.1342E+00, sys=0.3744E-02),
datapoint(x=0.25, xlow=0.2, binwidth=0.1, y=0.8116E-01, sys=0.1814E-02),
datapoint(x=0.35, xlow=0.3, binwidth=0.1, y=0.3410E-01, sys=0.5215E-03),
datapoint(x=0.45, xlow=0.4, binwidth=0.1, y=0.1429E-01, sys=0.4677E-03),
datapoint(x=0.55, xlow=0.5, binwidth=0.1, y=0.6504E-02, sys=0.2201E-03),
datapoint(x=0.65, xlow=0.6, binwidth=0.1, y=0.3151E-02, sys=0.8984E-04),
datapoint(x=0.75, xlow=0.7, binwidth=0.1, y=0.1283E-02, sys=0.4874E-04),
datapoint(x=0.85, xlow=0.8, binwidth=0.1, y=0.4691E-03, sys=0.1564E-04),
datapoint(x=0.95, xlow=0.9, binwidth=0.1, y=0.1156E-03, sys=0.5343E-05)
]

plus.dssv = [
datapoint(x=0.05, xlow=0.0, binwidth=0.1, y=0.6100E-04, sys=0.9855E-05),
datapoint(x=0.15, xlow=0.1, binwidth=0.1, y=0.6013E-03, sys=0.3802E-04),
datapoint(x=0.25, xlow=0.2, binwidth=0.1, y=0.3155E-03, sys=0.1138E-04),
datapoint(x=0.35, xlow=0.3, binwidth=0.1, y=0.1490E-03, sys=0.4420E-05),
datapoint(x=0.45, xlow=0.4, binwidth=0.1, y=0.7216E-04, sys=0.3986E-05),
datapoint(x=0.55, xlow=0.5, binwidth=0.1, y=0.3759E-04, sys=0.1239E-05),
datapoint(x=0.65, xlow=0.6, binwidth=0.1, y=0.2106E-04, sys=0.7834E-06),
datapoint(x=0.75, xlow=0.7, binwidth=0.1, y=0.1027E-04, sys=0.4140E-06),
datapoint(x=0.85, xlow=0.8, binwidth=0.1, y=0.4880E-05, sys=0.1765E-06),
datapoint(x=0.95, xlow=0.9, binwidth=0.1, y=0.1766E-05, sys=0.4686E-07)
]

plus.std = [
datapoint(x=0.05, xlow=0.0, binwidth=0.1, y=0.2195E-03, sys=0.5080E-04),
datapoint(x=0.15, xlow=0.1, binwidth=0.1, y=0.2784E-02, sys=0.9807E-04),
datapoint(x=0.25, xlow=0.2, binwidth=0.1, y=0.1370E-02, sys=0.3072E-04),
datapoint(x=0.35, xlow=0.3, binwidth=0.1, y=0.5996E-03, sys=0.1200E-04),
datapoint(x=0.45, xlow=0.4, binwidth=0.1, y=0.2724E-03, sys=0.1404E-04),
datapoint(x=0.55, xlow=0.5, binwidth=0.1, y=0.1328E-03, sys=0.4531E-05),
datapoint(x=0.65, xlow=0.6, binwidth=0.1, y=0.7032E-04, sys=0.3513E-05),
datapoint(x=0.75, xlow=0.7, binwidth=0.1, y=0.3053E-04, sys=0.1372E-05),
datapoint(x=0.85, xlow=0.8, binwidth=0.1, y=0.1156E-04, sys=0.6196E-06),
datapoint(x=0.95, xlow=0.9, binwidth=0.1, y=0.3644E-05, sys=0.1419E-06)
]

plus.max = [
datapoint(x=0.05, xlow=0.0, binwidth=0.1, y=0.9067E-03, sys=0.9024E-03),
datapoint(x=0.15, xlow=0.1, binwidth=0.1, y=0.1151E-01, sys=0.3106E-03),
datapoint(x=0.25, xlow=0.2, binwidth=0.1, y=0.6491E-02, sys=0.1460E-03),
datapoint(x=0.35, xlow=0.3, binwidth=0.1, y=0.2863E-02, sys=0.5373E-04),
datapoint(x=0.45, xlow=0.4, binwidth=0.1, y=0.1196E-02, sys=0.4938E-04),
datapoint(x=0.55, xlow=0.5, binwidth=0.1, y=0.5667E-03, sys=0.1677E-04),
datapoint(x=0.65, xlow=0.6, binwidth=0.1, y=0.2576E-03, sys=0.9428E-05),
datapoint(x=0.75, xlow=0.7, binwidth=0.1, y=0.1174E-03, sys=0.8087E-05),
datapoint(x=0.85, xlow=0.8, binwidth=0.1, y=0.3823E-04, sys=0.1976E-05),
datapoint(x=0.95, xlow=0.9, binwidth=0.1, y=0.9182E-05, sys=0.4135E-06)
]

plus.gsc = [
datapoint(x=0.05, xlow=0.0, binwidth=0.1, y=-0.4189E-05, sys=0.2336E-05),
datapoint(x=0.15, xlow=0.1, binwidth=0.1, y=-0.9303E-04, sys=0.2856E-05),
datapoint(x=0.25, xlow=0.2, binwidth=0.1, y=-0.2641E-04, sys=0.1604E-05),
datapoint(x=0.35, xlow=0.3, binwidth=0.1, y=-0.9422E-06, sys=0.1047E-05),
datapoint(x=0.45, xlow=0.4, binwidth=0.1, y= 0.1967E-05, sys=0.4163E-06),
datapoint(x=0.55, xlow=0.5, binwidth=0.1, y= 0.3812E-05, sys=0.3576E-06),
datapoint(x=0.65, xlow=0.6, binwidth=0.1, y= 0.2589E-05, sys=0.2137E-06),
datapoint(x=0.75, xlow=0.7, binwidth=0.1, y= 0.1562E-05, sys=0.1759E-06),
datapoint(x=0.85, xlow=0.8, binwidth=0.1, y= 0.1040E-05, sys=0.6509E-07),
datapoint(x=0.95, xlow=0.9, binwidth=0.1, y= 0.5663E-06, sys=0.2033E-07)
]
