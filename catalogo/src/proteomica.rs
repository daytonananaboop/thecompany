use crate::formato::{decimal, miles, porcentaje};
use crate::pizarra::{Marca, Narracion, Pizarra};
use crate::{Algoritmo, Bloque, Campo, Entrada, Paso, Punto, Salida, Tipo};
use nucleo::ecoli;
use nucleo::proteomica::{
    self, Busqueda, Digestion, Espectro, FDR_ESTANDAR, Indice, PeptidoDigerido, PeptidoIndexado,
    Psm,
};
use rayon::prelude::*;
use std::borrow::Cow;
use std::sync::OnceLock;
use std::time::Instant;

const ESPECTROS_DE_EJEMPLO: &str = "BEGIN IONS
TITLE=corrida_1.1200.1200.3
PEPMASS=609.99149
CHARGE=3+
RTINSECONDS=600
119.3243 2.9
134.1238 6.9
140.3129 29.0
187.1304 16.4
206.9171 1.6
318.1714 90.3
320.1206 8.0
417.2481 69.7
421.1767 43.7
643.4132 90.2
649.2867 42.4
707.2129 7.4
709.2402 8.7
748.3496 12.1
780.4736 73.8
837.4959 90.2
867.2114 12.5
876.4464 21.4
952.5252 53.8
991.4740 29.5
1080.6172 89.2
1096.1411 3.9
1101.1668 12.1
1164.2591 2.7
1177.1972 2.9
1179.6827 75.9
1200.2575 2.1
1205.8314 3.3
1220.8355 2.7
1258.3884 3.6
1308.7206 42.4
1407.4509 14.9
1407.7967 48.8
1411.7175 26.6
1434.3267 2.5
1476.3358 3.2
1480.3213 1.6
1508.8421 72.9
1512.9209 17.4
1603.5195 4.3
1624.8256 51.3
1655.0375 13.2
1680.8865 83.3
1681.7941 4.7
1704.6516 13.2
END IONS
BEGIN IONS
TITLE=corrida_1.1237.1237.2
PEPMASS=827.92155
CHARGE=2+
RTINSECONDS=623
123.2970 22.8
133.7385 3.4
175.1182 96.7
269.8764 34.4
276.1643 27.7
289.0316 9.5
334.1746 52.7
372.8573 1.6
376.4246 7.8
400.1687 6.8
439.2296 45.3
474.3966 4.4
544.3228 54.4
552.3256 43.1
748.2344 13.1
780.4297 78.5
808.2852 12.8
847.1491 17.1
875.4190 12.3
895.4598 77.6
931.6835 20.5
988.5002 50.0
1003.8357 15.6
1010.4765 75.5
1022.4529 14.8
1111.5221 32.3
1153.9542 4.3
1186.5976 4.7
1216.6186 19.2
1321.6619 25.4
1328.7944 1.8
1356.3587 1.8
1377.0217 3.1
1379.6752 45.2
1452.4487 11.4
1484.7289 51.3
1504.1973 5.0
1541.7460 88.9
1546.4677 17.0
1623.6507 3.9
END IONS
BEGIN IONS
TITLE=corrida_1.1274.1274.2
PEPMASS=693.40908
CHARGE=2+
RTINSECONDS=646
161.6337 5.4
178.0756 28.1
303.1787 62.3
316.1919 46.4
316.3926 9.4
320.9234 4.2
374.2047 58.2
379.2325 2.5
384.3862 15.5
417.2324 17.0
445.2582 52.8
511.6106 5.9
516.2862 91.1
530.3210 53.2
595.8081 3.3
607.7373 15.4
609.9549 7.2
613.2923 9.7
614.5459 14.2
658.3028 24.9
658.3800 8.4
670.4577 1.7
710.7246 21.8
711.9410 5.7
728.4428 22.5
757.4488 25.0
839.5764 27.0
856.5001 38.5
870.5356 12.6
887.6418 3.7
910.3996 2.9
986.8773 6.6
996.5330 2.1
1012.6041 10.6
1070.6344 24.1
1083.6422 39.8
1087.5454 3.8
1109.8763 23.2
1131.2509 4.7
1153.3366 8.4
1201.6310 19.4
1204.5587 4.0
1211.7039 19.7
1220.2242 12.4
1284.7675 92.1
END IONS
BEGIN IONS
TITLE=corrida_1.1311.1311.2
PEPMASS=746.94633
CHARGE=2+
RTINSECONDS=669
113.6150 10.6
117.7675 19.2
147.1154 74.2
171.1187 24.3
186.5567 1.0
204.1348 33.8
317.2183 27.5
353.8924 1.9
398.9250 7.0
412.2566 16.3
449.2643 31.8
511.6356 1.8
513.5919 5.1
524.7781 4.4
530.1378 4.4
530.4047 5.0
563.9730 17.2
583.3179 49.5
629.3964 62.0
653.7096 3.5
680.3684 32.9
700.4316 33.5
720.0426 12.9
734.1620 4.6
776.6972 10.5
813.5187 75.4
851.2173 3.6
862.0269 8.8
864.4885 14.6
878.0783 3.3
897.8968 2.9
910.5741 97.1
973.4618 2.7
1027.0175 4.4
1063.5864 51.1
1112.7142 28.9
1128.8215 8.3
1158.6685 7.7
1176.6780 17.4
1180.7057 61.1
1202.6667 1.9
1251.7455 66.4
1303.0744 15.4
1322.7830 68.1
1421.8505 71.5
1433.3933 11.4
1482.8690 6.1
1484.7532 1.3
END IONS
BEGIN IONS
TITLE=corrida_1.1348.1348.2
PEPMASS=493.28822
CHARGE=2+
RTINSECONDS=692
125.6943 20.6
136.3091 3.8
157.0980 52.9
185.8617 39.6
225.7931 13.8
246.1817 96.4
255.1689 1.9
295.0602 25.3
333.8771 7.4
361.1272 23.1
370.2132 46.9
371.2942 11.6
390.9156 9.7
407.7854 6.5
470.5785 34.2
492.5345 2.7
503.2856 75.8
546.6171 2.1
554.3310 39.8
616.3665 47.6
658.8044 17.6
669.3549 22.5
694.8680 20.9
740.3924 49.8
769.4544 18.9
770.9782 5.7
799.8840 12.0
924.9442 3.9
928.5528 52.2
END IONS
BEGIN IONS
TITLE=corrida_1.1385.1385.2
PEPMASS=706.39141
CHARGE=2+
RTINSECONDS=715
158.8819 15.9
174.6338 19.6
175.9925 23.3
196.8010 5.4
225.5450 9.1
250.7225 26.3
289.8126 19.1
313.9078 13.7
319.0836 21.8
465.2628 1.8
574.0028 1.9
612.0340 1.2
745.3278 1.8
829.0052 26.8
834.6498 31.9
887.7158 11.2
899.2176 8.6
914.4594 8.9
915.0741 24.2
1012.2718 8.3
1023.5443 51.2
1048.2313 25.9
1083.6677 36.2
1085.1269 1.5
1195.8926 3.9
1284.4143 9.7
1287.4670 25.6
1294.0644 2.2
1314.2205 8.0
1359.3622 32.2
1374.0546 8.2
1383.0773 1.8
END IONS
BEGIN IONS
TITLE=corrida_1.1422.1422.2
PEPMASS=993.52497
CHARGE=2+
RTINSECONDS=738
144.7436 14.6
158.1451 72.7
169.3997 38.6
237.2698 11.2
312.3441 3.0
383.0302 42.0
461.6861 5.4
620.5580 5.1
652.8521 14.1
676.9456 41.2
770.3123 39.3
820.8129 5.9
850.5911 29.2
908.6117 68.3
966.6966 13.1
1064.2943 27.1
1066.0623 9.4
1072.0664 29.5
1177.8145 44.0
1227.2497 43.2
1249.6114 11.9
1309.5496 63.0
1330.6446 10.6
1398.4105 14.3
1464.4595 10.9
1476.0286 38.7
1604.1736 7.9
1728.2667 8.1
1846.0687 3.5
1887.1782 36.1
1904.0363 1.4
1946.5257 7.9
END IONS
BEGIN IONS
TITLE=corrida_1.1459.1459.2
PEPMASS=845.41668
CHARGE=2+
RTINSECONDS=761
155.3990 2.3
171.2838 11.4
175.1159 79.8
288.2039 44.2
300.1600 36.7
313.4565 6.7
416.2554 74.1
552.0312 13.3
576.2978 22.6
643.3448 54.5
660.1070 5.1
705.3351 35.3
752.2265 7.2
753.5978 15.7
811.9926 5.0
833.3983 56.7
932.4605 71.7
985.4951 29.6
1020.2602 5.4
1047.4911 65.8
1094.0970 7.5
1114.5451 20.0
1170.2585 16.8
1194.5509 87.0
1244.1234 24.6
1274.5773 53.5
1293.6235 46.7
1320.2471 2.4
1328.6293 1.6
1356.7695 8.3
1361.7507 2.0
1390.6797 50.0
1495.9565 10.1
1519.7208 52.2
1564.6811 16.4
1655.1571 1.0
1659.1994 24.2
1670.6495 14.1
END IONS
BEGIN IONS
TITLE=corrida_1.1496.1496.2
PEPMASS=762.45351
CHARGE=2+
RTINSECONDS=784
121.4991 9.3
181.6404 31.0
191.6729 26.7
194.9362 8.7
250.2701 10.9
264.0156 26.0
271.2611 1.8
287.2582 48.9
309.4387 23.7
317.7016 16.9
374.9146 8.8
409.7477 13.6
410.4360 10.3
455.2676 14.8
482.6307 18.2
492.0173 36.0
540.5235 2.3
565.0910 48.5
603.7074 38.8
615.5186 36.8
639.7234 13.3
718.9152 24.5
720.0734 3.3
868.7201 19.5
982.5645 17.5
1013.0627 15.1
1027.5178 32.1
1030.1967 38.8
1034.0095 36.7
1065.5816 9.6
1085.0548 33.8
1094.6671 23.2
1096.6267 10.1
1142.0072 60.7
1142.5753 27.6
1203.5179 1.0
1210.5000 5.3
1278.9463 4.4
1323.7622 62.0
1406.6546 21.0
1426.9968 48.2
1430.4369 8.9
1445.8198 17.3
END IONS
BEGIN IONS
TITLE=corrida_1.1533.1533.2
PEPMASS=879.96849
CHARGE=2+
RTINSECONDS=807
115.7971 10.7
165.4516 1.8
171.5220 4.9
337.0326 35.6
342.9745 60.8
366.9572 9.5
377.2446 4.8
379.1743 9.1
445.6696 4.4
458.8676 14.1
502.4124 18.6
532.0844 68.0
567.7487 7.9
574.9331 2.1
594.0083 16.0
638.5843 17.9
655.3156 2.1
734.3790 8.2
826.2665 58.0
840.7350 5.5
866.1327 1.2
924.5672 2.3
925.5109 22.5
976.0581 28.1
1039.0227 6.0
1040.6283 4.7
1068.7795 23.3
1089.2297 21.0
1099.1196 3.7
1123.9975 15.3
1126.5014 52.2
1286.8960 35.3
1342.1134 11.8
1355.0189 11.9
1368.9151 20.7
1411.9602 4.8
1429.4454 2.3
1439.5998 19.2
1525.7542 18.7
1529.2964 40.3
1550.3018 29.4
1578.8148 7.0
1618.5524 15.0
1621.6483 46.7
1635.9954 2.0
END IONS
BEGIN IONS
TITLE=corrida_1.1570.1570.2
PEPMASS=523.28631
CHARGE=2+
RTINSECONDS=830
167.4615 14.8
175.1136 61.7
181.0187 3.7
239.7803 17.2
262.1554 90.5
268.8644 12.2
317.3859 5.5
327.8186 5.1
374.6325 4.4
376.1968 42.8
379.7270 5.8
385.2147 9.0
407.3471 4.6
456.2349 12.5
456.2451 10.8
489.2741 42.3
521.0053 3.2
541.3437 8.1
557.6189 10.3
577.5547 11.7
585.1180 2.0
590.3211 71.9
663.4426 1.6
665.9341 1.6
670.3810 20.5
758.4173 43.7
784.4183 48.4
808.5598 1.8
810.0560 2.4
845.4445 50.9
871.4512 50.2
901.0084 2.4
932.4714 80.5
951.3717 18.1
970.3691 11.0
1009.3141 1.2
END IONS
BEGIN IONS
TITLE=corrida_1.1607.1607.2
PEPMASS=902.88844
CHARGE=2+
RTINSECONDS=853
140.9980 2.5
195.4689 11.0
213.6193 39.4
278.8949 8.0
303.5074 43.0
392.0386 14.6
430.1499 35.8
433.0927 34.1
486.1655 17.5
498.4380 8.7
519.7870 18.9
551.0803 68.0
576.4666 11.0
713.4145 14.2
717.1898 44.1
749.8434 20.8
774.0318 49.4
801.4185 19.4
818.3065 14.9
862.9303 30.3
954.7401 16.3
975.2346 2.6
1015.2604 3.7
1050.9168 5.9
1092.9265 30.7
1097.2823 31.1
1130.7982 32.2
1175.0826 47.3
1260.0068 10.0
1345.7646 8.2
1353.4458 20.2
1356.0556 13.0
1420.8461 8.7
1538.9935 6.0
1554.1693 4.0
1571.3723 12.8
1603.0757 2.9
END IONS
BEGIN IONS
TITLE=corrida_1.1644.1644.2
PEPMASS=542.82450
CHARGE=2+
RTINSECONDS=876
152.0982 1.4
185.1317 15.3
208.3165 1.6
256.1680 26.9
300.4160 5.8
307.3260 10.0
332.8339 14.0
379.0183 5.0
389.2374 78.6
391.3753 31.4
407.0654 16.1
432.3364 2.8
443.8021 54.2
513.0298 9.2
529.6759 2.6
556.9139 1.4
568.3422 16.1
630.5232 6.7
653.7999 22.1
687.4088 68.5
696.4032 18.1
709.2024 3.4
719.8331 4.0
758.4382 57.9
788.9035 21.4
829.4777 94.8
837.2449 4.5
935.6802 7.4
1013.5963 65.8
1039.5600 12.9
1044.5547 1.2
END IONS
BEGIN IONS
TITLE=corrida_1.1681.1681.2
PEPMASS=435.56843
CHARGE=2+
RTINSECONDS=899
139.8814 7.0
174.1123 33.9
190.3428 8.3
195.9004 16.6
235.7046 13.4
251.4643 25.5
323.6027 34.2
323.9795 15.2
333.0203 6.6
391.9262 12.9
401.7081 2.6
404.3231 24.5
493.7261 52.8
517.7377 112.4
528.3582 12.6
529.8118 6.4
535.2626 44.9
544.7894 2.6
568.5676 4.2
595.6061 1.5
597.1158 3.6
598.8732 2.6
643.2838 30.0
659.3803 1.4
696.2468 6.0
708.7108 18.6
714.0148 9.4
737.4589 25.0
738.9776 8.2
743.5013 28.2
779.8825 22.1
781.3464 17.6
789.8436 4.8
817.9951 4.8
849.1688 5.1
END IONS
BEGIN IONS
TITLE=corrida_1.1718.1718.2
PEPMASS=862.46835
CHARGE=2+
RTINSECONDS=922
217.9113 32.8
232.2355 9.0
236.5427 4.5
279.5112 22.9
357.7309 11.5
399.2435 8.4
401.4853 19.1
402.6787 28.3
434.0523 9.7
460.3822 13.9
469.7116 2.3
486.2761 1.9
524.3557 28.0
585.4572 9.3
654.6055 5.8
656.9129 4.7
683.6987 4.1
692.5002 12.5
701.6626 83.4
780.6109 11.0
923.9665 17.5
979.4373 34.3
996.6154 1.0
1009.1639 8.7
1044.0842 2.1
1060.4573 19.1
1084.2833 1.4
1181.4644 8.5
1274.2056 77.0
1309.3470 2.1
1359.0301 22.3
1362.2303 12.7
1417.6724 2.6
1433.2871 11.8
1438.2596 30.8
1499.4495 13.8
1558.8511 16.1
1630.8198 25.7
1632.3822 13.4
1635.9112 16.6
1652.6729 3.8
1685.0225 11.5
1700.3340 9.1
1721.8072 6.8
END IONS
BEGIN IONS
TITLE=corrida_1.1755.1755.2
PEPMASS=605.03826
CHARGE=2+
RTINSECONDS=945
143.0314 55.9
217.2116 7.8
220.2091 9.4
292.6420 4.8
304.2939 7.1
318.7265 11.5
387.8209 22.4
393.5404 1.0
395.9434 19.0
410.0852 13.8
449.7269 15.8
455.8915 12.2
456.0918 7.9
619.7551 32.9
663.7563 6.2
685.7668 2.4
696.6644 12.7
710.4945 8.4
768.1963 9.2
779.2731 16.1
787.0694 39.7
864.6619 8.6
889.9593 36.3
900.4630 51.5
916.2580 30.8
919.2739 37.0
924.0944 3.4
996.6702 12.7
1005.9319 20.7
1028.3259 66.1
1067.4403 8.6
1081.6607 2.0
1123.9326 32.8
1204.8484 4.7
END IONS
BEGIN IONS
TITLE=corrida_1.1792.1792.2
PEPMASS=808.58865
CHARGE=2+
RTINSECONDS=968
145.5656 4.3
215.9454 12.7
225.7836 14.1
285.1846 26.5
327.7959 17.9
399.7574 14.4
399.8342 5.2
435.4311 11.9
455.5857 6.2
462.9233 46.1
476.5685 28.2
507.0191 8.1
511.4436 29.1
518.6506 88.6
539.8846 40.9
541.3250 41.4
557.0842 9.8
633.5876 22.7
643.4462 26.0
665.0320 43.0
701.2835 13.1
703.4487 53.5
707.1768 18.6
735.1185 15.1
761.7891 82.0
782.1332 5.0
836.7240 46.9
841.7072 27.3
848.0303 17.7
851.6157 34.1
890.1534 23.4
903.3500 20.5
945.6372 38.4
954.2729 4.6
982.5431 11.1
1055.6578 31.9
1065.1254 16.1
1105.0221 33.8
1184.9626 8.4
1211.3548 7.1
1255.6927 35.8
1289.2074 17.6
1305.5352 4.6
1329.5376 41.6
1559.7145 2.8
END IONS
BEGIN IONS
TITLE=corrida_1.1829.1829.2
PEPMASS=907.18171
CHARGE=2+
RTINSECONDS=991
144.7245 23.3
192.3996 35.9
331.5041 18.9
363.5407 3.5
372.6999 62.5
404.2600 4.4
419.1662 44.2
469.6641 10.4
503.3560 6.7
548.4282 18.5
557.0968 24.9
638.1275 16.1
659.5541 4.4
732.4724 55.3
735.4808 3.5
789.4865 1.7
825.8335 2.4
833.5625 5.0
834.3234 13.9
932.6925 168.2
1023.2614 25.4
1105.3540 11.9
1112.5342 5.6
1185.5620 7.1
1262.9459 5.7
1304.3898 50.3
1307.9246 2.5
1328.6215 3.8
1341.1729 14.1
1343.5639 32.0
1401.3788 14.2
1487.2772 1.1
1493.0907 2.2
1501.9760 13.9
1625.7378 4.1
1627.6823 38.4
1635.5258 10.0
1749.1437 2.8
1764.4079 34.6
END IONS
BEGIN IONS
TITLE=corrida_1.1866.1866.3
PEPMASS=654.98169
CHARGE=3+
RTINSECONDS=1014
134.4146 2.8
147.2238 27.9
175.1147 29.5
246.1544 76.5
283.0948 5.6
333.1229 11.7
376.5637 15.9
403.4676 8.9
409.2211 88.8
439.0301 21.9
522.3038 33.5
536.6267 4.6
619.3562 29.7
635.5247 5.8
666.6535 1.2
708.3010 24.3
734.3844 62.1
758.9921 9.7
779.3435 17.7
798.8962 8.7
847.4619 52.9
851.9946 16.3
918.5047 89.1
1027.1521 8.7
1045.4397 43.7
1055.5658 24.2
1219.3696 21.4
1229.5698 45.5
1255.6488 65.5
1270.7840 51.8
1316.6366 11.1
1321.6304 4.6
1344.5930 55.0
1412.9989 4.3
1441.6406 37.5
1483.7502 55.3
1490.0618 9.7
1604.9712 28.6
1630.8200 67.8
1646.6837 15.5
1649.4411 4.2
1774.9185 5.6
1828.2615 5.7
1843.7553 2.8
1846.9471 13.9
1861.8921 54.7
1903.5455 2.1
END IONS
BEGIN IONS
TITLE=corrida_1.1903.1903.2
PEPMASS=527.91297
CHARGE=2+
RTINSECONDS=1037
149.9198 21.7
172.1438 13.4
220.4502 30.7
228.8696 10.1
230.7056 41.1
281.7239 46.5
288.1192 5.3
330.2985 31.1
369.2048 10.1
384.5144 9.7
408.4659 4.9
410.2889 38.4
411.7872 32.0
456.7747 1.5
464.0100 15.1
518.5713 6.6
545.3205 13.0
572.5752 12.3
579.9069 1.5
597.3404 43.5
654.6176 14.1
668.1890 31.4
746.7103 29.9
784.5832 11.3
810.3097 22.0
827.0753 1.3
845.9382 42.0
846.2753 26.6
882.8970 15.2
910.0996 9.4
931.0774 63.5
931.2595 16.2
944.4538 33.9
949.9850 8.8
979.4708 14.8
989.6383 7.6
999.1778 55.9
1002.9988 5.2
END IONS
BEGIN IONS
TITLE=corrida_1.1940.1940.2
PEPMASS=665.36248
CHARGE=2+
RTINSECONDS=1060
118.6694 3.8
147.9227 7.3
175.1218 54.2
252.1681 3.7
255.3256 7.5
307.4014 14.4
334.0932 15.8
339.1056 25.2
349.8000 13.9
415.1994 41.4
441.4713 2.0
474.2979 31.5
551.7898 15.9
596.1438 7.8
603.3384 62.1
614.2941 31.0
707.9348 3.5
727.3803 11.3
787.4694 65.2
795.2357 1.7
915.5238 20.6
916.3867 2.2
931.8365 10.4
991.3000 28.2
1002.9704 25.2
1007.0639 9.3
1015.6196 16.0
1021.8005 18.0
1028.6074 65.1
1067.4039 12.0
1264.4369 2.7
END IONS
BEGIN IONS
TITLE=corrida_1.1977.1977.2
PEPMASS=772.92755
CHARGE=2+
RTINSECONDS=1083
118.3901 11.8
200.1028 47.6
222.4069 2.6
260.2066 92.3
299.1747 34.7
328.1253 8.6
331.2302 28.3
337.9369 8.2
344.7691 49.1
413.2184 4.6
441.7254 8.9
478.2957 31.0
571.1447 8.0
586.9927 1.2
591.3817 49.0
604.5050 9.2
611.3142 47.8
643.3055 12.8
692.4327 87.1
739.3731 44.0
751.5578 13.2
806.4793 37.5
853.4153 17.7
934.5362 65.9
954.4623 54.2
966.8625 1.5
977.6922 3.4
995.1660 4.9
1001.3315 2.3
1009.9476 3.0
1031.5858 49.1
1055.4147 8.6
1057.7472 5.3
1078.4388 1.6
1145.6287 95.5
1170.9669 8.0
1214.6196 54.5
1246.6725 37.0
1257.6022 15.6
1285.6557 52.7
1298.1872 1.7
1307.3271 3.5
1333.5797 7.1
1345.7566 35.3
1433.0924 12.3
1441.0417 17.3
1459.4823 28.1
1493.8162 1.3
END IONS
BEGIN IONS
TITLE=corrida_1.2014.2014.2
PEPMASS=920.36425
CHARGE=2+
RTINSECONDS=1106
209.0766 12.5
275.8243 6.1
333.5700 22.8
334.2052 39.4
340.0576 5.1
355.0007 12.7
380.3563 25.8
455.7441 17.2
575.7499 11.9
668.4332 14.5
701.6644 31.4
775.2008 4.1
782.0958 44.0
818.5644 14.5
841.9186 22.6
860.4653 24.4
936.5215 3.1
972.0636 7.8
985.1588 43.9
1015.0643 27.5
1023.4549 62.8
1086.2540 3.4
1213.5149 1.1
1233.7334 3.2
1277.4622 11.5
1361.1568 3.0
1406.8630 67.4
1411.2355 6.9
1435.6537 59.5
1465.7643 41.4
1546.4566 29.0
1588.0836 19.2
1592.3859 8.0
1607.7603 1.6
1697.0498 32.0
1741.3810 6.8
1772.3107 16.2
END IONS";

const CISTEINAS: &[&str] = &["carbamidometiladas", "libres"];
const MAXIMO_PROTEINA: usize = 100_000;
const MAXIMO_BYTES_MGF: u64 = 500_000_000;
const MAXIMO_ESPECTROS: usize = 1_000_000;
const MAXIMO_FILAS: usize = 2000;
const MAXIMO_PICOS_ANOTADOS: usize = 80;
const LARGO_EN_GRAFICA: usize = 40;
const SITIOS_DETALLADOS: usize = 12;
const COMBINACIONES_DETALLADAS: usize = 12;
const PEPTIDOS_RESUMIDOS_CADA: usize = 5;
const CANDIDATOS_DETALLADOS: usize = 8;
const PSM_DETALLADOS: usize = 40;
const SENUELOS_DETALLADOS: usize = 10;
const FILAS_EN_PIZARRA: usize = 12;
const LETRAS_EN_PIZARRA: usize = 60;
const ERROR_BIEN_CALIBRADO: f64 = 5.0;
const PSM_PARA_FDR_ESTABLE: usize = 100;
const COBERTURA_BUENA: f64 = 0.5;
const PPM_DE_ALTA_RESOLUCION: f64 = 20.0;

const L_DIGERIR: &str =
    "pub fn digerir(proteina: &[u8], digestion: &Digestion) -> Vec<PeptidoDigerido> {";
const L_ES_CORTE: &str = "(letra == b'K' || letra == b'R') && siguiente != Some(b'P')";
const L_FRONTERA: &str = "fronteras.push(posicion + 1);";
const L_FIN_DE_PROTEINA: &str = "fronteras.push(proteina.len());";
const L_LARGO_MAXIMO: &str = "if largo > digestion.largo_maximo {";
const L_LARGO_MINIMO: &str = "if largo < digestion.largo_minimo {";
const L_AGUA: &str = "let mut total = MASA_DEL_AGUA;";
const L_MASA_DE_RESIDUO: &str = "total += masa_de_residuo(letra, carbamidometilo)?;";
const L_GUARDAR: &str = "peptidos.push(PeptidoDigerido {";

const L_BUSCAR: &str = "pub fn buscar_espectros(";
const L_REVISAR: &str = "revisar_busqueda(busqueda)?;";
const L_SENUELOS: &str = "let proteinas = con_senuelos(objetivos);";
const L_DIGERIR_PROTEINA: &str = "digerir(&proteina.secuencia, digestion)";
const L_ORDENAR_INDICE: &str = "peptidos.par_sort_unstable_by(|a, b| {";
const L_CADA_ESPECTRO: &str =
    ".filter_map(|(numero, espectro)| mejor_psm(numero, espectro, indice, busqueda))";
const L_NORMALIZAR: &str = "let picos = normalizar(&espectro.picos);";
const L_MASA_NEUTRA: &str = "let masa = masa_neutra(espectro.mz_precursor, carga);";
const L_VENTANA: &str = "let (minima, maxima) = ventana(masa, busqueda.tolerancia_ppm);";
const L_EN_VENTANA: &str = "for peptido in indice.en_ventana(minima, maxima) {";
const L_IONES: &str = "llenar_iones(indice.secuencia(peptido), carbamidometilo, &mut b, &mut y);";
const L_PUNTUAR: &str = "let puntaje = puntuar(&b, &y, &picos, busqueda.tolerancia_fragmento);";
const L_RECORD: &str = "if puntaje.hiperpuntaje > record {";
const L_PICO: &str = "if let Some(indice) = pico_cercano(picos, ion, tolerancia) {";
const L_SUMA: &str = "intensidad += picos[indice].1;";
const L_HIPERPUNTAJE: &str =
    "(1.0 + intensidad).ln() + logaritmo_factorial(iones_b) + logaritmo_factorial(iones_y)";
const L_ORDENAR_PSM: &str = "ordenar_por_puntaje(psms);";
const L_CONTAR_SENUELO: &str = "if psms[fin].senuelo {";
const L_TASA: &str = "let tasa = tasa_de_falsos(senuelos, objetivos);";
const L_MINIMO: &str = "minimo = minimo.min(*tasa);";
const L_Q: &str = "psm.q = minimo;";

pub fn algoritmos() -> Vec<Algoritmo> {
    vec![
        Algoritmo {
            id: "antibioticos.digestion",
            capitulo: "antibioticos",
            nombre: "Digestión con tripsina (vanguardia)",
            resumen: "Corta una proteína después de K o R (salvo antes de P), admite cortes perdidos y da la masa monoisotópica exacta de cada péptido; a diferencia de las masas enteras del capítulo, es lo que usa un buscador real.",
            campos: &[
                Campo {
                    id: "proteina",
                    etiqueta: "Proteína",
                    tipo: Tipo::Proteina,
                    ejemplo: "@proteina:recA",
                },
                Campo {
                    id: "cortes_perdidos",
                    etiqueta: "Cortes perdidos permitidos (0 a 3)",
                    tipo: Tipo::Entero,
                    ejemplo: "1",
                },
                Campo {
                    id: "cisteinas",
                    etiqueta: "Cisteínas",
                    tipo: Tipo::Opcion(CISTEINAS),
                    ejemplo: "carbamidometiladas",
                },
            ],
            funcion: "proteomica::digerir",
            ejecutar: digestion,
            pasos: Some(pasos_digestion),
        },
        Algoritmo {
            id: "antibioticos.busqueda_de_espectros",
            capitulo: "antibioticos",
            nombre: "Identificar péptidos en una base de datos (vanguardia)",
            resumen: "Compara cada espectro MS/MS con todos los péptidos trípticos de E. coli de la misma masa, elige el mejor por hiperpuntaje y controla los falsos con señuelos al revés y q-valores; a diferencia de secuenciar de novo, no adivina: busca en el proteoma y dice cuánto azar hay.",
            campos: &[
                Campo {
                    id: "espectros",
                    etiqueta: "Espectros en formato MGF",
                    tipo: Tipo::Texto,
                    ejemplo: ESPECTROS_DE_EJEMPLO,
                },
                Campo {
                    id: "archivo",
                    etiqueta: "O un archivo .mgf (opcional)",
                    tipo: Tipo::Archivo,
                    ejemplo: "",
                },
                Campo {
                    id: "tolerancia_ppm",
                    etiqueta: "Tolerancia del precursor (ppm)",
                    tipo: Tipo::Decimal,
                    ejemplo: "10",
                },
                Campo {
                    id: "tolerancia_fragmento",
                    etiqueta: "Tolerancia de los fragmentos (Da)",
                    tipo: Tipo::Decimal,
                    ejemplo: "0.02",
                },
                Campo {
                    id: "cortes_perdidos",
                    etiqueta: "Cortes perdidos permitidos (0 a 3)",
                    tipo: Tipo::Entero,
                    ejemplo: "1",
                },
                Campo {
                    id: "cisteinas",
                    etiqueta: "Cisteínas",
                    tipo: Tipo::Opcion(CISTEINAS),
                    ejemplo: "carbamidometiladas",
                },
            ],
            funcion: "proteomica::buscar_espectros",
            ejecutar: busqueda_de_espectros,
            pasos: Some(pasos_busqueda),
        },
    ]
}

fn objetivos_de_ecoli() -> Vec<(String, Vec<u8>)> {
    ecoli::genes()
        .iter()
        .filter(|gen_de_ecoli| !gen_de_ecoli.proteina.is_empty())
        .map(|gen_de_ecoli| (gen_de_ecoli.nombre.clone(), gen_de_ecoli.proteina.clone()))
        .collect()
}

pub fn preparar() -> &'static Indice {
    static INDICE: OnceLock<Indice> = OnceLock::new();
    INDICE.get_or_init(|| {
        Indice::construir(&objetivos_de_ecoli(), Digestion::default())
            .expect("el proteoma incluido de E. coli se puede digerir")
    })
}

fn indice_para(digestion: Digestion) -> Result<Cow<'static, Indice>, String> {
    if digestion == Digestion::default() {
        return Ok(Cow::Borrowed(preparar()));
    }
    Indice::construir(&objetivos_de_ecoli(), digestion).map(Cow::Owned)
}

fn como_texto(letras: &[u8]) -> String {
    String::from_utf8_lossy(letras).to_string()
}

fn recortado(letras: &[u8], maximo: usize) -> String {
    if letras.len() <= maximo {
        return como_texto(letras);
    }
    format!("{}…", como_texto(&letras[..maximo]))
}

fn plural(numero: usize, singular: &str, varios: &str) -> String {
    if numero == 1 {
        format!("1 {singular}")
    } else {
        format!("{} {varios}", miles(numero))
    }
}

fn verbo(numero: usize, singular: &'static str, varios: &'static str) -> &'static str {
    if numero == 1 { singular } else { varios }
}

fn con_vecinos(proteina: &[u8], inicio: usize, fin: usize) -> String {
    let antes = if inicio == 0 {
        '-'
    } else {
        proteina[inicio - 1] as char
    };
    let despues = proteina.get(fin).map_or('-', |letra| *letra as char);
    format!("{antes}.{}.{despues}", como_texto(&proteina[inicio..fin]))
}

fn leer_cisteinas(entrada: &Entrada) -> Result<bool, String> {
    Ok(entrada.opcion("cisteinas", CISTEINAS)? == CISTEINAS[0])
}

fn leer_digestion(entrada: &Entrada) -> Result<Digestion, String> {
    let digestion = Digestion {
        cortes_perdidos: entrada.natural("cortes_perdidos")?,
        carbamidometilo: leer_cisteinas(entrada)?,
        ..Digestion::default()
    };
    proteomica::revisar_digestion(&digestion)?;
    Ok(digestion)
}

fn leer_proteina(entrada: &Entrada) -> Result<Vec<u8>, String> {
    let proteina: Vec<u8> = entrada
        .proteina("proteina")?
        .iter()
        .map(|letra| letra.to_ascii_uppercase())
        .filter(|letra| letra.is_ascii_uppercase())
        .collect();
    if proteina.is_empty() {
        return Err("La proteína no tiene letras.".to_string());
    }
    if proteina.len() > MAXIMO_PROTEINA {
        return Err(format!(
            "La proteína tiene {} letras; el máximo es {}.",
            miles(proteina.len()),
            miles(MAXIMO_PROTEINA)
        ));
    }
    Ok(proteina)
}

fn sitios_bloqueados(proteina: &[u8]) -> usize {
    proteina
        .windows(2)
        .filter(|par| matches!(par[0], b'K' | b'R') && par[1] == b'P')
        .count()
}

fn largos_de_pedazos(fronteras: &[usize]) -> Vec<usize> {
    fronteras.windows(2).map(|par| par[1] - par[0]).collect()
}

fn cobertura(proteina: &[u8], peptidos: &[PeptidoDigerido]) -> f64 {
    let mut cubiertas = vec![false; proteina.len()];
    for peptido in peptidos.iter().filter(|p| p.cortes_perdidos == 0) {
        for cubierta in &mut cubiertas[peptido.inicio..peptido.fin] {
            *cubierta = true;
        }
    }
    cubiertas.iter().filter(|&&cubierta| cubierta).count() as f64 / proteina.len().max(1) as f64
}

fn tiene_i_o_l(letras: &[u8]) -> bool {
    letras.iter().any(|letra| matches!(letra, b'I' | b'L'))
}

fn grafica_de_largos(largos: &[usize], digestion: &Digestion) -> Bloque {
    let mut cuentas = vec![0usize; LARGO_EN_GRAFICA + 1];
    for &largo in largos {
        cuentas[largo.min(LARGO_EN_GRAFICA + 1) - 1] += 1;
    }
    let mut categorias: Vec<String> = (1..=LARGO_EN_GRAFICA)
        .map(|largo| largo.to_string())
        .collect();
    categorias.push(format!("más de {LARGO_EN_GRAFICA}"));
    Bloque::barras(
        &format!(
            "Largo de los pedazos sin cortes perdidos (la búsqueda usa los de {} a {})",
            digestion.largo_minimo, digestion.largo_maximo
        ),
        "aminoácidos",
        "pedazos",
        categorias,
        cuentas.into_iter().map(|cuenta| cuenta as f64).collect(),
    )
}

fn digestion(entrada: &Entrada) -> Result<Salida, String> {
    let proteina = leer_proteina(entrada)?;
    let reglas = leer_digestion(entrada)?;
    let peptidos = proteomica::digerir(&proteina, &reglas);
    let fronteras = proteomica::fronteras(&proteina);
    let largos = largos_de_pedazos(&fronteras);
    let cortes = fronteras.len().saturating_sub(2);
    let bloqueados = sitios_bloqueados(&proteina);
    let cubierta = cobertura(&proteina, &peptidos);
    let completos = peptidos.iter().filter(|p| p.cortes_perdidos == 0).count();
    let cortos = largos
        .iter()
        .filter(|&&largo| largo < reglas.largo_minimo)
        .count();
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor(
        "Largo de la proteína",
        format!("{} aminoácidos", miles(proteina.len())),
    ));
    salida.agregar(Bloque::valor(
        "Sitios de corte (K o R sin P detrás)",
        miles(cortes),
    ));
    salida.agregar(Bloque::valor(
        "Sitios bloqueados por una prolina",
        miles(bloqueados),
    ));
    salida.agregar(Bloque::valor(
        &format!(
            "Péptidos de {} a {} aminoácidos",
            reglas.largo_minimo, reglas.largo_maximo
        ),
        format!(
            "{} ({} sin cortes perdidos)",
            miles(peptidos.len()),
            miles(completos)
        ),
    ));
    salida.agregar(Bloque::valor(
        "Secuencia cubierta por péptidos sin cortes perdidos",
        porcentaje(cubierta, 1),
    ));
    let filas: Vec<Vec<String>> = peptidos
        .iter()
        .take(MAXIMO_FILAS)
        .enumerate()
        .map(|(numero, peptido)| {
            vec![
                (numero + 1).to_string(),
                con_vecinos(&proteina, peptido.inicio, peptido.fin),
                format!("{}–{}", peptido.inicio + 1, peptido.fin),
                (peptido.fin - peptido.inicio).to_string(),
                peptido.cortes_perdidos.to_string(),
                decimal(peptido.masa, 5),
                decimal(proteomica::mz_con_carga(peptido.masa, 2), 5),
                decimal(proteomica::mz_con_carga(peptido.masa, 3), 5),
            ]
        })
        .collect();
    salida.agregar(Bloque::tabla_con_total(
        "Péptidos trípticos (el punto separa el residuo anterior y el siguiente)",
        &[
            "#",
            "Péptido",
            "Posición",
            "Largo",
            "Cortes perdidos",
            "Masa monoisotópica (Da)",
            "m/z 2+",
            "m/z 3+",
        ],
        filas,
        peptidos.len(),
    ));
    salida.agregar(grafica_de_largos(&largos, &reglas));
    explicar_digestion(
        &mut salida,
        &proteina,
        &peptidos,
        &reglas,
        cubierta,
        cortos,
        bloqueados,
    );
    Ok(salida)
}

fn explicar_digestion(
    salida: &mut Salida,
    proteina: &[u8],
    peptidos: &[PeptidoDigerido],
    reglas: &Digestion,
    cubierta: f64,
    cortos: usize,
    bloqueados: usize,
) {
    salida.explicar(format!(
        "La tripsina deja {} con hasta {}; esos son los únicos péptidos de esta proteína que una búsqueda en base de datos puede comparar con un espectro.",
        plural(peptidos.len(), "péptido de 6 a 40 aminoácidos", "péptidos de 6 a 40 aminoácidos"),
        plural(reglas.cortes_perdidos, "corte perdido", "cortes perdidos")
    ));
    if cubierta >= COBERTURA_BUENA {
        salida.explicar(format!(
            "Los péptidos sin cortes perdidos cubren {} de la secuencia: buena cobertura para identificar la proteína con varios péptidos.",
            porcentaje(cubierta, 1)
        ));
    } else {
        salida.explicar(format!(
            "Los péptidos sin cortes perdidos cubren solo {} de la secuencia: el resto cae en pedazos muy cortos o muy largos; una segunda enzima (Glu-C, quimotripsina) cubriría lo que falta.",
            porcentaje(cubierta, 1)
        ));
    }
    if cortos > 0 {
        salida.explicar(format!(
            "{} de menos de {} letras {}: su masa se repite en muchas proteínas y {} identificar nada por sí {}.",
            plural(cortos, "pedazo", "pedazos"),
            reglas.largo_minimo,
            verbo(cortos, "se descarta", "se descartan"),
            verbo(cortos, "casi nunca sirve para", "casi nunca sirven para"),
            verbo(cortos, "solo", "solos")
        ));
    }
    if bloqueados > 0 {
        salida.explicar(format!(
            "{} con K o R {} porque una prolina venía detrás (regla de Keil); los cortes perdidos reales salen de sitios lentos como ese o de K/R seguidas de D o E.",
            plural(bloqueados, "sitio", "sitios"),
            verbo(bloqueados, "no se cortó", "no se cortaron")
        ));
    }
    if tiene_i_o_l(proteina) {
        salida.explicar(
            "I y L pesan exactamente lo mismo (113.08406 Da): ningún espectro de masas las distingue, así que un péptido identificado con una de ellas vale igual con la otra.",
        );
    }
    salida.explicar(if reglas.carbamidometilo {
        "Las masas son monoisotópicas (con los isótopos más comunes, a 5 decimales), no las masas enteras de «Espectro de un péptido lineal»; las C llevan +57.02146 Da de carbamidometilo, como en casi todo experimento que reduce y alquila."
    } else {
        "Las masas son monoisotópicas (con los isótopos más comunes, a 5 decimales), no las masas enteras de «Espectro de un péptido lineal»; las C van sin modificar."
    });
}

struct PedidoDeBusqueda {
    espectros: Vec<Espectro>,
    busqueda: Busqueda,
    digestion: Digestion,
}

fn texto_de_espectros(entrada: &Entrada) -> Result<String, String> {
    if entrada.tiene("archivo") {
        let ruta = entrada.ruta("archivo")?;
        let bytes = std::fs::metadata(&ruta)
            .map(|datos| datos.len())
            .unwrap_or(0);
        if bytes > MAXIMO_BYTES_MGF {
            return Err(format!(
                "El archivo pesa más de {} MB.",
                MAXIMO_BYTES_MGF / 1_000_000
            ));
        }
        return std::fs::read_to_string(&ruta)
            .map_err(|error| format!("No se pudo leer {}: {error}", ruta.display()));
    }
    if entrada.tiene("espectros") {
        return entrada.texto("espectros");
    }
    Err("Pega espectros en formato MGF o elige un archivo .mgf.".to_string())
}

fn leer_busqueda(entrada: &Entrada) -> Result<PedidoDeBusqueda, String> {
    let busqueda = Busqueda {
        tolerancia_ppm: entrada.decimal("tolerancia_ppm")?,
        tolerancia_fragmento: entrada.decimal("tolerancia_fragmento")?,
    };
    proteomica::revisar_busqueda(&busqueda)?;
    let digestion = leer_digestion(entrada)?;
    let espectros = proteomica::leer_mgf(&texto_de_espectros(entrada)?)?;
    if espectros.len() > MAXIMO_ESPECTROS {
        return Err(format!(
            "El archivo tiene {} espectros; el máximo es {}.",
            miles(espectros.len()),
            miles(MAXIMO_ESPECTROS)
        ));
    }
    Ok(PedidoDeBusqueda {
        espectros,
        busqueda,
        digestion,
    })
}

fn contar_candidatos(espectros: &[Espectro], indice: &Indice, busqueda: &Busqueda) -> usize {
    espectros
        .par_iter()
        .map(|espectro| {
            espectro
                .cargas
                .iter()
                .map(|&carga| {
                    let masa = proteomica::masa_neutra(espectro.mz_precursor, carga);
                    let (minima, maxima) = proteomica::ventana(masa, busqueda.tolerancia_ppm);
                    indice.en_ventana(minima, maxima).len()
                })
                .sum::<usize>()
        })
        .sum()
}

fn nombre_de_proteina(indice: &Indice, peptido: &PeptidoIndexado) -> String {
    indice.proteinas[peptido.proteina as usize].nombre.clone()
}

fn peptido_con_vecinos(indice: &Indice, peptido: &PeptidoIndexado) -> String {
    con_vecinos(
        &indice.proteinas[peptido.proteina as usize].secuencia,
        peptido.inicio as usize,
        peptido.fin as usize,
    )
}

fn tipo_de_psm(psm: &Psm) -> &'static str {
    if psm.senuelo { "señuelo" } else { "objetivo" }
}

fn pasa(psm: &Psm) -> bool {
    !psm.senuelo && psm.q <= FDR_ESTANDAR
}

fn iones_de(indice: &Indice, peptido: &PeptidoIndexado) -> (Vec<f64>, Vec<f64>) {
    let mut b = Vec::new();
    let mut y = Vec::new();
    proteomica::llenar_iones(
        indice.secuencia(peptido),
        indice.digestion.carbamidometilo,
        &mut b,
        &mut y,
    );
    (b, y)
}

fn etiquetas_de_picos(
    picos: &[(f64, f64)],
    b: &[f64],
    y: &[f64],
    tolerancia: f64,
) -> Vec<Vec<String>> {
    let mut etiquetas = vec![Vec::new(); picos.len()];
    for (nombre, iones) in [("b", b), ("y", y)] {
        for (numero, &ion) in iones.iter().enumerate() {
            if let Some(cual) = proteomica::pico_cercano(picos, ion, tolerancia) {
                etiquetas[cual].push(format!("{nombre}{}", numero + 1));
            }
        }
    }
    etiquetas
}

fn espectro_anotado(
    espectro: &Espectro,
    psm: &Psm,
    indice: &Indice,
    busqueda: &Busqueda,
) -> Bloque {
    let picos = proteomica::normalizar(&espectro.picos);
    let (b, y) = iones_de(indice, &psm.peptido);
    let etiquetas = etiquetas_de_picos(&picos, &b, &y, busqueda.tolerancia_fragmento);
    let mut elegidos: Vec<usize> = (0..picos.len()).collect();
    if picos.len() > MAXIMO_PICOS_ANOTADOS {
        elegidos.sort_by(|&a, &c| {
            etiquetas[a]
                .is_empty()
                .cmp(&etiquetas[c].is_empty())
                .then(picos[c].1.total_cmp(&picos[a].1))
        });
        elegidos.truncate(MAXIMO_PICOS_ANOTADOS);
        elegidos.sort_unstable();
    }
    let categorias = elegidos
        .iter()
        .map(|&cual| {
            if etiquetas[cual].is_empty() {
                decimal(picos[cual].0, 2)
            } else {
                format!(
                    "{} {}",
                    decimal(picos[cual].0, 2),
                    etiquetas[cual].join("/")
                )
            }
        })
        .collect();
    let valores = elegidos.iter().map(|&cual| picos[cual].1).collect();
    Bloque::barras(
        &format!(
            "Espectro anotado: «{}» = {} (picos con nombre = iones b/y emparejados)",
            espectro.titulo,
            como_texto(indice.secuencia(&psm.peptido))
        ),
        "m/z",
        "intensidad relativa (máximo = 100)",
        categorias,
        valores,
    )
}

fn escalera_de_iones(
    espectro: &Espectro,
    psm: &Psm,
    indice: &Indice,
    busqueda: &Busqueda,
) -> Bloque {
    let picos = proteomica::normalizar(&espectro.picos);
    let peptido = indice.secuencia(&psm.peptido);
    let (b, y) = iones_de(indice, &psm.peptido);
    let n = peptido.len();
    let marca = |ion: f64| {
        if proteomica::pico_cercano(&picos, ion, busqueda.tolerancia_fragmento).is_some() {
            "sí"
        } else {
            "no"
        }
    };
    let filas = (1..=n)
        .map(|posicion| {
            let (texto_b, hallado_b) = if posicion < n {
                (decimal(b[posicion - 1], 4), marca(b[posicion - 1]))
            } else {
                ("—".to_string(), "")
            };
            let numero_y = n + 1 - posicion;
            let (texto_y, hallado_y) = if posicion > 1 {
                (decimal(y[numero_y - 1], 4), marca(y[numero_y - 1]))
            } else {
                ("—".to_string(), "")
            };
            vec![
                format!("b{posicion}"),
                (peptido[posicion - 1] as char).to_string(),
                texto_b,
                hallado_b.to_string(),
                texto_y,
                hallado_y.to_string(),
                format!("y{numero_y}"),
            ]
        })
        .collect();
    Bloque::tabla(
        "Escalera de iones del mejor PSM (iones de carga 1)",
        &["b", "Residuo", "m/z b", "¿Visto?", "m/z y", "¿Visto?", "y"],
        filas,
    )
}

fn busqueda_de_espectros(entrada: &Entrada) -> Result<Salida, String> {
    let pedido = leer_busqueda(entrada)?;
    let indice = indice_para(pedido.digestion)?;
    let reloj = Instant::now();
    let psms = proteomica::buscar_espectros(&pedido.espectros, &indice, &pedido.busqueda)?;
    let segundos = reloj.elapsed().as_secs_f64();
    let candidatos = contar_candidatos(&pedido.espectros, &indice, &pedido.busqueda);
    let aceptados = proteomica::aceptados(&psms, FDR_ESTANDAR);
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor(
        "Espectros leídos",
        miles(pedido.espectros.len()),
    ));
    salida.agregar(Bloque::valor(
        "PSM que pasan q ≤ 0.01 (1 % de FDR)",
        miles(aceptados.len()),
    ));
    salida.agregar(Bloque::valor(
        "Péptidos en el índice (objetivos + señuelos)",
        miles(indice.peptidos.len()),
    ));
    salida.agregar(Bloque::valor(
        "Candidatos puntuados (dentro de la ventana del precursor)",
        miles(candidatos),
    ));
    salida.agregar(Bloque::valor(
        "Tiempo de la búsqueda",
        format!("{} ms", decimal(segundos * 1000.0, 1)),
    ));
    let filas: Vec<Vec<String>> = psms
        .iter()
        .take(MAXIMO_FILAS)
        .map(|psm| {
            vec![
                pedido.espectros[psm.espectro].titulo.clone(),
                peptido_con_vecinos(&indice, &psm.peptido),
                nombre_de_proteina(&indice, &psm.peptido),
                format!("{}+", psm.carga),
                decimal(psm.error_ppm, 2),
                psm.puntaje.iones_b.to_string(),
                psm.puntaje.iones_y.to_string(),
                decimal(psm.puntaje.hiperpuntaje, 2),
                decimal(psm.q, 3),
                tipo_de_psm(psm).to_string(),
                if pasa(psm) { "sí" } else { "no" }.to_string(),
            ]
        })
        .collect();
    salida.agregar(Bloque::tabla_con_total(
        "Mejor PSM de cada espectro, de mayor a menor hiperpuntaje",
        &[
            "Espectro",
            "Péptido",
            "Proteína",
            "Carga",
            "Error (ppm)",
            "Iones b",
            "Iones y",
            "Hiperpuntaje",
            "q",
            "Tipo",
            "¿Pasa 1 %?",
        ],
        filas,
        psms.len(),
    ));
    let destacado = aceptados.first().or(psms.first());
    if let Some(psm) = destacado {
        let espectro = &pedido.espectros[psm.espectro];
        salida.agregar(espectro_anotado(espectro, psm, &indice, &pedido.busqueda));
        salida.agregar(escalera_de_iones(espectro, psm, &indice, &pedido.busqueda));
    }
    if !psms.is_empty() {
        let puntos = psms
            .iter()
            .enumerate()
            .map(|(lugar, psm)| Punto {
                x: (lugar + 1) as f64,
                y: psm.puntaje.hiperpuntaje,
                grupo: usize::from(psm.senuelo),
                nombre: format!(
                    "{} ({})",
                    como_texto(indice.secuencia(&psm.peptido)),
                    tipo_de_psm(psm)
                ),
            })
            .collect();
        salida.agregar(Bloque::dispersion(
            "Hiperpuntaje por lugar: objetivos (grupo 0) y señuelos (grupo 1)",
            "lugar",
            "hiperpuntaje",
            puntos,
        ));
    }
    explicar_busqueda(&mut salida, &pedido, &indice, &psms, &aceptados, candidatos);
    Ok(salida)
}

fn explicar_busqueda(
    salida: &mut Salida,
    pedido: &PedidoDeBusqueda,
    indice: &Indice,
    psms: &[Psm],
    aceptados: &[Psm],
    candidatos: usize,
) {
    let senuelos = psms.iter().filter(|psm| psm.senuelo).count();
    let objetivos_rechazados = psms.iter().filter(|psm| !psm.senuelo && !pasa(psm)).count();
    let aviso = if psms.len() < PSM_PARA_FDR_ESTABLE {
        format!(
            " Con menos de {PSM_PARA_FDR_ESTABLE} PSM ese cálculo es grueso: un solo señuelo mueve el FDR varios puntos; con las decenas de miles de espectros de una corrida real el 1 % es fino."
        )
    } else {
        String::new()
    };
    salida.explicar(format!(
        "De {}, {} {} el corte estándar q ≤ 0.01: entre los aceptados se espera como mucho 1 % de identificaciones falsas.{aviso}",
        plural(pedido.espectros.len(), "espectro", "espectros"),
        miles(aceptados.len()),
        verbo(aceptados.len(), "pasa", "pasan")
    ));
    if let Some(mejor) = aceptados.first() {
        salida.explicar(format!(
            "El mejor es {} de {}: explica {} b y {} y con {} ppm de error en el precursor.",
            como_texto(indice.secuencia(&mejor.peptido)),
            nombre_de_proteina(indice, &mejor.peptido),
            plural(mejor.puntaje.iones_b, "ion", "iones"),
            plural(mejor.puntaje.iones_y, "ion", "iones"),
            decimal(mejor.error_ppm, 2)
        ));
    }
    salida.explicar(format!(
        "{} {} como mejor pareja un señuelo (una proteína de E. coli leída al revés, que no puede ser la respuesta); cada señuelo estima un objetivo falso de puntaje parecido, y el q-valor es la menor razón señuelos/objetivos de un corte que incluye a ese PSM. {} {} sin pasar: {}.",
        plural(senuelos, "espectro", "espectros"),
        verbo(senuelos, "tuvo", "tuvieron"),
        plural(objetivos_rechazados, "objetivo", "objetivos"),
        verbo(objetivos_rechazados, "quedó", "quedaron"),
        verbo(
            objetivos_rechazados,
            "puede ser real con poca señal, pero no se distingue del azar",
            "pueden ser reales con poca señal, pero no se distinguen del azar"
        )
    ));
    let con_i_o_l = aceptados
        .iter()
        .filter(|psm| tiene_i_o_l(indice.secuencia(&psm.peptido)))
        .count();
    let con_cortes = aceptados
        .iter()
        .filter(|psm| psm.peptido.cortes_perdidos > 0)
        .count();
    if con_i_o_l > 0 || con_cortes > 0 {
        salida.explicar(format!(
            "{} de los aceptados {} I o L: la masa no las distingue (113.08406 Da las dos), así que esa letra es una convención. {} de los aceptados {} un corte perdido, algo normal en una digestión real.",
            miles(con_i_o_l),
            verbo(con_i_o_l, "tiene", "tienen"),
            miles(con_cortes),
            verbo(con_cortes, "tiene", "tienen")
        ));
    }
    if !aceptados.is_empty() {
        let error_medio =
            aceptados.iter().map(|psm| psm.error_ppm.abs()).sum::<f64>() / aceptados.len() as f64;
        if error_medio <= ERROR_BIEN_CALIBRADO
            && pedido.busqueda.tolerancia_ppm > PPM_DE_ALTA_RESOLUCION
        {
            salida.explicar(format!(
                "El error medio del precursor en los aceptados es {} ppm: el espectrómetro está bien calibrado y la ventana de {} ppm es más ancha de lo necesario; con {} ppm o menos habría menos candidatos y menos coincidencias por azar.",
                decimal(error_medio, 2),
                decimal(pedido.busqueda.tolerancia_ppm, 1),
                decimal(PPM_DE_ALTA_RESOLUCION, 0)
            ));
        } else if error_medio <= ERROR_BIEN_CALIBRADO {
            salida.explicar(format!(
                "El error medio del precursor en los aceptados es {} ppm: el espectrómetro está bien calibrado y la ventana de {} ppm es adecuada.",
                decimal(error_medio, 2),
                decimal(pedido.busqueda.tolerancia_ppm, 1)
            ));
        } else {
            salida.explicar(format!(
                "El error medio del precursor en los aceptados es {} ppm, más de {ERROR_BIEN_CALIBRADO}: revisa la calibración antes de estrechar la ventana.",
                decimal(error_medio, 2)
            ));
        }
    }
    let sin_indice = indice.peptidos.len() * pedido.espectros.len();
    salida.explicar(format!(
        "A diferencia de «Secuenciar un péptido cíclico (ramificación y poda)», que arma el péptido de novo con masas enteras, aquí solo se comparan péptidos que ya están en el proteoma: el índice por masa dejó {} candidatos en vez de {} comparaciones, pero un péptido que no está en la base (un contaminante, una variante) nunca puede identificarse bien.",
        miles(candidatos),
        miles(sin_indice)
    ));
}

fn pasos_digestion(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let proteina = leer_proteina(entrada)?;
    let reglas = leer_digestion(entrada)?;
    Ok(narrar_digestion(&proteina, &reglas).0)
}

fn pizarra_de_sitio(proteina: &[u8], posicion: usize, fronteras: &[usize]) -> Pizarra {
    let siguiente = proteina
        .get(posicion + 1)
        .map_or("(fin)".to_string(), |letra| (*letra as char).to_string());
    Pizarra::nueva()
        .variable("posicion", posicion + 1)
        .variable("letra", proteina[posicion] as char)
        .variable("siguiente", siguiente)
        .variable("fronteras", fronteras.len())
        .secuencia(proteina, &[(posicion, posicion + 2, Marca::Actual)])
}

fn narrar_sitios(proteina: &[u8], narracion: &mut Narracion) -> Vec<usize> {
    let mut fronteras = vec![0];
    let mut vistos = 0;
    for posicion in 0..proteina.len() {
        if !matches!(proteina[posicion], b'K' | b'R') {
            continue;
        }
        vistos += 1;
        let corta = proteomica::corta_despues(proteina, posicion);
        let ultima = posicion + 1 == proteina.len();
        if corta && !ultima {
            fronteras.push(posicion + 1);
        }
        if vistos > SITIOS_DETALLADOS || !narracion.hay_lugar() {
            continue;
        }
        let pizarra = pizarra_de_sitio(proteina, posicion, &fronteras);
        let letra = proteina[posicion] as char;
        let paso = if ultima {
            pizarra.paso(
                L_ES_CORTE,
                format!("La {letra} de la posición {} es la última letra: ahí la proteína ya termina, no hay nada que cortar.", posicion + 1),
            )
        } else if corta {
            pizarra.paso(
                L_FRONTERA,
                format!(
                    "{letra} en la posición {} y detrás {}: la tripsina corta aquí. Guardas la frontera {}.",
                    posicion + 1,
                    proteina[posicion + 1] as char,
                    posicion + 1
                ),
            )
        } else {
            pizarra.paso(
                L_ES_CORTE,
                format!(
                    "{letra} en la posición {}, pero detrás viene una P: la prolina bloquea a la tripsina y no hay corte.",
                    posicion + 1
                ),
            )
        };
        narracion.agregar(paso);
    }
    if vistos > SITIOS_DETALLADOS {
        narracion.agregar(
            Pizarra::nueva()
                .variable("fronteras", fronteras.len())
                .paso(
                    L_FRONTERA,
                    format!(
                        "Lo mismo con las otras {} K y R. Ya tienes {} fronteras.",
                        miles(vistos - SITIOS_DETALLADOS),
                        miles(fronteras.len())
                    ),
                ),
        );
    }
    fronteras.push(proteina.len());
    fronteras
}

fn lista_de_pedazos(proteina: &[u8], fronteras: &[usize]) -> Vec<String> {
    fronteras
        .windows(2)
        .take(FILAS_EN_PIZARRA * 2)
        .map(|par| {
            format!(
                "{}–{}: {}",
                par[0] + 1,
                par[1],
                recortado(&proteina[par[0]..par[1]], 30)
            )
        })
        .collect()
}

fn paso_de_masa(proteina: &[u8], inicio: usize, fin: usize, reglas: &Digestion) -> Option<Paso> {
    let mut total = proteomica::MASA_DEL_AGUA;
    let mut filas = Vec::new();
    for &letra in &proteina[inicio..fin] {
        let masa = proteomica::masa_de_residuo(letra, reglas.carbamidometilo)?;
        total += masa;
        filas.push(vec![
            (letra as char).to_string(),
            decimal(masa, 5),
            decimal(total, 5),
        ]);
    }
    Some(
        Pizarra::nueva()
            .variable("peptido", recortado(&proteina[inicio..fin], LETRAS_EN_PIZARRA))
            .variable("total", decimal(total, 5))
            .tabla("Masa residuo por residuo", &["Letra", "Masa (Da)", "total"], filas, &[])
            .paso(
                L_MASA_DE_RESIDUO,
                format!(
                    "Sumas la masa monoisotópica de cada letra a la del agua (18.01056 Da, el H del inicio y el OH del final). Total: {} Da.",
                    decimal(total, 5)
                ),
            ),
    )
}

fn pizarra_de_combinacion(
    proteina: &[u8],
    primera: usize,
    cortes_perdidos: usize,
    inicio: usize,
    fin: usize,
    cuantos: usize,
) -> Pizarra {
    Pizarra::nueva()
        .variable("primera", primera)
        .variable("cortes_perdidos", cortes_perdidos)
        .variable("inicio", inicio + 1)
        .variable("fin", fin)
        .variable("largo", fin - inicio)
        .variable("peptidos", cuantos)
        .secuencia(proteina, &[(inicio, fin, Marca::Actual)])
}

fn narrar_digestion(proteina: &[u8], reglas: &Digestion) -> (Vec<Paso>, Vec<PeptidoDigerido>) {
    let mut narracion = Narracion::nueva(400);
    narracion.agregar(
        Pizarra::nueva()
            .variable("proteina", recortado(proteina, LETRAS_EN_PIZARRA))
            .variable("cortes_perdidos", reglas.cortes_perdidos)
            .variable("largo_minimo", reglas.largo_minimo)
            .variable("largo_maximo", reglas.largo_maximo)
            .secuencia(proteina, &[])
            .paso(
                L_DIGERIR,
                format!(
                    "Vas a digerir en la computadora una proteína de {} aminoácidos con tripsina. Meta: todos los péptidos de {} a {} letras, con hasta {}, y la masa de cada uno.",
                    miles(proteina.len()),
                    reglas.largo_minimo,
                    reglas.largo_maximo,
                    plural(reglas.cortes_perdidos, "corte perdido", "cortes perdidos")
                ),
            ),
    );
    let fronteras = narrar_sitios(proteina, &mut narracion);
    narracion.agregar(
        Pizarra::nueva()
            .variable("fronteras", fronteras.len())
            .lista("Pedazos entre fronteras", lista_de_pedazos(proteina, &fronteras), &[])
            .paso(
                L_FIN_DE_PROTEINA,
                format!(
                    "Cierras con el final de la proteína: {} fronteras dejan {} pedazos sin cortes perdidos. Ahora juntas 1, 2… pedazos seguidos.",
                    miles(fronteras.len()),
                    miles(fronteras.len() - 1)
                ),
            ),
    );
    let mut peptidos = Vec::new();
    let mut combinacion = 0;
    for primera in 0..fronteras.len().saturating_sub(1) {
        for cortes_perdidos in 0..=reglas.cortes_perdidos {
            let ultima = primera + cortes_perdidos + 1;
            if ultima >= fronteras.len() {
                break;
            }
            let inicio = fronteras[primera];
            let fin = fronteras[ultima];
            let largo = fin - inicio;
            combinacion += 1;
            let detallar = combinacion <= COMBINACIONES_DETALLADAS && narracion.hay_lugar();
            if largo > reglas.largo_maximo {
                if detallar {
                    narracion.agregar(
                        pizarra_de_combinacion(proteina, primera, cortes_perdidos, inicio, fin, peptidos.len())
                            .paso(L_LARGO_MAXIMO, format!("Del pedazo {} con {}: {largo} letras, más de {}. Es demasiado largo y juntar más pedazos solo lo alarga: sales de este ciclo.", primera + 1, plural(cortes_perdidos, "corte perdido", "cortes perdidos"), reglas.largo_maximo)),
                    );
                }
                break;
            }
            if largo < reglas.largo_minimo {
                if detallar {
                    narracion.agregar(
                        pizarra_de_combinacion(proteina, primera, cortes_perdidos, inicio, fin, peptidos.len())
                            .paso(L_LARGO_MINIMO, format!("Del pedazo {} con {}: «{}», solo {largo} letras, menos de {}. Lo saltas y pruebas a juntar un pedazo más.", primera + 1, plural(cortes_perdidos, "corte perdido", "cortes perdidos"), como_texto(&proteina[inicio..fin]), reglas.largo_minimo)),
                    );
                }
                continue;
            }
            let Some(masa) =
                proteomica::masa_del_peptido(&proteina[inicio..fin], reglas.carbamidometilo)
            else {
                if detallar {
                    narracion.agregar(
                        pizarra_de_combinacion(proteina, primera, cortes_perdidos, inicio, fin, peptidos.len())
                            .paso(L_MASA_DE_RESIDUO, "Este pedazo tiene una letra que no es de los 20 aminoácidos (U, X…): masa_de_residuo devuelve None, el ? corta la suma y el péptido no se guarda."),
                    );
                }
                continue;
            };
            if peptidos.is_empty() && narracion.hay_lugar() {
                narracion.agregar(
                    pizarra_de_combinacion(proteina, primera, cortes_perdidos, inicio, fin, 0)
                        .paso(L_AGUA, format!("«{}» tiene {largo} letras: entra. Para su masa empiezas con la del agua.", recortado(&proteina[inicio..fin], LETRAS_EN_PIZARRA))),
                );
                if let Some(paso) = paso_de_masa(proteina, inicio, fin, reglas) {
                    narracion.agregar(paso);
                }
            }
            peptidos.push(PeptidoDigerido {
                inicio,
                fin,
                cortes_perdidos,
                masa,
            });
            let resumir = peptidos.len() % PEPTIDOS_RESUMIDOS_CADA == 0;
            if (detallar || resumir) && narracion.hay_lugar() {
                narracion.agregar(
                    pizarra_de_combinacion(
                        proteina,
                        primera,
                        cortes_perdidos,
                        inicio,
                        fin,
                        peptidos.len(),
                    )
                    .paso(
                        L_GUARDAR,
                        format!(
                            "Guardas «{}» ({largo} letras, {}, {} Da). Llevas {}.",
                            recortado(&proteina[inicio..fin], LETRAS_EN_PIZARRA),
                            plural(cortes_perdidos, "corte perdido", "cortes perdidos"),
                            decimal(masa, 4),
                            plural(peptidos.len(), "péptido", "péptidos")
                        ),
                    ),
                );
            }
        }
    }
    let lista = peptidos
        .iter()
        .take(FILAS_EN_PIZARRA * 2)
        .map(|peptido| {
            format!(
                "{} · {} Da",
                recortado(&proteina[peptido.inicio..peptido.fin], 30),
                decimal(peptido.masa, 4)
            )
        })
        .collect();
    let ultimo = Pizarra::nueva()
        .variable("peptidos", peptidos.len())
        .lista("Péptidos guardados", lista, &[])
        .paso(
            L_GUARDAR,
            format!(
                "Terminaste: {} de {} a {} letras. Mira «Qué significa» para la cobertura y los cortes perdidos.",
                plural(peptidos.len(), "péptido", "péptidos"),
                reglas.largo_minimo,
                reglas.largo_maximo
            ),
        );
    (narracion.terminar(ultimo), peptidos)
}

fn pasos_busqueda(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let pedido = leer_busqueda(entrada)?;
    let indice = indice_para(pedido.digestion)?;
    Ok(narrar_busqueda(&pedido, &indice).0)
}

fn filas_de_peptidos(proteina: &[u8], peptidos: &[PeptidoDigerido]) -> Vec<Vec<String>> {
    peptidos
        .iter()
        .take(FILAS_EN_PIZARRA)
        .map(|peptido| {
            vec![
                recortado(&proteina[peptido.inicio..peptido.fin], 30),
                peptido.cortes_perdidos.to_string(),
                decimal(peptido.masa, 4),
            ]
        })
        .collect()
}

fn narrar_indice(indice: &Indice, narracion: &mut Narracion) {
    let objetivos = indice.proteinas.len() / 2;
    let primera = &indice.proteinas[0];
    let su_senuelo = &indice.proteinas[objetivos];
    narracion.agregar(
        Pizarra::nueva()
            .variable("objetivos", objetivos)
            .variable("proteina", recortado(&primera.secuencia, LETRAS_EN_PIZARRA))
            .variable("senuelo", recortado(&su_senuelo.secuencia, LETRAS_EN_PIZARRA))
            .paso(
                L_SENUELOS,
                format!(
                    "Antes de ver espectros se arma el índice. Las {} proteínas de E. coli se copian y a cada una se le agrega su señuelo: la misma secuencia al revés. La primera es {} y su señuelo es «{}».",
                    miles(objetivos),
                    primera.nombre,
                    su_senuelo.nombre
                ),
            ),
    );
    for proteina in [primera, su_senuelo] {
        let peptidos = proteomica::digerir(&proteina.secuencia, &indice.digestion);
        narracion.agregar(
            Pizarra::nueva()
                .variable("proteina", recortado(&proteina.secuencia, LETRAS_EN_PIZARRA))
                .variable("peptidos", peptidos.len())
                .secuencia(
                    &proteina.secuencia,
                    &peptidos
                        .iter()
                        .map(|peptido| (peptido.inicio, peptido.fin, Marca::Hallado))
                        .collect::<Vec<_>>(),
                )
                .tabla(
                    &format!("Péptidos de {}", proteina.nombre),
                    &["Péptido", "Cortes perdidos", "Masa (Da)"],
                    filas_de_peptidos(&proteina.secuencia, &peptidos),
                    &[],
                )
                .paso(
                    L_DIGERIR_PROTEINA,
                    format!(
                        "Digieres {} con tripsina (después de K o R, no antes de P) y quedan {} de {} a {} letras, cada uno con su masa monoisotópica.",
                        proteina.nombre,
                        plural(peptidos.len(), "péptido", "péptidos"),
                        indice.digestion.largo_minimo,
                        indice.digestion.largo_maximo
                    ),
                ),
        );
    }
    let primera_masa = indice.peptidos.first().map_or(0.0, |p| p.masa);
    let ultima_masa = indice.peptidos.last().map_or(0.0, |p| p.masa);
    narracion.agregar(
        Pizarra::nueva()
            .variable("peptidos", indice.peptidos.len())
            .variable("masa_menor", decimal(primera_masa, 4))
            .variable("masa_mayor", decimal(ultima_masa, 4))
            .paso(
                L_ORDENAR_INDICE,
                format!(
                    "Lo mismo con las demás proteínas y sus señuelos, en paralelo: {} péptidos, ordenados por masa de {} a {} Da. Con la lista ordenada, encontrar los de una masa es una búsqueda binaria, no un recorrido.",
                    miles(indice.peptidos.len()),
                    decimal(primera_masa, 2),
                    decimal(ultima_masa, 2)
                ),
            ),
    );
}

fn filas_de_candidatos(
    indice: &Indice,
    candidatos: &[PeptidoIndexado],
    masa: f64,
) -> Vec<Vec<String>> {
    candidatos
        .iter()
        .take(FILAS_EN_PIZARRA)
        .map(|peptido| {
            vec![
                recortado(indice.secuencia(peptido), 30),
                nombre_de_proteina(indice, peptido),
                decimal(peptido.masa, 4),
                decimal(proteomica::error_ppm(masa, peptido.masa), 2),
            ]
        })
        .collect()
}

struct Elegido {
    peptido: PeptidoIndexado,
    carga: u8,
    puntaje: f64,
}

fn narrar_ventanas(
    espectro: &Espectro,
    picos: &[(f64, f64)],
    indice: &Indice,
    busqueda: &Busqueda,
    narracion: &mut Narracion,
) -> Option<Elegido> {
    let mut elegido: Option<Elegido> = None;
    for &carga in &espectro.cargas {
        let masa = proteomica::masa_neutra(espectro.mz_precursor, carga);
        narracion.agregar(
            Pizarra::nueva()
                .variable("carga", carga)
                .variable("masa", decimal(masa, 5))
                .paso(
                    L_MASA_NEUTRA,
                    format!(
                        "El precursor se midió en m/z {}. Si su carga es {carga}+, la masa del péptido sin protones es ({} − 1.00728) × {carga} = {} Da.",
                        decimal(espectro.mz_precursor, 5),
                        decimal(espectro.mz_precursor, 5),
                        decimal(masa, 5)
                    ),
                ),
        );
        let (minima, maxima) = proteomica::ventana(masa, busqueda.tolerancia_ppm);
        narracion.agregar(
            Pizarra::nueva()
                .variable("carga", carga)
                .variable("masa", decimal(masa, 5))
                .variable("minima", decimal(minima, 5))
                .variable("maxima", decimal(maxima, 5))
                .paso(
                    L_VENTANA,
                    format!(
                        "Con {} ppm de tolerancia la ventana va de {} a {} Da: apenas ±{} Da.",
                        decimal(busqueda.tolerancia_ppm, 1),
                        decimal(minima, 5),
                        decimal(maxima, 5),
                        decimal((maxima - minima) / 2.0, 5)
                    ),
                ),
        );
        let candidatos = indice.en_ventana(minima, maxima);
        narracion.agregar(
            Pizarra::nueva()
                .variable("carga", carga)
                .variable("candidatos", candidatos.len())
                .tabla(
                    "Péptidos del índice dentro de la ventana",
                    &["Péptido", "Proteína", "Masa (Da)", "Error (ppm)"],
                    filas_de_candidatos(indice, candidatos, masa),
                    &[],
                )
                .paso(
                    L_EN_VENTANA,
                    format!(
                        "Dos búsquedas binarias en el índice dan {} de {} péptidos. Solo esos se comparan con el espectro.",
                        plural(candidatos.len(), "candidato", "candidatos"),
                        miles(indice.peptidos.len())
                    ),
                ),
        );
        let mut b = Vec::new();
        let mut y = Vec::new();
        for (numero, peptido) in candidatos.iter().enumerate() {
            proteomica::llenar_iones(
                indice.secuencia(peptido),
                indice.digestion.carbamidometilo,
                &mut b,
                &mut y,
            );
            let puntaje = proteomica::puntuar(&b, &y, picos, busqueda.tolerancia_fragmento);
            let record = elegido.as_ref().map_or(0.0, |actual| actual.puntaje);
            let mejora = puntaje.hiperpuntaje > record;
            if mejora {
                elegido = Some(Elegido {
                    peptido: *peptido,
                    carga,
                    puntaje: puntaje.hiperpuntaje,
                });
            }
            if numero >= CANDIDATOS_DETALLADOS || !narracion.hay_lugar() {
                continue;
            }
            let texto = como_texto(indice.secuencia(peptido));
            let nota = if mejora {
                format!(
                    "«{texto}» ({}) explica {} b y {} y: hiperpuntaje {}. Supera al récord ({}), así que pasa a ser el mejor.",
                    nombre_de_proteina(indice, peptido),
                    puntaje.iones_b,
                    puntaje.iones_y,
                    decimal(puntaje.hiperpuntaje, 3),
                    decimal(record, 3)
                )
            } else {
                format!(
                    "«{texto}» ({}) explica {} b y {} y: hiperpuntaje {}. No supera al récord ({}).",
                    nombre_de_proteina(indice, peptido),
                    puntaje.iones_b,
                    puntaje.iones_y,
                    decimal(puntaje.hiperpuntaje, 3),
                    decimal(record, 3)
                )
            };
            narracion.agregar(
                Pizarra::nueva()
                    .variable("carga", carga)
                    .variable(
                        "peptido",
                        recortado(indice.secuencia(peptido), LETRAS_EN_PIZARRA),
                    )
                    .variable("iones_b", puntaje.iones_b)
                    .variable("iones_y", puntaje.iones_y)
                    .variable("hiperpuntaje", decimal(puntaje.hiperpuntaje, 3))
                    .variable("record", decimal(record.max(puntaje.hiperpuntaje), 3))
                    .paso(if mejora { L_RECORD } else { L_PUNTUAR }, nota),
            );
        }
        if candidatos.len() > CANDIDATOS_DETALLADOS {
            narracion.agregar(
                Pizarra::nueva()
                    .variable("carga", carga)
                    .variable("candidatos", candidatos.len())
                    .paso(
                        L_PUNTUAR,
                        format!(
                            "Lo mismo con los otros {} candidatos de esta carga.",
                            miles(candidatos.len() - CANDIDATOS_DETALLADOS)
                        ),
                    ),
            );
        }
    }
    elegido
}

fn narrar_iones(
    picos: &[(f64, f64)],
    indice: &Indice,
    elegido: &Elegido,
    busqueda: &Busqueda,
    narracion: &mut Narracion,
) -> f64 {
    let peptido = indice.secuencia(&elegido.peptido);
    let n = peptido.len();
    let (b, y) = iones_de(indice, &elegido.peptido);
    narracion.agregar(
        Pizarra::nueva()
            .variable("peptido", recortado(peptido, LETRAS_EN_PIZARRA))
            .variable("b", b.len())
            .variable("y", y.len())
            .secuencia(peptido, &[])
            .paso(
                L_IONES,
                format!(
                    "Ahora miras de cerca al ganador, «{}». Cuando el péptido se rompe en un enlace, el pedazo del inicio con un protón es un ion b y el del final con agua y un protón es un ion y: {} de cada uno.",
                    como_texto(peptido),
                    n - 1
                ),
            ),
    );
    let mut iones_b = 0;
    let mut iones_y = 0;
    let mut intensidad = 0.0;
    for (nombre, iones) in [("b", &b), ("y", &y)] {
        for (numero, &ion) in iones.iter().enumerate() {
            let largo = numero + 1;
            let tramo = if nombre == "b" {
                (0, largo)
            } else {
                (n - largo, n)
            };
            let cual = proteomica::pico_cercano(picos, ion, busqueda.tolerancia_fragmento);
            if let Some(indice_de_pico) = cual {
                intensidad += picos[indice_de_pico].1;
                if nombre == "b" {
                    iones_b += 1;
                } else {
                    iones_y += 1;
                }
            }
            if !narracion.hay_lugar() {
                continue;
            }
            let pizarra = Pizarra::nueva()
                .variable("ion", format!("{nombre}{largo}"))
                .variable("mz", decimal(ion, 4))
                .variable("iones_b", iones_b)
                .variable("iones_y", iones_y)
                .variable("intensidad", decimal(intensidad, 1))
                .secuencia(
                    peptido,
                    &[(
                        tramo.0,
                        tramo.1,
                        if cual.is_some() {
                            Marca::Hallado
                        } else {
                            Marca::Actual
                        },
                    )],
                );
            let paso = match cual {
                Some(indice_de_pico) => {
                    let (mz, altura) = picos[indice_de_pico];
                    pizarra.paso(L_SUMA, format!("{nombre}{largo} = «{}» debería estar en m/z {}. Hay un pico en {} (diferencia {} Da, dentro de {}): lo cuentas y sumas su intensidad {}.", como_texto(&peptido[tramo.0..tramo.1]), decimal(ion, 4), decimal(mz, 4), decimal(mz - ion, 4), decimal(busqueda.tolerancia_fragmento, 3), decimal(altura, 1)))
                }
                None => pizarra.paso(L_PICO, format!("{nombre}{largo} = «{}» debería estar en m/z {}, pero ningún pico cae a menos de {} Da: este ion no se vio.", como_texto(&peptido[tramo.0..tramo.1]), decimal(ion, 4), decimal(busqueda.tolerancia_fragmento, 3))),
            };
            narracion.agregar(paso);
        }
    }
    let puntaje = proteomica::hiperpuntaje(intensidad, iones_b, iones_y);
    narracion.agregar(
        Pizarra::nueva()
            .variable("iones_b", iones_b)
            .variable("iones_y", iones_y)
            .variable("intensidad", decimal(intensidad, 1))
            .variable("hiperpuntaje", decimal(puntaje, 3))
            .paso(
                L_HIPERPUNTAJE,
                format!(
                    "Hiperpuntaje = ln(1 + {}) + ln({iones_b}!) + ln({iones_y}!) = {}. Los factoriales premian las series largas de iones seguidos, que casi nunca salen por azar.",
                    decimal(intensidad, 1),
                    decimal(puntaje, 3)
                ),
            ),
    );
    puntaje
}

fn narrar_un_espectro(
    pedido: &PedidoDeBusqueda,
    indice: &Indice,
    numero: usize,
    narracion: &mut Narracion,
) {
    let espectro = &pedido.espectros[numero];
    let picos = proteomica::normalizar(&espectro.picos);
    let mas_altos: Vec<String> = {
        let mut orden: Vec<usize> = (0..picos.len()).collect();
        orden.sort_by(|&a, &c| picos[c].1.total_cmp(&picos[a].1));
        orden
            .iter()
            .take(FILAS_EN_PIZARRA)
            .map(|&cual| {
                format!(
                    "m/z {} · {}",
                    decimal(picos[cual].0, 4),
                    decimal(picos[cual].1, 1)
                )
            })
            .collect()
    };
    narracion.agregar(
        Pizarra::nueva()
            .variable("espectro", &espectro.titulo)
            .variable("picos", picos.len())
            .lista("Picos más altos (intensidad relativa)", mas_altos, &[])
            .paso(
                L_NORMALIZAR,
                format!(
                    "Miras en detalle el espectro «{}»: {}. Divides cada intensidad entre la mayor y la multiplicas por 100, para que los puntajes de espectros distintos se puedan comparar.",
                    espectro.titulo,
                    plural(picos.len(), "pico", "picos")
                ),
            ),
    );
    let Some(elegido) = narrar_ventanas(espectro, &picos, indice, &pedido.busqueda, narracion)
    else {
        narracion.agregar(Pizarra::nueva().paso(
            L_RECORD,
            "Ningún candidato explicó ni un ion: este espectro se queda sin pareja.",
        ));
        return;
    };
    let puntaje = narrar_iones(&picos, indice, &elegido, &pedido.busqueda, narracion);
    narracion.agregar(
        Pizarra::nueva()
            .variable("peptido", recortado(indice.secuencia(&elegido.peptido), LETRAS_EN_PIZARRA))
            .variable("carga", elegido.carga)
            .variable("hiperpuntaje", decimal(puntaje.max(elegido.puntaje), 3))
            .paso(
                L_RECORD,
                format!(
                    "Así queda el mejor PSM de este espectro: «{}» con carga {}+ y {} de proteína {}. Es la única pareja que este espectro lleva al cálculo del FDR.",
                    como_texto(indice.secuencia(&elegido.peptido)),
                    elegido.carga,
                    if indice.es_senuelo(&elegido.peptido) { "un señuelo" } else { "un objetivo" },
                    nombre_de_proteina(indice, &elegido.peptido)
                ),
            ),
    );
}

fn filas_de_psms(
    pedido: &PedidoDeBusqueda,
    indice: &Indice,
    psms: &[Psm],
    q: bool,
) -> Vec<Vec<String>> {
    psms.iter()
        .take(FILAS_EN_PIZARRA * 2)
        .map(|psm| {
            let mut fila = vec![
                recortado(pedido.espectros[psm.espectro].titulo.as_bytes(), 30),
                recortado(indice.secuencia(&psm.peptido), 30),
                decimal(psm.puntaje.hiperpuntaje, 2),
                tipo_de_psm(psm).to_string(),
            ];
            if q {
                fila.push(decimal(psm.q, 3));
            }
            fila
        })
        .collect()
}

fn narrar_valores_q(
    pedido: &PedidoDeBusqueda,
    indice: &Indice,
    psms: &[Psm],
    narracion: &mut Narracion,
) {
    narracion.agregar(
        Pizarra::nueva()
            .variable("psms", psms.len())
            .tabla(
                "PSM de mayor a menor puntaje",
                &["Espectro", "Péptido", "Hiperpuntaje", "Tipo"],
                filas_de_psms(pedido, indice, psms, false),
                &[],
            )
            .paso(
                L_ORDENAR_PSM,
                "Ordenas los PSM de mayor a menor hiperpuntaje. Ahora bajas un umbral por la lista: en cada punto cuentas cuántos objetivos y cuántos señuelos quedan arriba.",
            ),
    );
    let tasas = proteomica::tasas_por_umbral(psms);
    let mut objetivos = 0;
    let mut senuelos = 0;
    let cada = psms.len().div_ceil(PSM_DETALLADOS).max(1);
    for (lugar, psm) in psms.iter().enumerate() {
        if psm.senuelo {
            senuelos += 1;
        } else {
            objetivos += 1;
        }
        let interesante = (psm.senuelo && senuelos <= SENUELOS_DETALLADOS)
            || lugar < PSM_DETALLADOS / 2
            || lugar % cada == 0
            || lugar + 1 == psms.len();
        if !interesante || !narracion.hay_lugar() {
            continue;
        }
        let pizarra = Pizarra::nueva()
            .variable("lugar", lugar + 1)
            .variable("objetivos", objetivos)
            .variable("senuelos", senuelos)
            .variable("tasa", decimal(tasas[lugar], 3));
        let paso = if psm.senuelo {
            pizarra.paso(
                L_CONTAR_SENUELO,
                format!(
                    "Lugar {}: «{}» es un señuelo (puntaje {}). Ya van {} contra {}: FDR estimado = {senuelos}/{objetivos} = {}.",
                    lugar + 1,
                    como_texto(indice.secuencia(&psm.peptido)),
                    decimal(psm.puntaje.hiperpuntaje, 2),
                    plural(senuelos, "señuelo", "señuelos"),
                    plural(objetivos, "objetivo", "objetivos"),
                    decimal(tasas[lugar], 3)
                ),
            )
        } else {
            pizarra.paso(
                L_TASA,
                format!(
                    "Lugar {}: «{}» es un objetivo (puntaje {}). Con el corte aquí quedan {} y {}: FDR = {senuelos}/{objetivos} = {}.",
                    lugar + 1,
                    como_texto(indice.secuencia(&psm.peptido)),
                    decimal(psm.puntaje.hiperpuntaje, 2),
                    plural(objetivos, "objetivo", "objetivos"),
                    plural(senuelos, "señuelo", "señuelos"),
                    decimal(tasas[lugar], 3)
                ),
            )
        };
        narracion.agregar(paso);
    }
    let arrastrados = psms
        .iter()
        .zip(&tasas)
        .filter(|(psm, tasa)| psm.q < **tasa)
        .count();
    narracion.agregar(
        Pizarra::nueva()
            .variable("psms", psms.len())
            .tabla(
                "q-valor de cada PSM",
                &["Espectro", "Péptido", "Hiperpuntaje", "Tipo", "q"],
                filas_de_psms(pedido, indice, psms, true),
                &[],
            )
            .paso(
                L_MINIMO,
                format!(
                    "Ahora recorres la lista desde abajo guardando el menor FDR visto: ese mínimo es el q-valor, el FDR más bajo de cualquier corte que incluya al PSM. En {} el q-valor quedó por debajo de su FDR.",
                    plural(arrastrados, "PSM", "PSM")
                ),
            ),
    );
}

fn narrar_busqueda(pedido: &PedidoDeBusqueda, indice: &Indice) -> (Vec<Paso>, Vec<Psm>) {
    let mut narracion = Narracion::nueva(400);
    narracion.agregar(
        Pizarra::nueva()
            .variable("espectros", pedido.espectros.len())
            .variable("tolerancia_ppm", decimal(pedido.busqueda.tolerancia_ppm, 1))
            .variable("tolerancia_fragmento", decimal(pedido.busqueda.tolerancia_fragmento, 3))
            .paso(
                L_BUSCAR,
                format!(
                    "Tienes {}. Meta: para cada uno, el péptido de E. coli que mejor lo explica, y luego cuántas de esas parejas podrían ser azar.",
                    plural(pedido.espectros.len(), "espectro MS/MS", "espectros MS/MS")
                ),
            ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("tolerancia_ppm", decimal(pedido.busqueda.tolerancia_ppm, 1))
            .variable("tolerancia_fragmento", decimal(pedido.busqueda.tolerancia_fragmento, 3))
            .paso(
                L_REVISAR,
                format!(
                    "Las tolerancias son válidas: {} ppm para el precursor ({}) y {} Da para los fragmentos.",
                    decimal(pedido.busqueda.tolerancia_ppm, 1),
                    if pedido.busqueda.tolerancia_ppm <= PPM_DE_ALTA_RESOLUCION {
                        "lo que mide un espectrómetro de alta resolución"
                    } else {
                        "una ventana ancha: entran más candidatos y más azar"
                    },
                    decimal(pedido.busqueda.tolerancia_fragmento, 3)
                ),
            ),
    );
    narrar_indice(indice, &mut narracion);
    let mut psms: Vec<Psm> = pedido
        .espectros
        .par_iter()
        .enumerate()
        .filter_map(|(numero, espectro)| {
            proteomica::mejor_psm(numero, espectro, indice, &pedido.busqueda)
        })
        .collect();
    let detallado = psms
        .iter()
        .max_by(|a, c| {
            a.puntaje
                .hiperpuntaje
                .total_cmp(&c.puntaje.hiperpuntaje)
                .then(c.espectro.cmp(&a.espectro))
        })
        .map_or(0, |psm| psm.espectro);
    narracion.agregar(
        Pizarra::nueva()
            .variable("espectros", pedido.espectros.len())
            .paso(
                L_CADA_ESPECTRO,
                format!(
                    "Cada espectro busca su mejor pareja por separado (en paralelo). Te muestro uno completo: «{}».",
                    pedido.espectros[detallado].titulo
                ),
            ),
    );
    narrar_un_espectro(pedido, indice, detallado, &mut narracion);
    narracion.agregar(
        Pizarra::nueva()
            .variable("psms", psms.len())
            .tabla(
                "Mejor pareja de cada espectro",
                &["Espectro", "Péptido", "Hiperpuntaje", "Tipo"],
                filas_de_psms(pedido, indice, &psms, false),
                &[],
            )
            .paso(
                L_CADA_ESPECTRO,
                format!(
                    "Lo mismo con todos: {} de {} encontraron pareja ({} con objetivos, {} con señuelos).",
                    miles(psms.len()),
                    plural(pedido.espectros.len(), "espectro", "espectros"),
                    miles(psms.iter().filter(|psm| !psm.senuelo).count()),
                    miles(psms.iter().filter(|psm| psm.senuelo).count())
                ),
            ),
    );
    proteomica::asignar_valores_q(&mut psms);
    narrar_valores_q(pedido, indice, &psms, &mut narracion);
    let aceptados = proteomica::aceptados(&psms, FDR_ESTANDAR);
    let lista = aceptados
        .iter()
        .take(FILAS_EN_PIZARRA * 2)
        .map(|psm| {
            format!(
                "{} · {}",
                como_texto(indice.secuencia(&psm.peptido)),
                nombre_de_proteina(indice, &psm.peptido)
            )
        })
        .collect();
    let ultimo = Pizarra::nueva()
        .variable("aceptados", aceptados.len())
        .lista("Aceptados con q ≤ 0.01", lista, &[])
        .paso(
            L_Q,
            format!(
                "Te quedas con los objetivos de q ≤ 0.01: {}. Los señuelos nunca se reportan. Mira «Qué significa» para leer el resultado.",
                plural(aceptados.len(), "PSM aceptado", "PSM aceptados")
            ),
        );
    (narracion.terminar(ultimo), psms)
}

#[cfg(test)]
mod pruebas {
    use super::*;

    fn codigo_de(id: &str) -> String {
        crate::codigo(id)
            .unwrap()
            .iter()
            .map(|fragmento| fragmento.fuente.clone())
            .collect::<Vec<_>>()
            .join("\n")
    }

    fn revisar_lineas(id: &str, pasos: &[Paso]) {
        let codigo = codigo_de(id);
        for paso in pasos {
            assert!(
                codigo
                    .lines()
                    .any(|linea| linea.contains(paso.linea.as_str())),
                "{id}: {}",
                paso.linea
            );
        }
    }

    fn ejemplo_de_busqueda() -> Entrada {
        Entrada::nueva(&[
            ("espectros", ESPECTROS_DE_EJEMPLO),
            ("archivo", ""),
            ("tolerancia_ppm", "10"),
            ("tolerancia_fragmento", "0.02"),
            ("cortes_perdidos", "1"),
            ("cisteinas", "carbamidometiladas"),
        ])
    }

    #[test]
    fn la_digestion_narrada_da_lo_mismo_que_el_nucleo() {
        let mut azar = nucleo::azar::Azar::nuevo(5);
        for gen_de_ecoli in ["recA", "thrL", "lacZ", "tufA", "ompA"] {
            let proteina = ecoli::buscar_gen(gen_de_ecoli).unwrap().proteina.clone();
            for cortes in 0..=3 {
                let reglas = Digestion {
                    cortes_perdidos: cortes,
                    carbamidometilo: azar.entero(2) == 0,
                    ..Digestion::default()
                };
                let (pasos, peptidos) = narrar_digestion(&proteina, &reglas);
                assert_eq!(peptidos, proteomica::digerir(&proteina, &reglas));
                assert!(pasos.len() >= 5 && pasos.len() <= 400);
                revisar_lineas("antibioticos.digestion", &pasos);
            }
        }
        for _ in 0..50 {
            let largo = 1 + azar.entero(300);
            let proteina: Vec<u8> = (0..largo)
                .map(|_| b"KRPAGLUSTDE"[azar.entero(11)])
                .collect();
            let reglas = Digestion {
                cortes_perdidos: azar.entero(4),
                ..Digestion::default()
            };
            let (pasos, peptidos) = narrar_digestion(&proteina, &reglas);
            assert_eq!(peptidos, proteomica::digerir(&proteina, &reglas));
            revisar_lineas("antibioticos.digestion", &pasos);
        }
    }

    #[test]
    fn la_busqueda_narrada_da_lo_mismo_que_el_nucleo() {
        let pedido = leer_busqueda(&ejemplo_de_busqueda()).unwrap();
        let indice = indice_para(pedido.digestion).unwrap();
        let (pasos, psms) = narrar_busqueda(&pedido, &indice);
        let directos =
            proteomica::buscar_espectros(&pedido.espectros, &indice, &pedido.busqueda).unwrap();
        assert_eq!(psms, directos);
        assert!(pasos.len() >= 30 && pasos.len() <= 400, "{}", pasos.len());
        revisar_lineas("antibioticos.busqueda_de_espectros", &pasos);
    }

    #[test]
    fn el_ejemplo_identifica_los_diez_peptidos_y_rechaza_el_ajeno() {
        let pedido = leer_busqueda(&ejemplo_de_busqueda()).unwrap();
        let indice = indice_para(pedido.digestion).unwrap();
        let psms =
            proteomica::buscar_espectros(&pedido.espectros, &indice, &pedido.busqueda).unwrap();
        let mut aceptados: Vec<String> = proteomica::aceptados(&psms, FDR_ESTANDAR)
            .iter()
            .map(|psm| como_texto(indice.secuencia(&psm.peptido)))
            .collect();
        aceptados.sort();
        let mut esperados = vec![
            "ALAAALGQIEK",
            "AVAAVNGPIAQALIGK",
            "DWLQAIESVIR",
            "FDGTVEVKDGHLIVNGK",
            "GVNVLADAVK",
            "LGEPVFDVQECQIR",
            "LGYPITDDLDIYTR",
            "QAVTNPQNTLFAIK",
            "TCAFIDAEHALDPIYAR",
            "TTLTLQVIAAAQR",
        ];
        esperados.sort();
        assert_eq!(aceptados, esperados);
        let ajeno = pedido
            .espectros
            .iter()
            .position(|espectro| espectro.titulo == "corrida_1.1570.1570.2")
            .unwrap();
        assert!(
            psms.iter()
                .filter(|psm| psm.espectro == ajeno)
                .all(|psm| !pasa(psm))
        );
    }

    fn espectro_exacto(peptido: &[u8], carga: u8, titulo: &str) -> String {
        let (mut b, mut y) = (Vec::new(), Vec::new());
        proteomica::llenar_iones(peptido, true, &mut b, &mut y);
        let masa = proteomica::masa_del_peptido(peptido, true).unwrap();
        let mut texto = format!(
            "BEGIN IONS\nTITLE={titulo}\nPEPMASS={:.5}\nCHARGE={carga}+\n",
            proteomica::mz_con_carga(masa, carga)
        );
        for ion in b.iter().chain(&y) {
            texto.push_str(&format!("{ion:.4} 50\n"));
        }
        texto.push_str("END IONS\n");
        texto
    }

    fn busqueda_con(texto: &str, ppm: &str) -> Entrada {
        Entrada::nueva(&[
            ("espectros", texto),
            ("archivo", ""),
            ("tolerancia_ppm", ppm),
            ("tolerancia_fragmento", "0.02"),
            ("cortes_perdidos", "1"),
            ("cisteinas", "carbamidometiladas"),
        ])
    }

    #[test]
    fn el_texto_concuerda_en_singular_con_un_solo_espectro() {
        let texto = espectro_exacto(b"TCAFIDAEHALDPIYAR", 3, "unico");
        let salida = busqueda_de_espectros(&busqueda_con(&texto, "10")).unwrap();
        let todo = salida.interpretacion.join(" ");
        assert!(todo.contains("De 1 espectro, 1 pasa el corte"), "{todo}");
        assert!(todo.contains("0 espectros tuvieron"), "{todo}");
        assert!(todo.contains("0 objetivos quedaron"), "{todo}");
        assert!(!todo.contains("1 de los aceptados tienen"), "{todo}");
        assert!(todo.contains("1 de los aceptados tiene I o L"), "{todo}");
    }

    #[test]
    fn una_ventana_ancha_no_se_llama_adecuada() {
        let texto = espectro_exacto(b"TCAFIDAEHALDPIYAR", 3, "unico");
        let salida = busqueda_de_espectros(&busqueda_con(&texto, "200")).unwrap();
        let todo = salida.interpretacion.join(" ");
        assert!(todo.contains("más ancha de lo necesario"), "{todo}");
        assert!(!todo.contains("es adecuada"), "{todo}");
    }

    #[test]
    fn el_verbo_sigue_al_numero() {
        assert_eq!(verbo(1, "pasa", "pasan"), "pasa");
        assert_eq!(verbo(0, "pasa", "pasan"), "pasan");
        assert_eq!(verbo(7, "pasa", "pasan"), "pasan");
    }

    #[test]
    fn el_ejemplo_de_digestion_es_razonable() {
        let entrada = Entrada::nueva(&[
            ("proteina", "@proteina:recA"),
            ("cortes_perdidos", "1"),
            ("cisteinas", "carbamidometiladas"),
        ]);
        let salida = digestion(&entrada).unwrap();
        assert!(salida.interpretacion.len() >= 2 && salida.interpretacion.len() <= 6);
        let pasos = pasos_digestion(&entrada).unwrap();
        assert!(pasos.len() >= 15 && pasos.len() <= 400);
    }

    #[test]
    fn entradas_raras_dan_errores_en_espanol() {
        let malo = Entrada::nueva(&[
            ("espectros", "BEGIN IONS\nPEPMASS=500\n"),
            ("archivo", ""),
            ("tolerancia_ppm", "10"),
            ("tolerancia_fragmento", "0.02"),
            ("cortes_perdidos", "1"),
            ("cisteinas", "libres"),
        ]);
        assert!(busqueda_de_espectros(&malo).is_err());
        let sin_picos = Entrada::nueva(&[
            ("espectros", "BEGIN IONS\nPEPMASS=500\nEND IONS"),
            ("archivo", ""),
            ("tolerancia_ppm", "10"),
            ("tolerancia_fragmento", "0.02"),
            ("cortes_perdidos", "0"),
            ("cisteinas", "libres"),
        ]);
        let salida = busqueda_de_espectros(&sin_picos).unwrap();
        assert!(!salida.interpretacion.is_empty());
        let pasos = pasos_busqueda(&sin_picos).unwrap();
        assert!(!pasos.is_empty());
        let larga = Entrada::nueva(&[
            ("proteina", "@ecoli"),
            ("cortes_perdidos", "1"),
            ("cisteinas", "libres"),
        ]);
        assert!(digestion(&larga).unwrap_err().contains("máximo"));
    }
}
