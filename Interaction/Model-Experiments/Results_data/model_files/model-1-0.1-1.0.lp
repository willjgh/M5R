\ Model birth-death-regulation-capture-efficiency
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  k_reg
Subject To
 R0: k_deg_2 = 1
 B_lb_0,_0: p[0,0] >= 0.0904
 B_ub_0,_0: p[0,0] <= 0.1716
 B_lb_0,_1: p[0,1] >= 0.08
 B_ub_0,_1: p[0,1] <= 0.1606
 B_lb_0,_2: p[0,2] >= -0.0052
 B_ub_0,_2: p[0,2] <= 0.069
 B_lb_0,_3: p[0,3] >= -0.0436
 B_ub_0,_3: p[0,3] <= 0.0254
 B_lb_0,_4: p[0,4] >= -0.057
 B_ub_0,_4: p[0,4] <= 0.008
 B_lb_0,_5: p[0,5] >= -0.0608
 B_ub_0,_5: p[0,5] <= 1
 B_lb_0,_6: p[0,6] >= -0.0608
 B_ub_0,_6: p[0,6] <= 1
 B_lb_1,_0: p[1,0] >= 0.079
 B_ub_1,_0: p[1,0] <= 0.1592
 B_lb_1,_1: p[1,1] >= 0.0596
 B_ub_1,_1: p[1,1] <= 0.1388
 B_lb_1,_2: p[1,2] >= -0.0088
 B_ub_1,_2: p[1,2] <= 0.0646
 B_lb_1,_3: p[1,3] >= -0.044
 B_ub_1,_3: p[1,3] <= 0.0246
 B_lb_1,_4: p[1,4] >= -0.0588
 B_ub_1,_4: p[1,4] <= 0.0054
 B_lb_1,_5: p[1,5] >= -0.0608
 B_ub_1,_5: p[1,5] <= 1
 B_lb_1,_6: p[1,6] >= -0.0608
 B_ub_1,_6: p[1,6] <= 1
 B_lb_2,_0: p[2,0] >= -0.0002
 B_ub_2,_0: p[2,0] <= 0.0746
 B_lb_2,_1: p[2,1] >= -0.0044
 B_ub_2,_1: p[2,1] <= 0.07
 B_lb_2,_2: p[2,2] >= -0.0376
 B_ub_2,_2: p[2,2] <= 0.0324
 B_lb_2,_3: p[2,3] >= -0.0556
 B_ub_2,_3: p[2,3] <= 0.0098
 B_lb_2,_4: p[2,4] >= -0.0608
 B_ub_2,_4: p[2,4] <= 1
 B_lb_2,_5: p[2,5] >= -0.0608
 B_ub_2,_5: p[2,5] <= 1
 B_lb_2,_6: p[2,6] >= -0.0608
 B_ub_2,_6: p[2,6] <= 1
 B_lb_3,_0: p[3,0] >= -0.0414
 B_ub_3,_0: p[3,0] <= 0.0278
 B_lb_3,_1: p[3,1] >= -0.0456
 B_ub_3,_1: p[3,1] <= 0.0226
 B_lb_3,_2: p[3,2] >= -0.0548
 B_ub_3,_2: p[3,2] <= 0.011
 B_lb_3,_3: p[3,3] >= -0.0596
 B_ub_3,_3: p[3,3] <= 0.0038
 B_lb_3,_4: p[3,4] >= -0.0608
 B_ub_3,_4: p[3,4] <= 1
 B_lb_3,_5: p[3,5] >= -0.0608
 B_ub_3,_5: p[3,5] <= 1
 B_lb_3,_6: p[3,6] >= -0.0608
 B_ub_3,_6: p[3,6] <= 1
 B_lb_4,_0: p[4,0] >= -0.058
 B_ub_4,_0: p[4,0] <= 0.0066
 B_lb_4,_1: p[4,1] >= -0.059
 B_ub_4,_1: p[4,1] <= 0.0052
 B_lb_4,_2: p[4,2] >= -0.06
 B_ub_4,_2: p[4,2] <= 0.0034
 B_lb_4,_3: p[4,3] >= -0.0608
 B_ub_4,_3: p[4,3] <= 1
 B_lb_4,_4: p[4,4] >= -0.0608
 B_ub_4,_4: p[4,4] <= 1
 B_lb_4,_5: p[4,5] >= -0.0608
 B_ub_4,_5: p[4,5] <= 1
 B_lb_4,_6: p[4,6] >= -0.0608
 B_ub_4,_6: p[4,6] <= 1
 B_lb_5,_0: p[5,0] >= -0.0608
 B_ub_5,_0: p[5,0] <= 1
 B_lb_5,_1: p[5,1] >= -0.0608
 B_ub_5,_1: p[5,1] <= 1
 B_lb_5,_2: p[5,2] >= -0.0608
 B_ub_5,_2: p[5,2] <= 1
 B_lb_5,_3: p[5,3] >= -0.0608
 B_ub_5,_3: p[5,3] <= 1
 B_lb_5,_4: p[5,4] >= -0.0608
 B_ub_5,_4: p[5,4] <= 1
 B_lb_5,_5: p[5,5] >= -0.0608
 B_ub_5,_5: p[5,5] <= 1
 B_lb_5,_6: p[5,6] >= -0.0608
 B_ub_5,_6: p[5,6] <= 1
 B_marginal_lb_0: p1[0] >= 0.232495
 B_marginal_ub_0: p1[0] <= 0.416
 B_marginal_lb_1: p1[1] >= 0.1935
 B_marginal_ub_1: p1[1] <= 0.3762
 B_marginal_lb_2: p1[2] >= 0.0011
 B_marginal_ub_2: p1[2] <= 0.178
 B_marginal_lb_3: p1[3] >= -0.1095
 B_marginal_ub_3: p1[3] <= 0.0598
 B_marginal_lb_4: p1[4] >= -0.1489
 B_marginal_ub_4: p1[4] <= 0.0136
 B_marginal_lb_5: p1[5] >= -0.1559
 B_marginal_ub_5: p1[5] <= 0.0034
 B_marginal_lb_0: p2[0] >= 0.2612857142857143
 B_marginal_ub_0: p2[0] <= 0.4214
 B_marginal_lb_1: p2[1] >= 0.2194857142857143
 B_marginal_ub_1: p2[1] <= 0.3788
 B_marginal_lb_2: p2[2] >= 0.0156857142857143
 B_marginal_ub_2: p2[2] <= 0.169
 B_marginal_lb_3: p2[3] >= -0.0869142857142857
 B_marginal_ub_3: p2[3] <= 0.0588
 B_marginal_lb_4: p2[4] >= -0.1253142857142857
 B_marginal_ub_4: p2[4] <= 0.0138
 B_marginal_lb_5: p2[5] >= -0.1313142857142857
 B_marginal_ub_5: p2[5] <= 0.0054
 B_marginal_lb_6: p2[6] >= -0.1333142857142857
 B_marginal_ub_6: p2[6] <= 1
 frechet_lb_0,_0: - p[0,0] + p1[0] + p2[0] <= 1
 frechet_ub_x1_0: p[0,0] - p1[0] <= 0
 frechet_ub_x2_0: p[0,0] - p2[0] <= 0
 frechet_lb_0,_1: - p[0,1] + p1[0] + p2[1] <= 1
 frechet_ub_x1_0: p[0,1] - p1[0] <= 0
 frechet_ub_x2_1: p[0,1] - p2[1] <= 0
 frechet_lb_0,_2: - p[0,2] + p1[0] + p2[2] <= 1
 frechet_ub_x1_0: p[0,2] - p1[0] <= 0
 frechet_ub_x2_2: p[0,2] - p2[2] <= 0
 frechet_lb_0,_3: - p[0,3] + p1[0] + p2[3] <= 1
 frechet_ub_x1_0: p[0,3] - p1[0] <= 0
 frechet_ub_x2_3: p[0,3] - p2[3] <= 0
 frechet_lb_0,_4: - p[0,4] + p1[0] + p2[4] <= 1
 frechet_ub_x1_0: p[0,4] - p1[0] <= 0
 frechet_ub_x2_4: p[0,4] - p2[4] <= 0
 frechet_lb_0,_5: - p[0,5] + p1[0] + p2[5] <= 1
 frechet_ub_x1_0: p[0,5] - p1[0] <= 0
 frechet_ub_x2_5: p[0,5] - p2[5] <= 0
 frechet_lb_0,_6: - p[0,6] + p1[0] + p2[6] <= 1
 frechet_ub_x1_0: p[0,6] - p1[0] <= 0
 frechet_ub_x2_6: p[0,6] - p2[6] <= 0
 frechet_lb_1,_0: - p[1,0] + p1[1] + p2[0] <= 1
 frechet_ub_x1_1: p[1,0] - p1[1] <= 0
 frechet_ub_x2_0: p[1,0] - p2[0] <= 0
 frechet_lb_1,_1: - p[1,1] + p1[1] + p2[1] <= 1
 frechet_ub_x1_1: p[1,1] - p1[1] <= 0
 frechet_ub_x2_1: p[1,1] - p2[1] <= 0
 frechet_lb_1,_2: - p[1,2] + p1[1] + p2[2] <= 1
 frechet_ub_x1_1: p[1,2] - p1[1] <= 0
 frechet_ub_x2_2: p[1,2] - p2[2] <= 0
 frechet_lb_1,_3: - p[1,3] + p1[1] + p2[3] <= 1
 frechet_ub_x1_1: p[1,3] - p1[1] <= 0
 frechet_ub_x2_3: p[1,3] - p2[3] <= 0
 frechet_lb_1,_4: - p[1,4] + p1[1] + p2[4] <= 1
 frechet_ub_x1_1: p[1,4] - p1[1] <= 0
 frechet_ub_x2_4: p[1,4] - p2[4] <= 0
 frechet_lb_1,_5: - p[1,5] + p1[1] + p2[5] <= 1
 frechet_ub_x1_1: p[1,5] - p1[1] <= 0
 frechet_ub_x2_5: p[1,5] - p2[5] <= 0
 frechet_lb_1,_6: - p[1,6] + p1[1] + p2[6] <= 1
 frechet_ub_x1_1: p[1,6] - p1[1] <= 0
 frechet_ub_x2_6: p[1,6] - p2[6] <= 0
 frechet_lb_2,_0: - p[2,0] + p1[2] + p2[0] <= 1
 frechet_ub_x1_2: p[2,0] - p1[2] <= 0
 frechet_ub_x2_0: p[2,0] - p2[0] <= 0
 frechet_lb_2,_1: - p[2,1] + p1[2] + p2[1] <= 1
 frechet_ub_x1_2: p[2,1] - p1[2] <= 0
 frechet_ub_x2_1: p[2,1] - p2[1] <= 0
 frechet_lb_2,_2: - p[2,2] + p1[2] + p2[2] <= 1
 frechet_ub_x1_2: p[2,2] - p1[2] <= 0
 frechet_ub_x2_2: p[2,2] - p2[2] <= 0
 frechet_lb_2,_3: - p[2,3] + p1[2] + p2[3] <= 1
 frechet_ub_x1_2: p[2,3] - p1[2] <= 0
 frechet_ub_x2_3: p[2,3] - p2[3] <= 0
 frechet_lb_2,_4: - p[2,4] + p1[2] + p2[4] <= 1
 frechet_ub_x1_2: p[2,4] - p1[2] <= 0
 frechet_ub_x2_4: p[2,4] - p2[4] <= 0
 frechet_lb_2,_5: - p[2,5] + p1[2] + p2[5] <= 1
 frechet_ub_x1_2: p[2,5] - p1[2] <= 0
 frechet_ub_x2_5: p[2,5] - p2[5] <= 0
 frechet_lb_2,_6: - p[2,6] + p1[2] + p2[6] <= 1
 frechet_ub_x1_2: p[2,6] - p1[2] <= 0
 frechet_ub_x2_6: p[2,6] - p2[6] <= 0
 frechet_lb_3,_0: - p[3,0] + p1[3] + p2[0] <= 1
 frechet_ub_x1_3: p[3,0] - p1[3] <= 0
 frechet_ub_x2_0: p[3,0] - p2[0] <= 0
 frechet_lb_3,_1: - p[3,1] + p1[3] + p2[1] <= 1
 frechet_ub_x1_3: p[3,1] - p1[3] <= 0
 frechet_ub_x2_1: p[3,1] - p2[1] <= 0
 frechet_lb_3,_2: - p[3,2] + p1[3] + p2[2] <= 1
 frechet_ub_x1_3: p[3,2] - p1[3] <= 0
 frechet_ub_x2_2: p[3,2] - p2[2] <= 0
 frechet_lb_3,_3: - p[3,3] + p1[3] + p2[3] <= 1
 frechet_ub_x1_3: p[3,3] - p1[3] <= 0
 frechet_ub_x2_3: p[3,3] - p2[3] <= 0
 frechet_lb_3,_4: - p[3,4] + p1[3] + p2[4] <= 1
 frechet_ub_x1_3: p[3,4] - p1[3] <= 0
 frechet_ub_x2_4: p[3,4] - p2[4] <= 0
 frechet_lb_3,_5: - p[3,5] + p1[3] + p2[5] <= 1
 frechet_ub_x1_3: p[3,5] - p1[3] <= 0
 frechet_ub_x2_5: p[3,5] - p2[5] <= 0
 frechet_lb_3,_6: - p[3,6] + p1[3] + p2[6] <= 1
 frechet_ub_x1_3: p[3,6] - p1[3] <= 0
 frechet_ub_x2_6: p[3,6] - p2[6] <= 0
 frechet_lb_4,_0: - p[4,0] + p1[4] + p2[0] <= 1
 frechet_ub_x1_4: p[4,0] - p1[4] <= 0
 frechet_ub_x2_0: p[4,0] - p2[0] <= 0
 frechet_lb_4,_1: - p[4,1] + p1[4] + p2[1] <= 1
 frechet_ub_x1_4: p[4,1] - p1[4] <= 0
 frechet_ub_x2_1: p[4,1] - p2[1] <= 0
 frechet_lb_4,_2: - p[4,2] + p1[4] + p2[2] <= 1
 frechet_ub_x1_4: p[4,2] - p1[4] <= 0
 frechet_ub_x2_2: p[4,2] - p2[2] <= 0
 frechet_lb_4,_3: - p[4,3] + p1[4] + p2[3] <= 1
 frechet_ub_x1_4: p[4,3] - p1[4] <= 0
 frechet_ub_x2_3: p[4,3] - p2[3] <= 0
 frechet_lb_4,_4: - p[4,4] + p1[4] + p2[4] <= 1
 frechet_ub_x1_4: p[4,4] - p1[4] <= 0
 frechet_ub_x2_4: p[4,4] - p2[4] <= 0
 frechet_lb_4,_5: - p[4,5] + p1[4] + p2[5] <= 1
 frechet_ub_x1_4: p[4,5] - p1[4] <= 0
 frechet_ub_x2_5: p[4,5] - p2[5] <= 0
 frechet_lb_4,_6: - p[4,6] + p1[4] + p2[6] <= 1
 frechet_ub_x1_4: p[4,6] - p1[4] <= 0
 frechet_ub_x2_6: p[4,6] - p2[6] <= 0
 frechet_lb_5,_0: - p[5,0] + p1[5] + p2[0] <= 1
 frechet_ub_x1_5: p[5,0] - p1[5] <= 0
 frechet_ub_x2_0: p[5,0] - p2[0] <= 0
 frechet_lb_5,_1: - p[5,1] + p1[5] + p2[1] <= 1
 frechet_ub_x1_5: p[5,1] - p1[5] <= 0
 frechet_ub_x2_1: p[5,1] - p2[1] <= 0
 frechet_lb_5,_2: - p[5,2] + p1[5] + p2[2] <= 1
 frechet_ub_x1_5: p[5,2] - p1[5] <= 0
 frechet_ub_x2_2: p[5,2] - p2[2] <= 0
 frechet_lb_5,_3: - p[5,3] + p1[5] + p2[3] <= 1
 frechet_ub_x1_5: p[5,3] - p1[5] <= 0
 frechet_ub_x2_3: p[5,3] - p2[3] <= 0
 frechet_lb_5,_4: - p[5,4] + p1[5] + p2[4] <= 1
 frechet_ub_x1_5: p[5,4] - p1[5] <= 0
 frechet_ub_x2_4: p[5,4] - p2[4] <= 0
 frechet_lb_5,_5: - p[5,5] + p1[5] + p2[5] <= 1
 frechet_ub_x1_5: p[5,5] - p1[5] <= 0
 frechet_ub_x2_5: p[5,5] - p2[5] <= 0
 frechet_lb_5,_6: - p[5,6] + p1[5] + p2[6] <= 1
 frechet_ub_x1_5: p[5,6] - p1[5] <= 0
 frechet_ub_x2_6: p[5,6] - p2[6] <= 0
 Distribution: p[0,0] + p[0,1] + p[0,2] + p[0,3] + p[0,4] + p[0,5] + p[0,6]
   + p[1,0] + p[1,1] + p[1,2] + p[1,3] + p[1,4] + p[1,5] + p[1,6] + p[2,0]
   + p[2,1] + p[2,2] + p[2,3] + p[2,4] + p[2,5] + p[2,6] + p[3,0] + p[3,1]
   + p[3,2] + p[3,3] + p[3,4] + p[3,5] + p[3,6] + p[4,0] + p[4,1] + p[4,2]
   + p[4,3] + p[4,4] + p[4,5] + p[4,6] + p[5,0] + p[5,1] + p[5,2] + p[5,3]
   + p[5,4] + p[5,5] + p[5,6] <= 1
 Distribution_x1: p1[0] + p1[1] + p1[2] + p1[3] + p1[4] + p1[5] <= 1
 Distribution_x2: p2[0] + p2[1] + p2[2] + p2[3] + p2[4] + p2[5] + p2[6]
   <= 1
 Equation_0,_0: [ - p[0,0] * k_tx_1 - p[0,0] * k_tx_2 + p[0,1] * k_deg_2
   + p[1,0] * k_deg_1 + p[1,1] * k_reg ] = 0
 Equation_0,_1: [ p[0,0] * k_tx_2 - p[0,1] * k_tx_1 - p[0,1] * k_tx_2
   - p[0,1] * k_deg_2 + 2 p[0,2] * k_deg_2 + p[1,1] * k_deg_1
   + 2 p[1,2] * k_reg ] = 0
 Equation_0,_2: [ p[0,1] * k_tx_2 - p[0,2] * k_tx_1 - p[0,2] * k_tx_2
   - 2 p[0,2] * k_deg_2 + 3 p[0,3] * k_deg_2 + p[1,2] * k_deg_1
   + 3 p[1,3] * k_reg ] = 0
 Equation_0,_3: [ p[0,2] * k_tx_2 - p[0,3] * k_tx_1 - p[0,3] * k_tx_2
   - 3 p[0,3] * k_deg_2 + 4 p[0,4] * k_deg_2 + p[1,3] * k_deg_1
   + 4 p[1,4] * k_reg ] = 0
 Equation_0,_4: [ p[0,3] * k_tx_2 - p[0,4] * k_tx_1 - p[0,4] * k_tx_2
   - 4 p[0,4] * k_deg_2 + 5 p[0,5] * k_deg_2 + p[1,4] * k_deg_1
   + 5 p[1,5] * k_reg ] = 0
 Equation_0,_5: [ p[0,4] * k_tx_2 - p[0,5] * k_tx_1 - p[0,5] * k_tx_2
   - 5 p[0,5] * k_deg_2 + 6 p[0,6] * k_deg_2 + p[1,5] * k_deg_1
   + 6 p[1,6] * k_reg ] = 0
 Equation_1,_0: [ p[0,0] * k_tx_1 - p[1,0] * k_tx_1 - p[1,0] * k_tx_2
   - p[1,0] * k_deg_1 + p[1,1] * k_deg_2 + 2 p[2,0] * k_deg_1
   + 2 p[2,1] * k_reg ] = 0
 Equation_1,_1: [ p[0,1] * k_tx_1 + p[1,0] * k_tx_2 - p[1,1] * k_tx_1
   - p[1,1] * k_tx_2 - p[1,1] * k_deg_1 - p[1,1] * k_deg_2 - p[1,1] * k_reg
   + 2 p[1,2] * k_deg_2 + 2 p[2,1] * k_deg_1 + 4 p[2,2] * k_reg ] = 0
 Equation_1,_2: [ p[0,2] * k_tx_1 + p[1,1] * k_tx_2 - p[1,2] * k_tx_1
   - p[1,2] * k_tx_2 - p[1,2] * k_deg_1 - 2 p[1,2] * k_deg_2
   - 2 p[1,2] * k_reg + 3 p[1,3] * k_deg_2 + 2 p[2,2] * k_deg_1
   + 6 p[2,3] * k_reg ] = 0
 Equation_1,_3: [ p[0,3] * k_tx_1 + p[1,2] * k_tx_2 - p[1,3] * k_tx_1
   - p[1,3] * k_tx_2 - p[1,3] * k_deg_1 - 3 p[1,3] * k_deg_2
   - 3 p[1,3] * k_reg + 4 p[1,4] * k_deg_2 + 2 p[2,3] * k_deg_1
   + 8 p[2,4] * k_reg ] = 0
 Equation_1,_4: [ p[0,4] * k_tx_1 + p[1,3] * k_tx_2 - p[1,4] * k_tx_1
   - p[1,4] * k_tx_2 - p[1,4] * k_deg_1 - 4 p[1,4] * k_deg_2
   - 4 p[1,4] * k_reg + 5 p[1,5] * k_deg_2 + 2 p[2,4] * k_deg_1
   + 10 p[2,5] * k_reg ] = 0
 Equation_1,_5: [ p[0,5] * k_tx_1 + p[1,4] * k_tx_2 - p[1,5] * k_tx_1
   - p[1,5] * k_tx_2 - p[1,5] * k_deg_1 - 5 p[1,5] * k_deg_2
   - 5 p[1,5] * k_reg + 6 p[1,6] * k_deg_2 + 2 p[2,5] * k_deg_1
   + 12 p[2,6] * k_reg ] = 0
 Equation_2,_0: [ p[1,0] * k_tx_1 - p[2,0] * k_tx_1 - p[2,0] * k_tx_2
   - 2 p[2,0] * k_deg_1 + p[2,1] * k_deg_2 + 3 p[3,0] * k_deg_1
   + 3 p[3,1] * k_reg ] = 0
 Equation_2,_1: [ p[1,1] * k_tx_1 + p[2,0] * k_tx_2 - p[2,1] * k_tx_1
   - p[2,1] * k_tx_2 - 2 p[2,1] * k_deg_1 - p[2,1] * k_deg_2
   - 2 p[2,1] * k_reg + 2 p[2,2] * k_deg_2 + 3 p[3,1] * k_deg_1
   + 6 p[3,2] * k_reg ] = 0
 Equation_2,_2: [ p[1,2] * k_tx_1 + p[2,1] * k_tx_2 - p[2,2] * k_tx_1
   - p[2,2] * k_tx_2 - 2 p[2,2] * k_deg_1 - 2 p[2,2] * k_deg_2
   - 4 p[2,2] * k_reg + 3 p[2,3] * k_deg_2 + 3 p[3,2] * k_deg_1
   + 9 p[3,3] * k_reg ] = 0
 Equation_2,_3: [ p[1,3] * k_tx_1 + p[2,2] * k_tx_2 - p[2,3] * k_tx_1
   - p[2,3] * k_tx_2 - 2 p[2,3] * k_deg_1 - 3 p[2,3] * k_deg_2
   - 6 p[2,3] * k_reg + 4 p[2,4] * k_deg_2 + 3 p[3,3] * k_deg_1
   + 12 p[3,4] * k_reg ] = 0
 Equation_2,_4: [ p[1,4] * k_tx_1 + p[2,3] * k_tx_2 - p[2,4] * k_tx_1
   - p[2,4] * k_tx_2 - 2 p[2,4] * k_deg_1 - 4 p[2,4] * k_deg_2
   - 8 p[2,4] * k_reg + 5 p[2,5] * k_deg_2 + 3 p[3,4] * k_deg_1
   + 15 p[3,5] * k_reg ] = 0
 Equation_2,_5: [ p[1,5] * k_tx_1 + p[2,4] * k_tx_2 - p[2,5] * k_tx_1
   - p[2,5] * k_tx_2 - 2 p[2,5] * k_deg_1 - 5 p[2,5] * k_deg_2
   - 10 p[2,5] * k_reg + 6 p[2,6] * k_deg_2 + 3 p[3,5] * k_deg_1
   + 18 p[3,6] * k_reg ] = 0
 Equation_3,_0: [ p[2,0] * k_tx_1 - p[3,0] * k_tx_1 - p[3,0] * k_tx_2
   - 3 p[3,0] * k_deg_1 + p[3,1] * k_deg_2 + 4 p[4,0] * k_deg_1
   + 4 p[4,1] * k_reg ] = 0
 Equation_3,_1: [ p[2,1] * k_tx_1 + p[3,0] * k_tx_2 - p[3,1] * k_tx_1
   - p[3,1] * k_tx_2 - 3 p[3,1] * k_deg_1 - p[3,1] * k_deg_2
   - 3 p[3,1] * k_reg + 2 p[3,2] * k_deg_2 + 4 p[4,1] * k_deg_1
   + 8 p[4,2] * k_reg ] = 0
 Equation_3,_2: [ p[2,2] * k_tx_1 + p[3,1] * k_tx_2 - p[3,2] * k_tx_1
   - p[3,2] * k_tx_2 - 3 p[3,2] * k_deg_1 - 2 p[3,2] * k_deg_2
   - 6 p[3,2] * k_reg + 3 p[3,3] * k_deg_2 + 4 p[4,2] * k_deg_1
   + 12 p[4,3] * k_reg ] = 0
 Equation_3,_3: [ p[2,3] * k_tx_1 + p[3,2] * k_tx_2 - p[3,3] * k_tx_1
   - p[3,3] * k_tx_2 - 3 p[3,3] * k_deg_1 - 3 p[3,3] * k_deg_2
   - 9 p[3,3] * k_reg + 4 p[3,4] * k_deg_2 + 4 p[4,3] * k_deg_1
   + 16 p[4,4] * k_reg ] = 0
 Equation_3,_4: [ p[2,4] * k_tx_1 + p[3,3] * k_tx_2 - p[3,4] * k_tx_1
   - p[3,4] * k_tx_2 - 3 p[3,4] * k_deg_1 - 4 p[3,4] * k_deg_2
   - 12 p[3,4] * k_reg + 5 p[3,5] * k_deg_2 + 4 p[4,4] * k_deg_1
   + 20 p[4,5] * k_reg ] = 0
 Equation_3,_5: [ p[2,5] * k_tx_1 + p[3,4] * k_tx_2 - p[3,5] * k_tx_1
   - p[3,5] * k_tx_2 - 3 p[3,5] * k_deg_1 - 5 p[3,5] * k_deg_2
   - 15 p[3,5] * k_reg + 6 p[3,6] * k_deg_2 + 4 p[4,5] * k_deg_1
   + 24 p[4,6] * k_reg ] = 0
 Equation_4,_0: [ p[3,0] * k_tx_1 - p[4,0] * k_tx_1 - p[4,0] * k_tx_2
   - 4 p[4,0] * k_deg_1 + p[4,1] * k_deg_2 + 5 p[5,0] * k_deg_1
   + 5 p[5,1] * k_reg ] = 0
 Equation_4,_1: [ p[3,1] * k_tx_1 + p[4,0] * k_tx_2 - p[4,1] * k_tx_1
   - p[4,1] * k_tx_2 - 4 p[4,1] * k_deg_1 - p[4,1] * k_deg_2
   - 4 p[4,1] * k_reg + 2 p[4,2] * k_deg_2 + 5 p[5,1] * k_deg_1
   + 10 p[5,2] * k_reg ] = 0
 Equation_4,_2: [ p[3,2] * k_tx_1 + p[4,1] * k_tx_2 - p[4,2] * k_tx_1
   - p[4,2] * k_tx_2 - 4 p[4,2] * k_deg_1 - 2 p[4,2] * k_deg_2
   - 8 p[4,2] * k_reg + 3 p[4,3] * k_deg_2 + 5 p[5,2] * k_deg_1
   + 15 p[5,3] * k_reg ] = 0
 Equation_4,_3: [ p[3,3] * k_tx_1 + p[4,2] * k_tx_2 - p[4,3] * k_tx_1
   - p[4,3] * k_tx_2 - 4 p[4,3] * k_deg_1 - 3 p[4,3] * k_deg_2
   - 12 p[4,3] * k_reg + 4 p[4,4] * k_deg_2 + 5 p[5,3] * k_deg_1
   + 20 p[5,4] * k_reg ] = 0
 Equation_4,_4: [ p[3,4] * k_tx_1 + p[4,3] * k_tx_2 - p[4,4] * k_tx_1
   - p[4,4] * k_tx_2 - 4 p[4,4] * k_deg_1 - 4 p[4,4] * k_deg_2
   - 16 p[4,4] * k_reg + 5 p[4,5] * k_deg_2 + 5 p[5,4] * k_deg_1
   + 25 p[5,5] * k_reg ] = 0
 Equation_4,_5: [ p[3,5] * k_tx_1 + p[4,4] * k_tx_2 - p[4,5] * k_tx_1
   - p[4,5] * k_tx_2 - 4 p[4,5] * k_deg_1 - 5 p[4,5] * k_deg_2
   - 20 p[4,5] * k_reg + 6 p[4,6] * k_deg_2 + 5 p[5,5] * k_deg_1
   + 30 p[5,6] * k_reg ] = 0
Bounds
 p[0,0] <= 1
 p[0,1] <= 1
 p[0,2] <= 1
 p[0,3] <= 1
 p[0,4] <= 1
 p[0,5] <= 1
 p[0,6] <= 1
 p[1,0] <= 1
 p[1,1] <= 1
 p[1,2] <= 1
 p[1,3] <= 1
 p[1,4] <= 1
 p[1,5] <= 1
 p[1,6] <= 1
 p[2,0] <= 1
 p[2,1] <= 1
 p[2,2] <= 1
 p[2,3] <= 1
 p[2,4] <= 1
 p[2,5] <= 1
 p[2,6] <= 1
 p[3,0] <= 1
 p[3,1] <= 1
 p[3,2] <= 1
 p[3,3] <= 1
 p[3,4] <= 1
 p[3,5] <= 1
 p[3,6] <= 1
 p[4,0] <= 1
 p[4,1] <= 1
 p[4,2] <= 1
 p[4,3] <= 1
 p[4,4] <= 1
 p[4,5] <= 1
 p[4,6] <= 1
 p[5,0] <= 1
 p[5,1] <= 1
 p[5,2] <= 1
 p[5,3] <= 1
 p[5,4] <= 1
 p[5,5] <= 1
 p[5,6] <= 1
 p1[0] <= 1
 p1[1] <= 1
 p1[2] <= 1
 p1[3] <= 1
 p1[4] <= 1
 p1[5] <= 1
 p2[0] <= 1
 p2[1] <= 1
 p2[2] <= 1
 p2[3] <= 1
 p2[4] <= 1
 p2[5] <= 1
 p2[6] <= 1
 k_tx_1 <= 10
 k_tx_2 <= 10
 k_deg_1 <= 10
 k_deg_2 <= 10
 k_reg <= 10
End
