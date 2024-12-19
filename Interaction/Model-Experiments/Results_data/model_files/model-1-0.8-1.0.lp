\ Model birth-death-regulation-capture-efficiency
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  k_reg
Subject To
 R0: k_deg_2 = 1
 B_lb_0,_0: p[0,0] >= 0.177895
 B_ub_0,_0: p[0,0] <= 0.2414
 B_lb_0,_1: p[0,1] >= 0.1231
 B_ub_0,_1: p[0,1] <= 0.1846
 B_lb_0,_2: p[0,2] >= 0.0251
 B_ub_0,_2: p[0,2] <= 0.08
 B_lb_0,_3: p[0,3] >= -0.025105
 B_ub_0,_3: p[0,3] <= 0.023
 B_lb_0,_4: p[0,4] >= -0.037105
 B_ub_0,_4: p[0,4] <= 0.0072
 B_lb_0,_5: p[0,5] >= -0.0405
 B_ub_0,_5: p[0,5] <= 1
 B_lb_1,_0: p[1,0] >= 0.1333
 B_ub_1,_0: p[1,0] <= 0.1952
 B_lb_1,_1: p[1,1] >= 0.0567
 B_ub_1,_1: p[1,1] <= 0.1142
 B_lb_1,_2: p[1,2] >= -0.007705
 B_ub_1,_2: p[1,2] <= 0.0434
 B_lb_1,_3: p[1,3] >= -0.0351
 B_ub_1,_3: p[1,3] <= 0.0104
 B_lb_1,_4: p[1,4] >= -0.0397
 B_ub_1,_4: p[1,4] <= 0.0034
 B_lb_1,_5: p[1,5] >= -0.0405
 B_ub_1,_5: p[1,5] <= 1
 B_lb_2,_0: p[2,0] >= 0.0299
 B_ub_2,_0: p[2,0] <= 0.0854
 B_lb_2,_1: p[2,1] >= -0.0109
 B_ub_2,_1: p[2,1] <= 0.0398
 B_lb_2,_2: p[2,2] >= -0.0339
 B_ub_2,_2: p[2,2] <= 0.012
 B_lb_2,_3: p[2,3] >= -0.0405
 B_ub_2,_3: p[2,3] <= 1
 B_lb_2,_4: p[2,4] >= -0.0405
 B_ub_2,_4: p[2,4] <= 1
 B_lb_2,_5: p[2,5] >= -0.0405
 B_ub_2,_5: p[2,5] <= 1
 B_lb_3,_0: p[3,0] >= -0.0243
 B_ub_3,_0: p[3,0] <= 0.024
 B_lb_3,_1: p[3,1] >= -0.0333
 B_ub_3,_1: p[3,1] <= 0.0126
 B_lb_3,_2: p[3,2] >= -0.0405
 B_ub_3,_2: p[3,2] <= 1
 B_lb_3,_3: p[3,3] >= -0.0405
 B_ub_3,_3: p[3,3] <= 1
 B_lb_3,_4: p[3,4] >= -0.0405
 B_ub_3,_4: p[3,4] <= 1
 B_lb_3,_5: p[3,5] >= -0.0405
 B_ub_3,_5: p[3,5] <= 1
 B_lb_4,_0: p[4,0] >= -0.0389
 B_ub_4,_0: p[4,0] <= 0.0046
 B_lb_4,_1: p[4,1] >= -0.0405
 B_ub_4,_1: p[4,1] <= 1
 B_lb_4,_2: p[4,2] >= -0.0405
 B_ub_4,_2: p[4,2] <= 1
 B_lb_4,_3: p[4,3] >= -0.0405
 B_ub_4,_3: p[4,3] <= 1
 B_lb_4,_4: p[4,4] >= -0.0405
 B_ub_4,_4: p[4,4] <= 1
 B_lb_4,_5: p[4,5] >= -0.0405
 B_ub_4,_5: p[4,5] <= 1
 B_marginal_lb_0: p1[0] >= 0.342154
 B_marginal_ub_0: p1[0] <= 0.5152
 B_marginal_lb_1: p1[1] >= 0.180159
 B_marginal_ub_1: p1[1] <= 0.3518
 B_marginal_lb_2: p1[2] >= -0.030641
 B_marginal_ub_2: p1[2] <= 0.133
 B_marginal_lb_3: p1[3] >= -0.118441
 B_marginal_ub_3: p1[3] <= 0.0364
 B_marginal_lb_4: p1[4] >= -0.142441
 B_marginal_ub_4: p1[4] <= 0.0064
 B_marginal_lb_0: p2[0] >= 0.3828
 B_marginal_ub_0: p2[0] <= 0.5294
 B_marginal_lb_1: p2[1] >= 0.1938
 B_marginal_ub_1: p2[1] <= 0.3388
 B_marginal_lb_2: p2[2] >= -0.007
 B_marginal_ub_2: p2[2] <= 0.1308
 B_marginal_lb_3: p2[3] >= -0.0946
 B_marginal_ub_3: p2[3] <= 0.0336
 B_marginal_lb_4: p2[4] >= -0.1138
 B_marginal_ub_4: p2[4] <= 0.01
 B_marginal_lb_5: p2[5] >= -0.119
 B_marginal_ub_5: p2[5] <= 1
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
 Distribution: p[0,0] + p[0,1] + p[0,2] + p[0,3] + p[0,4] + p[0,5] + p[1,0]
   + p[1,1] + p[1,2] + p[1,3] + p[1,4] + p[1,5] + p[2,0] + p[2,1] + p[2,2]
   + p[2,3] + p[2,4] + p[2,5] + p[3,0] + p[3,1] + p[3,2] + p[3,3] + p[3,4]
   + p[3,5] + p[4,0] + p[4,1] + p[4,2] + p[4,3] + p[4,4] + p[4,5] <= 1
 Distribution_x1: p1[0] + p1[1] + p1[2] + p1[3] + p1[4] <= 1
 Distribution_x2: p2[0] + p2[1] + p2[2] + p2[3] + p2[4] + p2[5] <= 1
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
Bounds
 p[0,0] <= 1
 p[0,1] <= 1
 p[0,2] <= 1
 p[0,3] <= 1
 p[0,4] <= 1
 p[0,5] <= 1
 p[1,0] <= 1
 p[1,1] <= 1
 p[1,2] <= 1
 p[1,3] <= 1
 p[1,4] <= 1
 p[1,5] <= 1
 p[2,0] <= 1
 p[2,1] <= 1
 p[2,2] <= 1
 p[2,3] <= 1
 p[2,4] <= 1
 p[2,5] <= 1
 p[3,0] <= 1
 p[3,1] <= 1
 p[3,2] <= 1
 p[3,3] <= 1
 p[3,4] <= 1
 p[3,5] <= 1
 p[4,0] <= 1
 p[4,1] <= 1
 p[4,2] <= 1
 p[4,3] <= 1
 p[4,4] <= 1
 p[4,5] <= 1
 p1[0] <= 1
 p1[1] <= 1
 p1[2] <= 1
 p1[3] <= 1
 p1[4] <= 1
 p2[0] <= 1
 p2[1] <= 1
 p2[2] <= 1
 p2[3] <= 1
 p2[4] <= 1
 p2[5] <= 1
 k_tx_1 <= 10
 k_tx_2 <= 10
 k_deg_1 <= 10
 k_deg_2 <= 10
 k_reg <= 10
End
