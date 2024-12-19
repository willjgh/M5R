\ Model birth-death-regulation-capture-efficiency
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  k_reg
Subject To
 R0: k_deg_2 = 1
 B_lb_0,_0: p[0,0] >= 0.21736
 B_ub_0,_0: p[0,0] <= 0.2674
 B_lb_0,_1: p[0,1] >= 0.15976
 B_ub_0,_1: p[0,1] <= 0.2086
 B_lb_0,_2: p[0,2] >= 0.04216
 B_ub_0,_2: p[0,2] <= 0.0832
 B_lb_0,_3: p[0,3] >= -0.01164
 B_ub_0,_3: p[0,3] <= 0.0222
 B_lb_0,_4: p[0,4] >= -0.02384
 B_ub_0,_4: p[0,4] <= 0.006
 B_lb_1,_0: p[1,0] >= 0.15376
 B_ub_1,_0: p[1,0] <= 0.202
 B_lb_1,_1: p[1,1] >= 0.06116
 B_ub_1,_1: p[1,1] <= 0.1042
 B_lb_1,_2: p[1,2] >= -0.00604
 B_ub_1,_2: p[1,2] <= 0.029
 B_lb_1,_3: p[1,3] >= -0.02344
 B_ub_1,_3: p[1,3] <= 0.0068
 B_lb_1,_4: p[1,4] >= -0.02524
 B_ub_1,_4: p[1,4] <= 0.004
 B_lb_2,_0: p[2,0] >= 0.04156
 B_ub_2,_0: p[2,0] <= 0.0822
 B_lb_2,_1: p[2,1] >= -0.01084
 B_ub_2,_1: p[2,1] <= 0.0232
 B_lb_2,_2: p[2,2] >= -0.02484
 B_ub_2,_2: p[2,2] <= 0.0046
 B_lb_2,_3: p[2,3] >= -0.02644
 B_ub_2,_3: p[2,3] <= 1
 B_lb_2,_4: p[2,4] >= -0.02644
 B_ub_2,_4: p[2,4] <= 1
 B_lb_3,_0: p[3,0] >= -0.00964
 B_ub_3,_0: p[3,0] <= 0.0248
 B_lb_3,_1: p[3,1] >= -0.02504
 B_ub_3,_1: p[3,1] <= 0.0044
 B_lb_3,_2: p[3,2] >= -0.02644
 B_ub_3,_2: p[3,2] <= 1
 B_lb_3,_3: p[3,3] >= -0.02644
 B_ub_3,_3: p[3,3] <= 1
 B_lb_3,_4: p[3,4] >= -0.02644
 B_ub_3,_4: p[3,4] <= 1
 B_lb_4,_0: p[4,0] >= -0.02404
 B_ub_4,_0: p[4,0] <= 0.006
 B_lb_4,_1: p[4,1] >= -0.02644
 B_ub_4,_1: p[4,1] <= 1
 B_lb_4,_2: p[4,2] >= -0.02644
 B_ub_4,_2: p[4,2] <= 1
 B_lb_4,_3: p[4,3] >= -0.02644
 B_ub_4,_3: p[4,3] <= 1
 B_lb_4,_4: p[4,4] >= -0.02644
 B_ub_4,_4: p[4,4] <= 1
 B_lb_5,_0: p[5,0] >= -0.02644
 B_ub_5,_0: p[5,0] <= 1
 B_lb_5,_1: p[5,1] >= -0.02644
 B_ub_5,_1: p[5,1] <= 1
 B_lb_5,_2: p[5,2] >= -0.02644
 B_ub_5,_2: p[5,2] <= 1
 B_lb_5,_3: p[5,3] >= -0.02644
 B_ub_5,_3: p[5,3] <= 1
 B_lb_5,_4: p[5,4] >= -0.02644
 B_ub_5,_4: p[5,4] <= 1
 B_marginal_lb_0: p1[0] >= 0.4328333333333333
 B_marginal_ub_0: p1[0] <= 0.566
 B_marginal_lb_1: p1[1] >= 0.2004333333333334
 B_marginal_ub_1: p1[1] <= 0.3324
 B_marginal_lb_2: p1[2] >= -0.0151666666666667
 B_marginal_ub_2: p1[2] <= 0.1068
 B_marginal_lb_3: p1[3] >= -0.0857666666666667
 B_marginal_ub_3: p1[3] <= 0.0286
 B_marginal_lb_4: p1[4] >= -0.1027666666666667
 B_marginal_ub_4: p1[4] <= 0.0068
 B_marginal_lb_5: p1[5] >= -0.1057666666666667
 B_marginal_ub_5: p1[5] <= 1
 B_marginal_lb_0: p2[0] >= 0.402395
 B_marginal_ub_0: p2[0] <= 0.5608
 B_marginal_lb_1: p2[1] >= 0.1722
 B_marginal_ub_1: p2[1] <= 0.3292
 B_marginal_lb_2: p2[2] >= -0.0352
 B_marginal_ub_2: p2[2] <= 0.112405
 B_marginal_lb_3: p2[3] >= -0.1104
 B_marginal_ub_3: p2[3] <= 0.029
 B_marginal_lb_4: p2[4] >= -0.1262
 B_marginal_ub_4: p2[4] <= 0.0092
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
 Distribution: p[0,0] + p[0,1] + p[0,2] + p[0,3] + p[0,4] + p[1,0] + p[1,1]
   + p[1,2] + p[1,3] + p[1,4] + p[2,0] + p[2,1] + p[2,2] + p[2,3] + p[2,4]
   + p[3,0] + p[3,1] + p[3,2] + p[3,3] + p[3,4] + p[4,0] + p[4,1] + p[4,2]
   + p[4,3] + p[4,4] + p[5,0] + p[5,1] + p[5,2] + p[5,3] + p[5,4] <= 1
 Distribution_x1: p1[0] + p1[1] + p1[2] + p1[3] + p1[4] + p1[5] <= 1
 Distribution_x2: p2[0] + p2[1] + p2[2] + p2[3] + p2[4] <= 1
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
Bounds
 p[0,0] <= 1
 p[0,1] <= 1
 p[0,2] <= 1
 p[0,3] <= 1
 p[0,4] <= 1
 p[1,0] <= 1
 p[1,1] <= 1
 p[1,2] <= 1
 p[1,3] <= 1
 p[1,4] <= 1
 p[2,0] <= 1
 p[2,1] <= 1
 p[2,2] <= 1
 p[2,3] <= 1
 p[2,4] <= 1
 p[3,0] <= 1
 p[3,1] <= 1
 p[3,2] <= 1
 p[3,3] <= 1
 p[3,4] <= 1
 p[4,0] <= 1
 p[4,1] <= 1
 p[4,2] <= 1
 p[4,3] <= 1
 p[4,4] <= 1
 p[5,0] <= 1
 p[5,1] <= 1
 p[5,2] <= 1
 p[5,3] <= 1
 p[5,4] <= 1
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
 k_tx_1 <= 10
 k_tx_2 <= 10
 k_deg_1 <= 10
 k_deg_2 <= 10
 k_reg <= 10
End
