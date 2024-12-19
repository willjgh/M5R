\ Model birth-death-regulation-capture-efficiency
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  k_reg
Subject To
 R0: k_deg_2 = 1
 B_lb_0,_0: p[0,0] >= 0.18806
 B_ub_0,_0: p[0,0] <= 0.2462
 B_lb_0,_1: p[0,1] >= 0.14746
 B_ub_0,_1: p[0,1] <= 0.204
 B_lb_0,_2: p[0,2] >= 0.03186
 B_ub_0,_2: p[0,2] <= 0.0804
 B_lb_0,_3: p[0,3] >= -0.01974
 B_ub_0,_3: p[0,3] <= 0.0222
 B_lb_0,_4: p[0,4] >= -0.02994
 B_ub_0,_4: p[0,4] <= 0.0086
 B_lb_1,_0: p[1,0] >= 0.13706
 B_ub_1,_0: p[1,0] <= 0.193
 B_lb_1,_1: p[1,1] >= 0.05906
 B_ub_1,_1: p[1,1] <= 0.1098
 B_lb_1,_2: p[1,2] >= -0.00514
 B_ub_1,_2: p[1,2] <= 0.0392
 B_lb_1,_3: p[1,3] >= -0.02994
 B_ub_1,_3: p[1,3] <= 0.009
 B_lb_1,_4: p[1,4] >= -0.03434
 B_ub_1,_4: p[1,4] <= 1
 B_lb_2,_0: p[2,0] >= 0.03286
 B_ub_2,_0: p[2,0] <= 0.0818
 B_lb_2,_1: p[2,1] >= -0.01054
 B_ub_2,_1: p[2,1] <= 0.0332
 B_lb_2,_2: p[2,2] >= -0.02994
 B_ub_2,_2: p[2,2] <= 0.0088
 B_lb_2,_3: p[2,3] >= -0.03434
 B_ub_2,_3: p[2,3] <= 1
 B_lb_2,_4: p[2,4] >= -0.03434
 B_ub_2,_4: p[2,4] <= 1
 B_lb_3,_0: p[3,0] >= -0.01574
 B_ub_3,_0: p[3,0] <= 0.027
 B_lb_3,_1: p[3,1] >= -0.02934
 B_ub_3,_1: p[3,1] <= 0.0096
 B_lb_3,_2: p[3,2] >= -0.03434
 B_ub_3,_2: p[3,2] <= 1
 B_lb_3,_3: p[3,3] >= -0.03434
 B_ub_3,_3: p[3,3] <= 1
 B_lb_3,_4: p[3,4] >= -0.03434
 B_ub_3,_4: p[3,4] <= 1
 B_lb_4,_0: p[4,0] >= -0.03094
 B_ub_4,_0: p[4,0] <= 0.007205
 B_lb_4,_1: p[4,1] >= -0.03434
 B_ub_4,_1: p[4,1] <= 1
 B_lb_4,_2: p[4,2] >= -0.03434
 B_ub_4,_2: p[4,2] <= 1
 B_lb_4,_3: p[4,3] >= -0.03434
 B_ub_4,_3: p[4,3] <= 1
 B_lb_4,_4: p[4,4] >= -0.03434
 B_ub_4,_4: p[4,4] <= 1
 B_lb_5,_0: p[5,0] >= -0.03434
 B_ub_5,_0: p[5,0] <= 1
 B_lb_5,_1: p[5,1] >= -0.03434
 B_ub_5,_1: p[5,1] <= 1
 B_lb_5,_2: p[5,2] >= -0.03434
 B_ub_5,_2: p[5,2] <= 1
 B_lb_5,_3: p[5,3] >= -0.03434
 B_ub_5,_3: p[5,3] <= 1
 B_lb_5,_4: p[5,4] >= -0.03434
 B_ub_5,_4: p[5,4] <= 1
 B_marginal_lb_0: p1[0] >= 0.3974616666666666
 B_marginal_ub_0: p1[0] <= 0.540005
 B_marginal_lb_1: p1[1] >= 0.1974616666666666
 B_marginal_ub_1: p1[1] <= 0.3386
 B_marginal_lb_2: p1[2] >= -0.0129333333333333
 B_marginal_ub_2: p1[2] <= 0.1196
 B_marginal_lb_3: p1[3] >= -0.0887333333333333
 B_marginal_ub_3: p1[3] <= 0.036
 B_marginal_lb_4: p1[4] >= -0.1111333333333333
 B_marginal_ub_4: p1[4] <= 0.0082
 B_marginal_lb_5: p1[5] >= -0.1151333333333333
 B_marginal_ub_5: p1[5] <= 1
 B_marginal_lb_0: p2[0] >= 0.365039
 B_marginal_ub_0: p2[0] <= 0.5328
 B_marginal_lb_1: p2[1] >= 0.177639
 B_marginal_ub_1: p2[1] <= 0.3436
 B_marginal_lb_2: p2[2] >= -0.033761
 B_marginal_ub_2: p2[2] <= 0.1236
 B_marginal_lb_3: p2[3] >= -0.117961
 B_marginal_ub_3: p2[3] <= 0.0308
 B_marginal_lb_4: p2[4] >= -0.134561
 B_marginal_ub_4: p2[4] <= 0.01
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
