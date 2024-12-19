\ Model birth-death-regulation-capture-efficiency
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  k_reg
Subject To
 R0: k_deg_2 = 1
 B_lb_0,_0: p[0,0] >= 0.28965
 B_ub_0,_0: p[0,0] <= 0.325805
 B_lb_0,_1: p[0,1] >= 0.179245
 B_ub_0,_1: p[0,1] <= 0.2126
 B_lb_0,_2: p[0,2] >= 0.05705
 B_ub_0,_2: p[0,2] <= 0.0824
 B_lb_0,_3: p[0,3] >= 0.00605
 B_ub_0,_3: p[0,3] <= 0.0246
 B_lb_0,_4: p[0,4] >= -0.00775
 B_ub_0,_4: p[0,4] <= 0.0068
 B_lb_0,_5: p[0,5] >= -0.01075
 B_ub_0,_5: p[0,5] <= 1
 B_lb_1,_0: p[1,0] >= 0.18025
 B_ub_1,_0: p[1,0] <= 0.2134
 B_lb_1,_1: p[1,1] >= 0.04445
 B_ub_1,_1: p[1,1] <= 0.0686
 B_lb_1,_2: p[1,2] >= -0.00515
 B_ub_1,_2: p[1,2] <= 0.0104
 B_lb_1,_3: p[1,3] >= -0.00975
 B_ub_1,_3: p[1,3] <= 0.0036
 B_lb_1,_4: p[1,4] >= -0.01075
 B_ub_1,_4: p[1,4] <= 1
 B_lb_1,_5: p[1,5] >= -0.01075
 B_ub_1,_5: p[1,5] <= 1
 B_lb_2,_0: p[2,0] >= 0.05605
 B_ub_2,_0: p[2,0] <= 0.0816
 B_lb_2,_1: p[2,1] >= -0.00175
 B_ub_2,_1: p[2,1] <= 0.0148
 B_lb_2,_2: p[2,2] >= -0.01075
 B_ub_2,_2: p[2,2] <= 1
 B_lb_2,_3: p[2,3] >= -0.01075
 B_ub_2,_3: p[2,3] <= 1
 B_lb_2,_4: p[2,4] >= -0.01075
 B_ub_2,_4: p[2,4] <= 1
 B_lb_2,_5: p[2,5] >= -0.01075
 B_ub_2,_5: p[2,5] <= 1
 B_lb_3,_0: p[3,0] >= 0.00305
 B_ub_3,_0: p[3,0] <= 0.021
 B_lb_3,_1: p[3,1] >= -0.01075
 B_ub_3,_1: p[3,1] <= 1
 B_lb_3,_2: p[3,2] >= -0.01075
 B_ub_3,_2: p[3,2] <= 1
 B_lb_3,_3: p[3,3] >= -0.01075
 B_ub_3,_3: p[3,3] <= 1
 B_lb_3,_4: p[3,4] >= -0.01075
 B_ub_3,_4: p[3,4] <= 1
 B_lb_3,_5: p[3,5] >= -0.01075
 B_ub_3,_5: p[3,5] <= 1
 B_lb_4,_0: p[4,0] >= -0.00875
 B_ub_4,_0: p[4,0] <= 0.0054
 B_lb_4,_1: p[4,1] >= -0.01075
 B_ub_4,_1: p[4,1] <= 1
 B_lb_4,_2: p[4,2] >= -0.01075
 B_ub_4,_2: p[4,2] <= 1
 B_lb_4,_3: p[4,3] >= -0.01075
 B_ub_4,_3: p[4,3] <= 1
 B_lb_4,_4: p[4,4] >= -0.01075
 B_ub_4,_4: p[4,4] <= 1
 B_lb_4,_5: p[4,5] >= -0.01075
 B_ub_4,_5: p[4,5] <= 1
 B_lb_5,_0: p[5,0] >= -0.01075
 B_ub_5,_0: p[5,0] <= 1
 B_lb_5,_1: p[5,1] >= -0.01075
 B_ub_5,_1: p[5,1] <= 1
 B_lb_5,_2: p[5,2] >= -0.01075
 B_ub_5,_2: p[5,2] <= 1
 B_lb_5,_3: p[5,3] >= -0.01075
 B_ub_5,_3: p[5,3] <= 1
 B_lb_5,_4: p[5,4] >= -0.01075
 B_ub_5,_4: p[5,4] <= 1
 B_lb_5,_5: p[5,5] >= -0.01075
 B_ub_5,_5: p[5,5] <= 1
 B_marginal_lb_0: p1[0] >= 0.5113666666666666
 B_marginal_ub_0: p1[0] <= 0.629
 B_marginal_lb_1: p1[1] >= 0.1709666666666667
 B_marginal_ub_1: p1[1] <= 0.2868
 B_marginal_lb_2: p1[2] >= -0.0118333333333333
 B_marginal_ub_2: p1[2] <= 0.094605
 B_marginal_lb_3: p1[3] >= -0.0760333333333333
 B_marginal_ub_3: p1[3] <= 0.0222
 B_marginal_lb_4: p1[4] >= -0.0888333333333333
 B_marginal_ub_4: p1[4] <= 0.0054
 B_marginal_lb_5: p1[5] >= -0.0908333333333333
 B_marginal_ub_5: p1[5] <= 1
 B_marginal_lb_0: p2[0] >= 0.5046333333333332
 B_marginal_ub_0: p2[0] <= 0.625005
 B_marginal_lb_1: p2[1] >= 0.1704333333333333
 B_marginal_ub_1: p2[1] <= 0.2888
 B_marginal_lb_2: p2[2] >= -0.0173666666666667
 B_marginal_ub_2: p2[2] <= 0.0912
 B_marginal_lb_3: p2[3] >= -0.0745666666666667
 B_marginal_ub_3: p2[3] <= 0.027
 B_marginal_lb_4: p2[4] >= -0.0899666666666667
 B_marginal_ub_4: p2[4] <= 0.0072
 B_marginal_lb_5: p2[5] >= -0.0933666666666667
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
 Distribution: p[0,0] + p[0,1] + p[0,2] + p[0,3] + p[0,4] + p[0,5] + p[1,0]
   + p[1,1] + p[1,2] + p[1,3] + p[1,4] + p[1,5] + p[2,0] + p[2,1] + p[2,2]
   + p[2,3] + p[2,4] + p[2,5] + p[3,0] + p[3,1] + p[3,2] + p[3,3] + p[3,4]
   + p[3,5] + p[4,0] + p[4,1] + p[4,2] + p[4,3] + p[4,4] + p[4,5] + p[5,0]
   + p[5,1] + p[5,2] + p[5,3] + p[5,4] + p[5,5] <= 1
 Distribution_x1: p1[0] + p1[1] + p1[2] + p1[3] + p1[4] + p1[5] <= 1
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
 p[5,0] <= 1
 p[5,1] <= 1
 p[5,2] <= 1
 p[5,3] <= 1
 p[5,4] <= 1
 p[5,5] <= 1
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
 k_tx_1 <= 10
 k_tx_2 <= 10
 k_deg_1 <= 10
 k_deg_2 <= 10
 k_reg <= 10
End
