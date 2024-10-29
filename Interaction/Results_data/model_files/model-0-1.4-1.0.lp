\ Model birth-death-regulation-capture-efficiency
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  k_reg
Subject To
 R0: k_deg_2 = 1
 B_lb_0,_0: p[0,0] >= 0.2207
 B_ub_0,_0: p[0,0] <= 0.2622
 B_lb_0,_1: p[0,1] >= 0.1535
 B_ub_0,_1: p[0,1] <= 0.1932
 B_lb_0,_2: p[0,2] >= 0.0539
 B_ub_0,_2: p[0,2] <= 0.0866
 B_lb_0,_3: p[0,3] >= 0.0025
 B_ub_0,_3: p[0,3] <= 0.029
 B_lb_0,_4: p[0,4] >= -0.0155
 B_ub_0,_4: p[0,4] <= 0.0058
 B_lb_0,_5: p[0,5] >= -0.0179
 B_ub_0,_5: p[0,5] <= 1
 B_lb_0,_6: p[0,6] >= -0.0179
 B_ub_0,_6: p[0,6] <= 1
 B_lb_0,_7: p[0,7] >= -0.0179
 B_ub_0,_7: p[0,7] <= 1
 B_lb_1,_0: p[1,0] >= 0.1603
 B_ub_1,_0: p[1,0] <= 0.1992
 B_lb_1,_1: p[1,1] >= 0.0733
 B_ub_1,_1: p[1,1] <= 0.1076
 B_lb_1,_2: p[1,2] >= 0.0007
 B_ub_1,_2: p[1,2] <= 0.0268
 B_lb_1,_3: p[1,3] >= -0.0143
 B_ub_1,_3: p[1,3] <= 0.0076
 B_lb_1,_4: p[1,4] >= -0.0179
 B_ub_1,_4: p[1,4] <= 1
 B_lb_1,_5: p[1,5] >= -0.0179
 B_ub_1,_5: p[1,5] <= 1
 B_lb_1,_6: p[1,6] >= -0.0179
 B_ub_1,_6: p[1,6] <= 1
 B_lb_1,_7: p[1,7] >= -0.0179
 B_ub_1,_7: p[1,7] <= 1
 B_lb_2,_0: p[2,0] >= 0.0523
 B_ub_2,_0: p[2,0] <= 0.0852
 B_lb_2,_1: p[2,1] >= 0.0049
 B_ub_2,_1: p[2,1] <= 0.0318
 B_lb_2,_2: p[2,2] >= -0.0145
 B_ub_2,_2: p[2,2] <= 0.0072
 B_lb_2,_3: p[2,3] >= -0.0179
 B_ub_2,_3: p[2,3] <= 1
 B_lb_2,_4: p[2,4] >= -0.0179
 B_ub_2,_4: p[2,4] <= 1
 B_lb_2,_5: p[2,5] >= -0.0179
 B_ub_2,_5: p[2,5] <= 1
 B_lb_2,_6: p[2,6] >= -0.0179
 B_ub_2,_6: p[2,6] <= 1
 B_lb_2,_7: p[2,7] >= -0.0179
 B_ub_2,_7: p[2,7] <= 1
 B_lb_3,_0: p[3,0] >= -0.0021
 B_ub_3,_0: p[3,0] <= 0.0236
 B_lb_3,_1: p[3,1] >= -0.0145
 B_ub_3,_1: p[3,1] <= 0.0074
 B_lb_3,_2: p[3,2] >= -0.0179
 B_ub_3,_2: p[3,2] <= 1
 B_lb_3,_3: p[3,3] >= -0.0179
 B_ub_3,_3: p[3,3] <= 1
 B_lb_3,_4: p[3,4] >= -0.0179
 B_ub_3,_4: p[3,4] <= 1
 B_lb_3,_5: p[3,5] >= -0.0179
 B_ub_3,_5: p[3,5] <= 1
 B_lb_3,_6: p[3,6] >= -0.0179
 B_ub_3,_6: p[3,6] <= 1
 B_lb_3,_7: p[3,7] >= -0.0179
 B_ub_3,_7: p[3,7] <= 1
 B_lb_4,_0: p[4,0] >= -0.0153
 B_ub_4,_0: p[4,0] <= 0.006
 B_lb_4,_1: p[4,1] >= -0.0179
 B_ub_4,_1: p[4,1] <= 1
 B_lb_4,_2: p[4,2] >= -0.0179
 B_ub_4,_2: p[4,2] <= 1
 B_lb_4,_3: p[4,3] >= -0.0179
 B_ub_4,_3: p[4,3] <= 1
 B_lb_4,_4: p[4,4] >= -0.0179
 B_ub_4,_4: p[4,4] <= 1
 B_lb_4,_5: p[4,5] >= -0.0179
 B_ub_4,_5: p[4,5] <= 1
 B_lb_4,_6: p[4,6] >= -0.0179
 B_ub_4,_6: p[4,6] <= 1
 B_lb_4,_7: p[4,7] >= -0.0179
 B_ub_4,_7: p[4,7] <= 1
 B_lb_5,_0: p[5,0] >= -0.0179
 B_ub_5,_0: p[5,0] <= 1
 B_lb_5,_1: p[5,1] >= -0.0179
 B_ub_5,_1: p[5,1] <= 1
 B_lb_5,_2: p[5,2] >= -0.0179
 B_ub_5,_2: p[5,2] <= 1
 B_lb_5,_3: p[5,3] >= -0.0179
 B_ub_5,_3: p[5,3] <= 1
 B_lb_5,_4: p[5,4] >= -0.0179
 B_ub_5,_4: p[5,4] <= 1
 B_lb_5,_5: p[5,5] >= -0.0179
 B_ub_5,_5: p[5,5] <= 1
 B_lb_5,_6: p[5,6] >= -0.0179
 B_ub_5,_6: p[5,6] <= 1
 B_lb_5,_7: p[5,7] >= -0.0179
 B_ub_5,_7: p[5,7] <= 1
 B_marginal_lb_0: p1[0] >= 0.4162
 B_marginal_ub_0: p1[0] <= 0.5544
 B_marginal_lb_1: p1[1] >= 0.1926
 B_marginal_ub_1: p1[1] <= 0.3292
 B_marginal_lb_2: p1[2] >= -0.009
 B_marginal_ub_2: p1[2] <= 0.118805
 B_marginal_lb_3: p1[3] >= -0.0894
 B_marginal_ub_3: p1[3] <= 0.0302
 B_marginal_lb_4: p1[4] >= -0.1074
 B_marginal_ub_4: p1[4] <= 0.0074
 B_marginal_lb_5: p1[5] >= -0.1108
 B_marginal_ub_5: p1[5] <= 1
 B_marginal_lb_0: p2[0] >= 0.4438693750000001
 B_marginal_ub_0: p2[0] <= 0.555005
 B_marginal_lb_1: p2[1] >= 0.218074375
 B_marginal_ub_1: p2[1] <= 0.327605
 B_marginal_lb_2: p2[2] >= 0.015474375
 B_marginal_ub_2: p2[2] <= 0.1166
 B_marginal_lb_3: p2[3] >= -0.058125625
 B_marginal_ub_3: p2[3] <= 0.035
 B_marginal_lb_4: p2[4] >= -0.080525625
 B_marginal_ub_4: p2[4] <= 0.007
 B_marginal_lb_5: p2[5] >= -0.083725625
 B_marginal_ub_5: p2[5] <= 1
 B_marginal_lb_6: p2[6] >= -0.083725625
 B_marginal_ub_6: p2[6] <= 1
 B_marginal_lb_7: p2[7] >= -0.083725625
 B_marginal_ub_7: p2[7] <= 1
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
 frechet_lb_0,_7: - p[0,7] + p1[0] + p2[7] <= 1
 frechet_ub_x1_0: p[0,7] - p1[0] <= 0
 frechet_ub_x2_7: p[0,7] - p2[7] <= 0
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
 frechet_lb_1,_7: - p[1,7] + p1[1] + p2[7] <= 1
 frechet_ub_x1_1: p[1,7] - p1[1] <= 0
 frechet_ub_x2_7: p[1,7] - p2[7] <= 0
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
 frechet_lb_2,_7: - p[2,7] + p1[2] + p2[7] <= 1
 frechet_ub_x1_2: p[2,7] - p1[2] <= 0
 frechet_ub_x2_7: p[2,7] - p2[7] <= 0
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
 frechet_lb_3,_7: - p[3,7] + p1[3] + p2[7] <= 1
 frechet_ub_x1_3: p[3,7] - p1[3] <= 0
 frechet_ub_x2_7: p[3,7] - p2[7] <= 0
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
 frechet_lb_4,_7: - p[4,7] + p1[4] + p2[7] <= 1
 frechet_ub_x1_4: p[4,7] - p1[4] <= 0
 frechet_ub_x2_7: p[4,7] - p2[7] <= 0
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
 frechet_lb_5,_7: - p[5,7] + p1[5] + p2[7] <= 1
 frechet_ub_x1_5: p[5,7] - p1[5] <= 0
 frechet_ub_x2_7: p[5,7] - p2[7] <= 0
 Distribution: p[0,0] + p[0,1] + p[0,2] + p[0,3] + p[0,4] + p[0,5] + p[0,6]
   + p[0,7] + p[1,0] + p[1,1] + p[1,2] + p[1,3] + p[1,4] + p[1,5] + p[1,6]
   + p[1,7] + p[2,0] + p[2,1] + p[2,2] + p[2,3] + p[2,4] + p[2,5] + p[2,6]
   + p[2,7] + p[3,0] + p[3,1] + p[3,2] + p[3,3] + p[3,4] + p[3,5] + p[3,6]
   + p[3,7] + p[4,0] + p[4,1] + p[4,2] + p[4,3] + p[4,4] + p[4,5] + p[4,6]
   + p[4,7] + p[5,0] + p[5,1] + p[5,2] + p[5,3] + p[5,4] + p[5,5] + p[5,6]
   + p[5,7] <= 1
 Distribution_x1: p1[0] + p1[1] + p1[2] + p1[3] + p1[4] + p1[5] <= 1
 Distribution_x2: p2[0] + p2[1] + p2[2] + p2[3] + p2[4] + p2[5] + p2[6]
   + p2[7] <= 1
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
 Equation_0,_6: [ p[0,5] * k_tx_2 - p[0,6] * k_tx_1 - p[0,6] * k_tx_2
   - 6 p[0,6] * k_deg_2 + 7 p[0,7] * k_deg_2 + p[1,6] * k_deg_1
   + 7 p[1,7] * k_reg ] = 0
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
 Equation_1,_6: [ p[0,6] * k_tx_1 + p[1,5] * k_tx_2 - p[1,6] * k_tx_1
   - p[1,6] * k_tx_2 - p[1,6] * k_deg_1 - 6 p[1,6] * k_deg_2
   - 6 p[1,6] * k_reg + 7 p[1,7] * k_deg_2 + 2 p[2,6] * k_deg_1
   + 14 p[2,7] * k_reg ] = 0
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
 Equation_2,_6: [ p[1,6] * k_tx_1 + p[2,5] * k_tx_2 - p[2,6] * k_tx_1
   - p[2,6] * k_tx_2 - 2 p[2,6] * k_deg_1 - 6 p[2,6] * k_deg_2
   - 12 p[2,6] * k_reg + 7 p[2,7] * k_deg_2 + 3 p[3,6] * k_deg_1
   + 21 p[3,7] * k_reg ] = 0
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
 Equation_3,_6: [ p[2,6] * k_tx_1 + p[3,5] * k_tx_2 - p[3,6] * k_tx_1
   - p[3,6] * k_tx_2 - 3 p[3,6] * k_deg_1 - 6 p[3,6] * k_deg_2
   - 18 p[3,6] * k_reg + 7 p[3,7] * k_deg_2 + 4 p[4,6] * k_deg_1
   + 28 p[4,7] * k_reg ] = 0
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
 Equation_4,_6: [ p[3,6] * k_tx_1 + p[4,5] * k_tx_2 - p[4,6] * k_tx_1
   - p[4,6] * k_tx_2 - 4 p[4,6] * k_deg_1 - 6 p[4,6] * k_deg_2
   - 24 p[4,6] * k_reg + 7 p[4,7] * k_deg_2 + 5 p[5,6] * k_deg_1
   + 35 p[5,7] * k_reg ] = 0
Bounds
 p[0,0] <= 1
 p[0,1] <= 1
 p[0,2] <= 1
 p[0,3] <= 1
 p[0,4] <= 1
 p[0,5] <= 1
 p[0,6] <= 1
 p[0,7] <= 1
 p[1,0] <= 1
 p[1,1] <= 1
 p[1,2] <= 1
 p[1,3] <= 1
 p[1,4] <= 1
 p[1,5] <= 1
 p[1,6] <= 1
 p[1,7] <= 1
 p[2,0] <= 1
 p[2,1] <= 1
 p[2,2] <= 1
 p[2,3] <= 1
 p[2,4] <= 1
 p[2,5] <= 1
 p[2,6] <= 1
 p[2,7] <= 1
 p[3,0] <= 1
 p[3,1] <= 1
 p[3,2] <= 1
 p[3,3] <= 1
 p[3,4] <= 1
 p[3,5] <= 1
 p[3,6] <= 1
 p[3,7] <= 1
 p[4,0] <= 1
 p[4,1] <= 1
 p[4,2] <= 1
 p[4,3] <= 1
 p[4,4] <= 1
 p[4,5] <= 1
 p[4,6] <= 1
 p[4,7] <= 1
 p[5,0] <= 1
 p[5,1] <= 1
 p[5,2] <= 1
 p[5,3] <= 1
 p[5,4] <= 1
 p[5,5] <= 1
 p[5,6] <= 1
 p[5,7] <= 1
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
 p2[7] <= 1
 k_tx_1 <= 10
 k_tx_2 <= 10
 k_deg_1 <= 10
 k_deg_2 <= 10
 k_reg <= 10
End
