\ Model birth-death-regulation-capture-efficiency
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  k_reg
Subject To
 R0: k_deg_2 = 1
 B_lb_0,_0: p[0,0] + 0.4 p[0,1] + 0.16 p[0,2] + 0.064 p[0,3]
   + 0.0256 p[0,4] + 0.4 p[1,0] + 0.16 p[1,1] + 0.064 p[1,2]
   + 0.0256 p[1,3] + 0.01024 p[1,4] + 0.16 p[2,0] + 0.064 p[2,1]
   + 0.0256 p[2,2] + 0.01024 p[2,3] + 0.004096 p[2,4] + 0.064 p[3,0]
   + 0.0256 p[3,1] + 0.01024 p[3,2] + 0.004096 p[3,3] + 0.0016384 p[3,4]
   + 0.0256 p[4,0] + 0.01024 p[4,1] + 0.004096 p[4,2] + 0.0016384 p[4,3]
   + 0.00065536 p[4,4] >= 0.486
 B_ub_0,_0: p[0,0] + 0.4 p[0,1] + 0.16 p[0,2] + 0.064 p[0,3]
   + 0.0256 p[0,4] + 0.4 p[1,0] + 0.16 p[1,1] + 0.064 p[1,2]
   + 0.0256 p[1,3] + 0.01024 p[1,4] + 0.16 p[2,0] + 0.064 p[2,1]
   + 0.0256 p[2,2] + 0.01024 p[2,3] + 0.004096 p[2,4] + 0.064 p[3,0]
   + 0.0256 p[3,1] + 0.01024 p[3,2] + 0.004096 p[3,3] + 0.0016384 p[3,4]
   + 0.0256 p[4,0] + 0.01024 p[4,1] + 0.004096 p[4,2] + 0.0016384 p[4,3]
   + 0.00065536 p[4,4] <= 0.529605
 B_lb_0,_1: 0.6 p[0,1] + 0.48 p[0,2] + 0.288 p[0,3] + 0.1536 p[0,4]
   + 0.24 p[1,1] + 0.192 p[1,2] + 0.1152 p[1,3] + 0.06144 p[1,4]
   + 0.096 p[2,1] + 0.0768 p[2,2] + 0.04608 p[2,3] + 0.024576 p[2,4]
   + 0.0384 p[3,1] + 0.03072 p[3,2] + 0.018432 p[3,3] + 0.0098304 p[3,4]
   + 0.01536 p[4,1] + 0.012288 p[4,2] + 0.0073728 p[4,3]
   + 0.00393216 p[4,4] >= 0.1602
 B_ub_0,_1: 0.6 p[0,1] + 0.48 p[0,2] + 0.288 p[0,3] + 0.1536 p[0,4]
   + 0.24 p[1,1] + 0.192 p[1,2] + 0.1152 p[1,3] + 0.06144 p[1,4]
   + 0.096 p[2,1] + 0.0768 p[2,2] + 0.04608 p[2,3] + 0.024576 p[2,4]
   + 0.0384 p[3,1] + 0.03072 p[3,2] + 0.018432 p[3,3] + 0.0098304 p[3,4]
   + 0.01536 p[4,1] + 0.012288 p[4,2] + 0.0073728 p[4,3]
   + 0.00393216 p[4,4] <= 0.1978
 B_lb_0,_2: 0.36 p[0,2] + 0.432 p[0,3] + 0.3456 p[0,4] + 0.144 p[1,2]
   + 0.1728 p[1,3] + 0.13824 p[1,4] + 0.0576 p[2,2] + 0.06912 p[2,3]
   + 0.055296 p[2,4] + 0.02304 p[3,2] + 0.027648 p[3,3] + 0.0221184 p[3,4]
   + 0.009216 p[4,2] + 0.0110592 p[4,3] + 0.00884736 p[4,4] >= 0.0164
 B_ub_0,_2: 0.36 p[0,2] + 0.432 p[0,3] + 0.3456 p[0,4] + 0.144 p[1,2]
   + 0.1728 p[1,3] + 0.13824 p[1,4] + 0.0576 p[2,2] + 0.06912 p[2,3]
   + 0.055296 p[2,4] + 0.02304 p[3,2] + 0.027648 p[3,3] + 0.0221184 p[3,4]
   + 0.009216 p[4,2] + 0.0110592 p[4,3] + 0.00884736 p[4,4] <= 0.043
 B_lb_1,_0: 0.6 p[1,0] + 0.24 p[1,1] + 0.096 p[1,2] + 0.0384 p[1,3]
   + 0.01536 p[1,4] + 0.48 p[2,0] + 0.192 p[2,1] + 0.0768 p[2,2]
   + 0.03072 p[2,3] + 0.012288 p[2,4] + 0.288 p[3,0] + 0.1152 p[3,1]
   + 0.04608 p[3,2] + 0.018432 p[3,3] + 0.0073728 p[3,4] + 0.1536 p[4,0]
   + 0.06144 p[4,1] + 0.024576 p[4,2] + 0.0098304 p[4,3]
   + 0.00393216 p[4,4] >= 0.1504
 B_ub_1,_0: 0.6 p[1,0] + 0.24 p[1,1] + 0.096 p[1,2] + 0.0384 p[1,3]
   + 0.01536 p[1,4] + 0.48 p[2,0] + 0.192 p[2,1] + 0.0768 p[2,2]
   + 0.03072 p[2,3] + 0.012288 p[2,4] + 0.288 p[3,0] + 0.1152 p[3,1]
   + 0.04608 p[3,2] + 0.018432 p[3,3] + 0.0073728 p[3,4] + 0.1536 p[4,0]
   + 0.06144 p[4,1] + 0.024576 p[4,2] + 0.0098304 p[4,3]
   + 0.00393216 p[4,4] <= 0.1878
 B_lb_1,_1: 0.36 p[1,1] + 0.288 p[1,2] + 0.1728 p[1,3] + 0.09216 p[1,4]
   + 0.288 p[2,1] + 0.2304 p[2,2] + 0.13824 p[2,3] + 0.073728 p[2,4]
   + 0.1728 p[3,1] + 0.13824 p[3,2] + 0.082944 p[3,3] + 0.0442368 p[3,4]
   + 0.09216 p[4,1] + 0.073728 p[4,2] + 0.0442368 p[4,3]
   + 0.0235929600000013 p[4,4] >= 0.0068
 B_ub_1,_1: 0.36 p[1,1] + 0.288 p[1,2] + 0.1728 p[1,3] + 0.09216 p[1,4]
   + 0.288 p[2,1] + 0.2304 p[2,2] + 0.13824 p[2,3] + 0.073728 p[2,4]
   + 0.1728 p[3,1] + 0.13824 p[3,2] + 0.082944 p[3,3] + 0.0442368 p[3,4]
   + 0.09216 p[4,1] + 0.073728 p[4,2] + 0.0442368 p[4,3]
   + 0.0235929600000013 p[4,4] <= 0.0318
 B_lb_1,_2: 0.216 p[1,2] + 0.2592 p[1,3] + 0.20736 p[1,4] + 0.1728 p[2,2]
   + 0.20736 p[2,3] + 0.165888 p[2,4] + 0.10368 p[3,2] + 0.124416 p[3,3]
   + 0.0995328 p[3,4] + 0.055296 p[4,2] + 0.0663552 p[4,3]
   + 0.0530841599999977 p[4,4] >= -0.0146
 B_ub_1,_2: 0.216 p[1,2] + 0.2592 p[1,3] + 0.20736 p[1,4] + 0.1728 p[2,2]
   + 0.20736 p[2,3] + 0.165888 p[2,4] + 0.10368 p[3,2] + 0.124416 p[3,3]
   + 0.0995328 p[3,4] + 0.055296 p[4,2] + 0.0663552 p[4,3]
   + 0.0530841599999977 p[4,4] <= 0.0044
 B_lb_2,_0: 0.36 p[2,0] + 0.144 p[2,1] + 0.0576 p[2,2] + 0.02304 p[2,3]
   + 0.009216 p[2,4] + 0.432 p[3,0] + 0.1728 p[3,1] + 0.06912 p[3,2]
   + 0.027648 p[3,3] + 0.0110592 p[3,4] + 0.3456 p[4,0] + 0.13824 p[4,1]
   + 0.055296 p[4,2] + 0.0221184 p[4,3] + 0.00884736 p[4,4] >= 0.0166
 B_ub_2,_0: 0.36 p[2,0] + 0.144 p[2,1] + 0.0576 p[2,2] + 0.02304 p[2,3]
   + 0.009216 p[2,4] + 0.432 p[3,0] + 0.1728 p[3,1] + 0.06912 p[3,2]
   + 0.027648 p[3,3] + 0.0110592 p[3,4] + 0.3456 p[4,0] + 0.13824 p[4,1]
   + 0.055296 p[4,2] + 0.0221184 p[4,3] + 0.00884736 p[4,4] <= 0.0432
 B_lb_2,_1: 0.216 p[2,1] + 0.1728 p[2,2] + 0.10368 p[2,3] + 0.055296 p[2,4]
   + 0.2592 p[3,1] + 0.20736 p[3,2] + 0.124416 p[3,3] + 0.0663552 p[3,4]
   + 0.20736 p[4,1] + 0.165888 p[4,2] + 0.0995328 p[4,3]
   + 0.0530841599999977 p[4,4] >= -0.0144
 B_ub_2,_1: 0.216 p[2,1] + 0.1728 p[2,2] + 0.10368 p[2,3] + 0.055296 p[2,4]
   + 0.2592 p[3,1] + 0.20736 p[3,2] + 0.124416 p[3,3] + 0.0663552 p[3,4]
   + 0.20736 p[4,1] + 0.165888 p[4,2] + 0.0995328 p[4,3]
   + 0.0530841599999977 p[4,4] <= 0.0046
 B_lb_2,_2: 0.1296 p[2,2] + 0.15552 p[2,3] + 0.124416 p[2,4]
   + 0.15552 p[3,2] + 0.186624 p[3,3] + 0.1492991999999917 p[3,4]
   + 0.124416 p[4,2] + 0.1492991999999917 p[4,3]
   + 0.1194393599999972 p[4,4] >= -0.016
 B_ub_2,_2: 0.1296 p[2,2] + 0.15552 p[2,3] + 0.124416 p[2,4]
   + 0.15552 p[3,2] + 0.186624 p[3,3] + 0.1492991999999917 p[3,4]
   + 0.124416 p[4,2] + 0.1492991999999917 p[4,3]
   + 0.1194393599999972 p[4,4] <= 1
 B_marginal_lb_0: p1[0] + 0.4 p1[1] + 0.16 p1[2] + 0.064 p1[3]
   + 0.0256 p1[4] >= 0.6246
 B_marginal_ub_0: p1[0] + 0.4 p1[1] + 0.16 p1[2] + 0.064 p1[3]
   + 0.0256 p1[4] <= 0.7574
 B_marginal_lb_1: 0.6 p1[1] + 0.48 p1[2] + 0.288 p1[3] + 0.1536 p1[4]
   >= 0.0878
 B_marginal_ub_1: 0.6 p1[1] + 0.48 p1[2] + 0.288 p1[3] + 0.1536 p1[4]
   <= 0.219
 B_marginal_lb_2: 0.36 p1[2] + 0.432 p1[3] + 0.3456 p1[4] >= -0.0728
 B_marginal_ub_2: 0.36 p1[2] + 0.432 p1[3] + 0.3456 p1[4] <= 0.0464
 B_marginal_lb_0: p2[0] + 0.4 p2[1] + 0.16 p2[2] + 0.064 p2[3]
   + 0.0256 p2[4] >= 0.6153333333333334
 B_marginal_ub_0: p2[0] + 0.4 p2[1] + 0.16 p2[2] + 0.064 p2[3]
   + 0.0256 p2[4] <= 0.7488
 B_marginal_lb_1: 0.6 p2[1] + 0.48 p2[2] + 0.288 p2[3] + 0.1536 p2[4]
   >= 0.0967333333333334
 B_marginal_ub_1: 0.6 p2[1] + 0.48 p2[2] + 0.288 p2[3] + 0.1536 p2[4]
   <= 0.2292
 B_marginal_lb_2: 0.36 p2[2] + 0.432 p2[3] + 0.3456 p2[4]
   >= -0.0746716666666667
 B_marginal_ub_2: 0.36 p2[2] + 0.432 p2[3] + 0.3456 p2[4] <= 0.046
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
 Distribution: p[0,0] + p[0,1] + p[0,2] + p[0,3] + p[0,4] + p[1,0] + p[1,1]
   + p[1,2] + p[1,3] + p[1,4] + p[2,0] + p[2,1] + p[2,2] + p[2,3] + p[2,4]
   + p[3,0] + p[3,1] + p[3,2] + p[3,3] + p[3,4] + p[4,0] + p[4,1] + p[4,2]
   + p[4,3] + p[4,4] <= 1
 Distribution_x1: p1[0] + p1[1] + p1[2] + p1[3] + p1[4] <= 1
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
 k_tx_1 <= 10
 k_tx_2 <= 10
 k_deg_1 <= 10
 k_deg_2 <= 10
 k_reg <= 10
End
