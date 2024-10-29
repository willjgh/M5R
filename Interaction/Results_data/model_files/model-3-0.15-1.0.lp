\ Model birth-death-regulation-capture-efficiency
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  k_reg
Subject To
 R0: k_deg_2 = 1
 B_lb_0,_0: p[0,0] >= 0.0725477142857143
 B_ub_0,_0: p[0,0] <= 0.1614
 B_lb_0,_1: p[0,1] >= 0.0737477142857143
 B_ub_0,_1: p[0,1] <= 0.1624
 B_lb_0,_2: p[0,2] >= -0.0076522857142857
 B_ub_0,_2: p[0,2] <= 0.075605
 B_lb_0,_3: p[0,3] >= -0.0500522857142857
 B_ub_0,_3: p[0,3] <= 0.0272
 B_lb_0,_4: p[0,4] >= -0.0660522857142857
 B_ub_0,_4: p[0,4] <= 0.0068
 B_lb_0,_5: p[0,5] >= -0.0690522857142857
 B_ub_0,_5: p[0,5] <= 1
 B_lb_0,_6: p[0,6] >= -0.0690522857142857
 B_ub_0,_6: p[0,6] <= 1
 B_lb_1,_0: p[1,0] >= 0.0771427142857143
 B_ub_1,_0: p[1,0] <= 0.166205
 B_lb_1,_1: p[1,1] >= 0.0579477142857143
 B_ub_1,_1: p[1,1] <= 0.1458
 B_lb_1,_2: p[1,2] >= -0.0230522857142857
 B_ub_1,_2: p[1,2] <= 0.0584
 B_lb_1,_3: p[1,3] >= -0.0548522857142857
 B_ub_1,_3: p[1,3] <= 0.0218
 B_lb_1,_4: p[1,4] >= -0.0658522857142857
 B_ub_1,_4: p[1,4] <= 0.007
 B_lb_1,_5: p[1,5] >= -0.0690522857142857
 B_ub_1,_5: p[1,5] <= 1
 B_lb_1,_6: p[1,6] >= -0.0690522857142857
 B_ub_1,_6: p[1,6] <= 1
 B_lb_2,_0: p[2,0] >= -0.00685228571428571
 B_ub_2,_0: p[2,0] <= 0.0762
 B_lb_2,_1: p[2,1] >= -0.0102522857142857
 B_ub_2,_1: p[2,1] <= 0.0728
 B_lb_2,_2: p[2,2] >= -0.0474522857142857
 B_ub_2,_2: p[2,2] <= 0.0304
 B_lb_2,_3: p[2,3] >= -0.0638522857142857
 B_ub_2,_3: p[2,3] <= 0.01
 B_lb_2,_4: p[2,4] >= -0.0690522857142857
 B_ub_2,_4: p[2,4] <= 1
 B_lb_2,_5: p[2,5] >= -0.0690522857142857
 B_ub_2,_5: p[2,5] <= 1
 B_lb_2,_6: p[2,6] >= -0.0690522857142857
 B_ub_2,_6: p[2,6] <= 1
 B_lb_3,_0: p[3,0] >= -0.0528522857142857
 B_ub_3,_0: p[3,0] <= 0.0238
 B_lb_3,_1: p[3,1] >= -0.0538522857142857
 B_ub_3,_1: p[3,1] <= 0.0226
 B_lb_3,_2: p[3,2] >= -0.0648522857142857
 B_ub_3,_2: p[3,2] <= 0.0088
 B_lb_3,_3: p[3,3] >= -0.0690522857142857
 B_ub_3,_3: p[3,3] <= 1
 B_lb_3,_4: p[3,4] >= -0.0690522857142857
 B_ub_3,_4: p[3,4] <= 1
 B_lb_3,_5: p[3,5] >= -0.0690522857142857
 B_ub_3,_5: p[3,5] <= 1
 B_lb_3,_6: p[3,6] >= -0.0690522857142857
 B_ub_3,_6: p[3,6] <= 1
 B_lb_4,_0: p[4,0] >= -0.0656522857142857
 B_ub_4,_0: p[4,0] <= 0.0074
 B_lb_4,_1: p[4,1] >= -0.0666522857142857
 B_ub_4,_1: p[4,1] <= 0.0058
 B_lb_4,_2: p[4,2] >= -0.0680522857142857
 B_ub_4,_2: p[4,2] <= 0.0036
 B_lb_4,_3: p[4,3] >= -0.0690522857142857
 B_ub_4,_3: p[4,3] <= 1
 B_lb_4,_4: p[4,4] >= -0.0690522857142857
 B_ub_4,_4: p[4,4] <= 1
 B_lb_4,_5: p[4,5] >= -0.0690522857142857
 B_ub_4,_5: p[4,5] <= 1
 B_lb_4,_6: p[4,6] >= -0.0690522857142857
 B_ub_4,_6: p[4,6] <= 1
 B_marginal_lb_0: p1[0] >= 0.20148
 B_marginal_ub_0: p1[0] <= 0.4146
 B_marginal_lb_1: p1[1] >= 0.16808
 B_marginal_ub_1: p1[1] <= 0.3812
 B_marginal_lb_2: p1[2] >= -0.02632
 B_marginal_ub_2: p1[2] <= 0.1808
 B_marginal_lb_3: p1[3] >= -0.14592
 B_marginal_ub_3: p1[3] <= 0.0522
 B_marginal_lb_4: p1[4] >= -0.17692
 B_marginal_ub_4: p1[4] <= 0.0156
 B_marginal_lb_0: p2[0] >= 0.2579714285714285
 B_marginal_ub_0: p2[0] <= 0.4172
 B_marginal_lb_1: p2[1] >= 0.2327714285714285
 B_marginal_ub_1: p2[1] <= 0.3916
 B_marginal_lb_2: p2[2] >= 0.0127714285714286
 B_marginal_ub_2: p2[2] <= 0.1652
 B_marginal_lb_3: p2[3] >= -0.0874285714285714
 B_marginal_ub_3: p2[3] <= 0.0566
 B_marginal_lb_4: p2[4] >= -0.1232285714285714
 B_marginal_ub_4: p2[4] <= 0.0144
 B_marginal_lb_5: p2[5] >= -0.1318285714285714
 B_marginal_ub_5: p2[5] <= 1
 B_marginal_lb_6: p2[6] >= -0.1318285714285714
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
 Distribution: p[0,0] + p[0,1] + p[0,2] + p[0,3] + p[0,4] + p[0,5] + p[0,6]
   + p[1,0] + p[1,1] + p[1,2] + p[1,3] + p[1,4] + p[1,5] + p[1,6] + p[2,0]
   + p[2,1] + p[2,2] + p[2,3] + p[2,4] + p[2,5] + p[2,6] + p[3,0] + p[3,1]
   + p[3,2] + p[3,3] + p[3,4] + p[3,5] + p[3,6] + p[4,0] + p[4,1] + p[4,2]
   + p[4,3] + p[4,4] + p[4,5] + p[4,6] <= 1
 Distribution_x1: p1[0] + p1[1] + p1[2] + p1[3] + p1[4] <= 1
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
 p2[6] <= 1
 k_tx_1 <= 10
 k_tx_2 <= 10
 k_deg_1 <= 10
 k_deg_2 <= 10
 k_reg <= 10
End
