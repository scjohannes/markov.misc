# print and summary methods work

    Code
      print(result)
    Output
      Robust (Sandwich) Covariance Estimation for vglm
      ================================================
      
      Call:
      VGAM::vglm(formula = y ~ x, family = VGAM::binomialff, data = test_data)
      
      Number of observations: 100 
      
      Coefficients:
                (Intercept)       x
      Estimate       0.1030 -0.2566
      Robust SE      0.2029  0.2093
      
      (Use summary() for full coefficient table with z-values and p-values)

---

    Code
      summary(result)
    Output
      Robust (Sandwich) Covariance Estimation for vglm
      ================================================
      
      Call:
      VGAM::vglm(formula = y ~ x, family = VGAM::binomialff, data = test_data)
      
      Number of observations: 100 
      
      Coefficients (Robust SE):
                  Estimate Std. Error z value Pr(>|z|)  
      (Intercept)   0.1030     0.2029  0.5078   0.6116  
      x            -0.2566     0.2093 -1.2261   0.2202  
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      
      Model fit:
        deviance: 136.8
        loglikelihood: -68.38
      

---

    Code
      summary(result_cl)
    Output
      Robust (Sandwich) Covariance Estimation for vglm
      ================================================
      
      Call:
      VGAM::vglm(formula = y ~ x, family = VGAM::binomialff, data = test_data)
      
      Number of observations: 100 
      Number of clusters: 10 
      
      Coefficients (Robust SE):
                  Estimate Std. Error z value Pr(>|z|)  
      (Intercept)   0.1030     0.1642  0.6276   0.5303  
      x            -0.2566     0.1569 -1.6353   0.1020  
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      
      Model fit:
        deviance: 136.8
        loglikelihood: -68.38
      

---

    Code
      summary(result_adj)
    Output
      Robust (Sandwich) Covariance Estimation for vglm
      ================================================
      
      Call:
      VGAM::vglm(formula = y ~ x, family = VGAM::binomialff, data = test_data)
      
      Number of observations: 100 
      Number of clusters: 10 
      Small-sample adjustment: applied (G/(G-1))
      
      Coefficients (Robust SE):
                  Estimate Std. Error z value Pr(>|z|)  
      (Intercept)   0.1030     0.1731  0.5953   0.5516  
      x            -0.2566     0.1654 -1.5514   0.1208  
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      
      Model fit:
        deviance: 136.8
        loglikelihood: -68.38
      

