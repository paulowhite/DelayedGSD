### test-NonBindingHJ.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan  7 2022 (14:48) 
## Version: 
## Last-Updated: jan  7 2022 (14:56) 
##           By: Brice Ozenne
##     Update #: 4
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(testthat)

test_that("Boundary calculation with non-binding futility following Jennison method",{
    ## slide 106 from CJ_DSBS-v5.pdf
    
    bCJ <- NonBindingHJ(rho_alpha=1.345,
                        rho_beta=1.345,
                        alpha=0.025,
                        beta=0.1,
                        Kmax=3,
                        Info.max=12,
                        InfoR.i=c(3.5,6.75,12)/12,
                        InfoR.d=c(5.5,8.75)/12,
                        delta=1,  
                        abseps = 1e-06, 
                        direction="smaller",
                        sided=1
                        )
    expect_equal(round(bCJ$boundaries$l.k,3), c(-0.409,0.664,2.069))
    expect_equal(round(bCJ$boundaries$u.k,3), c(2.437,2.244,2.069))
    expect_equal(round(bCJ$boundaries$c.k,3), c(1.960,1.960,NA))
        
})
##----------------------------------------------------------------------
### test-NonBindingHJ.R ends here
