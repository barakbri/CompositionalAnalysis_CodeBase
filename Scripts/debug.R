
L1O_res = leave_one_out_validation(X =  data_generation$X,
                                   Y =  data_generation$Y,
                                   ref_obj_for_validation = ref_obj,
                                   test = DACOMP.TEST.NAME.WILCOXON,
                                   Q = 0.1,
                                   NR_PERMS = 1000,
                                   Verbose = F,disable_DSFDR = F)