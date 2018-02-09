#make clean
make 


# ./gctb       --bfile ../test2/ukb_subset_height_chr17_sub \
#              --make-full-ldm \
#              --out ../test2/ukb_subset_height_chr17_sub 


# ../bin/amber --bfile ~/Dropbox/Post_Doc_QBI/summary_stats_prob/amber/test2/ukb_subset_height_chr17_sub \
#  			 --make-ldm 1\
#   			 --wind 1  \
#  			 --out  ~/Dropbox/Post_Doc_QBI/summary_stats_prob/amber/test2/ukb_subset_height_chr17_sub 

./gctb       --bfile ../test2/ukb_subset_height_chr17_sub \
             --bayes C \
             --seed 2345 \
             --pheno ../test2/ukb_subset_height_chr17_sub.phen \
             --out ../test2/out/ukb_subset_height_chr17_sub_C_full \
             --chain-length 1000 --burn-in 100 --out-freq 1 2>&1 | tee ../test2/out/ukb_subset_height_chr17_sub_C_full.log

./gctb       --bfile ../test2/ukb_subset_height_chr17_sub \
             --bayes R \
             --seed 2345 \
             --pheno ../test2/ukb_subset_height_chr17_sub.phen \
             --out ../test2/out/ukb_subset_height_chr17_sub_R_full \
             --chain-length 1000 --burn-in 100 --out-freq 1 2>&1 | tee ../test2/out/ukb_subset_height_chr17_sub_R_full.log

              Mean            SD             
       Pi1    0.112507        0.066887       
       Pi2    0.887493        0.066887       
    NNZsnp    89.077774       6.736749       
   SigmaSq    0.857834        0.145436       
    ResVar    26.921997       0.066750       
    GenVar    26.877459       0.092768       
       hsq    0.499585        0.000992  

./gctb       --bfile ../test2/ukb_subset_height_chr17_sub \
             --bayes R \
             --seed 2345 \
             --set-pis 0.92,0.04,0.02,0.02 \
             --set-gammas 0,0.01,0.02,1 \
             --pheno ../test2/ukb_subset_height_chr17_sub.phen \
             --out ../test2/out/ukb_subset_height_chr17_sub_R_full \
             --chain-length 1000 --burn-in 1 --out-freq 1 2>&1 | tee ../test2/out/ukb_subset_height_chr17_sub_R_full.log

# Test 2 
./gctb       --bayes Cap \
             --bfile ../test2/ukb_subset_height_chr17_sub \
             --wind 100 \
             --seed 2345 \
             --pheno ../test2/ukb_subset_height_chr17_sub.phen \
             --out ../test2/out/ukb_subset_height_chr17_sub_Cap \
             --chain-length 1000 --burn-in 100 --out-freq 1 2>&1 | tee ../test2/out/ukb_subset_height_chr17_sub_Cap.log

./gctb       --sbayes C \
             --ldm ../test2/ukb_subset_height_chr17_sub.ldm.full \
             --gwas-summary ../test2/ukb_subset_height_chr17_sub.ma \
             --out ../test2/out/ukb_subset_height_chr17_sub_C \
             --chain-length 1000 --burn-in 100 --out-freq 1 2>&1 | tee ../test2/out/ukb_subset_height_chr17_sub_C.log

./gctb       --sbayes R \
             --ldm ../test2/ukb_subset_height_chr17_sub.ldm.full \
             --gwas-summary ../test2/ukb_subset_height_chr17_sub.ma \
             --set-pis 0.93,0.03,0.02,0.01,0.01 \
             --set-gammas 0,0.001,0.01,0.1,1.0 \
             --out ../test2/out/ukb_subset_height_chr17_sub_R \
             --seed 2345 \
             --chain-length 1000 --burn-in 100 --out-freq 1 2>&1 | tee ../test2/out/ukb_subset_height_chr17_sub_R.log
# ../bin/amber --sbayes C \
#              --ldmText ~/Dropbox/Post_Doc_QBI/summary_stats_prob/amber/test2/ukb_subset_height_chr17_sub.ldm.w1mb \
#              --gwas-summary ~/Dropbox/Post_Doc_QBI/summary_stats_prob/amber/test2/ukb_subset_height_chr17_sub.ma \
#              --out ~/Dropbox/Post_Doc_QBI/summary_stats_prob/amber/test2/out/ukb_subset_height_chr17_sub_C2 \
#              --chain-length 1000 --burn-in 100 --out-freq 1 2>&1 | tee ~/Desktop/file2.txt
# ../bin/amber --sbayes C \
#              --ldmText ~/Dropbox/Post_Doc_QBI/summary_stats_prob/amber/test2/ukb_subset_height_chr17_sub.ldm.w1000mb \
#              --gwas-summary ~/Dropbox/Post_Doc_QBI/summary_stats_prob/amber/test2/ukb_subset_height_chr17_sub.ma \
#              --out ~/Dropbox/Post_Doc_QBI/summary_stats_prob/amber/test2/out/ukb_subset_height_chr17_sub_C2 \
#              --chain-length 10000 --burn-in 1000 --out-freq 10 2>&1 | tee ~/Desktop/file2.txt
# ../bin/amber --sbayes R\
#              --ldm ~/Dropbox/Post_Doc_QBI/summary_stats_prob/amber/test2/ukb_subset_height_chr17_sub.ldm.w1000mb\
#              --gwas-summary ~/Dropbox/Post_Doc_QBI/summary_stats_prob/amber/test2/ukb_subset_height_chr17_sub.ma \
#              --out ~/Dropbox/Post_Doc_QBI/summary_stats_prob/amber/test2/out/ukb_subset_height_chr17_sub_R2 \
#              --chain-length 10000 --burn-in 1000 --out-freq 1  2>&1 | tee ~/Desktop/file1.txt
# ../bin/amber --sbayes R \
#             --ldmText ~/Dropbox/Post_Doc_QBI/summary_stats_prob/amber/test2/ukb_subset_height_chr17_sub.ldm.w1000mb \
#             --gwas-summary ~/Dropbox/Post_Doc_QBI/summary_stats_prob/amber/test2/ukb_subset_height_chr17_sub.ma \
#             --out ~/Dropbox/Post_Doc_QBI/summary_stats_prob/amber/test2/out/ukb_subset_height_chr17_sub_C2 \
#             --chain-length 10000 --burn-in 1000 --out-freq 1    2>&1 | tee ~/Desktop/file2.txt            
            
# Test 5

# ../bin/amber   --sbayes R \
#                --ldmText ../test5/ukbEURu_imp_chr22_v2_HM3_QC_sub_2_full_matlab_shrunk_deep.ldm.w100mb \
#                --gwas-summary ../test5/sim_1.ma \
#                --out ../test5/out/sim_1 \
#                --chain-length 10000 --burn-in 1000 --out-freq 100           
#               Mean            SD             
#        Pi1    0.242855        0.096095       
#        Pi2    0.022395        0.014326       
#        Pi3    0.333028        0.126290       
#        Pi4    0.401723        0.217374       
#     NNZsnp    1017.599976     554.014038     
#    SigmaSq    0.000183        0.000049       
#     ResVar    0.976920        0.003504       
#     GenVar    0.023334        0.002918       
#        hsq    0.023328        0.002914      
# 
#                   Mean            SD             
#         Pi    0.172657        0.031734       
#     NNZsnp    446.000000      92.474831      
#    SigmaSq    0.000241        0.000067       
#     ResVar    0.979676        0.002367       
#     GenVar    0.021447        0.001344       
#        hsq    0.021422        0.001333  


#               Mean            SD             
#        Pi1    0.047201        0.035529       
#        Pi2    0.048341        0.043436       
#        Pi3    0.092302        0.062154       
#        Pi4    0.812156        0.071035       
#     NNZsnp    82.800003       7.152599       
#    SigmaSq    1.031379        0.196127       
#     ResVar    25.912815       2.475392       
#     GenVar    27.862066       2.483257       
#        hsq    0.518117        0.046082     

# Test 2
#               Mean            SD             
#        Pi1    0.047201        0.035529       
#        Pi2    0.048341        0.043436       
#        Pi3    0.092302        0.062154       
#        Pi4    0.812156        0.071035       
#     NNZsnp    82.800003       7.152599       
#    SigmaSq    1.031379        0.196127       
#     ResVar    25.912815       2.475392       
#     GenVar    27.862066       2.483257       
#        hsq    0.518117        0.046082

#         Pi    0.639351        0.067780       
#     NNZsnp    68.800003       4.853820       
#    SigmaSq    1.224876        0.231907       
#     ResVar    25.964849       2.386140       
#     GenVar    27.815754       2.391799       
#        hsq    0.517203        0.044408    

#               Mean            SD             
#         Pi    0.897043        0.072625       
#     NNZsnp    90.537781       6.807624       
#    SigmaSq    1.289656        0.221877       
#     ResVar    26.720501       0.137109       
#     GenVar    27.056072       0.122783       
#        hsq    0.503121        0.002342       
# 
# Posterior summary:
# 
#               Mean            SD             
#         Pi    0.897043        0.072625       
#     NNZsnp    90.537781       6.807624       
#    SigmaSq    1.289656        0.221877       
#     ResVar    26.720501       0.137109       
#     GenVar    27.056072       0.122783       
#        hsq    0.503121        0.002342  

# LD matrix
# 
# Posterior summary:
# 
#               Mean            SD             
#         Pi    0.897043        0.072625       
#     NNZsnp    90.537781       6.807624       
#    SigmaSq    1.289656        0.221877       
#     ResVar    26.720501       0.137109       
#     GenVar    27.056072       0.122783       
#        hsq    0.503121        0.002342       



# Test 2 updated


#               Mean            SD             
#        Pi1    0.288090        0.056252       
#        Pi2    0.025436        0.024749       
#        Pi3    0.050920        0.046782       
#        Pi4    0.635554        0.064867       
#     NNZsnp    70.919998       4.018086       
#    SigmaSq    113.066963      24.477987      
#     ResVar    26.742041       0.110209       
#     GenVar    27.036072       0.089759       
#        hsq    0.502734        0.001777  

# When using the LD from text

#               Mean            SD             
#        Pi1    0.280067        0.065246       
#        Pi2    0.026531        0.026712       
#        Pi3    0.045929        0.043027       
#        Pi4    0.647474        0.070045       
#     NNZsnp    71.887779       5.025618       
#    SigmaSq    106.779716      21.951406      
#     ResVar    26.748009       0.112944       
#     GenVar    27.030994       0.092768       
#        hsq    0.502632        0.001819    
       
