####### Pioneer: Erica Baller
#### Date: 3/9/2021
####### Appropriator: Adam Pines
#### Date: 3/9/2022


### Will Need Spins in FS4
### Will need conversion of yeo7 to faces

##pre: right and left 10242 x 1000 matrices from matlab SpinPermuFS, yeo R & L assignments
##post: 2 7 x 1000 matrices (r & l) that contain the proportion of vertices within a network divided by the total number of vertices, and plots
## uses: Takes output of spin test, and calcualted the number of vertices within each of yeo's 7 networks out of the number of total possible vertices within the network
    #### 1) Read in the yeo network assignments and calculate total number of vertices per network
    #### 2) Multiply the yeo networks x the matrices (so every value is 1 -7 if they were within the mask, -1--7 if they were medial wall, and 0 otherwise)
    #### 3) Foreach permutation (r and l separately), and for each network, calculate the (# of vertices with a 1) divided(/) by the (number of total vertices within network minus number of negative vertices
    #### 4) Store
    #### 5) Plot

### dependencies: ggplot2, bigmemory, vroom


#library(bigmemory.sri)
library(ggplot2)
library(tidyr)

#################
### set home directory
source(paste0(homedir, "/baller/scripts/imco_functions.R"))

#initialize
hemis <- c("lh", "rh")
permNum <- 1000
yeo_num <- 7

### Read in matrices 
  
# real data
lh_t_fdr05_results <- t(read.table('/cbica/projects/pinesParcels/results/aggregated_data/AgeEfs_fs4_L.csv', sep = ","))
rh_t_fdr05_results <- t(read.table('/cbica/projects/pinesParcels/results/aggregated_data/AgeEfs_fs4_R.csv', sep = ","))
 
# spins                                         
lh_spin <-t(read.csv('~/data/lh_spin_test__AgeEf__output.csv',col.names=T))
rh_spin <-t(read.csv('~/data/rh_spin_test__AgeEf__output.csv',col.names=T))

# intiialize permutation matrix
permMat_L=matrix(0,2562,1000)
permMat_R=matrix(0,2562,1000)

# index out each of 1k permutations
for (p in seq(1,1000)){
	start=((p-1)*2562)+1
	finish=(p*2562)
	permMat_L[,p]=lh_spin[start:finish]
	permMat_R[,p]=rh_spin[start:finish]
}

#bring together, with original values as first column
lh_act_results_and_spin <- cbind(lh_t_fdr05_results, permMat_L)
rh_act_results_and_spin <- cbind(rh_t_fdr05_results, permMat_R)

# get ind of where 1000s are for setting to negative to match Erica's script later
lh_1k_inds=which(lh_act_results_and_spin==1000,arr.ind=T)
rh_1k_inds=which(rh_act_results_and_spin==1000,arr.ind=T)
# set 1000s back to 0 (1000s from spin script)
lh_act_results_and_spin[lh_act_results_and_spin==1000]=0
rh_act_results_and_spin[rh_act_results_and_spin==1000]=0
  
#grab list of yeo 7 networks in fsaverage4 space
lh_yeo_network <- read.csv('~/data/y7_R_3k.csv')[,1]
rh_yeo_network <- read.csv('~/data/y7_L_3k.csv')[,1]
  
#count up number of vertices per network
lh_yeo_network_count_table <- table(lh_yeo_network)
rh_yeo_network_count_table <- table(rh_yeo_network)
  
#multiply yeo network x spin test
# AP note: I think Erica made this a binarized "is there an effect here map" for it to work with yeo7 # multiplication
lh_act_results_and_spinBool=sapply(as.data.frame(lh_act_results_and_spin),as.logical)
rh_act_results_and_spinBool=sapply(as.data.frame(rh_act_results_and_spin),as.logical)

# set 1000s to negative to match Erica's script
lh_act_results_and_spinBool[lh_1k_inds]=-1
rh_act_results_and_spinBool[rh_1k_inds]=-1

# multiply them to turn boolean TRUE into network numbers
lh_spinxyeo <- lh_act_results_and_spinBool*lh_yeo_network
rh_spinxyeo <- rh_act_results_and_spinBool*rh_yeo_network

#proportions
#go through each hemisphere, go through each perm, and go through each network
  
lh_hemi_spin_proportions <- data.frame(matrix(nrow = yeo_num, ncol = (permNum + 1)))
rh_hemi_spin_proportions <- data.frame(matrix(nrow = yeo_num, ncol = (permNum + 1)))

for (hemi in hemis){
  for (perm in 1:(permNum + 1)){
    for (network in 1:yeo_num){
        
        #number of vertices within network that are fdr corrected
        num_pos_to_parse<- paste0("length(which(", hemi, "_spinxyeo[,", perm, "] == ", network, "))") 
        num_vertices_in_spin <- eval(parse(text = as.character(num_pos_to_parse)))
         
        #number of vertices within network that are negative (i.e., medial wall)
        num_neg_to_parse <- paste0("length(which(", hemi, "_spinxyeo[,", perm, "] == -", network, "))")
        num_neg <- eval(parse(text = as.character(num_neg_to_parse)))
        
        #total number of vertices in normal network : AP +1 because 0 hold first place in this table
        total_possible_to_parse <- paste0(hemi, "_yeo_network_count_table[", network, "+1]")
        total_possible <- eval(parse(text = as.character(total_possible_to_parse)))
        
        #proportion of vertices within network , with denominator being total possible by # in medial wall
        proportion_potential_vertices <- num_vertices_in_spin/(total_possible - num_neg)
        
        #store in matrix
        storing_to_parse <- paste0(hemi, "_hemi_spin_proportions[", network, ",", perm, "] = ", proportion_potential_vertices)
        eval(parse(text = as.character(storing_to_parse)))
    }
  }
}
  
write.table(lh_hemi_spin_proportions,'~/data/lh_spin_test_AgeEf_proportions.csv', sep = ",", col.names = F, row.names = F)
write.table(rh_hemi_spin_proportions,'~/data/rh_spin_test_AgeEf_proportions.csv', sep = ",", col.names = F, row.names = F)
#then plot
