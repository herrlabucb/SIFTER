########## 
#libraries
##########
library(Matrix)
library(MASS)
library(parallel)
library(kernlab)
library(ggplot2)
library(reshape2)
library(mclust)
library(flexclust)
library(gplots)
library(proxy)
library(grid)
library(gridExtra)
library(entropy)


########## 
#functions
##########
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

Cal_Similarity = function(criterion='none',gene_expression,method='pearson',
                          p_Lp=1,c_Lp=1,bin_width=0.1){
  #' calculate the similarity between columns of gene_expression
  #' pearson, spearman, kendall, Euclidean 
  
  if( !(method %in% c('pearson','spearman','kendall')) ){
    len = ncol(gene_expression)
    norm.len = nrow(gene_expression)
    tmp_pos = expand.grid(1:len,1:len)
    tmp_pos = as.list(as.data.frame(t(tmp_pos)))
  }
  
  if(method %in% c('pearson','spearman','kendall')){
    tmp = cor(gene_expression,method=method) 
    if(criterion == 'none'){
      return (exp((tmp-1)/c_Lp))
    } else if (criterion == 'abs'){
      return (abs(tmp))
    }
    
  } else if(method == 'Euclidean'){
    tmp = matrix(unlist(lapply(tmp_pos,FUN=function(x){
      sqrt(sum((gene_expression[,x[1]] - gene_expression[,x[2]])^2))
    })), nc=len,byrow = T)
    return (exp(-tmp^2/c_Lp/norm.len))
    
  }  
  
}


ncores = 1
n_specc = 2
pool_count = 2 
repeat.time = 3000


########## 
#read data
##########
dt_control = t(read.csv(paste0('./DMSO_ctrl_all_data', '.csv'), 
                     header = TRUE, #sep = ',',
                     stringsAsFactors = TRUE))
dt_treat = t(read.csv(paste0('./LatA_all_data', '.csv'), 
                    header = TRUE, #sep = ',',
                    stringsAsFactors = TRUE))
colnames(dt_control) = paste0('control_', 1:ncol(dt_control))
colnames(dt_treat) = paste0('treat_', 1:ncol(dt_treat))

# standardize together
dt_tmp = t(scale(t(cbind(dt_treat, dt_control))))[c('MT', 'F.actin', 'IF'), ]
dt_treat = dt_tmp[, 1:ncol(dt_treat)]
dt_control = dt_tmp[, (1 + ncol(dt_treat)) : (ncol(dt_treat) + ncol(dt_control))]

bait_file_wdir = './BAITS_nonoverlap/'
file_name = list.files(bait_file_wdir)


for(tmp_name in file_name){
  print(tmp_name)
  idx_check = read.csv(paste0(bait_file_wdir, tmp_name), 
                       header = FALSE, #sep = ',',
                       stringsAsFactors = TRUE)
  kk = 1
  BAIT_idx = idx_check[, kk]
  BAIT_idx = BAIT_idx[!is.na(BAIT_idx)]
  
  gene_bait = dt_treat[, BAIT_idx]
  gene_true = dt_treat[, -BAIT_idx]
  gene_false = dt_control

  n_true = ncol(gene_true)
  n_false = ncol(gene_false) 
  n_bait = ncol(gene_bait) 

  ##########
  ## CFR
  ##########
  Cal_CFR = function(data_fish){
    
    CF = mclapply(1:repeat.time, function(ii){
      CF = rep(0, ncol(data_fish))
      names(CF) = c(colnames(data_fish))
      
      idx.tmp = sample.int(ncol(data_fish), replace = F)

      tmp = floor(length(idx.tmp) * (0:pool_count)/pool_count)
      idx.tmp = lapply(1 : (length(tmp)-1), function(i){
        idx.tmp[ (tmp[i]+1) : tmp[i+1] ]
      })
      
      for (i in 1:length(idx.tmp)){
        gene_pool = data_fish[, idx.tmp[[i]]]

        W = Cal_Similarity(criterion = 'abs',method = 'Euclidean',
                           gene_expression = cbind(gene_bait, gene_pool))
        diag(W) = 0
        G.inverse.root = diag(sqrt(1 / rowSums(W)))
        L = diag(ncol(W)) - G.inverse.root %*% W %*% G.inverse.root
        rownames(L) = colnames(cbind(gene_bait, gene_pool)) 
        eigen_vec = eigen(L)$vectors
        rownames(eigen_vec) = rownames(L)
        
        eigen_val_idx = min(which(eigen(L)$values < 1e-6)) - 1 
        eigen_val_leftidx = max(1, (eigen_val_idx - n_specc + 1))
        
        BIC <- mclustBIC(data = eigen_vec[, eigen_val_leftidx:eigen_val_idx], G = 2)
        mod1 <- Mclust(data = eigen_vec[, eigen_val_leftidx:eigen_val_idx], x = BIC)
        pred_res = mod1$classification   

        data.tmp = as.data.frame(eigen_vec[, (eigen_val_idx-1):eigen_val_idx]) 
        data.tmp$group = factor(pred_res)

        fish_group = getmode(pred_res[names(pred_res) %in% colnames(gene_bait)])

        candidate.gene.tmp = pred_res[!(names(pred_res) %in% colnames(gene_bait))]
        if( sum(candidate.gene.tmp == fish_group) > 
            sum(candidate.gene.tmp != fish_group) ){
          CF[colnames(gene_pool)] = NA
          next
        }
        
        fish_gene = names(which(pred_res == fish_group))
        fish_gene = fish_gene[!(fish_gene %in% colnames(gene_bait))]
        
        CF[fish_gene] = CF[fish_gene] + 1
      }
      
      return (CF)
    }, mc.cores = ncores) 
    CI = do.call(rbind, CF)
    return(CI)
  }
  
  CI_all = vector(mode = 'list', length = 3)
  CI_all[[1]] = ncol(gene_true)
  CI_all[[2]] = Cal_CFR(data_fish = gene_false)
  CI_all[[3]] = ncol(cbind(gene_true, gene_false))
  
  save(CI_all, file = paste0('./Fishing_result/', 
                         strsplit(tmp_name, '[.]')[[1]][1], 
                         'CI',  '.Rda'))
}