#!/usr/bin/r


require(Rvcg)
require(Morpho)
require(mesheR)
require(RvtkStatismo)

if (!require(getopt))
    install.packages("getopt",repos="http://cran.rstudio.com/",quiet=TRUE)
require(getopt)
spec = matrix(c(
    'skull',     's', 1, "character",
    'landmarks', 'l', 1, "character",
    'output' ,   'o', 1, "character",
    'iterations', 'i', 2, "integer",
    'help'   , 'h', 0, "logical"
    ), byrow=TRUE, ncol=4);
opt <- getopt(spec)
if ( !is.null(opt$help) || is.null(opt$landmarks) || is.null(opt$skull) ) {
       cat(getopt(spec, usage=TRUE,,command="zygoForDave.r"));
       cat("   estimates a pointcloud of the soft-tissue nose\n")
       q(status=1);
   }

right <- FALSE
skulllm <- read.pts(opt$landmarks)
modlm <- read.pts("data/zlmr.pts")
modlmmirr <- mirror(modlm)
##get references
lmleft <- pcAlign(skulllm,modlm,iterations = 100)
lmright <-  pcAlign(skulllm,modlmmirr,iterations = 100)
##get associations
leftkd <- vcgKDtree(lmleft,modlm,k=1)
rightkd <- vcgKDtree(lmright,modlmmirr,k=1)
##use the config with the lower distance
if (mean(leftkd$distance) < mean(rightkd$distance)){
    skulllm <- skulllm[leftkd$index,]
} else {
    skulllm <- skulllm[rightkd$index,]
    right <- TRUE
    cat("using righthand model\n")
}


if (right) {
    cat("using right hand model\n")
    newmodRed <- statismoLoadModel("data/zrmrLowRes.h5")
    leftmean <- vcgImport("data/zlmr.ply")
    rightmean <- vcgImport("data/zrmr.ply")
    modlm <- transferPoints(modlm,leftmean,invertFaces(rightmean))
} else {
    newmodRed <- statismoLoadModel("data/zlmrLowRes.h5")
}
iterations <- opt$iterations
skull <- vcgQEdecim(vcgImport(opt$skull),edgeLength = 1)

myfun <- function(skull,skulllm,sdmax=NULL,iterations=10,icp=TRUE,matching="r",...) {
    initrot <- computeTransform(modlm,skulllm,type = "r")
    skullrot <- skullrot2 <- applyTransform(skull,initrot)
    skullrotlm <- skullrotlm2 <- applyTransform(skulllm,initrot)
    
    if (icp) {
        zygoicp <- icp(DrawMean(newmodRed),skullrot,iterations = 30,getTransform = T,type=matching,subsample=1000,rho=pi/2,uprange = 0.5)
        skullrot2 <- applyTransform(skullrot,zygoicp$transform,inverse = T)
        skullrotlm2 <- applyTransform(skullrotlm,zygoicp$transform,inverse = T)
    }
    cModRed <- statismoConstrainModel(newmodRed,skullrotlm2,modlm,ptValueNoise = 2)
                                        #wiredeform <- modelFitting(cMod,skullrot2,tardist = 4,refdist = 6,iterations = 10)
    wiredeformRed <- modelFitting(cModRed,skullrot2,tardist = 3,refdist = 6 ,iterations = iterations,sdmax = sdmax,mahaprob = "c",lbfgs.iter = 15,...)
    if (icp)
        finalRed <- applyTransform(wiredeformRed$mesh,zygoicp$transform)
    else
        finalRed <- wiredeformRed$mesh
    finalorigRed <- applyTransform(finalRed,initrot,inverse = T)
    return(list(estimate=finalorigRed,mod=cModRed,par=wiredeformRed$par,skullrot=skullrot2,estimrot=wiredeformRed$mesh))
}

out <- myfun(skull,skulllm,iterations=iterations,sdmax=4,orthantwise_c=0.01)
vcgPlyWrite(out$estimate,opt$output)
vtkMeshWrite(out$skullrot,paste0(opt$output,"_skull_rotated"))
vtkMeshWrite(out$estimrot,paste0(opt$output,"_estimation_rotated"))
statismoSaveModel(out$mod,modelname=paste0(opt$output,".h5"))
write.table(out$par,col.names=F,row.names=F,file=paste0(opt$output,"_params.csv"))
