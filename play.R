#!/usr/bin/r

#require(ANTsR)
#require(RANTs)
require(Rvcg)
require(Morpho)
require(mesheR)
require(rgl)
require(RvtkStatismo)
require(OpenMPController)

omp_set_num_threads(parallel::detectCores()) # use all the cores


landmarks = read.lmdta(file="lmply.dta") # $arr[,,filename]  ,  $idnames
representer = vcgImport(paste0("ply/",landmarks$idnames[1]))

modellm = landmarks$arr[,,landmarks$idnames[1]]

if(!exists("model")) {
    model = statismoModelFromPrepresenter(representer=representer,kernel=list(c(50,50)),ncomp=123)
} else {
    print("reusing existing model")
}

#Bayes = createBayes(model,sdmax = rep(4,50),wt=1.5,shrinkfun = function(x,i){ x = x*0.93^i})
#similarity = list(iterations=10,rhotol=pi/2)
#affine = list(iterations=10,rhotol=pi/2)

registerfun = function(i) {
 #   matchd = gaussMatch(Bayes,
 #                        target,
 #                        lm1=landmarks$arr[,,landmarks$idnames[1]],
 #                        lm2=landmarks$arr[,,landmarks$idnames[i]],
 #                        iterations=42,
 #                        sigma=30,
 #                        gamma=2,
 #                        toldist=30,
 #                        angtol=pi/2,
 #                        nh=100,
 #                        visualize=F,
 #                        similarity=similarity,
 #                        affine=affine)

    targetname = landmarks$idnames[i]
    target = vcgImport(paste0("ply/",targetname),silent=T);
    targetlm = landmarks$arr[,,targetname]
    transformation = computeTransform(modellm, targetlm, type="r")
    targetrot = applyTransform(target, transformation)
    


    print(paste("done w/",landmarks$idnames[i]))
}

if(F) {
    omp_set_num_threads(1) # needed for mclapply pattern
    out = parallel::mclapply(2:length(landmarks$idnames),registerfun,mc.cores=8)
    omp_set_num_threads(parallel::detectCores()) # use all the cores
} else {
    registerfun(2)
}


