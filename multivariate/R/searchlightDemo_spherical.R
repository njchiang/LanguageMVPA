make an error
###################################################################################################################################################################
# 13 May 2015.
# R code written by Joset A. Etzel (jetzel@artsci.wustl.edu, mvpa.blogspot.com) to perform an example single-subject searchlight analysis. This is an updated
# version of searchlight_demo.R. The differences are that .rds files are not used and the searchlights are shaped as in pyMVPA and the Princeton MVPA toolbox.
# This code may be adapted for personal use, provided the source is cited.
# The sample image (dataset.nii.gz) is at https://dl.dropboxusercontent.com/u/13098670/dataset.nii.gz
###################################################################################################################################################################
# make 3d.nii.gz: object with a 3d array (matching the size of the input image) of numbered voxels.
# the input image is 4d, with brain voxels nonzero. They have already been preprocessed.

library(oro.nifti);    # for readNIfTI and writeNIfTI
rm(list=ls());   # clear R's memory

img.path <- "d:/temp/searchlightDemo/";   # read in the image to get its size
out.path <- "d:/temp/searchlightDemo/withSpheres/";  # will write 3d.nii.gz into this directory

img <- readNIfTI(paste0(img.path, "dataset.nii.gz"), reorient=FALSE);   # read in one of the nifti images - they are all preprocessed and the same size.
dim(img);   # get the dimensions of the input images.
# [1] 41 48 35 20

inds <- which(img[,,,1] != 0, arr.ind=TRUE);  # find all the brain voxel coordinates.
# NOTE: usually the brain voxels should be read from an anatomical mask (such as the entire brain), not the functional data, but this is quick for the demo.
if (ncol(inds) != 3) { stop('wrong dims on inds'); }

out.array <- array(0, dim(img[,,,1]));  # make a blank 3d matrix 'brain', all zeros
out.array[inds] <- 1:nrow(inds);  # and put integers into the voxel places

# these NIfTI header parameters will match dataset.nii.gz and roi.nii.gz
out.img <- nifti(out.array, datatype=16, dim=dim(out.array), srow_x=c(4,0,0,-80), srow_y=c(0,4,0,-112), srow_z=c(0,0,4,-52),
                         qoffset_x=-80, qoffset_y=-112, qoffset_z=-52, xyzt_units=2);
out.img@sform_code <- 4; 
out.img@qform_code <- 4;
pixdim(out.img)[1:8] <- c(1, 4,4,4, 1, 0,0,0);  # qFactor, mm/voxel, number of timepoints
writeNIfTI(out.img, paste0(out.path, "3d"));    # write out as a NIfTI image

###################################################################################################################################################################
# make the look-up table; one row for each voxel in the brain. This only needs to be run once for each radius and mask (here, 3d_wustl333.nii.gz).
# The file this makes is read during the searchlight analysis to determine which voxels are in the "surround" for each voxel (ie where the searchlights are).
# this version makes pyMVPA/Princeton toolbox-shaped searchlights: http://mvpa.blogspot.com/2012/09/more-searchlights-princeton-mvpa-toolbox.html

# the code in the get.surround function is not particularly elegant; there are a lot of checks to make sure it doesn't try putting searchlights outside of the brain.

library(oro.nifti);
rm(list=ls());

# x,y,z are the coordinates (ijk space) of the voxel that's the searchlight center we're getting the surround for
# do.radius is the searchlight radius, in voxels
# mask.img is a 3d array (same ijk as the coordinates) with 0 for non-brain voxels and unique integers for brain voxels.
# get.surround returns an array of the surrounding voxels, by name in mask.img (not coordinates)
get.surround <- function(x,y,z, do.radius, mask.img) {    # do.radius <- 3; mask.img <- img3d; x <- 46; y <- 30; z <- 11;
  if (do.radius > 3 | do.radius < 1) { stop("radius can be 1, 2, or 3 voxels only"); }
  if (as.integer(do.radius) != do.radius) { stop("radius must be an integer"); }
  
  if (do.radius == 1) { surrounds <- rep(0, 6); }   # 6 voxels in the surround for a 1-voxel radius
  if (do.radius == 2) { surrounds <- rep(0, 32); } 
  if (do.radius == 3) { surrounds <- rep(0, 122); } 
  ctr <- 1;  # surround counter
  
  # all searchlights have the radius 1 voxels, so do that first
  if (y-1 > 0) { surrounds[ctr] <- mask.img[x, (y-1), z]; ctr <- ctr + 1; } # same slice as center voxel, sharing a face with the center
  if (y+1 <= dim(mask.img)[2]) { surrounds[ctr] <- mask.img[x, (y+1), z]; ctr <- ctr + 1; }
  if (x-1 > 0) { surrounds[ctr] <- mask.img[(x-1), y, z]; ctr <- ctr + 1; }
  if (x+1 <= dim(mask.img)[1]) { surrounds[ctr] <- mask.img[(x+1), y, z]; ctr <- ctr + 1; }
  if (z-1 > 0) { surrounds[ctr] <- mask.img[x, y, (z-1)];  ctr <- ctr + 1; } 
  if (z+1 <= dim(mask.img)[3]) { surrounds[ctr] <- mask.img[x, y, (z+1)]; ctr <- ctr + 1; } # slice above and below the center voxel, sharing a face with the center
 
  # now, voxels in the radius 2 surround
  if (do.radius > 1) {
    if (z-2 > 0) { surrounds[ctr] <- mask.img[x, y, (z-2)]; ctr <- ctr + 1; }  # two layers above and below the center voxel
    if (z+2 <= dim(mask.img)[3]) { surrounds[ctr] <- mask.img[x, y, (z+2)]; ctr <- ctr + 1; }
    if (y+2 <= dim(mask.img)[2]) { surrounds[ctr] <- mask.img[x, (y+2), z]; ctr <- ctr + 1; } # same slice as center voxel
    if (y-2 > 0) { surrounds[ctr] <- mask.img[x, (y-2), z]; ctr <- ctr + 1; }
    if (x-2 > 0) { surrounds[ctr] <- mask.img[x-2, y, z]; ctr <- ctr + 1; }
    if (x+2 <= dim(mask.img)[1]) { surrounds[ctr] <- mask.img[x+2, y, z]; ctr <- ctr + 1; }
    if (x+1 <= dim(mask.img)[1] & y+1 <= dim(mask.img)[2]) { surrounds[ctr] <- mask.img[x+1, y+1, z]; ctr <- ctr + 1; }
    if (x+1 <= dim(mask.img)[1] & y-1 > 0) { surrounds[ctr] <- mask.img[x+1, y-1, z]; ctr <- ctr + 1;}
    if (x-1 > 0 & y+1 <= dim(mask.img)[2]) { surrounds[ctr] <- mask.img[x-1, y+1, z]; ctr <- ctr + 1; }
    if (x-1 > 0 & y-1 > 0) { surrounds[ctr] <- mask.img[x-1, y-1, z]; ctr <- ctr + 1; }
    
    for (z.offset in c(-1, 1)) {  # layers above and below the center voxel
      if (z+z.offset <= dim(mask.img)[3] & z+z.offset > 0) { 
        if (x+1 <= dim(mask.img)[1]) { surrounds[ctr] <- mask.img[x+1, y, (z + z.offset)]; ctr <- ctr + 1; }
        if (x-1 > 0) { surrounds[ctr] <- mask.img[x-1, y, (z + z.offset)]; ctr <- ctr + 1; }
        if (x+1 <= dim(mask.img)[1] & y+1 <= dim(mask.img)[2]) { surrounds[ctr] <- mask.img[x+1, y+1, (z + z.offset)]; ctr <- ctr + 1; }
        if (y+1 <= dim(mask.img)[2]) { surrounds[ctr] <- mask.img[x, y+1, (z + z.offset)]; ctr <- ctr + 1; }
        if (x-1 > 0 & y+1 <= dim(mask.img)[2]) { surrounds[ctr] <- mask.img[x-1, y+1, (z + z.offset)]; ctr <- ctr + 1; }
        if (x+1 <= dim(mask.img)[1] & y-1 > 0) { surrounds[ctr] <- mask.img[x+1, y-1, (z + z.offset)]; ctr <- ctr + 1; }
        if (y-1 > 0) { surrounds[ctr] <- mask.img[x, y-1, (z + z.offset)]; ctr <- ctr + 1; }
        if (x-1 > 0 & y-1 > 0) { surrounds[ctr] <- mask.img[x-1, y-1, (z + z.offset)]; ctr <- ctr + 1; }
      } 
    }
  }


  # now, voxels in the radius 3 surround
  if (do.radius > 2) {
    if (z+3 <= dim(mask.img)[3]) { surrounds[ctr] <- mask.img[x, y, (z+3)]; ctr <- ctr + 1; }   # three layers above and below the center voxel
    if (z-3 > 0) { surrounds[ctr] <- mask.img[x, y, (z-3)]; ctr <- ctr + 1; }  # the brain is close to the edge, so this doesn't always work.
    if (y+3 <= dim(mask.img)[2]) { surrounds[ctr] <- mask.img[x, (y+3), z]; ctr <- ctr + 1; }    # same slice as center voxel
    if (y-3 > 0) { surrounds[ctr] <- mask.img[x, (y-3), z];  ctr <- ctr + 1; }
    if (x-3 > 0) { surrounds[ctr] <- mask.img[x-3, y, z]; ctr <- ctr + 1; }       # sharing faces
    if (x+3 <= dim(mask.img)[1]) { surrounds[ctr] <- mask.img[x+3, y, z]; ctr <- ctr + 1; }
     
    for (x.offset in c(-2, -1, 1, 2)) { 
      if (x+x.offset <= dim(mask.img)[1] & x+x.offset > 0) { 
        if (y-2 > 0) { surrounds[ctr] <- mask.img[x+x.offset, y-2, z]; ctr <- ctr + 1; }
        if (y+2 <= dim(mask.img)[2]) { surrounds[ctr] <- mask.img[x+x.offset, y+2, z]; ctr <- ctr + 1; }
      }
    }
    for (x.offset in c(-2, 2)) { 
      if (x+x.offset <= dim(mask.img)[1] & x+x.offset > 0) { 
        if (y-1 > 0) { surrounds[ctr] <- mask.img[x+x.offset, y-1, z]; ctr <- ctr + 1; }
        if (y+1 <= dim(mask.img)[2]) { surrounds[ctr] <- mask.img[x+x.offset, y+1, z]; ctr <- ctr + 1; }
      }
    }
     
    for (z.offset in c(-1, 1)) {  # layers above and below the center voxel
      if (z+z.offset <= dim(mask.img)[3] & z+z.offset > 0) { 
        for (x.offset in c(-2, -1, 0, 1, 2)) {
          if (x+x.offset <= dim(mask.img)[1] & x+x.offset > 0) { 
            if (y+2 <= dim(mask.img)[2]) { surrounds[ctr] <- mask.img[x+x.offset, y+2, (z + z.offset)]; ctr <- ctr + 1; }
            if (y-2 > 0) { surrounds[ctr] <- mask.img[x+x.offset, y-2, (z + z.offset)]; ctr <- ctr + 1; }
          }
        }
        for (x.offset in c(-2, 2)) {
          if (x+x.offset <= dim(mask.img)[1] & x+x.offset > 0) { 
            if (y+1 <= dim(mask.img)[2]) { surrounds[ctr] <- mask.img[x+x.offset, y+1, (z + z.offset)]; ctr <- ctr + 1; }
            surrounds[ctr] <- mask.img[x+x.offset, y, (z + z.offset)]; ctr <- ctr + 1; 
            if (y-1 > 0) { surrounds[ctr] <- mask.img[x+x.offset, y-1, (z + z.offset)]; ctr <- ctr + 1; }
          }
        }
      }
    } 
    
    for (z.offset in c(-2, 2)) {  # two layers above and below the center voxel
      if (z+z.offset <= dim(mask.img)[3] & z+z.offset > 0) { 
        for (x.offset in c(-2, -1, 0, 1, 2)) {
          if (x+x.offset <= dim(mask.img)[1] & x+x.offset > 0) { 
            if (y+1 <= dim(mask.img)[2]) { surrounds[ctr] <- mask.img[x+x.offset, y+1, (z + z.offset)]; ctr <- ctr + 1; }
            if (y-1 > 0) { surrounds[ctr] <- mask.img[x+x.offset, y-1, (z + z.offset)]; ctr <- ctr + 1; }
          }
        }
        for (x.offset in c(-2, -1, 1, 2)) { 
          if (x+x.offset <= dim(mask.img)[1] & x+x.offset > 0) { 
            surrounds[ctr] <- mask.img[x+x.offset, y, (z + z.offset)]; ctr <- ctr + 1;  
          }
        }
        for (x.offset in c(-1, 0, 1)) { 
          if (x+x.offset <= dim(mask.img)[1] & x+x.offset > 0) { 
            if (y+2 <= dim(mask.img)[2]) { surrounds[ctr] <- mask.img[x+x.offset, y+2, (z + z.offset)]; ctr <- ctr + 1; }
            if (y-2 > 0) { surrounds[ctr] <- mask.img[x+x.offset, y-2, (z + z.offset)]; ctr <- ctr + 1; }
          }
        }
      }
    } 
  }

  return(surrounds);
}


in.path <- "d:/temp/searchlightDemo/withSpheres/";   # location of 3d.nii.gz, as made in the previous block of code
out.path <- "d:/temp/searchlightDemo/withSpheres/";      # where the 'lookup' file will be written
need.radius <- 2;  # the voxel radius we want in the searchlights, as an integer. only 1, 2, and 3 voxel radii are possible.

# get the number of voxels and dimensionality from img3d.
img3d <- readNIfTI(paste0(in.path, "3d"), reorient=FALSE);  # unique integers for brain voxels, 0 elsewhere. MUST match the functional images.
num.roi.vox <- max(img3d);
inds <- which(img3d != 0, arr.ind=TRUE); 
if (dim(inds)[1] != num.roi.vox) { stop("dim(inds)[1] != num.roi.vox"); }

# start filling up the lookup list: the adjacent voxels for each non-zero voxel in the brain mask
if (need.radius == 1) { num.surround <- 6; }    # 6 voxels in the surround for a 1-voxel radius
if (need.radius == 2) { num.surround <- 32; }   # 32 voxels in the surround for a 2-voxel radius
if (need.radius == 3) { num.surround <- 122; }  # 122 voxels in the surround for a 3-voxel radius
  
# actually fill up the lookup table. This goes through every voxel, so can take some time.
lookup.tbl <- array(NA, c(nrow(inds), num.surround));  # list that will hold the surrounds for each center voxel
for (i in 1:nrow(inds)) {    # i <- 1;
  these.surrounds <- get.surround(inds[i,1], inds[i,2], inds[i,3], need.radius, img3d);
  if (length(these.surrounds) != num.surround) { stop("length(these.surrounds != num.surround)"); }
  vox.id <- img3d[inds[i,1], inds[i,2], inds[i,3]];
  lookup.tbl[vox.id,] <- sort(these.surrounds);
}

# replace the zeroes with NA so can sort better at runtime
inds <- which(lookup.tbl == 0, arr.ind=TRUE)
lookup.tbl[inds] <- NA;
write.table(lookup.tbl, gzfile(paste0(out.path, "lookup_radius", need.radius, ".txt.gz")));   # write the table, gzipped to save file size

###################################################################################################################################################################
# confirm that the searchlights are the expected shape by reading a few out of the new lookup table. This image can be overlaid on the functional or 3d image.

rm(list=ls());   # clear R's memory

in.path <- "d:/temp/searchlightDemo/withSpheres/";   # location of lookup_radius2.txt.gz and 3d.nii.gz
out.path <- "d:/temp/searchlightDemo/";  # where to write the test image
need.radius <- 2;  # searchlight radius for the lookup table we're checking

img.3d <- readNIfTI(paste0(in.path, "3d.nii.gz"), reorient=FALSE);
lookup.tbl <- read.table(gzfile(paste0(in.path, "lookup_radius", need.radius, ".txt.gz"))); 
# dim(lookup.tbl)   [1] 17582    32    # there are 17582 voxels in the demo dataset's lookup table. we can make searchlights out of any of them.
do.searchlights <- sample(1:nrow(lookup.tbl))[1:10];   # pick 10 searchlights to plot at random.

out.array <- array(NA, dim(img.3d));  # make a blank 3d matrix 'brain' to put the searchlights into
# put integers into out.array for each searchlight
for (i in 1:length(do.searchlights)) {    # i <- 1;
  voxs <- union(do.searchlights[i], unlist(lookup.tbl[do.searchlights[i],], use.names=FALSE)); # surrounding voxels for i; union to put center in the searchlight
  voxs <- voxs[which(!is.na(voxs))];  # get rid of NA.
  for (j in 1:length(voxs)) {    # j <- 1;
    coords <- which(img.3d == voxs[j], arr.ind=TRUE);   # location of this voxel
    if (ncol(coords) != 3 | nrow(coords) != 1) { stop("wrong sized coords"); }
    out.array[coords[1], coords[2], coords[3]] <- i;   # assign this voxel the searchlight number
  }
}

# these NIfTI header parameters will match dataset.nii.gz and roi.nii.gz
out.img <- nifti(out.array, datatype=16, dim=dim(img.3d), srow_x=c(4,0,0,-80), srow_y=c(0,4,0,-112), srow_z=c(0,0,4,-52),
                         qoffset_x=-80, qoffset_y=-112, qoffset_z=-52, xyzt_units=2);
out.img@sform_code <- 4; 
out.img@qform_code <- 4;
pixdim(out.img)[1:8] <- c(1, 4,4,4, 1, 0,0,0);  # qFactor, mm/voxel, number of timepoints
writeNIfTI(out.img, paste0(out.path, "searchlightTest"));    # write out as a NIfTI image

###################################################################################################################################################################
# now the preparations are finished, so we can actually run the searchlight analysis.
###################################################################################################################################################################
# do the searchlight analysis. Since this is a demo, it runs a very simple cross-validation scheme - not recommended for actual analyses!
# this code runs the searchlight analysis in four "chunks": each chunk is a separate part of the brain. I run each chunk as a separate job on a supercomputing
# cluster, then collect and combine the individual output files. Obviously, how the paths and arguments are set will vary with cluster computer.
# each chunk takes around an hour to run on my local machine; code speed-ups are probably possible.

library(oro.nifti);   # to read and write NIfTIs
library(e1071);   # for the linear svm

rm(list=ls()); on.cluster <- TRUE; 
# rm(list=ls()); on.cluster <- FALSE; 

if (on.cluster == TRUE) {   # get the chunk from the argument sent to R when the job was started.
  cA <- commandArgs();
  do.chunk  <- as.numeric(cA[5]);   # chunk to do in this job, sent as an argument when R was started

  img.path <- "/scratch/jetzel/searchlightDemo/";   # dataset.nii.gz
  spot.path <- "/scratch/jetzel/searchlightDemo/";  # 3d.nii.gz and lookup_radius2.txt.gz directory
  out.path <- "/scratch/jetzel/searchlightDemo/";   # write into here
} else {   # not on the cluster, so set the chunk here.
  img.path <- "d:/temp/searchlightDemo/";   # dataset.nii.gz
  spot.path <- "d:/temp/searchlightDemo/withSpheres/";  # 3d.nii.gz and lookup_radius2.txt.gz directory
  out.path <- "d:/temp/";   # write into here
  
  do.chunk <- 1;  # four total: which part of the brain to run
}

radius <- 2;   # searchlight radius
smallest.searchlight <- 2;  # how many voxels (counting the center) must be in a searchlight for it to be run; smaller ones are skipped.
lookup <- read.table(gzfile(paste0(spot.path, "lookup_radius", radius, ".txt.gz")));    # searchlight surrounds
vol3d <- readNIfTI(paste0(spot.path, "3d.nii.gz"), reorient=FALSE);    # for going from voxel numbers to 3d coordinates.
all.img <- readNIfTI(paste0(img.path, "dataset.nii.gz"), reorient=FALSE);  # first 10 class 'a', second 10 class 'b'
class.key <- c(rep("a", 10), rep("b", 10));  # vector of class labels. Known in this case; often read from a text file.

# figure out which voxels to run in this chunk
num.vox <- nrow(lookup);
if (max(vol3d) != num.vox) { stop("max(vol3d) != num.vox"); }
all.chunks <- seq(from=1, to=num.vox, by=4500); # for 4 chunks.
if (do.chunk > length(all.chunks)) { stop("do.chunk > length(all.chunks)"); }
if (do.chunk == length(all.chunks)) { do.centers <- all.chunks[do.chunk]:num.vox; }
if (do.chunk < length(all.chunks)) { do.centers <- all.chunks[do.chunk]:(all.chunks[do.chunk+1] - 1); }
if (length(class.key) != dim(all.img)[4]) { stop("number of labels doesn't match number of images"); }

# now do the classification for all the searchlights in this chunk.
out.img <- array(NA, dim(vol3d));
for (v in do.centers) {  # v <- do.centers[10];
  if (v%%500 == 0) { print(paste("at", v, "of", length(do.centers))); }   # print a message to show progress
  # find which voxels belong in this center voxel's searchlight
  voxs <- union(v, unlist(lookup[v,], use.names=FALSE));   # surrounding voxels for this center; union to put center in the searchlight
  voxs <- voxs[which(!is.na(voxs))];  # get rid of NAs. will be NA entries if some surrounding voxels not in the brain.
  if (length(voxs) > smallest.searchlight) {    # how many surrounding voxels must be in the searchlight? Smaller ones (edge of brain) will be skipped.
    # put the data into a matrix so can classify
    thisData <- array(NA, c(length(class.key), length(voxs)));  # images in the rows (voxels in this searchlight only), voxels in the columns
    for (i in 1:length(voxs)) {   # i <- 1;
      coords <- which(vol3d == voxs[i], arr.ind=TRUE);
      if (ncol(coords) != 3 | nrow(coords) != 1) { stop("wrong sized coords"); }
      thisOne <- all.img[coords[1], coords[2], coords[3],]
      if (sd(thisOne) > 0) { thisData[,i] <- thisOne; } else { stop("zero variance voxel"); } 
    }
    #thisData <- t(scale(t(thisData)));  # slow row-scaling ...
    
    # do the cross-validation and classification.
    # This is a stupidly simple cross-validation: leave out first of each class, then second, etc. Not recommended for real analyses!
    a.inds <- which(class.key == "a");
    b.inds <- which(class.key == "b");
    if (length(a.inds) != length(b.inds)) { stop("imbalance in a and b"); }   # need same number of examples of each class
    accs <- rep(NA, length(a.inds));  # leave-one-pair-out cross-validation
    for (i in 1:length(a.inds)) {   # i <- 1;
      test.inds <- c(a.inds[i], b.inds[i]);
      train.inds <- (1:length(class.key))[-test.inds]
      
      if (length(union(train.inds, test.inds)) != nrow(thisData)) { stop("length(union(train.inds, test.inds)) != nrow(thisData)"); }
      if (length(train.inds) < length(test.inds)) { stop("length(train.inds) < length(test.inds)"); }  # generally want more training than testing examples
      
      # do the linear svm
      fit <- svm(class.key[train.inds]~., data=thisData[train.inds,], type="C-classification", kernel="linear", cost=1, scale=FALSE);  # train
      c.matrix <- table(class.key[test.inds], predict(fit, thisData[test.inds,]));   # test, getting the confusion matrix
      accs[i] <- sum(diag(c.matrix))/sum(c.matrix);  # calculate proportion correct
    }
    coords <- which(vol3d == v, arr.ind=TRUE);   # find the coordinates of this searchlight center
    if (ncol(coords) != 3 | nrow(coords) != 1) { stop("wrong sized coords"); }
    out.img[coords[1], coords[2], coords[3]] <- mean(accs);   # store the mean accuracy in the correct place
  }
}

# these NIfTI header parameters will match dataset.nii.gz and roi.nii.gz
out.nii <- nifti(out.img, datatype=64, dim=dim(out.img), srow_x=c(4,0,0,-80), srow_y=c(0,4,0,-112), srow_z=c(0,0,4,-52),
                         qoffset_x=-80, qoffset_y=-112, qoffset_z=-52, xyzt_units=2);
out.nii@sform_code <- 4; 
out.nii@qform_code <- 4;
pixdim(out.nii)[1:8] <- c(1, 4,4,4, 1, 0,0,0);  # qFactor, mm/voxel, number of timepoints
writeNIfTI(out.nii, paste0(out.path, "chunk", do.chunk, "_rad", radius, "_slOut"));    # save the accuracies as a NIfTI image

###################################################################################################################################################################
# combine the chunk files into a single whole-brain output file

library(oro.nifti); 
rm(list=ls()); 

in.path <- "d:/temp/searchlightDemo/withSpheres/";
out.path <- "d:/temp/searchlightDemo/withSpheres/";

radius <- 2;   # searchlight radius
num.chunks <- 4;   # lazy hard-coding
img.dim <- c(41,48,35);

out.fname <- paste0(out.path, "searchlightAccuracies_rad", radius);
if (!file.exists(paste0(out.fname, ".nii.gz"))) {   # check so don't overwrite existing files
  ctr <- 0;
  allImg <- array(0, img.dim); 
  for (i in 1:num.chunks) {    # i <- 1;   # first, check to see if all chunk files are present
    fname <- paste0(in.path, "chunk", i, "_rad", radius, "_slOut.nii.gz")
    if (!file.exists(fname)) { 
      print(paste("missing:", fname)); 
    } else {
      inimg <- readNIfTI(fname, reorient=FALSE);
      if (dim(inimg)[1] != img.dim[1] | dim(inimg)[2] != img.dim[2] | dim(inimg)[3] != img.dim[3]) { stop("not expected dims"); }
      inds <- which(inimg != 0, arr.ind=TRUE);  # get the nonzero voxel locations: make sure don't have overlaps before adding!
      if (max(allImg[inds]) != 0 | min(allImg[inds])!= 0) { stop("overlapping chunks"); }  # stop if overlaps.
      allImg[inds] <- allImg[inds] + inimg[inds]; 
      ctr <- ctr + 1;
    }
  }
  if (ctr == num.chunks) {
    # complete, so write it out as a whole-brain NIfTI
    out.img <- nifti(allImg, datatype=64, dim=dim(allImg), srow_x=c(4,0,0,-80), srow_y=c(0,4,0,-112), srow_z=c(0,0,4,-52),
                         qoffset_x=-80, qoffset_y=-112, qoffset_z=-52, xyzt_units=2);
    pixdim(out.img)[1:8] <- c(1, 4,4,4, 1, 0,0,0);  # qFactor, mm/voxel, number of timepoints
    out.img@sform_code <- 4; 
    out.img@qform_code <- 4;
    writeNIfTI(out.img, out.fname);
    # this code will delete the chunk files if made the output file
    if (file.exists(paste0(out.fname, ".nii.gz"))) { 
      for (i in 1:num.chunks) {    # i <- 1;
        fname <- paste0(in.path, "chunk", i, "_rad", radius, "_slOut.nii.gz")
        msg <- file.remove(fname)
        if (msg == FALSE) { stop(paste("didn't remove", fname)); }
      }
    }
  } 
}

###################################################################################################################################################################
# now, we have a single 3d nifti image in which each voxel value is the accuracy of the voxel's searchlight.
# this sample dataset is for a single participant; if there were multiple participants we'd do group-level statistics with these images.
###################################################################################################################################################################
























































#
