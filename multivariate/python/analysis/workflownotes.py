Current workflow:

#general workflow
#read in attributes (labels and chunks)
attr=SampleAttributes('verb_attributes.txt')

#read in dataset with mask
fds=fmri_dataset(samples='LMVPA001_LSA_betaSeries.nii.gz', targets=attr.targets, chunks=attr.chunks, mask='LMVPA001_grayMatter.nii.gz')
fds=fmri_dataset(samples='LMVPA001_LSA_betaSeries.nii.gz', targets=attr.targets, chunks=attr.chunks, mask='LMVPA001_left_IFG_operc.nii.gz')

zscore(fds, chunks_attr=None)
#define classifier
clf=kNN(k=1, dfx=one_minus_correlation, voting='majority')
clf=LinearCSVMC()
#train
clf.train(fds)

#predict
clf.predict(fds.samples)

#to check accuracy.
np.mean(predictions== fds.sa.targets)

# disable post-processing again
clf.set_postproc(None)
# dataset generator
hpart = HalfPartitioner(attr='chunks')
# complete cross-validation facility
# cv is a class that you call on the dataset
cv = CrossValidation(clf, hpart) 

# two fold validation
cv = CrossValidation(clf, HalfPartitioner(attr='chunks'), errorfx=lambda p, t: np.mean(p == t))

#n fold validation
cv = CrossValidation(clf, NFoldPartitioner(), errorfx=lambda p, t: np.mean(p == t))

# print confusion matrix
cv = CrossValidation(clf, NFoldPartitioner(), errorfx=lambda p, t: np.mean(p == t), enable_ca=['stats'])

#do the cross validation
cv_results=cv(fds)

#check accuracy
np.mean(cv_results)
cv_results.samples
np.round(cvte.ca.stats.stats['ACC%'], 1)



#check confusing matrix
print cv.ca.stats.as_string(description=True)
print cv.ca.stats.matrix


####
#sphereical searchlight
sl = sphere_searchlight(cv, radius=6, postproc=mean_sample())

 #apply searchlight

 res=sl(fds)

#compare with error
sphere_accs = res.samples[0]
res_mean = np.mean(res)
res_std = np.std(res)
# we deal with errors here, hence 1.0 minus
err_chance_level = 1-(1.0 / len(fds.uniquetargets))
chance_level = (1.0 / len(fds.uniquetargets))

frac_higher = np.round(np.mean(sphere_accs > chance_level + 2 * res_std), 3)

map2nifti(fds, sphere_accs).to_filename('sl.nii.gz'))

»    ds = dataset.copy(deep=False,
                       sa=['targets', 'chunks'],
                       fa=['voxel_indices'],
                       a=['mapper'])

    map2nifti(ds, sphere_accs).to_filename(os.path.join(outPath, sub + '_' + mask + '_anim_' + con + '_' + t + '_cvsl.nii.gz'))
map2nifti(res, imghdr=fullSet.a.imghdr).to


##### 
#feature selection
fsel = SensitivityBasedFeatureSelection(OneWayAnova(), FixedNElementTailSelector(500, mode='select', tail='upper'))
fclf = FeatureSelectionClassifier(clf, fsel)
cv = CrossValidation(fclf, NFoldPartitioner(),  enable_ca=['stats'])

results=cv(fds)

print np.round(cv.ca.stats.stats['ACC%'],1)

#!/usr/bin/python

import threading
import time

exitFlag = 0

class myThread (threading.Thread):
    def __init__(self, threadID, name, counter):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = name
        self.counter = counter
    def run(self):
        print "Starting " + self.name
        print_time(self.name, self.counter, 5)
        print "Exiting " + self.name

def print_time(threadName, delay, counter):
    while counter:
        if exitFlag:
            threadName.exit()
        time.sleep(delay)
        print "%s: %s" % (threadName, time.ctime(time.time()))
        counter -= 1

# Create new threads
thread1 = myThread(1, "Thread-1", 1)
thread2 = myThread(2, "Thread-2", 2)

# Start new Threads
thread1.start()
thread2.start()

print "Exiting Main Thread"


# MULTISUBJECT

from mvpa2.suite import *

verbose.level = 2

verbose(1, "Loading data...")
filepath = os.path.join(cfg.get('location', 'tutorial data'),
                        'hyperalignment_tutorial_data.hdf5.gz')
ds_all = h5load(filepath)
# zscore all datasets individually
_ = [zscore(ds) for ds in ds_all]
# inject the subject ID into all datasets
for i,sd in enumerate(ds_all):
    sd.sa['subject'] = np.repeat(i, len(sd))
# number of subjects
nsubjs = len(ds_all)
# number of categories
ncats = len(ds_all[0].UT)
# number of run
nruns = len(ds_all[0].UC)
verbose(2, "%d subjects" % len(ds_all))
verbose(2, "Per-subject dataset: %i samples with %i features" % ds_all[0].shape)
verbose(2, "Stimulus categories: %s" % ', '.join(ds_all[0].UT))

# use same classifier
clf = LinearCSVMC()

# feature selection helpers
nf = 100
fselector = FixedNElementTailSelector(nf, tail='upper',
                                      mode='select',sort=False)
sbfs = SensitivityBasedFeatureSelection(OneWayAnova(), fselector,
                                        enable_ca=['sensitivities'])
# create classifier with automatic feature selection
fsclf = FeatureSelectionClassifier(clf, sbfs)

verbose(1, "Performing classification analyses...")
verbose(2, "within-subject...", cr=False, lf=False)
wsc_start_time = time.time()


# store results in a sequence
wsc_results = [cv(sd) for sd in ds_all]
wsc_results = vstack(wsc_results)
verbose(2, " done in %.1f seconds" % (time.time() - wsc_start_time,))



# leave-one-run-out for hyperalignment training
for test_run in range(nruns):
    # split in training and testing set
    ds_train = [sd[sd.sa.chunks != test_run,:] for sd in ds_all]
    ds_test = [sd[sd.sa.chunks == test_run,:] for sd in ds_all]

    # manual feature selection for every individual dataset in the list
    anova = OneWayAnova()
    fscores = [anova(sd) for sd in ds_train]
    featsels = [StaticFeatureSelection(fselector(fscore)) for fscore in fscores]
    ds_train_fs = [featsels[i].forward(sd) for i, sd in enumerate(ds_train)]


    # Perform hyperalignment on the training data with default parameters.
    # Computing hyperalignment parameters is as simple as calling the
    # hyperalignment object with a list of datasets. All datasets must have the
    # same number of samples and time-locked responses are assumed.
    # Hyperalignment returns a list of mappers corresponding to subjects in the
    # same order as the list of datasets we passed in.


    hyper = Hyperalignment()
    hypmaps = hyper(ds_train_fs)

    # Applying hyperalignment parameters is similar to applying any mapper in
    # PyMVPA. We start by selecting the voxels that we used to derive the
    # hyperalignment parameters. And then apply the hyperalignment parameters
    # by running the test dataset through the forward() function of the mapper.

    ds_test_fs = [featsels[i].forward(sd) for i, sd in enumerate(ds_test)]
    ds_hyper = [ hypmaps[i].forward(sd) for i, sd in enumerate(ds_test_fs)]

    # Now, we have a list of datasets with feature correspondence in a common
    # space derived from the training data. Just as in the between-subject
    # analyses of anatomically aligned data we can stack them all up and run the
    # crossvalidation analysis.

    ds_hyper = vstack(ds_hyper)
    # zscore each subject individually after transformation for optimal
    # performance
    zscore(ds_hyper, chunks_attr='subject')
    res_cv = cv(ds_hyper)
    bsc_hyper_results.append(res_cv)

bsc_hyper_results = hstack(bsc_hyper_results)
verbose(2, "done in %.1f seconds" % (time.time() - hyper_start_time,))

verbose(1, "Average classification accuracies:")
verbose(2, "within-subject: %.2f +/-%.3f"
        % (np.mean(wsc_results),
           np.std(wsc_results) / np.sqrt(nsubjs - 1)))
verbose(2, "between-subject (anatomically aligned): %.2f +/-%.3f"
        % (np.mean(bsc_mni_results),
           np.std(np.mean(bsc_mni_results, axis=1)) / np.sqrt(nsubjs - 1)))
verbose(2, "between-subject (hyperaligned): %.2f +/-%.3f" \
        % (np.mean(bsc_hyper_results),
           np.std(np.mean(bsc_hyper_results, axis=1)) / np.sqrt(nsubjs - 1)))