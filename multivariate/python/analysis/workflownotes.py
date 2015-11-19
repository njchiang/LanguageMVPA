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