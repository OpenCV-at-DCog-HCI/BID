## My goal is to create more programs that mirror Olaf Sporns' computations
## Main program is called BID for Behavioral InfoDynamics (or is it Behavior?)
## Jeremy believes in goals
## Ed believes in beliefs

## Created by Jeremy Karnowski August 1, 2011
## Updated November 10, 2011
## This version updated by Edwin Hutchins December 31, 2011

from scipy import *

def BID_norm(Mdata,lo_limit,hi_limit):
    '''
    Normalizes the input data 'Mdata' between the user-specified limits.

    INPUT
        Mdata = MxN matrix
        lo_limit = scalar         lower limit
        hi_limit = scalar         upper limit

    OUTPUT
        Mdata = MxN matrix         normalized output
    '''

    #I'm still not sure how to throw an error for too few arguments
    #try:
    scale = Mdata.max(1) - Mdata.min(1)                            # max and min of each data stream
    for i in range(0,Mdata.shape[0]):                            # for each data stream, subtract smallest
        Mdata[i,:] = (Mdata[i,:] - Mdata[i,:].min())/scale[i]    # element and scale appropriately
    return (lo_limit + Mdata*(hi_limit-lo_limit))
    #except:
    #    print "Must specify data matrix and both lower and upper limits"

def BID_discrete(Mdata,nst):
    '''
    Discretizes the input data into discrete states

    INPUT
        Mdata = MxN matrix
        nst = number of states         resolution

    OUTPUT
        Mstates = MxN matrix        every row ranges from 0 to nst-1
    '''
    M = zeros(Mdata.shape)
    final = zeros(Mdata.shape)
    scale = Mdata.max(1) - Mdata.min(1)                            # max and min of each data stream
    for i in range(0,Mdata.shape[0]):                            # for each data stream, subtract smallest
        M[i,:] = (Mdata[i,:] - Mdata[i,:].min())/scale[i]    # element and scale appropriately

    bins = arange(0.0,1.0,1.0/nst)
    bins = append(bins,1+1e-6)

    #loop over edges of bins and find elements in those bins
    for lower_limit,upper_limit in zip(bins[:-1],bins[1:]):
        final[(M>=lower_limit) & (M<upper_limit)] = round(upper_limit*nst)
    return final


def BID_jointH(Mdata):
    '''
    Computes the joint entropy of M data streams of length N.
    Data must be a numpy matrix with M rows and N columns and must
    consist of only binned data
    '''

    #Mdata has M data streams, all of which have been broken up into
    #a certain number of bins. In Matlab code, if you have 10 data streams,
    #you create a 10 dimensional matrix with each dimension having the
    #size of the number of bins. This blows up with huge amounts of bins
    #and dimension. Instead, we can create those same indices, but use
    #them as keys in a dictionary and create a sparse representation.

    #Each key is multidimensional matrix index and each value is the count
    #of how many times that configuration has occurred in our M-dim data.

    #Create dictionary and find length of our time series
    jointH = {}

    M = Mdata.shape[0]
    try:  # to treat Mdata as a an n > 1 dimensional array
        N = Mdata.shape[1]
        for x in range(0,N):
            try: # to add one to an already existing dictionary entry
                jointH[''.join([str(j) for j in int_(Mdata[:,x]).tolist()])] += 1.0
            except: # set the value of the new dictionary entry to 1
                jointH[''.join([str(j) for j in int_(Mdata[:,x]).tolist()])] = 1.0
    except: #  Mdata is one dimensional. Simplify indexing of jointH dictionary.
        #print "Input matrix has only one dimension."
        N = M
        for x in range(0,N):
            try:
                jointH[str(int(Mdata[x]))] += 1.0
            except:
                jointH[str(int(Mdata[x]))] = 1.0

    #Divide each entry by number of samples to normalize probabilities
    for key in jointH:
        jointH[key] /= N

    #Create new dictionary with p(x)*log(p(x))
    jointH_final = {}
    for key in jointH:
        jointH_final[key] = jointH[key]*log(jointH[key])

    #Compute entropy by summing them up
    jointH_final_sum = -sum(jointH_final.values())

    return jointH_final_sum

##      The following function in Jeremy's code called BID_joint(..), which threw an error.  I replaced
##      BID_joint with BID_jointH. I also changed the Mdata row indices to 0 and 1 instead of 1 and 2.
##      EH 12/31/2011

def BID_MI(Mdata):
    '''
    Computes the mutual information of M=2 data streams of length N.
    Data must be a numpy matrix with M rows and N columns and must
    consist of only binned data
    '''
    return BID_jointH(Mdata[0,:]) + BID_jointH(Mdata[1,:]) - BID_jointH(Mdata)

##    The following function in Jeremy's code included a reshape method, Mdata[i,:].reshape(1,N).
##    This had the effect of embedding the Mdata[i,:] vector as the only row in an array. This
##    prevented the BID_jointH function from finding the data.

def BID_integration(Mdata):
    '''
    Computes the integration of M data streams of length N.
    Data must be a numpy matrix with M rows and N columns and must
    consist of only binned data
    '''
    M = Mdata.shape[0]
    N = Mdata.shape[1]
    #return sum([BID_jointH(Mdata[i,:].reshape(1,N)) for i in range(0,M)]) - BID_jointH(Mdata)
    return sum([BID_jointH(Mdata[i,:]) for i in range(0,M)]) - BID_jointH(Mdata)

def BID_complexity(Mdata):
    '''
    Computes the complexity of M data streams of length N.
    Data must be a numpy matrix with M rows and N columns and must
    consist of only binned data
    '''
    M = Mdata.shape[0]
    N = Mdata.shape[1]
    return sum([BID_jointH(vstack((Mdata[:i],Mdata[i+1:]))) for i in range(0,M)]) - (M-1)*BID_jointH(Mdata)