## My goal is to create more programs that mirror Olaf Sporns' computations
## Main program is called BID for Behavioral InfoDynamics (or is it Behavior?)
## Jeremy believes in goals
## Ed believes in beliefs

## Created by Jeremy Karnowski August 1, 2011
## Updated by Jeremy Karnowski November 10, 2011
## Updated by Edwin Hutchins December 31, 2011
## Current version by Jeremy Karnowski June 9, 2012

from scipy import *

def BID_norm(Mdata,lo_limit,hi_limit):
    """
    Normalizes the input data 'Mdata' between the user-specified limits.

    INPUT
        Mdata = MxN matrix
        lo_limit = scalar         lower limit
        hi_limit = scalar         upper limit

    OUTPUT
        Mdata = MxN matrix         normalized output
    """
    # Possible error message to include later:
    # print "Must specify data matrix and both lower and upper limits"
    scale = Mdata.max(1) - Mdata.min(1)                          # max and min of each data stream
    for i in range(0,Mdata.shape[0]):                            # for each data stream, subtract smallest
        Mdata[i,:] = (Mdata[i,:] - Mdata[i,:].min())/scale[i]    # element and scale appropriately
    return lo_limit + Mdata*(hi_limit-lo_limit)


def BID_discrete(Mdata,nst):
    """
    Discretizes the input data into discrete states

    INPUT
        Mdata = MxN matrix
        nst = number of states      resolution

    OUTPUT
        Mstates = MxN matrix        every row ranges from 0 to nst-1
    """
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


def BID_jointH(Mdata,nst):
    """
    Computes the joint entropy of M data streams of length N.
    Data must be a numpy matrix with M rows and N columns and must
    consist of only binned data
    """

    #Mdata has M data streams, all of which have been broken up into
    #a certain number of bins. In Matlab code, if you have 10 data streams,
    #you create a 10 dimensional matrix with each dimension having the
    #size of the number of bins. This blows up with huge amounts of bins
    #and dimension. Instead, we can create those same indices, but use
    #them as keys in a dictionary and create a sparse representation.

    #Each key is multidimensional matrix index and each value is the count
    #of how many times that configuration has occurred in our M-dim data.

    #To fix data vectors not properly formatted
    try:
        N = Mdata.shape[1]
    except:
        Mdata = Mdata.reshape(1,Mdata.shape[0])

    #Create dictionary and find length of our time series
    jointH = {}
    M = Mdata.shape[0]
    N = Mdata.shape[1]

    if N < 3*(nst**M):
        print "Warning: Number of samples < 3*number of states. Potentially not meaningful estimate!"

    # If M = 1, Mdata is one dimensional. This means just compute entropy.
    if M==1:
        for x in range(0,N):
            try: # to add one to an already existing dictionary entry
                jointH[str(int(Mdata[:,x]))] += 1.0
            except: # set the value of the new dictionary entry to 1
                jointH[str(int(Mdata[:,x]))] = 1.0
    # M > 1, so compute joint entropy.
    else:
        for x in range(0,N):
            try: # to add one to an already existing dictionary entry
                jointH[tuple(int_(Mdata[:,x]).tolist())] += 1.0
            except: # set the value of the new dictionary entry to 1
                jointH[tuple(int_(Mdata[:,x]).tolist())] = 1.0

    # Divide each entry by number of samples to normalize probabilities
    for key in jointH:
        jointH[key] /= N

    # Create new dictionary with p(x)*log(p(x))
    jointH_final = {}
    for key in jointH:
        jointH_final[key] = jointH[key]*log2(jointH[key])

    # Compute entropy by summing them up
    jointH_obs = -sum(jointH_final.values())

    # Compute H_inf based on Roulston (1999)
    B_star = size(jointH_final.values())
    jointH_inf = jointH_obs + (B_star - 1)/(2*N)

    # Compute standard deviation of entropy based on Roulston (1999)
    jointSigma_Htemp = {}
    for key in jointH:
        jointSigma_Htemp[key] = ((log2(jointH[key]) + jointH_obs)**2)*jointH[key]*(1-jointH[key])
    jointSigma_H = sqrt(sum(jointSigma_Htemp.values())/N)

    return jointH_obs, jointH_inf, jointSigma_H


def BID_MI(Mdata,nst):
    """
    Computes the mutual information of M=2 data streams of length N.
    Data must be a numpy matrix with M rows and N columns and must
    consist of only binned data
    """
    if Mdata.shape[0] == 2:
        N = Mdata.shape[1]
        Mdata0 = Mdata[0,:].reshape(1,N)
        Mdata1 = Mdata[1,:].reshape(1,N)
        return BID_jointH(Mdata0,nst)[1] + BID_jointH(Mdata1,nst)[1] - BID_jointH(Mdata,nst)[1]


def BID_integration(Mdata,nst):
    """
    Computes the integration of M data streams of length N.
    Data must be a numpy matrix with M rows and N columns and must
    consist of only binned data
    """
    M = Mdata.shape[0]
    N = Mdata.shape[1]
    return sum([BID_jointH(Mdata[i,:].reshape(1,N),nst)[1] for i in range(0,M)]) - BID_jointH(Mdata,nst)[1]
    #return sum([BID_jointH(Mdata[i,:]) for i in range(0,M)]) - BID_jointH(Mdata)  # an old fix from Ed

def BID_complexity(Mdata,nst):
    """
    Computes the complexity of M data streams of length N.
    Data must be a numpy matrix with M rows and N columns and must
    consist of only binned data
    """
    # The below uses a shortcut based on the definition of joint entropy
    # H(Y|X) = H(X,Y) - H(X)
    M = Mdata.shape[0]
    N = Mdata.shape[1]
    return sum([BID_jointH(vstack((Mdata[:i],Mdata[i+1:])),nst)[1] for i in range(0,M)]) - (M-1)*BID_jointH(Mdata,nst)[1]
