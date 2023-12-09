import numpy as np
from sklearn.cross_decomposition import PLSRegression
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import ShuffleSplit, KFold
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import cross_val_score
from sklearn.utils import shuffle

from tqdm import tqdm
from numpy import dot, zeros, cumsum, diag, mean, sqrt
from numpy.linalg import inv, norm
from numpy.linalg import matrix_rank as rank

# def function of pls without iteration (SIMPLS)
def simpls_algorithm(X, y, nlv):
    '''
    Partial Least Squares, SIMPLS
    Ref https://github.com/freesiemens/SpectralMultivariateCalibration/blob/master/pypls.py
    :param X:
    :param y:
    :param nlv:
    :return:
    '''
    n_samples, n_variables = X.shape
    if np.ndim(y) == 1:
        y = y[:, np.newaxis]
    if n_samples != y.shape[0]:
        raise ValueError('The number of independent and dependent variable are inconsistent')

    nlv = np.min((nlv, n_samples, n_variables))
    V = zeros((n_variables, nlv))
    x_scores = zeros((n_samples, nlv))  # X scores (standardized)
    x_weights = zeros((n_variables, nlv))  # X weights
    x_loadings = zeros((n_variables, nlv))  # X loadings
    y_loadings = zeros((1, nlv))  # Y loadings
    y_scores = zeros((n_samples, nlv))  # Y scores
    s = dot(X.T, y).ravel()  # cross-product matrix between the X and y_data
    for i in range(nlv):
        r = s
        t = dot(X, r)
        tt = norm(t)
        t = t / tt
        r = r / tt
        p = dot(X.T, t)
        q = dot(y.T, t)
        u = dot(y, q)
        v = p  # P的正交基
        if i > 0:
            v = v - dot(V, dot(V.T, p))  # Gram-Schimidt orthogonal
            u = u - dot(x_scores, dot(x_scores.T, u))
        v = v / norm(v)
        s = s - dot(v, dot(v.T, s))
        x_weights[:, i] = r
        x_scores[:, i] = t
        x_loadings[:, i] = p
        y_loadings[:, i] = q
        y_scores[:, i] = u
        V[:, i] = v
    B = cumsum(dot(x_weights, diag(y_loadings.ravel())), axis=1)
    return {'B': B, 'x_scores': x_scores, 'x_loadings': x_loadings, 'y_loadings': y_loadings, \
            'x_scores_weights': x_weights, 'x_weights': x_weights, 'y_scores':y_scores}

class pls:
    def __init__(self, nlv=2):
        self.nlv = nlv
    
    def fit(self, X, y):
        self.X = X
        self.y = y
        meanX = mean(X, axis=0)
        meany = mean(y)
        Xcentered = X - meanX
        ycentered = y - meany
        model = simpls_algorithm(Xcentered, ycentered, self.nlv)
        meanX_hat = -1 * dot(meanX, model['B']) + meany
        model['B'] = np.append(meanX_hat[np.newaxis, :], model['B'],axis=0)
        self.model = model
        return self
       
    def predict(self, Xnew):
        B = self.model['B']
        if Xnew.shape[1] != B.shape[0]-1:
            raise ValueError('The feature number of predictor is isconsistent with that of indepnentent.')
        Xnew = np.append(np.ones([Xnew.shape[0],1]), Xnew, axis=1)
        ynew_hat = dot(Xnew, B)[:, -1]
        return ynew_hat
    
    def cross_val_predict(self, X, y, cv):
        yhat = np.zeros(y.shape)
        plsmodel = pls(nlv = self.nlv)
        for train, test in KFold(n_splits=cv).split(X):
            plsmodel.fit(X[train,:], y[train])
            yhat[test] = plsmodel.predict(X[test,:])
        return yhat
    
        
    def project(self, Xnew):
        meanX = mean(self.X, axis = 0)
        Xnew_c = Xnew - meanX
        Tnew = dot(Xnew_c,self.model['x_weights'])
        return Tnew
    
    def calc_VIP(self):
        # Reference: https://www.sciencedirect.com/topics/engineering/variable-importance-in-projection
        b = self.model['B'][1:, -1]
        T = self.model['x_scores']
        W = self.model['x_weights']
        up = 0
        low = 0
        for i in range(T.shape[1]):
            up += b**2 * dot(T[:, i], T[:, i]) * ((W[:, i])/norm(W[:, i]))**2
            low += b**2 * dot(T[:, i], T[:, i])
        VIP = sqrt(len(b)*up/low)
        return VIP
            
        
# Single step feature selection method
class MCUVE:
    def __init__(self, x, y, ncomp=1, nrep=500, testSize=0.2):
        self.x = x
        self.y = y
        # The number of latent components should not be larger than any dimension size of independent matrix
        self.ncomp = min([ncomp, rank(x)])
        self.nrep = nrep
        self.testSize = testSize
        self.criteria = None
        self.featureIndex = None
        self.featureR2 = np.full(self.x.shape[1], np.nan)
        self.selFeature = None
        self.tolerence = 1e-6

    def calcCriteria(self):
        PLSCoef = np.zeros((self.nrep, self.x.shape[1]))
        ss = ShuffleSplit(n_splits=self.nrep, test_size=self.testSize)
        step = 0
        with tqdm(total = self.nrep, desc='MCUVE iter') as pb:
            for train, test in ss.split(self.x, self.y):
                xtrain = self.x[train, :]
                ytrain = self.y[train]
                nlvi = min([self.ncomp, rank(xtrain)])
                plsModel = pls(nlvi)
                plsModel.fit(xtrain, ytrain)
                PLSCoef[step, :] = plsModel.model['B'][1:,-1]
                step += 1
                pb.update()
        PLSCoef = abs(PLSCoef)
        meanCoef = np.mean(PLSCoef, axis=0)
        stdCoef = np.std(PLSCoef, axis=0)
        self.criteria = meanCoef / stdCoef
        for i in range(len(self.criteria)):
            if meanCoef[i] + stdCoef[i] <= 1e-6 * mean(meanCoef):
                self.criteria[i] = 0
                

    def evalCriteria(self, cv=3):
        self.featureIndex = np.argsort(-np.abs(self.criteria))
        for i in range(self.x.shape[1]):
            xi = self.x[:, self.featureIndex[:i + 1]]
            if i<self.ncomp:
                regModel = LinearRegression()
            else:
                regModel = PLSRegression(min([self.ncomp, rank(xi)]))

            cvScore = cross_val_score(regModel, xi, self.y, cv=cv, n_jobs = -1)
            self.featureR2[i] = np.mean(cvScore)

    def cutFeature(self, *args):
        cuti = np.argmax(self.featureR2)
        self.selFeature = self.featureIndex[:cuti+1]
        if len(args) != 0:
            returnx = list(args)
            i = 0
            for argi in args:
                if argi.shape[1] == self.x.shape[1]:
                    returnx[i] = argi[:, self.selFeature]
                i += 1
        return tuple(returnx)


class RT(MCUVE):
    def calcCriteria(self):
        # calculate normal pls regression coefficient
        plsModel0 = pls(self.ncomp)
        plsModel0.fit(self.x, self.y)
        # calculate noise reference regression coefficient
        plsCoef0=plsModel0.model['B'][1:,-1]
        PLSCoef = np.zeros((self.nrep, self.x.shape[1]))
        print('-------Calculating RT criteria-----------')
        with tqdm(total=self.nrep, desc='RT iter') as pd:
            for i in range(self.nrep):
                randomidx = list(range(self.x.shape[0]))
                np.random.shuffle(randomidx)
                ytrain = self.y[randomidx]
                plsModel = pls(self.ncomp)
                plsModel.fit(self.x, ytrain)
                PLSCoef[i, :] = plsModel.model['B'][1:,-1]
                pd.update()
        plsCoef0 = np.tile(np.reshape(plsCoef0, [1, -1]), [self.nrep, 1])
        criteria = np.sum(np.abs(PLSCoef) > np.abs(plsCoef0), axis=0)/self.nrep
        self.criteria = criteria

    def evalCriteria(self, cv=3):
        # Note: small P value indicating important feature
        self.featureIndex = np.argsort(self.criteria)
        print('-------Evaluating RT criteria-----------')
        with tqdm(total = self.x.shape[1], desc='RT2 iter') as pd:
            for i in range(self.x.shape[1]):
                xi = self.x[:, self.featureIndex[:i + 1]]
                if i<self.ncomp:
                    regModel = LinearRegression()
                else:
                    regModel = PLSRegression(min([self.ncomp, rank(xi)]))
                cvScore = cross_val_score(regModel, xi, self.y, cv=cv)
                self.featureR2[i] = np.mean(cvScore)
                pd.update()


class VC(RT):
    def calcCriteria(self, cv=3):
        # calculate normal pls regression coefficient
        nVar = self.x.shape[1]
        sampleMatrix = np.ndarray([self.nrep,self.x.shape[1]], dtype=int)
        sampleMatrix[:, :] = 0
        errVector = np.ndarray([self.nrep,1])
        # The number of variable in combination should less than the total variable number
        if nVar > self.ncomp:
            nSample = max([self.ncomp, nVar//10])
        else:
            nSample = max([1, nVar-1])
        sampleidx = range(self.x.shape[1])
        for i in range(self.nrep):
            sampleidx = shuffle(sampleidx)
            seli = sampleidx[:nSample]
            plsModel = PLSRegression(n_components=min([self.ncomp, rank(self.x[:, seli])]))
            plsModel.fit(self.x[:, seli], self.y)
            sampleMatrix[i, seli] = 1
            yhati=cross_val_predict(plsModel, self.x[:, seli], self.y, cv=cv)
            errVector[i] = np.sqrt(mean_squared_error(yhati, self.y))
        plsModel = PLSRegression(n_components=self.ncomp)
        plsModel.fit(sampleMatrix, errVector)
        self.criteria = plsModel.coef_.ravel()

# Recursive feature selection method
class MSVC:
    def __init__(self, x, y, ncomp=1, nrep=7000, ncut=50, testSize=0.2):
        self.x = x
        self.y = y
        # The number of latent components should not be larger than any dimension size of independent matrix
        self.ncomp = min([ncomp, rank(x)])
        self.nrep = nrep
        self.ncut = ncut
        self.testSize = testSize
        self.criteria = np.full([ncut, self.x.shape[1]], np.nan)
        self.featureR2 = np.empty(ncut)
        self.selFeature = None

    def calcCriteria(self):
        varidx = np.array(range(self.x.shape[1]))
        ncuti = np.logspace(np.log(self.x.shape[1]), np.log(1), self.ncut, base=np.e)
        ncuti = (np.round(ncuti)).astype(int)
        for i in range(self.ncut):
            vcModel = VC(self.x[:, varidx], self.y, self.ncomp, nrep=self.nrep)
            vcModel.calcCriteria()
            self.criteria[i, varidx] = vcModel.criteria
            var_ranked = np.argsort(vcModel.criteria)
            if i < self.ncut - 1:
                varidx = varidx[var_ranked[:ncuti[i+1]]]

    def evalCriteria(self, cv=3):
        for i in range(self.ncut):
            varSeli = ~np.isnan(self.criteria[i, :])
            xi = self.x[:,varSeli]
            if sum(varSeli) < self.ncomp:
                regModel = LinearRegression()
            else:
                regModel = PLSRegression(min([self.ncomp, rank(xi)]))
            cvScore = cross_val_score(regModel, xi, self.y, cv=cv)
            self.featureR2[i] = np.mean(cvScore)

    def cutFeature(self, *args):
        cuti = np.argmax(self.featureR2)
        self.selFeature = ~np.isnan(self.criteria[cuti, :])
        if len(args) != 0:
            returnx = list(args)
            i = 0
            for argi in args:
                if argi.shape[1] == self.x.shape[1]:
                    returnx[i] = argi[:, self.selFeature]
                i += 1
        return tuple(returnx)


if __name__ == "__main__":
    print("This is the Feature selection library")