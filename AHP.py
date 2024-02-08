from numpy import *

class Mpc:
    """
    Matrix for pairwise comparsion (for Analytic Hierarchy Process).

    Attributes
    __________

    mA : ndarray(n,n)
        Matrix for pairwise comparsion.
        Square matrix with n elements.

    n : int
        Number of elements of matrix mA.

    evallist : ndarray(n)
        Eigenvalue of mA.
        Array with n elements.

    eveclist : ndarray(n,n)
        Eigenvector of mA.
        Have a total of n vectors with n elements.

    evalmax :  float64
        Maximum eigenvalue of mA.

    evecmax : ndarray
        Eigenvector for maximum eigenvalue.

    ci : float64
        CI(Consistency Index) of matrix mA.
        The matrix mA is required to have a CI of 0.1 or less.
    """

    def __init__(self, numitem):
        """
        Parameters
        __________

        numitem : int
            Number of items in the same hierarchy.
        """
        self.n = numitem
        self.mA = identity(self.n)
        self.caleig()
        self.ci = (self.evalmax - self.n) / (self.n - 1)

    def caleig(self):
        """
        Calculate maximum eigenvalue and eigenvector for maximum eigenvalue.
        """
        self.evallist,self.eveclist = linalg.eig(self.mA)
        self.evalmax = 0
        for i in arange(self.n):
            if self.evalmax < self.evallist[i]:
                self.evalmax = self.evallist[i]
                self.evecmax = self.eveclist[i]

    def setval(self,i,j,x):
        """
        Change the value of i-th row and j-th column of matrix mA to x.

        Parameters
        __________

        i : int
            Index of row of matrix mA.

        j : int
            Index of column of matrix mA.

        x : float64
            Value of i-th row and j-th column of matrix mA.
        """
        self.mA[i][j] = x
        self.mA[j][i] = 1 / x
        self.mA.caleig()

    def cons(self):
        """
        Calculate Consistency Index of matrix mA.

        Returns
        __________

        judge : bool
            If matrix mA is consistent, judge is True.
        """
        self.ci = (self.evalmax - self.n) / (self.n - 1)
        print(self.mA)
        print('ci = {}'.format(self.ci))
        if self.ci < 0.1:
            judge = True
        else:
            judge = False
        return judge

class Hierarchy:
    """
    Hierarchy data for Analytic Hierarchy Process.
    Instances of this class have evaluetion standards, alternative proposals
    and matrix of pairwise comparsion.

    Attributes
    __________

    self.numhie : int
        Number of hierarchy layers.

    fuctor : 2D list of str
        Name of evaluetion standards or alternative proposals.

    numfuc : ndarray
        Number of evaluetion standards of alternative proposals in each hierarchy.

    lismA : 2D list of Mpc
        List of matrix for pairwise comparsion.

    Notes
    __________

    The smaller self.numhie is the higher hierarchy.
    fuctor [self.numhie][] is alternative proposals.
    """
    def __init__(self, number = 2):
        """
        Parameters
        __________

        number : int
            Number of hierarchy layers.
        """
        self.numhie  = number - 1
        self.fuctor = [[] for i in arange(self.numhie)]
        self.numfuc = zeros(self.numhie, dtype = int)
        self.lismA = [[] for i in arange(self.numhie)]

    def addfuc(self, layer, name):
        """
        Add evaluetion standards of alternative proposals.

        Parameters
        __________

        layer : int
            Hierarchy of elements to be added.

        name : str
            Name of evaluetion standards of alternative proposals.
        """
        self.fuctor[layer].append(name)
        self.numfuc[layer] = self.numfuc[layer] + 1

    def makemat(self):
        """
        Making Matrix for pairwise comparsion.
        """
        self.lismA[0].append(Mpc(self.numfuc[0]))
        for i in range(1, self.numhie):
            for j in range(self.numfuc[i]):
                self.lismA[i].append(Mpc(self.numfuc[j]))

    def run(self):
        """
        Run Analytic Hierarchy Process.

        Parameters
        __________
        layer : int
            Number of hierarchy which is calculated importance.

        matev : ndarray
            The set of eigenvector for maximum eigenvalue of that hierarchy.

        Returns
        __________
        importance : ndarray
            Importance calculated based on Analytic Hierarchy Process.
            If layer equal 0, return eigenvector for maximum eigenvalue of matrix mA.

        Notes
        __________
        This is recursive call.
        """
        if layer == 0:
            return self.lismA[0][0].evecmax
        else:
            matvec = self.lismA[layer][0].evecmax
            for i in arange(1,self.numfuc[layer]):
                matvec = vstack((matvec,self.lismA[layer][i].evecmax))
            return dot(self.run(layer - 1),matvec)

if __name__ == '__main__':
    print('Pleace execute "GUI_AHP"')
