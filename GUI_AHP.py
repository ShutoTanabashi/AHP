from numpy import *
import tkinter
from tkinter import ttk
from matplotlib import *
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from random import *
#preace change font
font = {"family":"yumin"}
rc('font', **font)

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
        self.eveclist = self.eveclist.T
        for i in arange(self.n):

            self.eveclist[i] /= linalg.norm(self.eveclist[i])
            self.eveclist[i] /= self.eveclist[i].sum()
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
        self.caleig()

    def pristates(self):
        print(self.mA)
        print('ci = {}'.format(self.ci))

    def cons(self):
        """
        Calculate Consistency Index of matrix mA.

        Returns
        __________

        judge : bool
            If matrix mA is consistent, judge is True.
        """
        self.ci = (self.evalmax - self.n) / (self.n - 1)
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
        self.numhie  = number
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
            for j in range(self.numfuc[i - 1]):
                self.lismA[i].append(Mpc(self.numfuc[i]))

    def run(self, layer):
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
            for i in arange(1,self.numfuc[layer - 1]):
                matvec = vstack((matvec,self.lismA[layer][i].evecmax))
            return dot(self.run(layer - 1),matvec)

class Appahp(tkinter.Frame):
    """
    Making popup window.

    Attributes
    __________

    trghie : Mpc
        Target AHP system.

    fratxt : tkinter.Frame
        Window to input word.

    fraval : tkinter.Frame
        Window to input evaluetion values.

    frares : tkinter.Frame

    evatxt : list of str
        Words for buttons.

    evaval :  dictionary
        Correspondence table of button words and evaluetion values.

    evabtn : list of ttk.Button
        Buttons to decide evaluetion value.

    elelbl : list of ttk.Button
        Labels to show evaluetion standards or alternative proposals.

    numclick : int
        Number of times the button was pressed.

    """
    def __init__(self, trghie, master=None):
        super().__init__(master)
        self.pack()
        self.trghie = trghie
        self.evatxt = ['同程度', '少し重要', '重要', 'かなり重要', '絶対に重要']
        self.evaval = {self.evatxt[0]:1, self.evatxt[1]:2, self.evatxt[2]:3,\
            self.evatxt[3]:4, self.evatxt[4]:5}
        self.evabtn = []
        self.elelbl = []
        self.numclick = 0
        self.fraele = ttk.Frame(self)
        self.fraele.grid(row = 0, column = 0, sticky = "snew")
        self.create_widgets_fraele()
        self.fraval = ttk.Frame(self)
        self.fraval.grid(row = 0, column = 0, sticky = "snew")
        self.frares = tkinter.Frame(self)
        self.frares.grid(row = 0, column = 0, sticky = "snew")
        self.fraele.tkraise()

    def create_widgets_fraele(self):
        """
        Setting character entry box, buttons and labels to Frame "fraele".

        Parameters
        __________
        inpbox : ttk.Entry
            Entry box to add elements as evaluetion standards and
            alternative proposals.

        setbtn : ttk.Button
            Button to add new elements.

        bgnbtn : ttk.Button
            Button for finishing entry and moving to Frame"fraval".

        explbl : ttk.Label
            Label for showing explanatory text.

        ccklbl : ttk.Label
            Label for showing previous input text.

        nowreg : int
            Registration destination of input.
            If self.nowreg is 0, register entry to alternative proposal.
            If self.nowreg is more than 1, register entry to evaluetion standards.

        trglst : list of str
            List for entering text as evaluetion standards or
            alternative proposals.

        trgnum : int
            Number of entering as evaluetion standards or
            alternative proposals.
        """
        self.nowreg = 1
        self.trgnum = zeros(self.trghie.numhie, dtype = int)
        ttk.Style().configure("title.TLabel", foreground = '#228899', \
            anchor = "center")
        self.explbl = ttk.Label(self.fraele, text = '選択肢を入力してください',\
            width = 45, style = "title.TLabel")
        self.explbl.grid(column = 0, columnspan = 3, row = 0)

        self.ccklbl = ttk.Label(self.fraele, text = '', width = 45,\
            style = "title.TLabel")
        self.ccklbl.grid(column = 0, columnspan = 3, row = 2)

        self.inpbox = ttk.Entry(self.fraele ,width = 45)
        self.inpbox.bind("<Return>", self.addelement)
        self.inpbox.grid(column = 0, columnspan = 2, row = 1)

        ttk.Style().configure("setbtn.TButton", relief = 'flat',\
            background = '#116677', foreground = '#228899')
        self.setbtn = ttk.Button(self.fraele, text = "追加", width = 45,\
            style = "setbtn.TButton")
        self.setbtn.bind("<1>", self.addelement)
        self.setbtn.bind("<3>", self.showfuctor)
        self.setbtn.grid(column = 2, row = 1)

        ttk.Style().configure("bgnbtn.TButton", relief = 'flat',\
            background = '#771155', foreground = '#994477')
        self.bgnbtn = ttk.Button(self.fraele, text = "選択肢の追加を終了",\
            width = 45, style ="bgnbtn.TButton")
        self.bgnbtn.bind("<1>",self.changemode)
        self.bgnbtn.grid(column = 0, columnspan = 3, row = 3)

    def addelement(self, event):
        """
        Add new elements to list of evaluetion standards
        and alternative proposals.

        Parameters
        __________
        ent : str
            Text in the entry box.
        """
        self.ent = self.inpbox.get()
        if self.ent != '':
            self.trghie.addfuc(self.nowreg, self.ent)
            self.trgnum[self.nowreg] += 1
            self.ccklbl["text"] = self.ent + 'が追加されました'
            self.inpbox.delete(0, tkinter.END)

    def showfuctor(self, *event):
        """
        Show elements of target list.
        """
        print('Number of this hierarchy【{}】'.format(self.nowreg))
        print('Number of element in this hierarchy【{}】'.format(self.trgnum[self.nowreg]))
        j = 0
        for i in arange(self.trgnum[self.nowreg]):
            print(str(i) + ' : ' + self.trghie.fuctor[self.nowreg][i])

    def changemode(self, event):
        """
        Change input type or move to evaluetion.
        """
        self.showfuctor()
        if self.nowreg != 0:
            if self.nowreg == 1:
                self.explbl["text"] = '評価基準を入力してください'
            else:
                self.explbl["text"] = '{}層目の判断要素を入力してください'\
                    .format(self.nowreg + 1)
            self.nowreg -= 1
            self.ccklbl["text"] = ''
        else:
            self.trghie.makemat()
            self.create_widgets_fraval()
            self.fraval.tkraise()

    def create_widgets_fraval(self):
        """
        Setting buttons and texts to Frame "fraval".

        Parameters
        __________
        results : list of float
            Evaluetion results.

        button : ttk.Button
            Button for selecting evaluetion value

        titlbl : ttk.Label
            Label for explanatory text, evaluetion standards
            and alternative proposals.

        atelbl : ttk.Label
            Label for atention text.

        lbcol : str
            Text color for Label or Button

        tmp1 : str
            Temporary value for color calculation.

        tmp2 : str
            Temporary value for color calculation.

        tmpr : int
            Temporary value for label row.

        tmp : int
            Temporary value for color calculation.

        pos : str
            Character posision within label.

        eletop : int
            Number of evaluetion standards.

        ele1 : int
            Number of the comparing side.

        ele2 : int
            Number to be compered.

        txt : str
            elelbl text.
        """
        self.nowreg = 0
        self.results = []
        self.eletop = 0
        self.ele1 = len(self.trghie.fuctor[0]) - 1
        self.ele2 = self.ele1 - 1

        ttk.Style().configure("title.TLabel", foreground = '#228899',\
          anchor = "center")
        self.titlbl = ttk.Button(self.fraval, text = \
        '判断基準', width = 45, style = "title.TLabel")
        self.titlbl.grid(column = 3, columnspan = 3, row = 0)

        ttk.Style().configure("atention.TLabel", foreground = '#dd0011',\
            anchor = "center")
        self.atelbl = ttk.Button(self.fraval, text = '', width = 45, \
            style = "atention.TLabel")
        self.atelbl.grid(column = 3, columnspan = 3, row = 3)

        for i in range(2):
            if i==0:
                tmp1 = '00'
                tmp2 = '99'
                pos = "w"
                txt = self.trghie.fuctor[0][self.ele1]
            else:
                tmp1 = '99'
                tmp2 = '00'
                pos = "e"
                txt = self.trghie.fuctor[0][self.ele2]
            lbcol = '#33' + tmp1 + tmp2
            ttk.Style().configure(str(i) + "ele.TLabel",\
                background = lbcol, foreground = '#ffffff', anchor = pos)
            label = ttk.Label(self.fraval, text = txt,\
                width = 45, style = str(i) + "ele.TLabel")
            self.elelbl.append(label)
            tmpr = i * 6
            self.elelbl[i].grid(column = tmpr, columnspan = 3, row = 1)

        for i in arange(9):
            tmp = 8 - i
            lbcol = '#00{}{}{}{}'.format(i, i, tmp, tmp)
            ttk.Style().configure(str(i) + "eva.TButton",\
                relief = 'flat', background = lbcol, foreground = lbcol)
            button = ttk.Button(self.fraval, text = self.evatxt[abs(4-i)],\
                width = 15, style = str(i) + "eva.TButton")
            if i < 5:
                button.bind("<1>",self.judgeleft)
            else:
                button.bind("<1>",self.judgeright)
            self.evabtn.append(button)
            self.evabtn[i].grid(column = i, row = 2)

    def judgeleft(self, event):
        """
        Call judge from left buttons(When left element better than right one).

        Parameters
        __________

        val : float
            Evaluetion value.
        """
        val = self.judge(event)
        self.trghie.lismA[self.nowreg][self.eletop].setval(self.ele1, self.ele2, val)
        hoge = self.trghie.lismA[self.nowreg][self.eletop].cons()
        self.setnextelelbl()

    def judgeright(self, event):
        """
        Call judge from right buttons(When right element better than left one).

        Parameters
        __________
        val : float
            Evaluetion value.

        """
        val = self.judge(event)
        self.trghie.lismA[self.nowreg][self.eletop].setval(self.ele2, self.ele1, val)
        hoge = self.trghie.lismA[self.nowreg][self.eletop].cons()
        self.setnextelelbl()

    def judge(self, event):
        """
        Change the evaluetion value.
        """
        self.results.append(float(self.evaval[event.widget["text"]]))
        self.numclick += 1
        return float(self.evaval[event.widget["text"]])

    def setnextelelbl(self):
        if self.ele2 <= 0:
            self.ele1 -= 1
            if self.ele1 <= 0:
                self.trghie.lismA[self.nowreg][self.eletop].pristates()
                if self.trghie.lismA[self.nowreg][self.eletop].cons():
                    self.eletop -= 1
                    self.atelbl["text"] = ''
                    if self.eletop < 0:
                        if self.nowreg >= self.trghie.numhie - 1:
                            self.create_widgets_frares()
                            self.frares.tkraise()
                            #Avoid bug
                            self.nowreg = 0
                        else:
                            self.nowreg += 1
                            self.eletop = len(self.trghie.fuctor[self.nowreg - 1])\
                                - 1
                else:
                    self.atelbl["text"] = '判断の整合性が少ないです．もう一度行ってください'
                self.ele1 = len(self.trghie.fuctor[self.nowreg]) - 1
            self.ele2 = self.ele1
        self.ele2 -= 1
        if self.nowreg != 0:
            self.titlbl["text"] = self.trghie.fuctor[self.nowreg - 1][self.eletop]
        self.elelbl[0]["text"] = self.trghie.fuctor[self.nowreg][self.ele1]
        self.elelbl[1]["text"] = self.trghie.fuctor[self.nowreg][self.ele2]

    def create_widgets_frares(self):
        """
        Setting buttons and texts to Frame "fraval".

        Parameters
        __________

        reslbl : ttk.Label
            Label to desplay result.

        ahpres : ndarray
            Result value of AHP.

        tmpval : float
            Temporary variable to put result value for bubblesort.

        tmptxt : str
            Temporary variable to put elements for bubblesort.
        """
        self.ahpres = self.trghie.run(self.trghie.numhie - 1)
        for i in range(self.trghie.numfuc[self.trghie.numhie - 1]):
            print(self.trghie.fuctor[self.trghie.numhie - 1][i] )
            print(self.ahpres[i])
        #Bubblesort
        for i in range(len(self.trghie.fuctor[self.trghie.numhie - 1])-1, 0, -1):
            for j in range(i):
                if self.ahpres[j] < self.ahpres[j+1]:
                    tmpval = self.ahpres[j+1]
                    tmptxt = self.trghie.fuctor[self.trghie.numhie - 1][j+1]
                    self.ahpres[j+1] = self.ahpres[j]
                    self.trghie.fuctor[self.trghie.numhie - 1][j+1] =\
                        self.trghie.fuctor[self.trghie.numhie - 1][j]
                    self.ahpres[j] = tmpval
                    self.trghie.fuctor[self.trghie.numhie - 1][j] = tmptxt

        ttk.Style().configure("resl.TLabel", background = '#5ff')
        self.reslbl = ttk.Label(self.frares, text = "Best choice : {}"\
            .format(self.trghie.fuctor[self.trghie.numhie - 1][0]),\
            style = "resl.TLabel")
        self.reslbl.pack()

        fig = Figure(figsize=(6,4), dpi=100)
        ax = fig.add_subplot(111)
        ax.bar(range(len(self.trghie.fuctor[self.trghie.numhie - 1])),\
            self.ahpres,\
            tick_label = self.trghie.fuctor[self.trghie.numhie - 1],\
            color = "#228899", linewidth = 0, align = "center")
        ax.set_ylabel("Priority [%]")
        ax.grid(True)
        canvas = FigureCanvasTkAgg(fig, master=self.frares)
        canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH,\
            expand=1)
        canvas._tkcanvas.pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
        print('Preace close window')


if __name__ == '__main__':
    mat_ahp = Hierarchy()
    root = tkinter.Tk()
    root.title('Analytic Hierarchy Process')
    root.geometry("950x450+100+100")
    app = Appahp(mat_ahp, master=root)
    app.mainloop()
