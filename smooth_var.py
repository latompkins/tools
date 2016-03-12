import ROOT,array

##
##       *  Smoothing of a histogram using a Gaussian kernel
##       *  Author: Dag Gillberg, with input from
##       *  Bogdan Malaescu and Caterina Doglioni
##       *  Python-ized by Lauren Tompkins
##


class histSmoothing:
    def __init__(self,rel_width=0.08, Nsteps=200, Nsigma=5.0):
        self.m_width=rel_width
        self.m_Nsteps=Nsteps
        self.m_Nsigma=Nsigma
        self.m_interpolate = False
        self.m_useError=False
        self.m_userMin= -99
        self.m_userMax= -99
        self.m_Nbins= -99
        self.m_min=-99
        self.m_binning=[]

        return
    #Parameters:
    #  - Relative widht of Gaussina kernel: DeltaPt/Pt
    #  - how many steps to use in numberical integration
    #  - how wide range to integrate over

    def UseInterpolation(self, doIt=True):
        self.m_interpolate=doIt
        return
    
    def UseBinUncertainty(self, doIt=True):
        self.m_useError=doIt
        return
                                          
    def  SetRange(self,min, max):
        self.m_userMin=min
        self.m_userMax=max

        return
    
    def  SetNbins(self, Nbins):
        self.m_Nbins=Nbins
        return
    
    def  SetBinning(self, binning):
        self.m_binning = binning
        return
	

    def SetMinMax(self, histo):
	self.m_min = histo.GetBinLowEdge(1)
	self.m_max = histo.GetBinLowEdge(histo.GetNbinsX()+1)
	if self.m_interpolate:
            self.m_min = histo.GetBinCenter(1)
            self.m_max = histo.GetBinCenter(histo.GetNbinsX())
	if self.m_userMin!=-99:
            if self.m_userMin>=self.m_min:
                self.m_min=histo.GetBinCenter(histo.FindBin(self.m_userMin+1e-5))
	    else:
                print "JES_Smoothing::SetMinMax: cannot evaluate histogram in %.1f<pT<%.1f (ignoring this range)\n" %(self.m_userMin,self.m_min)
                
	    if self.m_userMax<=self.m_max:
                self.m_max=self.m_userMax
	    else:
                print "JES_Smoothing::SetMinMax: cannot evaluate histogram in %.1f<pT<%.1f (ignoring this range)\n" %(self.m_max,self.m_userMax)
	    
        return
 
    def GetSmoothedValue(self, histo, my_pT):
	if self.m_min==-99:
            self.SetMinMax(histo)
	if my_pT<=self.m_min:
            my_pT=self.m_min+1e-5
	if my_pT>self.m_max:
            my_pT=self.m_max-1e-5

	sumw=0
        sumwy=0
	pTwidth=self.m_width*my_pT
        
	if my_pT<5:
            pTwidth=self.m_width*5
	for istep in range(self.m_Nsteps):
	    pT = my_pT + pTwidth*(2.0*istep/self.m_Nsteps - 1.0)*self.m_Nsigma
	    if (pT<=self.m_min or pT>=self.m_max):
                continue
	    w = ROOT.TMath.Gaus(my_pT,pT,pTwidth)
            
	    y=0
	    if self.m_interpolate:
                y = histo.Interpolate(pT)
	    else:
                y = histo.GetBinContent(histo.FindBin(pT))
                
	    if self.m_useError:
		bin=histo.FindBin(pT)
                err = histo.GetBinError(bin)
		w *= 1.0/err/err
		
	    sumw+=w
            sumwy+=w*y
	
	return sumwy/sumw


    def SmoothHisto(self,histo):
	# 1. Create output histogram
	min=40
        max=450
        binWidth=1 # bin edges in GeV
        
	if self.m_userMin!=-99 :
	    min=self.m_userMin
            max=self.m_userMax
	
	Nbins=int((max-min)/binWidth)
	if self.m_Nbins!=-99:
            Nbins=self.m_Nbins
	smoothi=0 #give each new histo unique name
	
	if len(self.m_binning)==0:
	    out = ROOT.TH1D("SmoothedSpectrum%d" %(++smoothi),"",Nbins,min,max)
	else:
	    Nbins = self.m_binning.size() - 1;
	    out = ROOT.TH1D("SmoothedSpectrum%d" %(++smoothi),"", Nbins, array.arry('d',self.m_binning))
            

	# 2. Find the range to work with
	self.SetMinMax(histo);
	if self.m_useError:
	    if histo.GetSumw2N()==0:
		print "WARNNG in GetSmoothHisto:"
		print "  input histogram has no bin uncertainties.\n"
		print "  UseBinUncertainty therefore disabled\n"
		self.m_useError=false
	

	    # 3. Do the work
	for bini in range(Nbins):
	    pT=out.GetBinCenter(bini);
	    y=self.GetSmoothedValue(histo,pT)
	    out.SetBinContent(bini,y)
	
        return out

