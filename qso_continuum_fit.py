import pydl
import numpy as np
from pydl.pydlutils import bspline
import matplotlib.pyplot as plt
np.set_printoptions(threshold=np.inf)
def qso_continuum_fit(log_lam, flux, invar):
    '''
    performs bspline iteration fitting on quasar spectrum using python and pyidl
    install pyidl and pyidlutils via pip before using this
    input - log_lam - log of wavelength (numpy array)
            flux - quasar_unormalized_flux (numpy array)
            invar - Inverse of variance (numpy array)
    output - spline_fitted model (numpy array of length same as flux array)     
    '''
    firstpass = 30
    smoothscale = 20
    weightsmooth = 20
    firstupper = 10.0
    firstlower = 2.0
    firstneighbour = int(smoothscale / 2)
    bkptdist = 10000
    npix = len(log_lam)
    mask = invar <= 0
    tempivar = invar
    maxiter = 70
    
    for iiter in range(maxiter):
        weight = np.where(flux > 0, flux**2, 0) * tempivar + 1.0
        weight = pydl.smooth(weight**1.2, weightsmooth)
        maxbkptdist = (np.sum(weight) / 50.0)
        minbkptdist = (np.sum(weight) / 170.0)
        tmpbkptdist = max(min(bkptdist, maxbkptdist), minbkptdist)
        weight = weight/tmpbkptdist
        sumweight = weight
        for i in range(1, len(weight)):
            sumweight[i] = sumweight[i - 1] + weight[i]
        bkpt_indices = pydl.uniq(sumweight.astype(int))
        # Remove consecutive indices with a difference of 1
        cleaned_indices = []
        i = 0
        while i < len(bkpt_indices):
            if i < len(bkpt_indices) - 1 and (bkpt_indices[i + 1] - bkpt_indices[i] == 1):
                cleaned_indices.append(bkpt_indices[i])
                i += 2  # Skip the next index if the difference is 1
            else:
                cleaned_indices.append(bkpt_indices[i])
                i += 1

        bkpt_indices = np.array(cleaned_indices)
        bkpt = log_lam[bkpt_indices]
        firstset, _ = bspline.iterfit(log_lam, flux, invvar=tempivar, bkpt=bkpt, maxiter=70, upper=8., lower=firstlower)
        yfit1,mask_bspline = firstset.value(log_lam)
        diff = (flux - yfit1) * np.sqrt(tempivar)
        diffsmooth = pydl.smooth(diff, smoothscale) * np.sqrt(smoothscale)
        diffmask = diffsmooth < (-1.0 * firstlower)
        smoothmask = pydl.smooth(diffmask * smoothscale, smoothscale) > 0
        firstmask = smoothmask & ((diffsmooth * diff) > (firstlower ** 2))
        if np.sum(firstmask) == 0:
            break
        
        mask = mask | (firstmask > 0)
        recoverdiff = (flux - yfit1) * np.sqrt(invar)
        recover = np.where(recoverdiff > -0.5 * firstlower)[0]
        
        if recover.size > 0:
            mask[recover] = 0
        
        tempivar = invar * (mask == 0)
        #print(len(tempivar))
    model1 = yfit1
    return model1
