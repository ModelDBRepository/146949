# neuroplot.py -- Python script with routines for making pylab plots of neural
#   data.

# Script with pylab routines for plotting and processing neural data.
# 
# Last update: 10/31/11 (georgec)

# Import pylab so the functions and classes don't have to.
import pylab

#
# Data file I/O functions
#

# Load numpy arrays from a file. 
def ldnarrs(fstem):
   fname = 'savedpydata/%s.txt' % fstem
   return pylab.loadtxt(fname)

# Save numpy arrays to a file.
def svnarrs(fstem, X):
   fname = 'savedpydata/%s.txt' % fstem
   pylab.savetxt(fname, X)

#
# Data plotting functions
#

# Plot numpy arrays from a file.
# Note that tvec needs to be the first saved vector.  The other vectors
# will plot as separate traces. 
def plotnarrs(fstem, numtraces=1):
   vecsstr = 'tvec,vec'
   for ii in range(numtraces-1):
      vecsstr += ',vec%d' % (ii+2)
   exec('%s = ldnarrs("%s")' % (vecsstr, fstem))
   
   # Plot first plot.
   pylab.plot(tvec, vec)
   
   # Plot the remaining plots.
   for ii in range(numtraces-1):
      exec('pylab.plot(tvec,vec%d)' % (ii+2))

def plot_spike_times(tvec):
   xmarg = 10   # x display margin
   ymarg = 0.1  # y display margin
   if (tvec == []):
      return
   avec = pylab.ones(len(tvec))
   pylab.stem(tvec, avec)
#   pylab.xlim(-xmarg, sim_duration + xmarg)
   pylab.ylim(-ymarg, 1 + ymarg)
   pylab.title('Spike Times')
   pylab.xlabel('Time')

def plot_spike_raster(tvec, vec):
   xmarg = 10  # x display margin
   ymarg = 10  # y display margin
   marks = pylab.scatter(tvec, vec, marker='+')
   pylab.setp(marks, color='b')
#   pylab.xlim(-xmarg + simview_min_t, simview_max_t + xmarg)
   pylab.ylim(-ymarg, max(vec) + ymarg)
   pylab.title('Spike Events')
   pylab.xlabel('Time')
   pylab.ylabel('Unit ID #')

def plot_power_spect(x_fvec,x_pvec,logx=True,logy=False,\
   putupbands=True,putuplabels=True):
   freqs = pylab.copy(x_fvec)
   if ((not logx) and (not logy)):
      pylab.plot(freqs,x_pvec)
      if putuplabels:
         pylab.ylabel('Power')
   elif (logx and (not logy)):
      if (freqs[0] == 0.0):
         freqs[0] = 0.1
      pylab.semilogx(freqs,x_pvec)
      if putuplabels:
         pylab.ylabel('Power')
   elif ((not logx) and logy):
      pylab.semilogy(freqs,x_pvec)
      if putuplabels:
         pylab.ylabel('Log Power')
   else:
      if (freqs[0] == 0.0):
         freqs[0] = 0.1
      pylab.loglog(freqs,x_pvec)
      if putuplabels:
         pylab.ylabel('Log Power')
   if putuplabels:
      pylab.title('Spectrum')
   pylab.xlim(0,freqs[-1])
   if putuplabels:
      pylab.xlabel('Frequency (Hz)')
   if (putupbands):
      put_up_band_bounds(logx,logy)

def put_up_band_bounds(logx=True,logy=False):
   line = pylab.axvline(1, color='k')
   ax = line.get_axes()
   ymax = ax.get_ybound()[-1]
   pylab.axvline(4, color='k')
   pylab.axvline(7.5, color='k')
   pylab.axvline(14, color='k')
   pylab.axvline(30, color='k')
   pylab.axvline(70, color='k')
   if logx:
      if (not logy):
         pylab.text(1.75,ymax * 0.75,'delta')
         pylab.text(4.5,ymax * 0.75,'theta')
         pylab.text(8,ymax * 0.75,'alpha')
         pylab.text(17,ymax * 0.75,'beta')
         pylab.text(33,ymax * 0.75,'gamma')
      else:
         pylab.text(1.75,ymax * 0.10,'delta')
         pylab.text(4.5,ymax * 0.05,'theta')
         pylab.text(8,ymax * 0.10,'alpha')
         pylab.text(17,ymax * 0.05,'beta')
         pylab.text(33,ymax * 0.10,'gamma')
   else:
      if (not logy):
         pylab.text(1,ymax * 0.75,'delta')
         pylab.text(2,ymax * 0.7,'theta')
         pylab.text(8,ymax * 0.75,'alpha')
         pylab.text(20,ymax * 0.7,'beta')
         pylab.text(45,ymax * 0.75,'gamma')
      else:
         pylab.text(1,ymax * 0.10,'delta')
         pylab.text(2,ymax * 0.05,'theta')
         pylab.text(8,ymax * 0.10,'alpha')
         pylab.text(20,ymax * 0.05,'beta')
         pylab.text(45,ymax * 0.10,'gamma')

def plot_specgram(freqs, ts, Pxx, logz=False):
   if logz:
      pylab.imshow(pylab.log10(Pxx), origin='lower', aspect='auto', \
         extent=(0,ts[-1]*1000.0+self.ts[0]*1000.0,0,freqs[-1]))
      pylab.title('Log Spectrogram')
   else:
      pylab.imshow(Pxx, origin='lower', aspect='auto', \
         extent=(0,ts[-1]*1000.0+self.ts[0]*1000.0,0,freqs[-1]))
      pylab.title('Spectrogram')
   pylab.xlabel('Time (ms)')
   pylab.ylabel('Frequency (Hz)')
   pylab.spectral()
   pylab.colorbar()

def plot_band_specgram(bandnames, bts, Ptxx):
   pylab.imshow(Ptxx, origin='lower', aspect='auto', \
      extent=(0,bts[-1]*1000.0+self.bts[0]*1000.0,0,len(bandnames)), \
      interpolation='nearest')
   pylab.title('Band Spectrogram')
   pylab.xlabel('Time (ms)')
   pylab.ylabel('Spectral Band')
   pylab.yticks(pylab.arange(len(bandnames))+0.5, bandnames)
   pylab.spectral()
   pylab.colorbar()

def plot_specgram_spectdist(freqs,ts,Pxx,show='quartiles',logx=True,logy=False,putupbands=True):
   import scipy.stats

   # Get the mean, median, range, and standard deviation scores.
   mean_vec = Pxx.mean(axis=1)      # average over time
   std_vec = Pxx.std(axis=1)        # stdev over time
   upperstd_vec = mean_vec + std_vec
   lowerstd_vec = mean_vec - std_vec
   upperrng_vec = Pxx.max(axis=1)
   lowerrng_vec = Pxx.min(axis=1)

   # Calculate the 25th, 50th (median) and 75th percentile scores.
   p25_vec = pylab.zeros(len(mean_vec))
   median_vec = pylab.zeros(len(mean_vec))
   p75_vec = pylab.zeros(len(mean_vec))
   for ii in range(len(freqs)):
      p25_vec[ii] = scipy.stats.scoreatpercentile(Pxx[ii,:], 25)
      median_vec[ii] = scipy.stats.scoreatpercentile(Pxx[ii,:], 50)
      p75_vec[ii] = scipy.stats.scoreatpercentile(Pxx[ii,:], 75)

   if ((not logx) and (not logy)):
      if (show in ['meanstd','meanmed']):
         pylab.plot(freqs,mean_vec)
      if (show in ['meanmed','quartiles']):
         pylab.plot(freqs,median_vec)
      if (show in ['meanstd']):
         pylab.plot(freqs,upperstd_vec,'g')
         pylab.plot(freqs,lowerstd_vec,'g')
      if (show in ['quartiles']):
         pylab.plot(freqs,p75_vec,'g')
         pylab.plot(freqs,p25_vec,'g') 
      if (show in ['meanstd','meanmed','quartiles']):
         pylab.plot(freqs,upperrng_vec,'r')
         pylab.plot(freqs,lowerrng_vec,'r')
      pylab.ylabel('Power')
   elif (logx and (not logy)):
      if (freqs[0] == 0.0):
         freqs[0] = 0.1
      if (show in ['meanstd','meanmed']):
         pylab.semilogx(freqs,mean_vec)
      if (show in ['meanmed','quartiles']):
         pylab.semilogx(freqs,median_vec)
      if (show in ['meanstd']):
         pylab.semilogx(freqs,upperstd_vec,'g')
         pylab.semilogx(freqs,lowerstd_vec,'g')
      if (show in ['quartiles']):
         pylab.semilogx(freqs,p75_vec,'g')
         pylab.semilogx(freqs,p25_vec,'g') 
      if (show in ['meanstd','meanmed','quartiles']):
         pylab.semilogx(freqs,upperrng_vec,'r')
         pylab.semilogx(freqs,lowerrng_vec,'r')
      pylab.ylabel('Power')
   elif ((not logx) and logy):
      if (show in ['meanstd','meanmed']):
         pylab.semilogy(freqs,mean_vec)
      if (show in ['meanmed','quartiles']):
         pylab.semilogy(freqs,median_vec)
      if (show in ['meanstd']):
         pylab.semilogy(freqs,upperstd_vec,'g')
         pylab.semilogy(freqs,lowerstd_vec,'g')
      if (show in ['quartiles']):
         pylab.semilogy(freqs,p75_vec,'g')
         pylab.semilogy(freqs,p25_vec,'g') 
      if (show in ['meanstd','meanmed','quartiles']):
         pylab.semilogy(freqs,upperrng_vec,'r')
         pylab.semilogy(freqs,lowerrng_vec,'r')
      pylab.ylabel('Log Power')
   else:
      if (freqs[0] == 0.0):
         freqs[0] = 0.1
      if (show in ['meanstd','meanmed']):
         pylab.loglog(freqs,mean_vec)
      if (show in ['meanmed','quartiles']):
         pylab.loglog(freqs,median_vec)
      if (show in ['meanstd']):
         pylab.loglog(freqs,upperstd_vec,'g')
         pylab.loglog(freqs,lowerstd_vec,'g')
      if (show in ['quartiles']):
         pylab.loglog(freqs,p75_vec,'g')
         pylab.loglog(freqs,p25_vec,'g') 
      if (show in ['meanstd','meanmed','quartiles']):
         pylab.loglog(freqs,upperrng_vec,'r')
         pylab.loglog(freqs,lowerrng_vec,'r')
      pylab.ylabel('Log Power')
   pylab.title('Spectrogram Spectrum Distribution')
   pylab.xlim(0,freqs[-1])
   pylab.xlabel('Frequency (Hz)')
   if (putupbands):
      put_up_band_bounds(logx,logy)

#
# Data processing / analysis functions
#

def get_vec_subset(tvec, vec, mintvec=-1, maxtvec=-1, minvec=-1, maxvec=-1):
   if (mintvec == -1):
      mintvec = tvec.min()
   if (maxtvec == -1):
      maxtvec = tvec.max() + 0.0005  # allow inclusion of the max
   if (minvec == -1):
      minvec = vec.min()
   if (maxvec == -1):
      maxvec = vec.max() + 0.0005    # allow inclusion of the max
   filt_tvec = []
   filt_vec = []
   for ii in range(len(tvec)):
      # If the tvec and vec values are in the window, then include them.  Note
      # that the min boundaries include the minimum value, but the max bounds
      # do not include the max value, only values up to the max.
      if ((tvec[ii] >= mintvec) and (tvec[ii] < maxtvec) and \
         (vec[ii] >= minvec) and (vec[ii] < maxvec)):
         filt_tvec.append(tvec[ii])
         filt_vec.append(vec[ii])
   filt_tvec = pylab.array(filt_tvec)
   filt_vec = pylab.array(filt_vec)
   return filt_tvec, filt_vec

# Downsample a time sequence.
#   x_tvec - time stamps for sequence (units assumed to be in ms)
#   x_vec - amplitudes for sequence
#   newfs - new sample frequency (in Hz, must evenly divide the old fs)  
def downsample (x_tvec, x_vec, newfs):
   fsratio = (1000.0 / (x_tvec[1] - x_tvec[0])) / float(newfs)
   
   # If the fs ratio is invalid, give an error.
   if ((float(int(fsratio)) != fsratio) or (fsratio < 1)):
      print 'Error: The ratio of oldfs and newfs must be a positive integer'
      return None
      
   # Integerize the ratio.
   fsratio = int(fsratio)
   
   # Do the downsampling by averaging x_vec blockwise by blocks of size fsratio.
   npts = len(x_vec) / fsratio
   tvec = pylab.zeros(npts)
   vec = pylab.zeros(npts)
   for ii in range(npts):
      tvec[ii] = x_tvec[ii * fsratio]
      vec[ii] = x_vec[ii * fsratio:((ii+1) * fsratio)].mean()
      
   return tvec, vec

# Use convolution with a Hanning window to smooth a vector.
def smooth(x, win_size=11):
   y = pylab.zeros(x.size)
   smwin = pylab.hanning(win_size)
   smwin /= sum(smwin)  # scale the window so the total area is 1.0
   filtout = pylab.convolve(x,smwin,mode='valid')
   y[0:filtout.size] = filtout
   return y

def get_ave_fire_freqs(spks_tvec,spks_vec,num_units,fire_dur=0.0):
   if (fire_dur == 0.0):
      fire_dur = spks_tvec.max() - spks_tvec.min()
   fire_counts = pylab.zeros(num_units)
   for ii in range(len(spks_tvec)):
      fire_counts[spks_vec[ii]] += 1
   return fire_counts * 1000.0 / fire_dur

def get_unit_fire_bins(spks_tvec,spks_vec,num_units,fire_dur=0.0,bin_dur=100.0):
   if (fire_dur == 0.0):
      fire_dur = spks_tvec.max() - spks_tvec.min()
   num_bins = int(fire_dur / bin_dur)
   fire_bin_counts = pylab.zeros((num_units,num_bins))
   for ii in range(num_units):
      tvec, vec = get_vec_subset(spks_tvec, spks_vec, minvec=ii, maxvec=ii+1)
      n, bins = pylab.histogram(tvec,bins=num_bins,range=(0.0,fire_dur))
      fire_bin_counts[ii,:] = n[:]
   return bins,fire_bin_counts

# Given full spectrogram information, return a matrix where columns are all of 
# the time slices and rows are frequency band and the values are the sum, within 
# the time slice, of all of the power in that band.
def get_band_waterfall(freqs,ts,Pxx,normlz=False):
   bandnames = ['<Delta','Delta','Theta','Alpha','Beta','Gamma','>Gamma']
   bts = pylab.copy(ts)
   Ptxx = pylab.zeros((len(bandnames),len(ts)))
   for tt in range(len(ts)):
      for ff in range(len(freqs)):
         if (freqs[ff] < 1.0):
            Ptxx[0,tt] += Pxx[ff,tt]  # sub-delta
         elif ((freqs[ff] >= 1.0) and (freqs[ff] < 4.0)):
            Ptxx[1,tt] += Pxx[ff,tt]  # delta
         elif ((freqs[ff] >= 4.0) and (freqs[ff] < 7.5)):
            Ptxx[2,tt] += Pxx[ff,tt]  # theta
         elif ((freqs[ff] >= 7.5) and (freqs[ff] < 14.0)):
            Ptxx[3,tt] += Pxx[ff,tt]  # alpha
         elif ((freqs[ff] >= 14.0) and (freqs[ff] < 30.0)):
            Ptxx[4,tt] += Pxx[ff,tt]  # beta
         elif ((freqs[ff] >= 30.0) and (freqs[ff] < 70.0)):
            Ptxx[5,tt] += Pxx[ff,tt]  # gamma
         else:
            Ptxx[6,tt] += Pxx[ff,tt]  # super-gamma
      if (normlz):
         Ptxx[:,tt] /= sum(Ptxx[:,tt])
   return bandnames,bts,Ptxx

#
# Data plotting / browsing classes
#

class SpecgramBrowser:
   """
   Spectrogram browser for looking at spectra time-slices in a waterfall plot.

   Methods:
     constructor(xxx):

   Attributes:
     count (static): number of instances

   Usage:
   >>> sb = SpecgramBrowser(freqs,ts,Pxx)
   """

   count = 0

   def __init__(self, freqs, ts, Pxx, Pxx2=None, logz=False, logx=True, logy=False, \
      putupbands=True, nointerp=False):
      # Set up browser parameters.
      self.freqs = freqs
      self.ts = ts
      self.Pxx = Pxx
      self.Pxx2 = Pxx2
      self.logz = logz
      self.logx = logx
      self.logy = logy
      self.ymaxlock = False
      self.putupbands = putupbands
      self.nointerp = nointerp
      self.specgram_fig = pylab.figure()  # handle to specgram_browser spectrogram figure
      if (self.Pxx2 == None):
         self.specgram_fig2 = 0   # handle to specgram_browser spectrogram 2 figure
      else:
         self.specgram_fig2 = pylab.figure()  # handle to specgram_browser spectrogram figure
      self.spect_fig = pylab.figure()     # handle to specgram_browser spectrum figure
      self.tindex = 0           # index to the time slice to be shown in the spectrum figure/s
      self.specgram_mark = 0    # handle to spectrogram time slice marker
      self.specgram_mark2 = 0   # handle to spectrogram 2 time slice marker

      # Increment the instance count.
      SpecgramBrowser.count += 1 

      # Draw the initial spectrogram/s.
      self._draw_specgram()

      # Draw the initial spectrum/a.
      self._draw_spect()

      # Attach an onclick event to grab mouse clicks for spectrogram window.
      cid = self.specgram_fig.canvas.mpl_connect('button_press_event', self._onclick)

      # Attach an onkey event to grab key presses for spectrogram window.
      cid2 = self.specgram_fig.canvas.mpl_connect('key_press_event', self._onkey)

      if (self.Pxx2 != None):
         # Attach an onclick event to grab mouse clicks for spectrogram window.
         cid3 = self.specgram_fig2.canvas.mpl_connect('button_press_event', self._onclick)

         # Attach an onkey event to grab key presses for spectrogram window.
         cid4 = self.specgram_fig2.canvas.mpl_connect('key_press_event', self._onkey)

      # Attach an onkey event to grab key presses for spectrum window.
      cid5 = self.spect_fig.canvas.mpl_connect('key_press_event', self._onkey2)

   def _draw_specgram(self):
      def figdrawer(self, Pxx):
         pylab.clf()
         if self.logz:
            if self.nointerp:
               pylab.imshow(pylab.log10(Pxx), origin='lower', aspect='auto', \
                  extent=(0,self.ts[-1]*1000.0+self.ts[0]*1000.0,0,self.freqs[-1]), \
                  interpolation='nearest')
            else:
               pylab.imshow(pylab.log10(Pxx), origin='lower', aspect='auto', \
                  extent=(0,self.ts[-1]*1000.0+self.ts[0]*1000.0,0,self.freqs[-1]))
            pylab.title('Log Spectrogram')
         else:
            if self.nointerp:
               pylab.imshow(Pxx, origin='lower', aspect='auto', \
                  extent=(0,self.ts[-1]*1000.0+self.ts[0]*1000.0,0,self.freqs[-1]), \
                  interpolation='nearest')
            else:
               pylab.imshow(Pxx, origin='lower', aspect='auto', \
                  extent=(0,self.ts[-1]*1000.0+self.ts[0]*1000.0,0,self.freqs[-1]))
            pylab.title('Spectrogram')
         pylab.xlabel('Time (ms)')
         pylab.ylabel('Frequency (Hz)')
         pylab.spectral()
         pylab.colorbar()

      pylab.figure(self.specgram_fig.number)
      figdrawer(self, self.Pxx)
      self.specgram_mark = pylab.axvline(self.ts[self.tindex] * 1000.0,color='r')
      if (self.Pxx2 != None):
         pylab.figure(self.specgram_fig2.number)
         figdrawer(self, self.Pxx2)
         self.specgram_mark2 = pylab.axvline(self.ts[self.tindex] * 1000.0,color='r')

   def _move_specgram_mark(self):
      # Reset the line marker according to the new time slice index.  
      self.specgram_mark.set_xdata(pylab.array([self.ts[self.tindex] * 1000.0,self.ts[self.tindex] * 1000.0]))

      # Redraw the spectrogram canvas.
      self.specgram_fig.canvas.draw()

      # If we have another spectrogram...
      if (self.Pxx2 != None):
         # Reset the line marker according to the new time slice index.  
         self.specgram_mark2.set_xdata(pylab.array([self.ts[self.tindex] * \
            1000.0,self.ts[self.tindex] * 1000.0]))

         # Redraw the spectrogram canvas.
         self.specgram_fig2.canvas.draw()

   def _draw_spect(self):
      pylab.figure(self.spect_fig.number)
      pylab.clf()
      freqs2 = pylab.copy(self.freqs)
      x_pvec = pylab.copy(self.Pxx[:,self.tindex])
      if (self.Pxx2 != None):
         x_pvec2 = pylab.copy(self.Pxx2[:,self.tindex])
      if ((not self.logx) and (not self.logy)):
         pylab.plot(freqs2,x_pvec)
         if (self.Pxx2 != None):
            pylab.plot(freqs2,x_pvec2,'r')
         pylab.ylabel('Power')
      elif (self.logx and (not self.logy)):
         if (freqs2[0] == 0.0):
            freqs2[0] = 0.1
         pylab.semilogx(freqs2,x_pvec)
         if (self.Pxx2 != None):
            pylab.semilogx(freqs2,x_pvec2,'r')
         pylab.ylabel('Power')
      elif ((not self.logx) and self.logy):
         pylab.semilogy(freqs2,x_pvec)
         if (self.Pxx2 != None):
            pylab.semilogy(freqs2,x_pvec2,'r')
         pylab.ylabel('Log Power')
      else:
         if (freqs2[0] == 0.0):
            freqs2[0] = 0.1
         pylab.loglog(freqs2,x_pvec)
         if (self.Pxx2 != None):
            pylab.loglog(freqs2,x_pvec2,'r')
         pylab.ylabel('Log Power')
      if (self.Pxx2 == None):
         pylab.title('Spectrum at %.2f s' % self.ts[self.tindex])
      else:
         pylab.title('Spectra at %.2f s' % self.ts[self.tindex])
      pylab.xlim(0,freqs2[-1])
      if (self.ymaxlock):
         maxlock = self.Pxx.max()
         if (self.Pxx2 != None):
            maxlock = max(maxlock, self.Pxx2.max())
         maxlock *= 1.0
         pylab.ylim(0,maxlock)
      pylab.xlabel('Frequency (Hz)')
      if (self.putupbands):
         put_up_band_bounds(self.logx,self.logy)

   def _onclick(self,event):
      print 'Selected...'
#      print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f' % \
#         (event.button, event.x, event.y, event.xdata, event.ydata)

      # Create freq bin upper boundaries and array of indices where value 
      # falls in bin.
      fbinrad = self.freqs[1] / 2.0
      fbinubs = self.freqs + fbinrad
      fis = pylab.find(fbinubs >= event.ydata)

      # Create time bin upper boundaries and array of indices where value
      # falls in bin.
      tbinrad = (self.ts[1] - self.ts[0]) / 2.0
      tbinubs = (self.ts + tbinrad) * 1000.0
      tis = pylab.find(tbinubs >= event.xdata)

      # Show the bin and value.
      print 'freq bin=%f, time bin=%f, value=%f' % (self.freqs[fis[0]], self.ts[tis[0]], self.Pxx[fis[0]][tis[0]])

      # Switch the spectrum view to the corresponding time slice.
      self.tindex = tis[0]
      self._draw_spect()

      # Redraw the spectrogram (moving the marker).
      self._draw_specgram()

   def _onkey(self,event):
      # Note: for reasons I don't understand, Tkinter callbacks trap the 
      # following keystrokes, so I cannot use them:
      #   's'
      #   'l'   
#      print 'you pressed', event.key, event.xdata, event.ydata
      if (event.key == 'left'):
         if (self.tindex > 0):
            self.tindex -= 1
         self._draw_spect()
         self._draw_specgram()
      elif (event.key == 'right'):
         if (self.tindex < len(self.ts)-1):
            self.tindex += 1
         self._draw_spect()
         self._move_specgram_mark()
      elif (event.key == 'i'):
         self.nointerp = not self.nointerp
         self._draw_specgram()

   def _onkey2(self,event):
      # Note: for reasons I don't understand, Tkinter callbacks trap the 
      # following keystrokes, so I cannot use them:
      #   's'
      #   'l'  
#      print 'you pressed', event.key, event.xdata, event.ydata
      if (event.key == 'left'):
         if (self.tindex > 0):
            self.tindex -= 1
         self._draw_spect()
         self._move_specgram_mark()
      elif (event.key == 'right'):
         if (self.tindex < len(self.ts)-1):
            self.tindex += 1
         self._draw_spect()
         self._move_specgram_mark()
      elif (event.key == 'x'):
         self.logx = not self.logx
         self._draw_spect()
      elif (event.key == 'y'):
         self.logy = not self.logy
         self._draw_spect()
      elif (event.key == 'b'):
         self.putupbands = not self.putupbands
         self._draw_spect()
      elif (event.key == 'm'):
         self.ymaxlock = not self.ymaxlock
         self._draw_spect()
