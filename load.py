# load.py -- Script for loading and analyzing data for the Gordon-cell
#   M1 model.
#
# Updated: 11/1/12 (georgec)

#
# Global parameters
#

# Use the NEURON GUI?
use_NEURON_GUI = True

# Set up directory for paper figs.
paperfigsdir = '.'

# Index to grvec panel object for latest loaded sim.
curr_grvec_pind = 1

#
# Misc functions
#

def loadsimbyfstem (filestem,grvec_load=False):
   global num_cells
   global cell_types
   global cell_zlocs
   global curr_grvec_pind

   # Get LFPs NQS (nqLFP).
   h('if (nqLFP != nil) nqsdel(nqLFP)')
   h('nqLFP = new NQS()')
   fname = '%s-lfp.nqs' % filestem
   h.nqLFP.rd(fname)

   # Get spikes NQS (snq).
   h('if (snq != nil) nqsdel(snq)')
   h('snq = new NQS()')
   fname = '%s-spk.nqs' % filestem
   h.snq.rd(fname)
   spks_tvec = nqscol2narr('snq','t')
   spks_idvec = nqscol2narr('snq','id')

   # Get cells info NQS (cellsnq).
   h('if (cellsnq != nil) nqsdel(cellsnq)')
   h('cellsnq = new NQS()')
   fname = '%s-loc.nqs' % filestem
   h.cellsnq.rd(fname)

   # Get connections info NQS (connsnq).
   h('if (connsnq != nil) nqsdel(connsnq)')
   h('connsnq = new NQS()')
   fname = '%s-con.nqs' % filestem
   h.connsnq.rd(fname)

   # Set up second connsnq variable for info extraction.
   h('if (connsnq2 != nil) nqsdel(connsnq2)')
   h('connsnq2 = new NQS()')

   # Load information from the cellsnq NQS table.
   num_cells = int(h.cellsnq.size())
   cell_types = nqscol2narr('cellsnq','ty')
   cell_zlocs = nqscol2narr('cellsnq','zloc')

   # If the NEURON GUI is up, load the grvec printlist from the saved data.
   if use_NEURON_GUI and grvec_load:
      ldgrvec('%s-prlist.gvc' % outfilestem)
      curr_grvec_pind += 1

def ctypandind2cind (ctyp,ctind):
   cell_inds = find(cell_types == get_ctyp_num(ctyp))
   return cell_inds[ctind]

def cind2ctypandind (cell_num):
   # Get the cell type number.
   ctypnum = cell_types[cell_num]

   # Get the first cell index for that type.
   ctypfind = find(cell_types == ctypnum)[0]

   # Show the cell type string and the relative cell index.
   return (get_ctyp_str(ctypnum), cell_num - ctypfind)

#
# Data analysis functions
#

def get_ctyp_fire_rate (ctyp,min_t=0.0,max_t=-1.0,allcellrates=False):
   if (max_t == -1.0):
      max_t = sim_duration
   tvec = nqscol2narr('snq','t')
   idvec = nqscol2narr('snq','id')
   tvec,idvec = get_vec_subset(tvec,idvec,mintvec=min_t,maxtvec=max_t)
   freqs = get_ave_fire_freqs(tvec,idvec,num_cells,max_t - min_t)
   cell_inds = find(cell_types == get_ctyp_num(ctyp))
   if (allcellrates):
      return freqs[cell_inds]
   else:
      return freqs[cell_inds].mean()

# get_curr_fire_rate() -- get an estimate of a firing rate for a target cell at a 
#    particular time.  By default, the counting window is 100 ms and starts at time t.
#    Requires that snq table is loaded.
def get_curr_fire_rate (targtype='',targind=0,t=0.0,spkwindwid=100.0,spkwindoffset=0.0):
   # Remember the old snq table verbosity.
   oldverbose = h.snq.verbose

   # Turn off the table verbosity.
   h.snq.verbose = 0

   # Default the cell index to targind.
   cell_ind = targind

   # If we have a cell type, get the cell IDs and get the indexed value.
   if (targtype != ''):
      targ_cells = find(cell_types == get_ctyp_num(targtype)) 
      cell_ind = targ_cells[targind]     
   
   # Get number all of the desired cell spikes in the desired time window.
   num_spks = h.snq.select('id',cell_ind,'t','[]',t+spkwindoffset, \
      t+spkwindoffset+spkwindwid)

   # Get the firing rate from the spike count and window width.
   fire_rate = num_spks * (1000.0 / spkwindwid)

   # Reset the snq table.
   h.snq.tog()
   h.snq.verbose = oldverbose

   return fire_rate

def get_avg_ctyp_rates (min_t=0.0,max_t=-1.0):
   ctyps = ['E2','I2','I2L','E5R','E5B','I5','I5L','E6','I6','I6L']
   rates = zeros(len(ctyps))
   for rind in range(len(ctyps)):
      rates[rind] = get_ctyp_fire_rate(ctyps[rind],min_t,max_t)
   return ctyps,rates

def get_ctyp_spk_count (ctyp,min_t=0.0,max_t=-1.0):
   if (max_t == -1.0):
      max_t = sim_duration

   # Get all of the spikes.
   tvec = nqscol2narr('snq','t')
   idvec = nqscol2narr('snq','id')

   # Get all of the cell indices for the desire cell type.
   cell_inds = find(cell_types == get_ctyp_num(ctyp))

   # Get only the spikes of the right cell type and in the time window.
   tvec,idvec = get_vec_subset(tvec,idvec,mintvec=min_t,maxtvec=max_t, \
      minvec=cell_inds[0],maxvec=cell_inds[-1])

   # Return the numbrer of spikes.
   return len(tvec)

def get_bin_integral (F,xs,xstart,xend): 
   # Get the relevant bin indices using xs, xstart, and xend.
   bin_inds = find((xs >= xstart) & (xs <= xend))

   # Start with the sum being the sum all of the bins but the first and last.
   sumval = F[bin_inds[1:-1]].sum()

   # Add the portion of the left bin we have.
   sumval += (F[bin_inds[0]] * (xs[bin_inds[0]] - xstart) / (xs[1] - xs[0]))

   # Add the portion of the right bin we have.
   sumval += (F[bin_inds[-1]] * (xend - xs[bin_inds[-1]]) / (xs[1] - xs[0]))

   return sumval

# get_conv_conns() -- load the (static) convergence connections into a unit into 
#    the NQS table connsnq2.  Requires the connection table connsnq to be loaded.
def get_conv_conns (srctype='E2',targtype='E5R',targind=0):
   # Remember the old connsnq table verbosity, and turn it off.
   oldverbose = h.connsnq.verbose
   h.connsnq.verbose = 0

   # Remember the old connsnq2 table verbosity, and turn it off.
   oldverbose2 = h.connsnq2.verbose
   h.connsnq2.verbose = 0

   # Get source cell IDs and save in v1.
   src_cells = find(cell_types == get_ctyp_num(srctype))
   narr2hv('v1',src_cells)

   # Get target cell IDs.
   targ_cells = find(cell_types == get_ctyp_num(targtype))

   # Select the entries with the desired weights and return the number of connections found.
   nc = h.connsnq.select("id1","EQW",h.v1,"id2",targ_cells[targind])
   
   # If we found connections...
   if (nc > 0):
      # Copy those weights over to connsnq2.
      h.connsnq2.cp(h.connsnq.out)

      # Sort the list.
      h.connsnq2.sort("id1")

   # Reset the connsnq table.
   h.connsnq.tog("db")
   h.connsnq.verbose = oldverbose

   # Reset the connsnq2 table.
   h.connsnq2.verbose = oldverbose2

   # Return the number of connections found.
   return nc

# get_div_conns() -- load the (static) divergence connections out of a unit into 
#    the NQS table connsnq2.  Requires the connection table connsnq to be loaded.
def get_div_conns (srctype='E2',targtype='E5R',srcind=0):
   # Remember the old connsnq table verbosity, and turn it off.
   oldverbose = h.connsnq.verbose
   h.connsnq.verbose = 0

   # Remember the old connsnq2 table verbosity, and turn it off.
   oldverbose2 = h.connsnq2.verbose
   h.connsnq2.verbose = 0

   # Get source cell IDs.
   src_cells = find(cell_types == get_ctyp_num(srctype))

   # Get target cell IDs and save in v1.
   targ_cells = find(cell_types == get_ctyp_num(targtype))
   narr2hv('v1',targ_cells)

   # Select the entries with the desired weights and return the number of connections found.
   nc = h.connsnq.select("id1",src_cells[srcind],"id2","EQW",h.v1)

   # If we found connections...
   if (nc > 0):   
      # Copy those weights over to connsnq2.
      h.connsnq2.cp(h.connsnq.out)

      # Sort the list.
      h.connsnq2.sort("id2")

   # Reset the connsnq table.
   h.connsnq.tog("db")
   h.connsnq.verbose = oldverbose

   # Reset the connsnq2 table.
   h.connsnq2.verbose = oldverbose2

   # Return the number of connections found.
   return nc

# get_conv_wts() -- for the convergent weights into a unit (and potentially at a 
#    particular time) return an array of the source cell IDs and the weights (static 
#    or plastic).
def get_conv_wts (srctype='E2',targtype='E5R',targind=0):
   # Get source cell IDs.
   src_cells = find(cell_types == get_ctyp_num(srctype))

   # Get target cell IDs.
   targ_cells = find(cell_types == get_ctyp_num(targtype))

   # Get the appropriate connections into connsnq2 table.
   nc = get_conv_conns(srctype,targtype,targind)

   # If we actually have some connections...
   if (nc > 0):
      # Get IDs and subtract first index for type.
      nqids = nqscol2narr('connsnq2','id1')
      nqids -= src_cells[0]
      nqids = array(nqids,int)

      # Get the primary weights from the table.
      nqwts = nqscol2narr('connsnq2','wt1')
   else:
      nqids = []
      nqwts = []

   return [nqids, nqwts]

def get_avg_conv_wts (srctype='E2',targtype='E5R'):
   # Get target cell IDs.
   targ_cells = find(cell_types == get_ctyp_num(targtype))

   # Set up an array for the weights.
   avg_conv_wts = zeros(len(targ_cells))

   # For all of the target cells...
   for ii in range(len(targ_cells)):
      # Read the convergent weight IDs and weight values.
      nqids, nqwts = get_conv_wts(srctype, targtype, ii)

      # Save the average weight from this.
      if (len(nqwts) > 0):
         avg_conv_wts[ii] = nqwts.mean()
      else:
         avg_conv_wts[ii] = 0

   # Return the average of averages.
   return avg_conv_wts.mean()

#
# Display functions
#

def pravgrates (min_t=0.0,max_t=-1.0):
   print 'Average Firing Rates by Cell Type:'
   print '   E2: %.2f Hz' % get_ctyp_fire_rate('E2',min_t,max_t)
   print '   I2: %.2f Hz' % get_ctyp_fire_rate('I2',min_t,max_t)
   print '  I2L: %.2f Hz' % get_ctyp_fire_rate('I2L',min_t,max_t)
   print '  E5R: %.2f Hz' % get_ctyp_fire_rate('E5R',min_t,max_t)
   print '  E5B: %.2f Hz' % get_ctyp_fire_rate('E5B',min_t,max_t)
   print '   I5: %.2f Hz' % get_ctyp_fire_rate('I5',min_t,max_t)
   print '  I5L: %.2f Hz' % get_ctyp_fire_rate('I5L',min_t,max_t)
   print '   E6: %.2f Hz' % get_ctyp_fire_rate('E6',min_t,max_t)
   print '   I6: %.2f Hz' % get_ctyp_fire_rate('I6',min_t,max_t)
   print '  I6L: %.2f Hz' % get_ctyp_fire_rate('I6L',min_t,max_t)

#
# Plotting functions
#

def plot_raster (startt=0,endt=1000,xmarg=10,ymarg=10,drawpopbounds=True,bwmode=False):
   def plotcells (ctyp,pcol,drawpopbounds=True,drawpoplabels=True):
      # Select only the cells in the right time window and of the right type.
      spk_count = h.snq.select("t","[]",startt,endt,"type",get_ctyp_num(ctyp))

      # If there are spikes...
      if (spk_count):
         # Get the spikes and show them.
         ts = nqscol2narr('snq','t')
         ts /= 1000.0
         cs = nqscol2narr('snq','id')
         marks = scatter(ts, cs, marker='+')
         setp(marks, color=pcol)

      # Reset the snq table.
      h.snq.tog("DB")

      # If we want to draw population boundaries...
      if (drawpopbounds):
         cell_inds = find(cell_types == get_ctyp_num(ctyp))
         axhline(cell_inds[0] - 0.5,color='k')
         axhline(cell_inds[-1] + 0.5,color='k') 

      # If we want to draw population labels...
      if (drawpoplabels):
         cell_inds = find(cell_types == get_ctyp_num(ctyp))
         # If E5B, convert cell type to 'SPI'.
         if (ctyp == 'E5B'):
            ctyp = 'SPI'
         elif (ctyp == 'E5R'):
            ctyp = 'STR'
         text((xmarg + endt) / 1000.0 + ((endt - startt) / 1000.0 / 80.0), (cell_inds[0] + cell_inds[-1]) / 2.0 - 12, ctyp)

   # If not black-and-white mode...
   if (not bwmode):
      # Plot the cells with color-coding.
#      plotcells('E2','b',drawpopbounds)
#      plotcells('I2','r',drawpopbounds)  
#      plotcells('I2L','orange',drawpopbounds)
##      plotcells('E4','b',drawpopbounds) 
##      plotcells('I4','r',drawpopbounds)  
##      plotcells('I4L','orange',drawpopbounds)   
#      plotcells('E5R','b',drawpopbounds) 
#      plotcells('E5B','g',drawpopbounds)    
#      plotcells('I5','r',drawpopbounds)  
#      plotcells('I5L','orange',drawpopbounds)
#      plotcells('E6','b',drawpopbounds)  
#      plotcells('I6','r',drawpopbounds) 
#      plotcells('I6L','orange',drawpopbounds)

      plotcells('E2','b',drawpopbounds)
      plotcells('I2','orange',drawpopbounds)  
      plotcells('I2L','m',drawpopbounds)
#      plotcells('E4','b',drawpopbounds) 
#      plotcells('I4','orange',drawpopbounds)  
#      plotcells('I4L','m',drawpopbounds)   
      plotcells('E5R','g',drawpopbounds) 
      plotcells('E5B','r',drawpopbounds)    
      plotcells('I5','orange',drawpopbounds)  
      plotcells('I5L','m',drawpopbounds)
      plotcells('E6','c',drawpopbounds)  
      plotcells('I6','orange',drawpopbounds) 
      plotcells('I6L','m',drawpopbounds)
   else:
      # Plot the cells with black.
      plotcells('E2','k',drawpopbounds)
      plotcells('I2','k',drawpopbounds)  
      plotcells('I2L','k',drawpopbounds)
#      plotcells('E4','k',drawpopbounds) 
#      plotcells('I4','k',drawpopbounds)  
#      plotcells('I4L','k',drawpopbounds)   
      plotcells('E5R','k',drawpopbounds) 
      plotcells('E5B','k',drawpopbounds)    
      plotcells('I5','k',drawpopbounds)  
      plotcells('I5L','k',drawpopbounds)
      plotcells('E6','k',drawpopbounds)  
      plotcells('I6','k',drawpopbounds) 
      plotcells('I6L','k',drawpopbounds)

   # Set the x limits according to the min and max time.
   xlim((startt - xmarg) / 1000.0, (endt + xmarg) / 1000.0)

   # Set the y limits according to which cells exist.
   ylim(-ymarg, num_cells + ymarg - 1)

   # Label the axes.
   xlabel('Time (s)')
   ylabel('Cell #')

def plot_cell_densities ():
   plot(lyrbins,M1dens)
   title('Cell Density vs. Depth')
   xlabel('Distance from pia (microns)')
   ylabel('Cell Density (E and I)')
   axvline(143.0,color='k') # L1 / L2/3 boundary
   axvline(451.8,color='k') # L2/3 / L5A boundary
   axvline(663.6,color='k') # L5A / L5B boundary
   axvline(1059.0,color='k') # L5B / L6 boundary
   text(70,0.11,'L1')
   text(300,0.11,'L2/3')
   text(525,0.11,'L5A')
   text(850,0.11,'L5B')
   text(1250,0.11,'L6')

def plot_vtrace (tvec,vec,ctyp='E2'):
   plot(tvec,vec)
   ylim(-80,60)
   if (ctyp in ('E2','E5R','E6')):
      axhline(-25,color='r')
      axhline(-40,color='k')
      axhline(-65,color='k',linestyle='--')
   elif (ctyp == 'E5B'):
      axhline(-25,color='r')
      axhline(-37.5,color='k')
      axhline(-65,color='k',linestyle='--')
   elif (ctyp in ('I2','I5','I6')):
      axhline(-40,color='k')
      axhline(-63,color='k',linestyle='--')
   elif (ctyp in ('I2L','I5L','I6L')):
      axhline(-47,color='k')
      axhline(-65,color='k',linestyle='--')
   title('%s Voltage Trace' % ctyp)
   xlabel('Time (ms)')
   ylabel('Memb. Pot. (mV)')

def plot_inter_ctyp_matrix (ctyp_names, vals):
   # Make the reverse list of names.
   rnames=[]
   rnames.extend(ctyp_names)
   rnames.reverse() 

   npops = size(ctyp_names)
   figh = figure(figsize=(9,8)) # Open figure and store its handle
   axh = subplot(111)
   figh.subplots_adjust(left=0.1) # Less space on right
   figh.subplots_adjust(right=1.1) # Less space on right
   ex = array(range(npops+1))
   why = npops - ex
   tmp = axh.pcolor(ex,why,vals.transpose(),linestyle='solid',linewidth=0.5,color='k')
   xlim([0,npops])
   ylim([0,npops])
   axh.xaxis.set_ticks(array(range(npops))+0.5)
   axh.xaxis.set_ticklabels(ctyp_names)
   axh.xaxis.set_ticks_position('top')
   axh.yaxis.set_ticks(array(range(npops))+0.5)
   axh.yaxis.set_ticklabels(rnames)
   tmp.cmap = cm.bone_r # Blues_r # summer_r#cm.gist_earth # The _r reverses the order -- beautiful solution!
   colorbar(tmp, fraction=0.2)
   xlabel('Pre-Synaptic Cell Type')
   ylabel('Post-Synaptic Cell Type')

def plot_inter_ctyp_pmat ():
   def cellpops ():
      I2L=0; I2=1;  E2=2; I5L=3; I5=4; E5R=5; E5B=6; I6L=7 ;I6=8; E6=9;
      return I2L, I2, E2, I5L, I5, E5R, E5B, I6L, I6, E6

   # Define connectivity matrix (copied out of network.py)
   def definepmat (npops):
      I2L, I2, E2, I5L, I5, E5R, E5B, I6L, I6, E6 = cellpops()
      pmat = zeros((npops,npops))

      pmat[E2][E2]=0.2     # weak by wiring matrix in (Weiler et al., 2008)
      pmat[E2][E5B]=0.8    # strong by wiring matrix in (Weiler et al., 2008)
      pmat[E2][E5R]=0.8    # strong by wiring matrix in (Weiler et al., 2008)
      pmat[E2][I5L]=0.51   # L2/3 E -> L5 LTS I (justified by (Apicella et al., 2012)
      pmat[E2][E6]=0       # none by wiring matrix in (Weiler et al., 2008)
      pmat[E2][I2L]=0.51
      pmat[E2][I2]=0.43

      pmat[E5B][E2]=0          # none by wiring matrix in (Weiler et al., 2008)
      pmat[E5B][E5B]=0.04 * 4  # set using (Kiritani et al., 2012) Fig. 6D, Table 1, value x 4 
      pmat[E5B][E5R]=0         # set using (Kiritani et al., 2012) Fig. 6D, Table 1, value x 4 
      pmat[E5B][E6]=0          # none by suggestion of Ben and Gordon over phone
      pmat[E5B][I5L]=0         # ruled out by (Apicella et al., 2012) Fig. 7
      pmat[E5B][I5]=0.43 

      pmat[E5R][E2]=0.2        # weak by wiring matrix in (Weiler et al., 2008)
      pmat[E5R][E5B]=0.21 * 4  # set using (Kiritani et al., 2012) Fig. 6D, Table 1, value x 4
      pmat[E5R][E5R]=0.11 * 4  # set using (Kiritani et al., 2012) Fig. 6D, Table 1, value x 4 
      pmat[E5R][E6]=0.2        # weak by wiring matrix in (Weiler et al., 2008)
      pmat[E5R][I5L]=0         # ruled out by (Apicella et al., 2012) Fig. 7
      pmat[E5R][I5]=0.43

      pmat[E6][E2]=0           # none by wiring matrix in (Weiler et al., 2008)
      pmat[E6][E5B]=0.2        # weak by wiring matrix in (Weiler et al., 2008)
      pmat[E6][E5R]=0.2        # weak by wiring matrix in (Weiler et al., 2008)
      pmat[E6][E6]=0.2         # weak by wiring matrix in (Weiler et al., 2008)
      pmat[E6][I6L]=0.51
      pmat[E6][I6]=0.43

      pmat[I2L][E2]=0.35
      pmat[I2L][I2L]=0.09
      pmat[I2L][I2]=0.53
      pmat[I2][E2]=0.44
      pmat[I2][I2L]=0.34
      pmat[I2][I2]=0.62
      pmat[I5L][E5B]=0.35
      pmat[I5L][E5R]=0.35
      pmat[I5L][I5L]=0.09
      pmat[I5L][I5]=0.53
      pmat[I5][E5B]=0.44
      pmat[I5][E5R]=0.44
      pmat[I5][I5L]=0.34
      pmat[I5][I5]=0.62
      pmat[I6L][E6]=0.35
      pmat[I6L][I6L]=0.09
      pmat[I6L][I6]=0.53
      pmat[I6][E6]=0.44
      pmat[I6][I6L]=0.34
      pmat[I6][I6]=0.62

      return pmat

   ctyp_strs = ["I2L", "I2", "E2", "I5L", "I5", "STR", "SPI", "I6L", "I6", "E6"]
   pmat = definepmat(len(ctyp_strs))
   plot_inter_ctyp_matrix(ctyp_strs, pmat)

def plot_inter_ctyp_avg_wts ():
   def cellpops ():
      I2L=0; I2=1;  E2=2; I5L=3; I5=4; E5R=5; E5B=6; I6L=7 ;I6=8; E6=9;
      return I2L, I2, E2, I5L, I5, E5R, E5B, I6L, I6, E6

   I2L, I2, E2, I5L, I5, E5R, E5B, I6L, I6, E6 = cellpops()
   ctyp_labels = ["I2L", "I2", "E2", "I5L", "I5", "STR", "SPI", "I6L", "I6", "E6"]
   ctyp_strs = ["I2L", "I2", "E2", "I5L", "I5", "E5R", "E5B", "I6L", "I6", "E6"]
   wmat = zeros((10,10))

   # still under construction...
   for ii in range(10):
      for jj in range(10):
         wmat[ii][jj] = get_avg_conv_wts(ctyp_strs[ii], ctyp_strs[jj])
   
   plot_inter_ctyp_matrix(ctyp_labels, wmat)
   
#
# Batch functions
#

def batchA (simdatadir='data/12oct11_sims1',savebatch=True):
   # Set up ctyp range to use (useful to make smaller sometimes due to memory problems).
   ctyp_list = ('E2','I2','I2L','E5R','E5B','I5','I5L','E6','I6','I6L')
#   ctyp_list = ('E2','I2','I2L')
#   ctyp_list = ('E5R','E5B','I5','I5L')
#   ctyp_list = ('E6','I6','I6L')
#   ctyp_list = ('E5R','E5B','I5','I5L','E6','I6','I6L')

   # Set up upper bound for the firing response rate calculations.
   short_resp_dur = 10   # ms
   long_resp_dur = 1000  # ms
   short_resp_end = 1002 + short_resp_dur  # start is always at 1002 ms so excluding shock
   long_resp_end = 1002 + long_resp_dur
   if (short_resp_end > 2000):
      short_resp_end = 2000
   if (long_resp_end > 2000):
      long_resp_end = 2000

   for ctyp in ctyp_list:
      exec('global shock_%s_spk_counts' % ctyp)
      exec('global shock_%s_resp_rates' % ctyp)
      exec('global shock_%s_resp_rates2' % ctyp)

   # For each of the shock percentages...
   for ii in range(len(shock_pcts)):
      # Load L2/3 shock data.
      outfilestem = '%s/L2shock%dp' % (simdatadir, shock_pcts[ii])
#      outfilestem = '%s/L2isolshock%dp' % (simdatadir, shock_pcts[ii])
      loadsimbyfstem(outfilestem)
      for ctyp in ctyp_list:
         exec("shock_%s_spk_counts[ii,0] = get_ctyp_spk_count('%s',1000,1001)" % \
            (ctyp,ctyp))
         exec("shock_%s_resp_rates[ii,0] = get_ctyp_fire_rate('%s',1002,%d)" % \
            (ctyp,ctyp,short_resp_end))
         exec("shock_%s_resp_rates2[ii,0] = get_ctyp_fire_rate('%s',1002,%d)" % \
            (ctyp,ctyp,long_resp_end))

      # Load L5A shock data.
      outfilestem = '%s/L5Ashock%dp' % (simdatadir, shock_pcts[ii])
#      outfilestem = '%s/L5Aisolshock%dp' % (simdatadir, shock_pcts[ii])
      loadsimbyfstem(outfilestem)
      for ctyp in ctyp_list:
         exec("shock_%s_spk_counts[ii,1] = get_ctyp_spk_count('%s',1000,1001)" % \
            (ctyp,ctyp))
         exec("shock_%s_resp_rates[ii,1] = get_ctyp_fire_rate('%s',1002,%d)" % \
            (ctyp,ctyp,short_resp_end))
         exec("shock_%s_resp_rates2[ii,1] = get_ctyp_fire_rate('%s',1002,%d)" % \
            (ctyp,ctyp,long_resp_end))

      # Load L5B shock data.
      outfilestem = '%s/L5Bshock%dp' % (simdatadir, shock_pcts[ii])
#      outfilestem = '%s/L5Bisolshock%dp' % (simdatadir, shock_pcts[ii])
      loadsimbyfstem(outfilestem)
      for ctyp in ctyp_list:
         exec("shock_%s_spk_counts[ii,2] = get_ctyp_spk_count('%s',1000,1001)" % \
            (ctyp,ctyp))
         exec("shock_%s_resp_rates[ii,2] = get_ctyp_fire_rate('%s',1002,%d)" % \
            (ctyp,ctyp,short_resp_end))
         exec("shock_%s_resp_rates2[ii,2] = get_ctyp_fire_rate('%s',1002,%d)" % \
            (ctyp,ctyp,long_resp_end))

      # Load L6 shock data.
      outfilestem = '%s/L6shock%dp' % (simdatadir, shock_pcts[ii])
#      outfilestem = '%s/L6isolshock%dp' % (simdatadir, shock_pcts[ii])
      loadsimbyfstem(outfilestem)
      for ctyp in ctyp_list:
         exec("shock_%s_spk_counts[ii,3] = get_ctyp_spk_count('%s',1000,1001)" % \
            (ctyp,ctyp))
         exec("shock_%s_resp_rates[ii,3] = get_ctyp_fire_rate('%s',1002,%d)" % \
            (ctyp,ctyp,short_resp_end))
         exec("shock_%s_resp_rates2[ii,3] = get_ctyp_fire_rate('%s',1002,%d)" % \
            (ctyp,ctyp,long_resp_end))

      # Save the data if we want to.
      if (savebatch):
         dirname = '%s/batch_results' % simdatadir

         # If the directory doesn't exist, make it.
         if (not os.path.exists(dirname)):
            os.system('mkdir %s' % dirname)

         # Save all of the numpy arrays.
         for ctyp in ctyp_list:
            exec("fname = '%s/%s-spk_counts.txt'" % (dirname, ctyp))
            exec("savetxt('%s', shock_%s_spk_counts)" % (fname, ctyp))
            exec("fname = '%s/%s-resp_rates.txt'" % (dirname, ctyp))
            exec("savetxt('%s', shock_%s_resp_rates)" % (fname, ctyp))
            exec("fname = '%s/%s-resp_rates2.txt'" % (dirname, ctyp))
            exec("savetxt('%s', shock_%s_resp_rates2)" % (fname, ctyp))            

def loadbatchA (simdatadir='data/12oct11_sims1'):
   # This block doesn't seem to work for this function.
#   for ctyp in ('E2','I2','I2L','E5R','E5B','I5','I5L','E6','I6','I6L'):
#      exec('global shock_%s_spk_counts' % ctyp)
#      exec('global shock_%s_resp_rates' % ctyp)
#      exec('global shock_%s_resp_rates2' % ctyp)
   global shock_E2_spk_counts
   global shock_E2_resp_rates
   global shock_E2_resp_rates2
   global shock_I2_spk_counts
   global shock_I2_resp_rates
   global shock_I2_resp_rates2
   global shock_I2L_spk_counts
   global shock_I2L_resp_rates
   global shock_I2L_resp_rates2
   global shock_E5R_spk_counts
   global shock_E5R_resp_rates
   global shock_E5R_resp_rates2
   global shock_E5B_spk_counts
   global shock_E5B_resp_rates
   global shock_E5B_resp_rates2
   global shock_I5_spk_counts
   global shock_I5_resp_rates
   global shock_I5_resp_rates2
   global shock_I5L_spk_counts
   global shock_I5L_resp_rates
   global shock_I5L_resp_rates2
   global shock_E6_spk_counts
   global shock_E6_resp_rates
   global shock_E6_resp_rates2
   global shock_I6_spk_counts
   global shock_I6_resp_rates
   global shock_I6_resp_rates2
   global shock_I6L_spk_counts
   global shock_I6L_resp_rates
   global shock_I6L_resp_rates2

   dirname = '%s/batch_results' % simdatadir

   # Load all of the numpy arrays.
   # This is not working in this function for some reason.  Seems like the for loop
   # is its own scope and the information loaded in it is lost.
#   for ctyp in ('E2','I2','I2L','E5R','E5B','I5','I5L','E6','I6','I6L'):
#      exec("fname = '%s/%s-spk_counts.txt'" % (dirname, ctyp))
#      exec("shock_%s_spk_counts = loadtxt('%s')" % (ctyp, fname))
#      exec("fname = '%s/%s-resp_rates.txt'" % (dirname, ctyp))
#      exec("shock_%s_resp_rates = loadtxt('%s')" % (ctyp, fname))
#      exec("fname = '%s/%s-resp_rates2.txt'" % (dirname, ctyp))
#      exec("shock_%s_resp_rates2 = loadtxt('%s')" % (ctyp, fname))  
   fname = '%s/E2-spk_counts.txt' % dirname
   shock_E2_spk_counts = loadtxt(fname)
   fname = '%s/E2-resp_rates.txt' % dirname
   shock_E2_resp_rates = loadtxt(fname)
   fname = '%s/E2-resp_rates2.txt' % dirname
   shock_E2_resp_rates2 = loadtxt(fname)

   fname = '%s/I2-spk_counts.txt' % dirname
   shock_I2_spk_counts = loadtxt(fname)
   fname = '%s/I2-resp_rates.txt' % dirname
   shock_I2_resp_rates = loadtxt(fname)
   fname = '%s/I2-resp_rates2.txt' % dirname
   shock_I2_resp_rates2 = loadtxt(fname)

   fname = '%s/I2L-spk_counts.txt' % dirname
   shock_I2L_spk_counts = loadtxt(fname)
   fname = '%s/I2L-resp_rates.txt' % dirname
   shock_I2L_resp_rates = loadtxt(fname)
   fname = '%s/I2L-resp_rates2.txt' % dirname
   shock_I2L_resp_rates2 = loadtxt(fname)

   fname = '%s/E5R-spk_counts.txt' % dirname
   shock_E5R_spk_counts = loadtxt(fname)
   fname = '%s/E5R-resp_rates.txt' % dirname
   shock_E5R_resp_rates = loadtxt(fname)
   fname = '%s/E5R-resp_rates2.txt' % dirname
   shock_E5R_resp_rates2 = loadtxt(fname)

   fname = '%s/E5B-spk_counts.txt' % dirname
   shock_E5B_spk_counts = loadtxt(fname)
   fname = '%s/E5B-resp_rates.txt' % dirname
   shock_E5B_resp_rates = loadtxt(fname)
   fname = '%s/E5B-resp_rates2.txt' % dirname
   shock_E5B_resp_rates2 = loadtxt(fname)

   fname = '%s/I5-spk_counts.txt' % dirname
   shock_I5_spk_counts = loadtxt(fname)
   fname = '%s/I5-resp_rates.txt' % dirname
   shock_I5_resp_rates = loadtxt(fname)
   fname = '%s/I5-resp_rates2.txt' % dirname
   shock_I5_resp_rates2 = loadtxt(fname)

   fname = '%s/I5L-spk_counts.txt' % dirname
   shock_I5L_spk_counts = loadtxt(fname)
   fname = '%s/I5L-resp_rates.txt' % dirname
   shock_I5L_resp_rates = loadtxt(fname)
   fname = '%s/I5L-resp_rates2.txt' % dirname
   shock_I5L_resp_rates2 = loadtxt(fname)

   fname = '%s/E6-spk_counts.txt' % dirname
   shock_E6_spk_counts = loadtxt(fname)
   fname = '%s/E6-resp_rates.txt' % dirname
   shock_E6_resp_rates = loadtxt(fname)
   fname = '%s/E6-resp_rates2.txt' % dirname
   shock_E6_resp_rates2 = loadtxt(fname)

   fname = '%s/I6-spk_counts.txt' % dirname
   shock_I6_spk_counts = loadtxt(fname)
   fname = '%s/I6-resp_rates.txt' % dirname
   shock_I6_resp_rates = loadtxt(fname)
   fname = '%s/I6-resp_rates2.txt' % dirname
   shock_I6_resp_rates2 = loadtxt(fname)

   fname = '%s/I6L-spk_counts.txt' % dirname
   shock_I6L_spk_counts = loadtxt(fname)
   fname = '%s/I6L-resp_rates.txt' % dirname
   shock_I6L_resp_rates = loadtxt(fname)
   fname = '%s/I6L-resp_rates2.txt' % dirname
   shock_I6L_resp_rates2 = loadtxt(fname)

def batchB (simdatadir='data/12oct10_sims2',savebatch=True):
   # Set up ctyp range to use (useful to make smaller sometimes due to memory problems).
   ctyp_list = ('E2','I2','I2L','E5R','E5B','I5','I5L','E6','I6','I6L')

   for ctyp in ctyp_list:
      exec('global shock_%s_spk_counts' % ctyp)
      exec('global shock_%s_resp_rates' % ctyp)
      exec('global shock_%s_resp_rates2' % ctyp)
      exec('global shock_%s_resp_rates_mean' % ctyp)
      exec('global shock_%s_resp_rates_std' % ctyp)
      exec('global shock_%s_resp_rates2_mean' % ctyp)
      exec('global shock_%s_resp_rates2_std' % ctyp)

   # For each of the shock percentages...
   for ii in range(len(shock_pcts)):
      # Load L2/3 shock data.
      outfilestem = '%s/L2shock%dp' % (simdatadir, shock_pcts[ii])
      loadsimbyfstem(outfilestem)
      for ctyp in ctyp_list:
         exec("shock_%s_spk_counts[ii,0] = get_ctyp_spk_count('%s',1000,1001)" % \
            (ctyp,ctyp))
         exec("shock_%s_resp_rates[ii,0] = get_ctyp_fire_rate('%s',1002,1012)" % \
            (ctyp,ctyp))
         exec("shock_%s_resp_rates2[ii,0] = get_ctyp_fire_rate('%s',1002,2000)" % \
            (ctyp,ctyp))

      # Load L5A shock data.
      outfilestem = '%s/L5Ashock%dp' % (simdatadir, shock_pcts[ii])
      loadsimbyfstem(outfilestem)
      for ctyp in ctyp_list:
         exec("shock_%s_spk_counts[ii,1] = get_ctyp_spk_count('%s',1000,1001)" % \
            (ctyp,ctyp))
         exec("shock_%s_resp_rates[ii,1] = get_ctyp_fire_rate('%s',1002,1012)" % \
            (ctyp,ctyp))
         exec("shock_%s_resp_rates2[ii,1] = get_ctyp_fire_rate('%s',1002,2000)" % \
            (ctyp,ctyp))

      # Load L5B shock data.
      outfilestem = '%s/L5Bshock%dp' % (simdatadir, shock_pcts[ii])
      loadsimbyfstem(outfilestem)
      for ctyp in ctyp_list:
         exec("shock_%s_spk_counts[ii,2] = get_ctyp_spk_count('%s',1000,1001)" % \
            (ctyp,ctyp))
         exec("shock_%s_resp_rates[ii,2] = get_ctyp_fire_rate('%s',1002,1012)" % \
            (ctyp,ctyp))
         exec("shock_%s_resp_rates2[ii,2] = get_ctyp_fire_rate('%s',1002,2000)" % \
            (ctyp,ctyp))

      # Load L6 shock data.
      outfilestem = '%s/L6shock%dp' % (simdatadir, shock_pcts[ii])
      loadsimbyfstem(outfilestem)
      for ctyp in ctyp_list:
         exec("shock_%s_spk_counts[ii,3] = get_ctyp_spk_count('%s',1000,1001)" % \
            (ctyp,ctyp))
         exec("shock_%s_resp_rates[ii,3] = get_ctyp_fire_rate('%s',1002,1012)" % \
            (ctyp,ctyp))
         exec("shock_%s_resp_rates2[ii,3] = get_ctyp_fire_rate('%s',1002,2000)" % \
            (ctyp,ctyp))

def rasterdemo ():
   plot_raster(980,1060,xmarg=0,ymarg=0)
   xticks((0.96,0.98,1.0,1.02,1.04,1.06,1.08),('','-20','0','20','40','60',''))
   xlim(0.98,1.06)
   xlabel('Peri-stimulus Time (ms)')
   ylabel('')
   yticks([])

#
# Script code (run always)
#

# Load the sys and os packages.
import sys, os

# Set up the file stem pointing to the simulation output files.
# NOTE: When "ipython -pylab load.py" is called, only 1 arg is loaded, the last.
if (len(sys.argv) > 1):
   outfilestem = sys.argv[1]
else:
   outfilestem = './runsim'

# Load the interpreter object.
from neuron import h

# Load the NEURON GUI.
if use_NEURON_GUI:
   from neuron import gui

# Load the numpy namespace.
from numpy import *

# Load the neuroplot namespace.
from neuroplot import *

# Set up neurosim run-neuron code.
h.xopen("setup.hoc")
h.xopen("nrnoc.hoc")
h.load_file("init.hoc")

# Load functions for expediting NQS processing in Python interpreter.
execfile("miscfuncs.py")

# If the NEURON GUI is up, load the E5R and E5B cell data from the saved grvec file.
#if use_NEURON_GUI:
#   ldgrvec('celldata.gvc')

# Set up hoc objects for NQS tables for loadsimbyfstem()
h('objref nqLFP')
h('objref snq')
h('objref cellsnq')
h('objref connsnq')
h('objref connsnq2')

# Load the first simulation.
loadsimbyfstem(outfilestem,grvec_load=False)

#
# Set layer information.
#

# Set the M1 cell density array.
M1dens = [0.0044,0.0297,0.0253,0.0121,0.0141,0.0106,0.0195,0.0372,0.0361,0.0591,0.0748,0.0835,0.0947,0.0793,0.0923,0.0796,0.0702,0.0809,0.0749,0.0827,0.0792,0.0893,0.0475,0.0811,0.0755,0.0754,0.083,0.0688,0.0892,0.0765,0.1058,0.0791,0.091,0.1018,0.0931,0.0847,0.0919,0.0661,0.0896,0.0686,0.073,0.0774,0.0588,0.0661,0.0707,0.0496,0.0714,0.065,0.0463,0.0766,0.0546,0.0619,0.0584,0.0541,0.0612,0.047,0.069,0.0611,0.0519,0.0601,0.0646,0.0459,0.0634,0.0348,0.0705,0.0546,0.0485,0.0434,0.0659,0.0585,0.0655,0.0638,0.0922,0.0753,0.0833,0.0798,0.0751,0.0831,0.0893,0.1004,0.0951,0.0734,0.0845,0.0918,0.093,0.0867,0.0828,0.0971,0.088,0.0822,0.0925,0.0884,0.0944,0.0892,0.0659,0.0875,0.1046,0.123,0.1187,0.0656]
M1dens = array(M1dens)

# Set up the layer left bin boundary array.
lyrbins = arange(0,1,0.01)
lyrbins *= 1412.0

# Remember the simulation duration (in ms).
sim_duration = 2000.0

#
# Batch trend collection
#

# data for usual (old batch)
#shock_pcts = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,30,35,\
#   40,45,50,51,52,53,54,55]
# data for new batch
shock_pcts = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,30,35,\
   40,45,50,51,52,53,54,55,60,65,70,75,80,85,90,95,100]
# data for isolated layer models batch
#shock_pcts = [5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100]
inseeds = [1,2,3,4,5]
wseeds = [1,2,3,4,5]

# Set up the counts for shock spikes and responses.
for ctyp in ('E2','I2','I2L','E5R','E5B','I5','I5L','E6','I6','I6L'):
   # Batch A (single-seed version)
   exec('shock_%s_spk_counts = zeros((%d,4))' % (ctyp, len(shock_pcts)))
   exec('shock_%s_resp_rates = zeros((%d,4))' % (ctyp, len(shock_pcts)))
   exec('shock_%s_resp_rates2 = zeros((%d,4))' % (ctyp, len(shock_pcts)))

   # Batch B (multi-seed version)
#   exec('shock_%s_spk_counts = zeros((%d,4,%d,%d))' % \
#      (ctyp, len(shock_pcts), len(wseeds), len(inseeds)))
#   exec('shock_%s_resp_rates = zeros((%d,4,%d,%d))' % \
#      (ctyp, len(shock_pcts), len(wseeds), len(inseeds)))
#   exec('shock_%s_resp_rates2 = zeros((%d,4,%d,%d))' % \
#      (ctyp, len(shock_pcts), len(wseeds), len(inseeds)))
#   exec('shock_%s_resp_rates_mean = zeros((%d,4))' % (ctyp, len(shock_pcts)))
#   exec('shock_%s_resp_rates_std = zeros((%d,4))' % (ctyp, len(shock_pcts)))
#   exec('shock_%s_resp_rates2_mean = zeros((%d,4))' % (ctyp, len(shock_pcts)))
#   exec('shock_%s_resp_rates2_std = zeros((%d,4))' % (ctyp, len(shock_pcts)))

from numpy import *
from matplotlib.pyplot import *
from matplotlib.mlab import *
rasterdemo()
