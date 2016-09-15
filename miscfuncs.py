# miscfuncs.py -- miscellaneous Python functions to be callable from the 
#   running interpreter

# Various functions to be callable from the Python interpreter.
# This is kind of the dumping place for new function to be more properly 
# place later or functions that don't fit elsewhere, but are helpful.
#
# Note: This code should not make calls to Python functions in the scope 
# of the outer Python scope this Python / hoc code mixture uses.
#
# Also, this needs to run in the context of numpy namespace being loaded.
#
# Last update: 10/8/12 (georgec)

#
# Functions to ease interactive interface with hoc code
#

# Convert a hoc variable (passed as string) into a numpy array.
def hv2narr (hocvar='vec'):
   exec('%s = h.%s.to_python()' % (hocvar, hocvar))
   exec('%s = array(%s)' % (hocvar, hocvar))
   exec('x = %s' % hocvar)
   return x

# Convert a numpy array into a hoc variable (passed as string).
def narr2hv (hocvar,narr):
   h('objref %s' % hocvar)
   h('%s = new Vector()' % hocvar)
   h('objref tmp')
   h.tmp = narr
   h('{%s.from_python(tmp)}' % hocvar)

# Do a gethdrs() for an NQS table.
def shownqshdr (nqsvar='col[0].cellsnq'):
   print nqsvar
   h('%s.gethdrs()' % nqsvar)
 
# Do a pr(numrows) for an NQS table.
def nqspr (nqsvar='col[0].cellsnq', numrows=10):
   print nqsvar    
   h('%s.pr(%d)' % (nqsvar, numrows))
   
# Convert a hoc NQS column into a numpy array.
def nqscol2narr (nqsvar='col[0].cellsnq', colstr='col'):
   h('objref tmpv')
   h('tmpv = new Vector()')
   h('tmpv = %s.getcol("%s")' % (nqsvar, colstr))
   tmpv = h.tmpv.to_python()
   tmpv = array(tmpv)
   return tmpv
   
# Get the CTYP number from the string (e.g. 'E2').
def get_ctyp_num (ctyp_str):
   ctyp_ind = -1
   for ii in range(int(h.CTYPi)):
      if (h.CTYP.o(ii).s == ctyp_str):
         ctyp_ind = ii
   return ctyp_ind

# Get the CTYP string from the CTYP number.
def get_ctyp_str (ctyp_num):
   return h.CTYP.o(ctyp_num).s

# Get the STYP number from the string (e.g. 'AM2').
def get_styp_num (styp_str):
   styp_ind = -1
   for ii in range(int(h.STYPi)):
      if (h.STYP.o(ii).s == styp_str):
         styp_ind = ii
   return styp_ind

# Get the STYP string from the STYP number.
def get_styp_str (styp_num):
   return h.STYP.o(styp_num).s

# Get the average intracolumnar weight of a certain type between two
# cell types.
def get_ave_conn_wt (fromtype, totype, syntype, thecol=0):
   print 'wmat entry: %f' % \
      h.wmat[get_ctyp_num(fromtype)][get_ctyp_num(totype)][get_styp_num(syntype)][0]
   if (h.wsetting_INTF6):
      h.col[thecol].connsnq.verbose = 0
      h.col[thecol].selectconns2(get_ctyp_num(fromtype),get_ctyp_num(totype))
      if (syntype in ['AM','AM2','GA','GA2']):
         h('vec = col[%d].connsnq.getcol("wt1")' % thecol)
      else:
         h('vec = col[%d].connsnq.getcol("wt2")' % thecol)
      print 'actual weight average from connsnq: %f' % hv2narr('vec').mean()
      h('col[%d].connsnq.tog("db")' % thecol)
      h.col[thecol].connsnq.verbose = 1
    
#
# Parameter setting / reading functions (not in params.hoc)
#

#
#    vq input reading / setting functions (not usable for NetStim method)
#

# Show each of the input weights for each used cell type / synapse type
def show_ctyp_input_wts (thecol=0):
   ctyplist = ['E2','I2','I2L','E5R','E5B','I5','I5L','E6','I6','I6L']
#   ctyplist = ['E2','I2','I2L','E4','I4','I4L','E5R','E5B','I5','I5L', \
#      'E6','I6','I6L']
   styplist = ['AM2','NM2','GA','GA2']

   # Set up a hoc vector. 
   h('objref tmpv')
   h('tmpv = new Vector()')
   
   # Turn off verbosity of table.
   h('lcstim.o(%d).vq.verbose = 0' % thecol)
   
   # Print the table header.
   print 'ctyp\tstyp\twt'
   
   # For each of the cell types...
   for ctyp in ctyplist:
      # Find the first cell ID of this type. 
      h('tmpv = col[%d].cellsnq.getrow("celltype",%s)' % (thecol, ctyp))
      cell_id = int(h.tmpv[1])
      
      # For all of the used synapse types...
      for styp in styplist:
         h('{lcstim.o(%d).vq.select("ind",%d,"sy",%d)}' % \
            (thecol, cell_id, get_styp_num(styp)))
         h('tmpv = lcstim.o(%d).vq.getrow(0)' % thecol)         
         print ctyp, '\t', styp, '\t',
         print h.tmpv[2]
         
   # Reset the vq table
   h('lcstim.o(%d).vq.tog("DB")' % thecol)
   h('lcstim.o(%d).vq.verbose = 1' % thecol)

# Set all weights of a certain cell type / synapse type to newwt.
def set_ctyp_input_wts (ctyp='E2', styp='AM2', newwt=3.75):
   # Set up a hoc vector. 
   h('objref tmpv')
   h('tmpv = new Vector()')    

   # Turn off verbosity of table.
   thecol = 0
   h('lcstim.o(%d).vq.verbose = 0' % thecol)
   
   # Find the first cell ID of this type. 
   h('tmpv = col[%d].cellsnq.getrow("celltype",%s)' % (thecol, ctyp))
   cell_id = int(h.tmpv[1])

   # Get the old weight for the type.
   h('{lcstim.o(%d).vq.select("ind",%d,"sy",%d)}' % \
      (thecol, cell_id, get_styp_num(styp)))
   h('tmpv = lcstim.o(%d).vq.getrow(0)' % thecol)  
   oldwt = h.tmpv[2]
   
   # For all of the columns...
   for ii in range(int(h.numcols)):
      h('lcstim.o(%d).mulwts(%s, %s, %f)' % (ii, ctyp, styp, (newwt / oldwt)))
         
   # Reset the vq table
   h('lcstim.o(%d).vq.tog("DB")' % thecol)
   h('lcstim.o(%d).vq.verbose = 1' % thecol)
   
# Reset all weights of a certain cell type / synapse type to a multiplied 
# version of their original value.
def resetmult_ctyp_input_wts (ctyp='E2', styp='AM2', wtmult=1.0):
   # Reset the CSTIMs to their original value. 
   h.setcstim()

   # For all of the columns...
   for ii in range(int(h.numcols)):
      h('lcstim.o(%d).mulwts(%s, %s, %f)' % (ii, ctyp, styp, wtmult))
 
# Set all weights of a certain cell type / synapse type to a multiplied 
# version of their original value. 
def mult_ctyp_input_wts (ctyp='E2', styp='AM2', wtmult=1.0):
   # For all of the columns...
   for ii in range(int(h.numcols)):
      h('lcstim.o(%d).mulwts(%s, %s, %f)' % (ii, ctyp, styp, wtmult))
         
#
# grvec functions
#

# Load saved grvec simulation info.
# NOTE: Both the file and its dot-prefixed equivalent need to be present for gvnew() to succeed.
def ldgrvec (grvec_fname):
   h.gvnew(grvec_fname)

# Look at the grvec printlist for the current sim or a saved grvec file.
def lookgrveclist (gvcobjnum=0): 
   if (gvcobjnum == 0):
      gvcstr = 'current simulation'
   else:
      gvcstr = h.panobjl.o(gvcobjnum).filename
   print 'GRVEC List #%d (%s)' % (gvcobjnum, gvcstr)
   print '--------------------------------------------'
   if (gvcobjnum == 0):
      for ii in range(int(h.printlist.count())):
         prstr = '%d %s %d ' % (ii, h.printlist.o(ii).name, \
            h.printlist.o(ii).vec.size())
         print prstr,
         if (h.printlist.o(ii).tvec == None):
            print '(vec only)'
         else:
            print '\n',
   else:
      for ii in range(int(h.panobjl.o(gvcobjnum).printlist.count())):
         print '%d %s %d' % (ii, h.panobjl.o(gvcobjnum).printlist.o(ii).name, \
            h.panobjl.o(gvcobjnum).printlist.o(ii).size)

# Get tvec and vec (in numpy form) from grvec (gvcobjnum=0 means current 
# sim; >0 means saved grvec file)
def getgrvecdat (gvcobjnum=0, vecname='C0_X0_Y0_SPKS'):
   found_ind = -1
   for ii in range(int(h.panobjl.o(gvcobjnum).printlist.count())):
       if (h.panobjl.o(gvcobjnum).printlist.o(ii).name == vecname):
          found_ind = ii
   if (found_ind == -1):
      print 'No such array is on the printlist.'
      return None, None
   elif (gvcobjnum == 0):
      vec = h.printlist.o(found_ind).vec.to_python()    
      tvec = h.printlist.o(found_ind).tvec
      if (tvec == None):
         tvec = linspace(0,len(vec)-1,len(vec))
      else:
         tvec = array(tvec.to_python())
   else:
      h('goodread = panobjl.o(%d).rv_readvec(%d,tvec,vec)' % (gvcobjnum, found_ind))
      if (not h.goodread):
         h('panobjl.o(%d).rv_readvec(%d,vec)' % (gvcobjnum, found_ind))
         vec = h.vec.to_python()
         tvec = linspace(0,len(vec)-1,len(vec))
      else:         
         tvec = array(h.tvec.to_python())
         vec = h.vec.to_python()
   vec = array(vec)
   return tvec, vec

# Plot tvec and vec (in numpy form) from grvec (gvcobjnum=0 means current 
# sim; >0 means saved grvec file)
def plotgrvecdat (gvcobjnum=0, vecname='C0_X0_Y0_SPKS'):
   tvec,vec = getgrvecdat(gvcobjnum, vecname)
   if (tvec != None):
      plot(tvec,vec)

#
# Miscellaneous functions
#
