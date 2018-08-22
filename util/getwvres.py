# This spice wv reader for *.sw* only [support binary and acii]
from ctypes import *
from collections import Iterable, OrderedDict
import os,sys,glob,pickle,json
from copy import deepcopy
import numpy as np
import platform
import pdb

POST2001_WIDTH=11
relTolSP=0.005
relTolSP_meas=0.01
DELTA=1e-50
minVal=1e-10
HSPICE_END_SYMB='0.10000E+31'
HSPICE_END_SYMB_WIN='0.1000E+031'
DEBUG= not True
NdataOneSweep = 61  #data points in one sweep curve
VSweep=[]
MAXCHARONELINE = 80 #spice wv max char in header line 78 for data block
MAXERROR=1e9



class SpiceReader(object):
    @classmethod
    def get_array_fix_width(cls, wvfile):
        """
        read a ascii wv files (sw0/ac0/tr0) to an array whose each element has the fixed width
        :param wvfile:
        :return:
        """
        global HSPICE_END_SYMB
        sys=platform.system()
        if sys[:6] == 'CYGWIN' or sys[:7] == 'Windows':
            HSPICE_END_SYMB = HSPICE_END_SYMB_WIN
        _is_data_block = False
        all_data_array = []
        all_keys = []
        all_keys_content = ''
        find_VOLTS = False
        find_HERTZ = False
        find_TIME = False
        insert_space = True
        with open(wvfile) as fin:
            for idx,line in enumerate(fin):
                insert_space = True
                if (not _is_data_block) and line[MAXCHARONELINE-1] != ' ':
                    insert_space = False
                line = line.rstrip()  #hspice wv file has placed space between signals
                if idx==2:
                    nGroup = int(line.split()[-1])
                    print 'nGroup=',nGroup,wvfile
                    assert nGroup >=0
                if 'VOLTS' in line:
                    find_VOLTS = True
                if 'HERTZ' in line:
                    find_HERTZ = True
                if 'TIME' in line:
                    find_TIME = True
                if (find_VOLTS or find_HERTZ or find_TIME) and not _is_data_block:
                    all_keys_content += line+' ' if insert_space else line
                    #print '**dbg line len=',len(line) ,line
                    #print '**dbg toto-line:', all_keys_content, insert_space
                if '$&%#' in line:
                    _is_data_block = True
                    #print '**dbg found end symbol!'
                    continue
                if _is_data_block:
                    num = len(line)/POST2001_WIDTH
                    vec=[line[i*POST2001_WIDTH:(i+1)*POST2001_WIDTH] for i in range(num)]
                    if HSPICE_END_SYMB in vec:
                        #remove the end symbol of data block in wv files
                        vec.pop(vec.index(HSPICE_END_SYMB))
                    all_data_array.extend(vec)
        all_keys = all_keys_content.split()[:-1]
        if find_VOLTS:
            all_keys = all_keys[all_keys.index('VOLTS'):]
        elif find_HERTZ:
            all_keys = all_keys[all_keys.index('HERTZ'):]
        elif find_TIME:
            all_keys = all_keys[all_keys.index('TIME'):]
        else:
            print 'Error: No DC/AC/TR symbol found!'
            sys.exit(1)
        #print '**dbg all_keys=', all_keys,all_data_array
        Nblock=len(all_keys)
        NallData=len(all_data_array)
        #pdb.set_trace()
        if nGroup>0:
            # To handle double sweep wv files (only for DC now)
            # Tips: the data block is stored as below
            # vsweep2 vsweep1 <Y item list>
            # here we move vsweep2 into VSweep from data block vec in `all_data_array
            # and make sure NallData%Nblock == 0 to check the grouping method is valid
            assert find_VOLTS and not find_HERTZ, 'only support DC with double sweep!'
            all_keys = all_keys[:-1]
            Nblock -= 1
            if DEBUG:
                print 'N=',Nblock,NallData,nGroup,NdataOneSweep
            for idx in range(nGroup):
                idx_rm = idx*(NdataOneSweep*Nblock+1)
                if DEBUG:
                    print 'vseep-i',all_data_array[idx_rm]
                VSweep.append(all_data_array[idx_rm])
            for idx in range(nGroup):
                idx_rm = idx*(NdataOneSweep*Nblock+1)
                if DEBUG:
                    print 'idx_rm=',idx_rm
                all_data_array.pop(idx_rm)
            NallData = len(all_data_array)
        if find_HERTZ:
            # Handle real/img part of AC signal data in *.ac# file
            all_keys_expend = []
            for it in all_keys:
                if 'HERTZ' == it or it.lower().startswith('vdb(') \
                    or it.lower().startswith('vr(') or it.lower().startswith('vi(') \
                    or it.lower().startswith('ir(') or it.lower().startswith('ii('):
                    all_keys_expend.append(it)
                elif it !='':
                #remove it to be compitible with hspicev2013 with no `v(1' but `1'
                #elif it.lower().startswith('v(') or it.lower().startswith('i('):
                    all_keys_expend.extend([it+'_real', it+'_imag'])
            all_keys = all_keys_expend
            Nblock=len(all_keys)
        print 'all_keys=', all_keys, Nblock, NallData, VSweep
        assert NallData%Nblock == 0, 'Data grouping Error!'
        if DEBUG:  # print out the data grouping info
            for idx in range(0, NallData, Nblock):
                print 'data range: %d -- %d'%(idx, idx+Nblock)
                print ' '.join(all_data_array[idx:idx+Nblock]) + '\n'
        return all_data_array, all_keys

    @classmethod
    def spdiff_ascii(cls, sw0t, sw0g, NsweepPoints=-1):
        """
        diff two sw? file (post=2|ASCII), only compare the values, return True if all the same
        also False
        """
        global NdataOneSweep
        matched = True
        if NsweepPoints != -1:
            NdataOneSweep = NsweepPoints
        wv_vectt,wv_keyst = cls.get_array_fix_width(sw0t)
        wv_vectg,wv_keysg = cls.get_array_fix_width(sw0g)
        if wv_keyst != wv_keysg:
            print wv_keyst
            print wv_keysg
        for idx, it in enumerate(wv_vectg):
            if it != wv_vectt[idx]:
                cir_basename = os.path.basename(sw0t)
                err = (float(it) - float(wv_vectt[idx]))/(float(it)+DELTA)
                if abs(err) > relTolSP and (abs(float(it)) > minVal and abs(float(wv_vectt[idx])) > minVal):
                    signal = wv_keysg[idx%len(wv_keysg)]
                    vx = wv_vectt[idx-(idx%len(wv_keysg))]
                    vp = '0'
                    if VSweep:
                        vp = VSweep[idx/(len(wv_keysg)*NdataOneSweep)]
                    print '**Err: %s signal=%s vx=%s vp=%s mismatch: %s  %s  diff: %f(%%) idx=%d'%(cir_basename, signal, vx, vp, wv_vectt[idx], it, err*100, idx)
                    matched = False
        return matched

    @classmethod
    def getSpiceMeas(cls,measfile):
        '''read hspice *.ms? *.mt? etc meas file into a dict with single index(=1)
        Example:
        $DATA1 SOURCE='HSPICE' VERSION='H-2013.03 64-BIT'
        .TITLE '##for calc idlin/idsat/vtlin/vtsat to simu kop trend'
         idsat_10_10      idsat_5_10       idsat_3_10       idsat_0.3_10
         alter#
          5.247e-06        2.469e-06        1.418e-06        1.833e-07
         1
         '''
        if isinstance(measfile,file):
            ifp=measfile
        else:
            ifp=open(measfile,'rU')
        par={}
        par_name=[]
        par_val=[]
        meas_start=0
        for line in ifp:
            line=line.strip(os.linesep).strip()

            line=line.split()
            if not len(line):
                continue #next please...
            #pass through comments, they stay asis
            if line[0][0]=='*' or line[0][0]=='$':
                continue #next please...
            elif '.TITLE' in line[0]: #beginning of meas data block of hspice meas file
                meas_start+=1
            elif meas_start>0:
                if line[0][0].isalpha() and ('failed' not in line):
                    for it in line:
                        if len(it)>0:
                            par_name.append(it)
                else:
                    for it in line:
                        if len(it)>0:
                            if it == 'failed':
                                it='%e'%MAXERROR
                            par_val.append(it)
        assert(len(par_name)==len(par_val))  #par, val should be matched
        for idx in range(len(par_name)):
            if not par_name[idx] in ['temper','alter#']:
                par[par_name[idx]]=float(par_val[idx])
        return par  ##a dict

    @classmethod
    def getSpiceMultMeas(cls,fname,keys=('w','l'),targkeys=None):
        '''read hspice *.ms? *.mt? etc meas file into a dict with multiple index and instance=keys
        if targkeys is given, only output the list in targkeys which is case insensitive
        Example:
    $DATA1 SOURCE='HSPICE' VERSION='H-2013.03 64-BIT'
    .TITLE '##for calc idlin/idsat/vtlin/vtsat to simu kop trend'
     index            width            length           idlin
                      idsat            id1              id2
                      ideff            idoff            iboff
                      vtlin            vtlin_vb         vtsat
                      vtsat_vb         idswing1         idswing2
                      idswingsat1      idswingsat2      swing_lin
                      swing_sat        body_lin         body_sat
                      dibl             max_gm           vth_gmi
                      id_max           vth_maxgm_macro  vth_maxgm_compact
                      gdsm             gmm              temper
                      alter#
     1                 9.000e-06        1.800e-07       -5.120e-05
                      -4.679e-04       -3.802e-04       -1.551e-05
                      -1.979e-04       -1.621e-11        1.050e-12
                      '''
        #print 'reading meas:',fname
        if isinstance(fname,file):
            ifp=fname
        else:
            ifp=open(fname,'rU')
        if not ifp:
            print 'Cannot Open *.ms0 file!'
            sys.exit()
        #specialkey=['index','width','length','temper','alter#']
        specialkey=['index','temper','alter#'] + list(keys)
        par_name=[]
        par_val=[]
        keyv=[]
        meas_start=0
        in_header=1
        header_end = False
        idx_meas=0
        idx_key=0
        diw_meas={}
        idxspkey={}
        idxtrgkey={}
        targkeys=[_key.upper() for _key in targkeys]
        for line in ifp:
            line=line.strip(os.linesep).rstrip()

            line=line.lower().split()
            #print 'line',line
            if not len(line):
                continue #next please...
            #pass through comments, they stay asis
            if line[0][0]=='*' or line[0][0]=='$':
                continue #next please...
            elif '.title' in line[0]: #beginning of meas data block of hspice meas file
                meas_start+=1
            elif meas_start>0:
                if not header_end and not line[0].isdigit() and ('failed' not in line):
                    in_header += 1
                    for it in line:
                        if len(it)>0:
                            #print 'it=',it
                            if it in specialkey:
                                idxspkey[it]=idx_key
                            elif targkeys and it.upper() in targkeys:
                                #print 'add into targkey',it
                                idxtrgkey[it]=idx_key
                                par_name.append(it)
                            else:
                                if targkeys is None:
                                    idxtrgkey[it] = idx_key
                                par_name.append(it)
                            idx_key+=1
                elif in_header >0:
                    #print '-'*5 + 'end of header' + '-'*5
                    header_end = True
                    val0=eval(line[0]) if not line[0].startswith('failed') else None
                    #print 'idx_meas=',idx_meas
                    if len(line) == 1 and isinstance(val0,int): #skip the one 'alter#' line with no data before it
                        continue
                    if isinstance(val0,int) and idx_meas != val0 : # tips: ms0 index should starts from 1 (sweep monte=1)
                        ##new index's value starts
                        if idx_meas>0:
                            keyv = [FMT6G%eval(s) for s in keyv]
                            diw_meas[tuple(keyv)]=[float(_x) if isinstance(_x, (str, unicode)) else _x for _x in par_val]
                            #diw_meas[tuple(keyv)]=np.array(par_val).astype(np.float)
                            idx_key=0
                        else:
                            idx_key=0
                        keyv=[]
                        par_val=[]
                        idx_meas+=1
                    for it in line:
                        if len(it)>0:
                            #print 'it=',it
                            if it == 'failed':
                                it=None
                                #it='%e'%MAXERROR
                            for kk in keys:
                                if idx_key == idxspkey[kk]:
                                    keyv.append(it)
                                    break
                            if idx_key not in idxspkey.values():
                                if idx_key in idxtrgkey.values():
                                    par_val.append(it)
                            idx_key+=1
        keyv = [FMT6G%eval(s) for s in keyv]
        diw_meas[tuple(keyv)]=[float(_x) if isinstance(_x, (str, unicode)) else _x for _x in par_val]
        #diw_meas[tuple(keyv)]=np.array(par_val).astype(np.float)
        # arrange par_val with targkeys' order
        targkeys_in_measfile = [s.upper() for s in par_name if s.upper() in targkeys]
        for key in diw_meas:
            val_tmp = diw_meas[key]
            diw_meas[key] = [val_tmp[targkeys_in_measfile.index(s)] for s in targkeys]
        #assert(len(par_name)==len(par_val))  #par, val should be matched
        return par_name,diw_meas

    @classmethod
    def getSpiceDC(cls,wvfile,**kwargs):
        '''read back hspice DC sw0 file results'''
        pass

    @classmethod
    def getSpiceACLis(cls,lisfile,**kwargs):
        '''read back hspice AC infor from *.lis file results'''
        cmd="grep '^y' %s -B 1 | sed '/^--/d'|sed '/^y/d' | awk '{print $2}'" % lisfile
        #grep "^x" -A 16 hspiceCkt_flip_nrd.lis | sed '/^--/d' | sed '/^x/d' | sed '/^y/d' | sed '/current/d' | sed '/v_d/d' | sed '/^$/d' >nrd.ref
        ret=os.popen(cmd).read().strip()
        if ret:
            ret=ret.split()
            ret = [float(_x) for _x in ret]
        else:
            ret=None
        return ret

def test_getMultTargMs0(fms0, targs=None):
    if isinstance(targs, str):
        targs=targs.split(',')
    par_name, diw=SpiceReader.getSpiceMultMeas(fms0, targkeys=targs)
    print par_name
    print diw

def getWvBinCmpDC(topdir, spbasename1, spbasename2, N):
    for idx in range(int(N)):
        wvbinf1=os.path.join(topdir,spbasename1+'.sw%d'%(idx))
        wvbinf2=os.path.join(topdir,spbasename2+'.sw%d'%(idx))
        SpiceReader.spdiff_ascii(wvbinf1, wvbinf2)

def compTwoWv(wvbinf1,wvbinf2,NsweepPoints=-1):
    if os.path.exists(wvbinf1) and os.path.exists(wvbinf2):
        matched=SpiceReader.spdiff_ascii(wvbinf1, wvbinf2,int(NsweepPoints))
    else:
        print '**Err: %s or %s not found!' % (wvbinf1,wvbinf2)
        matched = False

    if matched:
        print '#'*20 + 'Test %s vs. %s PASS' % (wvbinf1, wvbinf2) + '#'*20
    else:
        print '#'*20 + 'Test %s vs. %s FAIL' % (wvbinf1, wvbinf2) + '#'*20

def compTwoMeas(measf1,measf2):
    if os.path.exists(measf1) and os.path.exists(measf2):
        diw1=SpiceReader.getSpiceMeas(measf1)
        diw2=SpiceReader.getSpiceMeas(measf2)
        matched = True
        assert len(diw1) == len(diw2)
        for key, val in diw1.iteritems():
            valg=diw2[key]
            if val != valg:
                err = (float(val) - float(valg))/(float(val)+DELTA)
                if abs(err) > relTolSP_meas and (abs(float(val)) > minVal and abs(float(valg)) > minVal):
                    signal = key
                    cir_basename = os.path.basename(measf1)
                    print '**Err: %s signal=%s mismatch: %s  %s  diff: %f(%%)'%(cir_basename, signal, val, valg, err*100)
                    matched = False
    else:
        print '**Err: %s or %s not found!' % (measf1, measf2)
        matched = False

    if matched:
        print '#'*20 + 'Test %s vs. %s PASS' % (measf1, measf2) + '#'*20
    else:
        print '#'*20 + 'Test %s vs. %s FAIL' % (measf1, measf2) + '#'*20

def main(simu_mode, f1, f2, NsweepPoints=-1):
    if simu_mode == 'tran':
        compTwoMeas(f1, f2)
    elif simu_mode == 'ac' or simu_mode == 'dc':
        compTwoWv(f1, f2, NsweepPoints)
    else:
        print 'Warn: No this type simulaion %s'%simu_mode
    return 0

if __name__=='__main__':
    if len(sys.argv) >2:
        #getWvBinCmpDC(*sys.argv[1:])
        #compTwoWv(*sys.argv[1:])
        #compTwoMeas(*sys.argv[1:])
        sys.exit(main(*sys.argv[1:]))


