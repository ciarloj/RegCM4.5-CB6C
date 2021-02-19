#!/usr/bin/python

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#    This file is part of ICTP RegCM.
#
#    ICTP RegCM is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    ICTP RegCM is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
#
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Script for running RegCM regression tests by M. Scarcia
#

# import system modules
import os,sys,shutil,subprocess,time

# import own modules
import parsing_editing,nc_stuff

def main(argv):

    if (len(sys.argv) < 2):
        print "Please specify a configuration file!"
        sys.exit(1)
        
    cfg=sys.argv[1]

    # get all the options from cfg file
    options = parsing_editing.parse_config(cfg)

    datadir = options["DATADIR"]
    bindir = options["BINDIR"]
    testdir = options["TESTDIR"]
    namelistdir = options["NLDIR"]
    referencedir = options["REFDIR"]
    teststodo = options["TESTSTODO"]
    mpistring=options["MPISTRING"]
    run_serial=int(options["SERIAL"])
    run_preproc=int(options["PREPROC"])
    run_clm=int(options["USECLM"])
    run_band=int(options["USEBAND"])
    run_diff=int(options["DIFF"])
    simdays=int(options["SIMDAYS"])

    datadir = os.path.abspath(datadir)
    testdir = os.path.abspath(testdir)
    bindir = os.path.abspath(bindir)
    referencedir = os.path.abspath(referencedir)

    # check running options compatibility + choose binary

    main_bin = "regcmMPI"

    if (run_clm + run_band) > 1:
        print "Running with BAND and CLM enabled is currently not supported!"
        os.sys.exit(1)

    if run_serial == 1 :
        main_bin = "regcmSerial"
    elif run_clm == 1 :
        main_bin = "regcm_clM"
    elif run_band == 1 :
        main_bin = "regcm_band"

    # now check if bins exist
    
    #print "bin = ",main_bin
    
    # will put the diff variables here?
    
    # check what tests to do
    if teststodo.rfind(",") > -1 :
        tests=teststodo.split(",")
        listtype=0
    elif teststodo.rfind("-") > -1:
        tests=teststodo.split("-")
        imin=int(tests[0])
        imax=int(tests[1])
        listtype=1
    else :
        tests=int(teststodo)
        listtype=2
    
    TOT_TESTS = 12 # number of total tests present

    if not os.path.isdir(testdir):
        os.mkdir(testdir)

    if listtype == 2 :
        if tests == 0 :
            imin = 1
            imax = TOT_TESTS
        else :
            imin = int(tests)
            imax = int(tests)
    elif listtype == 1 :
        imin = int(tests[0])
        imax = int(tests[1])
    else :
        imin = 0
        imax = len(tests)-1

    #print "imin =",imin
    #print "imax =",imax

    # main loop over tests        
    for i in range(imin,imax+1):

        if listtype == 0 :
            testname="test_"+str(tests[i]).zfill(3)
        else :
            testname="test_"+str(i).zfill(3)

        # create simulation directory tree
        simdir=testdir+"/"+testname
        testrefdir=referencedir+"/"+testname

        if not os.path.isdir(simdir):
            os.mkdir(simdir)

        if not os.path.isdir(simdir+"/input"):
            os.mkdir(simdir+"/input")
        if not os.path.isdir(simdir+"/output"):
            os.mkdir(simdir+"/output")

        namelist = simdir+"/regcm.in"
        shutil.copy(namelistdir+"/"+testname+".in",namelist)

        if (run_clm == 1) :
            try :
                shutil.copy(datadir+"/CLM/pft-physiology.c070207",simdir+"/input") # hardcoded for now
            except IOError :
                print "File",datadir+"/CLM/pft-physiology.c070207","not found. Stopping execution."
                os.sys.exit(1)
            
        # find idate0 and edit namelist for desired sim length
        idate0 = parsing_editing.parse_dates(namelist,simdays)
        
        #edit the namelist here
        parsing_editing.edit_namelist(namelist,datadir,simdir)

        # open log file
        writelog=True
        try:
            log = open(testname+".log","w")
        except :
            print "Unable to write log!"
            writelog=False

        exit_status = 0 # won't run Main if PreProc crashes...
        
        # run preproc
        if (run_preproc == 1):

            # check if binaries actually exist
            if not (os.path.isfile(bindir+"/terrain")) :
                print "Terrain binary not found! Skipping further steps."
                os.sys.exit(exit_status)
        
            p_terrain = subprocess.Popen(bindir+"/terrain "+namelist,stdout=log,stderr=log,shell=True)
            while not p_terrain.poll():
               if p_terrain.returncode is not None: break
               print "   ...Terrain: still alive..."
               time.sleep(10)

            if p_terrain.wait() != 0:
                print "\nError: Terrain in",testname,"crashed!!\n"
                exit_status = 1
            else:
                print "Terrain in",testname,"passed."

            if not (os.path.isfile(bindir+"/sst")) :
                print "SST binary not found! Skipping further steps."
                os.sys.exit(exit_status)
    
            p_sst=subprocess.Popen(bindir+"/sst "+namelist,stdout=log,stderr=log,shell=True)
            while not p_sst.poll():
               if p_sst.returncode is not None: break
               print "   ...SST: still alive..."
               time.sleep(10)
            if p_sst.wait() != 0:
                print "\nError: SST in",testname,"crashed!!\n"
                exit_status = 1
            else :
                print "SST in",testname,"passed."

            if not (os.path.isfile(bindir+"/icbc")) :
                print "ICBC binary not found! Skipping further steps."
                os.sys.exit(exit_status)
            
            p_icbc=subprocess.Popen(bindir+"/icbc "+namelist,stdout=log,stderr=log,shell=True)
            while not p_icbc.poll():
               if p_icbc.returncode is not None: break
               print "   ...ICBC: still alive..."
               time.sleep(10)
            if p_icbc.wait() != 0:
                print "\nError: ICBC in",testname,"crashed!!\n"
                exit_status = 1
            else :
                print "ICBC in",testname,"passed."
            
            if run_clm == 1:
                if not (os.path.isfile(bindir+"/clm2rcm")) :
                    print "clm2rcm binary not found! Skipping further steps."
                    os.sys.exit(exit_status)
                    
                p_clmpre=subprocess.Popen(bindir+"/clm2rcm "+namelist,stdout=log,stderr=log,shell=True)
                while not p_clmpre.poll():
                    if p_clmpre.returncode is not None: break
                    print "   ...CLMPRE: still alive..."
                    time.sleep(10)
                if p_clmpre.wait() != 0:
                    print "\nError: clm2rcm in",testname,"crashed!!\n"
                    exit_status = 1
                else :
                    print "clm2rcm in",testname,"passed."

            # compare preproc output only if everything went ok
            # and diff selected
            if (exit_status == 0) and (run_diff == 1) :

                dom_diff={}
                icbc_diff={}

                domain_file = "/input/"+testname+"_DOMAIN000.nc"
                icbc_file = "/input/"+testname+"_ICBC."+idate0+".nc"

                domain_vars = ["topo","landuse"]
                icbc_vars = ["u","v","t","ts"]

                # domain
                for var in domain_vars :    
                    dom_diff[var] = nc_stuff.compare_nc_file(simdir+domain_file,testrefdir+domain_file,var).rstrip("\n")
                    print var+" =",dom_diff[var]

                # icbc
                for var in icbc_vars :
                    icbc_diff[var] = nc_stuff.compare_nc_file(simdir+icbc_file,testrefdir+icbc_file,var).rstrip("\n")
                    print var+" =",icbc_diff[var]
                    
            sys.stdout.flush()

        # if preproc is ok, run main
        if exit_status == 0 :

            if not os.path.isfile(bindir+"/"+main_bin) :
                print "Main RegCM binary not found! Skipping further steps."
                log.close()
                os.sys.exit(exit_status)
            
            p_regcm=subprocess.Popen(mpistring+" "+bindir+"/"+main_bin+" "+namelist,stdout=log,stderr=log,shell=True)
            while not p_regcm.poll():
               if p_regcm.returncode is not None: break
               print "   ...RegCM: still alive..."
               time.sleep(10)
            if p_regcm.wait() != 0:
                print "\nError: RegCM",testname,"crashed!!\n"
                exit_status = 1
            else :
                print "RegCM",testname,"passed."
        else :
            print "Preprocessing did not complete correctly, RegCM main skipped."

        log.close()

        if exit_status == 1:
            outlog = open(testname+".log","r")
            stdouterr = outlog.read()
            print stdouterr

        # if everything ok and diff enabled compare output      
        if (exit_status == 0) and (run_diff == 1):

            srf_diff={}
            srf_file="/output/"+testname+"_SRF."+idate0+".nc"

            srf_vars = ["t2m"]

            for var in srf_vars :
                srf_diff[var] = nc_stuff.compare_nc_file(simdir+srf_file,testrefdir+srf_file,var).rstrip("\n")
                print var+" =",srf_diff[var]

    sys.stdout.flush()

    # end of the big loop
    if exit_status == 1:
        print "Warning! Some tests failed!"

    print "\n****  Test script terminated.  ****"

    sys.exit(exit_status)
    
if __name__ == "__main__":
    main(sys.argv[1:])

