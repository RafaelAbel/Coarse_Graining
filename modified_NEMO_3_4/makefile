# DCM NEMODRAK_3.3  (support for FCM)
# ------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MMS
DCMCCFG=dcmc.cfg
include $(DCMCCFG)
#:MMS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# --- Variables to be customized by users ---

# CONF = generic name of the configuration
# CASE   = case id for this variant
# CASEREF = reference case for this configuration (a former case ID or 'none' )


CONF = ORCA025.L46
CASE = KRA020
CONFREF = $(CONF)
#CASEREF = 'KREF3.4.1'
CASEREF = 'KIEL3.4.1'
CASEREF_PREFIX = $(HOMEDCM)/CONFIGS/$(CONFREF)-$(CASEREF)

# PREV_CONFIG  (type the full path or 'none' )
PREV_CONFIG = '/gpfs/fs6/home-geomar/smomw184/CONFIG_ORCA025.L46/ORCA025.L46-KRA010/'

# REFONLY ( yes or no : if yes, only reference is used, not the drakkar customs)
REFONLY = 'no'

# MACHINE : name of architecture for compilation. Assume that a arch_${MACHINE}.fcm exists in ARCH 
# to know available arch_xxx.fcm, use dcm_lsarch.ksh. You can edit/add your own arch_xxx.fcm in
# CONF-CASE/ARCH/ directory
#MACHINE = crayxc30
MACHINE = sxf90_sxace_nesh

# NCOMPIL_PROC : number of procs to use for the compilation of the code.
NCOMPIL_PROC = 6

# set to 'use' the NEMO component required for this config, to anything else
# in the actual state of NEMO, OPA and LIM must be set to use for correct compilation
#           if notused
OPA = 'use'
LIM2 = 'use'
LIM3 = 'notused'
TOP = 'notused'
OFF = 'notused'
TAM = 'notused'
AGRIF = 'notused'
# if set to 'check' maintain an 'install_history' file giving 
#  the date of install and the svn info of $HOMEDCM
SVN = 'notused'

# --- All customizations by users is done above this line ... Only experts touch things below ... ---
#----------------------------------------------------------------------------------------------------

# REF is the directory name holding the reference NEMOGCM tree from a NEMO-TEAM download
# CUS is the directory name holding the permanent customized files not yet under the reference, but valid
#     for all the users. Will overwrite the REF during install

REF = $(REFDIR)
CUS = $(CUSDIR)

CONFIG_DIR = $(UDIR)/CONFIG_$(CONF)/$(CONF)-$(CASE)
COMPIL_DIR = $(CDIR)/W$(CONF)-$(CASE)
WK =  $(COMPIL_DIR)/CONFIG/$(CONF)-$(CASE)/WORK

CFILE := $(CONF)-$(CASE)

# CFILE is used below for setting the cfg.txt file in TOOLS/COMPILE
ifeq ($(OPA),'use')
  CFILE:=$(CFILE) OPA_SRC
endif

ifeq ($(LIM2),'use')
  CFILE:=$(CFILE) LIM_SRC_2
endif

ifeq ($(LIM3),'use')
  CFILE:=$(CFILE) LIM_SRC_3
endif

ifeq ($(TOP),'use')
  CFILE:=$(CFILE) TOP_SRC
endif

ifeq ($(TAM),'use')
  CFILE:=$(CFILE) OPATAM_SRC
endif

ifeq ($(OFF),'use')
  CFILE:=$(CFILE) OFF_SRC
endif

ifeq ($(AGRIF),'use')
  CFILE:=$(CFILE) NST_SRC
endif


#####
# define the target that must always be done, in any condition, never up to date !
.PHONY: restore 

logcompile: compile modulelist

compile: copyfile copynamelist
	@echo ' '
	@echo '... compile opa  with fcm ...'
	@cp CPP.keys $(PDIR)/RUN_$(CONF)/$(CONF)-$(CASE)/EXE/
	@cpp2fcm.ksh > $(COMPIL_DIR)/CONFIG/$(CONF)-$(CASE)/cpp_$(CONF)-$(CASE).fcm
	@cd $(COMPIL_DIR)/CONFIG ; ./makenemo -n $(CONF)-$(CASE)  -j $(NCOMPIL_PROC)
	@dcm_chkmsh.ksh $(COMPIL_DIR)/CONFIG/$(CONF)-$(CASE)/BLD/ppsrc/nemo/par_oce.f90 >> $(PDIR)/RUN_$(CONF)/$(CONF)-$(CASE)/EXE/CPP.keys
	@ln -sf $(COMPIL_DIR)/CONFIG/$(CONF)-$(CASE)/BLD/bin/nemo.exe $(PDIR)/RUN_$(CONF)/$(CONF)-$(CASE)/EXE/opa


copyfile:
ifeq ($(CASEREF),'none')
	@echo "no CASEREF to take into account"
else
	@echo " CASEREF=$(CASEREF_PREFIX)"
  ifeq ($(OPA),'use') 
	# MMS:
	#@echo " retrieving CPP.keys from $(CASEREF_PREFIX)"
	#-cp -f -p -R $(CASEREF_PREFIX)/CPP.keys ./CPP.keys
	# :MMS
	@-cp -f -p -R $(CASEREF_PREFIX)/OPA_SRC/  $(COMPIL_DIR)/NEMO
	@-cp -f -p -R $(CASEREF_PREFIX)/IOIPSL/*  $(COMPIL_DIR)/EXTERNAL/IOIPSL/src/.
  endif
  ifeq ($(TAM),'use') 
	@-cp -f -p -R $(CASEREF_PREFIX)/OPATAM_SRC/  $(COMPIL_DIR)/NEMO
  endif
  ifeq ($(LIM2),'use')
	@-cp -f -p -R $(CASEREF_PREFIX)/LIM_SRC_2/  $(COMPIL_DIR)/NEMO
  endif
  ifeq ($(LIM3),'use')
	@-cp -f -p -R $(CASEREF_PREFIX)/LIM_SRC_3/  $(COMPIL_DIR)/NEMO
  endif
  ifeq ($(TOP),'use')
	@-cp -f -p -R $(CASEREF_PREFIX)/TOP_SRC/  $(COMPIL_DIR)/NEMO
  endif
  ifeq ($(OFF),'use')
	@-cp -f -p -R $(CASEREF_PREFIX)/OFF_SRC/  $(COMPIL_DIR)/NEMO
  endif
  ifeq ($(AGRIF),'use')
	@-cp -f -p -R $(CASEREF_PREFIX)/NST_SRC/  $(COMPIL_DIR)/NEMO
  endif
endif
ifeq ($(OPA),'use') 
	@-cp -f -p -R  $(CONFIG_DIR)/OPA_SRC/ $(COMPIL_DIR)/NEMO
	@-cp -f -p -R  $(CONFIG_DIR)/IOIPSL/*  $(COMPIL_DIR)/EXTERNAL/IOIPSL/src/.
endif
ifeq ($(TAM),'use') 
	@-cp -f -p -R $(CONFIG_DIR)/OPATAM_SRC/  $(COMPIL_DIR)/NEMO
endif
ifeq ($(LIM2),'use')
	@-cp -f -p -R  $(CONFIG_DIR)/LIM_SRC_2/ $(COMPIL_DIR)/NEMO
endif
ifeq ($(LIM3),'use')
	@-cp -f -p -R  $(CONFIG_DIR)/LIM_SRC_3/ $(COMPIL_DIR)/NEMO
endif
ifeq ($(TOP),'use')
	@-cp -f -p -R  $(CONFIG_DIR)/TOP_SRC/ $(COMPIL_DIR)/NEMO
endif
ifeq ($(OFF),'use')
	@-cp -f -p -R  $(CONFIG_DIR)/OFF_SRC/ $(COMPIL_DIR)/NEMO
endif
ifeq ($(AGRIF),'use')
	@-cp -f -p -R  $(CONFIG_DIR)/NST_SRC/ $(COMPIL_DIR)/NEMO
endif
	@-cp -f -p -R $(CONFIG_DIR)/ARCH $(COMPIL_DIR)
	@gmake links

copynamelist:
#	@rcp $(CONFIG_DIR)/CTL/namelist_ice.NATL4 uqbar:/u/rech/cli/rcli099/$(CONF)-$(CASE)/namelist_ice.NATL4
#	@rcp $(CONFIG_DIR)/CTL/namelist.NATL4 uqbar:/u/rech/cli/rcli099/$(CONF)-$(CASE)/namelist.NATL4


##############

install:   copyref copycpp update_history  WORK cleansvn copyctl

reinstall: copyref  

# Target copyconfig to duplicate an existant configuration to local directory
copyconfig:
ifeq ($(PREV_CONFIG),'none')
	@echo "no previous configuration defined"
else
	@echo "   copy a previous configuration $(PREV_CONFIG) to local directory"
	-cp -R $(PREV_CONFIG)/OPA_SRC .
	-cp -R $(PREV_CONFIG)/LIM_SRC_2 . 
	-cp -R $(PREV_CONFIG)/LIM_SRC_3 . 
	-cp -R $(PREV_CONFIG)/TOP_SRC . 
	-cp -R $(PREV_CONFIG)/NST_SRC . 
	-cp -R $(PREV_CONFIG)/IOIPSL . 
	-cp -R $(PREV_CONFIG)/OPATAM_SRC . 
endif
	@gmake links


# Target copyconfigall do copyconfig + copy CPP.keys and template namelists
copyconfigall:  copyconfig 
ifeq ($(PREV_CONFIG),'none')
	@echo "no previous configuration defined"
else
	@echo "   copy a previous configuration $(PREV_CONFIG) to local directory"
	-cp  $(PREV_CONFIG)/ARCH/* ARCH/
	-cp  $(PREV_CONFIG)/CPP.keys .
	-cp  $(PREV_CONFIG)/makefile makefile.prev
	-cp  $(PREV_CONFIG)/nam* $(PDIR)/RUN_$(CONF)/$(CONF)-$(CASE)/CTL/
endif
	@gmake links


# Target copyref is to be done each time we need to rebuild the
# source code from the reference sources.
copyref:
	@echo "    copy reference version to $(COMPIL_DIR) "
	@cp -f -R $(REF)/NEMOGCM/NEMO $(COMPIL_DIR)/.
	@cp -f -R $(REF)/NEMOGCM/EXTERNAL $(COMPIL_DIR)/.
	@if [ ! -d $(COMPIL_DIR)/CONFIG ] ; then mkdir $(COMPIL_DIR)/CONFIG ; fi
	@if [ ! -d $(COMPIL_DIR)/CONFIG/$(CONF)-$(CASE) ] ; then mkdir $(COMPIL_DIR)/CONFIG/$(CONF)-$(CASE) ; fi
	@if [ ! -d $(COMPIL_DIR)/CONFIG/$(CONF)-$(CASE)/EXP00 ] ; then mkdir $(COMPIL_DIR)/CONFIG/$(CONF)-$(CASE)/EXP00 ; fi
	@cp -f -R $(REF)/NEMOGCM/CONFIG/makenemo $(COMPIL_DIR)/CONFIG/.
	@cp -f -R $(REF)/NEMOGCM/ARCH $(COMPIL_DIR)/.
	@if [ ! -d $(COMPIL_DIR)/TOOLS ] ; then mkdir $(COMPIL_DIR)/TOOLS ; fi
	@cp -f -R $(REF)/NEMOGCM/TOOLS/COMPILE $(COMPIL_DIR)/TOOLS/.
ifeq ($(REFONLY),'yes')
	@echo "    WARNING : DRAKKAR version is not used "
else
	@echo "    copy DRAKKAR version to $(COMPIL_DIR) "
	@cp -f -R $(CUS)/NEMOGCM/NEMO $(COMPIL_DIR)/.
	@cp -f -R $(CUS)/NEMOGCM/EXTERNAL $(COMPIL_DIR)/.
	@cp -f -R $(CUS)/NEMOGCM/CONFIG/makenemo $(COMPIL_DIR)/CONFIG/.
	@cp -f -R $(CUS)/NEMOGCM/ARCH $(COMPIL_DIR)/.
	@cp -f -R $(CUS)/NEMOGCM/TOOLS/COMPILE $(COMPIL_DIR)/TOOLS/.
endif
	@gmake copyfile
	@if [ ! -d $(COMPIL_DIR)/CONFIG/$(CONF)-$(CASE) ] ; \
	then mkdir $(COMPIL_DIR)/CONFIG/$(CONF)-$(CASE) ; fi
	@echo $(CFILE) > $(COMPIL_DIR)/CONFIG/cfg.txt
	@cpp2fcm.ksh > $(COMPIL_DIR)/CONFIG/$(CONF)-$(CASE)/cpp_$(CONF)-$(CASE).fcm
	@cd $(COMPIL_DIR)/CONFIG/ ; ./makenemo -n $(CONF)-$(CASE) -m $(MACHINE) -j0


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MMS:
$(DCMCCFG):
	@echo "No $(DCMCCFG) avaiable. Using shell environment variables for Drakkar Config Manager..."
	@echo "DCMC=$(DCMC)">tmp.cfg
	@echo "HOMEDCM=$(HOMEDCM)">>tmp.cfg
	@echo "REFDIR=$(REFDIR)">>tmp.cfg
	@echo "CUSDIR=$(CUSDIR)">>tmp.cfg
	@echo "UDIR=$(UDIR)">>tmp.cfg
	@echo "CDIR=$(CDIR)">>tmp.cfg
	@echo "SDIR=$(SDIR)">>tmp.cfg
	@echo "PDIR=$(PDIR)">>tmp.cfg
	@echo "TDIR=$(TDIR)">>tmp.cfg
	@-mv -f tmp.cfg $(DCMCCFG)

createcfg: $(DCMCCFG)

copyctl:
	@if [ ! -d $(PDIR)/RUN_$(CONF)/$(CONF)-$(CASE)/CTL/EXAMPLE ] ; then        \
		echo create EXAMPLE directory in $(PDIR)/RUN_$(CONF)/$(CONF)-$(CASE)/CTL ;\
		mkdir $(PDIR)/RUN_$(CONF)/$(CONF)-$(CASE)/CTL/EXAMPLE                    ;\
	fi
	@echo Copying Example Job files into $(PDIR)/RUN_$(CONF)/$(CONF)-$(CASE)/CTL/EXAMPLE .
	@sed -e "s/:CONFIG:/$(CONF)/;s/:CASE:/$(CASE)/" $(CASEREF_PREFIX)/includefile.ksh > $(PDIR)/RUN_$(CONF)/$(CONF)-$(CASE)/CTL/EXAMPLE/includefile.ksh
	@sed -e "s/:USEREMAIL:/$(USEREMAIL)/;s/:CASE:/$(CASE)/" $(CASEREF_PREFIX)/CONF-CASE.ksh > $(PDIR)/RUN_$(CONF)/$(CONF)-$(CASE)/CTL/EXAMPLE/$(CONF)-$(CASE).ksh
	@-cp $(CASEREF_PREFIX)/CONF-CASE.db $(PDIR)/RUN_$(CONF)/$(CONF)-$(CASE)/CTL/EXAMPLE/$(CONF)-$(CASE).db
	@sed -e "s/:CONFIG:/$(CONF)/;s/:CASE:/$(CASE)/" $(CASEREF_PREFIX)/namelist.CONF-CASE > $(PDIR)/RUN_$(CONF)/$(CONF)-$(CASE)/CTL/EXAMPLE/namelist.$(CONF)-$(CASE)
	@-if [ -f $(CASEREF_PREFIX)/namelist_ice_lim2 ] ; then cp $(CASEREF_PREFIX)/namelist_ice_lim2 $(PDIR)/RUN_$(CONF)/$(CONF)-$(CASE)/CTL/EXAMPLE/namelist_ice_lim2 ; fi
	@-if [ -f $(CASEREF_PREFIX)/namelist_ice_lim3 ] ; then cp $(CASEREF_PREFIX)/namelist_ice_lim3 $(PDIR)/RUN_$(CONF)/$(CONF)-$(CASE)/CTL/EXAMPLE/namelist_ice_lim3 ; fi
	@-if [ -f $(CASEREF_PREFIX)/convert.ksh ] ; then sed -e "s/:USEREMAIL:/$(USEREMAIL)/" $(CASEREF_PREFIX)/convert.ksh > $(PDIR)/RUN_$(CONF)/$(CONF)-$(CASE)/CTL/EXAMPLE/convert.ksh ; fi
	@-cp $(CASEREF_PREFIX)/iodef.xml $(PDIR)/RUN_$(CONF)/$(CONF)-$(CASE)/CTL/EXAMPLE/iodef.xml
	@-cp $(CASEREF_PREFIX)/xmlio_server.def $(PDIR)/RUN_$(CONF)/$(CONF)-$(CASE)/CTL/EXAMPLE/xmlio_server.def
	@#-cp $(CASEREF_PREFIX)/mpi.conf $(PDIR)/RUN_$(CONF)/$(CONF)-$(CASE)/CTL/EXAMPLE/mpi.conf
ifeq ($(AGRIF),'use')
	@echo Copying AGRIF files, too....
	@-cp $(CASEREF_PREFIX)/1_CONF-CASE.db $(PDIR)/RUN_$(CONF)/$(CONF)-$(CASE)/CTL/EXAMPLE/1_$(CONF)-$(CASE).db
	@-sed -e "s/:CONFIG:/$(CONF)/;s/:CASE:/$(CASE)/" $(CASEREF_PREFIX)/1_namelist.CONF-CASE > $(PDIR)/RUN_$(CONF)/$(CONF)-$(CASE)/CTL/EXAMPLE/1_namelist.$(CONF)-$(CASE)
	@-if [ -f $(CASEREF_PREFIX)/1_namelist_ice_lim2 ] ; then cp $(CASEREF_PREFIX)/1_namelist_ice_lim2 $(PDIR)/RUN_$(CONF)/$(CONF)-$(CASE)/CTL/EXAMPLE/1_namelist_ice_lim2 ; fi
	@-if [ -f $(CASEREF_PREFIX)/1_namelist_ice_lim3 ] ; then cp $(CASEREF_PREFIX)/1_namelist_ice_lim3 $(PDIR)/RUN_$(CONF)/$(CONF)-$(CASE)/CTL/EXAMPLE/1_namelist_ice_lim3 ; fi
	@-cp $(CASEREF_PREFIX)/AGRIF_FixedGrids.in $(PDIR)/RUN_$(CONF)/$(CONF)-$(CASE)/CTL/EXAMPLE/AGRIF_FixedGrids.in
endif


svninit:
	@svn import --username $(SVNUSER) https://$(SVNSERVER)/EXP/$(CONF)-$(CASE)/Code
	@svn co  --username $(SVNUSER) https://$(SVNSERVER)/EXP/$(CONF)-$(CASE)/Code

modulelist:
	@printf '~%.0s' {1..100} >> log.make
	@printf '\n%s ' 'Compilation at' >> log.make
	@date >> log.make
	@printf '\n' >> log.make
	@bash -c 'module() { eval `/usr/bin/modulecmd bash $$*` ; }; module list -l >> log.make 2>&1; '

# :MMS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




WORK:
	@ln -s $(WK) WORK

cleaninst:
	@-\rm -rf $(CDIR)/W$(CONF)-$(CASE)/*

clean:
	#cd $(WK) ; gmake clean
	@cd $(COMPIL_DIR)/CONFIG ; ./makenemo clean

# restore target is used to restore a SINGLE file from the 3levels references
# must be invoked like : make restore RESTORE=OPA_SRC/SBC/ocesbc.F90  (for instance, with the full path relative to NEMO)
# it is used automatically with rmfile.

restore:
ifeq ($(RESTORE),)
	@echo "    Nothing to restore, guy !"
else
	@echo "    Restoring  $(RESTORE) from REFERENCES "
	@-cp  $(REF)/NEMOGCM/NEMO/$(RESTORE)  $(COMPIL_DIR)/NEMO/$(RESTORE)
  ifeq ($(REFONLY),'yes') 
	@echo "    WARNING : nothing restored from DRAKKAR "
  else
	@echo "    Restoring  $(RESTORE) from DRAKKAR "
	@-cp  $(CUS)/NEMOGCM/NEMO/$(RESTORE)  $(COMPIL_DIR)/NEMO/$(RESTORE)
  endif
  ifeq ($(CASEREF),'none')
	@echo "no CASEREF to take into account"
  else
	@echo " CASEREF=$(CASEREF) "
	@-cp  $(CASEREF_PREFIX)/$(RESTORE)  $(COMPIL_DIR)/NEMO/$(RESTORE)
  endif
endif


restoreioipsl:
ifeq ($(RESTORE),)
	@echo "    Nothing to restore, guy !"
else
	@echo "    Restoring  $(RESTORE) from REFERENCES "
	@-cp  $(REF)/NEMOGCM/EXTERNAL/IOIPSL/src/$(RESTORE)  $(COMPIL_DIR)/EXTERNAL/IOIPSL/src/$(RESTORE)
  ifeq ($(REFONLY),'yes') 
	@echo "    WARNING : nothing restored from DRAKKAR "
  else
	@echo "    Restoring  $(RESTORE) from DRAKKAR "
	@-cp  $(CUS)/NEMOGCM/EXTERNAL/IOIPSL/src/$(RESTORE)  $(COMPIL_DIR)/EXTERNAL/IOIPSL/src/$(RESTORE)
  endif
  ifeq ($(CASEREF),'none')
	@echo "no CASEREF to take into account"
  else
	@echo " CASEREF=$(CASEREF) "
	@-cp  $(CASEREF_PREFIX)/IOIPSL/$(RESTORE)  $(COMPIL_DIR)/EXTERNAL/IOIPSL/src/$(RESTORE)
  endif
endif
links:
	@-ln -sf OPA_SRC/*.*0 .
	@-ln -sf OPA_SRC/*/*.*0 .
	@-ln -sf OPATAM_SRC/*.*0 .
	@-ln -sf OPATAM_SRC/*/*.*0 .
	@-ln -sf LIM_SRC_2/*.*0 .
	@-ln -sf LIM_SRC_3/*.*0 .
	@-ln -sf TOP_SRC/*.*0 .
	@-ln -sf TOP_SRC/*/*.*0 .
	@-ln -sf NST_SRC/*.*0 .
	@-ln -sf IOIPSL/*.*0 .
	@-ln -sf ARCH/*.fcm .

update_history:
ifeq ($(SVN),'check')
	@date >> install_history
	@svn info $(HOMEDCM) >> install_history
	@echo install_history updated
else
	@echo SVN = $(SVN)
	@echo no tracks of install history
endif

cleansvn:
	@-find . -name ".svn" -exec \rm -rf {} \;
   
copycpp: 
	-cp -f $(CASEREF_PREFIX)/CPP.keys ./CPP.keys

