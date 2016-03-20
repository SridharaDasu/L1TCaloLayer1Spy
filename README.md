#L1TCaloLayer1Spy

This package contains an unusual CMSSW producer, which uses TCP/IP
client to connect to the server running on CTP7 embedded linux systems
and creates raw data objects for input and output of CMS Level-1
Calorimeter Layer-1 System, i.e.,
```bash
  produces<HcalTrigPrimDigiCollection>();
  produces<EcalTrigPrimDigiCollection>();
  produces<CaloTowerBxCollection>();
```

When run on machines with P5 online environment it can obtain a contiguous
set of zerobias crossings being processed by the Layer-1 system.  These
are particularly useful for Layer-1 system monitoring and pileup level
determination.

Before using this package one needs to obtain the TCP/IP client software
using the online cactus SVN repository and install it somewhere.  Note 
that this depends on google protobuf-lite, which needs to have been 
installed on the build machine.

```bash
svn co svn+ssh://svn.cern.ch/reps/cactus/trunk/cactusprojects/calol1/extern/UCT2016Layer1CTP7Client UCT2016Layer1CTP7Client
```

If you don't have access to online cactus svn, you can take stuff
from my area.  For your convenience the client is currently installed 
on CERN AFS in /afs/cern.ch/user/d/dasu/public and the directory trees 
are accessible to all system:anyuser.  So, you should be able to build 
whereever you have access to AFS with CERN tokens.

Typically this package is run along with the emulator so that we
can make Spy digis to emulated digis comparison to validate the
Layer-1 hardware.  Check that out from L1Trigger/L1TCaloLayer1
following Mike Mulhearn's recipe documented in cms-l1t-offline
repository.  Currently slightly modified version of it is:

```bash
cmsrel CMSSW_8_0_2
cd CMSSW_8_0_2
cmsenv
cd $CMSSW_RELEASE/extern
svn co svn+ssh://svn.cern.ch/reps/cactus/trunk/cactusprojects/calol1/extern/UCT2016Layer1CTP7Client UCT2016Layer1CTP7Client
cd UCT2016Layer1CTP7Client
make install
cd $CMSSW_RELEASE/src
git cms-init
git remote add cms-l1t-offline git@github.com:cms-l1t-offline/cmssw.git
git fetch cms-l1t-offline
git cms-merge-topic cms-l1t-offline:dasu-updates-$CMSSW_VERSION
git cms-addpkg L1Trigger/L1TCommon
cd $CMSSW_RELEASE/src/L1Trigger
git clone git@github.com:SridharaDasu/L1TCaloLayer1Spy.git
cd $CMSSW_RELEASE/src
scram b -j 8
```

To debug the Spy unpacker (L1TCaloLayer1Spy/testL1TCaloLayer1Spy.py) and emulator
(L1TCaloLayer1/testL1TCaloLayer1EmulatorWithSpy.py), there are several options
for getting input.  You can use pre-captured files, pre-determined
patterns or ethernet captures.  You can run using the syntax:

```bash
cmsRun testL1TCaloLayer1EmulatorWithSpy.py setupString=patterns:DecreasingEM,IncreasingHD
cmsRun testL1TCaloLayer1EmulatorWithSpy.py setupString=patterns:ConstantEM=0xA7,DecreasingHD
```

etc.  The keywords IncreasingXX, DecreasingXX, ConstantXX, ConstantXX=0xAB 
where XX is EM or HD and 0xAB can be any byte-sized hex-number, can be
mixed to get various patterns.  Pattern string is case sensitive!

The syntax for pre-captured files is:

```bash
cmsRun testL1TCaloLayer1EmulatorWithSpy.py setupString=captures:/data/dasu/Layer1ZeroBiasCaptureData/TT_Feb22_252Events_v1
```

The syntax for ethernet captures is:

```bash
cmsRun L1Trigger/L1TCaloLayer1/test/testL1TCaloLayer1EmulatorWithSpy.py setupString=ethernet:PhiMapFile.xml maxEvents=16200
```

where PhiMapFile.xml is path to an XML file specifying CTP7 card to TCP/IP address mapping.

In case of problems, report to : dasu@hep.wisc.edu
