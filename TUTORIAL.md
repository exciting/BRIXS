Tutorial: F K edge in LiF
============================================================================
First, we create a directory for the BSE calculations. All subsequent
calculations with the **exciting** code will be performed in that directory

```
mkdir lif-exciting
cd lif-exciting/
```

F K Edge BSE in LiF
------------------------------------------------------------------------------
First, we need to calculate the DFT electronic structure and F K edge absorption
spectrum. We do this with the following input file:

```xml
<input>
  <title>LiF F K Edge</title>
  <structure speciespath="./">
    <crystal scale="1.8897543768634038">
      <basevect>      0.00000000       2.01297060       2.01297060</basevect>
      <basevect>      2.01297060       0.00000000       2.01297060</basevect>
      <basevect>      2.01297060       2.01297060       0.00000000</basevect>
    </crystal>
    <species speciesfile="Li.xml">
      <atom coord="      0.00000000       0.00000000       0.00000000" />
    </species>
    <species speciesfile="F.xml">
      <atom coord="      0.50000000       0.50000000       0.50000000" />
    </species>
  </structure>
  <groundstate
    do="fromscratch"
    ngridk="2 2 2"
    xctype="GGA_PBE_SOL"/>
 <xs
      xstype="BSE"
      ngridk="2 2 2"
      ngridq="2 2 2"
      rgkmax="4.0" 
      vkloff="0.05 0.15 0.25"
      nempty="30" 
      gqmax="2.0" 
      broad="0.0276"
      tappinfo="true"
      tevout="true">

      <energywindow 
         intv="23.6 24.6" 
         points="1500" />
 
      <screening 
         screentype="full" 
         nempty="100"/>
 
      <BSE 
        chibar0="false"
        chibar0comp="1"
	      xas="true"
        xasspecies="2"
	      xasatom="1"
        xasedge="K"
        bsetype="singlet" 
	      nstlxas="1 2"/>
 
      <qpointset>
      <qpoint>0.0 0.0 0.0</qpoint>
      </qpointset>
       <storeexcitons MinNumberExcitons="1" MaxNumberExcitons="32"/>
   </xs>
</input>
```

We see that we write all 32 excitonic eigenstates to file. In a converged BSE
calculation, the number of eigenstates is huge and one might need to consider
carefully how many are necessary for a converged RIXS spectrum. Generally, this
will depend on the range of excitation energies one is interested in. The higher
the excitation energy is beyond the absorption edge, the more eigenstates are
required.

Executing the **exciting** code with this input.xml will generate a output file
named **bse_output.h5** which we rename by

```
mv bse_output.h5 core_output.h5
```

Renaming is important, because otherwise it will be overwritten when we
calculate the valence excitation spectrum.

Momentum Matrix Elements
------------------------------------------------------------------------------
An important ingredient of the RIXS calculation are the momentum matrix elements
between the core and conduction states. These can be calculated with low
computational cost by adjusting the input files such that it looks like:

```xml
<input>
  <title>LiF F K Edge</title>
  <structure speciespath="./">
    <crystal scale="1.8897543768634038">
      <basevect>      0.00000000       2.01297060       2.01297060</basevect>
      <basevect>      2.01297060       0.00000000       2.01297060</basevect>
      <basevect>      2.01297060       2.01297060       0.00000000</basevect>
    </crystal>
    <species speciesfile="Li.xml">
      <atom coord="      0.00000000       0.00000000       0.00000000" />
    </species>
    <species speciesfile="F.xml">
      <atom coord="      0.50000000       0.50000000       0.50000000" />
    </species>
  </structure>
  <groundstate
    do="fromscratch"
    ngridk="2 2 2"
    rgkmax="4.0"
    xctype="GGA_PBE_SOL"/>
 <xs
      xstype="BSE"
      ngridk="2 2 2"
      ngridq="2 2 2"
      rgkmax="4.0" 
      vkloff="0.05 0.15 0.25"
      nempty="30" 
      gqmax="2.0" 
      broad="0.0276"
      tappinfo="true"
      tevout="true">

      <energywindow 
         intv="23.6 24.6" 
         points="1500" />
 
      <screening 
         screentype="full" 
         nempty="100"/>
 
      <BSE 
        chibar0="false"
        chibar0comp="1"
	      xas="true"
        xasspecies="2"
	      xasatom="1"
        xasedge="K"
        bsetype="singlet" 
	      nstlxas="1 2"
        distribute="true"/>
 
      <qpointset>
      <qpoint>0.0 0.0 0.0</qpoint>
      </qpointset>
      <plan>
        <doonly task="writepmatxs"/>
        <doonly task="writepmatasc"/>
      </plan>
   </xs>
</input>
```

The input file is nearly identical to the one used for the calculation of the F K edge,
with the only change being the subelement `plan` which triggers the calculation
of the momentum matrix elements (`writepmatxs`) and the output to HDF5 file
(`writepmatasc`). The matrix elements are written to the file **bse_output.h5**.
We rename the file again

```
mv bse_output.h5 pmat.h5
```

Valence Spectrum in LiF
------------------------------------------------------------------------------
Following the F K edge spectrum, we now calculate the valence excitation
spectrum. While many parameters (such as `nstlbse` and `gqmax`) can be different
between the two BSE calculations, the parameters `ngridk`, `ngridq`, and
`vkloff` **HAVE TO BE IDENTICAL**. To make sure that the wavefunctions are
identical in both calculations, we do not calculate the DFT electronic structure
again (a repeated diagonalization of the Kohn-Sham Hamiltonian can introduce a
different global phase, which alters the results). As such, the input file looks
like this

```xml

<input>
  <title>LiF Valence Spectrum</title>
  <structure speciespath="./">
    <crystal scale="1.8897543768634038">
      <basevect>      0.00000000       2.01297060       2.01297060</basevect>
      <basevect>      2.01297060       0.00000000       2.01297060</basevect>
      <basevect>      2.01297060       2.01297060       0.00000000</basevect>
    </crystal>
    <species speciesfile="Li.xml">
      <atom coord="      0.00000000       0.00000000       0.00000000" />
    </species>
    <species speciesfile="F.xml">
      <atom coord="      0.50000000       0.50000000       0.50000000" />
    </species>
  </structure>
  <groundstate
    do="fromscratch"
    ngridk="2 2 2"
    rgkmax="4.0"
    xctype="GGA_PBE_SOL"/>
 <xs
      xstype="BSE"
      ngridk="2 2 2"
      ngridq="2 2 2"
      rgkmax="4.0" 
      vkloff="0.05 0.15 0.25"
      nempty="30" 
      gqmax="2.0" 
      broad="0.0276"
      tappinfo="true"
      tevout="true">

      <energywindow 
             intv="0 1.0" 
             points="1500" />

      <screening 
              screentype="full" 
              nempty="100" />

      <BSE
        chibar0="false"
        chibar0comp="1"
        bsetype="singlet" 
        nstlbse="2 5 1 2"
        distribute="true" />

      <qpointset>
                  <qpoint>0.0 0.0 0.0</qpoint>
       </qpointset>
       <storeexcitons MinNumberExcitons="1" MaxNumberExcitons="64"/>
       <plan>
         <doonly task="writepmatxs"/>
         <doonly task="scrcoulint"/>
         <doonly task="exccoulint"/>
         <doonly task="bse"/>
       </plan>
   </xs>
</input>

```

Once again, we write all 64 eigenstates of the BSE to file. With the `plan`
subelement, we specify which tasks need to be performed. We re-calculate neither
the KS electronic structure nor the screening. These quantities are supossed to
be identical in the two BSE calculations. Finally, we also rename this file

```
mv bse_output.h5 optical_output.h5
```
From this point on, only the files **core_output.h5**, **pmat.h5**, and
**optical_output.h5** are required. We will move them to a different directory

```
mkdir ../lif-rixs
mv *.h5 ../lif-rixs
cd ../lif-rixs
```

RIXS Spectrum LiF
------------------------------------------------------------------------------
Now we have to create the file **input.cfg** that describes the input parameters
for the **BRIXS** calculation. The file contains the following:

```
omega= 653 654 655 656
pol=1.0 0.0 0.0
broad=0.15
nblocks=2
eigstates_optical=64
eigstates_core=32
```
Here, we specify the excitation energies (`omega`) and the polarization (`pol`)
of the incoming photon. Note that the polarization should be identical to the
one specified by the parameter `chibar0comp` in the **exciting** calculations.
We furthermore specify the lifetime broadening (`broad`) of the intermediate
state. With `eigstates_optical` and `eigstates_core` we specify how many
eigenstates of the corresponding BSE calculations are used. Finally, `nblocks`
defines in how many blocks the eigenvectors should be separated. More blocks
will reduce the memory requirements of the calculations, but might increase the
runtime. A good estimate is to use the number of k-points along a specific
direction.

We then execute the first **BRIXS** executable with

```
/path/to/BRIXS/bin/rixs-pathway-serial
```

which creates the file **data.h5**. The final result is then obtained by
executing

```
/path/to/BRIXS/bin/rixs-oscstr-serial
```

The final result is stored in the file **rixs.h5**.
