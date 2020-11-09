BRIXS: Calculating Resonant Inelastic X-Ray Scattering Spectra from
Bethe-Salpeter Equation Calculations
============================================================================

This code determines resonant inelastic x-ray scattering (RIXS) spectra for
solids, starting from Bethe-Salpeter Equation (BSE) calculations in an
all-electron many-body framework. 

Theoretical Background
------------------------------------------------------------------------------
In RIXS, an x-ray photon with energy $`\omega_1`$ is absorbed, which leads to the
excitation of a core electron to a conduction band. Coherently, an x-ray photon
with energy $`\omega_2`$ is emitted as a valence electron fills the core hole. The
final state of the scattering process contains a valence hole and an excited
electron in a conduction band. The energy loss $`\omega_2-\omega_1`$ is
transferred to the electronic system. In RIXS experiments, the
double-differential cross section (DDCS) $`\mathrm{d}^2 \sigma/\mathrm{d}\Omega_2
\mathrm{d}\omega_2`$ is measured, which describes the rate of scattering for
photons with final energy $`[\omega_2,\omega_2+\mathrm{d}\omega_2]`$ in a solid
angle range $`[\Omega_2,\Omega_2+\mathrm{d}\Omega_2]`$.

Following the derivation presented in Phys. Rev. 2, 042003(R) (2020), the RIXS
double-differential cross section (DDCS) is defined as

```math
\frac{\mathrm{d}^2 \sigma}{\mathrm{d} \Omega_2 \mathrm{d}\omega_2}= \alpha^4
\left( \frac{\omega_2}{\omega_1} \right) \mathrm{Im} \sum_{\lambda_o}
\frac{|t^{(3)}_{\lambda_o}(\omega_1)|^2}{(\omega_1-\omega_2)-E^{\lambda_o}+\mathrm{i}\eta},
```

where $`E^{\lambda_o}`$ are the valence excitation energies of the system, and
$\alpha$ is the fine-structure constant. The DDCS is an explicit function of the
energy loss $`\omega_2-\omega_1`$, while the oscillator strength
$`t^{(3)}_{\lambda_o}(\omega_1)`$ depends solely on the excitation energy
$`\omega_1`$. The oscillator strenght is defined as

```math
t^{(3)}_{\lambda_o}(\omega_1) = \sum_{\lambda_c} \frac{t^{(2)}_{\lambda_o,
\lambda_c}t^{(1)}_{\lambda_c}}{\omega_1-E^{\lambda_c}+\mathrm{i}\eta},
```
where $`E^{\lambda_c}`$ are the core excitation energies of the system. The *core
excitation oscillator strength* $`t^{(1)}_{\lambda_c}`$ is defined as

```math
t^{(1)}_{\lambda_c}=\sum_{c \mu \mathbf{k}} X_{c \mu \mathbf{k},
\lambda_c}\left[ \mathbf{e}_1 \cdot \mathbf{P}_{c \mu \mathbf{k}} \right],
```
where $`X_{c \mu \mathbf{k}, \lambda_c}`$ is the BSE eigenstate corresponding to
the excitation energy $`E^{\lambda_c}`$, $`\mathbf{e}_1`$ is the polarization vector
of the incoming photon, and $`\mathbf{P}_{c \mu \mathbf{k}}=\langle c \mathbf{k}
| \hat{\mathbf{p}} | \mu \mathbf{k} \rangle`$ are the momentum matrix elements.

The *excitation pathway* $`t^{(2)}_{\lambda_o,\lambda_c}`$ is defined as\
```math
t^{(2)}_{\lambda_o,\lambda_c}=\sum_{c v \mathbf{k}} \sum_{\mu} X_{cv \mathbf{k},
\lambda_o}\left[ \mathbf{e}_2^* \cdot \mathbf{P}_{\mu v \mathbf{k}} \right]
\left[ X_{c \mu \mathbf{k}, \lambda_c} \right]^*
```

Code Usage
-----------------------------------------------------------------------------
To perform calculations with **BRIXS**, you first need to calculate two BSE
calculations to obtain the core excitation eigenstates $`X_{c \mu \mathbf{k},
\lambda_c}`$ and the valence excitation ones $`X_{c v \mathbf{k}, \lambda_o}`$.
Both calculations have to be performed on the same $`\mathbf{k}`$-grid, but the
number of empty states in the calculations can differ. 

**NOTE**: Here, we are interested in the eigenstates of the polarizability. In a
typical BSE calculation, one is only interested in the dielectric function. As
such, the eigenvectors from these calculations are typically obtained from a
reduced BSE Hamiltonian, using the adjusted Coulomb potential $`\bar{v}`$
[compare Eqs. 6.57 and 6.58 in Christian Vorwerk's PhD thesis and the further
discussion there]. Make sure that you calculate the eigenstates of the actual
BSE Hamiltonian when you create eigenstates for the RIXS spectra.

To use the **BRIXS** code, first perform a valence- and core-level BSE
calculation with the **exciting** code. To calculate the matrix elements of the
polarizability, add the following part to the *input.xml* 

```xml
chibar="true"
chibar0comp="1"

```
The first attribute triggers the calculation of the full polarizability at
vanishing momentum transfer, while the second one defines for which the BSE is
solved. The additional dependence occurs since we calculate the analytical limit
for the exchange interaction, *i.e.*

```math
V_{c \mu \mathbf{k}, c' \mu' \mathbf{k}'}(\mathbf{q}=0)=
V_{c \mu \mathbf{k}, c' \mu' \mathbf{k}'}(\mathbf{q}=0, \mathbf{G}=0)
+ \sum_{\mathbf{G}} v_{\mathbf{G}}(\mathbf{q})M^*_{c \mu
  \mathbf{k}}(\mathbf{G},\mathbf{q}=0) M_{c' \mu'
\mathbf{k}'}(\mathbf{G},\mathbf{q}=0),
```

where the first term in the sum is needed to tame the divergence of the Coulomb
potential as $`\lim_{q \rightarrow 0} v_0(\mathbf{q}) \propto 1/\mathbf{q}^2`$,
since

```math
\lim_{q \rightarrow 0} M_{c \mu \mathbf{k}}(\mathbf{q},\mathbf{G}=0)=
\lim_{q \rightarrow 0} \langle c \mathbf{k} |
\mathrm{e}^{-\mathrm{i}\mathrm{q}\mathbf{r}} | \mu \mathbf{k} \rangle=
-\mathrm{i}\mathbf{q}\langle c \mathbf{k} |\mathbf{r}| \mu \mathbf{k} \rangle=
-\mathrm{i}\mathbf{q}\frac{\langle c \mathbf{k}|\mathbf{p}| \mu
\mathbf{k}}{\epsilon_{c \mathbf{k}}-\epsilon_{\mu \mathbf{k}}}
```

This yields

```math
V_{c \mu \mathbf{k}, c' \mu' \mathbf{k}'}(\mathbf{q}=0, \mathbf{G}=0)=
\frac{1}{V} \mathrm{e}_q \frac{\langle c \mathbf{k} |\mathbf{p}| \mu \mathbf{k}
\rangle^* \langle c' \mathbf{k}' |\mathbf{p}| \mu'
\mathbf{k}'}{\left(\epsilon_{c \mathbf{k}}-\epsilon_{\mu \mathbf{k}}
\right)\left( \epsilon_{c' \mathbf{k}'}-\epsilon_{v' \mathbf{k}'} \right)}
```
As we can see, this matrix element depends on the direction $`\mathbf{e}_q`$ of
the momentum transfer as it goes to zero. Therefore, the BSE calculation depends
on the component of the dielectric tensor that one wants to calculate. 

**NOTE**: When one calculates the dielectric function as described above, only
the component of the dielectric tensor specified by `chibar0comp` will agree
with the corresponding component obtained from a calculation with
`chibar="false"`. The other components will differ and will be incorrect.

Once the BSE calculations are performed, the HDF5 files *optical_output.h5*,
*core_output.h5*, and *pmat.h5* have to be gathered in a directory together with
the input file *input.cfg*. Then, the execution of *rixs-pathway* creates the
file *data.h5*. Once  this file is present, the execution of *rixs-oscstr*
creates the final output in the file *rixs.h5*. A more detailed tutorial is
provided in the file *TUTORIAL.md*

