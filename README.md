# RegCM 4.5 - CB6C

RegCM v4.5 - Gas Phase Module CB6r2vC Configured to run as an option instead of CMBZ

## Usage

config and make as you would with your standard RegCM4 (https://www.ictp.it/research/esp/models/regcm4.aspx)

To use simply input CB6C instead of CBMZ as 'chemsimtype' in the 'chemparam' stanza of your namelist.

You can modify emission and boundary condition inputs in the files below, but you will need to run config and make again.

emissions: Main/chemlib/mod_che_ncio.F90

CHBC: PreProc/ICBC/mod_ch_icbc.F90 or (mod_ch_icbc_clim.F90)

More detailed information and a tutorial can be found in the presentations directory. Please cite the paper below when using this model.

## Reference
Ciarlo`, J.M., Aquilina, N.J., Strada, S. et al. A modified gas-phase scheme for advanced regional climate modelling with RegCM4. Clim Dyn (2021). https://doi.org/10.1007/s00382-021-05722-y


## Contact

For any queries please contact jciarlo[at]ictp.it
