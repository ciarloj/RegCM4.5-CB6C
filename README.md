# RegCM 4.5 - CB6C

RegCM v4.5 - Gas Phase Module CB6vC Configured to run as an option instead of CMBZ

## Usage

To use simply input CB6C instead of CBMZ as 'chemsimtype' in the 'chemparam' stanza of your namelist.

You can modify emission and boundary condition inputs in the files below, but you will need to run config and make again.

-emissions: Main/chemlib/mod_che_ncio.F90

-CH BC    : PreProc/ICBC/mod_ch_icbc.F90 or (mod_ch_icbc_clim.F90)

## Changelog



## Contact

For any queries please contact jciarlo[at]ictp.it
