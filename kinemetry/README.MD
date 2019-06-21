This distribution contains an IDL (http://www.exelisvis.com/ProductsServices/IDL.aspx)
implementation of the KINEMETRY method developed by
Krajnovic et al. 2006 MNRAS, 366, 787.

-------------
SAURON velocity and velocity dispersion measurements and associated errors for NGC2974 are
given in file NGC2974_SAURON_kinematics.dat. There are seven columns: number of bin, X and Y bin coordinates, velocity, velocity error, velocity dispersion and velocity dispersion error.

These data correspond to data published in Emsellem et al. (2004
MNRAS, 352, 271). These data are not part of the official SAURON data
release and may differ slightly from the final products.


To learn how to use KINEMETRY on velocity maps, run the example program
kinemetry_example_velocity.pro or read the detailed documentation at the top of
the file kinemetry.pro

To run kinemetry examples:
1. put all files in the directory 'kinemetry'
2. cd to that directory (within IDL: CD, kinemetry')
3. compile kinemetry_example_velocity.pro routine (.r kinemetry_example_velocity.pro)
4. type kinemetry_example_velocity at the IDL prompt

To learn how to use KINEMETRY on photometric images, run the example program kinemetry_example_phot.pro. The program will run kinemetry on a SLOAN r-band image of NGC2974 obtained at the Isaac Newton Telescope on La Palma, as part of the ATLAS3D Project (Cappellari et al. 2011, MNRAS, 413, 813), and published as part of Scott et al. 2013, MNRAS, 432, 1894.

For an example of running KINEMETRY on velocity dispersion maps, execute kinemetry_example_sigma.pro. 


The procedure KINEMETRY needs the following extra IDL routines:
- MPFIT: by C.B. Markwardt -> http://cow.physics.wisc.edu/~craigm/idl/

The IDL Astronomy User's Library is assumed to be installed
(http://idlastro.gfsc.nasa.gov/) and in the IDL path, as well as Coyote IDL Program Libraries (http://www.idlcoyote.com/documents/programs.php). In particular
KINEMETRY_EXAMPLE_ programs use routines RDFLOAT.PRO, PLOTERROR.PRO, MRDFITS.PRO and DIST_CIRCLE.PRO.


Next to the actual routines there are also .sav versions of the example files:

kinemetry_example_velocity_input.sav
kinemetry_example_sigma_input.sav
kinemetry_example_phot_input.sav

These .sav files allow you to execute example routine if you have an access to the IDL virtual machine only. You will need to pass the name of the input file (including the full path) with the data kinemetry will analyse. To run the .sav files, follow the procedure at the IDL prompt: 

restore, ‘kinemetry_example_velocity_input.sav’
kinemetry_example_velocity_input, ‘NGC2974_SAURON_kinematics.dat’

restore, ‘kinemetry_example_sigma_input.sav’
kinemetry_example_sigma_input, ‘NGC2974_SAURON_kinematics.dat’

restore, ‘kinemetry_example_phot_input.sav’
kinemetry_example_phot_input, ‘NGC2974r.fits’



The program was tested on IDL 7.1.

-------------

Davor Krajnovic
Potsdam, 23 June 2016


