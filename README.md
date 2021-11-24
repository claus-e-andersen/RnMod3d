RnMod3d
=======

Finite-volume radon transport model (porous media)

This numerical model is written in the programming language Pascal and it was originally designed for
the compilers Turbo Pascal 7.0 or Delphi. However, the code runs without problems (and modifications) with Free Pascal. 

To run it with Free Pascal I have tried the folowing on a Linux machine. I have tried two installations. In 2014, I used: openSUSE 13.1 running Free Pascal vers. 2.6.4 and the Lazarus IDE version 1.2.4. In 2021, I used Ubuntu 20.04 with Lazarus 2.0.6 (2019-12-15). Lazarus could not be installed using the Ubuntu software center, so I used: 
sudo apt-get install Lazarus.  


(1) Put the RNMOD3D files in a single folder. (2) Open Lazarus. (3) Go to Tools / "Convert Delphi project to Lazarus project" and select one of the .pdr-files (e.g. F0100prg.dpr). Lazarus will then automatically make a Lazarus .lpr-file. (4) From the Run menu item select Compile/build/run. (5) To see the output, open the terminal window using: View / Debug windows / Terminal output. (6) The output (such as f0100LOG.dat) is in the a folder such as RnMod3d/lib/x86_64-linux/. Use the linux command
find . -name f0100LOG.dat to find such files.

It's really easy.

- Claus
- August 10, 2014 + November 24, 2021


The manual for the code is the report from Risø National Laboratory: Risø-R-1201(EN)

Title:
    Radon transport modelling: User's guide to RnMod3d

Authors:
    Andersen, Claus Erik

Abstract:
    RnMod3d is a numerical computer model of soil-gas and radon transport in porous media. It can be used, for example, to study radon entry from soil into houses in response to indoor-outdoor pressure differences or changes in atmospheric pressure. It canalso be used for flux calculations of radon from the soil surface or to model radon exhalation from building materials such as concrete. The finite-volume model is a technical research tool, and it cannot be used meaningfully without good understandingof the involved physical equations. Some understanding of numerical mathematics and the programming language Pascal is also required. Originally, the code was developed for internal use at Risø only. With this guide, however, it should be possible forothers to use the model. Three-dimensional steady-state or transient problems with Darcy flow of soil gas and combined generation, radioactive decay, diffusion and advection of radon can be solved. Moisture is included in the model, and partitioning ofradon between air, water and soil grains (adsorption) is taken into account. Most parameters can change in time and space, and transport parameters (diffusivity and permeability) may be anisotropic. This guide includes benchmark tests based on simpleproblems with known solutions. RnMod3d has also been part of an international model intercomparison exercise based on more complicated problems without known solutions. All tests show that RnMod3d gives results of good quality.

Risø-R-1201, Risø-R-1201(EN)

