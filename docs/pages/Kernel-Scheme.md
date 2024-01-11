<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
<!-- SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder -->

# Kernel Scheme

five-part naming structure:

- first field: type of collision operator - one letter
- second field: year of release - two numbers
- third field: compressible/ incompressible/ selective
- fourth field: model
- fifth field (optional): Other

## First field - one letter:

- C Central Moments/ Cascade
- K Cumulant
- B BGK
- F Factorized
- M MRT

## Third field:
  
- Incompressible
- Compressible
- Selective

## Fourth field:

- NavierStokes
- AdvectionDiffusion
- Bingham

## Fifth field - optional with any length

e.g. Smagorinski, SharpInterfaceImmersedBoundary, PorousMedia, etc.


## Examples

| Koll.Operator | Year 	| compressible?  | Model                | MISC                           |
|--------------	|------	|----------------|----------------------|--------------------------------|
| K            	| 17   	| Compressible   | NavierStokes         | Smagorinski                    |
| C            	| 15   	| Selective      | AdvectionDiffusion   | SharpInterfaceImmersedBoundary |
| B            	| 15   	| Incompressible | Bingham              |                                |
| F            	| 23   	| Compressible   | NavierStokes         |                                |
| M            	| 02   	| Compressible   | NavierStokes         |                                |


- K17IncompressibleNavierStokes
- C23IncompressibleBinghamSharpInterfaceImmersedBoundary
- F99CompressibleNavierStokesSmagorinsky