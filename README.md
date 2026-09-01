# Mechanobiological FEBio Plugin

A custom FEBio material and reaction framework for the mechanobiological simulation of intervertebral disc (IVD) degeneration.

This repository contains a C++ implementation of a custom multiphasic material extending the FEBio framework to couple **mechanics, solute transport, cellular metabolism, extracellular matrix turnover, fixed charge density evolution, and degeneration-dependent tissue properties**.

The implementation was developed for research purposes in the context of computational mechanobiology of the intervertebral disc.

---

## Overview

Intervertebral disc degeneration is a long-term biological process resulting from interactions between mechanical loading, nutrient availability, cellular metabolism, extracellular matrix turnover, and changes in tissue composition.

The present implementation extends FEBio's multiphasic framework to represent these coupled processes within a finite element formulation.

The model includes:

* multiphasic material behavior;
* solute transport;
* glucose and oxygen metabolism;
* aerobic and anaerobic cellular metabolism;
* degeneration-dependent cell density;
* extracellular matrix biosynthesis and degradation;
* evolution of the referential fixed charge density (FCD);
* degeneration-dependent solute diffusivity;
* degeneration-dependent oriented collagen fiber (OCF) stiffness;
* a global degeneration grade derived from the local degeneration state of the nucleus pulposus.

The implementation is intended to provide a mechanistic framework for long-term simulation of IVD degeneration.

---

## Model architecture

The mechanobiological formulation couples several physical and biological processes:

```text
Mechanical deformation
        │
        ▼
Fluid and solute transport
        │
        ├───────────────┐
        ▼               ▼
Glucose availability   Oxygen availability
        │               │
        └───────┬───────┘
                ▼
        Cellular metabolism
                │
                ▼
           ATP production
                │
        ┌───────┴────────┐
        ▼                ▼
   Biosynthesis      Maintenance /
        │             degradation
        └───────┬────────┘
                ▼
      Referential FCD evolution
                │
                ▼
       Degeneration state
                │
        ┌───────┴────────┐
        ▼                ▼
 Diffusivity          OCF stiffness
```

The degeneration grade is treated as a global IVD-level variable. It is determined from the most advanced local degeneration state within the nucleus pulposus and subsequently used by the mechanobiological constitutive relationships.

---

## Main components

### 1. Mechanobiological multiphasic material

`FEMechanobiological` extends FEBio's `FEMultiphasic` material.

The material introduces a long-term evolution law for the referential fixed charge density:

$$
\frac{d FCD_r}{dt} = -S_{\mathrm{bio}} + D_{\mathrm{deg}}
$$

where:

* \(FCD_r\) is the referential fixed charge density;
* \(S_{\mathrm{bio}}\) is the biosynthesis contribution;
* \(D_{\mathrm{deg}}\) is the degradation contribution.

The material also stores history variables through the custom material-point class `FEMechanobiologicalPoint`.

---

### 2. Biosynthesis

The biosynthesis formulation accounts for:

* glucose availability;
* oxygen availability;
* cellular density;
* tissue region;
* metabolic activity;
* degeneration grade;
* ATP availability.

The implementation uses an ATP-based formulation to link cellular metabolism to extracellular matrix biosynthesis.

The parameter `delta` controls the degeneration-dependent reduction in biosynthetic efficiency.

---

### 3. Degradation

The degradation contribution is modeled as a degeneration-dependent loss of referential fixed charge density.

The implementation includes time-dependent and grade-dependent terms controlling the degradation rate.

The degeneration grade is bounded between:

```text
G0 = 1.5
Gmax = 5.0
```

---

### 4. Aerobic metabolism

`Aerobic_Reaction` implements oxygen consumption associated with aerobic cellular metabolism.

The reaction rate depends on:

* glucose concentration;
* oxygen concentration;
* effective cell density;
* tissue region;
* degeneration grade.

Analytical reaction tangents are implemented with respect to:

* deformation;
* effective solute concentration;
* pressure.

The deformation and concentration derivatives are provided to FEBio's nonlinear solver to improve the consistency of the Jacobian.

---

### 5. Anaerobic metabolism

`Anaerobic_Reaction` implements glucose consumption associated with anaerobic metabolism.

The model accounts for:

* basal lactate production;
* oxygen-dependent metabolic partitioning;
* maintenance metabolism;
* glucose-dependent effective cell density;
* tissue-specific metabolic parameters.

The glucose consumption rate is related to lactate production through the implemented stoichiometric relationship:

$$
Q_{\mathrm{glucose}} = \frac{Q_{\mathrm{lactate}}}{2}
$$

Analytical tangents are provided with respect to deformation and effective solute concentrations.

---

### 6. Degeneration-dependent diffusivity

`FEDiffusivity_Calcification` extends FEBio's `FESoluteDiffusivity`.

The diffusivity is modified according to the local degeneration grade, which is used as a surrogate for cartilage endplate calcification.

The implementation provides:

* free diffusivity;
* degeneration-dependent diffusivity;
* diffusivity validation;
* analytical tangent interfaces required by FEBio.

---

### 7. Degeneration grade controller

`Grade` is a FEBio load controller that computes the global degeneration grade from the local referential fixed charge density.

The controller:

1. identifies mechanobiological material domains;
2. restricts the search to the nucleus pulposus;
3. evaluates the local referential FCD at integration points;
4. identifies the maximum local degeneration;
5. converts the local FCD reduction into a degeneration grade;
6. returns the resulting global grade.

The grade is bounded between 1.5 and 5.

---

### 8. OCF stiffness controller

`OCF_stiffness_eff` computes an effective oriented collagen fiber (OCF) stiffness as a function of:

* lamellar position;
* healthy inner and outer AF stiffness;
* degeneration grade.

This allows the OCF mechanical contribution to evolve with degeneration.

---

## Tissue regions

The implementation distinguishes three tissue regions:

| Region | Tissue                   |
| -----: | ------------------------ |
|    `0` | Nucleus pulposus (NP)    |
|    `1` | Annulus fibrosus (AF)    |
|    `2` | Cartilage endplate (CEP) |

Region-dependent parameters are used for cellular density and metabolic activity.

---

## Source files

Currently, the repository contains the main source and header files:

```text
.
├── Mechanobiological.cpp
├── Mechanobiological.h
└── README.md
```

### `Mechanobiological.h`

Contains:

* material-point definitions;
* mechanobiological material class;
* degeneration controllers;
* aerobic reaction class;
* anaerobic reaction class;
* degeneration-dependent diffusivity class.

### `Mechanobiological.cpp`

Contains the implementation of:

* material initialization;
* solid-bound molecule updates;
* referential FCD evolution;
* biosynthesis;
* degradation;
* degeneration controllers;
* aerobic metabolism;
* anaerobic metabolism;
* analytical reaction tangents;
* degeneration-dependent diffusivity.

---

## Requirements

The implementation was developed and tested with:

* **FEBio 4.4**
* **Visual Studio 2022 Community**
* **C++**
* **FEBio SDK / FEBio source libraries**

The implementation is intended for the FEBio 4.4 API and may require modifications to compile against other FEBio versions.

FEBio is an open-source nonlinear finite element solver developed specifically for biomechanical applications. The official FEBio repository and releases are available on GitHub.


## Citation

The source code is archived on Zenodo and can be cited using the following DOI:

https://doi.org/10.5281/zenodo.22228993

The DOI corresponds to version 1.0.0 of the code used in the present study.
