/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "rhoLaunderGibsonRSTM.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(rhoLaunderGibsonRSTM, 0);
addToRunTimeSelectionTable(RASModel, rhoLaunderGibsonRSTM, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

rhoLaunderGibsonRSTM::rhoLaunderGibsonRSTM
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, U, phi, transport, turbulenceModelName),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            coeffDict_,
            0.09
        )
    ),
    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            coeffDict_,
            0.41
        )
    ),
    Clg1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Clg1",
            coeffDict_,
            1.8
        )
    ),
    Clg2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Clg2",
            coeffDict_,
            0.6
        )
    ),
    Clg3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Clg3",
            coeffDict_,
            0.75
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            coeffDict_,
            1.44
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            coeffDict_,
            1.92
        )
    ),
    C3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3",
            coeffDict_,
            0.00
        )
    ),
    C4_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C4",
            coeffDict_,
            1.00
        )
    ),
    C5_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C5",
            coeffDict_,
            -0.33
        )
    ),
    Cs_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cs",
            coeffDict_,
            0.25
        )
    ),
    Ceps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps",
            coeffDict_,
            0.15
        )
    ),
    sigmaR_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaR",
            coeffDict_,
            0.81967
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            coeffDict_,
            1.3
        )
    ),
    C1Ref_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1Ref",
            coeffDict_,
            0.5
        )
    ),
    C2Ref_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2Ref",
            coeffDict_,
            0.3
        )
    ),
    couplingFactor_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "couplingFactor",
            coeffDict_,
            0.0
        )
    ),

    yr_(mesh_),

    R_
    (
        IOobject
        (
            "R",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateR("R", mesh_)
    ),
    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateK("k", mesh_)
    ),
    epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateEpsilon("epsilon", mesh_)
    ),
    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateNut("nut", mesh_)
    )
{
    bound(k_, kMin_);
    dimensionedScalar kMin2_("kMin2_", k_.dimensions(), 1.0E-14);
    bound(k_, kMin2_);
    bound(epsilon_, epsilonMin_);
    dimensionedScalar epsilonMin2_("epsilonMin2_", epsilon_.dimensions(), 1.0E-10);
    bound(epsilon_, epsilonMin2_);

    nut_ = Cmu_*sqr(k_)/epsilon_;
    nut_.correctBoundaryConditions();

    if (couplingFactor_.value() < 0.0 || couplingFactor_.value() > 1.0)
    {
        FatalErrorIn
        (
            "rhoLaunderGibsonRSTM::rhoLaunderGibsonRSTM"
            "(const volVectorField& U, const surfaceScalarField& phi,"
            "transportModel& transport)"
        )   << "couplingFactor = " << couplingFactor_
            << " is not in range 0 - 1" << nl
            << exit(FatalError);
    }

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> rhoLaunderGibsonRSTM::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            R_ - nu()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> rhoLaunderGibsonRSTM::divDevReff(volVectorField& U) const
{
    if (couplingFactor_.value() > 0.0)
    {
        return
        (
            fvc::div(R_ + couplingFactor_*nut_*fvc::grad(U), "div(R)")
          + fvc::laplacian
            (
                (1.0 - couplingFactor_)*nut_,
                U,
                "laplacian(nuEff,U)"
            )
          - fvm::laplacian(nuEff(), U)
        );
    }
    else
    {
        return
        (
            fvc::div(R_)
          + fvc::laplacian(nut_, U, "laplacian(nuEff,U)")
          - fvm::laplacian(nuEff(), U)
        );
    }
}


tmp<fvVectorMatrix> rhoLaunderGibsonRSTM::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    volScalarField muEff("muEff", rho*nuEff());

    if (couplingFactor_.value() > 0.0)
    {
        return
        (
            fvc::div
            (
                rho*R_ + couplingFactor_*(rho*nut_)*fvc::grad(U),
                "div((rho*R))"
            )
          + fvc::laplacian
            (
                (1.0 - couplingFactor_)*rho*nut_,
                U,
                "laplacian(muEff,U)"
            )
          - fvm::laplacian(muEff, U)
          - fvc::div(rho*nu()*dev2(T(fvc::grad(U))))
        );
    }
    else
    {
        return
        (
            fvc::div(rho*R_)
          + fvc::laplacian(rho*nut_, U, "laplacian(muEff,U)")
          - fvm::laplacian(muEff, U)
          - fvc::div(rho*nu()*dev2(T(fvc::grad(U))))
        );
    }
}


bool rhoLaunderGibsonRSTM::read()
{
    if (RASModel::read())
    {
        Cmu_.readIfPresent(coeffDict());
        kappa_.readIfPresent(coeffDict());
        Clg1_.readIfPresent(coeffDict());
        Clg2_.readIfPresent(coeffDict());
        Clg3_.readIfPresent(coeffDict());
        C1_.readIfPresent(coeffDict());
        C2_.readIfPresent(coeffDict());
        C3_.readIfPresent(coeffDict());
        C4_.readIfPresent(coeffDict());
        C5_.readIfPresent(coeffDict());
        Cs_.readIfPresent(coeffDict());
        Ceps_.readIfPresent(coeffDict());
        sigmaR_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());
        C1Ref_.readIfPresent(coeffDict());
        C2Ref_.readIfPresent(coeffDict());

        couplingFactor_.readIfPresent(coeffDict());

        if (couplingFactor_.value() < 0.0 || couplingFactor_.value() > 1.0)
        {
            FatalErrorIn("rhoLaunderGibsonRSTM::read()")
                << "couplingFactor = " << couplingFactor_
                << " is not in range 0 - 1"
                << exit(FatalError);
        }

        return true;
    }
    else
    {
        return false;
    }
}


void rhoLaunderGibsonRSTM::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    if (mesh_.changing())
    {
        yr_.correct();
    }

    volScalarField divU(fvc::div(phi_));

    const volScalarField& rho_ = mesh_.lookupObject<volScalarField>("rho");
    const surfaceScalarField& rhoPhi_ = mesh_.lookupObject<surfaceScalarField>("rhoPhi");
    const volScalarField& io_up2gradp_ = mesh_.lookupObject<volScalarField>("io_sigma");
    const volSymmTensorField& io_up2gradpt_ = mesh_.lookupObject<volSymmTensorField>("io_sigmat");

    volSymmTensorField P(-twoSymm(R_ & fvc::grad(U_)));
    volScalarField G(GName(), 0.5*mag(tr(P)));

    // Update epsilon and G at the wall
    epsilon_.boundaryField().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(rho_, epsilon_)
      + fvm::div(rhoPhi_, epsilon_)
    //- fvm::laplacian(Ceps*(K/epsilon_)*R, epsilon_)
      - fvm::laplacian(rho_*DepsilonEff(), epsilon_)
      ==
        fvm::SuSp(C1_*rho_*G/k_, epsilon_)
      + fvm::SuSp(-C5_*rho_*divU, epsilon_)
      - fvm::Sp(C2_*rho_*epsilon_/k_, epsilon_)
      + fvm::SuSp(-C4_*io_up2gradp_/k_, epsilon_)
    );

    epsEqn().relax();

    epsEqn().boundaryManipulate(epsilon_.boundaryField());

    solve(epsEqn);
    bound(epsilon_, epsilonMin_);
    dimensionedScalar epsilonMin2_("epsilonMin2_", epsilon_.dimensions(), 1.0E-10);
    bound(epsilon_, epsilonMin2_);


    // Reynolds stress equation

    const fvPatchList& patches = mesh_.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isA<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];
                P[faceCelli] *=
                    min(G[faceCelli]/(0.5*mag(tr(P[faceCelli])) + SMALL), 1.0);
            }
        }
    }

    const volSymmTensorField reflect
    (
        C1Ref_*epsilon_/k_*R_ - C2Ref_*Clg2_*dev(P)
    );

    tmp<fvSymmTensorMatrix> REqn
    (
        fvm::ddt(rho_, R_)
      + fvm::div(rhoPhi_, R_)
      - fvm::laplacian(rho_*Cs_*(k_/epsilon_)*R_, R_)
      //- fvm::laplacian(rho_*DREff(), R_)
      + fvm::Sp(Clg1_*rho_*epsilon_/k_, R_)
      ==
        rho_*P
      - (2.0/3.0*(1 - Clg1_)*I)*rho_*epsilon_
      - Clg2_*rho_*dev(P)
      - Clg3_*dev(-io_up2gradpt_)
      - io_up2gradpt_

        // wall reflection terms
      + symm
        (
            I*((yr_.n() & reflect) & yr_.n())
          - 1.5*(yr_.n()*(reflect & yr_.n())
          + (yr_.n() & reflect)*yr_.n())
        )*pow(Cmu_, 0.75)*rho_*pow(k_, 1.5)/(kappa_*yr_*epsilon_)
    );

    REqn().relax();

    REqn().boundaryManipulate(R_.boundaryField());
    solve(REqn);

    R_.max
    (
        dimensionedSymmTensor
        (
            "zero",
            R_.dimensions(),
            symmTensor
            (
                kMin_.value(), -GREAT, -GREAT,
                kMin_.value(), -GREAT,
                kMin_.value()
            )
        )
    );

    k_ == 0.5*tr(R_);
    bound(k_, kMin_);
    dimensionedScalar kMin2_("kMin2_", k_.dimensions(), 1.0E-14);
    bound(k_, kMin2_);


    // Re-calculate turbulent viscosity
    nut_ = Cmu_*sqr(k_)/epsilon_;
    nut_.correctBoundaryConditions();


    // Correct wall shear stresses

    /*forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isA<wallFvPatch>(curPatch))
        {
            symmTensorField& Rw = R_.boundaryField()[patchi];

            const scalarField& nutw = nut_.boundaryField()[patchi];

            const vectorField snGradU(U_.boundaryField()[patchi].snGrad());

            const vectorField& faceAreas
                = mesh_.Sf().boundaryField()[patchi];

            const scalarField& magFaceAreas
                = mesh_.magSf().boundaryField()[patchi];

            forAll(curPatch, facei)
            {
                // Calculate near-wall velocity gradient
                tensor gradUw
                    = (faceAreas[facei]/magFaceAreas[facei])*snGradU[facei];

                // Calculate near-wall shear-stress tensor
                tensor tauw = -nutw[facei]*2*symm(gradUw);

                // Reset the shear components of the stress tensor
                Rw[facei].xy() = tauw.xy();
                Rw[facei].xz() = tauw.xz();
                Rw[facei].yz() = tauw.yz();
            }
        }
    }*/
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
