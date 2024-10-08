/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "sixDoFRigidBodyMotionSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"
#include "pointPatchDist.H"
#include "pointConstraints.H"
#include "uniformDimensionedFields.H"
#include "forces.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sixDoFRigidBodyMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        sixDoFRigidBodyMotionSolver,
        dictionary
    );
}

//---ModMorph{
// * * * * * Private member functions * * * * * //
void Foam::sixDoFRigidBodyMotionSolver::cosineTransition(pointScalarField& scale)
{
    // Convert the scale function to a cosine
        scale.primitiveFieldRef() =
            min
            (
                max
                (
                    0.5
                  - 0.5
                   *cos(scale.primitiveField()
                   *Foam::constant::mathematical::pi),
                    scalar(0)
                ),
                scalar(1)
            );
        return;
};
//---ModMorph}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionSolver::sixDoFRigidBodyMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    motion_
    (
        coeffDict(),
        IOobject
        (
            "sixDoFRigidBodyMotionState",
            mesh.time().timeName(),
            "uniform",
            mesh
        ).typeHeaderOk<IOdictionary>(true)
      ? IOdictionary
        (
            IOobject
            (
                "sixDoFRigidBodyMotionState",
                mesh.time().timeName(),
                "uniform",
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        )
      : coeffDict(),
        mesh.time()
    ),
    patches_(coeffDict().get<wordRes>("patches")),
    patchSet_(mesh.boundaryMesh().patchSet(patches_)),
    di_(coeffDict().get<scalar>("innerDistance")),
    do_(coeffDict().get<scalar>("outerDistance")),
    outerDistanceDomainCheck_(coeffDict().getOrDefault("containOuterDistance",false)),
    xdist_(coeffDict().getOrDefault<scalar>("xDistance",-1)),   // ModMorph: -1 means not used
    ydist_(coeffDict().getOrDefault<scalar>("yDistance",-1)),   // ModMorph: -1 means not used
    test_(coeffDict().getOrDefault("test", false)),
    rhoInf_(1.0),
    rhoName_(coeffDict().getOrDefault<word>("rho", "rho")),
    scale_
    (
        IOobject
        (
            "motionScale",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        pointMesh::New(mesh),
        dimensionedScalar(dimless, Zero)
    ),
    xscale_
    (
        IOobject
        (
            "xmotionScale",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        pointMesh::New(mesh),
        dimensionedScalar(dimless, Zero)
    ),
    yscale_
    (
        IOobject
        (
            "ymotionScale",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        pointMesh::New(mesh),
        dimensionedScalar(dimless, Zero)
    ),
    curTimeIndex_(-1),
    cOfGdisplacement_(coeffDict().getOrDefault<word>("cOfGdisplacement", "none"))
{
    if (rhoName_ == "rhoInf")
    {
        coeffDict().readEntry("rhoInf", rhoInf_);
    }

    // Calculate scaling factor everywhere
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        const pointMesh& pMesh = pointMesh::New(mesh);
        //---ModMorph{
        // Ensure that do_ does not reach beyond the domain. 
        // Note: This is only applicable for rectangular domains aligned with the main coordinate system. 
        tensor scales = tensor::I;
        if (outerDistanceDomainCheck_)
        {
            // Modification to avoid outer distance to go beyond the domain. 
            // Uses boundBox for minimal distance between rectangular domain and body patches. 
            point minPoint = mesh.bounds().min();
            point maxPoint = mesh.bounds().max();
            
            forAll(patches_,pp) // for all patches of the body, check bounding box and compute scale factor. 
            {   
                const label patchi = mesh.boundaryMesh().findPatchID(patches_[pp]);
                const pointField& pf = pMesh.boundary()[patchi].localPoints(); // get local patch points
                boundBox bb(pf); // compute bounding box of current patch pointField 
                point deltaMax = maxPoint-bb.max(); 
                point deltaMin = bb.min()-minPoint;
                // Compute scale factor for different coordinate directions 
                scales.xx() = max(scales.xx(), do_/min(deltaMax.x(), deltaMin.x()));
                scales.yy() = max(scales.yy(), do_/min(deltaMax.y(), deltaMin.y()));
                scales.zz() = max(scales.zz(), do_/min(deltaMax.z(), deltaMin.z()));
            }
        }
        pointPatchDist pDist(pMesh, patchSet_, (scales & points0()));
        //pointPatchDist pDist(pMesh, patchSet_, points0());
        //---ModMorph}
        
        // Scaling: 1 up to di then linear down to 0 at do away from patches
        scale_.primitiveFieldRef() =
            min
            (
                max
                (
                    (do_ - pDist.primitiveField())/(do_ - di_),
                    scalar(0)
                ),
                scalar(1)
            );

        // Convert the scale function to a cosine
        cosineTransition(scale_);   // ModMorph adaptation to shorten code. 
        pointConstraints::New(pMesh).constrain(scale_);
        scale_.write();

        //---ModMorph{
        // Compute x or y scale transitions
        // Note: Put implementation in motion_ to allow for parallel computing.
        if (xdist_>0 || ydist_>0)
        {
            motion_.updateXYScale(points0(), xdist_, ydist_, scale_, xscale_, yscale_);
            if (xdist_>0)
            {   cosineTransition(xscale_);      
                pointConstraints::New(pMesh).constrain(xscale_);
                xscale_.write();
            }
            if (ydist_>0)
            {   cosineTransition(yscale_);     
                pointConstraints::New(pMesh).constrain(yscale_);
                yscale_.write();
            }
        } 
        //---ModMorph}
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::sixDoFRigidBodyMotion&
Foam::sixDoFRigidBodyMotionSolver::motion() const
{
    return motion_;
}


Foam::tmp<Foam::pointField>
Foam::sixDoFRigidBodyMotionSolver::curPoints() const
{
    tmp<pointField> newPoints
    (
        points0() + pointDisplacement_.primitiveField()
    );

    if (!moveAllCells())
    {
        auto ttransformedPts = tmp<pointField>::New(mesh().points());
        auto& transformedPts = ttransformedPts.ref();

        UIndirectList<point>(transformedPts, pointIDs()) =
            pointField(newPoints.ref(), pointIDs());

        return ttransformedPts;
    }

    return newPoints;
}


void Foam::sixDoFRigidBodyMotionSolver::solve()
{
    const Time& t = mesh().time();

    if (mesh().nPoints() != points0().size())
    {
        FatalErrorInFunction
            << "The number of points in the mesh seems to have changed." << endl
            << "In constant/polyMesh there are " << points0().size()
            << " points; in the current mesh there are " << mesh().nPoints()
            << " points." << exit(FatalError);
    }

    // Store the motion state at the beginning of the time-step
    bool firstIter = false;
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        motion_.newTime();
        curTimeIndex_ = this->db().time().timeIndex();
        firstIter = true;
    }

    dimensionedVector g("g", dimAcceleration, Zero);

    if (db().time().foundObject<uniformDimensionedVectorField>("g"))
    {
        g = db().time().lookupObject<uniformDimensionedVectorField>("g");
    }
    else
    {
        coeffDict().readIfPresent("g", g);
    }

    // const scalar ramp = clamp((this->db().time().value() - 5)/10, 0, 1);
    const scalar ramp = 1.0;

    if (test_)
    {
        motion_.update
        (
            firstIter,
            ramp*(motion_.mass()*g.value()),
            ramp*(motion_.mass()*(motion_.momentArm() ^ g.value())),
            t.deltaTValue(),
            t.deltaT0Value()
        );
    }
    else
    {
        dictionary forcesDict;

        forcesDict.add("type", functionObjects::forces::typeName);
        forcesDict.add("patches", patches_);
        forcesDict.add("rhoInf", rhoInf_);
        forcesDict.add("rho", rhoName_);
        forcesDict.add("CofR", motion_.centreOfRotation());

        vector oldPos = motion_.centreOfRotation();

        functionObjects::forces f("forces", db(), forcesDict);

        f.calcForcesMoments();

        motion_.update
        (
            firstIter,
            ramp*(f.forceEff() + motion_.mass()*g.value()),
            ramp
           *(
               f.momentEff()
             + motion_.mass()*(motion_.momentArm() ^ g.value())
            ),
            t.deltaTValue(),
            t.deltaT0Value()
        );

        if (cOfGdisplacement_ != "none")
        {
            if
            (
                db().time().foundObject<uniformDimensionedVectorField>
                (
                    cOfGdisplacement_
                )
            )
            {
                auto& disp =
                    db().time().lookupObjectRef<uniformDimensionedVectorField>
                    (
                        cOfGdisplacement_
                    );

                disp += (motion_.centreOfRotation() - oldPos);
            }
        }
    }

    // Update the displacements
    //--- ModMorph{
    // Update x- or y-compensated morphed positions
    if (xdist_ >0 || ydist_>0)
    {
        pointDisplacement_.primitiveFieldRef() =
            motion_.transform(points0(), xdist_, ydist_, scale_, xscale_, yscale_) - points0();   
    }
    else // do normal morphing operation. 
    {
        pointDisplacement_.primitiveFieldRef() =
            motion_.transform(points0(), scale_) - points0();
    }
    //--- ModMorph}

    // Displacement has changed. Update boundary conditions
    pointConstraints::New
    (
        pointDisplacement_.mesh()
    ).constrainDisplacement(pointDisplacement_);
}


bool Foam::sixDoFRigidBodyMotionSolver::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    IOdictionary dict
    (
        IOobject
        (
            "sixDoFRigidBodyMotionState",
            mesh().time().timeName(),
            "uniform",
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    motion_.state().write(dict);
    return dict.regIOobject::writeObject(streamOpt, writeOnProc);
}


bool Foam::sixDoFRigidBodyMotionSolver::read()
{
    if (displacementMotionSolver::read())
    {
        motion_.read(coeffDict());

        return true;
    }

    return false;
}


// ************************************************************************* //
