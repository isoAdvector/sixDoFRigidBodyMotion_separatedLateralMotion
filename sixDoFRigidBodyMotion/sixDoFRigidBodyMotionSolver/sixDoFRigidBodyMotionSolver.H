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

Class
    Foam::sixDoFRigidBodyMotionSolver

Group
   grpMeshMotionSolvers

Description
    6-DoF solid-body mesh motion solver for an fvMesh.

    Applies SLERP interpolation of movement as function of distance to
    the object surface.

SourceFiles
    sixDoFRigidBodyMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef sixDoFRigidBodyMotionSolver_H
#define sixDoFRigidBodyMotionSolver_H

#include "displacementMotionSolver.H"
#include "sixDoFRigidBodyMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class sixDoFRigidBodyMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class sixDoFRigidBodyMotionSolver
:
    public displacementMotionSolver
{
    // Private data

        //- Six DoF motion object
        sixDoFRigidBodyMotion motion_;

        wordRes patches_;

        //- Patches to integrate forces
        const labelHashSet patchSet_;

        //- Inner morphing distance (limit of solid-body region)
        const scalar di_;

        //- Outer morphing distance (limit of linear interpolation region)
        const scalar do_;

        //---MorphMod{ 
        //- Switch to use the domain check for the outer distance
        bool outerDistanceDomainCheck_; 
        
        //- Morphing distances of x- and y-directions respectively (length of linear interpolation region)
        const scalar xdist_;

        const scalar ydist_;
        //---MorphMod}

        //- Test-mode in which only the gravitational body-force is applied
        bool test_;

        //- Reference density required by the forces object for
        //  incompressible calculations, required if rho == rhoInf
        scalar rhoInf_;

        //- Name of density field, optional unless used for an
        //  incompressible simulation, when this needs to be specified
        //  as rhoInf
        word rhoName_;

        //- Current interpolation scale (1 at patches, 0 at distance_)
        pointScalarField scale_;

        //---MorphMod{
        //- X-coordinate scale field (1 is rigid body transform, 0 is static mesh) 
        pointScalarField xscale_; 

        //- Y-coordinate scale field (1 is rigid body transform, 0 is static mesh)
        pointScalarField yscale_;
        //---MorphMod}

        //- Current time index (used for updating)
        label curTimeIndex_;

        //- Name of the uniformVectorField for CofG displacement
        word cOfGdisplacement_;


    // Private Member Functions

        //- No copy construct
        sixDoFRigidBodyMotionSolver
        (
            const sixDoFRigidBodyMotionSolver&
        ) = delete;

        //- No copy assignment
        void operator=(const sixDoFRigidBodyMotionSolver&) = delete;

        //---MorphMod{
        //- 
        void cosineTransition(pointScalarField& scale);
        
        //---MorphMod}

public:

    //- Runtime type information
    TypeName("sixDoFRigidBodyMotion");


    // Constructors

        //- Construct from polyMesh and IOdictionary
        sixDoFRigidBodyMotionSolver
        (
            const polyMesh&,
            const IOdictionary& dict
        );


    //- Destructor
    ~sixDoFRigidBodyMotionSolver() = default;


    // Member Functions

        //- Return the six DoF motion object
        const sixDoFRigidBodyMotion& motion() const;

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Write state using stream options
        virtual bool writeObject
        (
            IOstreamOption streamOpt,
            const bool writeOnProc
        ) const;

        //- Read dynamicMeshDict dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
