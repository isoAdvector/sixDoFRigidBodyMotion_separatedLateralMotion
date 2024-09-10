/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

#include "sixDoFRigidBodyMotion.H"
#include "sixDoFSolver.H"
#include "septernion.H"
//---ModMorph{
#include "boundBox.H"
//---ModMorph}
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotion::applyRestraints()
{
    if (restraints_.empty())
    {
        return;
    }

    if (Pstream::master())
    {
        forAll(restraints_, rI)
        {
            if (report_)
            {
                Info<< "Restraint " << restraints_[rI].name() << ": ";
            }

            // Restraint position
            point rP = Zero;

            // Restraint force
            vector rF = Zero;

            // Restraint moment
            vector rM = Zero;

            // Accumulate the restraints
            restraints_[rI].restrain(*this, rP, rF, rM);

            // Update the acceleration
            a() += rF/mass_;

            // Moments are returned in global axes, transforming to
            // body local to add to torque.
            tau() += Q().T() & (rM + ((rP - centreOfRotation()) ^ rF));
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotion::sixDoFRigidBodyMotion(const Time& time)
:
    time_(time),
    motionState_(),
    motionState0_(),
    restraints_(),
    constraints_(),
    tConstraints_(tensor::I),
    rConstraints_(tensor::I),
    initialCentreOfMass_(Zero),
    initialCentreOfRotation_(Zero),
    initialQ_(I),
    mass_(VSMALL),
    momentOfInertia_(diagTensor::one*VSMALL),
    aRelax_(1.0),
    aDamp_(1.0),
    report_(false),
    updateConstraints_(false),
    solver_(nullptr)
{}


Foam::sixDoFRigidBodyMotion::sixDoFRigidBodyMotion
(
    const dictionary& dict,
    const dictionary& stateDict,
    const Time& time
)
:
    time_(time),
    motionState_(stateDict),
    motionState0_(),
    restraints_(),
    constraints_(),
    tConstraints_(tensor::I),
    rConstraints_(tensor::I),
    initialCentreOfMass_
    (
        dict.getOrDefault
        (
            "initialCentreOfMass",
            dict.get<vector>("centreOfMass")
        )
    ),
    initialCentreOfRotation_(initialCentreOfMass_),
    initialQ_
    (
        dict.getOrDefault
        (
            "initialOrientation",
            dict.getOrDefault("orientation", tensor::I)
        )
    ),
    mass_(dict.get<scalar>("mass")),
    momentOfInertia_(dict.get<diagTensor>("momentOfInertia")),
    aRelax_(dict.getOrDefault<scalar>("accelerationRelaxation", 1)),
    aDamp_(dict.getOrDefault<scalar>("accelerationDamping", 1)),
    report_(dict.getOrDefault("report", false)),
    updateConstraints_(dict.getOrDefault("updateConstraints", false)),
    solver_(sixDoFSolver::New(dict.subDict("solver"), *this))
{
    addRestraints(dict);

    // Set constraints and initial centre of rotation
    // if different to the centre of mass
    addConstraints(dict);

    // If the centres of mass and rotation are different ...
    vector R(initialCentreOfMass_ - initialCentreOfRotation_);
    if (magSqr(R) > VSMALL)
    {
        // ... correct the moment of inertia tensor using parallel axes theorem
        momentOfInertia_ += mass_*diag(I*magSqr(R) - sqr(R));

        // ... and if the centre of rotation is not specified for motion state
        // update it
        if (!stateDict.found("centreOfRotation"))
        {
            motionState_.centreOfRotation() = initialCentreOfRotation_;
        }
    }

    // Save the old-time motion state
    motionState0_ = motionState_;
}


Foam::sixDoFRigidBodyMotion::sixDoFRigidBodyMotion
(
    const sixDoFRigidBodyMotion& sDoFRBM
)
:
    time_(sDoFRBM.time_),
    motionState_(sDoFRBM.motionState_),
    motionState0_(sDoFRBM.motionState0_),
    restraints_(sDoFRBM.restraints_),
    constraints_(sDoFRBM.constraints_),
    tConstraints_(sDoFRBM.tConstraints_),
    rConstraints_(sDoFRBM.rConstraints_),
    initialCentreOfMass_(sDoFRBM.initialCentreOfMass_),
    initialCentreOfRotation_(sDoFRBM.initialCentreOfRotation_),
    initialQ_(sDoFRBM.initialQ_),
    mass_(sDoFRBM.mass_),
    momentOfInertia_(sDoFRBM.momentOfInertia_),
    aRelax_(sDoFRBM.aRelax_),
    aDamp_(sDoFRBM.aDamp_),
    report_(sDoFRBM.report_),
    updateConstraints_(sDoFRBM.updateConstraints_),
    solver_(sDoFRBM.solver_.clone())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotion::~sixDoFRigidBodyMotion()
{} // Define here (incomplete type in header)


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotion::addRestraints
(
    const dictionary& dict
)
{
    if (dict.found("restraints"))
    {
        const dictionary& restraintDict = dict.subDict("restraints");

        label i = 0;

        restraints_.setSize(restraintDict.size());

        for (const entry& dEntry : restraintDict)
        {
            if (dEntry.isDict())
            {
                restraints_.set
                (
                    i++,
                    sixDoFRigidBodyMotionRestraint::New
                    (
                        dEntry.keyword(),
                        dEntry.dict()
                    )
                );
            }
        }

        restraints_.setSize(i);
    }
}


void Foam::sixDoFRigidBodyMotion::addConstraints
(
    const dictionary& dict
)
{
    if (dict.found("constraints"))
    {
        const dictionary& constraintDict = dict.subDict("constraints");

        label i = 0;

        constraints_.setSize(constraintDict.size());

        pointConstraint pct;
        pointConstraint pcr;

        for (const entry& dEntry : constraintDict)
        {
            if (dEntry.isDict())
            {
                constraints_.set
                (
                    i,
                    sixDoFRigidBodyMotionConstraint::New
                    (
                        dEntry.keyword(),
                        dEntry.dict(),
                        *this
                    )
                );

                constraints_[i].setCentreOfRotation(initialCentreOfRotation_);
                constraints_[i].constrainTranslation(pct);
                constraints_[i].constrainRotation(pcr);

                i++;
            }
        }

        constraints_.setSize(i);

        tConstraints_ = pct.constraintTransformation();
        rConstraints_ = pcr.constraintTransformation();

        Info<< "Translational constraint tensor " << tConstraints_ << nl
            << "Rotational constraint tensor " << rConstraints_ << endl;
    }
}


void Foam::sixDoFRigidBodyMotion::updateAcceleration
(
    const vector& fGlobal,
    const vector& tauGlobal
)
{
    static bool first = true;

    // Save the previous iteration accelerations for relaxation
    vector aPrevIter = a();
    vector tauPrevIter = tau();

    // Calculate new accelerations
    a() = fGlobal/mass_;
    tau() = (Q().T() & tauGlobal);
    applyRestraints();

    // Relax accelerations on all but first iteration
    if (!first)
    {
        a() = aRelax_*a() + (1 - aRelax_)*aPrevIter;
        tau() = aRelax_*tau() + (1 - aRelax_)*tauPrevIter;
    }
    else
    {
        first = false;
    }
}


void Foam::sixDoFRigidBodyMotion::updateConstraints()
{
    if (!updateConstraints_)
    {
        return;
    }

    pointConstraint pct;
    pointConstraint pcr;

    forAll(constraints_, i)
    {
        constraints_[i].setCentreOfRotation(initialCentreOfRotation_);
        constraints_[i].constrainTranslation(pct);
        constraints_[i].constrainRotation(pcr);
    }

    tConstraints_ = pct.constraintTransformation();
    rConstraints_ = pcr.constraintTransformation();

    Info<< "Translational constraint tensor " << tConstraints_ << nl
        << "Rotational constraint tensor " << rConstraints_ << endl;
}


void Foam::sixDoFRigidBodyMotion::update
(
    bool firstIter,
    const vector& fGlobal,
    const vector& tauGlobal,
    scalar deltaT,
    scalar deltaT0
)
{
    if (Pstream::master())
    {
        solver_->solve(firstIter, fGlobal, tauGlobal, deltaT, deltaT0);

        if (report_)
        {
            status();
        }
    }

    Pstream::broadcast(motionState_);
}


void Foam::sixDoFRigidBodyMotion::status() const
{
    Info<< "6-DoF rigid body motion" << nl
        << "    Centre of rotation: " << centreOfRotation() << nl
        << "    Centre of mass: " << centreOfMass() << nl
        << "    Orientation: " << orientation() << nl
        << "    Linear velocity: " << v() << nl
        << "    Angular velocity: " << omega()
        << endl;
}


Foam::tmp<Foam::pointField> Foam::sixDoFRigidBodyMotion::transform
(
    const pointField& initialPoints
) const
{
    return
    (
        centreOfRotation()
      + (Q() & initialQ_.T() & (initialPoints - initialCentreOfRotation_))
    );
}


Foam::tmp<Foam::pointField> Foam::sixDoFRigidBodyMotion::transform
(
    const pointField& initialPoints,
    const scalarField& scale
) const
{
    // Calculate the transformation septerion from the initial state
    septernion s
    (
        centreOfRotation() - initialCentreOfRotation(),
        quaternion(Q().T() & initialQ())
    );

    auto tpoints = tmp<pointField>::New(initialPoints);
    auto& points = tpoints.ref();

    forAll(points, pointi)
    {
        // Move non-stationary points
        if (scale[pointi] > SMALL)
        {
            // Use solid-body motion where scale = 1
            if (scale[pointi] > 1 - SMALL)
            {
                points[pointi] = transform(initialPoints[pointi]);
            }
            // Slerp septernion interpolation
            else
            {
                septernion ss(slerp(septernion::I, s, scale[pointi]));

                points[pointi] =
                    initialCentreOfRotation()
                  + ss.invTransformPoint
                    (
                        initialPoints[pointi]
                      - initialCentreOfRotation()
                    );
            }
        }
    }

    return tpoints;
}


//--- Moody{ 
// New transform method based on two scales. One for surge, and one for 5 dof motion.
Foam::tmp<Foam::pointField> Foam::sixDoFRigidBodyMotion::transform
(
    const pointField& initialPoints,
    const scalar& xdist, 
    const scalar& ydist,
    const scalarField& scale,
    const scalarField& xScale,
    const scalarField& yScale   
) const
{

    //- Get switches for different directions. True means active translation region -//
    const bool isXScale = xdist>0;
    const bool isYScale = ydist>0;

    // Compute translation displacement. 
    // tPoint- point to morph in translation regions of tDist
    // slerpPoint - point offset to morph with slerp region
    point tPoint = centreOfRotation() - initialCentreOfRotation();
    // Create slerp point as origin of rotation slerp. Remove translation if present as switches. 
    point slerpPoint(tPoint);
	if (isXScale)
        slerpPoint.x() = 0;
    if (isYScale)
        slerpPoint.y() = 0;

    // Calculate the transformation septerion from the initial state for the slerp dofs.  
    septernion s
    (
        slerpPoint,
        quaternion(Q().T() & initialQ())
    );
       
    // New pointfield for results
    tmp<pointField> tpoints(new pointField(initialPoints));
    pointField& points = tpoints.ref();


    forAll(points, pointi)
    {
        // Move non-stationary points in centre region
        if (scale[pointi] > SMALL)
        {
	    septernion ss(s);
	    if (scale[pointi] <= 1 - SMALL)			
		ss= slerp(septernion::I, s, scale[pointi]);

            points[pointi] =
                    initialCentreOfRotation()
                  + ss.invTransformPoint
                    (
                        initialPoints[pointi]
                      - initialCentreOfRotation()
                    );        
        }
        // Add x- and y-scale translations to the point location 
        if (xScale[pointi]>SMALL)
        {   points[pointi].x() += xScale[pointi]*tPoint.x();   }
        if (yScale[pointi]>SMALL)
        {   points[pointi].y() += yScale[pointi]*tPoint.y();   }
    
    }

    return tpoints;
}

// New method to update x and y-scale.
void Foam::sixDoFRigidBodyMotion::updateXYScale
(
    const pointField& initialPoints,
    const scalar& xdist,
    const scalar& ydist,
    const scalarField& scale,
    scalarField& xScale,
    scalarField& yScale
) const
{
	//- Collect initial points as copy of points
	tmp<pointField> tpoints(new pointField(initialPoints));
	pointField& points = tpoints.ref();

    // Find the bounding box of the inner (or outer) distance from scale	
	//- Move all points to inside the slerp scale domain. Then compute boundBox
	forAll(points,pointi)
	{
		if (scale[pointi] <= 1-SMALL)
			points[pointi] = initialCentreOfRotation(); 
	}
	
	//- Use boundBox to get min and max x-value of the deformation sphere of scale.
	boundBox box(points);
	vector minVal = box.min();
	vector maxVal = box.max();

	// Compute bound box of whole domain, stable for rectangular box domains. 
	boundBox boxi(initialPoints);
	vector minDomain = boxi.min();
	vector maxDomain = boxi.max();

	// Compute domain-adjusted interpolation lengths dx and dy
	scalar dx = min( min(minVal.x()-minDomain.x(),	maxDomain.x()-maxVal.x() )
			, xdist);	
	scalar dy = min( min(minVal.y()-minDomain.y(),	maxDomain.y()-maxVal.y() )
			, ydist);

	if (dx > SMALL)
		//- Set xScale values based on minVal and maxVal 
		forAll(initialPoints,pointi)
		{
			// Shorthand notation:
			const scalar& xVal = initialPoints[pointi].x();
						
			// Compute x-scale on right side of bound-box.
		    	if ( xVal >= maxVal.x() )
				xScale[pointi]= max( 1.0 - (xVal-maxVal.x())/dx, 0.0 );
					
			else if ( xVal <= minVal.x() ) // left side of bound box
					
				xScale[pointi]= max( 1.0 - (minVal.x()-xVal)/dx, 0.0 );
					
			else // inside the bound box (x-wise), use rigid body x-motion
				xScale[pointi]= 1.0;
	}

	// Repeat for y-scale
	if ( dy > SMALL ) 
    	{
		forAll(initialPoints,pointi)
		{
			const scalar& yVal = initialPoints[pointi].y();
			if ( yVal >= maxVal.y() ) 
				yScale[pointi]= max( 1.0 - (yVal-maxVal.y())/dy, 0.0 );
					
			else if ( yVal <= minVal.y() )				
				yScale[pointi]= max( 1.0 - (minVal.y()-yVal)/dy, 0.0 );
					
			else 
				yScale[pointi]= 1.0;			

		}
		
		// If x-scale is used, multiply to avoid moving relaxation zones at x=start and end
		if (dx > SMALL)
			yScale *= xScale;
	}
		
	return;
}

// Update method ends
//--- Moody}

// ************************************************************************* //
