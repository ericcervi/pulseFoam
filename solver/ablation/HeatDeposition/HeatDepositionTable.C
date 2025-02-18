/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "HeatDepositionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// template<class Type>
Foam::Function1Types::HeatDepositionTable::HeatDepositionTable
(
    const word& entryName,
    const dictionary& dict,
    scalar q0,
    const objectRegistry* obrPtr
)
:
    Table<scalar>(entryName, dict.subDict(entryName), obrPtr),
    explosiveBoilingThreshold_(dict.getOrDefault<scalar>("qThreshold" , 1e10)),
    isIncreasing_(true)
{
    // scale table
    scale(q0);

    // I want to create a twin reversed table to perform reverse lookup (if needed)
    // If the function is not increasing I allocate memory for another table, otherwise
    // I just use the original table.

    // Check the first two elements only assuming the function in the table is monotonic
    isIncreasing_ = this->table_[0].second() > this->table_[1].second() ? false : true;  

    // If not increasing (assuming it's monotonic) create a reversed table
    if (!isIncreasing_)
    {
        revTablePtr_.reset(new List<Tuple2<scalar, scalar>>(this->table_.size()));
        List<Tuple2<scalar, scalar>>& revTable = revTablePtr_();
        
        // forward label
        label j = 0;
        forAllReverse(this->table_, i)
        {
            revTable[j].first() = this->table_[i].first();
            revTable[j].second() = this->table_[i].second();
            j++;
        }
    }
}


// template<class Type>
// Foam::Function1Types::Table<Type>::Table(const Table<Type>& tbl)
// :
//     Table<Type>(tbl),
//     fName_(tbl.fName_)
// {}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::interpolationWeights&
Foam::Function1Types::HeatDepositionTable::interpolator() const
{
    // I need to check if the function is monotonically increasing or decreasing
    // If the dependent variable is decreasing I need to reverse the order of the table_
    // so that samples (aka *tableSamplePtr_) is populated 

    // When the interpolator is called, the existence of the pointer is checked first
    // If it doesn't exist the interpolator is defined with a selector
    // The interpolator needs to be feeded with a list of the dependent variable, in this case f(x)
    // If the function is not decreasing the samples are populated using the flipped table

    if (!revInterpolatorPtr_)
    { 
        // Re-work table into linear list

        // If function is not increasing populate the reverse lookup samples list
        // with elements from the flipped table

        revTableSamplesPtr_.reset(new scalarField(table_.size()));
        auto& revSamples = *revTableSamplesPtr_;

        if (isIncreasing_)
        {
            forAll(this->table_, i)
            {
                revSamples[i] = this->table_[i].second();
            }
            revInterpolatorPtr_ = interpolationWeights::New
            (
                interpolationScheme_,
                revSamples
            );

        }
        else
        {
            forAll(revTablePtr_(), i)
            {
                revSamples[i] = revTablePtr_()[i].second();
            }
            revInterpolatorPtr_ = interpolationWeights::New
            (
                interpolationScheme_,
                revSamples
            );
        }
        
    }

    return *revInterpolatorPtr_;
}


Foam::scalar Foam::Function1Types::HeatDepositionTable::reverseValue(const scalar f) const
{
    //scalar yDash = f;

    // if (checkMinBounds(f, yDash))
    // {
    //     return revTablePtr_().first().first();
    // }

    // if (checkMaxBounds(yDash, yDash))
    // {
    //     return revTablePtr_().last().first();
    // }

    // Use interpolator
    //interpolator().valueWeights(yDash, revCurrentIndices_, revCurrentWeights_);

    // The interpolator knows whether or not the function in the table is increasing or decreasing
    // so the weights should be correctly set
    interpolator().valueWeights(f, revCurrentIndices_, revCurrentWeights_);

    if (isIncreasing_)
    {

        scalar t(revCurrentWeights_[0]*this->table_[revCurrentIndices_[0]].first());
        for (label i = 1; i < revCurrentIndices_.size(); i++)
        {
            t += revCurrentWeights_[i]*this->table_[revCurrentIndices_[i]].first();
        }

        return t;
    }
    else
    {
        scalar t(revCurrentWeights_[0]*revTablePtr_()[revCurrentIndices_[0]].first());
        for (label i = 1; i < revCurrentIndices_.size(); i++)
        {
            t += revCurrentWeights_[i]*revTablePtr_()[revCurrentIndices_[i]].first();
        }
        
        return t;
    }
}

void Foam::Function1Types::HeatDepositionTable::scale(const scalar qVol)
{
    forAll(this->table_, i)
    {
        this->table_[i].second() *= qVol;
    }
}
// ************************************************************************* //
