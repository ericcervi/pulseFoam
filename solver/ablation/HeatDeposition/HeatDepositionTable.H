/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Class
    Foam::Function1Types::Table

Description
    Table class for heat deposition profile.
    Derived from the OpenFOAM Table class with additional
    support for reverse lookup.
    

    Dictionary specification, embedded content
    Dictionary form.
    \verbatim
        <entryName>
        {
            type table;
            values
            (
                (0.0 (1 2 3))
                (1.0 (4 5 6))
            );
        }
    \endverbatim

Note
    The external data reference (using the \c file keyword) is
    used in preference to the \c values specification.

See Also
    Foam::Function1Types::TableFile

SourceFiles
    HeatDepositionTable.C

\*---------------------------------------------------------------------------*/

#ifndef HeatDepositionTable_H
#define HeatDepositionTable_H

#include "Table.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                           Class Table Declaration
\*---------------------------------------------------------------------------*/


// No need to template the class
class HeatDepositionTable
:
    public Table<scalar>
{
    // Private Data

        //- Input name for file-based input (optional)
        //fileName fName_;

    const scalar explosiveBoilingThreshold_;


    // Private Member Functions

        //- No copy assignment
        //void operator=(const Table<Type>&) = delete;

protected:

    // One way to handle the increasing/decreasing nature of the function 
    // is to create a duplicate table only if the dependent variable is decreasing
    // I could define a pointer (or better a smart pointer) to a table that I set as nullptr unless the function
    // is decreasing. If the function is decreasing I create a reverse table reading the original one
    // and constructing from the bottom. I could do this in the constructor.

    //- Check if the table function is increasing. If not, a reversed table needs to be created
        bool isIncreasing_;

        autoPtr<List<Tuple2<scalar, scalar>>>  revTablePtr_;
     
     //- Extracted values from table. Only dependent variable
        mutable autoPtr<scalarField> revTableSamplesPtr_;

     //- Interpolator method for reverse interpolation
        mutable autoPtr<interpolationWeights> revInterpolatorPtr_;

         //- Cached indices and weights
        mutable labelList revCurrentIndices_;

        mutable scalarField revCurrentWeights_;


    // Protected Member Functions

         //- Return (demand driven) interpolator
        const interpolationWeights& interpolator() const;

public:

    //- Runtime type information
    //TypeName("heatDepositionTable");


    // Constructors

        //- Construct from entry name, dictionary and optional registry
        // q0 is the scaling factor for the heat deposition profile.
        HeatDepositionTable
        (
            const word& entryName,
            const dictionary& dict,
            scalar q0 = 1.0,
            const objectRegistry* obrPtr = nullptr
        );

        //- Copy construct
        //explicit HeatDepositionTable(const Table<Type>& tbl);

        //- Construct and return a clone
        virtual tmp<Function1<scalar>> clone() const
        {
            return tmp<Function1<scalar>>(new HeatDepositionTable(*this));
        }


    //- Destructor
    virtual ~HeatDepositionTable() = default;


    // Member Functions

     //- Return Table value
        scalar reverseValue(const scalar f) const;


    //- Scale table
        void scale(const scalar qVol);
        //- Write coefficients in dictionary format
        //virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
//     #include "Table.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
