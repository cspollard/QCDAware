//  QCDAware Package
//  Questions/Comments?  abuckley@cern.ch, cpollard@cern.ch
//
//  Copyright (c) 2014
//  Andy Buckley, Chris Pollard, Donatas Zaripovas
//
// $Id$
//
//----------------------------------------------------------------------
// This file is part of FastJet contrib.
//
// It is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// It is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this code. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------

#ifndef __FASTJET_CONTRIB_QCDAWARE_HH__
#define __FASTJET_CONTRIB_QCDAWARE_HH__



#include "fastjet/internal/base.hh"
#include "fastjet/JetDefinition.hh"
#include <queue>
#include <string>
#include <vector>
#include "DistanceMeasure.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


namespace contrib{

    struct PJDist {
        double dist;
        int pj1;
        int pj2;
    };


    //------------------------------------------------------------------------
    class QCDAware : public JetDefinition::Plugin {
        private:
            const DistanceMeasure *_dm;

        void insert_pj(ClusterSequence &cs,
                std::priority_queue<PJDist, std::vector<PJDist>, std::greater<PJDist> >& pjds,
                int iJet,
                std::vector<bool>& ismerged) const;

        void merge_iB(ClusterSequence &cs,
                const PJDist& dist,
                std::vector<bool>& ismerged) const;

        void merge_ij(ClusterSequence &cs,
                std::priority_queue<PJDist, std::vector<PJDist>, std::greater<PJDist> >& pjds,
                const PJDist& dist,
                std::vector<bool>& ismerged) const;

        public:
            /// default constructor
            QCDAware(DistanceMeasure *dm)
                : _dm(dm) {}

            /// default destructor
            ~QCDAware() {}

            void run_clustering(fastjet::ClusterSequence& cs) const;

            std::string description() const;

            double R() const;
    };


} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_QCDAWARE_HH__
