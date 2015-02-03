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

#include "QCDAware.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

using namespace std;
using namespace fastjet;


namespace contrib {

    bool operator > (const PJDist& pj1, const PJDist& pj2) {
        return pj1.dist > pj2.dist;
    }


    void QCDAware::insert_pj(ClusterSequence &cs,
            priority_queue<PJDist, vector<PJDist>, greater<PJDist> >& pjds,
            unsigned int iJet,
            vector<bool>& ismerged) const {

        const PseudoJet& ijet = cs.jets()[iJet];

        /*
        cout << "--------" << endl
            << "inserting new pseudojet " << iJet << endl
            << "with flavor " << cs.jets()[iJet].user_index() << endl
            << "pt eta phi " << ijet.pt() << " "
            << ijet.eta() << " " << ijet.phi() << endl;
        */

        for (unsigned int jJet = 0; jJet < iJet; jJet++) {
            // don't calculate distances for already-merged pjs
            if (ismerged[jJet])
                continue;

            const PseudoJet& jjet = cs.jets()[jJet];

            PJDist pjd;
            pjd.pj1 = iJet;
            pjd.pj2 = jJet;
            pjd.dist = _dm->dij(ijet, jjet);
            pjds.push(pjd);

            /*
            cout << "distance to pseudojet " << jJet << " with flavor " <<
                jjet.user_index() << endl
                << "pt eta phi " << jjet.pt() << " "
                << jjet.eta() << " " << jjet.phi()
                << ":" << endl << pjd.dist << endl;
            */
        }

        // calculate the beam distance
        PJDist pjd;
        pjd.pj1 = iJet;
        pjd.pj2 = -1;
        pjd.dist = _dm->diB(ijet);
        pjds.push(pjd);

        /*
        cout << "distance to beam:" << endl
            << pjd.dist << endl;

        cout << "--------" << endl;
        */

        ismerged.push_back(false);

        return;
    }


    void QCDAware::merge_iB(ClusterSequence &cs,
            const PJDist& pjd,
            std::vector<bool>& ismerged) const {

        /*
        cout << "--------" << endl
            << "merging pseudojet " << pjd.pj1 << endl
            << "into the beam" << endl
            << "--------" << endl;
            */

        cs.plugin_record_iB_recombination(pjd.pj1, pjd.dist);

        ismerged[pjd.pj1] = true;

        return;
    }

    void QCDAware::merge_ij(ClusterSequence &cs,
            std::priority_queue<PJDist, std::vector<PJDist>, std::greater<PJDist> >& pjds,
            const PJDist& pjd,
            std::vector<bool>& ismerged) const {

        // mark both old pjs as merged
        ismerged[pjd.pj1] = true;
        ismerged[pjd.pj2] = true;

        const PseudoJet& pj1 = cs.jets()[pjd.pj1];
        const PseudoJet& pj2 = cs.jets()[pjd.pj2];
        PseudoJet pj3 = pj1 + pj2;

        int labi = pj1.user_index();
        int labj = pj2.user_index();

        int abslabi = abs(labi);
        int abslabj = abs(labj);

        // qqbar -> g
        if (abslabi <= 6 && labi + labj == 0)
            pj3.set_user_index(21);
        // gg -> g
        else if (labi == 21 && labj == 21)
            pj3.set_user_index(21);
        // lgamma -> l
        else if ((abslabi == 11 || abslabi == 13) && labj == 22)
            pj3.set_user_index(labi);
        // gammal -> l
        else if ((abslabj == 11 || abslabj == 13) && labi == 22)
            pj3.set_user_index(labj);
        // qg and qgamma -> q
        else if (abslabi <= 6 && (labj == 21 || labj == 22))
            pj3.set_user_index(labi);
        // gq and gammaq -> q
        else if (abslabj <= 6 && (labi == 21 || labi == 22))
            pj3.set_user_index(labj);
        else {
            cout << "ERROR: attempting to merge pseudojets with pdgids "
                << labi << " and " << labj
                << ", which is not allowed: this will probably break." << endl;
            pj3.set_user_index(-999);
        }

        int newidx;
        cs.plugin_record_ij_recombination(pjd.pj1, pjd.pj2, pjd.dist, pj3, newidx);

        insert_pj(cs, pjds, newidx, ismerged);

        /*
        cout << "--------" << endl
            << "merging pseudojets " << pjd.pj1 << " and " << pjd.pj2 << endl
            << "with flavors " << labi << " " << labj << endl
            << "and distance " << pjd.dist << endl
            << "into pseudojet " << newidx << endl
            << "with flavor " << pj3.user_index() << endl
            << "--------" << endl;
        */


        return;
    }


    void QCDAware::run_clustering(ClusterSequence& cs) const {

        vector<bool> ismerged;

        priority_queue<PJDist, vector<PJDist>, greater<PJDist> > pjds;
        for (unsigned int iJet = 0; iJet < cs.jets().size(); iJet++)
            insert_pj(cs, pjds, iJet, ismerged);

        while (!pjds.empty()) {
            PJDist pjd = pjds.top();
            pjds.pop();

            // check for already merged pj1
            if (ismerged[pjd.pj1])
                continue;

            // check for the beam
            if (pjd.pj2 < 0) {
                merge_iB(cs, pjd, ismerged);
                continue;
            }

            // check for already merged pj2
            if (ismerged[pjd.pj2])
                continue;

            merge_ij(cs, pjds, pjd, ismerged);
        }

        return;
    }

    string QCDAware::description() const {
        return "QCDAware Jet Algorithm";
    }

    double QCDAware::R() const {
        return _dm->R();
    }


} // namespace contrib

FASTJET_END_NAMESPACE
