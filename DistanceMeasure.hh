#ifndef __DISTANCEMEASURE_HH__
#define __DISTANCEMEASURE_HH__

#include <float.h>
#include <algorithm>
#include "fastjet/PseudoJet.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {

    class DistanceMeasure {
        public:
            virtual double dij(const fastjet::PseudoJet& pji, const fastjet::PseudoJet& pjj) const = 0;
            virtual double diB(const fastjet::PseudoJet& pji) const = 0;
            virtual double R() const = 0;
    };


    class KtMeasure : public DistanceMeasure {
        private:
            double _R;

        public:
            KtMeasure(double R) : _R(R) {}

            inline double dij(const fastjet::PseudoJet& pji, const fastjet::PseudoJet& pjj) const {
                double drbyR = pji.delta_R(pjj) / _R;
                return std::min(pji.perp2(), pjj.perp2()) * drbyR * drbyR;
            }

            inline double diB(const fastjet::PseudoJet& pji) const {
                return pji.perp2();
            }

            virtual double R() const {
                return _R;
            }
    };


    class AntiKtMeasure : public DistanceMeasure {
        private:
            double _R;

        public:
            AntiKtMeasure(double R) : _R(R) {}

            inline double dij(const fastjet::PseudoJet& pji, const fastjet::PseudoJet& pjj) const {
                double drbyR = pji.delta_R(pjj) / _R;
                return 1.0 / std::max(pji.perp2(), pjj.perp2()) * drbyR * drbyR;
            }

            inline double diB(const fastjet::PseudoJet& pji) const {
                return 1.0 / pji.perp2();
            }

            virtual double R() const {
                return _R;
            }
    };


    class CAMeasure : public DistanceMeasure {
        private:
            double _R;

        public:
            CAMeasure(double R) : _R(R) {}

            inline double dij(const fastjet::PseudoJet& pji, const fastjet::PseudoJet& pjj) const {
                double drbyR = pji.delta_R(pjj) / _R;
                return drbyR * drbyR;
            }


            inline double diB(const fastjet::PseudoJet& pji) const {
                return 1.0;
            }

            virtual double R() const {
                return _R;
            }
    };

    // helper functions
    namespace {
        inline bool isQuark(const fastjet::PseudoJet& p) {
            return abs(p.user_index()) <= 6;
        }

        inline bool isGluon(const fastjet::PseudoJet& p) {
            return p.user_index() == 21;
        }

        inline bool isPhoton(const fastjet::PseudoJet& p) {
            return p.user_index() == 22;
        }

        inline bool canCluster(const fastjet::PseudoJet& p, const fastjet::PseudoJet& q) {
            // a quark can cluster with a photon or gluon.
            if (    // kind of complicated, so indented
                    (isQuark(p) && (isGluon(q) || isPhoton(q))) ||
                    ((isGluon(p) || isPhoton(p)) && isQuark(q))
                )
                return true;

            // gluons can cluster.
            else if (isGluon(p) && isGluon(q))
                return true;

            // same-flavor quark and anti-quark can cluster.
            else if (isQuark(p) && isQuark(q) &&
                    (p.user_index() + q.user_index() == 0))
                return true;

            // nothing else allowed. (for now... muahahaha!)
            return false;
        }
    }


    template <class T>
    class QCDAwareDistanceMeasure : public DistanceMeasure {
        private:
            const T t;

        public:
            QCDAwareDistanceMeasure(double R) : t(R) {}

            inline double dij(const fastjet::PseudoJet& pji, const fastjet::PseudoJet& pjj) const {

                /*
                std::cout << "calculating qcdaware distance: " << std::endl;
                std::cout << "original distance: " << t.dij(pji, pjj) << std::endl;
                std::cout << "pseudojet i pt eta phi lab : "
                    << pji.pt() << " "
                    << pji.eta() << " "
                    << pji.phi() << " "
                    << pji.user_index() << std::endl;

                std::cout << "pseudojet j pt eta phi lab : "
                    << pjj.pt() << " "
                    << pjj.eta() << " "
                    << pjj.phi() << " "
                    << pjj.user_index() << std::endl;
                */


                if (canCluster(pji, pjj))
                    return t.dij(pji, pjj);
                else
                    return DBL_MAX;
            }

            inline double diB(const fastjet::PseudoJet& pji) const {
                return t.diB(pji);
            }

            virtual double R() const {
                return t.R();
            }
    };


}

FASTJET_END_NAMESPACE

#endif
