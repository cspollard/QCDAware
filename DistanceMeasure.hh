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
}

FASTJET_END_NAMESPACE

#endif
