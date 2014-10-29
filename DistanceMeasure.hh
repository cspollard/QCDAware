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


    template <class T>
    class QCDAwareDistanceMeasure : public DistanceMeasure {
        private:
            const T t;

        public:
            QCDAwareDistanceMeasure(double R) : t(R) {}

            inline double dij(const fastjet::PseudoJet& pji, const fastjet::PseudoJet& pjj) const {

                // label = pdgid
                int labi = pji.user_index();
                int labj = pjj.user_index();

                if (abs(labi) > 6)
                    std::cout << "found a QCDAware jet constituent with label = " <<
                        labi << ".\nExiting." << std::endl;

                if (abs(labj) > 6)
                    std::cout << "found a QCDAware jet constituent with label = " <<
                        labj << ".\nExiting." << std::endl;


                double dist = t.dij(pji, pjj);
                /*
                std::cout << "calculating qcdaware distance: " << std::endl;
                std::cout << "original distance: " << dist << std::endl;
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


                // not allowed: qq or qbarqbar clustering
                if (labi * labj > 0) {
                    // std::cout << "attempting qq or qbarqbar clustering..." << std::endl;
                    dist = 999*t.diB(pji);
                // not allowed: different flavor clustering
                } else if (labi && labj && (labi + labj)) {
                    // std::cout << "attempting to cluster quarks with different flavors..." << std::endl;
                    dist = 999*t.diB(pji);
                }

                // std::cout << "final distance: " << dist << std::endl;

                return dist;
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
