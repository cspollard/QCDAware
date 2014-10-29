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

                if (labi == 21 || labi == 22)
                    labi = 0;

                if (labj == 21 || labj == 22)
                    labj = 0;

                if (abs(labi) > 6)
                    std::cout << "found a QCDAware jet constituent with pid = " <<
                        labi << ". Exiting." << std::endl;

                if (abs(labj) > 6)
                    std::cout << "found a QCDAware jet constituent with pid = " <<
                        labj << ". Exiting." << std::endl;


                // not allowed: qq or qbarqbar clustering
                if (labi * labj > 0)
                    return 999*t.diB(pji);

                // not allowed: different flavor clustering
                if (labi && labj && (labi + labj))
                    return 999*t.diB(pji);

                return t.dij(pji, pjj);
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
